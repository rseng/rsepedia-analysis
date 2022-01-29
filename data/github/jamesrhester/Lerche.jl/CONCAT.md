[![DOI](https://zenodo.org/badge/226250636.svg)](https://zenodo.org/badge/latestdoi/226250636)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.03497/status.svg)](https://doi.org/10.21105/joss.03497)
![Testing](https://github.com/jamesrhester/Lerche.jl/workflows/CI/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/jamesrhester/Lerche.jl/badge.svg?branch=master)](https://coveralls.io/github/jamesrhester/Lerche.jl?branch=master)
# Introduction

Lerche (German for Lark) is a partial port of the [Lark grammar processor](https://github.com/erezsh/lark-parser/lark) from
Python to Julia.  Lark grammars should work unchanged in Lerche.

**Installation**: at the Julia REPL, `using Pkg; Pkg.add("Lerche")`

**Documentation**: [![][docs-latest-img]][docs-latest-url]

# Quick start

See also 'Notes for Lark users' below.

Lerche reads Lark EBNF grammars to produce a parser. This parser, when
provided with text conforming to the grammar, produces a parse
tree. This tree can be visited and transformed using "rules". A rule is
a function named after the production whose arguments it should be called on, and
the first argument of a rule is an object which is a subtype of
``Visitor`` or ``Transformer``.

Given an EBNF grammar, it can be used to parse text into your data
structure as follows:
1. Define one or more subtypes of ``Transformer`` or ``Visitor`` instances of which will be
passed as the first argument to the appropriate rule. The instance can also be used to
hold information during transformation if you wish, in which case it must have a concrete type.
1. Define `visit_tokens(t::MyNewType) = false` if you will not be processing token values. This
is about 25% faster than leaving the default `true`.
1. For every production in your grammar that you wish to process,
write a rule with identical name to the production
1. The rule should be prefixed with macro ``@rule`` if the second argument
is an array containing all of the arguments to the grammar production
1. The rule should be prefixed with macro ``@inline_rule`` if the second
and following arguments refer to each argument in the grammar production
1. For every token which you wish to process, define an identically-named method
as for rules, but precede it with a ``@terminal`` macro instead of ``@rule``.

If your grammar is in ``String`` variable ``mygrammar``, your text to be parsed and transformed
is in ``String`` variable ``mytext``, and your ``Transformer`` subtype is ``MyTransformer``, the
following commands will produce a data structure from the text:

```julia
using Lerche
p = Lark(mygrammar,parser="lalr",lexer="contextual") #create parser
t = Lerche.parse(p,mytext)     #Create parse tree
x = Lerche.transform(MyTransformer(),t)  #transform parse tree
```

For a real-world example of usage, see [this file](https://github.com/jamesrhester/DrelTools.jl/blob/master/src/jl_transformer.jl).

## Citation

If you are publishing work where Lerche has been useful, please consider citing [the Lerche paper](https://doi.org/10.21105/joss.03497).

# Issues

Please raise any issues or problems with using Lerche in the [Github
issue tracker](https://github.com/jamesrhester/Lerche.jl/issues).

# Contributions

Contributions of all types are welcome. Examples include:
* Improvements to processing speed
* Improved documentation
* Links to projects using Lerche
* Commenting and triaging issues

The most straightforward way to make a contribution is to fork the
repository, make your changes, and create a pull request.

# Notes for Lark users

Please read the Lark documentation.  When converting from Lark
programs written in Python to Lerche programs written in Julia, the
changes outlined below are necessary.

1. All Transformer and Visitor classes become subtypes of Transformer/Visitor
1. All class method calls become Julia method calls with an instance of the type as the first argument
(i.e. replacing ``self``)
1. Transformation or visitor rules should be preceded by the ``@rule`` macro. Inline
rules use the ``@inline_rule`` macro and token processing methods use ``@terminal``. 
1. The first argument of transformer and visitor rules is a variable of the
desired transformer/visitor type.
1. Any grammars containing backslash-double quote sequences need to be fixed (see below).
1. Any grammars containing backslash-x to denote a byte value need to be fixed (see below).

## Inconsistencies with Lark

1. Earley and CYK grammars are not implemented. 
2. Dynamic lexer is not implemented. 
3. All errors with messages attached must be at the bottom of the
exception type hierarchy, as these are the only types that can have
contents. Thus an ``UnexpectedInput`` exception must become e.g 
an ``UnexpectedCharacter`` exception if a message is included.
4. The `PuppetParser` invoked when there is a parse error is not yet
functional
5. There may be issues with correctly interpreting import paths
to find imported grammars: please raise an issue if this happens.
6. No choice of ``regex`` engine, ``Tree`` structure or byte/string
choices are available as they make no sense for Julia.

# Implementation notes and hints

Lerche is currently based off Lark 0.11.1. The priority has been on
maintaining fidelity with Lark. For example, global `regex` flags
which are integers in Lark are still integers in Lerche, which means
you will need to look their values up. This may be changed to a more
Julian approach in future.

The ``@rule`` and ``@inline_rule`` macros define methods of Lerche function
`transformer_func`. Julia multiple dispatch is used to select the
appropriate method at runtime. ``@terminal`` similarly defines methods
of ``token_func``.

Parsing a large (500K) file suggest Lerche is about 3 times faster
than Lark with CPython for parsing. Parser generation is much slower as no
optimisation techniques have been applied (yet). Calculating and
storing your grammar in a Julia `const` variable at the top level 
of your package will allow it to be precompiled and thus avoid
grammar re-analysis each time your package is loaded.

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://jamesrhester.github.io/Lerche.jl/dev/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: http://jamesrhester.github.io/Lerche.jl/stable/
# Introduction: Lerche.jl

Lerche.jl creates a parser for a language specified in an EBNF-like
syntax. The resulting parse trees can be transformed using
easy-to-specify methods, where Lerche.jl takes care of the parse tree
traversal. Lerche.jl is a direct translation of the Python-language
Lark parser generator (Lerche is "Lark" in German).

Much of the extensive [Lark
documentation](https://lark-parser.readthedocs.io/) is also relevant.

# Quick start

If you are already familiar with Lark see 
[Notes for Lark users](#Notes-for-Lark-users) below.

Lerche reads EBNF grammars as recognised by Lark
("Lark grammars") to produce a parser. This parser, when provided with
text conforming to the grammar, produces a parse tree. This tree can
be visited and transformed using "rules". A rule is a function named
after the production whose arguments it should be called on, and the
first argument of a rule is an object which is a subtype of
[`Visitor`](@ref) or [`Transformer`](@ref).

Given an EBNF grammar, it can be used to parse text into your data
structure as follows:

  1. Define one or more subtypes of [`Transformer`](@ref) or
     [`Visitor`](@ref), instances of which will be passed as the first
     argument to the appropriate rule.  The instance can also be used
     to hold information during transformation if you wish, in which
     case it must have a concrete type.
  1. Define `visit_tokens(t::MyNewType) = false` if you will not be
     processing token values. This is about 25% faster than leaving the
     default `true`.
  1. For every production in your grammar that you wish to process,
     write a rule with identical name to the production
  1. The rule should be prefixed with macro [`@rule`](@ref) if the second argument
     is an array containing all of the arguments to the grammar production
  1. The rule should be prefixed with macro [`@inline_rule`](@ref) if the second
     and following arguments refer to each argument in the grammar production
  1. For every token which you wish to process, define an identically-named method
     as for rules, but precede it with a [`@terminal`](@ref) macro instead of `@rule`.


If your grammar is in `String` variable `mygrammar`, your text to be parsed and transformed
is in `String` variable `mytext`, and your `Transformer` subtype is `MyTransformer`, the
following commands will produce a data structure from the text:

```julia
p = Lark(mygrammar,parser="lalr",lexer="contextual") #create parser
t = Lerche.parse(p,mytext)     #Create parse tree
x = Lerche.transform(MyTransformer(),t)  #transform parse tree
```

For a real-world example of usage, see [this
file](https://github.com/jamesrhester/DrelTools.jl/blob/master/src/jl_transformer.jl).

!!! tip

    Fully qualify the `parse` call (i.e. write `Lerche.parse`) to avoid ambiguity with `parse` from
    other packages, including `Base.parse`
    
## Error handling

When the supplied text does not match the grammar, `parse` raises exceptions that
are subtypes of `UnexpectedInput`:[`UnexpectedToken`](@ref)
and [`UnexpectedCharacters`](@ref). `Base.show` for these types produces an informative
message regarding the position of the error and expected tokens.

# Defining a Grammar

The full Lark grammar is described [here](grammar.md).

# Example

The following example (condensed from [the JSON example in
Lark](https://lark-parser.readthedocs.io/en/latest/json_tutorial.html) )
how a simple JSON parser is implemented.

```@setup json
using Lerche
```

First, the grammar:
```@repl json
json_grammar = raw"""
    ?start: value

    ?value: object
          | array
          | string
          | SIGNED_NUMBER      -> number
          | "true"             -> t
          | "false"            -> f
          | "null"             -> null

    array  : "[" [value ("," value)*] "]"
    object : "{" [pair ("," pair)*] "}"
    pair   : string ":" value

    string : ESCAPED_STRING

    %import common.ESCAPED_STRING
    %import common.SIGNED_NUMBER
    %import common.WS

    %ignore WS
""";
```

Note that terminals are always uppercase, and common definitions
can be imported from definitions in the standard library supplied with
Lerche.

A method whose name matches the rule name (or alias) and whose first
argument has our subtype of [`Transformer`](@ref) will be called
whenever that rule is matched.

These methods are prefixed by the [`@rule`](@ref) macro (if all of the
parse tree children are collected into a single array argument) or
[`@inline_rule`](@ref) macro (if each parse tree child is assigned a
separate argument). 

```@example json

struct TreeToJson <: Transformer end

@inline_rule string(t::TreeToJson, s) = replace(s[2:end-1],"\\\""=>"\"")

@rule  array(t::TreeToJson,a) = Array(a)
@rule  pair(t::TreeToJson,p) = Tuple(p)
@rule  object(t::TreeToJson,o) = Dict(o)
@inline_rule number(t::TreeToJson,n) = Base.parse(Float64,n)

@rule  null(t::TreeToJson,_) = nothing
@rule  t(t::TreeToJson,_) = true
@rule  f(t::TreeToJson,_) = false
```

The above rules define a `TreeToJson` subtype of `Transformer`, and rules whose
names match the rule or alias names in the grammar. For example,
whenever the `string` rule is matched, the enclosing double quotes
are dropped and any `\"` sequences replaced by a double quote.

Finally, we create our parser by calling the `Lark` constructor:

```@repl json
json_parser = Lark(json_grammar, parser="lalr", lexer="standard", transformer=TreeToJson());
```

Passing the `transformer` argument at parser construction time avoids
a separate call to the [`transform`](@ref) method after parsing.

Now, we can parse JSON by calling the [`Lerche.parse`](@ref) method with
`json_parser` as the first argument and the text to parse as the
second argument:

```@repl json
text = raw"{\"key\": [\"item0\", \"item1\", 3.14]}"
j = Lerche.parse(json_parser,text)
```

The above example is available in the Examples directory for
study.

## Other examples

The `tests` directory contains many more very simple examples
of correctly-constructed grammars.

# Notes for Lark users

When converting from Lark programs written in Python to Lerche
programs written in Julia, make the following changes:

  1. All Transformer and Visitor classes become types  
  2. All class method calls become Julia method calls with an instance
     of the type as the first argument (i.e. replacing `self`)
  3. Transformation or visitor rules should be preceded by the
     [`@rule`](@ref) macro. Inline rules use the [`@inline_rule`](@ref)
     macro and token processing methods use [`@terminal`](@ref). 
  4. Any grammars containing backslash-double quote sequences need to be edited (see below).
  5. Any grammars containing backslash-x to denote a byte value need to be edited (see below).

## Grammars

Lark grammars work unchanged in Lerche, with the caveats below.  Note
that this guarantee applies to the sequence of characters after
interpretation by the Julia or Python language parser.  In particular
note the following differences:

  1. Raw strings in Julia are written `raw"<string contents>"` instead
     of Python's `r"<string contents>"`
  2. The sequence `\"` inside a Python raw, quote-delimited string
     encodes a two-character sequence.  However, it corresponds to a
     single quote in Julia. To obtain the two-character sequence in
     Julia, write `\\"`. Such a backslash-quote sequence is required
     in Lark grammars to represent a double quote, just like in
     Python; so these two characters must remain in the string after
     Julia has pre-processed it.
  3. While unicode escapes are recognised (`\uxxxx`), the Python `\x`
     combination to insert a particular byte value in the string is not.
     Simply replace with the appropriate Unicode character.

# API Documentation

```@meta
CurrentModule = Lerche
```
## Parsing

A parser is created from a Lark grammar by calling the Lark constructor. Parsing
is initiated by calling `parse`.

```@docs
Lark(grammar::String;options...)
Lark(grammar::IOStream,source;options...)
Lerche.open(grammar_filename;rel_to=nothing,options...)
Lerche.parse(l::Lark,text;start=nothing,on_error=nothing)
UnexpectedCharacters
UnexpectedToken
```

## Working with the parse tree

`transform` transforms the parse tree according to rules defined by
the user using `@rule` and `@inline_rule` macros. Tokens will also be
processed using methods defined using `@terminal` if `visit_tokens` 
returns `true` for that transformer type.
Token processing will slow down parse tree processing by around 20%.

```@docs
@rule s
@inline_rule s
@terminal s
Transformer
Transformer_InPlace
Transformer_InPlaceRecursive
visit_tokens(t::Transformer)
transform(tr::Transformer,tree)
Visitor
Visitor_Recursive
Interpreter
visit(v::Visitor,tree)
```
# Grammar Reference

Lerche.jl reads grammars in a syntax developed for the Python Lark
system.  These are referred to as "Lark grammars" below. The following
information is adapted from the grammar reference provided with Lark.

## Definitions

A **grammar** is a list of rules and terminals, that together define a language.

Terminals define the alphabet of the language, while rules define its structure.

In Lerche, a terminal may be a string, a regular expression, or a
concatenation of these and other terminals.

Each rule is a list of terminals and rules, whose location and nesting
define the structure of the resulting parse-tree.

A **parsing algorithm** is an algorithm that takes a grammar
definition and a sequence of symbols (members of the alphabet), and
matches the entirety of the sequence by searching for a structure that
is allowed by the grammar.

### General Syntax and notes

Lark grammars are based on
[EBNF](https://en.wikipedia.org/wiki/Extended_Backus–Naur_form)
syntax, with several enhancements.

EBNF is basically a short-hand for common BNF patterns.

Optionals are expanded:

```ebnf
  a b? c    ->    (a c | a b c)
```

Repetition is extracted into a recursion:

```ebnf
  a: b*    ->    a: _b_tag
                 _b_tag: (_b_tag b)?
```

And so on.

Lark grammars are composed of a list of definitions and directives,
each on its own line. A definition is either a named rule, or a named
terminal, with the following syntax, respectively:

```c
  rule: <EBNF EXPRESSION>
      | etc.

  TERM: <EBNF EXPRESSION>   // Rules aren't allowed
```


**Comments** start with `//` and last to the end of the line (C++ style)

Lerche begins the parse with the rule 'start', unless specified otherwise in the options.

Names of rules are always in lowercase, while names of terminals are
always in uppercase. This distinction has practical effects, for the
shape of the generated parse-tree, and the automatic construction of
the lexer (aka tokenizer, or scanner).


## Terminals

Terminals are used to match text into symbols. They can be defined as
a combination of literals and other terminals.

**Syntax:**

```html
<NAME> [. <priority>] : <literals-and-or-terminals>
```

Terminal names must be uppercase.

Literals can be one of:

* `"string"`
* `/regular expression+/`
* `"case-insensitive string"i`
* `/re with flags/imulx`
* Literal range: `"a".."z"`, `"1".."9"`, etc.

Terminals also support grammar operators, such as `|`, `+`, `*` and `?`.

Terminals are a linear construct, and therefore may not contain
themselves (recursion isn't allowed).

### Templates

Templates are expanded when preprocessing the grammar.

Definition syntax:

```ebnf
  my_template{param1, param2, ...}: <EBNF EXPRESSION>
```

Use syntax:

```ebnf
some_rule: my_template{arg1, arg2, ...}
```

Example:
```ebnf
_separated{x, sep}: x (sep x)*  // Define a sequence of 'x sep x sep x ...'

num_list: "[" _separated{NUMBER, ","} "]"   // Will match "[1, 2, 3]" etc.
```

### Priority

Terminals can be assigned priority only when using a lexer (future
versions may support Earley's dynamic lexing).

Priority can be either positive or negative. If not specified for a terminal, it defaults to 1.

Highest priority terminals are always matched first.

### Regexp Flags

You can use flags on regexps and strings. For example:

```perl
SELECT: "select"i     //# Will ignore case, and match SELECT or Select, etc.
MULTILINE_TEXT: /.+/s
SIGNED_INTEGER: /
    [+-]?  # the sign
    (0|[1-9][0-9]*)  # the digits
 /x
```

Supported flags are one of: `imslux`. See Julia's regex documentation for more details on each one.

Regexps/strings of different flags may not be concatenated.

#### Notes for using a lexer

When using a lexer (standard or contextual), it is the
grammar-author's responsibility to make sure the literals don't
collide, or that if they do, they are matched in the desired
order. Literals are matched according to the following precedence:

1. Highest priority first (priority is specified as: TERM.number: ...)
2. Length of match (for regexps, an estimate of longest theoretical match is used)
3. Length of literal / pattern definition
4. Name

**Examples:**
```perl
IF: "if"
INTEGER : /[0-9]+/
INTEGER2 : ("0".."9")+          //# Same as INTEGER
DECIMAL.2: INTEGER? "." INTEGER  //# Will be matched before INTEGER
WHITESPACE: (" " | /\t/ )+
SQL_SELECT: "select"i
```

### Regular expressions & Ambiguity

Each terminal is eventually compiled to a regular expression. All the
operators and references inside it are mapped to their respective
expressions.

For example, in the following grammar, `A1` and `A2`, are equivalent:
```perl
A1: "a" | "b"
A2: /a|b/
```

This means that inside terminals, Lerche.jl cannot detect or resolve
ambiguity.

For example, for this grammar:
```perl
start           : (A | B)+
A               : "a" | "ab"
B               : "b"
```

We get this behavior:

```bash
>>> parse(p, "ab")
Tree(start, [Token(A, "a"), Token(B, "b")])
```

This is happening because the regex engine always returns the first matching option.

## Rules

**Syntax:**
```html
<name> : <items-to-match>  [-> <alias> ]
       | ...
```

Names of rules and aliases are always in lowercase.

Rule definitions can be extended to the next line by using the OR
operator (signified by a pipe: `|` ).

An alias is a name for the specific rule alternative. It affects tree
construction.


Each item is one of:

* `rule`
* `TERMINAL`
* `"string literal"` or `/regexp literal/`
* `(item item ..)` - Group items
* `[item item ..]` - Maybe. Same as `(item item ..)?`, but when `maybe_placeholders=True`, generates `nothing` if there is no match.
* `item?` - Zero or one instances of item ("maybe")
* `item*` - Zero or more instances of item
* `item+` - One or more instances of item
* `item ~ n` - Exactly *n* instances of item
* `item ~ n..m` - Between *n* to *m* instances of item (not recommended for wide ranges, due to performance issues)

**Examples:**
```perl
hello_world: "hello" "world"
mul: (mul "*")? number     //# Left-recursion is allowed and encouraged!
expr: expr operator expr
    | value               //# Multi-line, belongs to expr

four_words: word ~ 4
```

### Priority

Rules can be assigned priority only when using Earley (future versions
may support LALR as well).

Priority can be either positive or negative. In not specified for a
terminal, it's assumed to be 1 (i.e. the default).

## Directives

### %ignore

All occurrences of the terminal will be ignored, and won't be part of the parse.

Using the `%ignore` directive results in a cleaner grammar.

It's especially important for the LALR(1) algorithm, because adding
whitespace (or comments, or other extraneous elements) explicitly in
the grammar, harms its predictive abilities, which are based on a
lookahead of 1.

**Syntax:**
```html
%ignore <TERMINAL>
```
**Examples:**
```perl
%ignore " "

COMMENT: "#" /[^\n]/*
%ignore COMMENT
```
### %import

Allows one to import terminals and rules from lark grammars.

When importing rules, all their dependencies will be imported into a
namespace, to avoid collisions. It's not possible to override their
dependencies (e.g. like you would when inheriting a class).

**Syntax:**
```html
%import <module>.<TERMINAL>
%import <module>.<rule>
%import <module>.<TERMINAL> -> <NEWTERMINAL>
%import <module>.<rule> -> <newrule>
%import <module> (<TERM1>, <TERM2>, <rule1>, <rule2>)
```

If the module path is absolute, Lerche will attempt to load it from
the built-in directory (which currently contains `common.lark`).

If the module path is relative, such as `.path.to.file`, Lark will
attempt to load it from the current working directory. Grammars must
have the `.lark` extension.

The rule or terminal can be imported under another name with the `->` syntax.

**Example:**
```perl
%import common.NUMBER

%import .terminals_file (A, B, C)

%import .rules_file.rulea -> ruleb
```

Note that `%ignore` directives cannot be imported. Imported rules will
abide by the `%ignore` directives declared in the main grammar.

### %declare

Declare a terminal without defining it. Useful for plugins.

