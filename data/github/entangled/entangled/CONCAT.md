---
title: Entangled
author: Johan Hidding
---

[![DOI](https://zenodo.org/badge/160842312.svg)](https://zenodo.org/badge/latestdoi/160842312)
[![Entangled badge](https://img.shields.io/badge/entangled-Use%20the%20source!-%2300aeff)](https://entangled.github.io/)
[![Haskell CI](https://github.com/entangled/entangled/actions/workflows/haskell-ci.yml/badge.svg)](https://github.com/entangled/entangled/actions/workflows/haskell-ci.yml)

> **Literate programming** [/ˈlɪtəɹət ˈpɹəʊɡɹæmɪŋ/]{.phonetic} (computing) Literate programming is a programming paradigm introduced by Donald Knuth in which a program is given as an explanation of the program logic in a natural language, such as English, interspersed with snippets of macros and traditional source code, from which a compilable source code can be generated. [(Wikipedia)](https://en.wikipedia.org/wiki/Literate_programming)

In short: you write Markdown containing code fragments. These code fragments are combined into working code in a process called **tangling**.

`Entangled` makes writing literate programs easier by keeping code blocks in markdown up-to-date with generated source files. By monitoring the tangled source files, any change in the master document or source files is reflected in the other. In practice this means:

*    Write well documented code using Markdown.
*    Use any programming language you like (or are forced to use).
*    Keep debugging and using other IDE features without change.
*    Generate a report in PDF or HTML from the same source (see examples at [Entangled homepage](https://entangled.github.io/)).

# Status

`Entangled` is approaching 1.0 release! It has been tested Linux, Windows and MacOS. Still, it is highly recommended to use version control and *commit often*. If you encounter unexpected behaviour, please post an issue and describe the steps to reproduce.

Features:

- live bi-directional updates
- (reasonably) robust against wrongly edited source files
- configurable with [Dhall](https://dhall-lang.org/)
- hackable through SQLite
- create PDF or HTML pages from literate source
- line directives to point compilers to markdown source

# Building

`Entangled` is written in [Haskell](https://www.haskell.org/), and uses the `cabal` build system. You can build an executable by running

    # (requires cabal >= 3.x)
    cabal build

Install the executable in your `~/.local/bin`

    cabal install

Run unit tests

    cabal test

# Using

`Entangled` should be run from the command-line. The idea is that you run it from the root folder of the project that you're working on. This folder should contain a `entangled.dhall` file that contains the configuration. You can get an example config file by running

    entangled config

This config asumes you have the markdown files in a folder named `./lit`, and stores information in a SQLite3 database located at `./.entangled/db`. To run the daemon,

    entangled daemon [files ...]

where the `[files ...]` bits is sequence of additional files that you want monitored.

# Syntax (markdown side)

The markdown syntax `Entangled` uses is compatible with `Pandoc`'s. This relies on the use of *fenced code attributes*. To tangle a code block to a file:

~~~markdown
``` {.bash file=src/count.sh}
   ...
```
~~~

Composing a file using multiple code blocks is done through *noweb* syntax. You can reference a named code block in another code block by putting something like `<<named-code-block>>` on a single line. This reference may be indented. Such an indentation is then prefixed to each line in the final result.

A named code block should have an identifier given:

~~~markdown
``` {.python #named-code-block}
   ...
```
~~~

If a name appears multiple times in the source, the code blocks are concatenated during tangling. When weaving, the first code block with a certain name will appear as `<<name>>=`, while consecutive code blocks with the same name will appear as `<<name>>+=`.

Please see the [Hello World](https://entangled.github.io/examples/hello-world.html) and [other examples](https://entangled.github.io/examples)!

## Syntax (source side)

In the source code we know exactly where the code came from, so there would be no strict need for extra syntax there. However, once we start to edit the source file it may not be clear where the extra code needs to end up. To make our life a little easier, named code blocks that were tangled into the file are marked with a comment at begin and end.

```cpp
// ~|~ begin <<lit/story.md|main-body>>[0]
std::cout << "Hello, World!" << std::endl;
// ~|~ end
```

These comments should not be tampered with!

## Running `entangled`

Assuming you have created a Markdown file, say `program.md`, you can start `entangled` by running

```
entangled daemon ./program.md
```

in the shell. You may run `entangled --help` to get help on options, or check out [the user manual](https://entangled.github.io/manual.html).

## Running `entangled` with Docker

Entangled is available as a [Docker image](https://hub.docker.com/r/nlesc/entangled).

Assuming you have created a Markdown file, say `program.md`, you can start `entangled` by running

```shell
docker run --rm --user $(id -u):$(id -g) --volume $PWD:/data nlesc/entangled daemon ./program.md
```

This command starts a Docker container with the current working directory mounted as /data and running with your user/group id so files are written with the correct ownership.

## Distribution

If you've written a literate code using Entangled and would like to distribute it, one way is to include the tangled source code in the tar ball. You may also wish to use the pandoc filters included in [`entangled/filters`](https://github.com/entangled/filters).

# Development

## Credits

The following persons have made contributions to Entangled:

- Michał J. Gajda (gh:mgajda), first implemented the line-directive feature
- Danny Wilson (gh:vizanto), first implemented the project annotation

## Generating manpage

```
pandoc lit/a2-manpage.md -s -t man | /usr/bin/man -l -
```

## License

Entangled is distributed under the Apache v2 license.

# Packaging
Everything related to packaging should go in this directory. 

# Static Linux binary
This contains singularity definition files. The following explains how to build in sandbox mode. From the project root, to build an image:

    singularity build -f --sandbox /tmp/alpine-entangled packaging/alpine-build.def

To make the static binary distribution:

    singularity run -f --no-home --bind .:/mnt --writable /tmp/alpine-entangled

The resulting tarball has the following structure:

    .
    ├── bin
    │   └── entangled
    ├── lib
    │   └── entangled
    │       ├── data
    │       │   ├── config-schema.dhall
    │       │   ├── example-config.dhall
    │       │   ├── minimal-config.dhall
    │       │   └── schema.sql
    │       └── entangled
    └── share
        └── doc
            └── entangled
                ├── CITATION.cff
                ├── LICENSE
                └── README.md

# Fedora

- [Packaging instructions](https://docs.fedoraproject.org/en-US/quick-docs/creating-rpm-packages/)
# Configuration

The configuration is written in Dhall and has the following schema:

``` {.dhall file=data/config-schema.dhall}
let Comment : Type = < Line : Text | Block : { start : Text, end : Text } >
let Language : Type = { name : Text, identifiers : List Text, comment : Comment }
let LineDirective : Type = { name : Text, format: Text }

<<config-comment-styles>>
<<config-languages>>

let Annotate = < Naked | Standard | Project >

let Syntax : Type =
    { matchCodeStart       : Text
    , extractLanguage      : Text
    , extractFileName      : Text
    , extractReferenceName : Text
    , matchCodeEnd         : Text
    , extractProperty      : Text -> Text }

let defaultSyntax : Syntax =
    { matchCodeStart       = "^[ ]*```[ ]*{[^{}]*}"
    , matchCodeEnd         = "^[ ]*```"
    , extractLanguage      = "^[ ]*```[ ]*{\\.([^{} \t]+)[^{}]*}"
    , extractReferenceName = "^[ ]*```[ ]*{[^{}]*#([^{} \t]*)[^{}]*}"
    , extractFileName      = "^[ ]*```[ ]*{[^{}]*file=([^{} \t]*)[^{}]*}"
    , extractProperty      = \(name : Text) -> "^[ ]*```[ ]*{[^{}]*${name}=([^{} \t]*)[^{}]*}" }

let Config =
    { Type =
        { version   : Text
        , languages : List Language
        , watchList : List Text
        , database  : Optional Text
        , syntax    : Syntax
        , annotate  : Annotate
        , lineDirectives : List LineDirective
        , useLineDirectives : Bool }
    , default =
        { version   = "1.3.0"
        , languages = languages
        , watchList = [] : List Text
        , database  = None Text
        , syntax    = defaultSyntax
        , annotate  = Annotate.Standard
        , lineDirectives = lineDirectives
        , useLineDirectives = False }
    }

in { Comment   = Comment
   , Language  = Language
   , LineDirective = LineDirective
   , Config    = Config
   , Annotate  = Annotate
   , Syntax    = Syntax
   , comments  = comments
   , languages = languages
   , lineDirectives = lineDirectives
   , defaultSyntax = defaultSyntax
   }
```

## Languages

``` {.dhall #config-comment-styles}
let comments =
    { hash         = Comment.Line "#"
    , lispStyle    = Comment.Line ";"
    , cStyle       = Comment.Block { start = "/*", end = "*/" }
    , cppStyle     = Comment.Line "//"
    , haskellStyle = Comment.Line "--"
    , mlStyle      = Comment.Block { start = "(*", end = "*)" }
    , xmlStyle     = Comment.Block { start = "<!--", end = "-->" }
    , texStyle     = Comment.Line "%"
    }
```

``` {.dhall #config-languages}
let languages =
    [ { name = "Awk",        identifiers = ["awk"],           comment = comments.hash }
    , { name = "C",          identifiers = ["c"],             comment = comments.cStyle }
    , { name = "C++",        identifiers = ["cpp", "c++"],    comment = comments.cppStyle }
    , { name = "Clojure",    identifiers = ["clojure"],       comment = comments.lispStyle }
    , { name = "CSS",        identifiers = ["css"],           comment = comments.cStyle }
    , { name = "D",          identifiers = ["d"],             comment = comments.cppStyle }
    , { name = "Dhall",      identifiers = ["dhall"],         comment = comments.haskellStyle }
    , { name = "Elm",        identifiers = ["elm"],           comment = comments.haskellStyle }
    , { name = "Gnuplot",    identifiers = ["gnuplot"],       comment = comments.hash }
    , { name = "Haskell",    identifiers = ["haskell"],       comment = comments.haskellStyle }
    , { name = "HTML",       identifiers = ["html"],          comment = comments.xmlStyle }
    , { name = "Idris",      identifiers = ["idris"],         comment = comments.haskellStyle }
    , { name = "Julia",      identifiers = ["julia"],         comment = comments.hash }
    , { name = "JavaScript", identifiers = ["js", "javascript", "ecma"],
                                                              comment = comments.cStyle }
    , { name = "LaTeX",      identifiers = ["latex"],         comment = comments.texStyle }
    , { name = "Lua",        identifiers = ["lua"],           comment = comments.haskellStyle }
    , { name = "Make",       identifiers = ["make", "makefile"],
                                                              comment = comments.hash }
    , { name = "OCaml",      identifiers = ["ocaml"],         comment = comments.mlStyle }
    , { name = "OpenCL",     identifiers = ["opencl"],        comment = comments.cStyle }
    , { name = "PureScript", identifiers = ["purs", "purescript"],
                                                              comment = comments.haskellStyle }
    , { name = "Python",     identifiers = ["py", "python"],  comment = comments.hash }
    , { name = "R",          identifiers = ["r"],             comment = comments.hash }
    , { name = "Rust",       identifiers = ["rust"],          comment = comments.cppStyle }
    , { name = "Scheme",     identifiers = ["scheme", "r6rs" ,"racket", "r7rs"],
                                                              comment = comments.lispStyle }
    , { name = "SQLite",     identifiers = ["sqlite"],        comment = comments.haskellStyle }
    , { name = "TOML",       identifiers = ["toml"],          comment = comments.hash }
    , { name = "TypeScript", identifiers = ["ts", "typescript"],
                                                              comment = comments.cppStyle }
    , { name = "YAML",       identifiers = ["yaml"],          comment = comments.hash }
    , { name = "<unknown>",  identifiers = [] : List Text,    comment = comments.hash }
    ]

let lineDirectives =
    [ { name = "C",          format = "#line {linenumber} \"{filename}\"" }
    , { name = "C++",        format = "#line {linenumber} \"{filename}\"" }
    , { name = "Haskell",    format = "{{-# LINE {linenumber} \"{filename}\" #-}}" }
    ]
```

## Older versions
We still would like to read older version of the schema.

``` {.haskell file=src/Config/Version_1_0_0.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Config.Version_1_0_0 where

import RIO

import Config.Record
import Format
import Dhall (auto, Decoder, record, field, setFromDistinctList)

<<config-1-0-0-record>>
<<config-1-0-0-decoder>>
```

``` {.haskell file=src/Config/Version_1_2_0.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Config.Version_1_2_0 where

import RIO

import Config.Record hiding (ConfigSyntax, configSyntaxDecoder)
import qualified Config.Version_1_0_0 as Version_1_0_0
import Format

import Dhall (auto, Decoder, record, field, setFromDistinctList )

<<config-1-2-0-record>>
<<config-1-2-0-decoder>>
```

``` {.haskell file=src/Config/Version_1_3_0.hs}
{-# LANGUAGE NoImplicitPrelude,UndecidableInstances #-}
module Config.Version_1_3_0 where

import RIO
import qualified RIO.Text as T

import Config.Record
import qualified Config.Version_1_2_0 as Version_1_2_0
import Format

import Paths_entangled
import Dhall (input, auto, Decoder, record, field, setFromDistinctList )

<<config-1-3-0-record>>
<<config-1-3-0-decoder>>
```

``` {.haskell file=src/Config/Record.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Config.Record where

import RIO
import qualified RIO.Map as M

<<config-imports>>
<<config-dhall-schema>>
```

## Reading config

``` {.haskell #config-imports}
import Dhall (FromDhall, ToDhall, auto, Decoder, record, list
             , field, constructor, unit, union)

import qualified Format
```

We need to match the Dhall schema with types in Haskell

### Language

``` {.haskell #config-dhall-schema}
data ConfigComment
    = Line  Text
    | Block { start :: Text, end :: Text }
    deriving (Generic, Show)

instance FromDhall ConfigComment
instance ToDhall ConfigComment

data ConfigLanguage = ConfigLanguage
    { languageName :: Text
    , languageIdentifiers :: [Text]
    , languageComment :: ConfigComment
    } deriving (Show)

configLanguage :: Decoder ConfigLanguage
configLanguage = record
    ( ConfigLanguage <$> field "name"        auto
                     <*> field "identifiers" auto
                     <*> field "comment"     auto
    )

instance Eq ConfigLanguage where
    a == b = languageName a == languageName b

instance Ord ConfigLanguage where
    compare a b = compare (languageName a) (languageName b)
```

### Line directives

``` {.haskell #config-dhall-schema}
decodeFormatSpec :: Decoder Format.Spec
decodeFormatSpec = fromMaybe [Format.Plain "illegal format spec"] . Format.spec <$> auto

lineDirectivesDecoder :: Decoder (Map Text Format.Spec)
lineDirectivesDecoder = M.fromList <$> list entry
    where entry = record ( pair <$> field "name" auto
                                <*> field "format" decodeFormatSpec )
          pair a b = (a, b)
```

### Annotation method

``` {.haskell #config-dhall-schema}
data AnnotateMethod = AnnotateNaked
                    | AnnotateStandard
                    | AnnotateProject
                    deriving (Show, Eq)

annotateDecoder :: Decoder AnnotateMethod
annotateDecoder = union
        (  ( AnnotateNaked          <$ constructor "Naked" unit )
        <> ( AnnotateStandard       <$ constructor "Standard" unit )
        <> ( AnnotateProject        <$ constructor "Project" unit ) )
```

### Syntax

``` {.haskell #config-dhall-schema}
data ConfigSyntax = ConfigSyntax
    { matchCodeStart       :: Text
    , matchCodeEnd         :: Text
    , extractLanguage      :: Text
    , extractReferenceName :: Text
    , extractFileName      :: Text
    , extractProperty      :: Text -> Text
    }

configSyntaxDecoder :: Decoder ConfigSyntax
configSyntaxDecoder = record
    ( ConfigSyntax <$> field "matchCodeStart" auto
                   <*> field "matchCodeEnd" auto
                   <*> field "extractLanguage" auto
                   <*> field "extractReferenceName" auto
                   <*> field "extractFileName" auto
                   <*> field "extractProperty" auto )
```

### Decoders

#### Version 1.0.0

``` {.haskell #config-1-0-0-record}
data Config = Config
    { configVersion   :: Text
    , configLanguages :: Set ConfigLanguage
    , configWatchList :: [Text]
    , configDatabase  :: Maybe Text
    , configAnnotate  :: AnnotateMethod
    , configLineDirectives :: Map Text Format.Spec
    , configUseLineDirectives :: Bool
    } deriving (Show)
```

``` {.haskell #config-1-0-0-decoder}
configDecoder :: Decoder Config
configDecoder = record
    ( Config <$> field "version" auto
             <*> field "languages" (setFromDistinctList configLanguage)
             <*> field "watchList" auto
             <*> field "database" auto
             <*> field "annotate" annotateDecoder
             <*> field "lineDirectives" lineDirectivesDecoder
             <*> field "useLineDirectives" auto
    )
```


#### Version 1.2.0

``` {.haskell #config-1-2-0-record}
data ConfigSyntax = ConfigSyntax
    { matchCodeStart       :: Text
    , matchCodeEnd         :: Text
    , extractLanguage      :: Text
    , extractReferenceName :: Text
    , extractFileName      :: Text
    }

configSyntaxDecoder :: Decoder ConfigSyntax
configSyntaxDecoder = record
    ( ConfigSyntax <$> field "matchCodeStart" auto
                   <*> field "matchCodeEnd" auto
                   <*> field "extractLanguage" auto
                   <*> field "extractReferenceName" auto
                   <*> field "extractFileName" auto )

data Config = Config
    { configVersion   :: Text
    , configLanguages :: Set ConfigLanguage
    , configWatchList :: [Text]
    , configDatabase  :: Maybe Text
    , configSyntax    :: ConfigSyntax
    , configAnnotate  :: AnnotateMethod
    , configLineDirectives :: Map Text Format.Spec
    , configUseLineDirectives :: Bool
    }

defaultSyntax :: ConfigSyntax
defaultSyntax = ConfigSyntax
    { matchCodeStart       = "^[ ]*```[ ]*{[^{}]*}"
    , matchCodeEnd         = "^[ ]*```"
    , extractLanguage      = "^[ ]*```[ ]*{\\.([^{} \t]+)[^{}]*}"
    , extractReferenceName = "^[ ]*```[ ]*{[^{}]*#([^{} \t]*)[^{}]*}"
    , extractFileName      = "^[ ]*```[ ]*{[^{}]*file=([^{} \t]*)[^{}]*}" }

class ToVersion_1_2_0 a where
    update :: a -> IO Config

instance ToVersion_1_2_0 Config where
    update = return

instance ToVersion_1_2_0 Version_1_0_0.Config where
    update old = do
        return Config
            { configVersion           = Version_1_0_0.configVersion           old
            , configLanguages         = Version_1_0_0.configLanguages         old
            , configWatchList         = Version_1_0_0.configWatchList         old
            , configDatabase          = Version_1_0_0.configDatabase          old
            , configSyntax            = defaultSyntax
            , configAnnotate          = Version_1_0_0.configAnnotate          old
            , configLineDirectives    = Version_1_0_0.configLineDirectives    old
            , configUseLineDirectives = Version_1_0_0.configUseLineDirectives old
            }
```

``` {.haskell #config-1-2-0-decoder}
configDecoder :: Decoder Config
configDecoder = record
    ( Config <$> field "version" auto
             <*> field "languages" (setFromDistinctList configLanguage)
             <*> field "watchList" auto
             <*> field "database" auto
             <*> field "syntax" configSyntaxDecoder
             <*> field "annotate" annotateDecoder
             <*> field "lineDirectives" lineDirectivesDecoder
             <*> field "useLineDirectives" auto
    )
```

#### Version 1.3.0

``` {.haskell #config-1-3-0-record}
data Config = Config
    { configVersion   :: Text
    , configLanguages :: Set ConfigLanguage
    , configWatchList :: [Text]
    , configDatabase  :: Maybe Text
    , configSyntax    :: ConfigSyntax
    , configAnnotate  :: AnnotateMethod
    , configLineDirectives :: Map Text Format.Spec
    , configUseLineDirectives :: Bool
    }

defaultSyntax :: IO ConfigSyntax
defaultSyntax = do
    path <- getDataFileName "data/config-schema.dhall"
    input configSyntaxDecoder $ "(" <> T.pack path <> ").defaultSyntax"

class ToVersion_1_3_0 a where
    update :: a -> IO Config

instance ToVersion_1_3_0 Config where
    update = return

syntax120to130 :: Version_1_2_0.ConfigSyntax -> ConfigSyntax
syntax120to130 Version_1_2_0.ConfigSyntax {..}
    = ConfigSyntax
        { extractProperty      = \name -> "^[ ]*```[ ]*{[^{}]*" <> name <> "=([^{} \t]*)[^{}]*}"
        , .. }

instance Version_1_2_0.ToVersion_1_2_0 a => ToVersion_1_3_0 a where
    update old' = do
        old <- Version_1_2_0.update old'
        return Config
            { configVersion           = Version_1_2_0.configVersion           old
            , configLanguages         = Version_1_2_0.configLanguages         old
            , configWatchList         = Version_1_2_0.configWatchList         old
            , configDatabase          = Version_1_2_0.configDatabase          old
            , configSyntax            = syntax120to130 $ Version_1_2_0.configSyntax old
            , configAnnotate          = Version_1_2_0.configAnnotate          old
            , configLineDirectives    = Version_1_2_0.configLineDirectives    old
            , configUseLineDirectives = Version_1_2_0.configUseLineDirectives old
            }
```

``` {.haskell #config-1-3-0-decoder}
configDecoder :: Decoder Config
configDecoder = record
    ( Config <$> field "version" auto
             <*> field "languages" (setFromDistinctList configLanguage)
             <*> field "watchList" auto
             <*> field "database" auto
             <*> field "syntax" configSyntaxDecoder
             <*> field "annotate" annotateDecoder
             <*> field "lineDirectives" lineDirectivesDecoder
             <*> field "useLineDirectives" auto
    )
```

``` {.haskell file=src/Config.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Config ( module Config
              , module Version_1_3_0
              , module Config.Record ) where

import RIO hiding (void)
import RIO.Directory
import RIO.FilePath
import RIO.List (scanl1, find)
import qualified RIO.Text as T

import Dhall (input, auto)
import System.FilePath.Glob

import qualified Config.Version_1_0_0 as Version_1_0_0
import qualified Config.Version_1_2_0 as Version_1_2_0
import qualified Config.Version_1_3_0 as Version_1_3_0
import Config.Version_1_3_0 (update, Config(..))
import Config.Record

import Errors
import Select

<<config-input>>
<<config-reader>>

class HasConfig env where
    config :: Lens' env Config

getDatabasePath :: (MonadIO m, MonadThrow m) => Config -> m FilePath
getDatabasePath cfg = do
    dbPath <- case configDatabase cfg of
        Nothing -> return ":memory:"
        Just db -> return $ T.unpack db
    liftIO $ createDirectoryIfMissing True (takeDirectory dbPath)
    return dbPath

getInputFiles :: (MonadIO m) => Config -> m [FilePath]
getInputFiles cfg = liftIO $ mconcat <$> mapM (glob . T.unpack) (configWatchList cfg)
```

> ~~Configuration can be stored in `${XDG_CONFIG_HOME}/entangled/config.dhall`. Also the local directory or its parents may contain a `.entangled.dhall` file. These override settings in the global configuration.~~
> This is a bad idea. Configurations should work the same for people cloning a repository. Without this idea, also stacking of configurations is not needed. Dhall should handle all that.

There are currently no customisation options for entangled, but I keep my options open: namespaces for references, enabling future features like git support, you name it.

### Reading config files

``` {.haskell #config-input}
findFileAscending :: String -> IO (Maybe FilePath)
findFileAscending filename = do
    path <- dropTrailingPathSeparator <$> getCurrentDirectory
    let parents = reverse $ scanl1 (</>) $ splitDirectories path
    findFile parents filename

getVersion :: FilePath -> IO Text
getVersion path =
    input auto $ "(" <> T.pack path <> ").entangled.version"

readLocalConfig :: IO Config
readLocalConfig = do
    cfg_path <- findFileAscending "entangled.dhall"
            >>= maybe (throwM $ SystemError "no config found") return
    version <- getVersion cfg_path
    decoder <- select (throwM $ SystemError $ "unrecognized version string '" <> version <> "'")
        [ ( version == "1.0.0", return $ update <=< input Version_1_0_0.configDecoder )
        , ( version == "1.2.0", return $ update <=< input Version_1_2_0.configDecoder )
        , ( version == "1.3.0", return $ input Version_1_3_0.configDecoder ) ]
    decoder $ "(" <> T.pack cfg_path <> ").entangled"
```

## Processing

``` {.haskell #config-reader}
lookupLanguage :: Config -> Text -> Maybe ConfigLanguage
lookupLanguage cfg x
    = find (elem x . languageIdentifiers)
    $ configLanguages cfg

languageFromName :: Config -> Text -> Maybe ConfigLanguage
languageFromName cfg x
    = find ((== x) . languageName)
    $ configLanguages cfg
```
# Main program
The main program runs the daemon, but also provides a number of commands to inspect and manipulate the database.

## Encoding
On linux consoles we use unicode bullet points (`•`). On Windows, those will just be asterisks (`*`). To facilitate this, we have to enable UTF-8 encoding.

``` {.haskell #main-imports}
import GHC.IO.Encoding
```

``` {.haskell #main-set-encoding}
setLocaleEncoding utf8
```

## Options
Options are parsed using `optparse-applicative`.

``` {.haskell #main-imports}
import Options.Applicative
```

All true options are left to the sub-commands. We're leaving `<<sub-commands>>` to be expanded.

``` {.haskell #main-options}
data Args = Args
    { versionFlag :: Bool
    , verboseFlag :: Bool
    , machineFlag :: Bool
    , checkFlag   :: Bool
    , preinsertFlag :: Bool
    , subCommand :: SubCommand }

data SubCommand
    = NoCommand
    <<sub-commands>>
    deriving (Show, Eq)
```

The same goes for the sub-command parsers, which are collected in `<<sub-parsers>>`.

``` {.haskell #main-options}
parseNoCommand :: Parser SubCommand
parseNoCommand = pure NoCommand

parseArgs :: Parser Args   {- HLINT ignore parseArgs -}
parseArgs = Args
    <$> switch (long "version" <> short 'v' <> help "Show version information.")
    <*> switch (long "verbose" <> short 'V' <> help "Be very verbose.")
    <*> switch (long "machine" <> short 'm' <> help "Machine readable output.")
    <*> switch (long "check"   <> short 'c' <> help "Don't do anything, returns 1 if changes would be made to file system.")
    <*> switch (long "preinsert" <> short 'p' <> help "Tangle everything as a first action, default when db is in-memory.")
    <*> ( subparser ( mempty
          <<sub-parsers>>
        ) <|> parseNoCommand )
```

And the runners.

``` {.haskell #main-imports}
import Config
```

``` {.haskell #main-run}
data Env = Env
    { connection' :: Connection
    , config'     :: Config
    , logFunc'    :: LogFunc }

instance HasConnection Env where
    connection = lens connection' (\ x y -> x { connection' = y })

instance HasConfig Env where
    config = lens config' (\x y -> x { config' = y })

instance HasLogFunc Env where
    logFuncL = lens logFunc' (\x y -> x { logFunc' = y })

run :: Args -> IO ()
run (Args True _ _ _ _ _)                           = putStrLn $ showVersion version
run (Args _ _ _ _ _ (CommandConfig ConfigArgs{..})) = printExampleConfig' minimalConfig
run Args{..}                                        = runWithEnv verboseFlag machineFlag checkFlag preinsertFlag (runSubCommand subCommand)

runWithEnv :: Bool -> Bool -> Bool -> Bool -> Entangled Env a -> IO a
runWithEnv verbose machineReadable dryRun preinsertFlag x = do
    cfg <- readLocalConfig
    dbPath <- getDatabasePath cfg
    logOptions <- setLogVerboseFormat True . setLogUseColor True
               <$> logOptionsHandle stderr verbose
    let preinsertFlag' = preinsertFlag || dbPath == ":memory:"
        x' = (if preinsertFlag' then preinsert else pure ()) >> x
    if dryRun
    then do
        todo <- withLogFunc logOptions (\logFunc
                -> withConnection dbPath (\conn
                    -> runRIO (Env conn cfg logFunc) (testEntangled x')))
        if todo then exitFailure else exitSuccess
    else withLogFunc logOptions (\logFunc
        -> withConnection dbPath (\conn
            -> runRIO (Env conn cfg logFunc) (runEntangled machineReadable Nothing x')))

preinsert :: (HasConfig env, HasLogFunc env, HasConnection env)
          => Entangled env ()
preinsert = do
    db createTables
    cfg <- view config
    abs_paths <- sort <$> getInputFiles cfg
    when (null abs_paths) $ throwM $ SystemError "No input files."
    rel_paths <- mapM makeRelativeToCurrentDirectory abs_paths
    insertSources rel_paths

runSubCommand :: (HasConfig env, HasLogFunc env, HasConnection env)
              => SubCommand -> Entangled env ()
runSubCommand sc = do
    db createTables
    case sc of
        NoCommand -> return ()
        <<sub-runners>>
```

This way we can add sub-commands independently in the following sections.

### Starting the daemon

``` {.haskell #main-imports}
import Daemon (runSession)
```

``` {.haskell #sub-commands}
| CommandDaemon DaemonArgs
```

``` {.haskell #sub-parsers}
<>  command "daemon" (info parseDaemonArgs ( progDesc "Run the entangled daemon." ))
```

``` {.haskell #main-options}
newtype DaemonArgs = DaemonArgs
    { inputFiles  :: [String]
    } deriving (Show, Eq)

parseDaemonArgs :: Parser SubCommand
parseDaemonArgs = CommandDaemon . DaemonArgs
    <$> many (argument str (metavar "FILES..."))
    <**> helper
```

``` {.haskell #sub-runners}
CommandDaemon DaemonArgs {..} -> runSession inputFiles
```

### Printing the config

``` {.haskell #main-options}
newtype ConfigArgs = ConfigArgs
    { minimalConfig :: Bool
    } deriving (Show, Eq)

parseConfigArgs :: Parser SubCommand
parseConfigArgs = CommandConfig . ConfigArgs
    <$> switch (long "minimal" <> short 'm' <> help "Print minimal config.")
    <**> helper

printExampleConfig' :: Bool -> IO ()
printExampleConfig' minimal = do
    let path = if minimal then "data/minimal-config.dhall"
               else "data/example-config.dhall"
    T.IO.putStr =<< T.IO.readFile =<< getDataFileName path
```

``` {.haskell #sub-commands}
| CommandConfig ConfigArgs
```

``` {.haskell #sub-parsers}
<> command "config" (info parseConfigArgs
                          (progDesc "Print an example configuration."))
```

``` {.haskell #sub-runners}
CommandConfig _ -> printExampleConfig
```

### Inserting files to the database

``` {.haskell #sub-commands}
| CommandInsert InsertArgs
```

``` {.haskell #sub-parsers}
<> command "insert" (info parseInsertArgs ( progDesc "Insert markdown files into database." ))
```

``` {.haskell #main-options}
data FileType = SourceFile | TargetFile deriving (Show, Eq)

data InsertArgs = InsertArgs
    { insertType :: FileType
    , insertFiles :: [FilePath] } deriving (Show, Eq)

parseFileType :: Parser FileType
parseFileType = flag' SourceFile (long "source" <> short 's' <> help "insert markdown source file")
            <|> flag' TargetFile (long "target" <> short 't' <> help "insert target code file")

parseInsertArgs :: Parser SubCommand
parseInsertArgs = CommandInsert <$> (InsertArgs
    <$> parseFileType
    <*> many (argument str (metavar "FILES..."))
    <**> helper)
```

``` {.haskell #sub-runners}
CommandInsert (InsertArgs SourceFile fs) -> insertSources fs
CommandInsert (InsertArgs TargetFile fs) -> insertTargets fs
```

### Tangling a single reference

``` {.haskell #sub-commands}
| CommandTangle TangleArgs
```

``` {.haskell #sub-parsers}
<> command "tangle" (info (CommandTangle <$> parseTangleArgs) ( progDesc "Retrieve tangled code." ))
```

``` {.haskell #main-options}
data TangleArgs = TangleArgs
    { tangleQuery :: TangleQuery
    , tangleDecorate :: Bool
    } deriving (Show, Eq)

parseTangleArgs :: Parser TangleArgs
parseTangleArgs = TangleArgs
    <$> (   (TangleFile <$> strOption ( long "file" <> short 'f'
                                      <> metavar "TARGET" <> help "file target" ))
        <|> (TangleRef  <$> strOption ( long "ref"  <> short 'r'
                                      <> metavar "TARGET" <> help "reference target" ))
        <|> flag' TangleAll (long "all" <> short 'a' <> help "tangle all and write to disk" ))
    <*> switch (long "decorate" <> short 'd' <> help "Decorate with stitching comments.")
    <**> helper
```

``` {.haskell #sub-runners}
CommandTangle TangleArgs {..} -> do
    cfg <- view config
    tangle tangleQuery (if tangleDecorate
                        then selectAnnotator cfg
                        else selectAnnotator (cfg {configAnnotate = AnnotateNaked}))
```

### Stitching a markdown source

``` {.haskell #sub-commands}
| CommandStitch StitchArgs
```

``` {.haskell #sub-parsers}
<> command "stitch" (info (CommandStitch <$> parseStitchArgs) ( progDesc "Retrieve stitched markdown." ))
```

``` {.haskell #main-options}
newtype StitchArgs = StitchArgs
    { stitchTarget :: FilePath
    } deriving (Show, Eq)

parseStitchArgs :: Parser StitchArgs
parseStitchArgs = StitchArgs
    <$> argument str ( metavar "TARGET" )
    <**> helper
```

``` {.haskell #sub-runners}
CommandStitch StitchArgs {..} -> stitch (StitchFile stitchTarget)
```

### Listing all target files

``` {.haskell #sub-commands}
| CommandList
```

``` {.haskell #sub-parsers}
<> command "list" (info (pure CommandList <**> helper) ( progDesc "List generated code files." ))
```

``` {.haskell #sub-runners}
CommandList -> listTargets
```

### Linter

``` {.haskell #main-options}
newtype LintArgs = LintArgs
    { lintFlags :: [Text]
    } deriving (Show, Eq)

parseLintArgs :: Parser LintArgs
parseLintArgs = LintArgs
    <$> many (argument str (metavar "LINTERS"))
    <**> helper
```

``` {.haskell #sub-commands}
| CommandLint LintArgs
```

``` {.haskell #sub-parsers}
<> command "lint" (info (CommandLint <$> parseLintArgs) ( progDesc ("Lint input on potential problems. Available linters: " <> RIO.Text.unpack (RIO.Text.unwords allLinters))))
```

``` {.haskell #sub-runners}
CommandLint LintArgs {..} -> lint lintFlags
```

### Cleaning orphan targets
This action deletes orphan targets from both the database and the file system.

``` {.haskell #sub-commands}
| CommandClearOrphans
```

``` {.haskell #sub-parsers}
<> command "clear-orphans" (info (pure CommandClearOrphans <**> helper) ( progDesc "Deletes orphan targets." ))
```

``` {.haskell #sub-runners}
CommandClearOrphans -> clearOrphans
```

## Main

``` {.haskell file=app/Main.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Main where

import RIO
import RIO.Text (unwords, unpack)
import RIO.Directory (makeRelativeToCurrentDirectory)
import RIO.List (sort)

import Prelude (putStrLn)
import qualified Data.Text.IO as T.IO
import Paths_entangled
import Data.Version (showVersion)

<<main-imports>>

import Tangle (selectAnnotator)
import Entangled
import Errors (EntangledError(..))
import Linters

<<main-options>>

main :: IO ()
main = do
    <<main-set-encoding>>
    run =<< execParser args
    where args = info (parseArgs <**> helper)
            (  fullDesc
            <> progDesc "Automatically tangles and untangles 'FILES...'."
            <> header   "Entangled -- daemonised literate programming"
            )

<<main-run>>
```

## Generics

### Create the empty database

``` {.haskell #main-imports}
import Database (HasConnection, connection, createTables, db)
-- import Comment (annotateNaked)
import Database.SQLite.Simple
```

```
dbPath <- getDatabasePath cfg
withSQL dbPath $ do 
```

## Wiring

``` {.haskell file=src/Entangled.hs}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE MultiParamTypeClasses #-}
module Entangled where

import RIO
import RIO.Writer (MonadWriter, WriterT, runWriterT, tell)
import qualified RIO.Text as T

import qualified Data.Map.Lazy as LM
import Control.Monad.Except ( MonadError(..) )

import FileIO
import Transaction

import Console (Doc, timeStamp)
import Paths_entangled
import Config (config, HasConfig, languageFromName)
import Database ( db, HasConnection, queryTargetRef, queryReferenceMap
                , listTargetFiles, insertDocument, stitchDocument, listSourceFiles
                , deduplicateRefs, updateTarget, listOrphanTargets, clearOrphanTargets
                , queryCodeAttr )
import Errors (EntangledError (..))

import Comment (headerComment)
import Document (ReferenceName(..))
import Tangle (ExpandedCode, Annotator, expandedCode, parseMarkdown')
import Stitch (untangle)

type FileTransaction env = Transaction (FileIO env)

newtype Entangled env a = Entangled { unEntangled :: WriterT (FileTransaction env) (RIO env) a }
    deriving ( Applicative, Functor, Monad, MonadIO, MonadThrow
             , MonadReader env, MonadWriter (FileTransaction env) )

instance MonadError EntangledError (Entangled env) where
    throwError = throwM
    catchError x _ = x

testEntangled :: (MonadIO m, MonadReader env m, HasLogFunc env)
              => Entangled env a -> m Bool
testEntangled (Entangled x) = do
    e <- ask
    (_, w) <- runRIO e (runWriterT x)
    runFileIO' $ testTransaction w

runEntangled :: (MonadIO m, MonadReader env m, HasLogFunc env)
             => Bool -> Maybe Doc -> Entangled env a -> m a
runEntangled True  _ = runEntangledMachine
runEntangled False h = runEntangledHuman h

runEntangledMachine :: (MonadIO m, MonadReader env m, HasLogFunc env)
             => Entangled env a -> m a
runEntangledMachine (Entangled x) = do
    e <- ask
    (r, w) <- runRIO e (runWriterT x)
    runFileIO' $ runTransactionMachine w
    return r

runEntangledHuman :: (MonadIO m, MonadReader env m, HasLogFunc env)
             => Maybe Doc -> Entangled env a -> m a
runEntangledHuman h (Entangled x) = do
    e <- ask
    (r, w) <- runRIO e (runWriterT x)
    ts <- timeStamp
    runFileIO' $ runTransaction (h >>= (\h' -> Just $ ts <> " " <> h')) w
    return r

instance (HasLogFunc env) => MonadFileIO (Entangled env) where
    readFile = readFile'
    dump = dump'

    writeFile path text = do
        old_content' <- liftRIO $ try $ runFileIO' $ readFile path
        case (old_content' :: Either IOException Text) of
            Right old_content | old_content == text -> return ()
                              | otherwise           -> actionw
            Left  _                                 -> actionc
        where actionw   = tell $ plan (WriteFile path) (writeFile path text)
              actionc   = tell $ plan (CreateFile path) (writeFile path text)

    deleteFile path     = tell $ plan (DeleteFile path) (deleteFile path)

data TangleQuery = TangleFile FilePath | TangleRef Text | TangleAll deriving (Show, Eq)

tangleRef :: (HasLogFunc env, HasConfig env)
    => ExpandedCode (Entangled env) -> ReferenceName -> Entangled env Text
tangleRef codes name =
    case codes LM.!? name of
        Nothing        -> throwM $ TangleError $ "Reference `" <> tshow name <> "` not found."
        Just t         -> t

toInt :: Text -> Maybe Int
toInt = readMaybe . T.unpack

takeLines :: Text -> Int -> [Text]
takeLines txt n = take n $ drop 1 $ T.lines txt

dropLines :: Text -> Int -> [Text]
dropLines txt n = take 1 lines_ <> drop (n+1) lines_
    where lines_ = T.lines txt

tangleFile :: (HasConnection env, HasLogFunc env, HasConfig env)
           => ExpandedCode (Entangled env) -> FilePath -> Entangled env Text
tangleFile codes path = do
    cfg <- view config
    db (queryTargetRef path) >>= \case
        Nothing              -> throwM $ TangleError $ "Target `" <> T.pack path <> "` not found."
        Just (ref, langName) -> do
            content <- tangleRef codes ref
            headerLen <- db (queryCodeAttr ref "header")
            case languageFromName cfg langName of
                Nothing -> throwM $ TangleError $ "Language unknown " <> langName
                Just lang -> return $ maybe (T.unlines [headerComment lang path, content])
                                            (\n -> T.unlines $ takeLines content n <> [headerComment lang path] <> dropLines content n)
                                            (toInt =<< headerLen)

tangle :: (HasConnection env, HasLogFunc env, HasConfig env)
       => TangleQuery -> Annotator (Entangled env) -> Entangled env ()
tangle query annotate = do
    cfg <- view config
    refs <- db (queryReferenceMap cfg)
    let codes = expandedCode annotate refs
    case query of
        TangleRef ref   -> dump =<< tangleRef codes (ReferenceName ref)
        TangleFile path -> dump =<< tangleFile codes path
        TangleAll       -> mapM_ (\f -> writeFile f =<< tangleFile codes f) =<< db listTargetFiles

data StitchQuery = StitchFile FilePath | StitchAll

stitchFile :: (HasConnection env, HasLogFunc env, HasConfig env)
       => FilePath -> Entangled env Text
stitchFile path = db (stitchDocument path)

stitch :: (HasConnection env, HasLogFunc env, HasConfig env)
       => StitchQuery -> Entangled env ()
stitch (StitchFile path) = dump =<< stitchFile path
stitch StitchAll = mapM_ (\f -> writeFile f =<< stitchFile f) =<< db listSourceFiles

listTargets :: (HasConnection env, HasLogFunc env, HasConfig env)
            => Entangled env ()
listTargets = dump . T.unlines . map T.pack =<< db listTargetFiles

insertSources :: (HasConnection env, HasLogFunc env, HasConfig env)
              => [FilePath] -> Entangled env ()
insertSources files = do
    logDebug $ display $ "inserting files: " <> tshow files
    mapM_ readDoc files
    where readDoc f = do
            document <- parseMarkdown' f =<< readFile f
            db (insertDocument f document)

insertTargets :: (HasConnection env, HasLogFunc env, HasConfig env)
              => [FilePath] -> Entangled env ()
insertTargets files = do
    logDebug $ display $ "inserting files: " <> tshow files
    mapM_ readTgt files
    where readTgt f = do
            refs <- untangle f =<< readFile f
            db (updateTarget =<< deduplicateRefs refs)

clearOrphans :: (HasConnection env, HasLogFunc env, HasConfig env)
             => Entangled env ()
clearOrphans = do
    files <- db $ do
        r <- listOrphanTargets
        clearOrphanTargets
        return r
    mapM_ deleteFile files

printExampleConfig :: (HasLogFunc env)
                   => Entangled env ()
printExampleConfig = dump =<< readFile =<< liftIO (getDataFileName "data/example-config.dhall")
```
# Database
We use an SQLite database to manage document content. Using SQL requires a remapping of the available data.

``` {.haskell file=src/Database.hs}
module Database where

import RIO
import RIO.List (initMaybe, sortOn, nub)
import qualified RIO.Text as T
import qualified RIO.Map as M

<<database-imports>>
<<database-types>>
<<database-create>>
<<database-insertion>>
<<database-update>>
<<database-queries>>
<<database-deduplicate>>
```

## Types

We wrap all SQL interaction in a `SQL` monad, which stores the `Connection` object and a logger.

``` {.haskell #database-imports}
import Paths_entangled

import Database.SQLite.Simple
import Database.SQLite.Simple.FromRow ()

import qualified Data.Text.IO as T.IO
```

``` {.haskell #database-types}
class HasConnection env where
    connection :: Lens' env Connection

data SQLEnv = SQLEnv
    { sqlLogFunc    :: LogFunc
    , sqlConnection :: Connection }

newtype SQL a = SQL { unSQL :: RIO SQLEnv a }  -- ReaderT Connection (LoggingT IO) a }
    deriving (Applicative, Functor, Monad, MonadIO, MonadThrow, MonadReader SQLEnv)

instance HasLogFunc SQLEnv where
    logFuncL = lens sqlLogFunc (\x y -> x { sqlLogFunc = y })

instance HasConnection SQLEnv where
    connection = lens sqlConnection (\x y -> x { sqlConnection = y })

getConnection :: (HasConnection env, MonadReader env m) => m Connection
getConnection = view connection

db :: (MonadIO m, MonadReader env m, HasConnection env, HasLogFunc env)
   => SQL a -> m a
db (SQL x) = do
    logFunc <- view logFuncL
    conn    <- view connection
    runRIO (SQLEnv logFunc conn) x

runSQL' :: (MonadIO m) => SQLEnv -> SQL a -> m a
runSQL' env (SQL x) = runRIO env x

expectUnique :: (MonadThrow m, Show a) => [a] -> m (Maybe a)
expectUnique []  = return Nothing
expectUnique [x] = return $ Just x
expectUnique lst = throwM $ DatabaseError $ "duplicate entry: " <> tshow lst

expectUnique' :: (MonadThrow m, Show a) => (a -> b) -> [a] -> m (Maybe b)
expectUnique' _ []  = return Nothing
expectUnique' f [x] = return $ Just (f x)
expectUnique' _ lst = throwM $ DatabaseError $ "duplicate entry: " <> tshow lst
```

The `SQLite.Simple` function `withTransaction` takes an `IO` action as argument. We somehow have to redirect logging information around the unpacking to `IO` and lifting back to `MonadSQL`. This is not the prettiest solution, and we see some repetition of the pattern where we unpack result and log, forward log to outer monad and return result pattern.

All the `RedirectLog` code is needed to aid the type checker, or it won't know what to do.

``` {.haskell #database-types}
withTransactionM :: SQL a -> SQL a
withTransactionM t = do
    env <- ask
    conn <- getConnection
    liftIO $ withTransaction conn (runSQL' env t)
```

``` {.haskell #database-imports}
import Document
import Config
import Select (select)
```

::: {.note}
use `pragma foreign_keys = true`.
:::

A diagram for this schema, generated with `sqleton`. It would be awesome if the following would work.

``` {.make #-knit- .-hidden-}
schema.svg: <<file|schema>>
    cat $< | sqlite3 database
    sqleton -o $@ database -e -L circo 
```

![Schema](schema.svg)

In `SQLite.Simple` the above schema becomes

``` {.sqlite file=data/schema.sql}
-- rules when editing this schema:
-- * always use double quotes for identifiers,
-- * align types for easy reading
pragma synchronous = off;
pragma journal_mode = memory;
pragma foreign_keys = on;

<<schema>>
-- vim:ft=mysql
```

::: {.TODO}
Implement the interface in Selda.
:::

## Create database

``` {.haskell #database-create}
schema :: IO [Query]
schema = do
    schema_path <- getDataFileName "data/schema.sql"
    qs <- initMaybe . T.split (== ';') <$> T.IO.readFile schema_path
    return $ maybe [] (map Query) qs

createTables :: (MonadIO m, MonadReader env m, HasConnection env) => m ()
createTables = do
    conn <- getConnection
    liftIO $ schema >>= mapM_ (execute_ conn)
```

## Insertion

### Documents

``` {.sqlite #schema}
-- this table should be sorted on order of inclusion
create table if not exists "documents"
    ( "id"        integer primary key autoincrement
    , "filename"  text not null
    , "time"      timestamp default current_timestamp not null
    );
```

### Codes

The code table maps to a `ReferencePair` containing both `CodeBlock` and `ReferenceId`. I use the `name`, `ordinal` pair as primary key.

``` {.sqlite #schema}
create table if not exists "codes"
    ( "id"         integer primary key autoincrement
    , "name"       text not null
    , "ordinal"    integer not null
    , "source"     text not null
    , "language"   text not null
    , "document"   integer not null
    , "linenumber" integer
    , foreign key ("document") references "documents"("id")
    );
```

``` {.haskell #database-insertion}
insertCode' :: (Text, Int, Text, Maybe Text, Int64, Maybe Int) -> [CodeProperty] -> SQL ()
insertCode' tup attrs =  do
    conn <- getConnection
    liftIO $ do
        execute conn "insert into `codes`(`name`,`ordinal`,`source`,`language`,`document`,`linenumber`) \
                     \ values (?,?,?,?,?,?)" tup
        codeId <- lastInsertRowId conn
        executeMany conn "insert into `classes` values (?,?)"
                    [(className, codeId) | className <- getClasses attrs]
        executeMany conn "insert into `attributes` values (?,?,?)"
                    [(attr, value, codeId) | (attr, value) <- getAttrs attrs]
    where getClasses [] = []
          getClasses (CodeClass c : cs) = c : getClasses cs
          getClasses (_ : cs) = getClasses cs

          getAttrs [] = []
          getAttrs (CodeAttribute attr val : cs) = (attr, val) : getAttrs cs
          getAttrs (_ : cs) = getAttrs cs

insertCode :: Int64 -> ReferencePair -> SQL ()
insertCode docId ( ReferenceId file (ReferenceName name) count
                 , CodeBlock lang attrs source linenum ) = do
    langName <- case lang of
        UnknownClass c  -> logWarn (display $ "unknown language `" <> c <> "` in " <> location)
                        >> return (Just "<unknown>")
        NoLanguage      -> logWarn (display $ "no language class in " <> location)
                        >> return (Just "<unknown>") -- cannot use null language
        KnownLanguage l -> return (Just l)
    insertCode' (name, count, source, langName, docId, linenum) attrs
  where
    location = T.pack file <> ":<<" <> name <> ">>" <> T.pack (maybeLineNo linenum)
    maybeLineNo Nothing    = ""
    maybeLineNo (Just lno) = " at #" <> show lno

insertCodes :: Int64 -> ReferenceMap -> SQL ()
insertCodes docId codes = mapM_ (insertCode docId) (M.toList codes)
```

A table that references specific code blocks should reference both the code name and the code ordinal.

``` {.sqlite #reference-code}
, "code"        text not null
, foreign key ("code") references "codes"("id") on delete cascade
```

### Classes and attributes

``` {.sqlite #schema}
create table if not exists "classes"
    ( "class"     text not null
    <<reference-code>>
    );

create table if not exists "attributes"
    ( "attribute"   text not null
    , "value"       text not null
    <<reference-code>>
    );
```



### Content

``` {.sqlite #schema}
create table if not exists "content"
    ( "id"          integer primary key autoincrement
    , "document"    integer not null
    , "plain"       text
    , "code"        integer
    , foreign key ("document") references "documents"("id")
    , foreign key ("code") references "codes"("id")
    );
    -- , check ("plain" is not null or ("codeName" is not null and "codeOrdinal" is not null)) )
```

``` {.haskell #database-insertion}
queryCodeId' :: Int64 -> ReferenceId -> SQL (Maybe Int64)
queryCodeId' docId (ReferenceId _ (ReferenceName name) count) = do
    conn <- getConnection
    expectUnique' fromOnly =<< liftIO (query conn codeQuery (docId, name, count))
    where codeQuery = "select `id` from `codes` where \
                      \ `document` is ? and `name` is ? and `ordinal` is ?"

queryCodeId :: ReferenceId -> SQL (Maybe Int64)
queryCodeId (ReferenceId doc (ReferenceName name) count) = do
    conn    <- getConnection
    expectUnique' fromOnly =<< liftIO (query conn codeQuery (doc, name, count))
    where codeQuery = "select `codes`.`id` \
                      \ from `codes` inner join `documents` \
                      \     on `documents`.`id` is `codes`.`document` \
                      \ where   `documents`.`filename` is ? \
                      \     and `name` is ? \
                      \     and `ordinal` is ?"

queryReferenceCount :: ReferenceName -> SQL Int
queryReferenceCount (ReferenceName name) = do
    conn <- getConnection
    result <- liftIO (query conn codeQuery (Only name))
    case result of
        [Only a] -> return a
        _        -> return 0
    where codeQuery = "select count(`id`) from `codes` where `name` is ?"

queryCodeSource :: ReferenceId -> SQL (Maybe Text)
queryCodeSource (ReferenceId doc (ReferenceName name) count) = do
    conn    <- getConnection
    expectUnique' fromOnly =<< liftIO (query conn codeQuery (doc, name, count))
    where codeQuery = "select `codes`.`source` \
                      \ from `codes` inner join `documents` \
                      \     on `documents`.`id` is `codes`.`document` \
                      \ where   `documents`.`filename` is ? \
                      \     and `name` is ? \
                      \     and `ordinal` is ?"

queryCodeAttr :: ReferenceName -> Text -> SQL (Maybe Text)
queryCodeAttr (ReferenceName name) attr = do
    conn    <- getConnection
    expectUnique' fromOnly =<< liftIO (query conn codeQuery (name, attr))
    where codeQuery = "select `attributes`.`value` \
                      \ from `attributes` inner join `codes` \
                      \     on `codes`.`id` is `attributes`.`code` \
                      \ where `name` is ? and `ordinal` is 0 \
                      \     and `attribute` is ?"

contentToRow :: Int64 -> Content -> SQL (Int64, Maybe Text, Maybe Int64)
contentToRow docId (PlainText text) = return (docId, Just text, Nothing)
contentToRow docId (Reference ref)  = do
    codeId' <- queryCodeId' docId ref
    case codeId' of
        Nothing     -> throwM $ DatabaseError $ "expected reference: " <> tshow ref
        Just codeId -> return (docId, Nothing, Just codeId)

insertContent :: Int64 -> [Content] -> SQL ()
insertContent docId content = do
    rows <- mapM (contentToRow docId) content
    conn <- getConnection
    liftIO $ executeMany conn "insert into `content`(`document`,`plain`,`code`) values (?,?,?)" rows
```

### Targets

``` {.sqlite #schema}
create table if not exists "targets"
    ( "filename"  text primary key
    , "codename"  text not null
    , "language"  text not null
    , "document"  integer          -- in case this is null, the target is orphaned
    , "time"      timestamp default current_timestamp not null
    -- , foreign key ("codename") references "codes"("name")
    , foreign key ("document") references "documents"("id")
    );
```

``` {.haskell #database-insertion}
insertTargets :: Int64 -> FileMap -> SQL ()
insertTargets docId files = do
        conn <- getConnection
        -- liftIO $ print rows
        liftIO $ executeMany conn "replace into `targets`(`filename`,`codename`,`language`,`document`) values (?, ?, ?, ?)" rows
    where targetRow (path, (ReferenceName name, lang)) = (path, name, lang, docId)
          rows = map targetRow (M.toList files)
```

## Updating

### Tangle

When a Markdown source is written to, we first remove all data that was inserted from that document and then reinsert the new data. This is to prevent that code blocks remain from older versions. This problem could also be solved by a separate garbage collection step.

``` {.haskell #database-update}
getDocumentId :: FilePath -> SQL (Maybe Int64)
getDocumentId rel_path = do
        conn <- getConnection
        expectUnique' fromOnly =<< liftIO (query conn documentQuery (Only rel_path))
    where documentQuery = "select `id` from `documents` where `filename` is ?"

orphanTargets :: Int64 -> SQL ()
orphanTargets docId = do
    conn <- getConnection
    liftIO $ execute conn "update `targets` set `document` = null where `document` is ?" (Only docId)

clearOrphanTargets :: SQL ()
clearOrphanTargets = do
    conn <- getConnection
    liftIO $ execute_ conn "delete from `targets` where `document` is null"

removeDocumentData :: Int64 -> SQL ()
removeDocumentData docId = do
    orphanTargets docId
    conn <- getConnection
    liftIO $ do
        execute conn "delete from `content` where `document` is ?" (Only docId)
        execute conn "delete from `codes` where `document` is ?" (Only docId)

insertDocument :: FilePath -> Document -> SQL ()
insertDocument rel_path Document{..} = do
    let refNames = nub $ map referenceName $ M.keys references
    conn <- getConnection
    docId' <- getDocumentId rel_path
    docId <- case docId' of
        Just docId -> do
            logDebug $ display $ "Replacing '" <> T.pack rel_path <> "'."
            removeDocumentData docId >> return docId
        Nothing    -> do
            logDebug $ display $ "Inserting new '" <> T.pack rel_path <> "'."
            liftIO $ execute conn "insert into `documents`(`filename`) values (?)" (Only rel_path)
            liftIO $ lastInsertRowId conn
    refCountMap <- M.fromList . zip refNames
                <$> mapM queryReferenceCount refNames
    let mapref r@ReferenceId{..}
            = maybe r (\c -> r {referenceCount=referenceCount+c})
                    (refCountMap M.!? referenceName)
        refs = M.mapKeys mapref references
        cont = map (\case { Reference rid -> Reference (mapref rid);
                            x             -> x }) documentContent
    insertCodes docId refs
    insertContent docId cont
    insertTargets docId documentTargets
```

### Stitch

::: {.TODO}
Make this update more efficient. See https://stackoverflow.com/questions/11563869/update-multiple-rows-with-different-values-in-a-single-sql-query for an example.
Other performance related post: https://stackoverflow.com/questions/1711631/improve-insert-per-second-performance-of-sqlite
:::

``` {.haskell #database-update}
fromMaybeM :: Monad m => m a -> m (Maybe a) -> m a
fromMaybeM whenNothing x = maybe whenNothing return =<< x

updateCode :: ReferencePair -> SQL ()
updateCode (ref, CodeBlock {codeSource}) = do
    codeId <- fromMaybeM (throwM $ DatabaseError $ "could not find code '" <> tshow ref <> "'")
                         (queryCodeId ref)
    conn   <- getConnection
    liftIO $ execute conn "update `codes` set `source` = ? where `id` is ?" (codeSource, codeId)

updateTarget :: [ReferencePair] -> SQL ()
updateTarget refs = withTransactionM $ mapM_ updateCode refs
```

## Queries

``` {.haskell #database-queries}
listOrphanTargets :: SQL [FilePath]
listOrphanTargets = do
    conn <- getConnection
    map fromOnly <$>
        liftIO (query_ conn "select `filename` from `targets` where `document` is null" :: IO [Only FilePath])

listTargetFiles :: SQL [FilePath]
listTargetFiles = do
    conn <- getConnection
    map fromOnly <$>
        liftIO (query_ conn "select `filename` from `targets` where `document` is not null" :: IO [Only FilePath])
```

``` {.haskell #database-queries}
listSourceFiles :: SQL [FilePath]
listSourceFiles = do
    conn <- getConnection
    map fromOnly <$>
        liftIO (query_ conn "select `filename` from `documents`" :: IO [Only FilePath])
```

``` {.haskell #database-queries}
queryTargetRef :: FilePath -> SQL (Maybe (ReferenceName, Text))
queryTargetRef rel_path = do
    conn <- getConnection
    val  <- expectUnique =<< liftIO (query conn targetQuery (Only rel_path))
    return $ case val of
        Nothing           -> Nothing
        Just (name, lang) -> Just (ReferenceName name, lang)
    where targetQuery = "select `codename`,`language` from `targets` where `filename` is ?"
```

``` {.haskell #database-queries}
stitchDocument :: FilePath -> SQL Text
stitchDocument rel_path = do
    conn <- getConnection
    docId <- getDocumentId rel_path >>= \case
        Nothing -> throwM $ StitchError $ "File `" <> T.pack rel_path <> "` not in database."
        Just x  -> return x
    result <- liftIO (query conn "select coalesce(`plain`,`codes`.`source`) from `content` \
                                 \  left outer join `codes` \
                                 \    on `code` is `codes`.`id` \
                                 \  where `content`.`document` is ?" (Only docId) :: IO [Only Text])
    return $ T.unlines $ map fromOnly result
```

### References

``` {.haskell #database-queries}
queryReferenceMap :: Config -> SQL ReferenceMap
queryReferenceMap cfg = do
        conn <- getConnection
        rows <- liftIO (query_ conn "select `documents`.`filename`, `name`, `ordinal`, `source`, `linenumber`, `language`\
                                    \from `codes` inner join `documents` on `codes`.`document` is `documents`.`id`"
                        :: IO [(FilePath, Text, Int, Text, Maybe Int, Text)])
        M.fromList <$> mapM refpair rows
    where refpair (rel_path, name, ordinal, source, linenum, lang) = do
            let lang' = maybe (UnknownClass lang) (const $ KnownLanguage lang)
                              (languageFromName cfg lang)
            return ( ReferenceId rel_path (ReferenceName name) ordinal
                   , CodeBlock lang' [] source linenum )
```

## Deduplication
When a target file that gets writen to has multiple copies of the same reference, it may not be obvious which of them got changed. The best we can do, is check that only one of the new versions is different from the one in the markdown source.

The function `deduplicateRefs` takes a list of reference pairs (read from a target file), and uses the information in the database to make sure that the resulting list of reference pairs is unique by `ReferenceId`. If no decision is reached, a `StitchError` is raised.

``` {.haskell #database-deduplicate}
deduplicateRefs :: [ReferencePair] -> SQL [ReferencePair]
deduplicateRefs refs = dedup $ sortOn fst refs
    where dedup [] = return []
          dedup [x1] = return [x1]
          dedup ((ref1, code1@CodeBlock{codeSource=s1}) : (ref2, code2@CodeBlock{codeSource=s2}) : xs)
                | ref1 /= ref2 = ((ref1, code1) :) <$> dedup ((ref2, code2) : xs)
                | s1 == s2     = dedup ((ref1, code1) : xs)
                | otherwise    = do
                    old_code <- queryCodeSource ref1
                    case old_code of
                        Nothing -> throwM $ StitchError $ "ambiguous update: " <> tshow ref1 <> " not in database."
                        Just c  -> select (throwM $ StitchError $ "ambiguous update to " <> tshow ref1)
                                    [(s1 == c && s2 /= c, dedup ((ref2, code2) : xs))
                                    ,(s1 /= c && s2 == c, dedup ((ref1, code1) : xs))]
```

### Document
# Text utilities

``` {.haskell #import-text}
import RIO (Text)
import qualified RIO.Text as T
```

We need a version of `unlines` that doesn't append a final newline. This also means changing `T.lines` a little.

``` {.haskell #unlines}
lines' :: Text -> [Text]
lines' text
    | text == ""               = [""]
    | "\n" `T.isSuffixOf` text = T.lines text <> [""]
    | otherwise                = T.lines text

unlines' :: [Text] -> Text
unlines' = T.intercalate "\n"
```

This way, `unlines'` has nicer properties. No single `Text` in Entangled ends in a newline, meaning the `unlines'` operation is associative:

``` {.haskell}
unlines' [unlines' a, unlines' b] = unlines' (a <> b)
```

With the notable exception of empty lists. The `[Text]` type is a proper monoid, however when we consider `Text` itself in the context of an element of an `unlines`ed text, the zero element for `Text`, being `""` changes meaning. An empty line is something different than the absence of text. To preserve both associativity and make the `unlines . lines = id` property more strict we need to wrap `Text` in a `Maybe`.

``` {.haskell #maybe-unlines}
mLines :: Maybe Text -> [Text]
mLines Nothing = []
mLines (Just text) = lines' text

mUnlines :: [Text] -> Maybe Text
mUnlines [] = Nothing
mUnlines ts = Just $ unlines' ts
```

Now we can test:

``` {.haskell #test-unlines-inverse}
t == mUnlines (mLines t)
```

and

``` {.haskell #test-unlines-associative}
mUnlines (catMaybes [mUnlines a, mUnlines b]) == mUnlines (a <> b)
```

### `TextUtils` module 

``` {.haskell file=src/TextUtil.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module TextUtil where

import RIO
import qualified RIO.Text as T

import Data.Char (isSpace)

<<indent>>
<<unindent>>
<<unlines>>
<<maybe-unlines>>
```

``` {.haskell #indent}
indent :: Text -> Text -> Text
indent pre text
    = unlines' $ map indentLine $ lines' text
    where indentLine line
            | line == "" = line
            | otherwise  = pre <> line
```

``` {.haskell #unindent}
unindent :: Text -> Text -> Maybe Text
unindent prefix s
    = unlines' <$> mapM unindentLine (lines' s)
    where unindentLine t
            | T.all isSpace t = Just ""
            | otherwise       = T.stripPrefix prefix t
```

#### Tests

``` {.haskell file=test/TextUtilSpec.hs}
module TextUtilSpec where

<<import-text>>
import Data.Maybe (catMaybes, isJust)
import Test.Hspec
import Test.QuickCheck
import Test.QuickCheck.Instances.Text ()

import TextUtil

propUnlines :: Maybe Text -> Bool
propUnlines t =
    <<test-unlines-inverse>>

propUnlineLists :: ([Text], [Text]) -> Bool
propUnlineLists (a, b) =
    <<test-unlines-associative>>

genLine :: Gen Text
genLine = T.pack <$> (listOf $ elements ['!'..'~'])

genText :: Gen Text
genText = unlines' <$> listOf genLine

genPair :: Gen a -> Gen b -> Gen (a, b)
genPair x y = do
    i <- x
    j <- y
    return (i, j)

propIndent :: (Text, Text) -> Bool
propIndent (a, b) = unindent a (indent a b) == Just b

propUnindentFail :: (Text, Text) -> Bool
propUnindentFail (a, b) = (isJust $ unindent a b) == a `T.isPrefixOf` b

textUtilSpec :: Spec
textUtilSpec = do
    describe "property check" $ do
        it "mUnlines inverses mLines" $
            property $ propUnlines
        it "mUnlines can be nested (associativity)" $
            property $ propUnlineLists
        it "unindent inverts indent" $
            property $ forAll (genPair genLine genText)  propIndent
        it "unindent fails on wrong indent" $
            property $ forAll (genPair genLine genLine) propIndent
```

# Linting

``` {.haskell file=src/Linters.hs}
{-# LANGUAGE NoImplicitPrelude #-}
{-# LANGUAGE ScopedTypeVariables #-}
module Linters where

import RIO
import qualified RIO.Text as T
import RIO.List (repeat, zip3)
import RIO.Set (elems)
import qualified RIO.Set as S
import qualified RIO.Map as M

import Data.Graph (Graph, Vertex)
import qualified Data.Graph as G

import Database.SQLite.Simple (query_, query, fromOnly, Only(..))

import Entangled
import Document (ReferenceMap, ReferenceName(..), codeBlocksByName, CodeBlock(..), referenceNames)
import Database (HasConnection, db, queryReferenceMap, getConnection)
import Daemon (Session)
import Config (HasConfig, config)
import Tangle (parseCode, CodeLine(..))
import FileIO

referencesInFragment :: ReferenceMap -> ReferenceName -> [ReferenceName]
referencesInFragment refs r =
    concatMap (mapMaybe getReference . parseCode r . codeSource)
              (codeBlocksByName refs r)
    where getReference (NowebReference s _) = Just s
          getReference _                    = Nothing

referenceGraph :: (HasConnection env, HasLogFunc env, HasConfig env)
               => Entangled env (Graph, Vertex -> ((), ReferenceName, [ReferenceName]), ReferenceName -> Maybe Vertex)
referenceGraph = do
    cfg <- view config
    refs <- db (queryReferenceMap cfg)
    let names = elems $ referenceNames refs
    return $ G.graphFromEdges $ zip3 (repeat ()) names (map (referencesInFragment refs) names)

dumpToGraphViz :: (HasConnection env, HasLogFunc env, HasConfig env)
               => Entangled env ()
dumpToGraphViz = do
    (graph, vertexToNode, _) <- referenceGraph
    let graphvizEdges = T.unlines $ map edgeToGV $ G.edges graph
        nodeToKey (_, ReferenceName b, _) = "\"" <> b <> "\""
        edgeToGV (a, b) = nodeToKey (vertexToNode a) <> " -> " <> nodeToKey (vertexToNode b)
    dump $ "digraph {\n" <> graphvizEdges <> "}"

listUnusedFragments :: (HasConnection env, HasLogFunc env, HasConfig env)
                    => Entangled env ()
listUnusedFragments = do
    (graph, vertexToNode, keyToVertex) <- referenceGraph
    conn <- getConnection
    roots <- mapMaybe (keyToVertex . ReferenceName . fromOnly)
          <$> liftIO (query_ conn "select `codename` from `targets` \
                                  \where `document` is not null" :: IO [Only Text])
    let reach = foldMap (S.fromList . G.reachable graph) roots
        unused = S.fromList (G.vertices graph) S.\\ reach
        nodeToKey (_, ReferenceName b, _) = b
        getSource name = liftIO $ query conn "select `codes`.`name`,`documents`.`filename` from `codes` \
                                             \inner join `documents` on `codes`.`document` is `documents`.`id`\
                                             \where `codes`.`name` is ? and `codes`.`ordinal` is 0" (Only name)
    codes <- foldMapM getSource (map (nodeToKey . vertexToNode) (S.toList unused))
    dump $ T.unlines $ map (\(name, src) -> name <> " in " <> src) codes

linters :: (HasConnection env, HasLogFunc env, HasConfig env)
        => Map Text (Entangled env ())
linters = M.fromList
    [ ("dumpToGraphViz", dumpToGraphViz)
    , ("listUnusedFragments", listUnusedFragments) ]

allLinters :: [Text]
allLinters  = M.keys (linters :: Map Text (Entangled Session ()))

lint :: (HasConnection env, HasLogFunc env, HasConfig env)
     => [Text] -> Entangled env ()
lint [] | allLinters /= [] = lint allLinters -- use all available linters
lint args = sequence_ $ mapMaybe (linters M.!?) args
```
# Untangling

```{.haskell file=src/Stitch.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Stitch where

import RIO hiding (some)
import qualified RIO.Text as T
import RIO.List.Partial (head, tail)

<<stitch-imports>>
<<source-parser>>
<<stitch>>
```

Untangling starts with reading the top line to identify the file and language. Following that should be one series of referenced code block items.

The result is a list of `ReferencePair` giving all the content of the code blocks as it is given in the source file. We don't yet make it a `ReferenceMap`, since there may be duplicate conflicting entries. In this case we want to get the entry that is different from the one that we already know of.

``` {.haskell #source-parser}
sourceDocument :: ( MonadParsec e (ListStream Text) m
                  , MonadFail m
                  , MonadReader Config m )
               => m [ReferencePair]
sourceDocument = do
    cfg <- ask
    (header, (prop, _)) <- manyTill_ anySingle (tokenP topHeader)
    -- (prop, _) <- tokenP topHeader
    lang <- maybe (fail "No valid language found in header.") return
                  $ getAttribute prop "language" >>= languageFromName cfg
    (_, refs) <- mconcat <$> some (sourceBlock lang)
    return (prependCode header (head refs) : tail refs)
    where prependCode h (ref, cb@CodeBlock{..})
            = (ref, cb {codeSource = unlines' (h <> [codeSource])})
```

A `sourceBlock` starts with a *begin* marker, then has many lines of plain source or nested `sourceBlock`s. Both `sourceBlock` and `sourceLine` return pairs of texts and references. The content of these pairs are concatenated. If a `sourceBlock` is the first in a series (index 0), the noweb reference is generated with the correct indentation.

``` {.haskell #source-parser}
sourceBlock :: ( MonadReader Config m, MonadParsec e (ListStream Text) m, MonadFail m )
            => ConfigLanguage -> m ([Text], [ReferencePair])
sourceBlock lang = do
    Config{..} <- ask
    ((ref, beginIndent), _) <- tokenP (commented lang beginBlock)
    when configUseLineDirectives (void anySingle)
    (ilines, refpairs) <- mconcat <$> manyTill
                (sourceBlock lang <|> sourceLine)
                (tokenP (commented lang endBlock))
    let unindentedLines = map (unindent beginIndent) ilines
    when (any isNothing unindentedLines) $ fail "Indentation error"
    let content = unlines' $ catMaybes unindentedLines
    return ( [ indent beginIndent $ showNowebReference $ referenceName ref
             | referenceCount ref == 0 ]
           , (ref, CodeBlock (KnownLanguage $ languageName lang) [] content Nothing):refpairs )

sourceLine :: ( MonadParsec e (ListStream Text) m )
           => m ([Text], [ReferencePair])
sourceLine = do
    x <- anySingle
    return ([x], [])
```

## The `stitch` function

In the `stitch` function we take out the `Config` from the anonymous `MonadReader` and put it in a `SourceParser` monad. This transformation is the `asks . runReaderT` combo. It seems silly that we can't "inherit" the outer monad here. I tried turning the transformers around like: `ParsecT Void (ListStream Text) m`, but type deduction fails on that one.

``` {.haskell #stitch}
type SourceParser = ReaderT Config (Parsec Void (ListStream Text))

untangle :: ( MonadReader env m, HasConfig env, MonadThrow m )
         => FilePath -> Text -> m [ReferencePair]
untangle file text = do
    doc <- runReaderT (sourceDocument :: SourceParser [ReferencePair]) <$> view config
    either (throwM . StitchError . T.pack . errorBundlePretty)
           return $ parse doc file $ ListStream (T.lines text)
```

## Imports

``` {.haskell #stitch-imports}
import ListStream (ListStream(..), tokenP)
import Document
import Config (config, HasConfig, Config(..), languageFromName, ConfigLanguage(..))
import Comment (topHeader, beginBlock, endBlock, commented)
import TextUtil (indent, unindent, unlines')

import Text.Megaparsec
    ( MonadParsec, Parsec, parse, anySingle, manyTill, some, errorBundlePretty, manyTill_ )
```
# Transactions
The idea here is that we'd like to create an abstract class to do all our file system interactions. For this we define a `Transaction` that encodes a monadic action, together with a description of that action and the possibility of needing human confirmation on the transaction.

``` {.haskell file=src/Transaction.hs}
{-# LANGUAGE NoImplicitPrelude,UndecidableInstances #-}
module Transaction where

import RIO
import qualified RIO.Text as T
<<transaction-imports>>

<<transaction>>
```

``` {.haskell #transaction}
data Description = Message Doc
                 | CreateFile FilePath
                 | WriteFile FilePath
                 | DeleteFile FilePath

data Transaction m = Transaction
  { action :: Maybe (m ())
  , description :: [Description]
  , needConfirm :: Bool }
```

Here `Doc` is pretty-printing class from `Console`.
Now, considering the different options for displaying a transaction, we need a human readable and a machine readable description

``` {.haskell #transaction}
showDescriptionHuman :: [Description] -> Doc
showDescriptionHuman = mconcat . map (\case
    Message x    -> x
    CreateFile f -> msgCreate f
    WriteFile f  -> msgWrite f
    DeleteFile f -> msgDelete f)

showDescriptionMachine :: [Description] -> Text
showDescriptionMachine = T.unlines . mapMaybe (\case
    Message _    -> Nothing
    CreateFile f -> Just $ "+ " <> T.pack f
    WriteFile f  -> Just $ "~ " <> T.pack f
    DeleteFile f -> Just $ "- " <> T.pack f)
```

When an event happened we need to respond, usually by writing out several files. Since `IO` is a `Monoid`, we can append `IO` actions and keep track of describing the gathered events in a `Transaction`. There are some things that we may need to ask the user permission for, like overwriting files in dubious circumstances. Messaging is done through pretty-printed `Doc`.

``` {.haskell #transaction-imports}
import Console (Doc, group, msgCreate, msgDelete, msgWrite)
import qualified Console
import qualified Data.Text.IO as T.IO
```

The `action` is wrapped in a `Maybe` so that we can tell if the `Transaction` does anything. A `Transaction` is a `Monoid`.

``` {.haskell #transaction}
instance (Semigroup (m ())) => Semigroup (Transaction m) where
    (Transaction al dl cl) <> (Transaction ar dr cr)
      = Transaction (al <> ar) (dl <> dr) (cl || cr)

instance (Monoid (m ())) => Monoid (Transaction m) where
    mempty = Transaction mempty mempty False
```

We can build `Transaction`s by appending elemental parts.

``` {.haskell #transaction}
plan :: (Monad m) => Description -> m () -> Transaction m
plan d action = Transaction (Just action) [d] False

doc :: Doc -> Transaction m
doc x = Transaction Nothing [Message x] False

confirm :: Transaction m
confirm = Transaction Nothing mempty True
```

In most of the program logic, `Transaction` will be available in terms of a `MonadWriter`.

``` {.haskell #transaction}
testTransaction :: (MonadIO m) => Transaction m -> m Bool
testTransaction (Transaction Nothing _ _)  = return False
testTransaction (Transaction (Just _) d _) = liftIO $ T.IO.putStr (showDescriptionMachine d) >> return True

runTransactionMachine :: (MonadIO m) => Transaction m -> m ()
runTransactionMachine (Transaction Nothing d _) = liftIO $ T.IO.putStr (showDescriptionMachine d)
runTransactionMachine (Transaction (Just x) d _) = do
    liftIO $ T.IO.putStr (showDescriptionMachine d)
    x

runTransaction :: (MonadIO m) => Maybe Doc -> Transaction m -> m ()
runTransaction (Just h) (Transaction Nothing d _) = liftIO $ Console.putTerminal $ group h (showDescriptionHuman d)
runTransaction Nothing (Transaction Nothing d _) = liftIO $ Console.putTerminal (showDescriptionHuman d)
runTransaction h (Transaction (Just x) d c) = do
    liftIO $ Console.putTerminal $ maybe (showDescriptionHuman d) (`group` showDescriptionHuman d) h
    if c then do
        reply <- liftIO $ do
            T.IO.putStr "confirm? (y/n) "
            hFlush stdout
            T.IO.getLine
        liftIO $ T.IO.putStrLn ""
        when (reply == "y") x
    else x
```

## File transactions

There is a limited set of IO file system actions that result from a tangle or stitch. We define a little language using a type class.

- When a target gets created or modified, we need to `writeFile`.
- When a target is removed or renamed, we need to `deleteFile`.
- For every tangle and stitch operation we need to `readFile`.

The behaviour will be such that, if a file is deleted and no other file remains in its containing directory, the directory is removed. If we write a file to a directory that does not exists, the directory is created.

``` {.haskell file=src/FileIO.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module FileIO where

import RIO
<<file-io-imports>>

class Monad m => MonadFileIO m where
    writeFile :: FilePath -> Text -> m ()
    deleteFile :: FilePath -> m ()
    readFile :: FilePath -> m Text
    dump :: Text -> m ()

<<file-io-prim>>
<<file-io-instance>>
```

## Primitives

``` {.haskell #file-io-imports}
import qualified RIO.Text as T

import Errors (EntangledError(SystemError))
import Select (selectM)
```

### Make/Remove dir

``` {.haskell #file-io-imports}
import RIO.Directory ( createDirectoryIfMissing, doesDirectoryExist
                     , listDirectory, removeFile, removeDirectory )
import RIO.FilePath  ( (</>), splitDirectories, takeDirectory )
import RIO.List      ( scanl1 )
```

``` {.haskell #file-io-prim}
ensurePath :: (MonadIO m, MonadReader env m, HasLogFunc env)
           => FilePath -> m ()
ensurePath path = selectM (return ())
    [ ( not <$> doesDirectoryExist path
      , logDebug (display $ "creating directory `" <> T.pack path <> "`")
        >> createDirectoryIfMissing True path ) ]
```

``` {.haskell #file-io-prim}
rmDirIfEmpty :: (MonadIO m, MonadThrow m, MonadReader env m, HasLogFunc env)
             => FilePath -> m ()
rmDirIfEmpty path = selectM (return ())
    [ ( not <$> doesDirectoryExist path
      , throwM $ SystemError $ "could not remove dir: `" <> T.pack path <> "`")
    , ( null <$> listDirectory path
      , logDebug (display $ "removing empty directory `" <> T.pack path <> "`")
        >> removeDirectory path ) ]

parents :: FilePath -> [FilePath]
parents = scanl1 (</>) . splitDirectories

rmPathIfEmpty :: (MonadIO m, MonadThrow m, MonadReader env m, HasLogFunc env)
              => FilePath -> m ()
rmPathIfEmpty = mapM_ rmDirIfEmpty . reverse . parents
```

### Write file

``` {.haskell #file-io-imports}
import RIO.File ( writeBinaryFileDurable )
import qualified RIO.ByteString as B
```

``` {.haskell #file-io-prim}
writeIfChanged :: (MonadIO m, MonadThrow m, MonadReader env m, HasLogFunc env)
               => FilePath -> Text -> m ()
writeIfChanged path text = do
    old_content' <- liftIO $ try $ B.readFile path
    case (old_content' :: Either IOException B.ByteString) of
        Right old_content | old_content == new_content -> return ()
                          | otherwise                  -> write
        Left  _                                        -> write
    where new_content = T.encodeUtf8 text
          write       = logDebug (display $ "writing `" <> T.pack path <> "`")
                      >> writeBinaryFileDurable path new_content

dump' :: (MonadIO m, MonadReader env m, HasLogFunc env)
      => Text -> m ()
dump' text = logDebug "dumping to stdio"
         >> B.hPutStr stdout (T.encodeUtf8 text)
```

## FileIO instance

``` {.haskell #file-io-instance}
newtype FileIO env a = FileIO { unFileIO :: RIO env a }
    deriving (Applicative, Functor, Semigroup, Monoid, Monad, MonadIO, MonadThrow, MonadReader env)

readFile' :: ( MonadIO m, HasLogFunc env, MonadReader env m, MonadThrow m )
          => FilePath -> m Text
readFile' path = decodeUtf8With lenientDecode
               <$> (logDebug (display $ "reading `" <> T.pack path <> "`")
                    >> B.readFile path)

runFileIO' :: ( MonadIO m, MonadReader env m, HasLogFunc env )
          => FileIO env a -> m a
runFileIO' (FileIO f) = do
    env <- ask
    runRIO env f

instance (HasLogFunc env) => MonadFileIO (FileIO env) where
    writeFile path text = ensurePath (takeDirectory path)
                        >> writeIfChanged path text

    deleteFile path     = logDebug (display $ "deleting `" <> T.pack path <> "`")
                        >> removeFile path
                        >> rmPathIfEmpty (takeDirectory path)

    readFile            = readFile'

    dump                = dump'
```

These are IO actions that need logging, possible confirmation by the user and execution. Also, using this we can do some mock testing.
# Logging

A logging type class.
# Document structure

``` {.haskell file=src/Document.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Document
    ( module Document
    , module Errors ) where

import RIO
import RIO.List (sort)
import qualified RIO.Set as S
import qualified RIO.Map as M
import Errors

<<document-structure>>
```

## Structure

A document is modelled as a series of text and code blocks. A code block is delimited by two lines starting with three backtics:

~~~
 Normal text.

 ```
 print("this is a code block")
 ```
~~~

Each code block that is of relevance to Entangled has a reference attached.

~~~
 ``` {.language #<reference>}
 ```
~~~

or a filename.

~~~
 ``` {.language file=<path>}
 ```
~~~

An identifier may be repeated, in which case the code block is concatenated to the previous instances. Each instance will have an associated integer marking the place in the sequence.

### Reference Id
We define the type `ReferenceId`.

``` {.haskell #document-structure}
newtype ReferenceName = ReferenceName
    { unReferenceName :: Text
    } deriving (Show, Eq, Ord, Display)

data ReferenceId = ReferenceId
    { referenceFile       :: FilePath
    , referenceName       :: ReferenceName
    , referenceCount      :: Int
    } deriving (Show, Eq, Ord)

showNowebReference :: ReferenceName -> Text
showNowebReference (ReferenceName x) = "<<" <> x <> ">>"
```

The type is deriving `Ord`, so that we can use it as an index for a `Map`.

### Document
Using the `ReferenceId` type we can represent a text (or `Document`) as a sequence of plain `Text` or `ReferenceId`, paired up with a map linking each `ReferenceId` with a code block.

``` {.haskell #document-structure}
data Content
    = PlainText Text
    | Reference ReferenceId
    deriving (Show, Eq)

type ReferencePair = (ReferenceId, CodeBlock)
type ReferenceMap = Map ReferenceId CodeBlock
type FileMap = Map FilePath (ReferenceName, Text)

data Document = Document
    { references      :: ReferenceMap
    , documentContent :: [Content]
    , documentTargets :: FileMap
    } deriving (Show)
```

### Accessors

``` {.haskell #document-structure}
referenceNames :: ReferenceMap -> Set ReferenceName
referenceNames = S.fromList . map referenceName . M.keys

referencesByName :: ReferenceMap -> ReferenceName -> [ReferenceId]
referencesByName refs name
    = (sort . filter ((== name) . referenceName) . M.keys) refs

codeBlocksByName :: ReferenceMap -> ReferenceName -> [CodeBlock]
codeBlocksByName refs name = mapMaybe (refs M.!?) $ referencesByName refs name

mapReferences :: (ReferenceId -> ReferenceId) -> Document -> Document
mapReferences f (Document {..}) =
    Document (M.mapKeys f references)
             (map mapContent documentContent)
             documentTargets
    where mapContent (Reference rid) = Reference (f rid)
          mapContent x               = x
```

### Code blocks

A code block can have all attributes that are normally associated in the context of CSS. The following markdown,

~~~
 ``` {.class #id attribute=value}
 source
 ```
~~~

can be mapped to a list of `CodeProperty`.

``` {.haskell #document-structure}
data CodeProperty
    = CodeId Text
    | CodeAttribute Text Text
    | CodeClass Text
    deriving (Eq, Show)

getAttribute :: [CodeProperty] -> Text -> Maybe Text
getAttribute [] _ = Nothing
getAttribute (CodeAttribute k v:ps) l
    | k == l    = Just v
    | otherwise = getAttribute ps l
getAttribute (_:ps) l = getAttribute ps l
```

In our case the *class* represents the code language or an abbreviation thereof.

``` {.haskell #document-structure}
data CodeBlock = CodeBlock
    { codeLanguage   :: ProgrammingLanguage
    , codeProperties :: [CodeProperty]
    , codeSource     :: Text
    , codeLineNumber :: Maybe Int
    } deriving (Show, Eq)
```

A `ProgrammingLanguage` is either a known language (for which we know how to generate comments) an unknown (or misspelled) class identifier, or empty.

``` {.haskell #document-structure}
data ProgrammingLanguage
    = KnownLanguage Text
    | UnknownClass Text
    | NoLanguage
    deriving (Show, Eq)
```

# The Daemon
The Entangled daemon runs a Milkshake main loop.

## Strategy

``` {.haskell file=src/Daemon.hs #daemon}
{-# LANGUAGE NoImplicitPrelude,ScopedTypeVariables #-}
module Daemon where

import RIO
import RIO.List (sort)
import RIO.Writer (tell)
import qualified RIO.Text as T

<<daemon-imports>>

<<daemon-events>>
<<daemon-session>>
<<daemon-watches>>
<<daemon-main-loop>>
<<daemon-start>>
```

Using lenses we build a `Session` record. In combination with the `State` monad (and a `Reader Config`, and `IO`) we can manage the Entangled Daemon.

### FSNotify

We use `System.FSNotify` to setup watches on the files.

``` {.haskell #daemon-imports}
import qualified System.FSNotify as FSNotify
```

The watches are set to trigger events of type `Event`.

``` {.haskell #daemon-events}
data DaemonState
    = Idle
    | Tangling
    | Stitching
    deriving (Show, Eq)

data Event
    = WriteSource FilePath
    | WriteTarget FilePath
    | DebugEvent Text
    deriving (Show)
```

### Session

``` {.haskell #daemon-imports}
import Database.SQLite.Simple (Connection)
import Database (db, connection, HasConnection, listSourceFiles, listTargetFiles)

import Transaction (doc)
-- import FileIO
import Tangle (Annotator, selectAnnotator)
import Entangled
import Config (Config(..), HasConfig, config, getInputFiles, configWatchList, AnnotateMethod(..))
import Errors (EntangledError(..), formatError)

import Console (Doc, putTerminal)
import qualified Console
import qualified Prettyprinter as P

import Control.Monad.Except ( MonadError )
-- import Control.Concurrent.Chan
-- import Control.Concurrent
-- import Control.Monad.Catch
```

``` {.haskell #daemon-session}
data Session = Session
    { watches       :: MVar [FSNotify.StopListening]
    , manager       :: FSNotify.WatchManager
    , eventChannel  :: Chan Event
    , daemonState   :: MVar DaemonState
    , connection'   :: Connection
    , config'       :: Config
    , logFunc'      :: LogFunc
    }

instance HasConfig Session where
    config = lens config' (\x y -> x { config' = y })

instance HasLogFunc Session where
    logFuncL = lens logFunc' (\x y -> x { logFunc' = y })

instance HasConnection Session where
    connection = lens connection' (\x y -> x { connection' = y })

newtype Daemon a = Daemon { unDaemon :: RIO Session a }
    deriving ( Applicative, Functor, Monad, MonadIO, MonadReader Session, MonadThrow, MonadUnliftIO )
```

Every time an event happens we send it to `_eventChannel`. When we tangle we change the `_daemonState` to `Tangling`. Write events to target files are then not triggering a stitch. The other way around, if the daemon is in `Stitching` state, write events to the markdown source do not trigger a tangle. Because these events will arrive asynchronously, we use an `MVar` to change the state.

Setting the state is a bit involved (a bit more than with an `IORef`), but it guarantees safe use in a multi-threaded environment. An `MVar` can only be set if it is first emptied, the combined action can be done with `modifyMVar_`.

``` {.haskell #daemon-session}
setDaemonState :: DaemonState -> Daemon ()
setDaemonState s = do
    state <- asks daemonState
    modifyMVar_ state (const $ return s)
```

## Watching

``` {.haskell #daemon-imports}
import Data.List (nub)
-- import Control.Monad (mapM)
```

The `passEvent` function acts as the call-back for the FSNotify watcher. There are two strategies in which editors save files:

    * remove and create: Vim (by default), gedit and many more editors do this. Actually, the content is written to a temporary file, which is then renamed to the existing file. These operations are atomic so that no data is lost if the system crashes.
    * modify: VS Code does this.

We have to be flexible in how we interpret the incoming events. The most important bit is that we need to ignore the `delete` events. 

`FSNotify` lets us watch directories containing the files we're interested in. In `passEvent` we need to check if the event is actually on an involved file.

``` {.haskell #daemon-watches}
passEvent :: MVar DaemonState -> Chan Event
          -> [FilePath] -> [FilePath] -> FSNotify.Event -> IO ()
passEvent _      _       _    _    FSNotify.Removed {} = return ()
passEvent state' channel srcs tgts fsEvent = do
    abs_path <- canonicalizePath $ FSNotify.eventPath fsEvent
    state    <- readMVar state'

    let isSourceFile = any (equalFilePath abs_path) srcs
        isTargetFile = any (equalFilePath abs_path) tgts
        pass         = case state of
                         Idle      -> isSourceFile || isTargetFile
                         Tangling  -> isSourceFile
                         Stitching -> isTargetFile

    when pass $ do
        let etype = if isSourceFile then WriteSource else WriteTarget
        writeChan channel (etype abs_path)
```

``` {.haskell #daemon-watches}
setWatch :: Daemon ()
setWatch = do
    srcs <- db listSourceFiles >>= (liftIO . mapM canonicalizePath)
    tgts <- db listTargetFiles >>= (liftIO . mapM canonicalizePath)
    fsnotify <- asks manager
    channel  <- asks eventChannel

    let abs_dirs = nub $ map takeDirectory (srcs <> tgts)
    rel_dirs <- mapM makeRelativeToCurrentDirectory abs_dirs

    state <- asks daemonState
    stopActions <- liftIO $ mapM
        (\dir -> FSNotify.watchDir fsnotify dir (const True)
                                   (passEvent state channel srcs tgts))
        abs_dirs
    watchesMVar <- asks watches
    putMVar watchesMVar stopActions

    logDebug $ display $ "watching: " <> tshow rel_dirs
```

``` {.haskell #daemon-watches}
closeWatch :: Daemon ()
closeWatch = do
    stopActions <- takeMVar =<< asks watches
    liftIO $ sequence_ stopActions
    logDebug "suspended watches"
```

## Main loop

``` {.haskell #daemon-imports}
import RIO.Directory (makeRelativeToCurrentDirectory, canonicalizePath)
import RIO.FilePath (equalFilePath, takeDirectory)
```

The `mainLoop` is fed events and handles them.

``` {.haskell #daemon-main-loop}
tryEntangled :: (MonadReader env m, MonadUnliftIO m, MonadIO m, HasLogFunc env)
             => Maybe Doc -> Entangled env a -> m ()
tryEntangled msg action = catch (void $ runEntangledHuman msg action)
                                (\(err :: EntangledError) -> logError $ display $ formatError err)

mainLoop :: Event -> Daemon ()
<<main-loop-cases>>
```

After the first event we need to wait a bit, there may be more coming.

``` {.haskell #main-loop-cases}
mainLoop (WriteSource abs_path) = do
    rel_path <- makeRelativeToCurrentDirectory abs_path
    logDebug $ display $ "tangle triggered on `" <> T.pack rel_path <> "`"

    setDaemonState Tangling
    closeWatch

    tryEntangled (Just $ "tangling on `" <> P.pretty rel_path <> "`") $ do
        insertSources [rel_path]
        tangle TangleAll =<< getAnnotator
        clearOrphans

    setWatch
    setDaemonState Idle

mainLoop (WriteTarget abs_path) = do
    rel_path <- liftIO $ makeRelativeToCurrentDirectory abs_path
    logDebug $ display $ "stitch triggered on `" <> T.pack rel_path <> "`"
    setDaemonState Stitching
    closeWatch

    tryEntangled (Just $ "stitching on `" <> P.pretty rel_path <> "`") $ do
        insertTargets [rel_path]
        stitch StitchAll

    tryEntangled (Just $ "tangling after stitch `" <> P.pretty rel_path <> "`") $
        tangle TangleAll =<< getAnnotator

    setWatch
    setDaemonState Idle

mainLoop _ = return ()
```

## Initialisation

``` {.haskell #daemon-start}
printMsg :: Doc -> Daemon ()
printMsg = liftIO . Console.putTerminal

initSession :: Daemon ()
initSession = do
    cfg <- view config
    abs_paths <- sort <$> getInputFiles cfg
    when (null abs_paths) $ throwM $ SystemError "No input files."
    rel_paths <- mapM makeRelativeToCurrentDirectory abs_paths

    printMsg Console.banner
    runEntangledHuman (Just "initializing") $ do
        tell $ doc $
             P.align (P.vsep
                   $ map (Console.bullet
                         . (P.pretty ("Monitoring " :: Text) <>)
                         . Console.fileRead)
                           rel_paths)
             <> P.line

        insertSources rel_paths
        tangle TangleAll =<< getAnnotator

    setWatch

getAnnotator :: (HasConfig env, MonadReader env m, MonadIO m, MonadThrow m, MonadError EntangledError n)
             => m (Annotator n)
getAnnotator = do
    cfg <- view config
    when (configAnnotate cfg == AnnotateNaked) $
        throwM $ ConfigError "cannot run daemon with `Naked` annotation"
    return $ selectAnnotator cfg

runSession :: (HasConfig env, HasLogFunc env, HasConnection env, MonadReader env m, MonadIO m)
           => [FilePath] -> m ()
runSession inputFiles = do
    hSetBuffering stdout LineBuffering
    cfg' <- view config

    let cfg = cfg' { configWatchList = configWatchList cfg' <> map T.pack inputFiles }
    conn <- view connection
    logFunc <- view logFuncL
    fsnotify <- liftIO FSNotify.startManager
    channel <- newChan
    daemonState' <- newMVar Idle
    watches' <- newEmptyMVar

    let session = Session watches' fsnotify channel daemonState' conn cfg logFunc
    runRIO session $ unDaemon $ do
        initSession
        mapM_ mainLoop =<< getChanContents channel

    liftIO $ FSNotify.stopManager fsnotify
```
# Tangling

``` {.haskell file=src/Tangle.hs}
{-# LANGUAGE NoImplicitPrelude,ScopedTypeVariables #-}
module Tangle where

import RIO hiding (try, some, many)
import qualified RIO.Text as T
import qualified RIO.Map as M
<<import-lazy-map>>

import Text.Megaparsec
    ( MonadParsec, Parsec, parse, getOffset
    , chunk, many, some, eof
    , manyTill, anySingle, try, lookAhead, takeWhile1P, takeWhileP
    , (<?>) )
import Text.Megaparsec.Char
    ( space )

import Control.Monad.State (MonadState, gets, modify, StateT, evalStateT)
import Control.Monad.Except ( MonadError )

import ListStream
import Document
import Config (config, HasConfig, Config(..), lookupLanguage, ConfigLanguage(..), AnnotateMethod(..), ConfigSyntax(..))

<<tangle-imports>>

<<parse-markdown>>
<<generate-code>>
```

The task of tangling means:

* Parse the markdown to a `Document`
* Generate annotated source files.

Parsing the markdown is done on a line per line basis. We don't try to parse the markdown itself, rather we try to detect lines that start and end code-blocks.

Remember the golden rule:

> Untangling from a generated source returns the same markdown **to the byte**.

We will be adding unit tests to check exactly this property.

## CSS Attributes

We'll reuse the attribute parser later on, so we put it in a separate module.

``` {.haskell #tangle-imports}
import Text.Regex.TDFA
```

``` {.haskell file=src/Attributes.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Attributes where

import RIO
import qualified RIO.Text as T

<<attributes-imports>>
<<parse-attributes>>
```

``` {.haskell #attributes-imports}
import Document (CodeProperty(..))
import Text.Megaparsec
    ( MonadParsec, takeWhileP, chunk, endBy )
import Text.Megaparsec.Char
    ( space, letterChar )
```

The one function that we will use is the `attributes` parser.

``` {.haskell #parse-attributes}
attributes :: (MonadParsec e Text m)
           => m [CodeProperty]
attributes = (  codeClass
            <|> codeId
            <|> codeAttribute
             ) `endBy` space
```

#### Identifiers and values

Anything goes, as long as it doesn't conflict with a space separated curly braces delimited list.

``` {.haskell #parse-attributes}
cssIdentifier :: (MonadParsec e Text m)
              => m Text
cssIdentifier = do
    firstLetter <- letterChar
    rest        <- takeWhileP (Just "identifier")
                        (`notElem` (" {}=<>|" :: String))
    return $ T.singleton firstLetter <> rest

cssValue :: (MonadParsec e Text m)
         => m Text
cssValue = takeWhileP (Just "value")
                      (`notElem` (" {}=<>|" :: String))
```

#### class
A class starts with a period (`.`), and cannot contain spaces or curly braces.

``` {.haskell #parse-attributes}
codeClass :: (MonadParsec e Text m)
          => m CodeProperty
codeClass = chunk "." >> CodeClass <$> cssIdentifier
```

#### id
An id starts with a hash (`#`), and cannot contain spaces or curly braces.

``` {.haskell #parse-attributes}
codeId :: (MonadParsec e Text m)
       => m CodeProperty
codeId = chunk "#" >> CodeId <$> cssIdentifier
```

#### attribute
An generic attribute is written as key-value-pair, separated by an equals sign (`=`). Again no spaces or curly braces are allowed inside.

``` {.haskell #parse-attributes}
codeAttribute :: (MonadParsec e Text m)
              => m CodeProperty
codeAttribute = do
    key <- cssIdentifier
    _ <- chunk "="
    CodeAttribute key <$> cssValue
```

## Quasi-parsing Markdown

Parsing the markdown using MegaParsec,

### Parsing code blocks

::: {.TODO}
Add Wirth Syntax Notation with all the parsers.
:::

Code blocks are delimited with three back-quotes, and again closed with three back-quotes. Pandoc allows for any number larger than three, and using other symbols like tildes. The only form Entangled reacts to is the following **unindented** template: 

~~~markdown
 ``` {[#<reference>|.<language>|<key>=<value>] ...}
 <code> ...
 ```
~~~

Any code block is one of the following:
    
* A **referable** block: has exactly one **reference id** (`#<reference>`) and exactly one class giving the **language** of the code block (`.<language>`). Example:

~~~markdown
 ``` {.rust #hello-rust}
 println!("Hello, Rust!");
 ```
~~~

* A **file** block: has a key-value pair giving the path to the file (`file=<path>`), absolute or relative to the directory in which `entangled` runs. Again, there should be exactly one class giving the **language** of the block (`.<language>`). Example:

~~~markdown
 ``` {.rust #main-function file=src/main.rs}
 fn main() {
    <<hello-rust>>
 }
 ```
~~~

* An **ignored** block: anything not matching the previous two.

### Running configure regexes

``` {.haskell #parse-markdown}
type Match = Maybe (Text, Text, Text, [Text])

matchCodeHeader :: ConfigSyntax -> Text -> Maybe ([CodeProperty], Text)
matchCodeHeader syntax line =
    if line =~ matchCodeStart syntax
    then Just (fromMaybe [] (getLanguage' <> getFileName <> getReferenceName <> getHeaderLen), line)
    else Nothing
    where
          getLanguage' :: Maybe [CodeProperty]
          getLanguage' = do
            (_, _, _, lang) <- line =~~ extractLanguage syntax :: Match
            return (map CodeClass lang)

          getFileName :: Maybe [CodeProperty]
          getFileName = do
            (_, _, _, file) <- line =~~ extractFileName syntax :: Match
            return (map (CodeAttribute "file") file)

          getReferenceName :: Maybe [CodeProperty]
          getReferenceName = do
            (_, _, _, ref)  <- line =~~ extractReferenceName syntax :: Match
            return (map CodeId ref)

          getHeaderLen :: Maybe [CodeProperty]
          getHeaderLen = do
            (_, _, _, headLen) <- line =~~ extractProperty syntax "header" :: Match
            return $ CodeAttribute "header" <$> headLen

matchCodeFooter :: ConfigSyntax -> Text -> Maybe ((), Text)
matchCodeFooter syntax line =
    if line =~ matchCodeEnd syntax
    then Just ((), line)
    else Nothing
```

### Extracting data from a list of `CodeProperty`

In case the language is not given, or is misspelled, the system should be aware of that, so we can give an error message or warning.

``` {.haskell #parse-markdown}
getLanguage :: ( MonadReader Config m )
            => [CodeProperty] -> m ProgrammingLanguage
getLanguage [] = return NoLanguage
getLanguage (CodeClass cls : _)
    = maybe (UnknownClass cls)
            KnownLanguage
            <$> asks (\cfg -> languageName <$> lookupLanguage cfg cls)
getLanguage (_ : xs) = getLanguage xs
```

On the other hand, if an identifyer is left out, the parser should fail, so that the code-block is included as normal markdown. If only a file attribute is given, the filename is also the identifyer. If an id is also given, that takes precedence. We use a `StateT ReferenceCount` to keep track of the number of references with the same name, so that we can concatenate these code blocks in order.

``` {.haskell #parse-markdown}
getReference :: ( MonadState ReferenceCount m )
             => [CodeProperty] -> m (Maybe ReferenceId)
getReference [] = return Nothing
getReference (CodeId x:_) = Just <$> newReference (ReferenceName x)
getReference (CodeAttribute k v:xs)
    | k == "file" = do
        a <- getReference xs
        b <- Just <$> newReference (ReferenceName v)
        return $ a <|> b
    | otherwise   = getReference xs
getReference (_:xs) = getReference xs
```

We keep a separate `Map` to link certain references with top-level files.

``` {.haskell #parse-markdown}
getFilePath :: [CodeProperty] -> Maybe FilePath
getFilePath [] = Nothing
getFilePath (CodeAttribute k v:xs)
    | k == "file" = Just $ T.unpack v
    | otherwise   = getFilePath xs
getFilePath (_:xs) = getFilePath xs

getFileMap :: [ReferencePair] -> FileMap
getFileMap = M.fromList . mapMaybe filePair
    where filePair (ref, CodeBlock{..}) = do
              path <- getFilePath codeProperties
              case codeLanguage of
                  KnownLanguage l -> return (path, (referenceName ref, l))
                  _               -> Nothing
```

### References

To build up the `ReferenceMap` we define `ReferenceCount`.

``` {.haskell #parse-markdown}
data ReferenceCount = ReferenceCount
    { currentDocument :: FilePath
    , refCounts       :: Map ReferenceName Int }

countReference :: ( MonadState ReferenceCount m )
               => ReferenceName -> m Int
countReference r = do
    x <- gets (M.findWithDefault 0 r . refCounts)
    modify $ \s -> s { refCounts = M.insert r (x + 1) (refCounts s) }
    return x

newReference :: ( MonadState ReferenceCount m )
             => ReferenceName -> m ReferenceId
newReference n = do
    doc <- gets currentDocument
    x   <- countReference n
    return $ ReferenceId doc n x
```

### Parsing the document

::: {.TODO}
Add Wirth Syntax Notation with all the parsers.
:::

Using these two parsers, we can create a larger parser that works on a line-by-line basis. We define several helpers to create a parser for `ListStream Text` using single line parsers for `Text`.

``` {.haskell #tangle-imports}
-- import ListStream (parseLine, parseLineNot, tokenLine)
```

To parse markdown, we first try to parse a code block (as given above), stored in `CodeBlock`. If that fails lines are interpreted as other markdown, and stored in `PlainText`.

``` {.haskell #parse-markdown}
type DocumentParser = ReaderT Config (StateT ReferenceCount (Parsec Text (ListStream Text)))

codeBlock :: ( MonadParsec e (ListStream Text) m
             , MonadReader Config m
             , MonadState ReferenceCount m )
          => m ([Content], [ReferencePair])
codeBlock = do
    Config{..} <- ask
    -- linenum        <- Just . unPos . sourceLine . pstateSourcePos . statePosState <$> getParserState
    linenum        <- Just . (+ 2) <$> getOffset
    (props, begin) <- tokenLine (matchCodeHeader configSyntax)
    code           <- unlines'
                   <$> manyTill (anySingle <?> "code line")
                                (try $ lookAhead $ tokenLine (matchCodeFooter configSyntax))
    (_, end)       <- tokenLine (matchCodeFooter configSyntax)
    language       <- getLanguage props
    ref'           <- getReference props
    return $ case ref' of
        Nothing  -> ( [ PlainText $ unlines' [begin, code, end] ], [] )
        Just ref -> ( [ PlainText begin, Reference ref, PlainText end ]
                    , [ ( ref, CodeBlock language props code linenum ) ] )

normalText :: ( MonadParsec e (ListStream Text) m
              , MonadReader Config m )
           => m ([Content], [ReferencePair])
normalText = do
    Config{..} <- ask
    text <- unlines' <$> some (tokenLineNot $ matchCodeHeader configSyntax)
    return ( [ PlainText text ], [] )

markdown :: DocumentParser ([Content], [ReferencePair])
markdown = mconcat <$> many (codeBlock <|> normalText)

parseMarkdown' :: ( MonadReader env m, HasConfig env, MonadThrow m )
               => FilePath -> Text -> m Document
parseMarkdown' f t = do
    cfg <- view config
    let result' = parse (evalStateT (runReaderT markdown cfg) (ReferenceCount f mempty))
                        f (ListStream $ T.lines t)
    case result' of
        Left err              -> throwM $ TangleError $ tshow err
        Right (content, refs) -> return $ Document (M.fromList refs) content (getFileMap refs)
```

## Generating output files


### Indentation
Tangling is indentation sensitive. Given two code blocks

~~~markdown
 ``` {.python #multiply}
 x *= i
 ```

 ``` {.python #factorial}
 x = 1
 for i in range(n):
    <<multiply>>
 ```
~~~

We should get the code

``` {.python}
x = 1
for i in range(n):
    x *= i
```

Indentation is done by `lines` $\to$ `map (append indent)` $\to$ `unlines`, with the distinction that empty lines are not indented, and that the inverse of `lines` is `unlines'`, which doesn't append a final newline.

``` {.haskell #tangle-imports}
import TextUtil (indent, unlines')
```

``` {.haskell #generate-code}
data CodeLine = PlainCode Text
              | NowebReference ReferenceName Text
              deriving (Eq, Show)

type CodeParser = Parsec Void Text
```

We try to parse every line in a code block as a *noweb* reference. If that fails, the line is accepted as normal source text. A *noweb* reference should stand alone on a line, maybe indented with white space. Space at the end of the line is ignored.

    +--- indentation ---+--- reference  ----+--- possible space ---+
                        <<no-spaces-alowed>>

The function `parseCode` takes the `Text` from a code block and generates a list of `CodeLine`, being either `PlainCode` or `NowebReference`.

``` {.haskell #generate-code}
nowebReference :: CodeParser CodeLine
nowebReference = do
    indent' <- takeWhileP Nothing (`elem` (" \t" :: String))
    _ <- chunk "<<"
    id' <- takeWhile1P Nothing (`notElem` (" \t<>" :: String))
    _ <- chunk ">>"
    space >> eof
    return $ NowebReference (ReferenceName id') indent'

parseCode :: ReferenceName -> Text -> [CodeLine]
parseCode name = map parseLine' . T.lines
    where parseLine' l = fromRight (PlainCode l)
                       $ parse nowebReference
                              (T.unpack $ unReferenceName name)
                              l
```

## Code expansion

We don't want code blocks refering to themselves. I used to keep a history of visited `ReferenceId`s to prevent circular dependencies. This part will be moved to a validation stage.
I now use a lazy map to provide the recursion, which becomes cumbersome if we also have to detect cycles in the dependency graph.

``` {.haskell #generate-code}
type ExpandedCode m = LM.Map ReferenceName (m Text)
type Annotator m = ReferenceMap -> ReferenceId -> m Text
```

The map of expanded code blocks is generated using an induction pattern here illustrated on lists. Suppose we already have the resulting `output` and a function `g :: [Output] -> Input -> Output` that generates any element in `output` given an input and the rest of all `output`, then `f` generates `output` as follows.

``` {.haskell}
f :: [Input] -> [Output]
f input = output
    where output = map generate input
          generate = g output
```

Lazy evaluation does the rest.

``` {.haskell #generate-code}
expandedCode :: (MonadIO m, MonadReader env m, HasLogFunc env)
    => Annotator m -> ReferenceMap -> ExpandedCode m
expandedCode annotate refs = result
    where result = LM.fromSet expand (referenceNames refs)
          expand name = unlines' <$> mapM
                        (annotate refs >=> expandCodeSource result name)
                        (referencesByName refs name)

expandCodeSource :: (MonadIO m, MonadReader env m, HasLogFunc env)
    => ExpandedCode m -> ReferenceName -> Text -> m Text
expandCodeSource result name t
    = unlines' <$> mapM codeLineToText (parseCode name t)
    where codeLineToText (PlainCode x) = return x
          codeLineToText (NowebReference name' i)
              = indent i <$> fromMaybe (logWarn ("unknown reference <<" <> display name' <> ">> in #" <> display name) >> return "")
                                       (result LM.!? name')
```

We have two types of annotators:

* Naked annotator: creates the output code without annotation. From such an expansion it is not possible to untangle.

``` {.haskell #generate-code}
selectAnnotator :: (MonadError EntangledError m) => Config -> Annotator m
selectAnnotator cfg@Config{..} = case configAnnotate of
    AnnotateNaked         -> \rmap rid -> runReaderT (annotateNaked rmap rid) cfg
    AnnotateStandard      -> \rmap rid -> runReaderT (annotateComment rmap rid) cfg
    AnnotateProject       -> \rmap rid -> runReaderT (annotateProject rmap rid) cfg
```

* Commenting annotator: adds annotations in comments, from which we can locate the original code block.

We put comments in a separate module, where we also address parsing back the generated comments.

``` {.haskell #tangle-imports}
import Comment (annotateComment, annotateProject, annotateNaked)
```

## Entangled comments

``` {.haskell file=src/Comment.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Comment where

import RIO
import qualified RIO.Text as T

<<comment-imports>>

<<generate-comment>>
<<parse-comment>>
```

``` {.haskell #comment-imports}
import Control.Monad.Except

```

Tangled files start with a single line header comment, followed by nested blocks of (possibly indented) code delimeted by lines:

``` {.c}
\\ -#- entangled -#- begin <<reference-name>>[3]
code content;
\\ -#- entangled -#- end
```

The `-#- entangled -#-` bit is stored in `delim`.

``` {.haskell #generate-comment}
delim :: Text
delim = " ~\\~ "
```

Given any content, the `comment` function generates a commented line following the above prescription.

``` {.haskell #comment-imports}
import Config
    ( Config(..), ConfigLanguage(..), ConfigComment(..), languageFromName )
```

``` {.haskell #generate-comment}
getLangName :: (MonadError EntangledError m)
            => ProgrammingLanguage -> m Text
getLangName (UnknownClass cls)       = throwError $ UnknownLanguageClass cls
getLangName NoLanguage               = throwError MissingLanguageClass
getLangName (KnownLanguage langName) = return langName

comment :: (MonadReader Config m, MonadError EntangledError m)
        => ProgrammingLanguage
        -> Text
        -> m Text
comment (UnknownClass cls) _ = throwError $ UnknownLanguageClass cls
comment NoLanguage         _ = throwError MissingLanguageClass
comment (KnownLanguage langName) text = do
    cfg <- ask
    maybe (throwError $ SystemError $ "language named " <> langName <> " is not in config.")
          (\lang -> return $ formatComment lang text)
          (languageFromName cfg langName)

commentStart :: ConfigComment -> Text
commentStart (Block x _) = x
commentStart (Line x) = x

commentEnd :: ConfigComment -> Maybe Text
commentEnd (Block _ x) = Just x
commentEnd (Line _) = Nothing

formatComment :: ConfigLanguage -> Text -> Text
formatComment lang text = pre <> text <> post
    where pre  = commentStart (languageComment lang) <> delim
          post = maybe "" (" " <>) $ commentEnd $ languageComment lang
```

Using this we can write the `annotateComment` function. Given a `ReferenceId` this retrieves the code text and annotates it with a begin and end comment line.

``` {.haskell #comment-imports}
import qualified RIO.Map as M

import Document
import TextUtil (unlines')
import qualified Format
```

``` {.haskell #generate-comment}
standardPreComment :: ReferenceId -> Text
standardPreComment (ReferenceId file (ReferenceName name) count) =
    "begin <<" <> T.pack file <> "|" <> name <> ">>[" <> tshow count <> "]"

getReference :: (MonadError EntangledError m) => ReferenceMap -> ReferenceId -> m CodeBlock
getReference refs ref = maybe (throwError $ ReferenceError $ "not found: " <> tshow ref)
                              return (refs M.!? ref)

lineDirective :: (MonadReader Config m, MonadError EntangledError m)
              => ReferenceId -> CodeBlock -> m Text
lineDirective ref code = do
    Config{..} <- ask
    lang <- getLangName $ codeLanguage code
    spec <- maybe (throwError $ ConfigError $ "line directives not configured for " <> lang)
                  return (configLineDirectives M.!? lang)
    maybe (throwError $ ConfigError $ "error formatting on " <> tshow spec)
          return (Format.formatMaybe spec $ M.fromList [ ("linenumber" :: Text, tshow (fromMaybe 0 $ codeLineNumber code))
                                                       , ("filename"          , T.pack (referenceFile ref))])

annotateNaked :: (MonadReader Config m, MonadError EntangledError m)
              => ReferenceMap -> ReferenceId -> m Text
annotateNaked refs ref = do
    Config{..} <- ask
    code <- getReference refs ref
    if configUseLineDirectives then do
        line <- lineDirective ref code
        return $ unlines' [line, codeSource code]
    else return $ codeSource code

annotateComment :: (MonadReader Config m, MonadError EntangledError m)
                => ReferenceMap -> ReferenceId -> m Text
annotateComment refs ref = do
    code <- getReference refs ref
    naked <- annotateNaked refs ref
    pre <- comment (codeLanguage code) $ standardPreComment ref
    post <- comment (codeLanguage code) "end"
    return $ unlines' [pre, naked, post]

annotateProject :: (MonadReader Config m, MonadError EntangledError m)
                => ReferenceMap -> ReferenceId -> m Text
annotateProject refs ref@(ReferenceId file _ _) = do
    code <- getReference refs ref
    naked <- annotateNaked refs ref
    let line = fromMaybe 0 (codeLineNumber code)
    pre  <- comment (codeLanguage code) (standardPreComment ref <> " project://" <> T.pack file <> "#" <> tshow line)
    post <- comment (codeLanguage code) "end"
    return $ unlines' [pre, naked, post]

headerComment :: ConfigLanguage -> FilePath -> Text
headerComment lang path = formatComment lang
    $ "language=" <> languageName lang <> " filename=" <> T.pack path
```

### Parsing comments

``` {.haskell #comment-imports}
import Text.Megaparsec            ( MonadParsec, chunk, skipManyTill
                                  , anySingle, (<?>), takeWhileP, takeRest )
import Text.Megaparsec.Char.Lexer ( decimal )

import Attributes (attributes, cssIdentifier, cssValue)
```

The same comment lines have to be parsed back when we untangle. The first line is the top header comment. We don't know yet what the language is going to be, so we `skipManyTill` we find the `delim` text.

``` {.haskell #parse-comment}
topHeader :: ( MonadParsec e Text m )
          => m [CodeProperty]
topHeader = do
    _    <- skipManyTill (anySingle <?> "open comment")
                         (chunk delim)
    attr <- attributes
    _    <- takeRest
    return attr
```

Other parsers will always be combined with `commented`, giving the value of the original parser and a `Text` that gives the indentation level of the parsed comment.

``` {.haskell #parse-comment}
commented :: (MonadParsec e Text m)
          => ConfigLanguage -> m a -> m (a, Text)
commented lang p = do
    indent <- takeWhileP (Just "initial indent") (`elem` (" \t" :: String))
    _ <- chunk $ commentStart (languageComment lang) <> delim
    x <- p
    -- _ <- chunk (fromMaybe "" $ commentEnd $ languageComment lang)
    -- space
    -- eof
    _ <- takeRest
    return (x, indent)
```

``` {.haskell #parse-comment}
beginBlock :: (MonadParsec e Text m)
           => m ReferenceId
beginBlock = do
    _ <- chunk "begin <<"
    doc  <- cssValue
    _ <- chunk "|"
    name <- cssIdentifier
    _ <- chunk ">>["
    count <- decimal
    _ <- chunk "]"
    _ <- takeRest
    return $ ReferenceId (T.unpack doc) (ReferenceName name) count

endBlock :: (MonadParsec e Text m)
         => m ()
endBlock = void $ chunk "end"
```

# Parsing anything with Megaparsec

Megaparsec is considered to be preferable over older Parsec. It is a bit more work to parse different things than strings though. We will be parsing lists of items.

``` {.haskell file=src/ListStream.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module ListStream where

import RIO
import RIO.List (splitAt, headMaybe)
import Text.Megaparsec ( Parsec, MonadParsec, token, parse, satisfy
                       , Stream (..), VisualStream (..), TraversableStream(..)
                       , PosState (..), SourcePos (..), mkPos, unPos )

<<instance-list-stream>>

<<list-stream-helpers>>
<<line-parser>>
```

The minimal definition of a stream must define `tokensToChunk`, `chunkToTokens`, `chunkLength`, `take1_`, `takeN_`, `takeWhile_`. We should alse like to implement `VisualStream`, for which we need `showTokens`, and  `TraversableStream` which needs `reachOffset`.

Because our instance will overlap with that of `[Char] = String` we'll have to wrap the `ListStream` in a `newtype`.

``` {.haskell #instance-list-stream}
newtype ListStream a = ListStream { unListStream :: [a] }
    deriving (Show, Foldable, Semigroup, Monoid, Functor, Applicative, Monad)

instance (Eq a, Ord a, Show a) => Stream (ListStream a) where
    <<list-stream-associated-types>>
    <<list-stream-methods>>

instance (Eq a, Ord a, Show a) => VisualStream (ListStream a) where
    <<list-visual-stream-methods>>

instance (Eq a, Ord a, Show a) => TraversableStream (ListStream a) where
    <<list-traversable-stream-methods>>
```

``` {.haskell #list-stream-associated-types}
type Token (ListStream a) = a
type Tokens (ListStream a) = [a]
```

Converting a token to a chunk is just `pure`.

``` {.haskell #list-stream-methods}
tokenToChunk Proxy = pure
```

Tokens to chunk or even chunk to tokens are both the identity function

``` {.haskell #list-stream-methods}
tokensToChunk Proxy = id
chunkToTokens Proxy = id
chunkLength Proxy = length
chunkEmpty Proxy = null
```

In fact, these definitions are identical to those in the instance definition for `Stream String` in Megaparsec.

``` {.haskell #list-stream-methods}
take1_ (ListStream [])     = Nothing
take1_ (ListStream (x:xs)) = Just (x, ListStream xs)
takeN_ n (ListStream xs)
    | n <= 0    = Just ([], ListStream xs)
    | null xs   = Nothing
    | otherwise = Just (h, ListStream t) where (h, t) = splitAt n xs
takeWhile_ p (ListStream xs) = (t, ListStream h) where (h, t) = break p xs
```

``` {.haskell #list-visual-stream-methods}
showTokens Proxy = show
```

Here's where things get different. Offsets are lines in a source file. Therefore, `sourceColumn` is always 1. The precise working of `reachOffset` is a bit guesswork.

``` {.haskell #list-stream-helpers}
offsetSourcePos :: Int -> SourcePos -> SourcePos
offsetSourcePos offset sourcePos@SourcePos{..}
    = sourcePos { sourceLine = mkPos (unPos sourceLine + offset) }
```

``` {.haskell #list-traversable-stream-methods}
reachOffset offset state@PosState{..} = (Just repr, state')
    where sourcePos = offsetSourcePos (offset - pstateOffset) pstateSourcePos
          input     = ListStream $ drop (offset - pstateOffset) (unListStream pstateInput)
          repr      = maybe "<end of stream>" show $ headMaybe $ unListStream input
          state'    = state
                    { pstateInput = input
                    , pstateOffset = offset
                    , pstateSourcePos = sourcePos }
```

## Text-parser tokens

``` {.haskell #line-parser}
type LineParser = Parsec Void Text

parseLine :: LineParser a -> Text -> Maybe (a, Text)
parseLine p t = either (const Nothing) (\x -> Just (x, t))
              $ parse p "" t

parseLineNot :: LineParser a -> Text -> Maybe Text
parseLineNot p t = either (const $ Just t) (const Nothing)
                 $ parse p "" t

tokenLine :: ( MonadParsec e (ListStream Text) m )
          => (Text -> Maybe a) -> m a
tokenLine f = token f mempty

tokenLineNot :: ( MonadParsec e (ListStream Text) m )
             => (Text -> Maybe a) -> m Text
tokenLineNot f = satisfy (isNothing . f)

tokenP :: ( MonadParsec e (ListStream Text) m )
       => LineParser a -> m (a, Text)
tokenP = tokenLine . parseLine
```


## Testing on lists

``` {.haskell file=test/ListStreamSpec.hs}
module ListStreamSpec where

import Test.Hspec
import Test.QuickCheck

import Text.Megaparsec
import Data.Void

import ListStream

<<list-stream-props>>

listStreamSpec :: Spec
listStreamSpec =
    <<list-stream-spec>>
```

### Lists of integers

``` {.haskell #list-stream-props}
prop_parser :: (Eq a) => Parsec e s a -> s -> a -> Bool
prop_parser p input expected = success (parse p "" input)
    where success (Left _) = False
          success (Right x) = x == expected

parseAny :: (Show a, Ord a, Eq a) => Parsec Void (ListStream a) [a]
parseAny = takeWhileP Nothing (const True)

prop_tokens :: [Int] -> Bool
prop_tokens xs = prop_parser parseAny (ListStream xs) xs
```

``` {.haskell #list-stream-spec}
describe "Parsing integers" $
    it "takeWhileP yield input" $
        property prop_tokens
```
# All you'll ever need is Markdown
Markdown is widely used for writing up documents. It is supported by many blog engines, GitHub readme, you name it. There are many varieties, dialects and extensions to Markdown. The standard is described on the Daring Fireball website. One place where many flavours of Markdown meet is in Pandoc. Pandoc can convert Markdown to many different formats like HTML, LaTeX, PDF (through LaTeX), RTF, DOC, EPUB and even back to Markdown of a different flavour. There are some little known features (or extensions) of Markdown that make it very versatile and suitable for any rich text content, especially if you use Pandoc.

The take-away message is: you will never again have to write a document in LaTeX or HTML directly, not for writing notes, reports or papers, no power-point for presentations, no awkward WYSIWYG editor for web content (I'm talking to you Medium!); just use Markdown.
You might ask "But Johan, how do I then ...?", shush! The answer is always going to be: **Pandoc**.

## Primer
Markdown is a way of writing up rich content (i.e. text, document hierarchy, images, lists, quotes, etc.) in plain text using a human readable format. The format is aimed to mimic in plain text the formatted result:

~~~markdown
Level 1 Header
==============
Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do 
eiusmod tempor incididunt ut labore et dolore magna aliqua.

Level 2 Header
--------------
Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris 
nisi ut aliquip ex ea commodo consequat.

* list item
* Another item

## Also a Level 2 Header
Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
> quote: Excepteur sint occaecat cupidatat non proident, sunt in 
> culpa qui officia deserunt mollit anim id est laborum
~~~

You can read Daring Fireball for more syntax. The basics of Markdown are highly intuitive, but in many instances standard Markdown does not suffice: there is no support for equations, citation management, cross references or numbered figures. However, Pandoc supports many extensions and flavours of Markdown that have evolved over the years. Markdown feels a bit like a natural language. This has lead to some [critique on the use of Markdown in technical documents](https://www.ericholscher.com/blog/2016/mar/15/dont-use-markdown-for-technical-docs/).

I will highlight some extensions, all clearly documented and supported by Pandoc, that transform Markdown into a versatile and extendable format, suitable for even the most technical and demanding documents, yet easy to use when no such demands arise.

## Flavours and extensions
There are many extensions to Markdown, some of which are more widely used than others. One of the most influential flavour of Markdown is the superset by Github, its major addition being that of delimited code blocks.

~~~markdown
```javascript
var factorial = n => (n < 2) ? 1 : n * factorial(n - 1);
```
~~~

### Attributes
Nothing good ever comes from PHP, except for PHP Markdown, which adds attributes to Markdown. Any element in a Markdown document can be adorned with (CSS style) attributes by appending them in curly braces: `{#id .class .class key=value}`. For example,

~~~markdown
## The proof {#proof .section color=red}
~~~

when converted to HTML looks like

~~~html
<h2 id="proof" class="section" color="red">The proof</h2>
~~~

or to LaTeX

~~~latex
\hypertarget{proof}{%
\subsection{The proof}\label{proof}}
~~~

In the case of LaTeX, the color attribute as well as the class is ignored, because Pandoc (by default) doesn't know what to do with it. The #proof id attached to this header can be used further on in the document to make cross-references.

The attribute syntax also applies to code blocks: a line starting with `\`\`\`python` is equivalent to `\`\`\` {.python}`.

### DIV elements
A `<div>` element can be added using three (or more) colons.

~~~markdown
::: {.warning}
This is a warning!
:::
~~~

This can be used to write down any non-standard element. So how is this rendered? If your output format is HTML, change the style sheet. For generating LaTeX however we need to do a bit more work. This is where Pandoc comes in. Pandoc has support for filtering elements and creating relevant output, but we'll get back to that.

### YAML metadata
Information that is usually contained in a HTML `<header>` region can be included in the YAML metadata block. This is a block delimited by hyphens at the top of the document.

~~~markdown
---
title: A theorem on right angled triangles
author: Pythagoras of Samos
date: ~520 BC, a hot summer night
---
~~~

This is also the place where you can configure options for Pandoc and its filters.

### Citations
Citations are managed using the pandoc-citeproc plugin. For those who have worked with BibTeX (or another bibliography database) before, you include the bibliography by adding a line to the YAML header block. Say you have stored your references in ref.bib , and want to create a dedicated section called "References" to list the citations:

~~~markdown
---
bibliography: ref.bib
reference-section-title: References
---
~~~

There are many ways to cite papers, the syntax for which is documented in the Pandoc manual. Basic citation looks like `[@Hidding2014]`, where `Hidding2014` is an entry in the bibliography. Several editors (VSCode and Emacs) even support autocompletion on the included bibliography.

### Equations
Equations can be entered using the famous LaTeX DSL for equations. Use single $ characters to delimit an inline equation and double $$ to delimit a full width equation.

~~~markdown
Given a right angled triangle of sides $a$, $b$ and hypotenuse $c$, we may state that,
$$a^2 + b^2 = c^2.$$
~~~

When converting to a LaTeX based output this will work trivially. For HTML it is probably best to use MathJax, enabled in Pandoc with the `--mathjax` and `--standalone` flags.

I promised to show you how Markdown can be extended with arbitrary functionality using Pandoc. It is inevitable that I will get a bit more technical here, but the rewards are high.

## Pandoc basics
Pandoc reads and writes documents of many formats. It does so by converting to and from a native intermediate representation. We'll create a document named last-theorem.md and enter the following:

~~~markdown
---
title: Last theorem
author: P. de Fermat
---

For integer values of $a$, $b$, $c$ and $n > 2$, the equation
$$a^n + b^n = c^n,$$
has no solution. I have a truly remarkable proof for this.
~~~

We can see how Pandoc reads a document by running
```
pandoc -f markdown+yaml_metadata_block -t native -s last-theorem.md
```
which will give (slightly edited for readability)

```haskell
Pandoc
 (Meta {unMeta = fromList
   [("author", MetaInlines [Str "P.",Space,Str "de",Space,Str "Fermat"])
   ,("title", MetaInlines [Str "Last",Space,Str "theorem"])]})
 [Para [Str "I",Space,Str "state",Space,Str "that,",Space
       ,Str "for",Space,Str "integer",Space,Str "values",Space
       ,Str "of",Space,Math InlineMath "a",Str ",",Space
       ,Math InlineMath "b",Str ",",Space
       ,Math InlineMath "c",Space,Str "and",Space
       ,Math InlineMath "n > 2"
       ,Str ",",Space,Str "the",Space,Str "equation"]
 ,Para [Math DisplayMath "a^n + b^n = c^n,"]
 ,Para [Str "has",Space,Str "no",Space,Str "solution.",Space
       ,Str "I",Space,Str "have",Space,Str "a",Space
       ,Str "truly",Space,Str "elegant",Space,Str "proof"
       ,Space,Str "for",Space,Str "this,",Space,Str "but",Space
       ,Str "this",Space,Str "is",Space,Str "not",Space
       ,Str "the",Space,Str "place",Space,Str "to",Space
       ,Str "give",Space,Str "it."]]
```

Yes, this is Haskell syntax. It just means what you think it means. Using native output in Pandoc will be very useful if you start developing your own filters. For now it just serves to illustrates how Pandoc works. If you're not a hacker, you'll never have to look at this again 🤓👍.

Let's create a PDF from this mathematics wizardry.

```
pandoc -f markdown+yaml_metadata_block -t latex -o last-theorem.pdf -s last-theorem.md
```

Resulting in a nice PDF rendering:

![Typeset paper on a famous theorem](last-theorem-1.png)

### Pandoc filters
There is a big problem with the above example. The equation is not numbered! Pandoc filters let you change the intermediate representation. Let's try the `pandoc-eqnos` filter. This filter is written in Python (using `pypandoc`); any language that can read the intermediate JSON representation works.

```
pip install --user pandoc-eqnos
```

We'll need to change the last-theorem.md document a bit. To get an equation numbered, add an `id` to the equation.

~~~markdown
$$a^n + b^n = c^n,$$ {#eq:fermat}
~~~

Also add a sentence to the end.

~~~markdown
The prove of inexistence of a solution for Equation @eq:fermat had eluded mathematicians for centuries.
~~~

Now to invoke the equation numbering filter, add the `--filter pandoc-eqnos` argument to the command line.

```
pandoc -f markdown+yaml_metadata_block -t latex --filter pandoc-eqnos -o last-theorem.pdf -s last-theorem.md
```

![Typeset paper on a famous theorem with numbered equations](last-theorem-2.png)

Pandoc command lines can grow out of hand rather quickly. It is advisable to manage your Pandoc settings in a Bash script or Make file, whatever you prefer.

### Lua filters
Pandoc has built-in support for filters written in Lua. Filters written in Lua are generally faster than those written in Python or other external languages. Lua filters forego documents being passed to an external program via JSON, but rather work directly on the abstract syntax tree as it is represented in Pandoc itself.

Let's add a feature. Add the following to our budding math paper:

~~~markdown
::: {.warning}
Modular forms and elliptic curves ahead!
:::
~~~

To parse this, Pandoc needs the extension `fenced_divs` enabled.

```
pandoc -f markdown+yaml_metadata_block+fenced_divs --filter pandoc-eqnos -t native -s last-theorem.md
```

At the end of the output will be the expression:

```haskell
Div 
 ("",["warning"],[])
 [Para [Str "Modular",Space,Str "forms",Space,Str "and",Space
 ,Str "elliptic",Space,Str "curves",Space,Str "ahead!"]]
```

Once we generate HTML from this example, we can add the proper CSS to the .warning class to change the looks of the warning. However, in the generated PDF there's no change. Let's make a filter that creates a nice coloured box in the LaTeX output. We define a filter that runs on all Div elements in the file `warning-div.lua`.

```lua
function Div(el)
  if el.classes[1] == "warning" then
    -- insert element in front
    table.insert(
      el.content, 1,
      pandoc.RawBlock("latex", "\\begin{warning}"))
    -- insert element at the back
    table.insert(
      el.content,
      pandoc.RawBlock("latex", "\\end{warning}"))
  end
  return el
end
```

The filter checks if the div has class `warning`, if so, it adds a LaTeX `RawBlock` at the start and end of the `div`.

Running Pandoc with `--lua-filter=warning-div.lua` now converts the div element to a LaTeX string

```latex
\begin{warning}
Modular forms and elliptic curves ahead!
\end{warning}
```

This is not standard LaTeX, so we'll need to define a macro in `warning.tex`:

```latex
\usepackage{tcolorbox}
\newenvironment{warning}
    {\begin{tcolorbox}[colbacktitle=red!50!white,
                       title=Warning,coltitle=black,
                       fonttitle=\bfseries]}
    {\end{tcolorbox}}
```

We can run Pandoc again to generate the PDF. The `-H` option can be used to include files into the generated output.

```
pandoc -f markdown+yaml_metadata_block+fenced_divs --filter pandoc-eqnos --lua-filter warning-div.lua -H warning.tex -t latex -o last-theorem.pdf -s last-theorem.md
```

Did I mention Pandoc command lines tend to grow out of hand?

![Typeset paper on a famous theorem with numbered equations and a warning](last-theorem-3.png)

## Skies and limits
Granted, to unlock the full power of Markdown for the web, you'll need to know some HTML and CSS, and to tweak PDF output to your impossibly high standards you need to grok LaTeX. All the more reason to create an ecosystem of scripts, themes and tutorials to ease the learning curve. Also, code editors should offer better support for more varieties of markdown. I don't mean cluttering up the editing experience with more distracting tool-tips, snippets etc. Just this: correct and efficient highlighting and proper outline support.

This was just a teaser of what is possible with Pandoc. Did I mention creating slide shows with `reveal.js`? Or how about doing some literate programming with `entangled`? The documentation of Pandoc is excellent, so just go ahead and write all your content in Markdown!
% ENTANGLED(1)
% Johan Hidding (j.hidding@esciencecenter.nl)
% 2019

# NAME

entangled - literate programming with markdown

# SYNOPSIS

`entangled` [*input files*] [*`--help`*] [*`--defaults`*] --- **Watch input files for changes**  
`entangled` *`insert`* [*input files*]                        --- **Insert input files into database**  

# DESCRIPTION

**entangled** is a tool to aid literate programming with Markdown. In literate programming, you do all programming in code blocks inside a natural language literary work. These code blocks are extracted to generate the source files that are suitable for the compiler, a process known as **tangling**. Using **entangled**, you can also do the inverse, patching changes to the generated files back to the markdown documents, also known as **stitching**. **entangled** runs as a daemon in the background, and is triggered on write events to any of the files involved.

Configuration of **entangled** is managed from a *`entangled.toml`* file, which is searched for in all parent folders and *`$HOME/.config/entangled`*.

# GENERAL OPTIONS

*`-h`*, *`--help`*
: Display a friendly help message.

*`--defaults`*
: Print an `entangled.toml` file with default settings.
---
title: Entangled, literate programming Swiss army knife
author: Johan Hidding
---

Entangled makes writing literate programs easier by keeping code blocks in markdown up-to-date with generated source files. By monitoring the tangled source files, any change in the master document or source files is reflected in the other. In practice this means:

* Write well documented code using Markdown.
* Use any programming language you like (or are forced to use).
* Keep debugging and using other IDE features without change.
* Generate a report in PDF or HTML from the same source (see examples on the right).

# Preliminaries

## Modules

Entangled uses `RIO`. This saves a lot of standard import statements and comes with all the perks of the Reader-IO pattern and `UnliftIO`. In one instance we will be using lazy maps though.

``` {.haskell #import-lazy-map}
import qualified Data.Map.Lazy as LM
```

We will be using the `Text` module everywhere. All parsing will be done through `Megaparsec`. We use `FSNotify` to watch events on the filesystem.

## Errors

``` {.haskell file=src/Errors.hs}
{-# LANGUAGE NoImplicitPrelude #-}
module Errors where

import RIO

data EntangledError
    = TangleError Text
    | StitchError Text
    | ReferenceError Text
    | CyclicReference Text
    | UnknownLanguageClass Text
    | DatabaseError Text
    | SystemError Text
    | MissingLanguageClass
    | NotYetImplemented Text
    | ConfigError Text
    | UnknownError
    deriving (Show, Ord, Eq, Typeable)

toEntangledError :: (Show e)
                 => (Text -> EntangledError) -> Either e a
                 -> Either EntangledError a
toEntangledError _ (Right x) = Right x
toEntangledError f (Left x) = Left $ f $ tshow x

instance Exception EntangledError

formatError :: EntangledError -> Text
formatError (TangleError t) = "tangling: " <> t
formatError (StitchError t) = "stitching: " <> t
formatError x = tshow x
```

# Design

Entangled is a command-line tool, structured around a SQLite database. The daemon ties up the different command-line facilities in a loop managed by `FSNotify`.

``` {.python #root}
Hello
<<doesnt-exists>>
World
```
# Target files with shared reference

``` {.scheme #a}
(display "a")
```

``` {.scheme file=a.scm}
<<a>>
```

``` {.scheme file=b.scm}
<<a>>
```

# Simple test

This tests basic one-shot behaviour of entangled.

``` {.scheme file=hello.scm}
<<hello-world>>
```

``` {.scheme #hello-world}
(display "Hello, World!") (newline)
```

## Factorials

``` {.scheme file=factorial.scm}
<<factorial>>

(display (factorial 10)) (newline)
```

``` {.scheme #factorial}
(define (factorial n)
  (let loop ((x 1)
             (n n))
    (if (zero? n)
      x
      (loop (* x n) (- n 1)))))
```

# Testing multi-layer

``` {.scheme file=case3.scm}
<<a>>
<<b>>
(display "c")
```

``` {.scheme #b}
(display "B")
```

Test if files and directories get removed when they are orphaned.

``` {.scheme file=sits/in/nested/directories.scm}
(let* ((yin
         ((lambda (cc) (display #\@) cc) (call/cc (lambda (c) c))))
       (yang
         ((lambda (cc) (display #\*) cc) (call/cc (lambda (c) c)))))
    (yin yang))
```

``` {.scheme file=sits/unchanged.scm}
(define Y
  (lambda (h)
    ((lambda (x) (h (lambda (a) ((x x) a))))
     (lambda (x) (h (lambda (a) ((x x) a)))))))
```

# Testing multiple outputs

``` {.scheme #a}
(display "a")
```

``` {.scheme #b}
(display "b")
<<a>>
```

## Case 1

``` {.scheme file=case1.scm}
<<a>>
<<b>>
```

## Case 2

``` {.scheme file=case2.scm}
<<a>>
<<a>>
```

# Editing one of multiple instances

``` {.scheme #a}
(display "a")
```

``` {.scheme file=a1.scm}
<<a>>
<<a>>
(display "1")
```

``` {.scheme file=a2.scm}
<<a>>
(display "2")
```

``` {.scheme file=b.scm}
(display "b")
```

# Command-line tests
These tests run in Bash. Each starts with a set of Markdown files that are sequentially patched to see if Entangled give the correct output. The different testing suites should have a `.test` extension. For each suite the contents of this folder (the `.test`, `.dhall`, `.md` and `.diff` files that is) are copied to a temporary folder, where the commands in the `.test` file are executed, after which the temporary folder is deleted. Run with:

    ./run-tests.sh

If things take a long time, this means that `cabal` is recompiling the code. If you are developing a test, you may want to run the tests locally to do a post-mortem. In that case, run with `-u` to select the test you're working on, `-d` to run locally and `-x` to quit after the first failure:

    ./run-tests.sh -xdu 01-basic

## Rationale
Entangled is a simple program, but its workings are sometimes subtle and very stateful. In stead of loading a debugger and checking the state of the memory when something goes wrong, we work with a Sqlite database. All the actions that the daemon does can be programmed from the command-line, allowing for easy inspection.

What happens if we try to declare a file twice?

``` {.scheme file=hello.scm}
(display "Hello, ")
```

``` {.scheme file=hello.scm}
(display "World!") (newline)
```

# Simple test with block comments

This tests basic one-shot behaviour of entangled.

``` {.c file=hello.c}
#include <stdio.h>
#include <stdlib.h>

int main() {
    <<hello-world>>
    return EXIT_SUCCESS;
}
```

``` {.c #hello-world}
printf("Hello, World!\n");
```

## Factorials

``` {.c file=factorial.c}
#include <stdio.h>
#include <stdlib.h>

<<factorial>>

int main() {
    printf("%u! = %u\n", 10, factorial(10));
    return EXIT_SUCCESS;
}
```

``` {.c #factorial}
unsigned factorial(unsigned n) {
    unsigned m = 1, i;
    for (i = 1; i <= n; ++i) {
        m *= i;
    }
    return m;
}
```

# This should give line directives

``` {.c++ file=test.cc}
<<includes>>
<<main>>
```

``` {.c++ #includes}
#include <iostream>
```

``` {.c++ #main}
int main() {
    <<hello>>
    return 0;
}
```

``` {.c++ #hello}
std::cout << "Hello, World!" << std::endl;
```

Try different syntax

```c file="main.c"
#include <stdio.h>
#include <stdlib.h>

int main() {
    printf("Hello, World!\n");
    return EXIT_SUCCESS;
}
```
---
title: test1
---

```
This will not be tangled.
```

We present our first program:

``` {.cpp file=hello.cc}
#include <iostream>
#include <cstdlib>

int main() {
    std::cout << "Hello, World!" << std::endl;
    return EXIT_SUCCESS;
}
```

Will print `Hello, World!` on the console.

Empty code blocks should not be a problem:

``` {.scm file=lala.scm}
```
---
title: test3
---

We present our first program:

``` {.cpp #imports}
#include <iostream>
#include <cstdlib>
```

``` {.cpp file=hello.cc}
<<imports>>

int main() {
    <<main-body>>
}
```

Will print `Hello, World!` on the console.

``` {.cpp #main-body}
std::cout << "Hello, World!" << std::endl;
```

And make sure we exit properly:

``` {.cpp #main-body}
return EXIT_SUCCESS;
```
---
title: test4
---

``` {.hs file=hello.hs}
main :: IO ()
main = putStrLn "Hello World!"
```

``` {.scm file=hello.scm}
(display "Hello World!")
(newline)
```
---
title: test2
---

We present our first program:

``` {.cpp file=hello.cc}
#include <iostream>
#include <cstdlib>

int main() {
    <<main-body>>
}
```

Will print `Hello, World!` on the console.

``` {.cpp #main-body}
std::cout << "Hello, World!" << std::endl;
return EXIT_SUCCESS;
```
---
title: Cyclic reference
---

Tangling this file should fail with a `CyclicReference` error.

``` {.cpp #test}
std::cout << "When time becomes a loop ...\n";
<<test>>
```

``` {.cpp file=loop.cc}
#include <iostream>
#include <cstdlib>

int main() {
    <<test>>
    return EXIT_SUCCESS;
}
```

``` {.python file=a.py}
pass
```

``` {.python file=a.py}
pass
```


``` {.python #a}
<<b>>
```

``` {.python #b}
<<a>>
```

