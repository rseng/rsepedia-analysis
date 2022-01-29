# ChangeLog

Version numbers follow [semantic versioning](https://semver.org)


## phonemizer-3.0.1

* **bug fixes**

  * The method `BaseBackend.phonemize` now raises a `RuntimeError` if the input
    text is a str instead of a list of of str (was only logging an error
    message).

  * Preserve punctuation alignement when using `--preserve-punctuation`, was
    inserting a space before each punctuation token, see [issue
    #97](https://github.com/bootphon/phonemizer/issues/97).


## phonemizer-3.0

* **breaking change**

  * Do not remove empty lines from output. For example:

  ```python
  # this is now
  phonemize(["hello", "!??"]) == ['həloʊ ', '']
  # this was
  phonemize(["hello", "!??"]) == ['həloʊ ']
  ```

  * Default backend in the `phonemize` function is now `espeak` (was
    `festival`).

  * `espeak-mbrola` backend now requires `espeak>=1.49`.

  * `--espeak-path` option renamed as `--espeak-library`and
    `PHONEMIZER_ESPEAK_PATH` environment variable renamed as
    `PHONEMIZER_ESPEAK_LIBRARY`.

  * `--festival-path` option renamed as `--festival-executable` and
    `PHONEMIZER_FESTIVAL_PATH` environment variable renamed as
    `PHONEMIZER_FESTIVAL_EXECUTABLE`.

  * The methods `backend.phonemize()` from the backend classes take only a list
    of str a input text (was either a str or a list of str).

  * The methods `backend.version()` from the backend classes returns a tuple of
    int instead of a str.

* **improvements**

  * `espeak` and `mbrola` backends now rely on the `espeak` shared library using
    the `ctypes` Python module, instead of reliying on the `espeak` executable
    through subprocesses. This implies **drastic speed improvments, up to 40 times
    faster**.

* **new features**

  * New option `--prepend-text` to prepend the input text to phonemized
    utterances, so as to have both orthographic and phonemized available at
    output.

  * New option `--tie` for the `espeak` backend to display a tie character
    within multi-letter phonemes. (see issue
    [#74](https://github.com/bootphon/phonemizer/issues/74)).

  * New option `--words-mismatch` for the `espeak` backend. This allows to
    detect when espeak merge consecutive words or drop a word from the
    orthographic text. Possible actions are to ignore those misatches, to issue
    a warning for each line where a mismatch is detectd, or to remove those
    lines from the output.

* **bugfixes**

  * phonemizer's logger no more conflicts with other loggers when imported from
    Python (see PR [#61](https://github.com/bootphon/phonemizer/pull/61)).


## phonemizer-2.2.2

* **bugfixes**

  * fixed installation from source (bug introduced in 2.2.1, see
    issue [#52](https://github.com/bootphon/phonemizer/issues/52)).

  * Fixed a bug when trying to restore punctuation on an empty text (see issue
    [#54](https://github.com/bootphon/phonemizer/issues/54)).

  * Fixed an edge case bug when using custom punctuation marks (see issue
    [#55](https://github.com/bootphon/phonemizer/issues/55)).

  * Fixed regex issue that causes digits to be considered punctuation (see
    issue [#60](https://github.com/bootphon/phonemizer/pull/60)).


## phonemizer-2.2.1

* **improvements**

  From Python import the phonemize function using `from phonemizer import
  phonemize` instead of `from phonemizer.phonemize import phonemize`. The
  second import is still available for compatibility.

* **bugfixes**

  * Fixed a minor bug in `utils.chunks`.

  * Fixed warnings on language switching for espeak backend when using parallel
    jobs (see issue [#50](https://github.com/bootphon/phonemizer/issues/50)).

  * Save file in utf-8 explicitly for Windows compat (see issue
    [#43](https://github.com/bootphon/phonemizer/issues/43)).

  * Fixed build and tests in Dockerfile (see issue
    [#45](https://github.com/bootphon/phonemizer/issues/45)).


## phonemizer-2.2

* **new features**

  * New option ``--list-languages`` to list the available languages for a given
    backend from the command line.

  * The ``--sampa`` option of the ``espeak`` backend has been replaced by a new
    backend ``espeak-mbrola``.

    * The former ``--sampa`` option (introduced in phonemizer-2.0) outputs
      phones that are not standard SAMPA but are adapted to the espeak TTS
      front-end.

    * On the other hand the ``espeak-mbrola`` backend allows espeak to output
      phones in standard SAMPA (adapted to the mbrola TTS front-end). This
      backend requires mbrola to be installed, as well as additional mbrola
      voices to support needed languages. **This backend does not support word
      separation nor punctuation preservation**.

* **bugfixes**

  * Fixed issues with punctuation processing on some corner cases, see issues
    [#39](https://github.com/bootphon/phonemizer/issues/39) and
    [#40](https://github.com/bootphon/phonemizer/issues/40).

  * Improvments and updates in the documentation (Readme, ``phonemize --help``
    and Python code).

  * Fixed a test when using ``espeak>=1.50``.

  * Empty lines are correctly ignored when reading text from a file.


## phonemizer-2.1

* **new features**

  * Possibility to preserve the punctuation (ignored and silently removed by
    default) in the phonemized output with the new option
    ``--preserve-punctuation`` from command line (or the equivalent
    ``preserve-punctuation`` from Python API). With the ``punctuation-marks``
    option, one can overload the default marls considered as punctuation.

  * It is now possible to specify the path to a custom ``espeak`` or
    ``festival`` executable (for instance to use a local installation or to test
    different versions). Either specify the ``PHONEMIZER_ESPEAK_PATH``
    environment variable, the ``--espeak-path`` option from command line or use
    the ``EspeakBackend.set_espeak_path`` method from the Python API. Similarly
    for festival use ``PHONEMIZER_FESTIVAL_PATH``, ``--festival-path`` or
    ``FestivalBackend.set_festival_path``.

  * The ``--sampa`` option is now available for espeak (was available only for
    espeak-ng).

  * When using ``espeak`` with SAMPA output, some SAMPA phones are corrected to
    correspond to the normalized SAMPA alphabet (espeak seems not to respect
    it). The corrections are language specific. A correction file must be placed
    in ``phonemizer/share/espeak``. This have been implemented only for French
    by now.

* **bugfixes**

  * parses correctly the version of ``espeak-ng`` even for dev versions (e.g.
    ``1.51-dev``).

  * fixed an issue with ``espeak`` backend, where multiple phone separators can be
    present at the end of a word, see
    [#31](https://github.com/bootphon/phonemizer/issues/31).

  * added an additional stress symbol ``-`` for ``espeak``.


## phonemizer-2.0.1

* **bugfixes**

  * ``keep-flags`` was not the default argument for ``language_switch`` in the
    class ``EspeakBackend``.

  * fixed an issue with punctuation processing in the espeak backend, see
    [#26](https://github.com/bootphon/phonemizer/issues/26)

* **improvements**

  * log a warning if using ``python2``.


## phonemizer-2.0

* **incompatible change**

  Starting with ``phonemizer-2.0`` only python3 is supported. **Compatibility
  with python2 is no more ensured nor tested.** https://pythonclock.org.

* **bugfixes**

  * new ``--language-switch`` option to use with ``espeak`` backend to deals
    with language switching on phonemized output. In previous version there was
    a bug in detection of the language switching flags (sometimes removed,
    sometimes not). Now you can choose to keep the flags, to remove them, or to
    delete the whole utterance.

  * bugfix in a test with `espeak>=1.49.3`.

  * bugfix using `NamedTemporaryFile` on windows, see
    [#21](https://github.com/bootphon/phonemizer/issues/21).

  * bugfix when calling *festival* or *espeak* subprocesses on Windows, see
    [#17](https://github.com/bootphon/phonemizer/issues/17).

  * bugfix in detecting recent versions of *espeak-ng*, see
    [#18](https://github.com/bootphon/phonemizer/issues/18).

  * bugfix when using utf8 input on *espeak* backend (python2), see
    [#19](https://github.com/bootphon/phonemizer/issues/19).


* **new features and improvements**

  * new `--sampa` option to output phonemes in SAMPA alphabet instead of IPA,
    available for espeak-ng only.

  * new ``--with-stress`` option to use with ``espeak`` backend to not remove the
    stresses on phonemized output. For instance:

        $ echo "hello world" | phonemize
        həloʊ wɜːld
        $ echo "hello world" | phonemize --with-stress
        həlˈoʊ wˈɜːld

  * improved logging: by default only warnings are displayed, use the new
    ``--quiet`` option to inhibate all log messages or ``--verbose`` to see all of
    them. Log messages now display level name (debug/info/warning).

  * improved code organization:

    * backends are now implemented in the ``backend`` submodule
      as separated source files.

    * improved version string (displays uninstalled backends, moved outside of
      main for use from Python).

    * improved logger implemented in its own module so as a call to phonemizer
      from CLI or API yields the same log messages.


## phonemizer-1.0

* **incompabile changes**

  The following changes break the compatibility with previous versions
  of phonemizer (0.X.Y):

  * command-line `phonemize` program: new `--backend
    <espeak|festival|segments>` option, default language is now
    *espeak en-us* (was *festival en-us*),

  * it is now illegal to have the same separator at different levels
    (for instance a space for both word and phone),

  * from Python, must import the phonemize function as `from
    phonemizer.phonemize import phonemize`, was `from phonemizer
    import phonemize`.

* New backend [segments](https://github.com/cldf/segments) for
  phonemization based on grapheme-to-phoneme mappings.

* Major refactoring of the backends implementation and separators (as
  Python classes).

* Input to phonemizer now supports utf8.

* Better handling of errors (display of a meaningful message).

* Fixed a bug in fetching espeak version on macos, see
  [#14](https://github.com/bootphon/phonemizer/issues/14).

## phonemizer-0.3.3

* Fix a bug introduced in phonemizer-0.3.2 (apostrophes in festival
  backend). See [#12](https://github.com/bootphon/phonemizer/issues/12).


## phonemizer-0.3.2

* Continuous integration with tracis-ci.

* Support for docker.

* Better support for different versions of espeak/festival.

* Minor bugfixes and improved tests.


## phonemizer-0.3.1

* New espeak or espeak-ng backend with more than 100 languages.

* Support for Python 2.7 and 3.5.

* Integration with zenodo for citation.

* Various bugfixes and minor improvments.


## phonemizer-0.2

* First public release.

* Support for festival backend, American English only.
| **Tests** | [![Linux][badge-test-linux]](https://github.com/bootphon/phonemizer/actions/workflows/linux.yaml) [![MacOS][badge-test-macos]](https://github.com/bootphon/phonemizer/actions/workflows/macos.yaml) [![Windows][badge-test-windows]](https://github.com/bootphon/phonemizer/actions/workflows/windows.yaml) [![Codecov][badge-codecov]](https://codecov.io/gh/bootphon/phonemizer) |
| ---: | --- |
| **Release** | [![GitHub release (latest SemVer)][badge-github-version]](https://github.com/bootphon/phonemizer/releases/latest) [![PyPI][badge-pypi-version]](https://pypi.python.org/pypi/phonemizer) [![downloads][badge-pypi-downloads]](https://pypi.python.org/pypi/phonemizer) |
| **Citation** | [![status][badge-joss]](https://joss.theoj.org/papers/08d1ffc14f233f56942f78f3742b266e) [![DOI][badge-zenodo]](https://doi.org/10.5281/zenodo.1045825) |

---

# Phonemizer -- *foʊnmaɪzɚ*

* The phonemizer allows simple phonemization of words and texts in many languages.

* Provides both the `phonemize` command-line tool and the Python function
  `phonemizer.phonemize`. See [function documentation][phonemize-function].

* It is based on four backends: **espeak**, **espeak-mbrola**, **festival** and
  **segments**. The backends have different properties and capabilities resumed
  in table below. The backend choice is let to the user.

  * [espeak-ng](https://github.com/espeak-ng/espeak-ng) is a Text-to-Speech
    software supporting a lot of languages and IPA (International Phonetic
    Alphabet) output.

  * [espeak-ng-mbrola](https://github.com/espeak-ng/espeak-ng/blob/master/docs/mbrola.md)
    uses the SAMPA phonetic alphabet instead of IPA but does not preserve word
    boundaries.

  * [festival](http://www.cstr.ed.ac.uk/projects/festival) is another
    Tex-to-Speech engine. Its phonemizer backend currently supports only
    American English. It uses a [custom phoneset][festival-phoneset], but it
    allows tokenization at the syllable level.

  * [segments](https://github.com/cldf/segments) is a Unicode tokenizer that
    build a phonemization from a grapheme to phoneme mapping provided as a file
    by the user.

  |                              | espeak                   | espeak-mbrola           | festival                    | segments           |
  | ---:                         | ---                      | ---                     | ---                         | ---                |
  | **phone set**                | [IPA]                    | [SAMPA]                 | [custom][festival-phoneset] | user defined       |
  | **supported languages**      | [100+][espeak-languages] | [35][mbrola-languages] | US English                  | user defined       |
  | **processing speed**         | fast                     | slow                    | very slow                   | fast               |
  | **phone tokens**             | :heavy_check_mark:       | :heavy_check_mark:      | :heavy_check_mark:          | :heavy_check_mark: |
  | **syllable tokens**          | :x:                      | :x:                     | :heavy_check_mark:          | :x:                |
  | **word tokens**              | :heavy_check_mark:       | :x:                     | :heavy_check_mark:          | :heavy_check_mark: |
  | **punctuation preservation** | :heavy_check_mark:       | :x:                     | :heavy_check_mark:          | :heavy_check_mark: |
  | **stressed phones**          | :heavy_check_mark:       | :x:                     | :x:                         | :x:                |
  | [**tie**][tie-IPA]           | :heavy_check_mark:       | :x:                     | :x:                         | :x:                |


## Installation

### Dependencies

**You need python>=3.6.** If you really need to use python2, use the
[phonemizer-1.0] release.

You need to install
[festival](http://www.festvox.org/docs/manual-2.4.0/festival_6.html#Installation),
[espeak-ng](https://github.com/espeak-ng/espeak-ng#espeak-ng-text-to-speech)
and/or [mbrola](https://github.com/numediart/MBROLA) in order to use the
corresponding `phonemizer` backends. Follow instructions for your system below.

<details><summary>on Debian/Unbuntu</summary>
<p>

To install dependencies, simply run `sudo apt-get install festival espeak-ng mbrola`.

When using the **espeak-mbrola** backend, additional mbrola voices must be
installed (see
[here](https://github.com/espeak-ng/espeak-ng/blob/master/docs/mbrola.md)). List
the installable voices with `apt search mbrola`.

</p>
</details>

<details><summary>on CentOS/Fedora</summary>
<p>

To install dependencies, simply run `sudo yum install festival espeak-ng`.

When using the **espeak-mbrola** backend, the mbrola binary and additional
mbrola voices must be installed (see
[here](https://github.com/espeak-ng/espeak-ng/blob/master/docs/mbrola.md)).

</p>
</details>

<details><summary>on MacOS</summary>
<p>

**espeak** is available on brew at version 1.48: `brew install espeak`. If you
want a more recent version you have to [compile it from
sources](https://github.com/espeak-ng/espeak-ng/blob/master/docs/building.md#linux-mac-bsd).
To install **festival**, **mbrola** and additional mbrola voices, use the
script provided [here](https://github.com/pettarin/setup-festival-mbrola).

</p>
</details>

<details><summary>on Windows</summary>
<p>

Install **espeak-ng** with the `.msi` Windows installer provided with [espeak
releases](https://github.com/espeak-ng/espeak-ng/releases). **festival** must be
compiled from sources (see
[here](https://github.com/festvox/festival/blob/master/INSTALL) and
[here](https://www.eguidedog.net/doc/doc_build_win_festival.php)). **mbrola** is
not available for Windows.

</p>
</details>


### Phonemizer

* The simplest way is using pip:

        pip install phonemizer

* **OR** install it from sources with:

  ```shell
  git clone https://github.com/bootphon/phonemizer
  cd phonemizer
  python setup.py install
  ```

  If you experiment an error such as `ImportError: No module named setuptools`
  during installation, refeer to [issue
  #11](https://github.com/bootphon/phonemizer/issues/11).


### Docker image

Alternatively you can run the phonemizer within docker, using the
provided `Dockerfile`. To build the docker image, have a:

```shell
git clone https://github.com/bootphon/phonemizer
cd phonemizer
sudo docker build -t phonemizer .
```

Then run an interactive session with:

```shell
sudo docker run -it phonemizer /bin/bash
```


### Testing

When installed from sources or whithin a Docker image, you can run the tests
suite from the root `phonemizer` folder (once you installed `pytest`):

```shell
pip install pytest
pytest
```

### Developers

The `phonemizer` project is open-source and is welcoming contributions from
everyone. Please look at the [contributors guidelines](CONTRIBUTING.md) if you
wish to contribute.


## Citation

To refenrece the `phonemizer` in your own work, please cite the following [JOSS
paper](https://joss.theoj.org/papers/10.21105/joss.03958).

```bibtex
@article{Bernard2021,
  doi = {10.21105/joss.03958},
  url = {https://doi.org/10.21105/joss.03958},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {68},
  pages = {3958},
  author = {Mathieu Bernard and Hadrien Titeux},
  title = {Phonemizer: Text to Phones Transcription for Multiple Languages in Python},
  journal = {Journal of Open Source Software}
}
```

## Python usage

In Python import the `phonemize` function with `from phonemizer import
phonemize`. See the [function documentation][phonemize-function].


### Advice for best performances

It is much more efficient to minimize the number of calls to the `phonemize`
function. Indeed the initialization of the phonemization backend can be
expensive, especially for espeak. In one exemple:

```python
from phonemizer import phonemize

text = [line1, line2, ...]

# Do this:
phonemized = phonemize(text, ...)

# Not this:
phonemized = [phonemize(line, ...) for line in text]

# An alternative is to directly instanciate the backend and to call the
# phonemize function from it:
from phonemizer.backend import EspeakBackend
backend = EspeakBackend('en-us', ...)
phonemized = [backend.phonemize(line, ...) for line in text]
```

### Exemple 1: phonemize a text with festival

The following exemple downloads a text and phonemizes it using the festival
backend, preserving punctuation and using 4 jobs in parallel. The phones are not
separated, words are separated by a space and syllables by `|`.

```python
# need to pip install requests
import requests
from phonemizer import phonemize
from phonemizer.separator import Separator

# text is a list of 190 English sentences downloaded from github
url = (
    'https://gist.githubusercontent.com/CorentinJ/'
    '0bc27814d93510ae8b6fe4516dc6981d/raw/'
    'bb6e852b05f5bc918a9a3cb439afe7e2de570312/small_corpus.txt')
text = requests.get(url).content.decode()
text = [line.strip() for line in text.split('\n') if line]

# phn is a list of 190 phonemized sentences
phn = phonemize(
    text,
    language='en-us',
    backend='festival',
    separator=Separator(phone=None, word=' ', syllable='|'),
    strip=True,
    preserve_punctuation=True,
    njobs=4)
```

### Exemple 2: build a lexicon with espeak

The following exemple extracts a list of words present in a text, ignoring
punctuation, and builds a dictionary `word: [phones]`, e.g. `{'students': 's t
uː d ə n t s', 'cobb': 'k ɑː b', 'its': 'ɪ t s', 'put': 'p ʊ t', ...}`. We
consider here the same text as in the previous exemple.

```python
from phonemizer.backend import EspeakBackend
from phonemizer.punctuation import Punctuation
from phonemizer.separator import Separator

# remove all the punctuation from the text, condidering only the specified
# punctuation marks
text = Punctuation(';:,.!"?()-').remove(text)

# build the set of all the words in the text
words = {w.lower() for line in text for w in line.strip().split(' ') if w}

# initialize the espeak backend for English
backend = EspeakBackend('en-us')

# separate phones by a space and ignoring words boundaries
separator = Separator(phone=' ', word=None)

# build the lexicon by phonemizing each word one by one. The backend.phonemize
# function expect a list as input and outputs a list.
lexicon = {
    word: backend.phonemize([word], separator=separator, strip=True)[0]
    for word in words}
```

## Command-line examples

**The above examples can be run from Python using the `phonemize` function**


For a complete list of available options, have a:

```shell
phonemize --help
```

See the installed backends with the `--version` option:

```shell
$ phonemize --version
phonemizer-3.0
available backends: espeak-ng-1.50, espeak-mbrola, festival-2.5.0, segments-2.1.3
```


### Input/output exemples

* from stdin to stdout:

  ```shell
  $ echo "hello world" | phonemize
  həloʊ wɜːld
  ```

* Prepend the input text to output:

  ```shell
  $ echo "hello world" | phonemize --prepend-text
  hello world | həloʊ wɜːld

  $ echo "hello world" | phonemize --prepend-text=';'
  hello world ; həloʊ wɜːld
  ```

* from file to stdout

  ```shell
  $ echo "hello world" > hello.txt
  $ phonemize hello.txt
  həloʊ wɜːld
  ```

* from file to file

  ```shell
  $ phonemize hello.txt -o hello.phon --strip
  $ cat hello.phon
  həloʊ wɜːld
  ```


### Backends

* The default is to use **espeak** us-english:

  ```shell
  $ echo "hello world" | phonemize
  həloʊ wɜːld

  $ echo "hello world" | phonemize -l en-us -b espeak
  həloʊ wɜːld

  $ echo 'hello world' | phonemize -l en-us -b espeak --tie
  həlo͡ʊ wɜːld
  ```

* Use **festival** US English instead

  ```shell
  $ echo "hello world" | phonemize -l en-us -b festival
  hhaxlow werld
  ```

* In French, using **espeak** and **espeak-mbrola**, with custom token
  separators (see below). **espeak-mbrola** does not support words separation.

  ```shell
  $ echo "bonjour le monde" | phonemize -b espeak -l fr-fr -p ' ' -w '/w '
  b ɔ̃ ʒ u ʁ /w l ə /w m ɔ̃ d /w

  $ echo "bonjour le monde" | phonemize -b espeak-mbrola -l mb-fr1 -p ' ' -w '/w '
  b o~ Z u R l @ m o~ d
  ```

* In Japanese, using **segments**

  ```shell
  $ echo 'konnichiwa' | phonemize -b segments -l japanese
  konnitʃiwa

  $ echo 'konnichiwa' | phonemize -b segments -l ./phonemizer/share/japanese.g2p
  konnitʃiwa
  ```

### Supported languages

The exhaustive list of supported languages is available with the command
`phonemize --list-languages [--backend <backend>]`.

* Languages supported by **espeak** are available [here][espeak-languages].

* Languages supported by **espeak-mbrola** are available
  [here][mbrola-languages]. Please note that the mbrola voices are not bundled
  with the phonemizer nor the mbrola binary and must be installed separately.

* Languages supported by **festival** are:

  ```
  en-us -> english-us
  ```

* Languages supported by the **segments** backend are:

  ```
  chintang  -> ./phonemizer/share/segments/chintang.g2p
  cree      -> ./phonemizer/share/segments/cree.g2p
  inuktitut -> ./phonemizer/share/segments/inuktitut.g2p
  japanese  -> ./phonemizer/share/segments/japanese.g2p
  sesotho   -> ./phonemizer/share/segments/sesotho.g2p
  yucatec   -> ./phonemizer/share/segments/yucatec.g2p
  ```

  Instead of a language you can also provide a file specifying a
  grapheme to phone mapping (see the files above for examples).


### Token separators

You can specify separators for phones, syllables (**festival** only) and
words (excepted **espeak-mbrola**).

```shell
$ echo "hello world" | phonemize -b festival -w ' ' -p ''
hhaxlow werld

$ echo "hello world" | phonemize -b festival -p ' ' -w ''
hh ax l ow w er l d

$ echo "hello world" | phonemize -b festival -p '-' -s '|'
hh-ax-l-|ow-| w-er-l-d-|

$ echo "hello world" | phonemize -b festival -p '-' -s '|' --strip
hh-ax-l|ow w-er-l-d

$ echo "hello world" | phonemize -b festival -p ' ' -s ';esyll ' -w ';eword '
hh ax l ;esyll ow ;esyll ;eword w er l d ;esyll ;eword
```

You cannot specify the same separator for several tokens (for instance
a space for both phones and words):

```shell
$ echo "hello world" | phonemize -b festival -p ' ' -w ' '
fatal error: illegal separator with word=" ", syllable="" and phone=" ",
must be all differents if not empty
```


### Punctuation

By default the punctuation is removed in the phonemized output. You can preserve
it using the ``--preserve-punctuation`` option (not supported by the
**espeak-mbrola** backend):

```shell
$ echo "hello, world!" | phonemize --strip
həloʊ wɜːld

$ echo "hello, world!" | phonemize --preserve-punctuation --strip
həloʊ, wɜːld!
```


### Espeak specific options

* The espeak backend can output the **stresses** on phones:

  ```shell
  $ echo "hello world" | phonemize -l en-us -b espeak --with-stress
  həlˈoʊ wˈɜːld
  ```

* The espeak backend can add **tie** on multi-characters phonemes:

  ```shell
  $ echo "hello world" | phonemize -l en-us -b espeak --tie
  həlo͡ʊ wɜːld
  ```

* :warning: The espeak backend can **switch languages** during phonemization (below from
  French to English), use the ``--language-switch`` option to deal with it:

  ```shell
  $ echo "j'aime le football" | phonemize -l fr-fr -b espeak --language-switch keep-flags
  [WARNING] fount 1 utterances containing language switches on lines 1
  [WARNING] extra phones may appear in the "fr-fr" phoneset
  [WARNING] language switch flags have been kept (applying "keep-flags" policy)
  ʒɛm lə- (en)fʊtbɔːl(fr)

  $ echo "j'aime le football" | phonemize -l fr-fr -b espeak --language-switch remove-flags
  [WARNING] fount 1 utterances containing language switches on lines 1
  [WARNING] extra phones may appear in the "fr-fr" phoneset
  [WARNING] language switch flags have been removed (applying "remove-flags" policy)
  ʒɛm lə- fʊtbɔːl

  $ echo "j'aime le football" | phonemize -l fr-fr -b espeak --language-switch remove-utterance
  [WARNING] removed 1 utterances containing language switches (applying "remove-utterance" policy)
  ```

* :warning: The espeak backend sometimes **merge words together** in the output, use the
  `--words-mismatch` option to deal with it:

  ```shell
  $ echo "that's it, words are merged" | phonemize -l en-us -b espeak
  [WARNING] words count mismatch on 100.0% of the lines (1/1)
  ðætsɪt wɜːdz ɑːɹ mɜːdʒd
  ```


## Licence

**Copyright 2015-2021 Mathieu Bernard**

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.


[badge-test-linux]: https://github.com/bootphon/phonemizer/actions/workflows/linux.yaml/badge.svg?branch=master
[badge-test-macos]: https://github.com/bootphon/phonemizer/actions/workflows/macos.yaml/badge.svg?branch=master
[badge-test-windows]: https://github.com/bootphon/phonemizer/actions/workflows/windows.yaml/badge.svg?branch=master
[badge-codecov]: https://img.shields.io/codecov/c/github/bootphon/phonemizer
[badge-github-version]: https://img.shields.io/github/v/release/bootphon/phonemizer
[badge-pypi-version]: https://img.shields.io/pypi/v/phonemizer
[badge-pypi-downloads]: https://img.shields.io/pypi/dm/phonemizer
[badge-joss]: https://joss.theoj.org/papers/08d1ffc14f233f56942f78f3742b266e/status.svg
[badge-zenodo]: https://zenodo.org/badge/56728069.svg
[phonemizer-1.0]: https://github.com/bootphon/phonemizer/releases/tag/v1.0
[festival-phoneset]: http://www.festvox.org/bsv/c4711.html
[IPA]: https://en.wikipedia.org/wiki/International_Phonetic_Alphabet
[SAMPA]: https://en.wikipedia.org/wiki/SAMPA
[phonemize-function]: https://github.com/bootphon/phonemizer/blob/c5e2f3878d6db391ec7253173f44e4a85cfe41e3/phonemizer/phonemize.py#L33-L156
[tie-IPA]: https://en.wikipedia.org/wiki/Tie_(typography)#International_Phonetic_Alphabet
[espeak-languages]: https://github.com/espeak-ng/espeak-ng/blob/master/docs/languages.md
[mbrola-languages]: https://github.com/numediart/MBROLA-voices
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, caste, color, religion, or sexual
identity and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
  community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or advances of
  any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email address,
  without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to Mathieu Bernard at the email address provided on it's [GitHub
profile page](https://github.com/mmmaat). All complaints will be reviewed and
investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.1, available at
[https://www.contributor-covenant.org/version/2/1/code_of_conduct.html][v2.1].

Community Impact Guidelines were inspired by
[Mozilla's code of conduct enforcement ladder][Mozilla CoC].

For answers to common questions about this code of conduct, see the FAQ at
[https://www.contributor-covenant.org/faq][FAQ]. Translations are available at
[https://www.contributor-covenant.org/translations][translations].

[homepage]: https://www.contributor-covenant.org
[v2.1]: https://www.contributor-covenant.org/version/2/1/code_of_conduct.html
[Mozilla CoC]: https://github.com/mozilla/diversity
[FAQ]: https://www.contributor-covenant.org/faq
[translations]: https://www.contributor-covenant.org/translations
# How to contribute?

We welcome and encourage every bug report, feature and pull request.


## You have an issue with `phonemizer` or a feature request

Please open an [issue on github](https://github.com/bootphon/phonemizer/issues)
and follow the template from there.


## You want to contribute code

If you're willing to take it upon yourself to improve `phonemizer`, via
bugfixes, improvements and new features, please follow these steps:

- Submit an issue explaining what you're willing to fix or add to this package.
  We can discuss with you on the the best way to do it, considering the current
  state of things.

- Fork the `phonemizer` repo, code away and open a pull-request. If you add some
  code or change significantly a function, please test it by adding more unit
  tests.

- Ensure the tests are passing on Linux, MacOS and Windows on the [github
  continuous integration page](https://github.com/bootphon/phonemizer/actions).
  This is important because `phonemizer` is used on a wide range of systems and
  configurations.

- Please conform to the following conventions:

    - Python code follows [PEP 8 style](https://pep8.org).
    - Docstrings follow [Google
      style](https://google.github.io/styleguide/pyguide.html#s3.8-comments-and-docstrings).
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**Phonemizer version**
The output of `phonemize --version` from command line, very helpfull!

**System**
Your OS (Linux distribution, Windows, ...), eventually Python version.

**To reproduce**
A short example (Python script or command) reproducing the bug.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
