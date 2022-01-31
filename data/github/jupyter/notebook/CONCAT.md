# Changelog

A summary of changes in the Jupyter notebook. For more detailed
information, see [GitHub](https://github.com/jupyter/notebook).

Use `pip install notebook --upgrade` or `conda upgrade notebook` to
upgrade to the latest release.

We strongly recommend that you upgrade pip to version 9+ of pip before
upgrading `notebook`.

Use `pip install pip --upgrade` to upgrade pip. Check pip version with
`pip --version`.

<!-- <START NEW CHANGELOG ENTRY> -->

## 6.4.8

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.7...479902d83a691253e0cff8439a33577e82408317))

### Bugs fixed

- Fix to remove potential memory leak on Jupyter Notebooks ZMQChannelHandler code [#6251](https://github.com/jupyter/notebook/pull/6251) ([@Vishwajeet0510](https://github.com/Vishwajeet0510))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2022-01-12&to=2022-01-25&type=c))

[@Vishwajeet0510](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AVishwajeet0510+updated%3A2022-01-12..2022-01-25&type=Issues)

<!-- <END NEW CHANGELOG ENTRY> -->

## 6.4.7

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.6...b77b5e38b8fa1a20150d7fa4d735dbf1c4f00418))

### Bugs fixed

- Fix Chinese punctuation [#6268](https://github.com/jupyter/notebook/pull/6268) ([@LiHua-Official](https://github.com/LiHua-Official))
- Add date field to kernel message header [#6265](https://github.com/jupyter/notebook/pull/6265) ([@kevin-bates](https://github.com/kevin-bates))
- Fix deprecation warning [#6253](https://github.com/jupyter/notebook/pull/6253) ([@tornaria](https://github.com/tornaria))

### Maintenance and upkeep improvements

- Enforce labels on PRs [#6235](https://github.com/jupyter/notebook/pull/6235) ([@blink1073](https://github.com/blink1073))
- Fix: CI error for python 3.6 & macOS [#6215](https://github.com/jupyter/notebook/pull/6215) ([@penguinolog](https://github.com/penguinolog))

### Other merged PRs

- handle KeyError when get session [#6245](https://github.com/jupyter/notebook/pull/6245) ([@ccw630](https://github.com/ccw630))
- Updated doc for passwd [#6209](https://github.com/jupyter/notebook/pull/6209) ([@antoinecarme](https://github.com/antoinecarme))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-11-16&to=2022-01-12&type=c))

[@antoinecarme](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aantoinecarme+updated%3A2021-11-16..2022-01-12&type=Issues) | [@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-11-16..2022-01-12&type=Issues) | [@ccw630](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Accw630+updated%3A2021-11-16..2022-01-12&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-11-16..2022-01-12&type=Issues) | [@LiHua-Official](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3ALiHua-Official+updated%3A2021-11-16..2022-01-12&type=Issues) | [@penguinolog](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Apenguinolog+updated%3A2021-11-16..2022-01-12&type=Issues) | [@tornaria](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Atornaria+updated%3A2021-11-16..2022-01-12&type=Issues)

## 6.4.6

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.5...160c27d3c23dafe8b42240571db21b0d5cbae2fe))

### Bugs fixed

- Fix `asyncio` error when opening notebooks [#6221](https://github.com/jupyter/notebook/pull/6221) ([@dleen](https://github.com/dleen))
- Change to use a universal Chinese translation on certain words [#6218](https://github.com/jupyter/notebook/pull/6218) ([@jackexu](https://github.com/jackexu))
- Fix Chinese translation typo [#6211](https://github.com/jupyter/notebook/pull/6211) ([@maliubiao](https://github.com/maliubiao)
- Fix `send2trash` tests failing on Windows [#6127](https://github.com/jupyter/notebook/pull/6127) ([@dolfinus](https://github.com/dolfinus))

### Maintenance and upkeep improvements

- TST: don't look in user site for serverextensions [#6233](https://github.com/jupyter/notebook/pull/6233) ([@bnavigator](https://github.com/bnavigator))
- Enable terminal tests as `pywinpty` is ported for python 3.9 [#6228](https://github.com/jupyter/notebook/pull/6228) ([@nsait-linaro](https://github.com/nsait-linaro))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-10-19&to=2021-11-16&type=c))

[@bnavigator](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Abnavigator+updated%3A2021-10-19..2021-11-16&type=Issues) | [@dleen](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Adleen+updated%3A2021-10-19..2021-11-16&type=Issues) | [@dolfinus](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Adolfinus+updated%3A2021-10-19..2021-11-16&type=Issues) | [@jackexu](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ajackexu+updated%3A2021-10-19..2021-11-16&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-10-19..2021-11-16&type=Issues) | [@maliubiao](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amaliubiao+updated%3A2021-10-19..2021-11-16&type=Issues) | [@nsait-linaro](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ansait-linaro+updated%3A2021-10-19..2021-11-16&type=Issues) | [@takluyver](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Atakluyver+updated%3A2021-10-19..2021-11-16&type=Issues) | [@Zsailer](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AZsailer+updated%3A2021-10-19..2021-11-16&type=Issues)

## 6.4.5

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.4...ccd9665571107e02a325a738b8baebd6532b2d3d))

### Bug fixes

- Recover from failure to render mimetype [#6181](https://github.com/jupyter/notebook/pull/6181) ([@martinRenou](https://github.com/martinRenou))

### Maintenance and upkeep improvements

- Fix crypto handling [#6197](https://github.com/jupyter/notebook/pull/6197) ([@blink1073](https://github.com/blink1073))
- Fix `jupyter_client` warning [#6178](https://github.com/jupyter/notebook/pull/6178) ([@martinRenou](https://github.com/martinRenou))

### Documentation improvements

- Fix nbsphinx settings [#6200](https://github.com/jupyter/notebook/pull/6200) ([@mgeier](https://github.com/mgeier))
- Fully revert the pinning of `nbsphinx` to 0.8.6 [#6201](https://github.com/jupyter/notebook/pull/6201) ([@kevin-bates](https://github.com/kevin-bates))
- Pin `nbsphinx` to 0.8.6, clean up orphaned resources [#6194](https://github.com/jupyter/notebook/pull/6194) ([@kevin-bates](https://github.com/kevin-bates))
- Fix typo in docstring [#6188](https://github.com/jupyter/notebook/pull/6188) ([@jgarte](https://github.com/jgarte))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-09-03&to=2021-10-19&type=c))

[@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-09-03..2021-10-19&type=Issues) | [@jgarte](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ajgarte+updated%3A2021-09-03..2021-10-19&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-09-03..2021-10-19&type=Issues) | [@martinRenou](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AmartinRenou+updated%3A2021-09-03..2021-10-19&type=Issues) | [@mgeier](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amgeier+updated%3A2021-09-03..2021-10-19&type=Issues)

## 6.4.4

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.3...c06c340574e1d2207940c5bd1190eb73d82ab945))

### Documentation improvements

- Update Manual Release Instructions [#6152](https://github.com/jupyter/notebook/pull/6152) ([@blink1073](https://github.com/blink1073))

### Other merged PRs

- Use default JupyterLab CSS sanitizer options for Markdown [#6160](https://github.com/jupyter/notebook/pull/6160) ([@krassowski](https://github.com/krassowski))
- Fix syntax highlight [#6128](https://github.com/jupyter/notebook/pull/6128) ([@massongit](https://github.com/massongit))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-08-11&to=2021-09-03&type=c))

[@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-08-11..2021-09-03&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-08-11..2021-09-03&type=Issues) | [@krassowski](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akrassowski+updated%3A2021-08-11..2021-09-03&type=Issues) | [@massongit](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amassongit+updated%3A2021-08-11..2021-09-03&type=Issues) | [@minrk](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aminrk+updated%3A2021-08-11..2021-09-03&type=Issues) | [@Zsailer](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AZsailer+updated%3A2021-08-11..2021-09-03&type=Issues)

## 6.4.3

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.2...c373bd89adaaddffbb71747ebbcfe8a749cae0a8))

### Bugs fixed

- Add @babel/core dependency [#6133](https://github.com/jupyter/notebook/pull/6133) ([@afshin](https://github.com/afshin))
- Switch webpack to production mode [#6131](https://github.com/jupyter/notebook/pull/6131) ([@afshin](https://github.com/afshin))

### Maintenance and upkeep improvements

- Clean up link checking [#6130](https://github.com/jupyter/notebook/pull/6130) ([@blink1073](https://github.com/blink1073))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-08-06&to=2021-08-10&type=c))

[@afshin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aafshin+updated%3A2021-08-06..2021-08-10&type=Issues) | [@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-08-06..2021-08-10&type=Issues) | [@Zsailer](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AZsailer+updated%3A2021-08-06..2021-08-10&type=Issues)

## 6.4.2

([Full Changelog](https://github.com/jupyter/notebook/compare/v6.4.0...999e8322bcd24e0ed62b180c19ec13db3f48165b))

### Bugs fixed

- Add missing file to manifest [#6122](https://github.com/jupyter/notebook/pull/6122) ([@afshin](https://github.com/afshin))
- Fix issue #3218 [#6108](https://github.com/jupyter/notebook/pull/6108) ([@Nazeeh21](https://github.com/Nazeeh21))
- Fix version of jupyter-packaging in pyproject.toml [#6101](https://github.com/jupyter/notebook/pull/6101) ([@frenzymadness](https://github.com/frenzymadness))
- "#element".tooltip is not a function on home page fixed. [#6070](https://github.com/jupyter/notebook/pull/6070) ([@ilayh123](https://github.com/ilayh123))

### Maintenance and upkeep improvements

- Enhancements to the desktop entry [#6099](https://github.com/jupyter/notebook/pull/6099) ([@Amr-Ibra](https://github.com/Amr-Ibra))
- Add missing spaces to help messages in config file [#6085](https://github.com/jupyter/notebook/pull/6085) ([@saiwing-yeung](https://github.com/saiwing-yeung))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-05-17&to=2021-08-06&type=c))

[@afshin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aafshin+updated%3A2021-05-17..2021-08-06&type=Issues) | [@Amr-Ibra](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AAmr-Ibra+updated%3A2021-05-17..2021-08-06&type=Issues) | [@frenzymadness](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Afrenzymadness+updated%3A2021-05-17..2021-08-06&type=Issues) | [@ilayh123](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ailayh123+updated%3A2021-05-17..2021-08-06&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-05-17..2021-08-06&type=Issues) | [@Nazeeh21](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3ANazeeh21+updated%3A2021-05-17..2021-08-06&type=Issues) | [@saiwing-yeung](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Asaiwing-yeung+updated%3A2021-05-17..2021-08-06&type=Issues)

## 6.4.0

([Full Changelog](https://github.com/jupyter/notebook/compare/6.3.0...80eb286f316838afc76a9a84b06f54e7dccb6c86))

### Bugs fixed

- Fix Handling of Encoded Paths in Save As Dialog [#6030](https://github.com/jupyter/notebook/pull/6030) ([@afshin](https://github.com/afshin))
- Fix: split_cell doesn't always split cell [#6017](https://github.com/jupyter/notebook/pull/6017) ([@gamestrRUS](https://github.com/gamestrRUS))
- Correct 'Content-Type' headers [#6026](https://github.com/jupyter/notebook/pull/6026) ([@faucct](https://github.com/faucct))
- Fix skipped tests & remove deprecation warnings [#6018](https://github.com/jupyter/notebook/pull/6018) ([@befeleme](https://github.com/befeleme))
- [Gateway] Track only this server's kernels [#5980](https://github.com/jupyter/notebook/pull/5980) ([@kevin-bates](https://github.com/kevin-bates))
- Bind the HTTPServer in start [#6061](https://github.com/jupyter/notebook/pull/6061)

### Maintenance and upkeep improvements

- Revert "do not apply asyncio patch for tornado >=6.1" [#6052](https://github.com/jupyter/notebook/pull/6052) ([@minrk](https://github.com/minrk))
- Use Jupyter Releaser [#6048](https://github.com/jupyter/notebook/pull/6048) ([@afshin](https://github.com/afshin))
- Add Workflow Permissions for Lock Bot [#6042](https://github.com/jupyter/notebook/pull/6042) ([@jtpio](https://github.com/jtpio))
- Fixes related to the recent changes in the documentation [#6021](https://github.com/jupyter/notebook/pull/6021) ([@frenzymadness](https://github.com/frenzymadness))
- Add maths checks in CSS reference test [#6035](https://github.com/jupyter/notebook/pull/6035) ([@stef4k](https://github.com/stef4k))
- Add Issue Lock and Answered Bots [#6019](https://github.com/jupyter/notebook/pull/6019) ([@afshin](https://github.com/afshin))

### Documentation improvements

- Spelling correction [#6045](https://github.com/jupyter/notebook/pull/6045) ([@wggillen](https://github.com/wggillen))
- Minor typographical and comment changes [#6025](https://github.com/jupyter/notebook/pull/6025) ([@misterhay](https://github.com/misterhay))
- Fixes related to the recent changes in the documentation [#6021](https://github.com/jupyter/notebook/pull/6021) ([@frenzymadness](https://github.com/frenzymadness))
- Fix readthedocs environment [#6020](https://github.com/jupyter/notebook/pull/6020) ([@blink1073](https://github.com/blink1073))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-03-22&to=2021-05-12&type=c))

[@afshin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aafshin+updated%3A2021-03-22..2021-05-12&type=Issues) | [@befeleme](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Abefeleme+updated%3A2021-03-22..2021-05-12&type=Issues) | [@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-03-22..2021-05-12&type=Issues) | [@faucct](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Afaucct+updated%3A2021-03-22..2021-05-12&type=Issues) | [@frenzymadness](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Afrenzymadness+updated%3A2021-03-22..2021-05-12&type=Issues) | [@gamestrRUS](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AgamestrRUS+updated%3A2021-03-22..2021-05-12&type=Issues) | [@jtpio](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ajtpio+updated%3A2021-03-22..2021-05-12&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-03-22..2021-05-12&type=Issues) | [@minrk](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aminrk+updated%3A2021-03-22..2021-05-12&type=Issues) | [@misterhay](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amisterhay+updated%3A2021-03-22..2021-05-12&type=Issues) | [@stef4k](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Astef4k+updated%3A2021-03-22..2021-05-12&type=Issues) | [@wggillen](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Awggillen+updated%3A2021-03-22..2021-05-12&type=Issues)

## 6.3.0

### Merged PRs

* Add square logo and desktop entry files [#6010](https://github.com/jupyter/notebook/pull/6010) ([@befeleme](https://github.com/befeleme))
* Modernize Changelog [#6008](https://github.com/jupyter/notebook/pull/6008) ([@afshin](https://github.com/afshin))
* Add missing "import inspect" [#5999](https://github.com/jupyter/notebook/pull/5999) ([@mgeier](https://github.com/mgeier))
* Add Codecov badge to README [#5989](https://github.com/jupyter/notebook/pull/5989) ([@thomasrockhu](https://github.com/thomasrockhu))
* Remove configuration for nosetests from setup.cfg [#5986](https://github.com/jupyter/notebook/pull/5986) ([@frenzymadness](https://github.com/frenzymadness))
* Update security.rst [#5978](https://github.com/jupyter/notebook/pull/5978) ([@dlrice](https://github.com/dlrice))
*  Docs-Translations: Updated Hindi and Chinese Readme.md [#5976](https://github.com/jupyter/notebook/pull/5976) ([@rjn01](https://github.com/rjn01))
* Allow /metrics by default if auth is off [#5974](https://github.com/jupyter/notebook/pull/5974) ([@blairdrummond](https://github.com/blairdrummond))
* Skip terminal tests on Windows 3.9+ (temporary) [#5968](https://github.com/jupyter/notebook/pull/5968) ([@kevin-bates](https://github.com/kevin-bates))
* Update GatewayKernelManager to derive from AsyncMappingKernelManager [#5966](https://github.com/jupyter/notebook/pull/5966) ([@kevin-bates](https://github.com/kevin-bates))
* Drop use of deprecated pyzmq.ioloop [#5965](https://github.com/jupyter/notebook/pull/5965) ([@kevin-bates](https://github.com/kevin-bates))
* Drop support for Python 3.5 [#5962](https://github.com/jupyter/notebook/pull/5962) ([@kevin-bates](https://github.com/kevin-bates))
* Allow jupyter_server-based contents managers in notebook [#5957](https://github.com/jupyter/notebook/pull/5957) ([@afshin](https://github.com/afshin))
* Russian translation fixes [#5954](https://github.com/jupyter/notebook/pull/5954) ([@insolor](https://github.com/insolor))
* Increase culling test idle timeout [#5952](https://github.com/jupyter/notebook/pull/5952) ([@kevin-bates](https://github.com/kevin-bates))
* Re-enable support for answer_yes flag [#5941](https://github.com/jupyter/notebook/pull/5941) ([@afshin](https://github.com/afshin))
* Replace Travis and Appveyor with Github Actions [#5938](https://github.com/jupyter/notebook/pull/5938) ([@kevin-bates](https://github.com/kevin-bates))
* DOC: Server extension, extra docs on configuration/authentication. [#5937](https://github.com/jupyter/notebook/pull/5937) ([@Carreau](https://github.com/Carreau))

### Contributors to this release

([GitHub contributors page for this release](https://github.com/jupyter/notebook/graphs/contributors?from=2021-01-13&to=2021-03-18&type=c))

[@abielhammonds](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aabielhammonds+updated%3A2021-01-13..2021-03-18&type=Issues) | [@afshin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aafshin+updated%3A2021-01-13..2021-03-18&type=Issues) | [@ajharry](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aajharry+updated%3A2021-01-13..2021-03-18&type=Issues) | [@Alokrar](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AAlokrar+updated%3A2021-01-13..2021-03-18&type=Issues) | [@befeleme](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Abefeleme+updated%3A2021-01-13..2021-03-18&type=Issues) | [@blairdrummond](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablairdrummond+updated%3A2021-01-13..2021-03-18&type=Issues) | [@blink1073](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ablink1073+updated%3A2021-01-13..2021-03-18&type=Issues) | [@bollwyvl](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Abollwyvl+updated%3A2021-01-13..2021-03-18&type=Issues) | [@Carreau](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3ACarreau+updated%3A2021-01-13..2021-03-18&type=Issues) | [@ChenChenDS](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AChenChenDS+updated%3A2021-01-13..2021-03-18&type=Issues) | [@cosmoscalibur](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Acosmoscalibur+updated%3A2021-01-13..2021-03-18&type=Issues) | [@dlrice](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Adlrice+updated%3A2021-01-13..2021-03-18&type=Issues) | [@dwanneruchi](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Adwanneruchi+updated%3A2021-01-13..2021-03-18&type=Issues) | [@ElisonSherton](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AElisonSherton+updated%3A2021-01-13..2021-03-18&type=Issues) | [@FazeelUsmani](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AFazeelUsmani+updated%3A2021-01-13..2021-03-18&type=Issues) | [@frenzymadness](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Afrenzymadness+updated%3A2021-01-13..2021-03-18&type=Issues) | [@goerz](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Agoerz+updated%3A2021-01-13..2021-03-18&type=Issues) | [@insolor](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ainsolor+updated%3A2021-01-13..2021-03-18&type=Issues) | [@jasongrout](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ajasongrout+updated%3A2021-01-13..2021-03-18&type=Issues) | [@JianghuiDu](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AJianghuiDu+updated%3A2021-01-13..2021-03-18&type=Issues) | [@JuzerShakir](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AJuzerShakir+updated%3A2021-01-13..2021-03-18&type=Issues) | [@kevin-bates](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Akevin-bates+updated%3A2021-01-13..2021-03-18&type=Issues) | [@Khalilsqu](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AKhalilsqu+updated%3A2021-01-13..2021-03-18&type=Issues) | [@meeseeksdev](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ameeseeksdev+updated%3A2021-01-13..2021-03-18&type=Issues) | [@mgeier](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amgeier+updated%3A2021-01-13..2021-03-18&type=Issues) | [@michaelpedota](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amichaelpedota+updated%3A2021-01-13..2021-03-18&type=Issues) | [@mjbright](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Amjbright+updated%3A2021-01-13..2021-03-18&type=Issues) | [@MSeal](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AMSeal+updated%3A2021-01-13..2021-03-18&type=Issues) | [@ncoughlin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ancoughlin+updated%3A2021-01-13..2021-03-18&type=Issues) | [@NTimmons](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3ANTimmons+updated%3A2021-01-13..2021-03-18&type=Issues) | [@ProsperousHeart](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AProsperousHeart+updated%3A2021-01-13..2021-03-18&type=Issues) | [@rjn01](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Arjn01+updated%3A2021-01-13..2021-03-18&type=Issues) | [@slw07g](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Aslw07g+updated%3A2021-01-13..2021-03-18&type=Issues) | [@stenivan](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Astenivan+updated%3A2021-01-13..2021-03-18&type=Issues) | [@takluyver](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Atakluyver+updated%3A2021-01-13..2021-03-18&type=Issues) | [@thomasrockhu](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Athomasrockhu+updated%3A2021-01-13..2021-03-18&type=Issues) | [@wgilpin](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Awgilpin+updated%3A2021-01-13..2021-03-18&type=Issues) | [@wxtt522](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Awxtt522+updated%3A2021-01-13..2021-03-18&type=Issues) | [@yuvipanda](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3Ayuvipanda+updated%3A2021-01-13..2021-03-18&type=Issues) | [@Zsailer](https://github.com/search?q=repo%3Ajupyter%2Fnotebook+involves%3AZsailer+updated%3A2021-01-13..2021-03-18&type=Issues)

## 6.2.0

## Merged PRs

- Increase minimum tornado version ([5933](https://github.com/jupyter/notebook/pull/5933))
- Adjust skip decorators to avoid remaining dependency on nose ([5932](https://github.com/jupyter/notebook/pull/5932))
- Ensure that cell ids persist after save ([5928](https://github.com/jupyter/notebook/pull/5928))
- Add reconnection to Gateway (form nb2kg) ([5924](https://github.com/jupyter/notebook/pull/5924))
- Fix some typos ([5917](https://github.com/jupyter/notebook/pull/5917))
- Handle TrashPermissionError, now that it exist ([5894](https://github.com/jupyter/notebook/pull/5894))

Thank you to all the contributors:

- @kevin-bates
- @mishaschwartz
- @oyvsyo
- @user202729
- @stefanor

## 6.1.6

## Merged PRs

-  do not require nose for testing ([5826](https://github.com/jupyter/notebook/pull/5826))
- [docs] Update Chinese and Hindi readme.md ([5823](https://github.com/jupyter/notebook/pull/5823))
- Add support for creating terminals via GET ([5813](https://github.com/jupyter/notebook/pull/5813))
- Made doc translations in Hindi and Chinese ([5787](https://github.com/jupyter/notebook/pull/5787))

Thank you to all the contributors:

- @pgajdos
- @rjn01
- @kevin-bates
- @virejdasani

## 6.1.5

6.1.5 is a security release, fixing one vulnerability:

- Fix open redirect vulnerability GHSA-c7vm-f5p4-8fqh (CVE to be assigned)

## 6.1.4

- Fix broken links to jupyter documentation ([5686](https://github.com/jupyter/notebook/pull/5686))
- Add additional entries to troubleshooting section ([5695](https://github.com/jupyter/notebook/pull/5695))
- Revert change in page alignment ([5703](https://github.com/jupyter/notebook/pull/5703))
- Bug fix: remove double encoding in download files ([5720](https://github.com/jupyter/notebook/pull/5720))
- Fix typo for Check in zh_CN ([5730](https://github.com/jupyter/notebook/pull/5730))
- Require a file name in the "Save As" dialog ([5733](https://github.com/jupyter/notebook/pull/5733))

Thank you to all the contributors:

- bdbai
- Jaipreet Singh
- Kevin Bates
- Pavel Panchekha
- Zach Sailer

## 6.1.3

- Title new buttons with label if action undefined ([5676](https://github.com/jupyter/notebook/pull/5676))

Thank you to all the contributors:

- Kyle Kelley

## 6.1.2

- Fix russian message format for delete/duplicate actions ([5662](https://github.com/jupyter/notebook/pull/5662))
- Remove unnecessary import of bind_unix_socket ([5666](https://github.com/jupyter/notebook/pull/5666))
- Tooltip style scope fix ([5672](https://github.com/jupyter/notebook/pull/5672))

Thank you to all the contributors:

- Dmitry Akatov
- Kevin Bates
- Magda Stenius

## 6.1.1

- Prevent inclusion of requests_unixsocket on Windows ([5650](https://github.com/jupyter/notebook/pull/5650))

Thank you to all the contributors:

- Kevin Bates

## 6.1.0

Please note that this repository is currently maintained by a skeleton
crew of maintainers from the Jupyter community. For our approach moving
forward, please see this
[notice](https://github.com/jupyter/notebook#notice) from the README.
Thank you.

Here is an enumeration of changes made since the last release and
included in 6.1.0.

- Remove deprecated encoding parameter for Python 3.9 compatibility. ([5174](https://github.com/jupyter/notebook/pull/5174))
- Add support for async kernel management ([4479](https://github.com/jupyter/notebook/pull/4479))
- Fix typo in password_required help message ([5320](https://github.com/jupyter/notebook/pull/5320))
- Gateway only: Ensure launch and request timeouts are in sync ([5317](https://github.com/jupyter/notebook/pull/5317))
- Update Markdown Cells example to HTML5 video tag ([5411](https://github.com/jupyter/notebook/pull/5411))
- Integrated LoginWidget into edit to enable users to logout from the t... ([5406](https://github.com/jupyter/notebook/pull/5406))
- Update message about minimum Tornado version ([5222](https://github.com/jupyter/notebook/pull/5222))
- Logged notebook type ([5425](https://github.com/jupyter/notebook/pull/5425))
- Added nl language ([5354](https://github.com/jupyter/notebook/pull/5354))
- Add UNIX socket support to notebook server. ([4835](https://github.com/jupyter/notebook/pull/4835))
- Update CodeMirror dependency ([5198](https://github.com/jupyter/notebook/pull/5198))
- Tree added download multiple files ([5351](https://github.com/jupyter/notebook/pull/5351))
- Toolbar buttons tooltip: show help instead of label ([5107](https://github.com/jupyter/notebook/pull/5107))
- Remove unnecessary import of requests_unixsocket ([5451](https://github.com/jupyter/notebook/pull/5451))
- Add ability to cull terminals and track last activity ([5372](https://github.com/jupyter/notebook/pull/5372))
- Code refactoring notebook.js ([5352](https://github.com/jupyter/notebook/pull/5352))
- Install terminado for docs build ([5462](https://github.com/jupyter/notebook/pull/5462))
- Convert notifications JS test to selenium ([5455](https://github.com/jupyter/notebook/pull/5455))
- Add cell attachments to markdown example ([5412](https://github.com/jupyter/notebook/pull/5412))
- Add Japanese document ([5231](https://github.com/jupyter/notebook/pull/5231))
- Migrate Move multiselection test to selenium ([5158](https://github.com/jupyter/notebook/pull/5158))
- Use `cmdtrl-enter` to run a cell ([5120](https://github.com/jupyter/notebook/pull/5120))
- Fix broken "Raw cell MIME type" dialog ([5385](https://github.com/jupyter/notebook/pull/5385))
- Make a notebook writable after successful save-as ([5296](https://github.com/jupyter/notebook/pull/5296))
- Add actual watch script ([4738](https://github.com/jupyter/notebook/pull/4738))
- Added `--autoreload` flag to `NotebookApp` ([4795](https://github.com/jupyter/notebook/pull/4795))
- Enable check_origin on gateway websocket communication ([5471](https://github.com/jupyter/notebook/pull/5471))
- Restore detection of missing terminado package ([5465](https://github.com/jupyter/notebook/pull/5465))
- Culling: ensure `last_activity` attr exists before use ([5355](https://github.com/jupyter/notebook/pull/5355))
- Added functionality to allow filter kernels by Jupyter Enterprise Gat... ([5484](https://github.com/jupyter/notebook/pull/5484))
- 'Play' icon for run-cell toolbar button ([2922](https://github.com/jupyter/notebook/pull/2922))
- Bump minimum version of jQuery to 3.5.0 ([5491](https://github.com/jupyter/notebook/pull/5491))
- Remove old JS markdown tests, add a new one in selenium ([5497](https://github.com/jupyter/notebook/pull/5497))
- Add support for more RTL languages ([5036](https://github.com/jupyter/notebook/pull/5036))
- Make markdown cells stay RTL in edit mode ([5037](https://github.com/jupyter/notebook/pull/5037))
- Unforce RTL output display ([5039](https://github.com/jupyter/notebook/pull/5039))
- Fixed multicursor backspacing ([4880](https://github.com/jupyter/notebook/pull/4880))
- Implemented Split Cell for multicursor ([4824](https://github.com/jupyter/notebook/pull/4824))
- Alignment issue \[FIXED\] ([3173](https://github.com/jupyter/notebook/pull/3173))
- MathJax: Support for `\gdef` ([4407](https://github.com/jupyter/notebook/pull/4407))
- Another (Minor) Duplicate Code Reduction ([5316](https://github.com/jupyter/notebook/pull/5316))
- Update readme regarding maintenance ([5500](https://github.com/jupyter/notebook/pull/5500))
- Document contents chunks ([5508](https://github.com/jupyter/notebook/pull/5508))
- Backspace deletes empty line ([5516](https://github.com/jupyter/notebook/pull/5516))
- The dropdown submenu at notebook page is not keyboard accessible ([4732](https://github.com/jupyter/notebook/pull/4732))
- Tooltips visible through keyboard navigation for specified buttons ([4729](https://github.com/jupyter/notebook/pull/4729))
- Fix for recursive symlink ([4670](https://github.com/jupyter/notebook/pull/4670))
- Fix for the terminal shutdown issue ([4180](https://github.com/jupyter/notebook/pull/4180))
- Add japanese translation files ([4490](https://github.com/jupyter/notebook/pull/4490))
- Workaround for socket permission errors on Cygwin ([4584](https://github.com/jupyter/notebook/pull/4584))
- Implement optional markdown header and footer files ([4043](https://github.com/jupyter/notebook/pull/4043))
- Remove double link when using `custom_display_url` ([5544](https://github.com/jupyter/notebook/pull/5544))
- Respect `cell.is_editable` during find-and-replace ([5545](https://github.com/jupyter/notebook/pull/5545))
- Fix exception causes all over the codebase ([5556](https://github.com/jupyter/notebook/pull/5556)
- Improve login shell heuristics ([5588](https://github.com/jupyter/notebook/pull/5588))
- Added support for `JUPYTER_TOKEN_FILE` ([5587](https://github.com/jupyter/notebook/pull/5587))
- Kill notebook itself when server cull idle kernel ([5593](https://github.com/jupyter/notebook/pull/5593))
- Implement password hashing with bcrypt ([3793](https://github.com/jupyter/notebook/pull/3793))
- Fix broken links ([5600](https://github.com/jupyter/notebook/pull/5600))
- Russian internationalization support ([5571](https://github.com/jupyter/notebook/pull/5571))
- Add a metadata tag to override notebook direction (ltr/rtl) ([5052](https://github.com/jupyter/notebook/pull/5052))
- Paste two images from clipboard in markdown cell ([5598](https://github.com/jupyter/notebook/pull/5598))
- Add keyboard shortcuts to menu dropdowns ([5525](https://github.com/jupyter/notebook/pull/5525))
- Update codemirror to `5.56.0+components1` ([5637](https://github.com/jupyter/notebook/pull/5637))

Thank you to all the contributors:

- Aaron Myatt
- Adam Blake
- Afshin Taylor Darian
- Aman Bansal
- Ben Thayer
- berendjan
- Bruno P. Kinoshita
- bzinberg
- Christophe Cadilhac
- Daiki Katsuragawa
- David Lukes
- Dmitriy Q
- dmpe
- dylanzjy
- dSchurch
- E. M. Bray
- ErwinRussel
- Felix Mönckemeyer
- Grant Nestor
- Jarrad Whitaker
- Jesus Panales Castillo
- Joshua Zeltser
- Karthikeyan Singaravelan
- Kenichi Ito
- Kevin Bates
- Koki Nishihara
- Kris Wilson
- Kyle Kelley
- Laura Merlo
- levinxo
- Luciano Resende
- Luis Cabezon Manchado
- Madhusudhan Srinivasa
- Matthias Geier
- mattn
- Max Klein
- Min RK
- Mingxuan Lin
- Mohammad Mostafa Farzan
- Niko Felger
- Norah Abanumay
- Onno Broekmans
- PierreMB
- pinarkavak
- Ram Rachum
- Reece Hart
- Remi Rampin
- Rohit Sanjay
- Shane Canon
- Simon Li
- Steinar Sturlaugsson
- Steven Silvester
- taohan16
- Thew Dhanat
- Thomas Kluyver
- Toon Baeyens
- Vidar Tonaas Fauske
- Zachary Sailer

## 6.0.3

- Dependency updates to fix startup issues on Windows platform
- Add support for nbconvert 6.x
- Creation of recent tab

Thanks for all the contributors:

- Luciano Resende
- Kevin Bates
- ahangsleben
- Zachary Sailer
- Pallavi Bharadwaj
- Thomas Kluyver
- Min RK
- forest0
- Bibo Hao
- Michal Charemza
- Sergey Shevelev
- Shuichiro MAKIGAKI
- krinsman
- TPartida
- Landen McDonald
- Tres DuBiel

## 6.0.2

- Update JQuery dependency to version 3.4.1 to fix security vulnerability (CVE-2019-11358)
- Update CodeMirror to version 5.48.4 to fix Python formatting issues
- Continue removing obsolete Python 2.x code/dependencies
- Multiple documentation updates

Thanks for all the contributors:

- David Robles
- Jason Grout
- Kerwin Sun
- Kevin Bates
- Kyle Kelley
- Luciano Resende
- Marcus D Sherman
- Sasaki Takeru
- Tom Jarosz
- Vidar Tonaas Fauske
- Wes Turner
- Zachary Sailer

## 6.0.1

- Attempt to re-establish websocket connection to Gateway ([4777](https://github.com/jupyter/notebook/pull/4777))
- Add missing react-dom js to package data ([4772](https://github.com/jupyter/notebook/pull/4772))

Thanks for all the contributors:

- Eunsoo Park
- Min RK

## 6.0

This is the first major release of the Jupyter Notebook since version
5.0 (March 2017).

We encourage users to start trying JupyterLab, which has just announced
it's 1.0 release in preparation for a future transition.

- Remove Python 2.x support in favor of Python 3.5 and higher.
- Multiple accessibility enhancements and bug-fixes.
- Multiple translation enhancements and bug-fixes.
- Remove deprecated ANSI CSS styles.
- Native support to forward requests to Jupyter Gateway(s) (Embedded
    NB2KG).
- Use JavaScript to redirect users to notebook homepage.
- Enhanced SSL/TLS security by using PROTOCOL_TLS which selects the
    highest ssl/tls protocol version available that both the client and
    server support. When PROTOCOL_TLS is not available use
    PROTOCOL_SSLv23.
- Add `?no_track_activity=1` argument to allow API requests. to not be
    registered as activity (e.g. API calls by external activity
    monitors).
- Kernels shutting down due to an idle timeout is no longer considered
    an activity-updating event.
- Further improve compatibility with tornado 6 with improved checks
    for when websockets are closed.
- Launch the browser with a local file which redirects to the server
    address including the authentication token. This prevents another
    logged-in user from stealing the token from command line arguments
    and authenticating to the server. The single-use token previously
    used to mitigate this has been removed. Thanks to Dr. Owain Kenway
    for suggesting the local file approach.
- Respect nbconvert entrypoints as sources for exporters
- Update to CodeMirror to 5.37, which includes f-string syntax for
    Python 3.6.
- Update jquery-ui to 1.12
- Execute cells by clicking icon in input prompt.
- New "Save as" menu option.
- When serving on a loopback interface, protect against DNS rebinding
    by checking the `Host` header from the browser. This check can be
    disabled if necessary by setting `NotebookApp.allow_remote_access`. (Disabled by default while we work out some Mac issues in
    [3754](https://github.com/jupyter/notebook/issues/3754)).
- Add kernel_info_timeout traitlet to enable restarting slow kernels.
- Add `custom_display_host` config option to override displayed URL.
- Add /metrics endpoint for Prometheus Metrics.
- Optimize large file uploads.
- Allow access control headers to be overriden in
    jupyter_notebook_config.py to support greater CORS and proxy
    configuration flexibility.
- Add support for terminals on windows.
- Add a "restart and run all" button to the toolbar.
- Frontend/extension-config: allow default json files in a .d
    directory.
- Allow setting token via jupyter_token env.
- Cull idle kernels using `--MappingKernelManager.cull_idle_timeout`.
- Allow read-only notebooks to be trusted.
- Convert JS tests to Selenium.

Security Fixes included in previous minor releases of Jupyter Notebook
and also included in version 6.0.

- Fix Open Redirect vulnerability (CVE-2019-10255) where certain
    malicious URLs could redirect from the Jupyter login page to a
    malicious site after a successful login.
- Contains a security fix for a cross-site inclusion (XSSI)
    vulnerability (CVE-2019--9644), where files at a known URL could be
    included in a page from an unauthorized website if the user is
    logged into a Jupyter server. The fix involves setting the
    `X-Content-Type-Options: nosniff` header, and applying CSRF checks
    previously on all non-GET API requests to GET requests to API
    endpoints and the /files/ endpoint.
- Check Host header to more securely protect localhost deployments
    from DNS rebinding. This is a pre-emptive measure, not fixing a
    known vulnerability. Use `.NotebookApp.allow_remote_access` and
    `.NotebookApp.local_hostnames` to configure access.
- Upgrade bootstrap to 3.4, fixing an XSS vulnerability, which has
    been assigned
    [CVE-2018-14041](https://nvd.nist.gov/vuln/detail/CVE-2018-14041).
- Contains a security fix preventing malicious directory names from
    being able to execute javascript.
- Contains a security fix preventing nbconvert endpoints from
    executing javascript with access to the server API. CVE request
    pending.

Thanks for all the contributors:

- AAYUSH SINHA
- Aaron Hall, MBA
- Abhinav Sagar
- Adam Rule
- Adeel Ahmad
- Alex Rothberg
- Amy Skerry-Ryan
- Anastasis Germanidis
- Andrés Sánchez
- Arjun Radhakrishna
- Arovit Narula
- Benda Xu
- Björn Grüning
- Brian E. Granger
- Carol Willing
- Celina Kilcrease
- Chris Holdgraf
- Chris Miller
- Ciaran Langton
- Damian Avila
- Dana Lee
- Daniel Farrell
- Daniel Nicolai
- Darío Hereñú
- Dave Aitken
- Dave Foster
- Dave Hirschfeld
- Denis Ledoux
- Dmitry Mikushin
- Dominic Kuang
- Douglas Hanley
- Elliott Sales de Andrade
- Emilio Talamante Lugo
- Eric Perry
- Ethan T. Hendrix
- Evan Van Dam
- Francesco Franchina
- Frédéric Chapoton
- Félix-Antoine Fortin
- Gabriel
- Gabriel Nützi
- Gabriel Ruiz
- Gestalt LUR
- Grant Nestor
- Gustavo Efeiche
- Harsh Vardhan
- Heng GAO
- Hisham Elsheshtawy
- Hong Xu
- Ian Rose
- Ivan Ogasawara
- J Forde
- Jason Grout
- Jessica B. Hamrick
- Jiaqi Liu
- John Emmons
- Josh Barnes
- Karthik Balakrishnan
- Kevin Bates
- Kirit Thadaka
- Kristian Gregorius Hustad
- Kyle Kelley
- Leo Gallucci
- Lilian Besson
- Lucas Seiki Oshiro
- Luciano Resende
- Luis Angel Rodriguez Guerrero
- M Pacer
- Maarten Breddels
- Mac Knight
- Madicken Munk
- Maitiú Ó Ciaráin
- Marc Udoff
- Mathis HAMMEL
- Mathis Rosenhauer
- Matthias Bussonnier
- Matthias Geier
- Max Vovshin
- Maxime Mouchet
- Michael Chirico
- Michael Droettboom
- Michael Heilman
- Michael Scott Cuthbert
- Michal Charemza
- Mike Boyle
- Milos Miljkovic
- Min RK
- Miro Hrončok
- Nicholas Bollweg
- Nitesh Sawant
- Ondrej Jariabka
- Park Hae Jin
- Paul Ivanov
- Paul Masson
- Peter Parente
- Pierre Tholoniat
- Remco Verhoef
- Roland Weber
- Roman Kornev
- Rosa Swaby
- Roy Hyunjin Han
- Sally
- Sam Lau
- Samar Sultan
- Shiti Saxena
- Simon Biggs
- Spencer Park
- Stephen Ward
- Steve (Gadget) Barnes
- Steven Silvester
- Surya Prakash Susarla
- Syed Shah
- Sylvain Corlay
- Thomas Aarholt
- Thomas Kluyver
- Tim
- Tim Head
- Tim Klever
- Tim Metzler
- Todd
- Tom Jorquera
- Tyler Makaro
- Vaibhav Sagar
- Victor
- Vidar Tonaas Fauske
- Vu Minh Tam
- Vít Tuček
- Will Costello
- Will Starms
- William Hosford
- Xiaohan Li
- Yuvi Panda
- ashley teoh
- nullptr

## 5.7.8

- Fix regression in restarting kernels in 5.7.5. The restart handler
    would return before restart was completed.
- Further improve compatibility with tornado 6 with improved checks
    for when websockets are closed.
- Fix regression in 5.7.6 on Windows where .js files could have the
    wrong mime-type.
- Fix Open Redirect vulnerability (CVE-2019-10255) where certain
    malicious URLs could redirect from the Jupyter login page to a
    malicious site after a successful login. 5.7.7 contained only a
    partial fix for this issue.

## 5.7.6

5.7.6 contains a security fix for a cross-site inclusion (XSSI)
vulnerability (CVE-2019--9644), where files at a known URL could be
included in a page from an unauthorized website if the user is logged
into a Jupyter server. The fix involves setting the
`X-Content-Type-Options: nosniff` header, and applying CSRF checks
previously on all non-GET API requests to GET requests to API endpoints
and the /files/ endpoint.

The attacking page is able to access some contents of files when using
Internet Explorer through script errors, but this has not been
demonstrated with other browsers.

## 5.7.5

- Fix compatibility with tornado 6 ([4392](https://github.com/jupyter/notebook/pull/4392), [4449](https://github.com/jupyter/notebook/pull/4449)).
- Fix opening integer filedescriptor during startup on Python 2 ([4349](https://github.com/jupyter/notebook/pull/4349))
- Fix compatibility with asynchronous
    [KernelManager.restart_kernel]{.title-ref} methods ([4412](https://github.com/jupyter/notebook/pull/4412))

## 5.7.4

5.7.4 fixes a bug introduced in 5.7.3, in which the
`list_running_servers()` function attempts to parse HTML files as JSON,
and consequently crashes ([4284](https://github.com/jupyter/notebook/pull/4284)).

## 5.7.3

5.7.3 contains one security improvement and one security fix:

- Launch the browser with a local file which redirects to the server
    address including the authentication token ([4260](https://github.com/jupyter/notebook/pull/4260)). This prevents another logged-in user from stealing
    the token from command line arguments and authenticating to the
    server. The single-use token previously used to mitigate this has
    been removed. Thanks to Dr. Owain Kenway for suggesting the local
    file approach.
- Upgrade bootstrap to 3.4, fixing an XSS vulnerability, which has
    been assigned
    [CVE-2018-14041](https://nvd.nist.gov/vuln/detail/CVE-2018-14041) ([4271](https://github.com/jupyter/notebook/pull/4271)).

## 5.7.2

5.7.2 contains a security fix preventing malicious directory names from
being able to execute javascript. CVE request pending.

## 5.7.1

5.7.1 contains a security fix preventing nbconvert endpoints from
executing javascript with access to the server API. CVE request pending.

## 5.7.0

New features:

- Update to CodeMirror to 5.37, which includes f-string syntax for
    Python 3.6 ([3816](https://github.com/jupyter/notebook/pull/3816))
- Update jquery-ui to 1.12 ([3836](https://github.com/jupyter/notebook/pull/3836))
- Check Host header to more securely protect localhost deployments
    from DNS rebinding. This is a pre-emptive measure, not fixing a
    known vulnerability ([3766](https://github.com/jupyter/notebook/pull/3766)). Use
    `.NotebookApp.allow_remote_access` and
    `.NotebookApp.local_hostnames` to configure access.
- Allow access-control-allow-headers to be overridden ([3886](https://github.com/jupyter/notebook/pull/3886))
- Allow configuring max_body_size and max_buffer_size ([3829](https://github.com/jupyter/notebook/pull/3829))
- Allow configuring get_secure_cookie keyword-args ([3778](https://github.com/jupyter/notebook/pull/3778))
- Respect nbconvert entrypoints as sources for exporters ([3879](https://github.com/jupyter/notebook/pull/3879))
- Include translation sources in source distributions ([3925](https://github.com/jupyter/notebook/pull/3925), [3931](https://github.com/jupyter/notebook/pull/3931))
- Various improvements to documentation ([3799](https://github.com/jupyter/notebook/pull/3799), [3800](https://github.com/jupyter/notebook/pull/3800),
    [3806](https://github.com/jupyter/notebook/pull/3806), [3883](https://github.com/jupyter/notebook/pull/3883), [3908](https://github.com/jupyter/notebook/pull/3908))

Fixing problems:

- Fix breadcrumb link when running with a base url ([3905](https://github.com/jupyter/notebook/pull/3905))
- Fix possible type error when closing activity stream ([3907](https://github.com/jupyter/notebook/pull/3907))
- Disable metadata editing for non-editable cells ([3744](https://github.com/jupyter/notebook/pull/3744))
- Fix some styling and alignment of prompts caused by regressions in
    5.6.0.
- Enter causing page reload in shortcuts editor ([3871](https://github.com/jupyter/notebook/pull/3871))
- Fix uploading to the same file twice ([3712](https://github.com/jupyter/notebook/pull/3712))

See the 5.7 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.7)
involved in this release.

Thanks to the following contributors:

- Aaron Hall
- Benjamin Ragan-Kelley
- Bill Major
- bxy007
- Dave Aitken
- Denis Ledoux
- Félix-Antoine Fortin
- Gabriel
- Grant Nestor
- Kevin Bates
- Kristian Gregorius Hustad
- M Pacer
- Madicken Munk
- Maitiu O Ciarain
- Matthias Bussonnier
- Michael Boyle
- Michael Chirico
- Mokkapati, Praneet(ES)
- Peter Parente
- Sally Wilsak
- Steven Silvester
- Thomas Kluyver
- Walter Martin

## 5.6.0

New features:

- Execute cells by clicking icon in input prompt ([3535](https://github.com/jupyter/notebook/pull/3535), [3687](https://github.com/jupyter/notebook/pull/3687))
- New "Save as" menu option ([3289](https://github.com/jupyter/notebook/pull/3289))
- When serving on a loopback interface, protect against DNS rebinding
    by checking the `Host` header from the browser ([3714](https://github.com/jupyter/notebook/pull/3714)). This check can be
    disabled if necessary by setting `NotebookApp.allow_remote_access`. (Disabled by default while we work out some Mac issues in
    [3754](https://github.com/jupyter/notebook/issues/3754)).
- Add kernel_info_timeout traitlet to enable restarting slow kernels ([3665](https://github.com/jupyter/notebook/pull/3665))
- Add `custom_display_host` config option to override displayed URL ([3668](https://github.com/jupyter/notebook/pull/3668))
- Add /metrics endpoint for Prometheus Metrics ([3490](https://github.com/jupyter/notebook/pull/3490))
- Update to MathJax 2.7.4 ([3751](https://github.com/jupyter/notebook/pull/3751))
- Update to jQuery 3.3 ([3655](https://github.com/jupyter/notebook/pull/3655))
- Update marked to 0.4 ([3686](https://github.com/jupyter/notebook/pull/3686))

Fixing problems:

- Don't duplicate token in displayed URL ([3656](https://github.com/jupyter/notebook/pull/3656))
- Clarify displayed URL when listening on all interfaces ([3703](https://github.com/jupyter/notebook/pull/3703))
- Don't trash non-empty directories on Windows ([3673](https://github.com/jupyter/notebook/pull/3673))
- Include LICENSE file in wheels ([3671](https://github.com/jupyter/notebook/pull/3671))
- Don't show "0 active kernels" when starting the notebook ([3696](https://github.com/jupyter/notebook/pull/3696))

Testing:

- Add find replace test ([3630](https://github.com/jupyter/notebook/pull/3630))
- Selenium test for deleting all cells ([3601](https://github.com/jupyter/notebook/pull/3601))
- Make creating a new notebook more robust ([3726](https://github.com/jupyter/notebook/pull/3726))

Thanks to the following contributors:

- Arovit Narula ([arovit](https://github.com/arovit))
- lucasoshiro ([lucasoshiro](https://github.com/lucasoshiro))
- M Pacer ([mpacer](https://github.com/mpacer))
- Thomas Kluyver ([takluyver](https://github.com/takluyver))
- Todd ([toddrme2178](https://github.com/toddrme2178))
- Yuvi Panda ([yuvipanda](https://github.com/yuvipanda))

See the 5.6 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.6)
involved in this release.

## 5.5.0

New features:

- The files list now shows file sizes ([3539](https://github.com/jupyter/notebook/pull/3539))
- Add a quit button in the dashboard ([3004](https://github.com/jupyter/notebook/pull/3004))
- Display hostname in the terminal when running remotely ([3356](https://github.com/jupyter/notebook/pull/3356), [3593](https://github.com/jupyter/notebook/pull/3593))
- Add slides exportation/download to the menu ([3287](https://github.com/jupyter/notebook/pull/3287))
- Add any extra installed nbconvert exporters to the "Download as"
    menu ([3323](https://github.com/jupyter/notebook/pull/3323))
- Editor: warning when overwriting a file that is modified on disk ([2783](https://github.com/jupyter/notebook/pull/2783))
- Display a warning message if cookies are not enabled ([3511](https://github.com/jupyter/notebook/pull/3511))
- Basic `__version__` reporting for extensions ([3541](https://github.com/jupyter/notebook/pull/3541))
- Add `NotebookApp.terminals_enabled` config option ([3478](https://github.com/jupyter/notebook/pull/3478))
- Make buffer time between last modified on disk and last modified on
    last save configurable ([3273](https://github.com/jupyter/notebook/pull/3273))
- Allow binding custom shortcuts for 'close and halt' ([3314](https://github.com/jupyter/notebook/pull/3314))
- Add description for 'Trusted' notification ([3386](https://github.com/jupyter/notebook/pull/3386))
- Add `settings['activity_sources']` ([3401](https://github.com/jupyter/notebook/pull/3401))
- Add an `output_updated.OutputArea` event ([3560](https://github.com/jupyter/notebook/pull/3560))

Fixing problems:

- Fixes to improve web accessibility ([3507](https://github.com/jupyter/notebook/pull/3507))
- Fixed color contrast issue in tree.less ([3336](https://github.com/jupyter/notebook/pull/3336))
- Allow cancelling upload of large files ([3373](https://github.com/jupyter/notebook/pull/3373))
- Don't clear login cookie on requests without cookie ([3380](https://github.com/jupyter/notebook/pull/3380))
- Don't trash files on different device to home dir on Linux ([3304](https://github.com/jupyter/notebook/pull/3304))
- Clear waiting asterisks when restarting kernel ([3494](https://github.com/jupyter/notebook/pull/3494))
- Fix output prompt when `execution_count` missing ([3236](https://github.com/jupyter/notebook/pull/3236))
- Make the 'changed on disk' dialog work when displayed twice ([3589](https://github.com/jupyter/notebook/pull/3589))
- Fix going back to root directory with history in notebook list ([3411](https://github.com/jupyter/notebook/pull/3411))
- Allow defining keyboard shortcuts for missing actions ([3561](https://github.com/jupyter/notebook/pull/3561))
- Prevent default on pageup/pagedown when completer is active ([3500](https://github.com/jupyter/notebook/pull/3500))
- Prevent default event handling on new terminal ([3497](https://github.com/jupyter/notebook/pull/3497))
- ConfigManager should not write out default values found in the .d
    directory ([3485](https://github.com/jupyter/notebook/pull/3485))
- Fix leak of iopub object in activity monitoring ([3424](https://github.com/jupyter/notebook/pull/3424))
- Javascript lint in notebooklist.js ([3409](https://github.com/jupyter/notebook/pull/3409))
- Some Javascript syntax fixes ([3294](https://github.com/jupyter/notebook/pull/3294))
- Convert native for loop to `Array.forEach()` ([3477](https://github.com/jupyter/notebook/pull/3477))
- Disable cache when downloading nbconvert output ([3484](https://github.com/jupyter/notebook/pull/3484))
- Add missing digestmod arg to HMAC ([3399](https://github.com/jupyter/notebook/pull/3399))
- Log OSErrors failing to create less-critical files during startup ([3384](https://github.com/jupyter/notebook/pull/3384))
- Use powershell on Windows ([3379](https://github.com/jupyter/notebook/pull/3379))
- API spec improvements, API handler improvements ([3368](https://github.com/jupyter/notebook/pull/3368))
- Set notebook to dirty state after change to kernel metadata ([3350](https://github.com/jupyter/notebook/pull/3350))
- Use CSP header to treat served files as belonging to a separate
    origin ([3341](https://github.com/jupyter/notebook/pull/3341))
- Don't install gettext into builtins ([3330](https://github.com/jupyter/notebook/pull/3330))
- Add missing `import _` ([3316](https://github.com/jupyter/notebook/pull/3316),
    [3326](https://github.com/jupyter/notebook/pull/3326))
- Write `notebook.json` file atomically ([3305](https://github.com/jupyter/notebook/pull/3305))
- Fix clicking with modifiers, page title updates ([3282](https://github.com/jupyter/notebook/pull/3282))
- Upgrade jQuery to version 2.2 ([3428](https://github.com/jupyter/notebook/pull/3428))
- Upgrade xterm.js to 3.1.0 ([3189](https://github.com/jupyter/notebook/pull/3189))
- Upgrade moment.js to 2.19.3 ([3562](https://github.com/jupyter/notebook/pull/3562))
- Upgrade CodeMirror to 5.35 ([3372](https://github.com/jupyter/notebook/pull/3372))
- "Require" pyzmq\>=17 ([3586](https://github.com/jupyter/notebook/pull/3586))

Documentation:

- Documentation updates and organisation ([3584](https://github.com/jupyter/notebook/pull/3584))
- Add section in docs about privacy ([3571](https://github.com/jupyter/notebook/pull/3571))
- Add explanation on how to change the type of a cell to Markdown ([3377](https://github.com/jupyter/notebook/pull/3377))
- Update docs with confd implementation details ([3520](https://github.com/jupyter/notebook/pull/3520))
- Add more information for where `jupyter_notebook_config.py` is
    located ([3346](https://github.com/jupyter/notebook/pull/3346))
- Document options to enable nbextensions in specific sections ([3525](https://github.com/jupyter/notebook/pull/3525))
- jQuery attribute selector value MUST be surrounded by quotes ([3527](https://github.com/jupyter/notebook/pull/3527))
- Do not execute special notebooks with nbsphinx ([3360](https://github.com/jupyter/notebook/pull/3360))
- Other minor fixes in [3288](https://github.com/jupyter/notebook/pull/3288),
    [3528](https://github.com/jupyter/notebook/pull/3528), [3293](https://github.com/jupyter/notebook/pull/3293), [3367](https://github.com/jupyter/notebook/pull/3367)

Testing:

- Testing with Selenium & Sauce labs ([3321](https://github.com/jupyter/notebook/pull/3321))
- Selenium utils + markdown rendering tests ([3458](https://github.com/jupyter/notebook/pull/3458))
- Convert insert cell tests to Selenium ([3508](https://github.com/jupyter/notebook/pull/3508))
- Convert prompt numbers tests to Selenium ([3554](https://github.com/jupyter/notebook/pull/3554))
- Convert delete cells tests to Selenium ([3465](https://github.com/jupyter/notebook/pull/3465))
- Convert undelete cell tests to Selenium ([3475](https://github.com/jupyter/notebook/pull/3475))
- More selenium testing utilities ([3412](https://github.com/jupyter/notebook/pull/3412))
- Only check links when build is trigger by Travis Cron job ([3493](https://github.com/jupyter/notebook/pull/3493))
- Fix Appveyor build errors ([3430](https://github.com/jupyter/notebook/pull/3430))
- Undo patches in teardown before attempting to delete files ([3459](https://github.com/jupyter/notebook/pull/3459))
- Get tests running with tornado 5 ([3398](https://github.com/jupyter/notebook/pull/3398))
- Unpin ipykernel version on Travis ([3223](https://github.com/jupyter/notebook/pull/3223))

Thanks to the following contributors:

- Arovit Narula ([arovit](https://github.com/arovit))
- Ashley Teoh ([ashleytqy](https://github.com/ashleytqy))
- Nicholas Bollweg ([bollwyvl](https://github.com/bollwyvl))
- Alex Rothberg ([cancan101](https://github.com/cancan101))
- Celina Kilcrease ([ckilcrease](https://github.com/ckilcrease))
- dabuside ([dabuside](https://github.com/dabuside))
- Damian Avila ([damianavila](https://github.com/damianavila))
- Dana Lee ([danagilliann](https://github.com/danagilliann))
- Dave Hirschfeld ([dhirschfeld](https://github.com/dhirschfeld))
- Heng GAO ([ehengao](https://github.com/ehengao))
- Leo Gallucci ([elgalu](https://github.com/elgalu))
- Evan Van Dam ([evandam](https://github.com/evandam))
- forbxy ([forbxy](https://github.com/forbxy))
- Grant Nestor ([gnestor](https://github.com/gnestor))
- Ethan T. Hendrix ([hendrixet](https://github.com/hendrixet))
- Miro Hrončok ([hroncok](https://github.com/hroncok))
- Paul Ivanov ([ivanov](https://github.com/ivanov))
- Darío Hereñú ([kant](https://github.com/kant))
- Kevin Bates ([kevin-bates](https://github.com/kevin-bates))
- Maarten Breddels ([maartenbreddels](https://github.com/maartenbreddels))
- Michael Droettboom ([mdboom](https://github.com/mdboom))
- Min RK ([minrk](https://github.com/minrk))
- M Pacer ([mpacer](https://github.com/mpacer))
- Peter Parente ([parente](https://github.com/parente))
- Paul Masson ([paulmasson](https://github.com/paulmasson))
- Philipp Rudiger ([philippjfr](https://github.com/philippjfr))
- Mac Knight ([Shels1909](https://github.com/Shels1909))
- Hisham Elsheshtawy ([Sheshtawy](https://github.com/Sheshtawy))
- Simon Biggs ([SimonBiggs](https://github.com/SimonBiggs))
- Sunil Hari (`@sunilhari`)
- Thomas Kluyver ([takluyver](https://github.com/takluyver))
- Tim Klever ([tklever](https://github.com/tklever))
- Gabriel Ruiz ([unnamedplay-r](https://github.com/unnamedplay-r))
- Vaibhav Sagar ([vaibhavsagar](https://github.com/vaibhavsagar))
- William Hosford ([whosford](https://github.com/whosford))
- Hong ([xuhdev](https://github.com/xuhdev))

See the 5.5 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.5)
involved in this release.

## 5.4.1

A security release to fix [CVE-2018-8768](http://cve.mitre.org/cgi-bin/cvename.cgi?name=CVE-2018-8768).

Thanks to [Alex](https://hackerone.com/pisarenko) for identifying this
bug, and Jonathan Kamens and Scott Sanderson at Quantopian for verifying
it and bringing it to our attention.

## 5.4.0

- Fix creating files and folders after navigating directories in the
    dashboard ([3264](https://github.com/jupyter/notebook/pull/3264)).
- Enable printing notebooks in colour, removing the CSS that made
    everything black and white ([3212](https://github.com/jupyter/notebook/pull/3212)).
- Limit the completion options displayed in the notebook to 1000, to
    avoid performance issues with very long lists ([3195](https://github.com/jupyter/notebook/pull/3195)).
- Accessibility improvements in `tree.html` ([3271](https://github.com/jupyter/notebook/pull/3271)).
- Added alt-text to the kernel logo image in the notebook UI ([3228](https://github.com/jupyter/notebook/pull/3228)).
- Added a test on Travis CI to flag if symlinks are accidentally
    introduced in the future. This should prevent the issue that
    necessitated `release-5.3.1`{.interpreted-text role="ref"} ([3227](https://github.com/jupyter/notebook/pull/3227)).
- Use lowercase letters for random IDs generated in our Javascript ([3264](https://github.com/jupyter/notebook/pull/3264)).
- Removed duplicate code setting `TextCell.notebook` ([3256](https://github.com/jupyter/notebook/pull/3256)).

Thanks to the following contributors:

- Alex Soderman ([asoderman](https://github.com/asoderman))
- Matthias Bussonnier ([Carreau](https://github.com/Carreau))
- Min RK ([minrk](https://github.com/minrk))
- Nitesh Sawant ([ns23](https://github.com/ns23))
- Thomas Kluyver ([takluyver](https://github.com/takluyver))
- Yuvi Panda ([yuvipanda](https://github.com/yuvipanda))

See the 5.4 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.4)
involved in this release.

## 5.3.1

Replaced a symlink in the repository with a copy, to fix issues
installing on Windows ([3220](https://github.com/jupyter/notebook/pull/3220)).

## 5.3.0

This release introduces a couple noteable improvements, such as terminal
support for Windows and support for OS trash (files deleted from the
notebook dashboard are moved to the OS trash vs. deleted permanently).

- Add support for terminals on windows ([3087](https://github.com/jupyter/notebook/pull/3087)).
- Add a "restart and run all" button to the toolbar ([2965](https://github.com/jupyter/notebook/pull/2965)).
- Send files to os trash mechanism on delete ([1968](https://github.com/jupyter/notebook/pull/1968)).
- Allow programmatic copy to clipboard ([3088](https://github.com/jupyter/notebook/pull/3088)).
- Use DOM History API for navigating between directories in the file
    browser ([3115](https://github.com/jupyter/notebook/pull/3115)).
- Add translated files to folder(docs-translations) ([3065](https://github.com/jupyter/notebook/pull/3065)).
- Allow non empty dirs to be deleted ([3108](https://github.com/jupyter/notebook/pull/3108)).
- Set cookie on base_url ([2959](https://github.com/jupyter/notebook/pull/2959)).
- Allow token-authenticated requests cross-origin by default ([2920](https://github.com/jupyter/notebook/pull/2920)).
- Change cull_idle_timeout_minimum to 1 from 300 ([2910](https://github.com/jupyter/notebook/pull/2910)).
- Config option to shut down server after n seconds with no kernels ([2963](https://github.com/jupyter/notebook/pull/2963)).
- Display a "close" button on load notebook error ([3176](https://github.com/jupyter/notebook/pull/3176)).
- Add action to command pallette to run CodeMirror's "indentAuto"
    on selection ([3175](https://github.com/jupyter/notebook/pull/3175)).
- Add option to specify extra services ([3158](https://github.com/jupyter/notebook/pull/3158)).
- Warn_bad_name should not use global name ([3160](https://github.com/jupyter/notebook/pull/3160)).
- Avoid overflow of hidden form ([3148](https://github.com/jupyter/notebook/pull/3148)).
- Fix shutdown trans loss ([3147](https://github.com/jupyter/notebook/pull/3147)).
- Find available kernelspecs more efficiently ([3136](https://github.com/jupyter/notebook/pull/3136)).
- Don't try to translate missing help strings ([3122](https://github.com/jupyter/notebook/pull/3122)).
- Frontend/extension-config: allow default json files in a .d
    directory ([3116](https://github.com/jupyter/notebook/pull/3116)).
- Use [requirejs]{.title-ref} vs. [require]{.title-ref} ([3097](https://github.com/jupyter/notebook/pull/3097)).
- Fixes some ui bugs in firefox \#3044 ([3058](https://github.com/jupyter/notebook/pull/3058)).
- Compare non-specific language code when choosing to use arabic
    numerals ([3055](https://github.com/jupyter/notebook/pull/3055)).
- Fix save-script deprecation ([3053](https://github.com/jupyter/notebook/pull/3053)).
- Include moment locales in package_data ([3051](https://github.com/jupyter/notebook/pull/3051)).
- Fix moment locale loading in bidi support ([3048](https://github.com/jupyter/notebook/pull/3048)).
- Tornado 5: periodiccallback loop arg will be removed ([3034](https://github.com/jupyter/notebook/pull/3034)).
- Use [/files]{.title-ref} prefix for pdf-like files ([3031](https://github.com/jupyter/notebook/pull/3031)).
- Add folder for document translation ([3022](https://github.com/jupyter/notebook/pull/3022)).
- When login-in via token, let a chance for user to set the password ([3008](https://github.com/jupyter/notebook/pull/3008)).
- Switch to jupyter_core implementation of ensure_dir_exists ([3002](https://github.com/jupyter/notebook/pull/3002)).
- Send http shutdown request on 'stop' subcommand ([3000](https://github.com/jupyter/notebook/pull/3000)).
- Work on loading ui translations ([2969](https://github.com/jupyter/notebook/pull/2969)).
- Fix ansi inverse ([2967](https://github.com/jupyter/notebook/pull/2967)).
- Add send2trash to requirements for building docs ([2964](https://github.com/jupyter/notebook/pull/2964)).
- I18n readme.md improvement ([2962](https://github.com/jupyter/notebook/pull/2962)).
- Add 'reason' field to json error responses ([2958](https://github.com/jupyter/notebook/pull/2958)).
- Add some padding for stream outputs ([3194](https://github.com/jupyter/notebook/pull/3194)).
- Always use setuptools in `setup.py` ([3206](https://github.com/jupyter/notebook/pull/3206)).
- Fix clearing cookies on logout when `base_url` is configured ([3207](https://github.com/jupyter/notebook/pull/3207)).

Thanks to the following contributors:

- bacboc ([bacboc](https://github.com/bacboc))
- Steven Silvester ([blink1073](https://github.com/blink1073))
- Matthias Bussonnier ([Carreau](https://github.com/Carreau))
- ChungJooHo ([ChungJooHo](https://github.com/ChungJooHo))
- edida ([edida](https://github.com/edida))
- Francesco Franchina (`ferdas`)
- forbxy ([forbxy](https://github.com/forbxy))
- Grant Nestor ([gnestor](https://github.com/gnestor))
- Josh Barnes ([jcb91](https://github.com/jcb91))
- JocelynDelalande ([JocelynDelalande](https://github.com/JocelynDelalande))
- Karthik Balakrishnan ([karthikb351](https://github.com/karthikb351))
- Kevin Bates ([kevin-bates](https://github.com/kevin-bates))
- Kirit Thadaka ([kirit93](https://github.com/kirit93))
- Lilian Besson ([Naereen](https://github.com/Naereen))
- Maarten Breddels ([maartenbreddels](https://github.com/maartenbreddels))
- Madhu94 ([Madhu94](https://github.com/Madhu94))
- Matthias Geier ([mgeier](https://github.com/mgeier))
- Michael Heilman ([mheilman](https://github.com/mheilman))
- Min RK ([minrk](https://github.com/minrk))
- PHaeJin ([PHaeJin](https://github.com/PHaeJin))
- Sukneet ([Sukneet](https://github.com/Sukneet))
- Thomas Kluyver ([takluyver](https://github.com/takluyver))

See the 5.3 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.3)
involved in this release.

## 5.2.1

- Fix invisible CodeMirror cursor at specific browser zoom levels ([2983](https://github.com/jupyter/notebook/pull/2983)).
- Fix nbconvert handler causing broken export to PDF ([2981](https://github.com/jupyter/notebook/pull/2981)).
- Fix the prompt_area argument of the output area constructor. ([2961](https://github.com/jupyter/notebook/pull/2961)).
- Handle a compound extension in new_untitled ([2949](https://github.com/jupyter/notebook/pull/2949)).
- Allow disabling offline message buffering ([2916](https://github.com/jupyter/notebook/pull/2916)).

Thanks to the following contributors:

- Steven Silvester ([blink1073](https://github.com/blink1073))
- Grant Nestor ([gnestor](https://github.com/gnestor))
- Jason Grout ([jasongrout](https://github.com/jasongrout))
- Min RK ([minrk](https://github.com/minrk))
- M Pacer ([mpacer](https://github.com/mpacer))

See the 5.2.1 milestone on GitHub for a complete list of [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.2.1)
involved in this release.

## 5.2.0

- Allow setting token via jupyter_token env ([2921](https://github.com/jupyter/notebook/pull/2921)).
- Fix some errors caused by raising 403 in get_current_user ([2919](https://github.com/jupyter/notebook/pull/2919)).
- Register contents_manager.files_handler_class directly ([2917](https://github.com/jupyter/notebook/pull/2917)).
- Update viewable_extensions ([2913](https://github.com/jupyter/notebook/pull/2913)).
- Show edit shortcuts modal after shortcuts modal is hidden ([2912](https://github.com/jupyter/notebook/pull/2912)).
- Improve edit/view behavior ([2911](https://github.com/jupyter/notebook/pull/2911)).
- The root directory of the notebook server should never be hidden ([2907](https://github.com/jupyter/notebook/pull/2907)).
- Fix notebook require config to match tools/build-main ([2888](https://github.com/jupyter/notebook/pull/2888)).
- Give page constructor default arguments ([2887](https://github.com/jupyter/notebook/pull/2887)).
- Fix codemirror.less to match codemirror's expected padding layout ([2880](https://github.com/jupyter/notebook/pull/2880)).
- Add x-xsrftoken to access-control-allow-headers ([2876](https://github.com/jupyter/notebook/pull/2876)).
- Buffer messages when websocket connection is interrupted ([2871](https://github.com/jupyter/notebook/pull/2871)).
- Load locale dynamically only when not en-us ([2866](https://github.com/jupyter/notebook/pull/2866)).
- Changed key strength to 2048 bits ([2861](https://github.com/jupyter/notebook/pull/2861)).
- Resync jsversion with python version ([2860](https://github.com/jupyter/notebook/pull/2860)).
- Allow copy operation on modified, read-only notebook ([2854](https://github.com/jupyter/notebook/pull/2854)).
- Update error handling on apihandlers ([2853](https://github.com/jupyter/notebook/pull/2853)).
- Test python 3.6 on travis, drop 3.3 ([2852](https://github.com/jupyter/notebook/pull/2852)).
- Avoid base64-literals in image tests ([2851](https://github.com/jupyter/notebook/pull/2851)).
- Upgrade xterm.js to 2.9.2 ([2849](https://github.com/jupyter/notebook/pull/2849)).
- Changed all python variables named file to file_name to not override
    built_in file ([2830](https://github.com/jupyter/notebook/pull/2830)).
- Add more doc tests ([2823](https://github.com/jupyter/notebook/pull/2823)).
- Typos fix ([2815](https://github.com/jupyter/notebook/pull/2815)).
- Rename and update license \[ci skip\] ([2810](https://github.com/jupyter/notebook/pull/2810)).
- Travis builds doc ([2808](https://github.com/jupyter/notebook/pull/2808)).
- Pull request i18n ([2804](https://github.com/jupyter/notebook/pull/2804)).
- Factor out output_prompt_function, as is done with input prompt ([2774](https://github.com/jupyter/notebook/pull/2774)).
- Use rfc5987 encoding for filenames ([2767](https://github.com/jupyter/notebook/pull/2767)).
- Added path to the resources metadata, the same as in
    from_filename(...) in nbconvert.exporters.py ([2753](https://github.com/jupyter/notebook/pull/2753)).
- Make "extrakeys" consistent for notebook and editor ([2745](https://github.com/jupyter/notebook/pull/2745)).
- Bidi support ([2357](https://github.com/jupyter/notebook/pull/2357)).

Special thanks to [samarsultan](https://github.com/samarsultan) and the
Arabic Competence and Globalization Center Team at IBM Egypt for adding
RTL (right-to-left) support to the notebook!

See the 5.2 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A5.2)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.2)
involved in this release.

## 5.1.0

- Preliminary i18n implementation ([2140](https://github.com/jupyter/notebook/pull/2140)).
- Expose URL with auth token in notebook UI ([2666](https://github.com/jupyter/notebook/pull/2666)).
- Fix search background style ([2387](https://github.com/jupyter/notebook/pull/2387)).
- List running notebooks without requiring `--allow-root` ([2421](https://github.com/jupyter/notebook/pull/2421)).
- Allow session of type other than notebook ([2559](https://github.com/jupyter/notebook/pull/2559)).
- Fix search background style ([2387](https://github.com/jupyter/notebook/pull/2387)).
- Fix some Markdown styling issues ([2571](https://github.com/jupyter/notebook/pull/2571)), ([2691](https://github.com/jupyter/notebook/pull/2691)) and ([2534](https://github.com/jupyter/notebook/pull/2534)).
- Remove keymaps that conflict with non-English keyboards ([2535](https://github.com/jupyter/notebook/pull/2535)).
- Add session-specific favicons (notebook, terminal, file) ([2452](https://github.com/jupyter/notebook/pull/2452)).
- Add /api/shutdown handler ([2507](https://github.com/jupyter/notebook/pull/2507)).
- Include metadata when copying a cell ([2349](https://github.com/jupyter/notebook/pull/2349)).
- Stop notebook server from command line ([2388](https://github.com/jupyter/notebook/pull/2388)).
- Improve "View" and "Edit" file handling in dashboard ([2449](https://github.com/jupyter/notebook/pull/2449)) and ([2402](https://github.com/jupyter/notebook/pull/2402)).
- Provide a promise to replace use of the
    `app_initialized.NotebookApp` event ([2710](https://github.com/jupyter/notebook/pull/2710)).
- Fix disabled collapse/expand output button ([2681](https://github.com/jupyter/notebook/pull/2681)).
- Cull idle kernels using `--MappingKernelManager.cull_idle_timeout` ([2215](https://github.com/jupyter/notebook/pull/2215)).
- Allow read-only notebooks to be trusted ([2718](https://github.com/jupyter/notebook/pull/2718)).

See the 5.1 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A5.1)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A5.1)
involved in this release.

## 5.0.0

This is the first major release of the Jupyter Notebook since version
4.0 was created by the "Big Split" of IPython and Jupyter.

We encourage users to start trying JupyterLab in preparation for a
future transition.

We have merged more than 300 pull requests since 4.0. Some of the major
user-facing changes are described here.

### File sorting in the dashboard

Files in the dashboard may now be sorted by last modified date or name
([943](https://github.com/jupyter/notebook/pull/943)):

![image](/_static/images/dashboard-sort.png)

### Cell tags

There is a new cell toolbar for adding *cell tags*
([2048](https://github.com/jupyter/notebook/pull/2048)):

![image](/_static/images/cell-tags-toolbar.png)

Cell tags are a lightweight way to customise the behaviour of tools
working with notebooks; we're working on building support for them into
tools like [nbconvert](https://nbconvert.readthedocs.io/en/latest/) and
[nbval](https://github.com/computationalmodelling/nbval). To start using
tags, select `Tags` in the `View > Cell Toolbar` menu in a notebook.

The UI for editing cell tags is basic for now; we hope to improve it in
future releases.

### Table style

The default styling for tables in the notebook has been updated
([1776](https://github.com/jupyter/notebook/pull/1776)).

Before:

![image](/_static/images/table-style-before.png)

After:

![image](/_static/images/table-style-after.png)

### Customise keyboard shortcuts

You can now edit keyboard shortcuts for *Command Mode* within the UI
([1347](https://github.com/jupyter/notebook/pull/1347)):

![image](/_static/images/shortcut-editor.png)

See the `Help > Edit Keyboard Shortcuts` menu item and follow the
instructions.

### Other additions

- You can copy and paste cells between notebooks, using
    `Ctrl-C`{.interpreted-text role="kbd"} and
    `Ctrl-V`{.interpreted-text role="kbd"} (`Cmd-C`{.interpreted-text
    role="kbd"} and `Cmd-V`{.interpreted-text role="kbd"} on Mac).
- It's easier to configure a password for the notebook with the new
    `jupyter notebook password` command ([2007](https://github.com/jupyter/notebook/pull/2007)).
- The file list can now be ordered by *last modified* or by *name* ([943](https://github.com/jupyter/notebook/pull/943)).
- Markdown cells now support attachments. Simply drag and drop an
    image from your desktop to a markdown cell to add it. Unlike
    relative links that you enter manually, attachments are embedded in
    the notebook itself. An unreferenced attachment will be
    automatically scrubbed from the notebook on save ([621](https://github.com/jupyter/notebook/pull/621)).
- Undoing cell deletion now supports undeleting multiple cells. Cells
    may not be in the same order as before their deletion, depending on
    the actions you did on the meantime, but this should should help
    reduce the impact of accidentally deleting code.
- The file browser now has *Edit* and *View* buttons.
- The file browser now supports moving multiple files at once ([1088](https://github.com/jupyter/notebook/pull/1088)).
- The Notebook will refuse to run as root unless the `--allow-root`
    flag is given ([1115](https://github.com/jupyter/notebook/pull/1115)).
- Keyboard shortcuts are now declarative ([1234](https://github.com/jupyter/notebook/pull/1234)).
- Toggling line numbers can now affect all cells ([1312](https://github.com/jupyter/notebook/pull/1312)).
- Add more visible *Trusted* and *Untrusted* notifications ([1658](https://github.com/jupyter/notebook/pull/1658)).
- The favicon (browser shortcut icon) now changes to indicate when the
    kernel is busy ([1837](https://github.com/jupyter/notebook/pull/1837)).
- Header and toolbar visibility is now persisted in nbconfig and
    across sessions ([1769](https://github.com/jupyter/notebook/pull/1769)).
- Load server extensions with ConfigManager so that merge happens
    recursively, unlike normal config values, to make it load more
    consistently with frontend extensions([2108](https://github.com/jupyter/notebook/pull/2108)).
- The notebook server now supports the [bundler
    API](https://jupyter-notebook.readthedocs.io/en/latest/extending/bundler_extensions.html)
    from the [jupyter_cms incubator
    project](https://github.com/jupyter-incubator/contentmanagement) ([1579](https://github.com/jupyter/notebook/pull/1579)).
- The notebook server now provides information about kernel activity
    in its kernel resource API ([1827](https://github.com/jupyter/notebook/pull/1827)).

Remember that upgrading `notebook` only affects the user interface.
Upgrading kernels and libraries may also provide new features, better
stability and integration with the notebook interface.

## 4.4.0

- Allow override of output callbacks to redirect output messages. This
    is used to implement the ipywidgets Output widget, for example.
- Fix an async bug in message handling by allowing comm message
    handlers to return a promise which halts message processing until
    the promise resolves.

See the 4.4 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A4.4)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A4.4)
involved in this release.

## 4.3.2

4.3.2 is a patch release with a bug fix for CodeMirror and improved
handling of the "editable" cell metadata field.

- Monkey-patch for CodeMirror that resolves
    [\#2037](https://github.com/jupyter/notebook/issues/2037) without
    breaking [\#1967](https://github.com/jupyter/notebook/issues/1967)
- Read-only (`"editable": false`) cells can be executed but cannot be
    split, merged, or deleted

See the 4.3.2 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A4.3.2)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A4.3.2)
involved in this release.

## 4.3.1

4.3.1 is a patch release with a security patch, a couple bug fixes, and
improvements to the newly-released token authentication.

**Security fix**:

- CVE-2016-9971. Fix CSRF vulnerability, where malicious forms could
    create untitled files and start kernels (no remote execution or
    modification of existing files) for users of certain browsers (Firefox, Internet Explorer / Edge). All previous notebook releases
    are affected.

Bug fixes:

- Fix carriage return handling
- Make the font size more robust against fickle browsers
- Ignore resize events that bubbled up and didn't come from window
- Add Authorization to allowed CORS headers
- Downgrade CodeMirror to 5.16 while we figure out issues in Safari

Other improvements:

- Better docs for token-based authentication
- Further highlight token info in log output when autogenerated

See the 4.3.1 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A4.3.1)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A4.3.1)
involved in this release.

## 4.3.0

4.3 is a minor release with many bug fixes and improvements. The biggest
user-facing change is the addition of token authentication, which is
enabled by default. A token is generated and used when your browser is
opened automatically, so you shouldn't have to enter anything in the
default circumstances. If you see a login page (e.g. by switching
browsers, or launching on a new port with `--no-browser`), you get a
login URL with the token from the command `jupyter notebook list`, which
you can paste into your browser.

Highlights:

- API for creating mime-type based renderer extensions using
    `OutputArea.register_mime_type` and `Notebook.render_cell_output`
    methods. See
    [mimerender-cookiecutter](https://github.com/jupyterlab/mimerender-cookiecutter)
    for reference implementations and cookiecutter.
- Enable token authentication by default. See
    `server_security`{.interpreted-text role="ref"} for more details.
- Update security docs to reflect new signature system
- Switched from term.js to xterm.js

Bug fixes:

- Ensure variable is set if exc_info is falsey
- Catch and log handler exceptions in `events.trigger`
- Add debug log for static file paths
- Don't check origin on token-authenticated requests
- Remove leftover print statement
- Fix highlighting of Python code blocks
- `json_errors` should be outermost decorator on API handlers
- Fix remove old nbserver info files
- Fix notebook mime type on download links
- Fix carriage symbol behavior
- Fix terminal styles
- Update dead links in docs
- If kernel is broken, start a new session
- Include cross-origin check when allowing login URL redirects

Other improvements:

- Allow JSON output data with mime type `application/*+json`
- Allow kernelspecs to have spaces in them for backward compat
- Allow websocket connections from scripts
- Allow `None` for post_save_hook
- Upgrade CodeMirror to 5.21
- Upgrade xterm to 2.1.0
- Docs for using comms
- Set `dirty` flag when output arrives
- Set `ws-url` data attribute when accessing a notebook terminal
- Add base aliases for nbextensions
- Include `@` operator in CodeMirror IPython mode
- Extend mathjax_url docstring
- Load nbextension in predictable order
- Improve the error messages for nbextensions
- Include cross-origin check when allowing login URL redirects

See the 4.3 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?utf8=%E2%9C%93&q=is%3Aissue%20milestone%3A4.3%20)
and [pull
requests](https://github.com/jupyter/notebook/pulls?utf8=%E2%9C%93&q=is%3Apr%20milestone%3A4.3%20)
involved in this release.

## 4.2.3

4.2.3 is a small bugfix release on 4.2.

> Highlights:

- Fix regression in 4.2.2 that delayed loading custom.js until after
    `notebook_loaded` and `app_initialized` events have fired.
- Fix some outdated docs and links.

## 4.2.2

4.2.2 is a small bugfix release on 4.2, with an important security fix.
All users are strongly encouraged to upgrade to 4.2.2.

> Highlights:

- **Security fix**: CVE-2016-6524, where untrusted latex output could
    be added to the page in a way that could execute javascript.
- Fix missing POST in OPTIONS responses.
- Fix for downloading non-ascii filenames.
- Avoid clobbering ssl_options, so that users can specify more
    detailed SSL configuration.
- Fix inverted load order in nbconfig, so user config has highest
    priority.
- Improved error messages here and there.

## 4.2.1

4.2.1 is a small bugfix release on 4.2. Highlights:

- Compatibility fixes for some versions of ipywidgets
- Fix for ignored CSS on Windows
- Fix specifying destination when installing nbextensions

## 4.2.0

Release 4.2 adds a new API for enabling and installing extensions.
Extensions can now be enabled at the system-level, rather than just
per-user. An API is defined for installing directly from a Python
package, as well.

Highlighted changes:

- Upgrade MathJax to 2.6 to fix vertical-bar appearing on some
    equations.
- Restore ability for notebook directory to be root (4.1 regression)
- Large outputs are now throttled, reducing the ability of output
    floods to kill the browser.
- Fix the notebook ignoring cell executions while a kernel is starting
    by queueing the messages.
- Fix handling of url prefixes (e.g. JupyterHub) in terminal and edit
    pages.
- Support nested SVGs in output.

And various other fixes and improvements.

## 4.1.0

Bug fixes:

- Properly reap zombie subprocesses
- Fix cross-origin problems
- Fix double-escaping of the base URL prefix
- Handle invalid unicode filenames more gracefully
- Fix ANSI color-processing
- Send keepalive messages for web terminals
- Fix bugs in the notebook tour

UI changes:

- Moved the cell toolbar selector into the *View* menu. Added a button
    that triggers a "hint" animation to the main toolbar so users can
    find the new location. (Click here to see a
    [screencast](https://cloud.githubusercontent.com/assets/335567/10711889/59665a5a-7a3e-11e5-970f-86b89592880c.gif)
    )

    > ![image](/_static/images/cell-toolbar-41.png)

- Added *Restart & Run All* to the *Kernel* menu. Users can also bind
    it to a keyboard shortcut on action
    `restart-kernel-and-run-all-cells`.

- Added multiple-cell selection. Users press `Shift-Up/Down` or
    `Shift-K/J` to extend selection in command mode. Various actions
    such as cut/copy/paste, execute, and cell type conversions apply to
    all selected cells.

    ![image](/_static/images/multi-select-41.png)

- Added a command palette for executing Jupyter actions by name. Users
    press `Cmd/Ctrl-Shift-P` or click the new command palette icon on
    the toolbar.

    ![image](/_static/images/command-palette-41.png)

- Added a *Find and Replace* dialog to the *Edit* menu. Users can also
    press `F` in command mode to show the dialog.

    ![image](/_static/images/find-replace-41.png)

Other improvements:

- Custom KernelManager methods can be Tornado coroutines, allowing
    async operations.
- Make clearing output optional when rewriting input with
    `set_next_input(replace=True)`.
- Added support for TLS client authentication via
    `--NotebookApp.client-ca`.
- Added tags to `jupyter/notebook` releases on DockerHub. `latest`
    continues to track the master branch.

See the 4.1 milestone on GitHub for a complete list of
[issues](https://github.com/jupyter/notebook/issues?page=3&q=milestone%3A4.1+is%3Aclosed+is%3Aissue&utf8=%E2%9C%93)
and [pull
requests](https://github.com/jupyter/notebook/pulls?q=milestone%3A4.1+is%3Aclosed+is%3Apr)
handled.

## 4.0.x

### 4.0.6

- fix installation of mathjax support files
- fix some double-escape regressions in 4.0.5
- fix a couple of cases where errors could prevent opening a notebook

### 4.0.5

Security fixes for maliciously crafted files.

- [CVE-2015-6938](http://www.openwall.com/lists/oss-security/2015/09/02/3):
    malicious filenames
- [CVE-2015-7337](http://www.openwall.com/lists/oss-security/2015/09/16/3):
    malicious binary files in text editor.

Thanks to Jonathan Kamens at Quantopian and Juan Broullón for the
reports.

### 4.0.4

- Fix inclusion of mathjax-safe extension

### 4.0.2

- Fix launching the notebook on Windows
- Fix the path searched for frontend config

### 4.0.0

First release of the notebook as a standalone package.
# Jupyter Notebook

[![Google Group](https://img.shields.io/badge/-Google%20Group-lightgrey.svg)](https://groups.google.com/forum/#!forum/jupyter)
[![Build Status](https://travis-ci.org/jupyter/notebook.svg?branch=master)](https://travis-ci.org/jupyter/notebook)
[![Documentation Status](https://readthedocs.org/projects/jupyter-notebook/badge/?version=latest)](https://jupyter-notebook.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/jupyter/notebook/branch/master/graph/badge.svg)](https://codecov.io/gh/jupyter/notebook)

The Jupyter notebook is a web-based notebook environment for interactive
computing.

![Jupyter notebook example](docs/resources/running_code_med.png "Jupyter notebook example")

### Notice
Please note that this repository is currently maintained by a skeleton crew of maintainers from the Jupyter community. We encourage users to transition to JupyterLab, where more immediate support can occur. Our approach moving forward will be:

1. To maintain the security of the Jupyter Notebook. That means security-related issues and pull requests are our highest priority.
2. To address JupyterLab [feature parity issues](https://github.com/jupyterlab/jupyterlab/issues?q=is%3Aopen+is%3Aissue+label%3A%22tag%3AFeature+Parity%22). As part of this effort, we are also working on a better [notebook-only experience](https://github.com/jupyterlab/jupyterlab/issues/8450) in JupyterLab for users who prefer the UI of the classic Jupyter Notebook.
3. To be responsive to the hard work of community members who have opened pull requests. We are triaging these PRs. We cannot support or maintain new features at this time, but we welcome security and other sustainability fixes.

If you have an open pull request with a new feature or if you were planning to open one, please consider shipping it as a [notebook extension](https://jupyter-notebook.readthedocs.io/en/stable/extending/) instead.

##### Alternatives to contributing to `notebook`
Additionally, please consider whether your contribution would be appropriate for either the underlying server for Jupyter front-ends, [jupyter_server](https://github.com/jupyter/jupyter_server) or in the [JupyterLab front-end](https://github.com/jupyterlab/jupyterlab).

### Jupyter notebook, the language-agnostic evolution of IPython notebook
Jupyter notebook is a language-agnostic HTML notebook application for
Project Jupyter. In 2015, Jupyter notebook was released as a part of
The Big Split™ of the IPython codebase. IPython 3 was the last major monolithic
release containing both language-agnostic code, such as the *IPython notebook*,
and language specific code, such as the *IPython kernel for Python*. As
computing spans across many languages, Project Jupyter will continue to develop the
language-agnostic **Jupyter notebook** in this repo and with the help of the
community develop language specific kernels which are found in their own
discrete repos.
[[The Big Split™ announcement](https://blog.jupyter.org/the-big-split-9d7b88a031a7)]
[[Jupyter Ascending blog post](https://blog.jupyter.org/jupyter-ascending-1bf5b362d97e)]

## Installation
You can find the installation documentation for the
[Jupyter platform, on ReadTheDocs](https://jupyter.readthedocs.io/en/latest/install.html).
The documentation for advanced usage of Jupyter notebook can be found
[here](https://jupyter-notebook.readthedocs.io/en/latest/).

For a local installation, make sure you have
[pip installed](https://pip.readthedocs.io/en/stable/installing/) and run:

    $ pip install notebook

## Usage - Running Jupyter notebook

### Running in a local installation

Launch with:

    $ jupyter notebook

### Running in a remote installation

You need some configuration before starting Jupyter notebook remotely. See [Running a notebook server](https://jupyter-notebook.readthedocs.io/en/stable/public_server.html).

## Development Installation

See [`CONTRIBUTING.rst`](CONTRIBUTING.rst) for how to set up a local development installation.

## Contributing

If you are interested in contributing to the project, see [`CONTRIBUTING.rst`](CONTRIBUTING.rst).

## Resources
- [Project Jupyter website](https://jupyter.org)
- [Online Demo at jupyter.org/try](https://jupyter.org/try)
- [Documentation for Jupyter notebook](https://jupyter-notebook.readthedocs.io/en/latest/) [[PDF](https://media.readthedocs.org/pdf/jupyter-notebook/latest/jupyter-notebook.pdf)]
- [Korean Version of Installation](https://github.com/ChungJooHo/Jupyter_Kor_doc/)
- [Documentation for Project Jupyter](https://jupyter.readthedocs.io/en/latest/index.html) [[PDF](https://media.readthedocs.org/pdf/jupyter/latest/jupyter.pdf)]
- [Issues](https://github.com/jupyter/notebook/issues)
- [Technical support - Jupyter Google Group](https://groups.google.com/forum/#!forum/jupyter)
# Making a Release of Notebook

## Using `jupyter_releaser`

The recommended way to make a release is to use [`jupyter_releaser`](https://github.com/jupyter-server/jupyter_releaser#checklist-for-adoption).

## Manual Release Process

### Start from a fresh git checkout and conda environment

#### Set the release branch

```bash
export release_branch=master
```

#### Create the git checkout

```bash
git clone git@github.com:jupyter/notebook.git
cd notebook
git checkout ${release_banch}
```

#### Create and activate the conda environment

```bash
conda create -n notebook-release -c conda-forge jupyter
conda activate notebook-release
```

### Perform a local dev install

```bash
pip install -ve .
```

### Install release dependencies

```bash
conda install -c conda-forge nodejs babel
npm install -g po2json
pip install jupyter_releaser  # used for build dependencies (build, twine, tbump)
```

### Update the version

```bash
tbump --only-patch <new_version> # set the new version
python setup.py jsversion
git commit -am "Release $(python setup.py --version)"
git tag $(python setup.py --version)
```

### Create the artifacts

```bash
rm -rf dist
python -m build .
```

### Upload the artifacts

```bash
twine check dist/* && twine upload dist/*
```

### Change back to dev version

```bash
tbump --only-patch <dev_version>  # Add the .dev suffix
python setup.py jsversion
git commit -am "Back to dev version"
```

### Push the commits and tags

```bash
git push origin ${release_branch} --tags
```
# IPython Notebook JavaScript Tests

This directory includes regression tests for the web notebook. These tests
depend on [CasperJS](https://github.com/casperjs/casperjs/), which in turn requires a recent
version of [PhantomJS](http://phantomjs.org/).

The JavaScript tests are organized into subdirectories that match those in
`static` (`base`, `notebook`, `services`, `tree`, etc.).

To run all of the JavaScript tests do:

```
python -m notebook.jstest 
```

To run the JavaScript tests for a specific file (`base/utils.js` in this case)
do:

```
python -m notebook.jstest base/utils.js
```

The file `jstest.py` will automatically launch a notebook server to run the
tests against. You can however specify the url of a running notebook server
by using `--url=http://localhost:8888`.
# Implementation Notes for Internationalization of Jupyter Notebook

The implementation of i18n features for jupyter notebook is still a work-in-progress:

- User interface strings are (mostly) handled
- Console messages are not handled (their usefulness in a translated environment is questionable)
- Tooling has to be refined

However…

## How the language is selected ?

1. `jupyter notebook` command reads the `LANG` environment variable at startup,
(`xx_XX` or just `xx` form, where `xx` is the language code you're wanting to
run in).

Hint: if running Windows, you can set it in PowerShell with `${Env:LANG} = "xx_XX"`.
      if running Ubuntu 14, you should set environment variable `LANGUAGE="xx_XX"`.

2. The preferred language for web pages in your browser settings (`xx`) is
   also used. At the moment, it has to be first in the list.

## Contributing and managing translations

Finding and translating the `.pot` files could be (partially) done with a translation API, see the repo [Jupyter Notebook Azure Translator](https://github.com/berendjan/Jupyter-Notebook-Azure-Translator.git) for a possible starting point. (Not affiliated with Jupyter)

### Requirements

- *pybabel* (could be installed `pip install babel`)
- *po2json* (could be installed with `npm install -g po2json`)

**All i18n-related commands are done from the related directory :**

    cd notebook/i18n/

### Message extraction

The translatable material for notebook is split into 3 `.pot` files, as follows:

- *notebook/i18n/notebook.pot* - Console and startup messages, basically anything that is
	produced by Python code.
- *notebook/i18n/nbui.pot* - User interface strings, as extracted from the Jinja2 templates
	in *notebook/templates/\*.html*
- *noteook/i18n/nbjs.pot* - JavaScript strings and dialogs, which contain much of the visible
	user interface for Jupyter notebook.

To extract the messages from the source code whenever new material is added, use the
`pybabel` command:

```shell
pybabel extract -F babel_notebook.cfg -o notebook.pot --no-wrap --project Jupyter .
pybabel extract -F babel_nbui.cfg -o nbui.pot --no-wrap --project Jupyter .
pybabel extract -F babel_nbjs.cfg -o nbjs.pot --no-wrap --project Jupyter .
```

After this is complete you have 3 `.pot` files that you can give to a translator for your favorite language.

Finding and translating the `.pot` files could be (partially done with a translation API, see the repo [Jupyter Notebook Azure Translator](https://github.com/berendjan/Jupyter-Notebook-Azure-Translator.git) for a possible starting point. (Not affiliated with Jupyter)

### Messages compilation

After the source material has been translated, you should have 3 `.po` files with the same base names
as the `.pot` files above.  Put them in `notebook/i18n/${LANG}/LC_MESSAGES`, where `${LANG}` is the language
code for your desired language ( i.e. German = "de", Japanese = "ja", etc. ).

*notebook.po* and *nbui.po* need to be converted from `.po` to `.mo` format for
use at runtime.

```shell
pybabel compile -D notebook -f -l ${LANG} -i ${LANG}/LC_MESSAGES/notebook.po -o ${LANG}/LC_MESSAGES/notebook.mo
pybabel compile -D nbui -f -l ${LANG} -i ${LANG}/LC_MESSAGES/nbui.po -o ${LANG}/LC_MESSAGES/nbui.mo
```

*nbjs.po* needs to be converted to JSON for use within the JavaScript code, with  *po2json*, as follows:

    po2json -p -F -f jed1.x -d nbjs ${LANG}/LC_MESSAGES/nbjs.po ${LANG}/LC_MESSAGES/nbjs.json

When new languages get added, their language codes should be added to *notebook/i18n/nbjs.json*
under the `supported_languages` element.

### Tips for Jupyter developers

The biggest "mistake" I found while doing i18n enablement was the habit of constructing UI messages
from English "piece parts".  For example, code like:

```javascript
var msg = "Enter a new " + type + "name:"
```

where `type` is either "file", "directory", or "notebook"....

is problematic when doing translations, because the surrounding text may need to vary
depending on the inserted word. In this case, you need to switch it and use complete phrases,
as follows:

```javascript
var rename_msg = function (type) {
    switch(type) {
        case 'file': return _("Enter a new file name:");
        case 'directory': return _("Enter a new directory name:");
        case 'notebook': return _("Enter a new notebook name:");
        default: return _("Enter a new name:");
    }
}
```

Also you need to remember that adding an "s" or "es" to an English word to
create the plural form doesn't translate well.  Some languages have as many as 5 or 6 different
plural forms for differing numbers, so using an API such as ngettext() is necessary in order
to handle these cases properly.

### Known issues and future evolutions

1. Right now there are two different places where the desired language is set.  At startup time, the Jupyter console's messages pay attention to the setting of the `${LANG}` environment variable
as set in the shell at startup time.  Unfortunately, this is also the time where the Jinja2
environment is set up, which means that the template stuff will always come from this setting.
We really want to be paying attention to the browser's settings for the stuff that happens in the
browser, so we need to be able to retrieve this information after the browser is started and somehow
communicate this back to Jinja2.  So far, I haven't yet figured out how to do this, which means that if the ${LANG} at startup doesn't match the browser's settings, you could potentially get a mix
of languages in the UI ( never a good thing ).

2. We will need to decide if console messages should be translatable, and enable them if desired.
3. The keyboard shortcut editor was implemented after the i18n work was completed, so that portion
does not have translation support at this time.
4. Babel's documentation has instructions on how to integrate messages extraction
into your *setup.py* so that eventually we can just do:

    ./setup.py extract_messages

I hope to get this working at some point in the near future.
5. The conversions from `.po` to `.mo` probably can and should be done using `setup.py install`.


Any questions or comments please let me know @JCEmmons on github (emmo@us.ibm.com)
---
name: Is this a bug in Notebook? Open an issue.
about: If you're not sure, feel free to post your question on Jupyter's Discourse channel.
title: ''
labels: ''
assignees: ''

---

<!--
BEFORE YOU OPEN AN ISSUE, PLEASE READ THIS

Hello! Thank you for using Jupyter Notebook. We're glad you're here.

Right now, you're opening an issue. Before you do, let's make sure this is the right place to post your question/issue.

First, it's important to know that Jupyter Notebook development has moved into a phase of maintenance-only. There are very few people with limited time maintaining this repository. This means, we won't likely accept new features here. Instead, we recommend that you check out JupyterLab (https://github.com/jupyterlab/jupyterlab)—Jupyter's next generation Notebook interface.

Here, we're looking for specific bugs in the Jupyter Notebook codebase. If you think you've identified such a bug, you can continue opening your issue here. We'd appreciate if you include as much detail as possible, i.e. links to the offending code, snapshots of the UI issue, code-blocks with your console logs, etc.

If you're having issues installing Jupyter Notebook, or you're having another issue and don't know how to proceed, try the following:

1. scan the "What to do when things go wrong" (https://jupyter-notebook.readthedocs.io/en/stable/troubleshooting.html#what-to-do-when-things-go-wrong) page in our documentation to see if your question has already been answered

2. post your question on the Jupyter Notebook discourse channel (https://discourse.jupyter.org/c/notebook/31). There are many more people in the Jupyter community that engage on that channel.
-->

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
git hooks for Jupyter

add these to your `.git/hooks`

For now, we just have `post-checkout` and `post-merge`,
both of which attempt to rebuild javascript and css sourcemaps,
so make sure that you have a fully synced repo whenever you checkout or pull.

To use these hooks, run `./install-hooks.sh`.  
# Jupyter Notebook

[![Google Group](https://img.shields.io/badge/-Google%20Group-lightgrey.svg)](https://groups.google.com/forum/#!forum/jupyter)
[![Build Status](https://travis-ci.org/jupyter/notebook.svg?branch=master)](https://travis-ci.org/jupyter/notebook)
[![Documentation Status](https://readthedocs.org/projects/jupyter-notebook/badge/?version=latest)](https://jupyter-notebook.readthedocs.io/en/latest/?badge=latest)
                
英語版のリンク : [[English Version](http://github.com/jupyter/notebook/)]

Jupyter Notebookは、インタラクティブなWebベースのノートブック形式の環境です。

![Jupyter notebook example](resources/running_code_med.png "Jupyter notebook example")

### Jupyter Notebook, 言語に依存しないIPython Notebookの進化

Jupyter Notebookは、Project Jupyter用の言語に依存しないHTMLノートブックアプリケーションです。
2015年、Jupyter NotebookはIPythonコードベースのThe Big Split™の一部としてリリースされました。IPython3はIPython Notebookなどのユーザーの言語に依存しないコードとIPython kernel for Pythonのような特定の言語ベースのコードの機能を持ってリリースしました。
コンピューティングは多くの言語にまたがるため、Project Jupyterはこのリポジトリで言語に依存しない**Jupyter Notebook**を継続的に開発します。そして、コミュニティの助けを借りて、独自のリポジトリにある言語固有のカーネルを開発します。
[[The Big Split™ announcement](https://blog.jupyter.org/the-big-split-9d7b88a031a7)]
[[Jupyter Ascending blog post](https://blog.jupyter.org/jupyter-ascending-1bf5b362d97e)]

## インストール

[Jupyter platform, on ReadTheDocs](https://jupyter.readthedocs.io/en/latest/install.html)から、インストールドキュメントをご覧になれます。
Jupyter Notebookの発展的な使用方法に関するドキュメントは、[こちら](https://jupyter-notebook.readthedocs.io/en/latest/)をご覧ください。

ローカルへのインストールの場合、[pip](https://pip.readthedocs.io/en/stable/installing/)をインストールしていることを確認し、以下のコマンドを実行してください。

    $ pip install notebook

## 使用方法 - Jupyter Notebookの実行

### ローカルへのインストールにおける実行

以下のコマンドをを実行してください：

    $ jupyter notebook

### リモートへのインストールにおける実行

Jupyter Notebookをリモートで起動する前に、いくつかの構成が必要です。 [Notebookサーバーの実行](https://jupyter-notebook.readthedocs.io/en/stable/public_server.html)を参照してください。

## 開発用インストール

開発用インストールのセットアップ方法については、[`CONTRIBUTING.rst`](https://github.com/jupyter/notebook/blob/master/CONTRIBUTING.rst)を参照してください。

## 貢献

プロジェクトへの貢献に興味がある場合は、[`CONTRIBUTING.rst`](https://github.com/jupyter/notebook/blob/master/CONTRIBUTING.rst)をご覧ください。

## 参考

- [Project Jupyter website](https://jupyter.org)
- [Online Demo at try.jupyter.org](https://try.jupyter.org)
- [Documentation for Jupyter notebook](https://jupyter-notebook.readthedocs.io/en/latest/) [[PDF](https://media.readthedocs.org/pdf/jupyter-notebook/latest/jupyter-notebook.pdf)]
- [Documentation for Project Jupyter](https://jupyter.readthedocs.io/en/latest/index.html) [[PDF](https://media.readthedocs.org/pdf/jupyter/latest/jupyter.pdf)]
- [Issues](https://github.com/jupyter/notebook/issues)
- [Technical support - Jupyter Google Group](https://groups.google.com/forum/#!forum/jupyter)
# Jupyter Notebook

[![Google Group](https://img.shields.io/badge/-Google%20Group-lightgrey.svg)](https://groups.google.com/forum/#!forum/jupyter)
[![Build Status](https://travis-ci.org/jupyter/notebook.svg?branch=master)](https://travis-ci.org/jupyter/notebook)
[![Documentation Status](https://readthedocs.org/projects/jupyter-notebook/badge/?version=latest)](http://jupyter-notebook.readthedocs.io/en/latest/?badge=latest)
                
English 버전 링크 : [[English Version](http://github.com/jupyter/notebook/)]

Jupyter notebook 은 상호 교환을 위한 웹 기반 환경입니다.

![Jupyter notebook example](resources/running_code_med.png "Jupyter notebook example")

### Jupyter notebook, 사용자의 언어에 독립적인 IPython notebook의 진화
Jupyter notebook은 Jupyter 프로젝트를 위한 사용자 언어에 독립적인 HTML 응용 프로그램입니다.
2015년에 Jupyter notebook은 IPython 코드 기반의 The Big Split™ 의 일부분으로 시작되었습니다.
IPython 3는 *IPython notebook* 과 같은 사용자 언어에 독립적인 코드와 *IPython kernel for Python* 과 같은 특정 언어 기반의 코드의 기능을 가지고 출시되었습니다.
컴퓨터에는 많은 언어가 사용되기 때문에, Jupyter 프로젝트는 사용자 언어에 독립적인 **Jupyter notebook** 을 이 저장소와 개인의 독립적인 저장소에 있는 특정 언어 중심의 커널의 도움으로 지속적으로 개발할 것입니다.
[[The Big Split™ announcement](https://blog.jupyter.org/2015/04/15/the-big-split/)]
[[Jupyter Ascending blog post](http://blog.jupyter.org/2015/08/12/first-release-of-jupyter/)]

## 설치
설치법 문서는 다음 주소에서 찾을 수 있습니다.
You can find the installation documentation for the
[Jupyter platform, on ReadTheDocs](https://jupyter.readthedocs.io/en/latest/install.html).
조금 더 심화된 Jupyter notebook의 사용은 다음 주소에서 볼 수 있습니다.
[here](https://jupyter-notebook.readthedocs.io/en/latest/).

설치를 위해서는 
[pip installed](https://pip.readthedocs.io/en/stable/installing/) 가 있는지 확인한 후 다음을 실행해주세요:

    $ pip install notebook

## 활용 - Jupyter notebook 실행하기

### 로컬에서 실행할 때

이와 같이 실행하세요:

    $ jupyter notebook

## 개발 설치

[`CONTRIBUTING.rst`](CONTRIBUTING.rst) 을 통해 설치법을 확인하세요.

## 기여하기

이 프로젝트에 기여를 하고 싶다면, [`CONTRIBUTING.rst`](CONTRIBUTING.rst) 을 참고해주세요.

## 자료
- [Project Jupyter website](https://jupyter.org)
- [Online Demo at try.jupyter.org](https://try.jupyter.org)
- [Documentation for Jupyter notebook](https://jupyter-notebook.readthedocs.io/en/latest/) [[PDF](https://media.readthedocs.org/pdf/jupyter-notebook/latest/jupyter-notebook.pdf)]
- [Documentation for Project Jupyter](https://jupyter.readthedocs.io/en/latest/index.html) [[PDF](https://media.readthedocs.org/pdf/jupyter/latest/jupyter.pdf)]
- [Issues](https://github.com/jupyter/notebook/issues)
- [Technical support - Jupyter Google Group](https://groups.google.com/forum/#!forum/jupyter)
# Notebook 실행하기

 ## 첫 걸음 
  1. 다음 명령어를 통해 Notebook 서버를 시작하세요 :

	$ jupyter notebook

  2. 브라우저에 Notebook이 실행된 것을 확인할 수 있습니다.


# Notebook 서버 시작하기
 
  Notebook을 컴퓨터에 설치하였으면 Notebook 서버를 시작할 수 있습니다. 다음 명령어를 이용하여 Notebook서버를 시작할 수 있습니다.

	$ jupyter notebook

 이 명령어를 실행하면, 터미널에 웹 응용프로그램의 주소와 서버에 대한 정보가 출력됩니다.

	$ jupyter notebook
	$ [I 08:58:24.417 NotebookApp] Serving notebooks from local directory: /Users/catherline
	$ [I 08:58:24.417 NotebookApp] 0 active kernels
	$ [I 08:58:24.417 NotebookApp] The Jupyter Notebook is running at: http://localhost:8888/ 
	$ [I 08:58:24.417 NotebookApp] Use Control-C to stop this server and shut down all kernels

 기본 브라우저를 통해 이 주소가 열립니다.

 Notebook이 브라우저에 열리면, Notebook의 목록을 보여주는 Notebook Dashboard를 볼 수 있습니다. 대체로 가장 상위의 디렉토리를 열어줄 것입니다. 

 **Notebook Dashboard**

![Notebook Dashboard example](resources/dashboard.GIF "Notebook Dashboard")

# Notebook 서버의 명령어 소개

 ## 커스텀 IP 나 포트를 이용하여 시작하려면 어떻게 해야할까?
  
 기본값으로, Notebook 서버는 포트 8888로 시작됩니다. 만약 포트8888이 사용할 수 없다면, Notebook 서버는 다른 가능한 포트를 찾습니다. 또한 임의로 포트를 설정해주는 것도 가능합니다. 예를 들어 포트 9999로 실행하면 :
 
	$ jupyter notebook --port 9999


 ## 브라우저를 열지않고 Notebook를 열기

 브라우저를 열지 않고 Notebook 서버를 시작하려면 :

	$ jupyter notebook --no-browser

 
 ## Notebook 서버 옵션 도움말 보기

 Notebook 서버는 --help 옵션을 통해 도움말 메시지를 제공합니다 :

	$ jupyter notebook --help


 

# UI 기능

버그 리포트나 Jupyter Mailing list에 메일을 보내려고 할 때, 개발자나 사용자들이 버그를 진단하거나 해결해줄 경우 다른 UI를 사용하면 시간이 단축된다.
이번 장에서는 Notebook과 Notebook의 다른 모드의 UI 요소를 알려줄 것이다.

 ## Notebook Dashboard

jupyter notebook 명령어를 실행하면 가장 먼저  Notebook Dashboard가 나타난다.

![Notebook Dashboard example](resources/dashboard.GIF "Notebook Dashboard")



 ## Notebook 편집기

편집을 위해 Notebook을 선택했다면, Notebook은 Notebook편집기를 열어준다.

![Notebook Editor example](resources/Notebook_Editor.GIF "Notebook Editor")
 ## Notebook 의 사용자 도움 인터페이스

만약 Notebook 편집기의 특정 요소를 더 배우고 싶다면, 도움 메뉴 - 사용자 인터페이스 를 선택함으로써 사용사 인터페이스 도움말을 볼 수 있습니다.

 ## 편집 모드와 Notebook편집기

셀이 편집모드에 있다면, 셀 모드 지시자는 셀의 상태를 반영합니다. 이 상태는 오른쪽 위의 작은 연필모양으로 선택가능합니다. 셀이 명령 모드에 있다면, 그 위치에 아이콘이 없습니다.

![Edit Mode example](resources/edit_mode.GIF "Edit Mode")

 ## 파일 편집기

이제 Notebook Dashboard 안의 Notebook 파일이 아닌 표시된 파일을 선택하여 열어야한다고 한다면, 파일은 파일 편집기로 열립니다.

![File Editor example](resources/file_editor.GIF "File Editor")





# Jupyter Notebook 설치하기

## 필요한 것 : Python

Jupyter Notebook 을 설치하기 위해선 Jupyter가 많은 프로그래밍 언어들로 동작되기 때문에, Python이 필요합니다. (Python 3.3이상, Python 2.7)

Python과 Jupyter를 설치할 때 Anaconda를 이용하는 것을 추천합니다. 밑에서 이를 이용하여 설치할 것입니다.

## Anaconda 와 conda 를 이용하여 Jupyter 설치하기

새로운 이용자들은 Anaconda를 설치하는 것을 강력하게 추천합니다. Anaconda는 Python과 Jupyter를 쉽게 설치하게 해주고, 과학적인 계산과 데이터를 위한 자주 사용되는 패키지들의 설치에도 유용합니다.

설치 순서 : 
   
   1. Anaconda를 다운받으세요. Anaconda의 가장 최신의 Python 3버전을 다운 받는 것을 추천합니다.
   2. 다운 받은 Anaconda 의 다운로드 페이지에 있는 설명을 읽고 설치해주세요.
   3. 축하합니다. Jupyter Notebook 을 설치하셨습니다. Jupyter Notebook을 실행하려면 :

	$ jupyter notebook

## 숙련된 Python 이용자 : pip을 통해 설치하기
 
 Python 이용자라면, Anaconda 대신에 Python의 패키지 매니저 pip을 이용하여 설치하세요.

 첫째로, 가장 최신의 pip인지를 확인하세요; 구 버전은 독립성에 문제가 있을 수 있습니다.
   
    $ pip install --upgrade pip

 이제 다음을 이용하여 Jupyter Notebook를 설치하세요 :

    $ pip install jupyter

 (축하합니다. Jupyter Notebook를 설치하셨습니다.)
# Jupyter Notebook

[![Google Group](https://img.shields.io/badge/-Google%20Group-lightgrey.svg)](https://groups.google.com/forum/#!forum/jupyter)
[![Build Status](https://travis-ci.org/jupyter/notebook.svg?branch=master)](https://travis-ci.org/jupyter/notebook)
[![Documentation Status](https://readthedocs.org/projects/jupyter-notebook/badge/?version=latest)](https://jupyter-notebook.readthedocs.io/en/latest/?badge=latest)
                


Jupyter नोटबुक इंटरैक्टिव के लिए एक वेब-आधारित नोटबुक वातावरण है
कंप्यूटिंग।

![Jupyter notebook example](resources/running_code_med.png "Jupyter notebook example")

### नोटिस
कृपया ध्यान दें कि इस भंडार का रखरखाव वर्तमान में जुपिटर समुदाय के एक कंकाल के दल द्वारा किया जाता है। हम उपयोगकर्ताओं को जुपिटरलैब में संक्रमण के लिए प्रोत्साहित करते हैं, जहां अधिक तत्काल समर्थन हो सकता है। हमारा दृष्टिकोण आगे बढ़ेगा:

1. जुपिटर नोटबुक की सुरक्षा बनाए रखने के लिए। इसका मतलब है कि सुरक्षा से संबंधित मुद्दे और पुल अनुरोध हमारी सर्वोच्च प्राथमिकता है।
2. JupyterLab को संबोधित करने के लिए [समता मुद्दों की सुविधा](https://github.com/jupyterlab/jupyterlab/issues?q=is%3Aopen+is%3Aissue+label%3A%22tag%3AFeature+Parity%22)| इस प्रयास के हिस्से के रूप में, हम एक बेहतर [नोटबुक-ओनली एक्सपीरियंस](https://github.com/jupyterlab/jupyterlab/issues/8450)JupyterLab में उन उपयोगकर्ताओं के लिए जो क्लासिक Jupyter नोटबुक के UI को पसंद करते हैं।
3. समुदाय के सदस्यों की कड़ी मेहनत के प्रति उत्तरदायी होना जिन्होंने पुल अनुरोधों को खोला है। हम इन पीआर को ट्राई कर रहे हैं। हम इस समय नई सुविधाओं का समर्थन या रखरखाव नहीं कर सकते हैं, लेकिन हम सुरक्षा और अन्य स्थिरता सुधारों का स्वागत करते हैं।

यदि आपके पास एक नई सुविधा के साथ एक खुला पुल अनुरोध है या यदि आप एक खोलने की योजना बना रहे हैं, तो कृपया इसे [नोटबुक एक्सटेंशन](https://jupyter-notebook.readthedocs.io/en/stable/extending/) के रूप में शिपिंग करने पर विचार करें। बजाय।

##### `नोटबुक` में योगदान करने के लिए विकल्प
इसके अतिरिक्त, कृपया विचार करें कि क्या आपका योगदान Jupyter फ्रंट-एंड के लिए अंतर्निहित सर्वर के लिए उपयुक्त होगा, [jupyter server](https://github.com/jupyter/jupyter_server) या में [JupyterLab फ़्रंट एंड](https://github.com/jupyterlab/jupyterlab).

### जुपिटर नोटबुक, आइपीथॉन नोटबुक की भाषा-अज्ञेय विकास
Jupyter नोटबुक एक भाषा-अज्ञेय HTML नोटबुक अनुप्रयोग है
प्रोजेक्ट जुपिटर। 2015 में, जुपिटर नोटबुक के एक भाग के रूप में जारी किया गया था
IPython कोडबेस का बिग स्प्लिट ™। IPython 3 अंतिम प्रमुख अखंड था
दोनों भाषा-अज्ञेयवादी कोड, जैसे *IPython नोटबुक*,
और भाषा विशिष्ट कोड, जैसे कि *अजगर के लिए आईपीथॉन कर्नेल*। जैसा
कई भाषाओं में कंप्यूटिंग स्पैन, प्रोजेक्ट जुपिटर विकसित करना जारी रखेगा
भाषा-अज्ञेय **जुपिटर नोटबुक** इस रेपो में और की मदद से
समुदाय भाषा विशिष्ट गुठली विकसित करते हैं जो अपने आप में पाए जाते हैं
असतत रेपो।
[[Big Split™ घोषणा](https://blog.jupyter.org/the-big-split-9d7b88a031a7)]
[[Jupyter आरोही ब्लॉग पोस्ट](https://blog.jupyter.org/jupyter-ascending-1bf5b362d97e)]

## स्थापना
आप के लिए स्थापना प्रलेखन पा सकते हैं
[बृहस्पति मंच, ReadTheDocs पर](https://jupyter.readthedocs.io/en/latest/install.html).
जुपिटर नोटबुक के उन्नत उपयोग के लिए दस्तावेज पाया जा सकता है
[यहाँ](https://jupyter-notebook.readthedocs.io/en/latest/).

स्थानीय स्थापना के लिए, सुनिश्चित करें कि आपके पास है
[pip स्थापित](https://pip.readthedocs.io/en/stable/installing/) और भाग खड़ा हुआ:

    $ pip install notebook

## उपयोग - जुपिटर नोटबुक चल रहा है

### स्थानीय स्थापना में चल रहा है

इसके साथ लॉन्च करें:

    $ jupyter notebook

### एक दूरस्थ स्थापना में चल रहा है

आपको बृहस्पति नोटबुक को दूरस्थ रूप से शुरू करने से पहले कुछ कॉन्फ़िगरेशन की आवश्यकता है। देखें [नोटबुक सर्वर चला रहा है](https://jupyter-notebook.readthedocs.io/en/stable/public_server.html).

## विकास स्थापना

स्थानीय विकास की स्थापना कैसे करें, इसके लिए [`CONTRIBUTING.rst`](CONTRIBUTING.rst) देखें।

## योगदान

यदि आप इस परियोजना में योगदान देने में रुचि रखते हैं, तो [`CONTRIBUTING.rst`](CONTRIBUTING.rst) देखें।

## साधन
- [Project Jupyter website](https://jupyter.org)
- [Online Demo at jupyter.org/try](https://jupyter.org/try)
- [Documentation for Jupyter notebook](https://jupyter-notebook.readthedocs.io/en/latest/) [[PDF](https://media.readthedocs.org/pdf/jupyter-notebook/latest/jupyter-notebook.pdf)]
- [Korean Version of Installation](https://github.com/ChungJooHo/Jupyter_Kor_doc/)
- [Documentation for Project Jupyter](https://jupyter.readthedocs.io/en/latest/index.html) [[PDF](https://media.readthedocs.org/pdf/jupyter/latest/jupyter.pdf)]
- [Issues](https://github.com/jupyter/notebook/issues)
- [Technical support - Jupyter Google Group](https://groups.google.com/forum/#!forum/jupyter) # Jupyter Notebook

[![Google Group](https://img.shields.io/badge/-Google%20Group-lightgrey.svg)](https://groups.google.com/forum/#!forum/jupyter)
[![Build Status](https://travis-ci.org/jupyter/notebook.svg?branch=master)](https://travis-ci.org/jupyter/notebook)
[![Documentation Status](https://readthedocs.org/projects/jupyter-notebook/badge/?version=latest)](https://jupyter-notebook.readthedocs.io/en/latest/?badge=latest)
                


Jupyter Notebook是用于交互的基于Web的笔记本环境
计算。

![Jupyter notebook example](resources/running_code_med.png "Jupyter notebook example")

### 注意
请注意，这家商店目前由木星社区的骨干团队维护。我们鼓励用户过渡到 JupyterLab，那里可能会立即提供更多支持。我们的方法将向前发展：

1.维护Jupiter笔记本电脑的安全性。这意味着与安全相关的问题和请求是我们的首要任务。
2.解决JupyterLab [促进平等问题](https://github.com/jupyterlab/jupyterlab/issues?q=is%3Aopen+is%3Aissue+label%3A%22tag%3AFeature+Parity%22)|作为这项工作的一部分，我们有更好的[仅限笔记本电脑的体验](https://github.com/jupyterlab/jupyterlab/issues/8450)在JupyterLab中，适合喜欢经典Jupyter笔记本UI的用户。
3.负责提出请求请求的社区成员的辛勤工作。我们正在尝试这些PR。我们目前无法支持或维护新设施，但是我们欢迎安全性和其他稳定性方面的改进。

如果您有一个具有新功能的打开请求请求，或者您打算打开一个请求，请将该请求命名为[notebook extension](https://jupyter-notebook.readthedocs.io/en/stable/extending/) 考虑运送为。代替。

##### 选择贡献“笔记本”
此外，请考虑您的贡献是否适合Jupyter前端的基础服务器， [jupyter server](https://github.com/jupyter/jupyter_server) 或在 [JupyterLab 前端](https://github.com/jupyterlab/jupyterlab).

### Jupyter笔记本，与IPython笔记本无关的语言开发
Jupyter Notebook是与语言无关的HTML Notebook应用程序
木星计划。 2015年，木星作为笔记本的一部分发布
IPython代码库的Big Split™。 IPython 3是最后一个主要的整体
两种与语言无关的代码，例如 *IPython notebook*，
以及特定语言的代码，例如 *用于Python的IPython内核* 。如
通过多种语言计算SPAN，Jupyter项目将继续发展
与语言无关 **Jupyter Notebook** 在此仓库中更多帮助下
社区开发自己发现的特定于语言的内核
离散回购。
[[Big Split™ 宣言](https://blog.jupyter.org/the-big-split-9d7b88a031a7)]
[[Jupyter 升序博客文章](https://blog.jupyter.org/jupyter-ascending-1bf5b362d97e)]

## 成立
您可以找到以下安装文件
[Jupiter论坛，在ReadTheDocs上](https://jupyter.readthedocs.io/en/latest/install.html).
可以找到有关Jupiter笔记本的高级使用的文档
[这里](https://jupyter-notebook.readthedocs.io/en/latest/).

对于本地安装，请确保您已经
[pip 成立时间](https://pip.readthedocs.io/en/stable/installing/) 并运行：

    $ pip install notebook

## 用法-运行木星笔记本

### 在本地安装中运行

与启动

    $ jupyter笔记本

### 在远程安装中运行

在远程启动Jupiter笔记本电脑之前，需要进行一些配置。请参阅 [运行笔记本服务器](https://jupyter-notebook.readthedocs.io/en/stable/public_server.html).

## 开发设置

有关如何建立本地发展 [`CONTRIBUTING.rst`](CONTRIBUTING.rst) 看到。

## 贡献

如果您有兴趣为这个项目做贡献，请参阅 [`CONTRIBUTING.rst`](CONTRIBUTING.rst).

## 资源
- [Project Jupyter website](https://jupyter.org)
- [Online Demo at jupyter.org/try](https://jupyter.org/try)
- [Documentation for Jupyter notebook](https://jupyter-notebook.readthedocs.io/en/latest/) [[PDF](https://media.readthedocs.org/pdf/jupyter-notebook/latest/jupyter-notebook.pdf)]
- [Korean Version of Installation](https://github.com/ChungJooHo/Jupyter_Kor_doc/)
- [Documentation for Project Jupyter](https://jupyter.readthedocs.io/en/latest/index.html) [[PDF](https://media.readthedocs.org/pdf/jupyter/latest/jupyter.pdf)]
- [Issues](https://github.com/jupyter/notebook/issues)
- [Technical support - Jupyter Google Group](https://groups.google.com/forum/#!forum/jupyter) Contributing to the Jupyter Notebook
====================================

If you're reading this section, you're probably interested in contributing to
Jupyter.  Welcome and thanks for your interest in contributing!

Please take a look at the Contributor documentation, familiarize yourself with
using the Jupyter Notebook, and introduce yourself on the mailing list and
share what area of the project you are interested in working on.

General Guidelines
------------------

For general documentation about contributing to Jupyter projects, see the
`Project Jupyter Contributor Documentation`__.

__ https://jupyter.readthedocs.io/en/latest/contributing/content-contributor.html


Setting Up a Development Environment
------------------------------------

Installing Node.js and npm
^^^^^^^^^^^^^^^^^^^^^^^^^^

Building the Notebook from its GitHub source code requires some tools to
create and minify JavaScript components and the CSS,
specifically Node.js and Node's package manager, ``npm``.
It should be node version ≥ 6.0.

If you use ``conda``, you can get them with::

    conda install -c conda-forge nodejs

If you use `Homebrew <https://brew.sh/>`_ on Mac OS X::

    brew install node

Installation on Linux may vary, but be aware that the `nodejs` or `npm` packages
included in the system package repository may be too old to work properly.

You can also use the installer from the `Node.js website <https://nodejs.org>`_.


Installing the Jupyter Notebook
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have installed the dependencies mentioned above, use the following
steps::

    pip install --upgrade setuptools pip
    git clone https://github.com/jupyter/notebook
    cd notebook
    pip install -e .

If you are using a system-wide Python installation and you only want to install the notebook for you,
you can add ``--user`` to the install commands.

Once you have done this, you can launch the master branch of Jupyter notebook
from any directory in your system with::

    jupyter notebook

Verification
^^^^^^^^^^^^

While running the notebook, select one of your notebook files (the file will have the extension ``.ipynb``).
In the top tab you will click on "Help" and then click on "About". In the pop window you will see information about the version of Jupyter that you are running. You will see "The version of the notebook server is:".
If you are working in development mode, you will see that your version of Jupyter notebook will include the word "dev". If it does not include the word "dev", you are currently not working in development mode and should follow the steps below to uninstall and reinstall Jupyter.

Troubleshooting the Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you do not see that your Jupyter Notebook is not running on dev mode, it's possible that you are
running other instances of Jupyter Notebook. You can try the following steps:

1. Uninstall all instances of the notebook package. These include any installations you made using
   pip or conda.
2. Run ``python3 -m pip install -e .`` in the notebook repository to install the notebook from there.
3. Run ``npm run build`` to make sure the Javascript and CSS are updated and compiled.
4. Launch with ``python3 -m notebook --port 8989``, and check that the browser is pointing to ``localhost:8989``
   (rather than the default 8888). You don't necessarily have to launch with port 8989, as long as you use
   a port that is neither the default nor in use, then it should be fine.
5. Verify the installation with the steps in the previous section.


Rebuilding JavaScript and CSS
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

There is a build step for the JavaScript and CSS in the notebook.
To make sure that you are working with up-to-date code, you will need to run
this command whenever there are changes to JavaScript or LESS sources::

    npm run build

**IMPORTANT:** Don't forget to run ``npm run build`` after switching branches.
When switching between branches of different versions (e.g. ``4.x`` and
``master``), run ``pip install -e .``. If you have tried the above and still
find that the notebook is not reflecting the current source code, try cleaning
the repo with ``git clean -xfd`` and reinstalling with ``pip install -e .``.

Development Tip
"""""""""""""""

When doing development, you can use this command to automatically rebuild
JavaScript and LESS sources as they are modified::

    npm run build:watch

Git Hooks
"""""""""

If you want to automatically update dependencies and recompile JavaScript and
CSS after checking out a new commit, you can install post-checkout and
post-merge hooks which will do it for you::

    git-hooks/install-hooks.sh

See ``git-hooks/README.md`` for more details.


Running Tests
-------------

Python Tests
^^^^^^^^^^^^

Install dependencies::

    pip install -e '.[test]'

To run the Python tests, use::

    pytest

If you want coverage statistics as well, you can run::

    py.test --cov notebook -v --pyargs notebook

JavaScript Tests
^^^^^^^^^^^^^^^^

To run the JavaScript tests, you will need to have PhantomJS and CasperJS
installed::

    npm install -g casperjs phantomjs-prebuilt

Then, to run the JavaScript tests::

    python -m notebook.jstest [group]

where ``[group]`` is an optional argument that is a path relative to
``notebook/tests/``.
For example, to run all tests in ``notebook/tests/notebook``::

    python -m notebook.jstest notebook

or to run just ``notebook/tests/notebook/deletecell.js``::

    python -m notebook.jstest notebook/deletecell.js


Building the Documentation
--------------------------

To build the documentation you'll need `Sphinx <http://www.sphinx-doc.org/>`_,
`pandoc <http://pandoc.org/>`_ and a few other packages.

To install (and activate) a conda environment named ``notebook_docs``
containing all the necessary packages (except pandoc), use::

    conda create -n notebook_docs pip
    conda activate notebook_docs  # Linux and OS X
    activate notebook_docs        # Windows
    pip install .[docs]

If you want to install the necessary packages with ``pip``, use the following instead::

    pip install .[docs]

Once you have installed the required packages, you can build the docs with::

    cd docs
    make html

After that, the generated HTML files will be available at
``build/html/index.html``. You may view the docs in your browser.

You can automatically check if all hyperlinks are still valid::

    make linkcheck

Windows users can find ``make.bat`` in the ``docs`` folder.

You should also have a look at the `Project Jupyter Documentation Guide`__.

__ https://jupyter.readthedocs.io/en/latest/contributing/docs-contributions/index.html
.. highlight:: sh

.. include:: ../../CONTRIBUTING.rst

.. _server_security:

Security in the Jupyter notebook server
=======================================

Since access to the Jupyter notebook server means access to running arbitrary code,
it is important to restrict access to the notebook server.
For this reason, notebook 4.3 introduces token-based authentication that is **on by default**.

.. note::

    If you enable a password for your notebook server,
    token authentication is not enabled by default,
    and the behavior of the notebook server is unchanged from versions earlier than 4.3.

When token authentication is enabled, the notebook uses a token to authenticate requests.
This token can be provided to login to the notebook server in three ways:

- in the ``Authorization`` header, e.g.::

    Authorization: token abcdef...

- In a URL parameter, e.g.::

    https://my-notebook/tree/?token=abcdef...

- In the password field of the login form that will be shown to you if you are not logged in.

When you start a notebook server with token authentication enabled (default),
a token is generated to use for authentication.
This token is logged to the terminal, so that you can copy/paste the URL into your browser::

    [I 11:59:16.597 NotebookApp] The Jupyter Notebook is running at:
    http://localhost:8888/?token=c8de56fa4deed24899803e93c227592aef6538f93025fe01


If the notebook server is going to open your browser automatically
(the default, unless ``--no-browser`` has been passed),
an *additional* token is generated for launching the browser.
This additional token can be used only once,
and is used to set a cookie for your browser once it connects.
After your browser has made its first request with this one-time-token,
the token is discarded and a cookie is set in your browser.

At any later time, you can see the tokens and URLs for all of your running servers with :command:`jupyter notebook list`::

    $ jupyter notebook list
    Currently running servers:
    http://localhost:8888/?token=abc... :: /home/you/notebooks
    https://0.0.0.0:9999/?token=123... :: /tmp/public
    http://localhost:8889/ :: /tmp/has-password

For servers with token-authentication enabled, the URL in the above listing will include the token,
so you can copy and paste that URL into your browser to login.
If a server has no token (e.g. it has a password or has authentication disabled),
the URL will not include the token argument.
Once you have visited this URL,
a cookie will be set in your browser and you won't need to use the token again,
unless you switch browsers, clear your cookies, or start a notebook server on a new port.

Alternatives to token authentication
------------------------------------

If a generated token doesn't work well for you,
you can set a password for your notebook.
:command:`jupyter notebook password` will prompt you for a password,
and store the hashed password in your :file:`jupyter_notebook_config.json`.

.. versionadded:: 5.0

    :command:`jupyter notebook password` command is added.


It is possible to disable authentication altogether by setting the token and password to empty strings,
but this is **NOT RECOMMENDED**, unless authentication or access restrictions are handled at a different layer in your web application:

.. sourcecode:: python

    c.NotebookApp.token = ''
    c.NotebookApp.password = ''


.. _notebook_security:

Security in notebook documents
==============================

As Jupyter notebooks become more popular for sharing and collaboration,
the potential for malicious people to attempt to exploit the notebook
for their nefarious purposes increases. IPython 2.0 introduced a
security model to prevent execution of untrusted code without explicit
user input.

The problem
-----------

The whole point of Jupyter is arbitrary code execution. We have no
desire to limit what can be done with a notebook, which would negatively
impact its utility.

Unlike other programs, a Jupyter notebook document includes output.
Unlike other documents, that output exists in a context that can execute
code (via Javascript).

The security problem we need to solve is that no code should execute
just because a user has **opened** a notebook that **they did not
write**. Like any other program, once a user decides to execute code in
a notebook, it is considered trusted, and should be allowed to do
anything.

Our security model
------------------

-  Untrusted HTML is always sanitized
-  Untrusted Javascript is never executed
-  HTML and Javascript in Markdown cells are never trusted
-  **Outputs** generated by the user are trusted
-  Any other HTML or Javascript (in Markdown cells, output generated by
   others) is never trusted
-  The central question of trust is "Did the current user do this?"

The details of trust
--------------------

When a notebook is executed and saved, a signature is computed from a
digest of the notebook's contents plus a secret key. This is stored in a
database, writable only by the current user. By default, this is located at::

    ~/.local/share/jupyter/nbsignatures.db  # Linux
    ~/Library/Jupyter/nbsignatures.db       # OS X
    %APPDATA%/jupyter/nbsignatures.db       # Windows

Each signature represents a series of outputs which were produced by code the
current user executed, and are therefore trusted.

When you open a notebook, the server computes its signature, and checks if it's
in the database. If a match is found, HTML and Javascript
output in the notebook will be trusted at load, otherwise it will be
untrusted.

Any output generated during an interactive session is trusted.

Updating trust
**************

A notebook's trust is updated when the notebook is saved. If there are
any untrusted outputs still in the notebook, the notebook will not be
trusted, and no signature will be stored. If all untrusted outputs have
been removed (either via ``Clear Output`` or re-execution), then the
notebook will become trusted.

While trust is updated per output, this is only for the duration of a
single session. A newly loaded notebook file is either trusted or not in its
entirety.

Explicit trust
**************

Sometimes re-executing a notebook to generate trusted output is not an
option, either because dependencies are unavailable, or it would take a
long time. Users can explicitly trust a notebook in two ways:

-  At the command-line, with::

    jupyter trust /path/to/notebook.ipynb

-  After loading the untrusted notebook, with ``File / Trust Notebook``

These two methods simply load the notebook, compute a new signature, and add
that signature to the user's database.

Reporting security issues
-------------------------

If you find a security vulnerability in Jupyter, either a failure of the
code to properly implement the model described here, or a failure of the
model itself, please report it to security@ipython.org.

If you prefer to encrypt your security reports,
you can use :download:`this PGP public key <ipython_security.asc>`.

Affected use cases
------------------

Some use cases that work in Jupyter 1.0 became less convenient in
2.0 as a result of the security changes. We do our best to minimize
these annoyances, but security is always at odds with convenience.

Javascript and CSS in Markdown cells
************************************

While never officially supported, it had become common practice to put
hidden Javascript or CSS styling in Markdown cells, so that they would
not be visible on the page. Since Markdown cells are now sanitized (by
`Google Caja <https://developers.google.com/caja>`__), all Javascript
(including click event handlers, etc.) and CSS will be stripped.

We plan to provide a mechanism for notebook themes, but in the meantime
styling the notebook can only be done via either ``custom.css`` or CSS
in HTML output. The latter only have an effect if the notebook is
trusted, because otherwise the output will be sanitized just like
Markdown.

Collaboration
*************

When collaborating on a notebook, people probably want to see the
outputs produced by their colleagues' most recent executions. Since each
collaborator's key will differ, this will result in each share starting
in an untrusted state. There are three basic approaches to this:

-  re-run notebooks when you get them (not always viable)
-  explicitly trust notebooks via ``jupyter trust`` or the notebook menu
   (annoying, but easy)
-  share a notebook signatures database, and use configuration dedicated to the
   collaboration while working on the project.

To share a signatures database among users, you can configure:

.. code-block:: python

    c.NotebookNotary.data_dir = "/path/to/signature_dir"

to specify a non-default path to the SQLite database (of notebook hashes,
essentially). We are aware that SQLite doesn't work well on NFS and we are
`working out better ways to do this <https://github.com/jupyter/notebook/issues/1782>`_.
.. _frontend_config:

Configuring the notebook frontend
=================================

.. note::

    The ability to configure the notebook frontend UI and preferences is
    still a work in progress.

This document is a rough explanation on how you can persist some configuration
options for the notebook JavaScript.

There is no exhaustive list of all the configuration options as most options
are passed down to other libraries, which means that non valid
configuration can be ignored without any error messages.


How front end configuration works
---------------------------------
The frontend configuration system works as follows:

  - get a handle of a configurable JavaScript object.
  - access its configuration attribute.
  - update its configuration attribute with a JSON patch.


Example - Changing the notebook's default indentation
-----------------------------------------------------
This example explains how to change the default setting ``indentUnit``
for CodeMirror Code Cells::

    var cell = Jupyter.notebook.get_selected_cell();
    var config = cell.config;
    var patch = {
          CodeCell:{
            cm_config:{indentUnit:2}
          }
        }
    config.update(patch)

You can enter the previous snippet in your browser's JavaScript console once.
Then reload the notebook page in your browser. Now, the preferred indent unit
should be equal to two spaces. The custom setting persists and you do not need
to reissue the patch on new notebooks.

``indentUnit``, used in this example, is one of the many `CodeMirror options
<https://codemirror.net/doc/manual.html#option_indentUnit>`_ which are available
for configuration.

You can similarly change the options of the file editor by entering the following
snippet in the browser's Javascript console once (from a file editing page).::

   var config = Jupyter.editor.config
   var patch = {
         Editor: {
           codemirror_options: {
             indentUnit: 2
           }
         }
       }
   config.update(patch)

Example - Restoring the notebook's default indentation
------------------------------------------------------
If you want to restore a notebook frontend preference to its default value,
you will enter a JSON patch with a ``null`` value for the preference setting.

For example, let's restore the indent setting ``indentUnit`` to its default of
four spaces. Enter the following code snippet in your JavaScript console::

    var cell = Jupyter.notebook.get_selected_cell();
    var config = cell.config;
    var patch = {
          CodeCell:{
            cm_config:{indentUnit: null} // only change here.
          }
        }
    config.update(patch)

Reload the notebook in your browser and the default indent should again be two
spaces.

Persisting configuration settings
---------------------------------
Under the hood, Jupyter will persist the preferred configuration settings in
``~/.jupyter/nbconfig/<section>.json``, with ``<section>``
taking various value depending on the page where the configuration is issued.
``<section>`` can take various values like ``notebook``, ``tree``, and
``editor``. A ``common`` section contains configuration settings shared by all
pages.
.. _htmlnotebook:

The Jupyter Notebook
====================

Introduction
------------

The notebook extends the console-based approach to interactive computing in
a qualitatively new direction, providing a web-based application suitable for
capturing the whole computation process: developing, documenting, and
executing code, as well as communicating the results.  The Jupyter notebook
combines two components:

**A web application**: a browser-based tool for interactive authoring of
documents which combine explanatory text, mathematics, computations and their
rich media output.

**Notebook documents**: a representation of all content visible in the web
application, including inputs and outputs of the computations, explanatory
text, mathematics, images, and rich media representations of objects.

.. seealso::

    See the :ref:`installation guide <jupyter:install>` on how to install the
    notebook and its dependencies.


Main features of the web application
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* In-browser editing for code, with automatic syntax highlighting,
  indentation, and tab completion/introspection.

* The ability to execute code from the browser, with the results of
  computations attached to the code which generated them.

* Displaying the result of computation using rich media representations, such
  as HTML, LaTeX, PNG, SVG, etc. For example, publication-quality figures
  rendered by the matplotlib_ library, can be included inline.

* In-browser editing for rich text using the Markdown_ markup language, which
  can provide commentary for the code, is not limited to plain text.

* The ability to easily include mathematical notation within markdown cells
  using LaTeX, and rendered natively by MathJax_.



.. _MathJax: https://www.mathjax.org/


Notebook documents
~~~~~~~~~~~~~~~~~~
Notebook documents contains the inputs and outputs of a interactive session as
well as additional text that accompanies the code but is not meant for
execution.  In this way, notebook files can serve as a complete computational
record of a session, interleaving executable code with explanatory text,
mathematics, and rich representations of resulting objects. These documents
are internally JSON_ files and are saved with the ``.ipynb`` extension. Since
JSON is a plain text format, they can be version-controlled and shared with
colleagues.

.. _JSON: https://en.wikipedia.org/wiki/JSON

Notebooks may be exported to a range of static formats, including HTML (for
example, for blog posts), reStructuredText, LaTeX, PDF, and slide shows, via
the nbconvert_ command.

Furthermore, any  ``.ipynb`` notebook document available from a public
URL can be shared via the Jupyter Notebook Viewer <nbviewer>.
This service loads the notebook document from the URL and renders it as a
static web page.  The results may thus be shared with a colleague, or as a
public blog post, without other users needing to install the Jupyter notebook
themselves.  In effect, nbviewer is simply nbconvert_ as
a web service, so you can do your own static conversions with nbconvert,
without relying on nbviewer.



.. seealso::

    :ref:`Details on the notebook JSON file format <nbformat:notebook_file_format>`


Notebooks and privacy
~~~~~~~~~~~~~~~~~~~~~

Because you use Jupyter in a web browser, some people are understandably
concerned about using it with sensitive data.
However, if you followed the standard
`install instructions <https://jupyter.readthedocs.io/en/latest/install.html>`_,
Jupyter is actually running on your own computer.
If the URL in the address bar starts with ``http://localhost:`` or
``http://127.0.0.1:``, it's your computer acting as the server.
Jupyter doesn't send your data anywhere else—and as it's open source,
other people can check that we're being honest about this.

You can also use Jupyter remotely:
your company or university might run the server for you, for instance.
If you want to work with sensitive data in those cases,
talk to your IT or data protection staff about it.

We aim to ensure that other pages in your browser or other users on the same
computer can't access your notebook server. See :ref:`server_security` for
more about this.


Starting the notebook server
----------------------------

You can start running a notebook server from the command line using the
following command::

    jupyter notebook

This will print some information about the notebook server in your console,
and open a web browser to the URL of the web application (by default,
``http://127.0.0.1:8888``).

The landing page of the Jupyter notebook web application, the **dashboard**,
shows the notebooks currently available in the notebook directory (by default,
the directory from which the notebook server was started).

You can create new notebooks from the dashboard with the ``New Notebook``
button, or open existing ones by clicking on their name.  You can also drag
and drop ``.ipynb`` notebooks and standard ``.py`` Python source code files
into the notebook list area.

When starting a notebook server from the command line, you can also open a
particular notebook directly, bypassing the dashboard, with ``jupyter notebook
my_notebook.ipynb``. The ``.ipynb`` extension is assumed if no extension is
given.

When you are inside an open notebook, the `File | Open...` menu option will
open the dashboard in a new browser tab, to allow you to open another notebook
from the notebook directory or to create a new notebook.


.. note::

   You can start more than one notebook server at the same time, if you want
   to work on notebooks in different directories.  By default the first
   notebook server starts on port 8888, and later notebook servers search for
   ports near that one.  You can also manually specify the port with the
   ``--port`` option.

Creating a new notebook document
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A new notebook may be created at any time, either from the dashboard, or using
the :menuselection:`File --> New` menu option from within an active notebook.
The new notebook is created within the same directory and will open in a new
browser tab. It will also be reflected as a new entry in the notebook list on
the dashboard.

.. image:: _static/images/new-notebook.gif


Opening notebooks
~~~~~~~~~~~~~~~~~
An open notebook has **exactly one** interactive session connected to a
kernel, which will execute code sent by the user
and communicate back results.  This kernel remains active if the web browser
window is closed, and reopening the same notebook from the dashboard will
reconnect the web application to the same kernel. In the dashboard, notebooks
with an active kernel have a ``Shutdown`` button next to them, whereas
notebooks without an active kernel have a ``Delete`` button in its place.

Other clients may connect to the same kernel.
When each kernel is started, the notebook server prints to the terminal a
message like this::

    [NotebookApp] Kernel started: 87f7d2c0-13e3-43df-8bb8-1bd37aaf3373

This long string is the kernel's ID which is sufficient for getting the
information necessary to connect to the kernel. If the notebook uses the IPython
kernel, you can also see this
connection data by running the ``%connect_info`` :ref:`magic
<magics_explained>`, which will print the same ID information along with other
details.

You can then, for example, manually start a Qt console connected to the *same*
kernel from the command line, by passing a portion of the ID::

    $ jupyter qtconsole --existing 87f7d2c0

Without an ID, ``--existing`` will  connect to the most recently
started kernel.

With the IPython kernel, you can also run the ``%qtconsole``
:ref:`magic <magics_explained>` in the notebook to open a Qt console connected
to the same kernel.

.. seealso::

    :ref:`ipythonzmq`

Notebook user interface
-----------------------

When you create a new notebook document, you will be presented with the
**notebook name**, a **menu bar**, a **toolbar** and an empty **code cell**.

.. image:: ./_static/images/blank-notebook-ui.png

**Notebook name**: The name displayed at the top of the page,
next to the Jupyter logo, reflects the name of the ``.ipynb`` file.
Clicking on the notebook name brings up a dialog which allows you to rename it.
Thus, renaming a notebook
from "Untitled0" to "My first notebook" in the browser, renames the
``Untitled0.ipynb`` file to ``My first notebook.ipynb``.

**Menu bar**: The menu bar presents different options that may be used to
manipulate the way the notebook functions.

**Toolbar**: The tool bar gives a quick way of performing the most-used
operations within the notebook, by clicking on an icon.

**Code cell**: the default type of cell; read on for an explanation of cells.


Structure of a notebook document
--------------------------------

The notebook consists of a sequence of cells.  A cell is a multiline text input
field, and its contents can be executed by using :kbd:`Shift-Enter`, or by
clicking either the "Play" button the toolbar, or :guilabel:`Cell`, :guilabel:`Run` in the menu bar.
The execution behavior of a cell is determined by the cell's type.  There are three
types of cells: **code cells**, **markdown cells**, and **raw cells**.  Every
cell starts off being a **code cell**, but its type can be changed by using a
drop-down on the toolbar (which will be "Code", initially), or via
:ref:`keyboard shortcuts <keyboard-shortcuts>`.

For more information on the different things you can do in a notebook,
see the `collection of examples
<https://nbviewer.jupyter.org/github/jupyter/notebook/tree/master/docs/source/examples/Notebook/>`_.

Code cells
~~~~~~~~~~
A *code cell* allows you to edit and write new code, with full syntax
highlighting and tab completion. The programming language you use depends
on the *kernel*, and the default kernel (IPython) runs Python code.

When a code cell is executed, code that it contains is sent to the kernel
associated with the notebook.  The results that are returned from this
computation  are then displayed in the notebook as the cell's *output*. The
output is not limited to text, with many other possible forms of output are
also possible, including ``matplotlib`` figures and HTML tables (as used, for
example, in the ``pandas`` data analysis package). This is known as IPython's
*rich display* capability.

.. seealso::

   `Rich Output`_  example notebook

Markdown cells
~~~~~~~~~~~~~~
You can document the computational process in a literate way, alternating
descriptive text with code, using *rich text*. In IPython this is accomplished
by marking up text with the Markdown language. The corresponding cells are
called *Markdown cells*. The Markdown language provides a simple way to
perform this text markup, that is, to specify which parts of the text should
be emphasized (italics), bold, form lists, etc.

If you want to provide structure for your document, you can use markdown
headings. Markdown headings consist of 1 to 6 hash # signs ``#`` followed by a
space and the title of your section. The markdown heading will be converted
to a clickable link for a section of the notebook. It is also used as a hint
when exporting to other document formats, like PDF.

When a Markdown cell is executed, the Markdown code is converted into
the corresponding formatted rich text. Markdown allows arbitrary HTML code for
formatting.

Within Markdown cells, you can also include *mathematics* in a straightforward
way, using standard LaTeX notation: ``$...$`` for inline mathematics and
``$$...$$`` for displayed mathematics. When the Markdown cell is executed,
the LaTeX portions are automatically rendered in the HTML output as equations
with high quality typography. This is made possible by MathJax_, which
supports a `large subset <https://docs.mathjax.org/en/latest/input/tex/index.html>`__ of LaTeX functionality

Standard mathematics environments defined by LaTeX and AMS-LaTeX (the
``amsmath`` package) also work, such as
``\begin{equation}...\end{equation}``, and ``\begin{align}...\end{align}``.
New LaTeX macros may be defined using standard methods,
such as ``\newcommand``, by placing them anywhere *between math delimiters* in
a Markdown cell. These definitions are then available throughout the rest of
the IPython session.

.. seealso::

    `Working with Markdown Cells`_ example notebook

Raw cells
~~~~~~~~~

*Raw* cells provide a place in which you can write *output* directly.
Raw cells are not evaluated by the notebook.
When passed through nbconvert_, raw cells arrive in the
destination format unmodified. For example, you can type full LaTeX
into a raw cell, which will only be rendered by LaTeX after conversion by
nbconvert.


Basic workflow
--------------

The normal workflow in a notebook is, then, quite similar to a standard
IPython session, with the difference that you can edit cells in-place multiple
times until you obtain the desired results, rather than having to
rerun separate scripts with the ``%run`` magic command.


Typically, you will work on a computational problem in pieces, organizing
related ideas into cells and moving forward once previous parts work
correctly. This is much more convenient for interactive exploration than
breaking up a computation into scripts that must be executed together, as was
previously necessary, especially if parts of them take a long time to run.

To interrupt a calculation which is taking too long, use the :guilabel:`Kernel`,
:guilabel:`Interrupt` menu option, or the :kbd:`i,i` keyboard shortcut.
Similarly, to restart the whole computational process,
use the :guilabel:`Kernel`, :guilabel:`Restart` menu option or :kbd:`0,0`
shortcut.

A notebook may be downloaded as a ``.ipynb`` file or converted to a number of
other formats using the menu option :guilabel:`File`, :guilabel:`Download as`.

.. seealso::

    `Running Code in the Jupyter Notebook`_ example notebook

    `Notebook Basics`_ example notebook

.. _keyboard-shortcuts:

Keyboard shortcuts
~~~~~~~~~~~~~~~~~~
All actions in the notebook can be performed with the mouse, but keyboard
shortcuts are also available for the most common ones. The essential shortcuts
to remember are the following:

* :kbd:`Shift-Enter`:  run cell
    Execute the current cell, show any output, and jump to the next cell below.
    If :kbd:`Shift-Enter` is invoked on the last cell, it makes a new cell below.
    This is equivalent to clicking the :guilabel:`Cell`, :guilabel:`Run` menu
    item, or the Play button in the toolbar.

* :kbd:`Esc`: Command mode
    In command mode, you can navigate around the notebook using keyboard shortcuts.

* :kbd:`Enter`: Edit mode
    In edit mode, you can edit text in cells.

For the full list of available shortcuts, click :guilabel:`Help`,
:guilabel:`Keyboard Shortcuts` in the notebook menus.

Plotting
--------
One major feature of the Jupyter notebook is the ability to display plots that
are the output of running code cells. The IPython kernel is designed to work
seamlessly with the matplotlib_ plotting library to provide this functionality.
Specific plotting library integration is a feature of the kernel.

Installing kernels
------------------

For information on how to install a Python kernel, refer to the
`IPython install page <https://ipython.org/install.html>`__.

The Jupyter wiki has a long list of `Kernels for other languages
<https://github.com/jupyter/jupyter/wiki/Jupyter-kernels>`_.
They usually come with instructions on how to make the kernel available
in the notebook.


.. _signing_notebooks:

Trusting Notebooks
------------------

To prevent untrusted code from executing on users' behalf when notebooks open,
we store a signature of each trusted notebook.
The notebook server verifies this signature when a notebook is opened.
If no matching signature is found,
Javascript and HTML output will not be displayed
until they are regenerated by re-executing the cells.

Any notebook that you have fully executed yourself will be
considered trusted, and its HTML and Javascript output will be displayed on
load.

If you need to see HTML or Javascript output without re-executing,
and you are sure the notebook is not malicious, you can tell Jupyter to trust it
at the command-line with::

    $ jupyter trust mynotebook.ipynb

See :ref:`notebook_security` for more details about the trust mechanism.

Browser Compatibility
---------------------

The Jupyter Notebook aims to support the latest versions of these browsers:

* Chrome
* Safari
* Firefox

Up to date versions of Opera and Edge may also work, but if they don't, please
use one of the supported browsers.

Using Safari with HTTPS and an untrusted certificate is known to not work
(websockets will fail).

.. include:: links.txt
What to do when things go wrong
===============================

First, have a look at the common problems listed below. If you can figure it out
from these notes, it will be quicker than asking for help.

Check that you have the latest version of any packages that look relevant.
Unfortunately it's not always easy to figure out what packages are relevant,
but if there was a bug that's already been fixed,
it's easy to upgrade and get on with what you wanted to do.

Jupyter fails to start
----------------------

* Have you `installed it <https://jupyter.org/install.html>`__? ;-)
* If you're using a menu shortcut or Anaconda launcher to start it, try
  opening a terminal or command prompt and running the command ``jupyter notebook``.
* If it can't find ``jupyter``,
  you may need to configure your ``PATH`` environment variable.
  If you don't know what that means, and don't want to find out,
  just (re)install Anaconda with the default settings,
  and it should set up PATH correctly.
* If Jupyter gives an error that it can't find ``notebook``,
  check with pip or conda that the ``notebook`` package is installed.
* Try running ``jupyter-notebook`` (with a hyphen). This should normally be the
  same as ``jupyter notebook`` (with a space), but if there's any difference,
  the version with the hyphen is the 'real' launcher, and the other one wraps
  that.

Jupyter doesn't load or doesn't work in the browser
---------------------------------------------------

* Try in another browser (e.g. if you normally use Firefox, try with Chrome).
  This helps pin down where the problem is.
* Try disabling any browser extensions and/or any Jupyter extensions you have
  installed.
* Some internet security software can interfere with Jupyter.
  If you have security software, try turning it off temporarily,
  and look in the settings for a more long-term solution.
* In the address bar, try changing between ``localhost`` and ``127.0.0.1``.
  They should be the same, but in some cases it makes a difference.

Jupyter can't start a kernel
----------------------------

Files called *kernel specs* tell Jupyter how to start different kinds of kernels.
To see where these are on your system, run ``jupyter kernelspec list``::

    $ jupyter kernelspec list
    Available kernels:
      python3      /home/takluyver/.local/lib/python3.6/site-packages/ipykernel/resources
      bash         /home/takluyver/.local/share/jupyter/kernels/bash
      ir           /home/takluyver/.local/share/jupyter/kernels/ir

There's a special fallback for the Python kernel:
if it doesn't find a real kernelspec, but it can import the ``ipykernel`` package,
it provides a kernel which will run in the same Python environment as the notebook server.
A path ending in ``ipykernel/resources``, like in the example above,
is this default kernel.
The default often does what you want,
so if the ``python3`` kernelspec points somewhere else
and you can't start a Python kernel,
try deleting or renaming that kernelspec folder to expose the default.

If your problem is with another kernel, not the Python one we maintain,
you may need to look for support about that kernel.

Python Environments
-------------------
Multiple python environments, whether based on Anaconda or Python Virtual environments,
are often the source of reported issues.  In many cases, these issues stem from the
Notebook server running in one environment, while the kernel and/or its resources,
derive from another environment.  Indicators of this scenario include:

* ``import`` statements within code cells producing ``ImportError`` or ``ModuleNotFound`` exceptions.
* General kernel startup failures exhibited by nothing happening when attempting
  to execute a cell.

In these situations, take a close look at your environment structure and ensure all
packages required by your notebook's code are installed in the correct environment.
If you need to run the kernel from different environments than your Notebook
server, check out `IPython's documentation <https://ipython.readthedocs.io/en/stable/install/kernel_install.html#kernels-for-different-environments>`_
for using kernels from different environments as this is the recommended approach.
Anaconda's `nb_conda_kernels <https://github.com/Anaconda-Platform/nb_conda_kernels>`_
package might also be an option for you in these scenarios.

Another thing to check is the ``kernel.json`` file that will be located in the
aforementioned *kernel specs* directory identified by running ``jupyter kernelspec list``.
This file will contain an ``argv`` stanza that includes the actual command to run
when launching the kernel.  Oftentimes, when reinstalling python environments, a previous
``kernel.json`` will reference an python executable from an old or non-existent location.
As a result, it's always a good idea when encountering kernel startup issues to validate
the ``argv`` stanza to ensure all file references exist and are appropriate.

Windows Systems
---------------
Although Jupyter Notebook is primarily developed on the various flavors of the Unix
operating system it also supports Microsoft
Windows - which introduces its own set of commonly encountered issues,
particularly in the areas of security, process management and lower-level libraries.

pywin32 Issues
^^^^^^^^^^^^^^^^^^
The primary package for interacting with Windows' primitives is ``pywin32``.

* Issues surrounding the creation of the kernel's communication file utilize
  ``jupyter_core``'s ``secure_write()`` function.  This function ensures a file is
  created in which only the owner of the file has access.  If libraries like ``pywin32``
  are not properly installed, issues can arise when it's necessary to use the native
  Windows libraries.

  Here's a portion of such a traceback::

    File "c:\users\jovyan\python\myenv.venv\lib\site-packages\jupyter_core\paths.py", line 424, in secure_write
    win32_restrict_file_to_user(fname)
    File "c:\users\jovyan\python\myenv.venv\lib\site-packages\jupyter_core\paths.py", line 359, in win32_restrict_file_to_user
    import win32api
    ImportError: DLL load failed: The specified module could not be found.

* As noted earlier, the installation of ``pywin32`` can be problematic on Windows
  configurations.  When such an issue occurs, you may need to revisit how the environment
  was setup.  Pay careful attention to whether you're running the 32 or 64 bit versions
  of Windows and be sure to install appropriate packages for that environment.

  Here's a portion of such a traceback::

    File "C:\Users\jovyan\AppData\Roaming\Python\Python37\site-packages\jupyter_core\paths.py", line 435, in secure_write
    win32_restrict_file_to_user(fname)
    File "C:\Users\jovyan\AppData\Roaming\Python\Python37\site-packages\jupyter_core\paths.py", line 361, in win32_restrict_file_to_user
    import win32api
    ImportError: DLL load failed: %1 is not a valid Win32 application

Resolving pywin32 Issues
""""""""""""""""""""""""""""
  In this case, your ``pywin32`` module may not be installed correctly and the following
  should be attempted:
  ::

    pip install --upgrade pywin32

  or::

    conda install --force-reinstall pywin32

  followed by::

    python.exe Scripts/pywin32_postinstall.py -install

  where ``Scripts`` is located in the active Python's installation location.

* Another common failure specific to Windows environments is the location of various
  python commands.  On ``*nix`` systems, these typically reside in the ``bin`` directory
  of the active Python environment.  However, on Windows, these tend to reside in the
  ``Scripts`` folder - which is a sibling to ``bin``.  As a result, when encountering
  kernel startup issues, again, check the ``argv`` stanza and verify it's pointing to a
  valid file.  You may find that it's pointing in ``bin`` when ``Scripts`` is correct, or
  the referenced file does not include its ``.exe`` extension - typically resulting in
  ``FileNotFoundError`` exceptions.

This Worked An Hour Ago
-----------------------
The Jupyter stack is very complex and rightfully so, there's a lot going on.  On occassion
you might find the system working perfectly well, then, suddenly, you can't get past a
certain cell due to ``import`` failures.  In these situations, it's best to ask yourself
if any new python files were added to your notebook development area.

These issues are usually evident by carefully analyzing the traceback produced in
the notebook error or the Notebook server's command window.  In these cases, you'll typically
find the Python kernel code (from ``IPython`` and ``ipykernel``) performing *its* imports
and notice a file from your Notebook development error included in that traceback followed
by an ``AttributeError``::

    File "C:\Users\jovyan\anaconda3\lib\site-packages\ipykernel\connect.py", line 13, in
    from IPython.core.profiledir import ProfileDir
    File "C:\Users\jovyan\anaconda3\lib\site-packages\IPython_init.py", line 55, in
    from .core.application import Application
    ...
    File "C:\Users\jovyan\anaconda3\lib\site-packages\ipython_genutils\path.py", line 13, in
    import random
    File "C:\Users\jovyan\Desktop\Notebooks\random.py", line 4, in
    rand_set = random.sample(english_words_lower_set, 12)
    AttributeError: module 'random' has no attribute 'sample'

What has happened is that you have named a file that conflicts with an installed package
that is used by the kernel software and now introduces a conflict preventing the
kernel's startup.

**Resolution**: You'll need to rename your file.  A best practice would be to prefix or
*namespace* your files so as not to conflict with any python package.


Asking for help
---------------

As with any problem, try searching to see if someone has already found an answer.
If you can't find an existing answer, you can ask questions at:

* The `Jupyter Discourse Forum <https://discourse.jupyter.org/>`_
* The `jupyter-notebook tag on Stackoverflow <https://stackoverflow.com/questions/tagged/jupyter-notebook>`_
* Peruse the `jupyter/help repository on Github <https://github.com/jupyter/help>`_ (read-only)
* Or in an issue on another repository, if it's clear which component is
  responsible.  Typical repositories include:

    * `jupyter_core <https://github.com/jupyter/jupyter_core>`_ - ``secure_write()``
      and file path issues
    * `jupyter_client <https://github.com/jupyter/jupyter_core>`_ - kernel management
      issues found in Notebook server's command window.
    * `IPython <https://github.com/ipython/ipython>`_ and
      `ipykernel <https://github.com/ipython/ipykernel>`_ - kernel runtime issues
      typically found in Notebook server's command window and/or Notebook cell execution.

Gathering Information
^^^^^^^^^^^^^^^^^^^^^
Should you find that your problem warrants that an issue be opened in
`notebook <https://github.com/jupyter/notebook>`_ please don't forget to provide details
like the following:

* What error messages do you see (within your notebook and, more importantly, in
  the Notebook server's command window)?
* What platform are you on?
* How did you install Jupyter?
* What have you tried already?

The ``jupyter troubleshoot`` command collects a lot of information
about your installation, which can also be useful.

When providing textual information, it's most helpful if you can *scrape* the contents
into the issue rather than providing a screenshot.  This enables others to select
pieces of that content so they can search more efficiently and try to help.

Remember that it's not anyone's job to help you.
We want Jupyter to work for you,
but we can't always help everyone individually.
.. _development_faq:

Developer FAQ
=============

1. How do I install a prerelease version such as a beta or release candidate?

.. code-block:: bash

    python -m pip install notebook --pre --upgrade
User interface components
=========================

When opening bug reports or sending emails to the Jupyter mailing list, it is
useful to know the names of different UI components so that other developers
and users have an easier time helping you diagnose your problems. This section
will familiarize you with the names of UI elements within the Notebook and the
different Notebook modes.

Notebook Dashboard
-------------------

When you launch ``jupyter notebook`` the first page that you encounter is the
Notebook Dashboard.

.. image:: ./_static/images/jupyter-notebook-dashboard.png

Notebook Editor
---------------

Once you've selected a Notebook to edit, the Notebook will open in the Notebook
Editor.

.. image:: ./_static/images/jupyter-notebook-default.png

Interactive User Interface Tour of the Notebook
-----------------------------------------------

If you would like to learn more about the specific elements within the Notebook
Editor, you can go through the user interface tour by selecting *Help* in the
menubar then selecting *User Interface Tour*.

Edit Mode and Notebook Editor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When a cell is in edit mode, the Cell Mode Indicator will change to reflect
the cell's state. This state is indicated by a small pencil icon on the
top right of the interface. When the cell is in command mode, there is no
icon in that location.

.. image:: ./_static/images/jupyter-notebook-edit.png

File Editor
-----------

Now let's say that you've chosen to open a Markdown file instead of a Notebook
file whilst in the Notebook Dashboard. If so, the file will be opened in the
File Editor.

.. image:: ./_static/images/jupyter-file-editor.png
====================
The Jupyter Notebook
====================

* `Installation <https://jupyter.readthedocs.io/en/latest/install.html>`_
* `Starting the Notebook <https://jupyter.readthedocs.io/en/latest/running.html>`_

.. toctree::
   :maxdepth: 1
   :caption: User Documentation

   notebook
   ui_components
   examples/Notebook/examples_index.rst
   troubleshooting
   changelog
   comms

.. toctree::
   :maxdepth: 1
   :caption: Configuration

   config_overview
   config
   public_server
   security
   frontend_config
   examples/Notebook/Distributing Jupyter Extensions as Python Packages
   extending/index.rst

.. toctree::
   :maxdepth: 1
   :caption: Contributor Documentation

   contributing
   development_faq

.. toctree::
    :hidden:

    examples/Notebook/nbpackage/mynotebook.ipynb
    examples/Notebook/nbpackage/nbs/other.ipynb
.. _configuration-overview:

Configuration Overview
======================

Beyond the default configuration settings, you can configure a rich array of
options to suit your workflow. Here are areas that are commonly configured
when using Jupyter Notebook:

    - :ref:`Jupyter's common configuration system <configure_common>`
    - :ref:`Notebook server <configure_nbserver>`
    - :ref:`Notebook front-end client <configure_nbclient>`
    - :ref:`Notebook extensions <configure_nbextensions>`

Let's look at highlights of each area.

.. _configure_common:

Jupyter's Common Configuration system
-------------------------------------
Jupyter applications, from the Notebook to JupyterHub to nbgrader, share a
common configuration system. The process for creating a configuration file
and editing settings is similar for all the Jupyter applications.

    - `Jupyter’s Common Configuration Approach <https://jupyter.readthedocs.io/en/latest/use/config.html>`_
    - `Common Directories and File Locations <https://jupyter.readthedocs.io/en/latest/use/jupyter-directories.html>`_
    - `Language kernels <https://jupyter.readthedocs.io/en/latest/projects/kernels.html>`_
    - `traitlets <https://traitlets.readthedocs.io/en/latest/config.html#module-traitlets.config>`_
      provide a low-level architecture for configuration.

.. _configure_nbserver:

Notebook server
---------------
The Notebook server runs the language kernel and communicates with the
front-end Notebook client (i.e. the familiar notebook interface).

  - Configuring the Notebook server

      To create a ``jupyter_notebook_config.py`` file in the ``.jupyter``
      directory, with all the defaults commented out, use the following
      command::

            $ jupyter notebook --generate-config

        :ref:`Command line arguments for configuration <config>` settings are
        documented in the configuration file and the user documentation.

  - :ref:`Running a Notebook server <working_remotely>`
  - Related: `Configuring a language kernel <https://ipython.readthedocs.io/en/latest/install/kernel_install.html>`_
    to run in the Notebook server enables your server to run other languages, like R or Julia.

.. _configure_nbclient:

Notebook front-end client
-------------------------

.. toctree::
   :maxdepth: 2

   frontend_config

.. _configure_nbextensions:

Notebook extensions
-------------------
- `Distributing Jupyter Extensions as Python Packages <https://jupyter-notebook.readthedocs.io/en/latest/examples/Notebook/Distributing%20Jupyter%20Extensions%20as%20Python%20Packages.html#Distributing-Jupyter-Extensions-as-Python-Packages>`_
- `Extending the Notebook <https://jupyter-notebook.readthedocs.io/en/latest/extending/index.html>`_


:ref:`Security in Jupyter notebooks:  <notebook_security>` Since security
policies vary from organization to organization, we encourage you to
consult with your security team on settings that would be best for your use
cases. Our documentation offers some responsible security practices, and we
recommend becoming familiar with the practices.
Comms
=====

*Comms* allow custom messages between the frontend and the kernel. They are used,
for instance, in `ipywidgets <https://ipywidgets.readthedocs.io/en/latest/>`__ to
update widget state.

A comm consists of a pair of objects, in the kernel and the frontend, with an
automatically assigned unique ID. When one side sends a message, a callback on
the other side is triggered with that message data. Either side, the frontend
or kernel, can open or close the comm.

.. seealso::

    `Custom Messages <https://jupyter-client.readthedocs.io/en/latest/messaging.html#custom-messages>`__
      The messaging specification section on comms

Opening a comm from the kernel
------------------------------

First, the function to accept the comm must be available on the frontend. This
can either be specified in a `requirejs` module, or registered in a registry, for
example when an :doc:`extension <extending/frontend_extensions>` is loaded.
This example shows a frontend comm target registered in a registry:

.. code-block:: javascript

    Jupyter.notebook.kernel.comm_manager.register_target('my_comm_target',
        function(comm, msg) {
            // comm is the frontend comm instance
            // msg is the comm_open message, which can carry data

            // Register handlers for later messages:
            comm.on_msg(function(msg) {...});
            comm.on_close(function(msg) {...});
            comm.send({'foo': 0});
        });

Now that the frontend comm is registered, you can open the comm from the kernel:

.. code-block:: python

    from ipykernel.comm import Comm

    # Use comm to send a message from the kernel
    my_comm = Comm(target_name='my_comm_target', data={'foo': 1})
    my_comm.send({'foo': 2})

    # Add a callback for received messages.
    @my_comm.on_msg
    def _recv(msg):
        # Use msg['content']['data'] for the data in the message


This example uses the IPython kernel; it's up to each language kernel what API,
if any, it offers for using comms.

Opening a comm from the frontend
--------------------------------

This is very similar to above, but in reverse. First, a comm target must be
registered in the kernel. For instance, this may be done by code displaying
output: it will register a target in the kernel, and then display output
containing Javascript to connect to it.

.. code-block:: python

    def target_func(comm, open_msg):
        # comm is the kernel Comm instance
        # msg is the comm_open message

        # Register handler for later messages
        @comm.on_msg
        def _recv(msg):
            # Use msg['content']['data'] for the data in the message
            comm.send({'echo': msg['content']['data']})

        # Send data to the frontend on creation
        comm.send({'foo': 5})

    get_ipython().kernel.comm_manager.register_target('my_comm_target', target_func)

This example uses the IPython kernel again; this example will be different in
other kernels that support comms. Refer to the specific language kernel's
documentation for comms support.

And then open the comm from the frontend:

.. code-block:: javascript

    const comm = Jupyter.notebook.kernel.comm_manager.new_comm('my_comm_target', {'foo': 6})
    // Send data
    comm.send({'foo': 7})

    // Register a handler
    comm.on_msg(function(msg) {
        console.log(msg.content.data.foo);
    });
.. _working_remotely:

Running a notebook server
=========================


The :doc:`Jupyter notebook <notebook>` web application is based on a
server-client structure.  The notebook server uses a :ref:`two-process kernel
architecture <ipython:ipythonzmq>` based on ZeroMQ_, as well as Tornado_ for
serving HTTP requests.

.. note::
   By default, a notebook server runs locally at 127.0.0.1:8888
   and is accessible only from `localhost`. You may access the
   notebook server from the browser using `http://127.0.0.1:8888`.

This document describes how you can
:ref:`secure a notebook server <notebook_server_security>` and how to
:ref:`run it on a public interface <notebook_public_server>`.

.. important::

    **This is not the multi-user server you are looking for**. This document
    describes how you can run a public server with a single user. This should
    only be done by someone who wants remote access to their personal machine.
    Even so, doing this requires a thorough understanding of the set-ups
    limitations and security implications. If you allow multiple users to
    access a notebook server as it is described in this document, their
    commands may collide, clobber and overwrite each other.

    If you want a multi-user server, the official solution is  JupyterHub_.
    To use JupyterHub, you need a Unix server (typically Linux) running
    somewhere that is accessible to your users on a network. This may run over
    the public internet, but doing so introduces additional
    `security concerns <https://jupyterhub.readthedocs.io/en/latest/getting-started/security-basics.html>`_.



.. _ZeroMQ: http://zeromq.org

.. _Tornado: http://www.tornadoweb.org

.. _JupyterHub: https://jupyterhub.readthedocs.io/en/latest/

.. _notebook_server_security:

Securing a notebook server
--------------------------

You can protect your notebook server with a simple single password. As of notebook
5.0 this can be done automatically. To set up a password manually you can configure the
:attr:`NotebookApp.password` setting in :file:`jupyter_notebook_config.py`.


Prerequisite: A notebook configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Check to see if you have a notebook configuration file,
:file:`jupyter_notebook_config.py`. The default location for this file
is your Jupyter folder located in your home directory:

    - Windows: :file:`C:\\Users\\USERNAME\\.jupyter\\jupyter_notebook_config.py`
    - OS X: :file:`/Users/USERNAME/.jupyter/jupyter_notebook_config.py`
    - Linux: :file:`/home/USERNAME/.jupyter/jupyter_notebook_config.py`

If you don't already have a Jupyter folder, or if your Jupyter folder doesn't contain
a notebook configuration file, run the following command::

  $ jupyter notebook --generate-config

This command will create the Jupyter folder if necessary, and create notebook
configuration file, :file:`jupyter_notebook_config.py`, in this folder.


Automatic Password setup
~~~~~~~~~~~~~~~~~~~~~~~~

As of notebook 5.3, the first time you log-in using a token, the notebook server
should give you the opportunity to setup a password from the user interface.

You will be presented with a form asking for the current _token_, as well as
your _new_ _password_ ; enter both and click on ``Login and setup new password``.

Next time you need to log in you'll be able to use the new password instead of
the login token, otherwise follow the procedure to set a password from the
command line.

The ability to change the password at first login time may be disabled by
integrations by setting the ``--NotebookApp.allow_password_change=False``


Starting at notebook version 5.0, you can enter and store a password for your
notebook server with a single command. :command:`jupyter notebook password` will
prompt you for your password and record the hashed password in your
:file:`jupyter_notebook_config.json`.

.. code-block:: bash

    $ jupyter notebook password
    Enter password:  ****
    Verify password: ****
    [NotebookPasswordApp] Wrote hashed password to /Users/you/.jupyter/jupyter_notebook_config.json

This can be used to reset a lost password; or if you believe your credentials
have been leaked and desire to change your password. Changing your password will
invalidate all logged-in sessions after a server restart.

.. _hashed-pw:

Preparing a hashed password
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can prepare a hashed password manually, using the function
:func:`notebook.auth.security.passwd`:

.. code-block:: ipython

    In [1]: from notebook.auth import passwd
    In [2]: passwd()
    Enter password:
    Verify password:
    Out[2]: 'sha1:67c9e60bb8b6:9ffede0825894254b2e042ea597d771089e11aed'

.. caution::

  :func:`~notebook.auth.security.passwd` when called with no arguments
  will prompt you to enter and verify your password such as
  in the above code snippet. Although the function can also
  be passed a string as an argument such as ``passwd('mypassword')``, please
  **do not** pass a string as an argument inside an IPython session, as it
  will be saved in your input history.

Adding hashed password to your notebook configuration file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
You can then add the hashed password to your
:file:`jupyter_notebook_config.py`. The default location for this file
:file:`jupyter_notebook_config.py` is in your Jupyter folder in your home
directory, ``~/.jupyter``, e.g.::

    c.NotebookApp.password = u'sha1:67c9e60bb8b6:9ffede0825894254b2e042ea597d771089e11aed'

Automatic password setup will store the hash in ``jupyter_notebook_config.json``
while this method stores the hash in ``jupyter_notebook_config.py``. The ``.json``
configuration options take precedence over the ``.py`` one, thus the manual
password may not take effect if the Json file has a password set.


Using SSL for encrypted communication
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
When using a password, it is a good idea to also use SSL with a web
certificate, so that your hashed password is not sent unencrypted by your
browser.

.. important::
   Web security is rapidly changing and evolving. We provide this document
   as a convenience to the user, and recommend that the user keep current on
   changes that may impact security, such as new releases of OpenSSL.
   The Open Web Application Security Project (`OWASP`_) website is a good resource
   on general security issues and web practices.

You can start the notebook to communicate via a secure protocol mode by setting
the ``certfile`` option to your self-signed certificate, i.e. ``mycert.pem``,
with the command::

    $ jupyter notebook --certfile=mycert.pem --keyfile mykey.key

.. tip::

    A self-signed certificate can be generated with ``openssl``.  For example,
    the following command will create a certificate valid for 365 days with
    both the key and certificate data written to the same file::

        $ openssl req -x509 -nodes -days 365 -newkey rsa:2048 -keyout mykey.key -out mycert.pem

When starting the notebook server, your browser may warn that your self-signed
certificate is insecure or unrecognized.  If you wish to have a fully
compliant self-signed certificate that will not raise warnings, it is possible
(but rather involved) to create one, as explained in detail in this
`tutorial`_. Alternatively, you may use `Let's Encrypt`_ to acquire a free SSL
certificate and follow the steps in :ref:`using-lets-encrypt` to set up a
public server.

.. _OWASP: https://www.owasp.org
.. _tutorial: https://arstechnica.com/information-technology/2009/12/how-to-get-set-with-a-secure-sertificate-for-free/

.. _notebook_public_server:

Running a public notebook server
--------------------------------

If you want to access your notebook server remotely via a web browser,
you can do so by running a public notebook server. For optimal security
when running a public notebook server, you should first secure the
server with a password and SSL/HTTPS as described in
:ref:`notebook_server_security`.

Start by creating a certificate file and a hashed password, as explained in
:ref:`notebook_server_security`.

If you don't already have one, create a
config file for the notebook using the following command line::

  $ jupyter notebook --generate-config

In the ``~/.jupyter`` directory, edit the notebook config file,
``jupyter_notebook_config.py``.  By default, the notebook config file has
all fields commented out. The minimum set of configuration options that
you should uncomment and edit in :file:`jupyter_notebook_config.py` is the
following::

     # Set options for certfile, ip, password, and toggle off
     # browser auto-opening
     c.NotebookApp.certfile = u'/absolute/path/to/your/certificate/mycert.pem'
     c.NotebookApp.keyfile = u'/absolute/path/to/your/certificate/mykey.key'
     # Set ip to '*' to bind on all interfaces (ips) for the public server
     c.NotebookApp.ip = '*'
     c.NotebookApp.password = u'sha1:bcd259ccf...<your hashed password here>'
     c.NotebookApp.open_browser = False

     # It is a good idea to set a known, fixed port for server access
     c.NotebookApp.port = 9999

You can then start the notebook using the ``jupyter notebook`` command.

.. _using-lets-encrypt:

Using Let's Encrypt
~~~~~~~~~~~~~~~~~~~
`Let's Encrypt`_ provides free SSL/TLS certificates. You can also set up a
public server using a `Let's Encrypt`_ certificate.

:ref:`notebook_public_server` will be similar when using a Let's Encrypt
certificate with a few configuration changes. Here are the steps:

1. Create a `Let's Encrypt certificate <https://letsencrypt.org/getting-started/>`_.
2. Use :ref:`hashed-pw` to create one.
3. If you don't already have config file for the notebook, create one
   using the following command:

   .. code-block:: bash

       $ jupyter notebook --generate-config

4. In the ``~/.jupyter`` directory, edit the notebook config file,
``jupyter_notebook_config.py``.  By default, the notebook config file has
all fields commented out. The minimum set of configuration options that
you should to uncomment and edit in :file:`jupyter_notebook_config.py` is the
following::

     # Set options for certfile, ip, password, and toggle off
     # browser auto-opening
     c.NotebookApp.certfile = u'/absolute/path/to/your/certificate/fullchain.pem'
     c.NotebookApp.keyfile = u'/absolute/path/to/your/certificate/privkey.pem'
     # Set ip to '*' to bind on all interfaces (ips) for the public server
     c.NotebookApp.ip = '*'
     c.NotebookApp.password = u'sha1:bcd259ccf...<your hashed password here>'
     c.NotebookApp.open_browser = False

     # It is a good idea to set a known, fixed port for server access
     c.NotebookApp.port = 9999

You can then start the notebook using the ``jupyter notebook`` command.

.. important::

    **Use 'https'.**
    Keep in mind that when you enable SSL support, you must access the
    notebook server over ``https://``, not over plain ``http://``.  The startup
    message from the server prints a reminder in the console, but *it is easy
    to overlook this detail and think the server is for some reason
    non-responsive*.

    **When using SSL, always access the notebook server with 'https://'.**

You may now access the public server by pointing your browser to
``https://your.host.com:9999`` where ``your.host.com`` is your public server's
domain.

.. _`Let's Encrypt`: https://letsencrypt.org


Firewall Setup
~~~~~~~~~~~~~~

To function correctly, the firewall on the computer running the jupyter
notebook server must be configured to allow connections from client
machines on the access port ``c.NotebookApp.port`` set in
:file:`jupyter_notebook_config.py` to allow connections to the
web interface.  The firewall must also allow connections from
127.0.0.1 (localhost) on ports from 49152 to 65535.
These ports are used by the server to communicate with the notebook kernels.
The kernel communication ports are chosen randomly by ZeroMQ, and may require
multiple connections per kernel, so a large range of ports must be accessible.

Running the notebook with a customized URL prefix
-------------------------------------------------

The notebook dashboard, which is the landing page with an overview
of the notebooks in your working directory, is typically found and accessed
at the default URL ``http://localhost:8888/``.

If you prefer to customize the URL prefix for the notebook dashboard, you can
do so through modifying ``jupyter_notebook_config.py``. For example, if you
prefer that the notebook dashboard be located with a sub-directory that
contains other ipython files, e.g. ``http://localhost:8888/ipython/``,
you can do so with configuration options like the following (see above for
instructions about modifying ``jupyter_notebook_config.py``):

.. code-block:: python

    c.NotebookApp.base_url = '/ipython/'


Embedding the notebook in another website
-----------------------------------------

Sometimes you may want to embed the notebook somewhere on your website,
e.g. in an IFrame. To do this, you may need to override the
Content-Security-Policy to allow embedding. Assuming your website is at
`https://mywebsite.example.com`, you can embed the notebook on your website
with the following configuration setting in
:file:`jupyter_notebook_config.py`:

.. code-block:: python

    c.NotebookApp.tornado_settings = {
        'headers': {
            'Content-Security-Policy': "frame-ancestors https://mywebsite.example.com 'self' "
        }
    }

When embedding the notebook in a website using an iframe,
consider putting the notebook in single-tab mode.
Since the notebook opens some links in new tabs by default,
single-tab mode keeps the notebook from opening additional tabs.
Adding the following to :file:`~/.jupyter/custom/custom.js` will enable
single-tab mode:

.. code-block:: javascript

    define(['base/js/namespace'], function(Jupyter){
        Jupyter._target = '_self';
    });


Using a gateway server for kernel management
--------------------------------------------

You are now able to redirect the management of your kernels to a Gateway Server
(i.e., `Jupyter Kernel Gateway <https://jupyter-kernel-gateway.readthedocs.io/en/latest/>`_ or
`Jupyter Enterprise Gateway <https://jupyter-enterprise-gateway.readthedocs.io/en/latest/>`_)
simply by specifying a Gateway url via the following command-line option:

    .. code-block:: bash

        $ jupyter notebook --gateway-url=http://my-gateway-server:8888

the environment:

    .. code-block:: bash

        JUPYTER_GATEWAY_URL=http://my-gateway-server:8888

or in :file:`jupyter_notebook_config.py`:

   .. code-block:: python

      c.GatewayClient.url = http://my-gateway-server:8888

When provided, all kernel specifications will be retrieved from the specified Gateway server and all
kernels will be managed by that server.  This option enables the ability to target kernel processes
against managed clusters while allowing for the notebook's management to remain local to the Notebook
server.

Known issues
------------

Proxies
~~~~~~~

When behind a proxy, especially if your system or browser is set to autodetect
the proxy, the notebook web application might fail to connect to the server's
websockets, and present you with a warning at startup. In this case, you need
to configure your system not to use the proxy for the server's address.

For example, in Firefox, go to the Preferences panel, Advanced section,
Network tab, click 'Settings...', and add the address of the notebook server
to the 'No proxy for' field.

Content-Security-Policy (CSP)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Certain `security guidelines
<https://infosec.mozilla.org/guidelines/web_security.html#content-security-policy>`_
recommend that servers use a Content-Security-Policy (CSP) header to prevent
cross-site scripting vulnerabilities, specifically limiting to ``default-src:
https:`` when possible.  This directive causes two problems with Jupyter.
First, it disables execution of inline javascript code, which is used
extensively by Jupyter.  Second, it limits communication to the https scheme,
and prevents WebSockets from working because they communicate via the wss
scheme (or ws for insecure communication).  Jupyter uses WebSockets for
interacting with kernels, so when you visit a server with such a CSP, your
browser will block attempts to use wss, which will cause you to see
"Connection failed" messages from jupyter notebooks, or simply no response
from jupyter terminals.  By looking in your browser's javascript console, you
can see any error messages that will explain what is failing.

To avoid these problem, you need to add ``'unsafe-inline'`` and ``connect-src
https: wss:`` to your CSP header, at least for pages served by jupyter.  (That
is, you can leave your CSP unchanged for other parts of your website.)  Note
that multiple CSP headers are allowed, but successive CSP headers can only
restrict the policy; they cannot loosen it.  For example, if your server sends
both of these headers

    Content-Security-Policy "default-src https: 'unsafe-inline'"
    Content-Security-Policy "connect-src https: wss:"

the first policy will already eliminate wss connections, so the second has no
effect.  Therefore, you can't simply add the second header; you have to
actually modify your CSP header to look more like this:

    Content-Security-Policy "default-src https: 'unsafe-inline'; connect-src https: wss:"



Docker CMD
~~~~~~~~~~

Using ``jupyter notebook`` as a
`Docker CMD <https://docs.docker.com/engine/reference/builder/#cmd>`_ results in
kernels repeatedly crashing, likely due to a lack of `PID reaping
<https://blog.phusion.nl/2015/01/20/docker-and-the-pid-1-zombie-reaping-problem/>`_.
To avoid this, use the `tini <https://github.com/krallin/tini>`_ ``init`` as your
Dockerfile `ENTRYPOINT`::

  # Add Tini. Tini operates as a process subreaper for jupyter. This prevents
  # kernel crashes.
  ENV TINI_VERSION v0.6.0
  ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
  RUN chmod +x /usr/bin/tini
  ENTRYPOINT ["/usr/bin/tini", "--"]

  EXPOSE 8888
  CMD ["jupyter", "notebook", "--port=8888", "--no-browser", "--ip=0.0.0.0"]
File save hooks
===============

You can configure functions that are run whenever a file is saved. There are
two hooks available:

* ``ContentsManager.pre_save_hook`` runs on the API path and model with
  content. This can be used for things like stripping output that people don't
  like adding to VCS noise.
* ``FileContentsManager.post_save_hook`` runs on the filesystem path and model
  without content. This could be used to commit changes after every save, for
  instance.

They are both called with keyword arguments:

.. code-block:: python

    pre_save_hook(model=model, path=path, contents_manager=cm)
    post_save_hook(model=model, os_path=os_path, contents_manager=cm)

Examples
--------

These can both be added to :file:`jupyter_notebook_config.py`.

A pre-save hook for stripping output:

.. code-block:: python

    def scrub_output_pre_save(model, **kwargs):
        """scrub output before saving notebooks"""
        # only run on notebooks
        if model['type'] != 'notebook':
            return
        # only run on nbformat v4
        if model['content']['nbformat'] != 4:
            return

        for cell in model['content']['cells']:
            if cell['cell_type'] != 'code':
                continue
            cell['outputs'] = []
            cell['execution_count'] = None

    c.FileContentsManager.pre_save_hook = scrub_output_pre_save

A post-save hook to make a script equivalent whenever the notebook is saved
(replacing the ``--script`` option in older versions of the notebook):

.. code-block:: python

    import io
    import os
    from notebook.utils import to_api_path

    _script_exporter = None

    def script_post_save(model, os_path, contents_manager, **kwargs):
        """convert notebooks to Python script after save with nbconvert

        replaces `jupyter notebook --script`
        """
        from nbconvert.exporters.script import ScriptExporter

        if model['type'] != 'notebook':
            return

        global _script_exporter

        if _script_exporter is None:
            _script_exporter = ScriptExporter(parent=contents_manager)

        log = contents_manager.log

        base, ext = os.path.splitext(os_path)
        script, resources = _script_exporter.from_filename(os_path)
        script_fname = base + resources.get('output_extension', '.txt')
        log.info("Saving script /%s", to_api_path(script_fname, contents_manager.root_dir))

        with io.open(script_fname, 'w', encoding='utf-8') as f:
            f.write(script)

    c.FileContentsManager.post_save_hook = script_post_save


This could be a simple call to ``jupyter nbconvert --to script``, but spawning
the subprocess every time is quite slow.
Custom request handlers
=======================

The notebook webserver can be interacted with using a well `defined
RESTful
API <http://petstore.swagger.io/?url=https://raw.githubusercontent.com/jupyter/notebook/master/notebook/services/api/api.yaml>`__.
You can define custom RESTful API handlers in addition to the ones
provided by the notebook. As described below, to define a custom handler
you need to first write a notebook server extension. Then, in the
extension, you can register the custom handler.

Writing a notebook server extension
-----------------------------------

The notebook webserver is written in Python, hence your server extension
should be written in Python too. Server extensions, like IPython
extensions, are Python modules that define a specially named load
function, ``load_jupyter_server_extension``. This function is called
when the extension is loaded.

.. code:: python

    def load_jupyter_server_extension(nb_server_app):
        """
        Called when the extension is loaded.

        Args:
            nb_server_app (NotebookWebApplication): handle to the Notebook webserver instance.
        """
        pass

To get the notebook server to load your custom extension, you'll need to
add it to the list of extensions to be loaded. You can do this using the
config system. ``NotebookApp.nbserver_extensions`` is a config variable
which is a dictionary of strings, each a Python module to be imported, mapping
to ``True`` to enable or ``False`` to disable each extension.
Because this variable is notebook config, you can set it two different
ways, using config files or via the command line.

For example, to get your extension to load via the command line add a
double dash before the variable name, and put the Python dictionary in
double quotes. If your package is "mypackage" and module is
"mymodule", this would look like
``jupyter notebook --NotebookApp.nbserver_extensions="{'mypackage.mymodule':True}"``
.
Basically the string should be Python importable.

Alternatively, you can have your extension loaded regardless of the
command line args by setting the variable in the Jupyter config file.
The default location of the Jupyter config file is
``~/.jupyter/jupyter_notebook_config.py`` (see :doc:`/config_overview`). Inside
the config file, you can use Python to set the variable. For example,
the following config does the same as the previous command line example.

.. code:: python

    c = get_config()
    c.NotebookApp.nbserver_extensions = {
        'mypackage.mymodule': True,
    }

Before continuing, it's a good idea to verify that your extension is
being loaded. Use a print statement to print something unique. Launch
the notebook server and you should see your statement printed to the
console.

Registering custom handlers
---------------------------

Once you've defined a server extension, you can register custom handlers
because you have a handle to the Notebook server app instance
(``nb_server_app`` above). However, you first need to define your custom
handler. To declare a custom handler, inherit from
``notebook.base.handlers.IPythonHandler``. The example below[1] is a
Hello World handler:

.. code:: python

    from notebook.base.handlers import IPythonHandler

    class HelloWorldHandler(IPythonHandler):
        def get(self):
            self.finish('Hello, world!')

The Jupyter Notebook server use
`Tornado <http://www.tornadoweb.org/en/stable/>`__ as its web framework.
For more information on how to implement request handlers, refer to the
`Tornado documentation on the
matter <http://www.tornadoweb.org/en/stable/web.html#request-handlers>`__.

After defining the handler, you need to register the handler with the
Notebook server. See the following example:

.. code:: python

    web_app = nb_server_app.web_app
    host_pattern = '.*$'
    route_pattern = url_path_join(web_app.settings['base_url'], '/hello')
    web_app.add_handlers(host_pattern, [(route_pattern, HelloWorldHandler)])

Putting this together with the extension code, the example looks like the
following:

.. code:: python

    from notebook.utils import url_path_join
    from notebook.base.handlers import IPythonHandler

    class HelloWorldHandler(IPythonHandler):
        def get(self):
            self.finish('Hello, world!')

    def load_jupyter_server_extension(nb_server_app):
        """
        Called when the extension is loaded.

        Args:
            nb_server_app (NotebookWebApplication): handle to the Notebook webserver instance.
        """
        web_app = nb_server_app.web_app
        host_pattern = '.*$'
        route_pattern = url_path_join(web_app.settings['base_url'], '/hello')
        web_app.add_handlers(host_pattern, [(route_pattern, HelloWorldHandler)])


Extra Parameters and authentication
===================================

Here is a quick rundown of what you need to know to pass extra parameters to the handler and enable authentication:

 - extra arguments to the ``__init__`` constructor are given in a dictionary after the  handler class in ``add_handlers``:

.. code:: python


    class HelloWorldHandler(IPythonHandler):

        def __init__(self, *args, **kwargs):
            self.extra = kwargs.pop('extra')
            ...

    def load_jupyter_server_extension(nb_server_app):

        ...

        web_app.add_handlers(host_pattern,
            [
               (route_pattern, HelloWorldHandler, {"extra": nb_server_app.extra})
            ])


All handler methods that require authentication _MUST_ be decorated with ``@tornado.web.authenticated``:


.. code:: python

    from tornado import web

    class HelloWorldHandler(IPythonHandler):

        ...

        @web.authenticated
        def  get(self, *args, **kwargs):
             ...

        @web.authenticated
        def  post(self, *args, **kwargs):
             ...


References:

1. `Peter Parente's Mindtrove <https://mindtrove.info/4-ways-to-extend-jupyter-notebook/#nb-server-exts>`__
Custom front-end extensions
===========================

This describes the basic steps to write a JavaScript extension for the Jupyter
notebook front-end. This allows you to customize the behaviour of the various
pages like the dashboard, the notebook, or the text editor.

The structure of a front-end extension
--------------------------------------

.. note::

    The notebook front-end and Javascript API are not stable, and are subject
    to a lot of changes. Any extension written for the current notebook is
    almost guaranteed to break in the next release.

.. _AMD module: https://en.wikipedia.org/wiki/Asynchronous_module_definition

A front-end extension is a JavaScript file that defines an `AMD module`_
which exposes at least a function called ``load_ipython_extension``, which
takes no arguments. We will not get into the details of what each of these
terms consists of yet, but here is the minimal code needed for a working
extension:

.. code:: javascript

    // file my_extension/main.js

    define(function(){

        function load_ipython_extension(){
            console.info('this is my first extension');
        }

        return {
            load_ipython_extension: load_ipython_extension
        };
    });

.. note::

    Although for historical reasons the function is called
    ``load_ipython_extension``, it does apply to the Jupyter notebook in
    general, and will work regardless of the kernel in use.

If you are familiar with JavaScript, you can use this template to require any
Jupyter module and modify its configuration, or do anything else in client-side
Javascript. Your extension will be loaded at the right time during the notebook
page initialisation for you to set up a listener for the various events that
the page can trigger.

You might want access to the current instances of the various Jupyter notebook
components on the page, as opposed to the classes defined in the modules. The
current instances are exposed by a module named ``base/js/namespace``. If you
plan on accessing instances on the page, you should ``require`` this module
rather than accessing the global variable ``Jupyter``, which will be removed in
future. The following example demonstrates how to access the current notebook
instance:

.. code:: javascript

    // file my_extension/main.js

    define([
        'base/js/namespace'
    ], function(
        Jupyter
    ) {
        function load_ipython_extension() {
            console.log(
                'This is the current notebook application instance:',
                Jupyter.notebook
            );
        }

        return {
            load_ipython_extension: load_ipython_extension
        };
    });


Modifying key bindings
----------------------

One of the abilities of extensions is to modify key bindings, although once
again this is an API which is not guaranteed to be stable. However, custom key
bindings are frequently requested, and are helpful to increase accessibility,
so in the following we show how to access them.

Here is an example of an extension that will unbind the shortcut ``0,0`` in
command mode, which normally restarts the kernel, and bind ``0,0,0`` in its
place:

.. code:: javascript

    // file my_extension/main.js

    define([
        'base/js/namespace'
    ], function(
        Jupyter
    ) {

        function load_ipython_extension() {
            Jupyter.keyboard_manager.command_shortcuts.remove_shortcut('0,0');
            Jupyter.keyboard_manager.command_shortcuts.add_shortcut('0,0,0', 'jupyter-notebook:restart-kernel');
        }

        return {
            load_ipython_extension: load_ipython_extension
        };
    });

.. note::

    The standard keybindings might not work correctly on non-US keyboards.
    Unfortunately, this is a limitation of browser implementations and the
    status of keyboard event handling on the web in general. We appreciate your
    feedback if you have issues binding keys, or have any ideas to help improve
    the situation.

You can see that I have used the **action name**
``jupyter-notebook:restart-kernel`` to bind the new shortcut. There is no API
yet to access the list of all available *actions*, though the following in the
JavaScript console of your browser on a notebook page should give you an idea
of what is available:

.. code:: javascript

    Object.keys(require('base/js/namespace').actions._actions);

In this example, we changed a keyboard shortcut in **command mode**; you
can also customize keyboard shortcuts in **edit mode**.
However, most of the keyboard shortcuts in edit mode are handled by CodeMirror,
which supports custom key bindings via a completely different API.


Defining and registering your own actions
-----------------------------------------

As part of your front-end extension, you may wish to define actions, which can
be attached to toolbar buttons, or called from the command palette. Here is an
example of an extension that defines an (not very useful!) action to show an
alert, and adds a toolbar button using the full action name:

.. code:: javascript

    // file my_extension/main.js

    define([
        'base/js/namespace'
    ], function(
        Jupyter
    ) {
        function load_ipython_extension() {

            var handler = function () {
                alert('this is an alert from my_extension!');
            };

            var action = {
                icon: 'fa-comment-o', // a font-awesome class used on buttons, etc
                help    : 'Show an alert',
                help_index : 'zz',
                handler : handler
            };
            var prefix = 'my_extension';
            var action_name = 'show-alert';

            var full_action_name = Jupyter.actions.register(action, action_name, prefix); // returns 'my_extension:show-alert'
            Jupyter.toolbar.add_buttons_group([full_action_name]);
        }

        return {
            load_ipython_extension: load_ipython_extension
        };
    });

Every action needs a name, which, when joined with its prefix to make the full
action name, should be unique. Built-in actions, like the
``jupyter-notebook:restart-kernel`` we bound in the earlier
`Modifying key bindings`_ example, use the prefix ``jupyter-notebook``. For
actions defined in an extension, it makes sense to use the extension name as
the prefix. For the action name, the following guidelines should be considered:

.. adapted from notebook/static/notebook/js/actions.js

* First pick a noun and a verb for the action. For example, if the action is
  "restart kernel," the verb is "restart" and the noun is "kernel".
* Omit terms like "selected" and "active" by default, so "delete-cell", rather
  than "delete-selected-cell". Only provide a scope like "-all-" if it is other
  than the default "selected" or "active" scope.
* If an action has a secondary action, separate the secondary action with
  "-and-", so "restart-kernel-and-clear-output".
* Use above/below or previous/next to indicate spatial and sequential
  relationships.
* Don't ever use before/after as they have a temporal connotation that is
  confusing when used in a spatial context.
* For dialogs, use a verb that indicates what the dialog will accomplish, such
  as "confirm-restart-kernel".


Installing and enabling extensions
----------------------------------

You can install your nbextension with the command::

    jupyter nbextension install path/to/my_extension/ [--user|--sys-prefix]

The default installation is system-wide. You can use ``--user`` to do a
per-user installation, or ``--sys-prefix`` to install to Python's prefix (e.g.
in a virtual or conda environment). Where my_extension is the directory
containing the Javascript files. This will copy it to a Jupyter data directory
(the exact location is platform dependent - see :ref:`jupyter_path`).

For development, you can use the ``--symlink`` flag to symlink your extension
rather than copying it, so there's no need to reinstall after changes.

To use your extension, you'll also need to **enable** it, which tells the
notebook interface to load it. You can do that with another command::

    jupyter nbextension enable my_extension/main [--sys-prefix][--section='common']

The argument refers to the Javascript module containing your
``load_ipython_extension`` function, which is ``my_extension/main.js`` in this
example. The ``--section='common'`` argument will affect all pages, by default 
it will be loaded on the notebook view only. 
There is a corresponding ``disable`` command to stop using an
extension without uninstalling it.

.. versionchanged:: 4.2

    Added ``--sys-prefix`` argument


Kernel Specific extensions
--------------------------

.. warning::

  This feature serves as a stopgap for kernel developers who need specific
  JavaScript injected onto the page. The availability and API are subject to
  change at anytime.


It is possible to load some JavaScript on the page on a per kernel basis. Be
aware that doing so will make the browser page reload without warning as
soon as the user switches the kernel without notice.

If you, a kernel developer, need a particular piece of JavaScript to be loaded
on a "per kernel" basis, such as:

* if you are developing a CodeMirror mode for your language
* if you need to enable some specific debugging options

your ``kernelspecs`` are allowed to contain a ``kernel.js`` file that defines
an AMD module. The AMD module should define an `onload` function that will be
called when the kernelspec loads, such as:

* when you load a notebook that uses your kernelspec
* change the active kernelspec of a notebook to your kernelspec.

Note that adding a `kernel.js` to your kernelspec will add an unexpected side
effect to changing a kernel in the notebook. As it is impossible to "unload"
JavaScript, any attempt to change the kernelspec again will save the current
notebook and reload the page without confirmations.

Here is an example of ``kernel.js``:

.. code:: javascript

    define(function(){
      return {onload: function(){
        console.info('Kernel specific javascript loaded');

        // do more things here, like define a codemirror mode

      }}

    });
.. _contents_api:

Contents API
============

.. currentmodule:: notebook.services.contents

The Jupyter Notebook web application provides a graphical interface for
creating, opening, renaming, and deleting files in a virtual filesystem.

The :class:`~manager.ContentsManager` class defines an abstract
API for translating these interactions into operations on a particular storage
medium. The default implementation,
:class:`~filemanager.FileContentsManager`, uses the local
filesystem of the server for storage and straightforwardly serializes notebooks
into JSON.  Users can override these behaviors by supplying custom subclasses
of ContentsManager.

This section describes the interface implemented by ContentsManager subclasses.
We refer to this interface as the **Contents API**.

Data Model
----------

.. currentmodule:: notebook.services.contents.manager

Filesystem Entities
~~~~~~~~~~~~~~~~~~~
.. _notebook models:

ContentsManager methods represent virtual filesystem entities as dictionaries,
which we refer to as **models**.

Models may contain the following entries:

+--------------------+-----------+------------------------------+
| Key                | Type      |Info                          |
+====================+===========+==============================+
|**name**            |unicode    |Basename of the entity.       |
+--------------------+-----------+------------------------------+
|**path**            |unicode    |Full                          |
|                    |           |(:ref:`API-style<apipaths>`)  |
|                    |           |path to the entity.           |
+--------------------+-----------+------------------------------+
|**type**            |unicode    |The entity type. One of       |
|                    |           |``"notebook"``, ``"file"`` or |
|                    |           |``"directory"``.              |
+--------------------+-----------+------------------------------+
|**created**         |datetime   |Creation date of the entity.  |
+--------------------+-----------+------------------------------+
|**last_modified**   |datetime   |Last modified date of the     |
|                    |           |entity.                       |
+--------------------+-----------+------------------------------+
|**content**         |variable   |The "content" of the entity.  |
|                    |           |(:ref:`See                    |
|                    |           |Below<modelcontent>`)         |
+--------------------+-----------+------------------------------+
|**mimetype**        |unicode or |The mimetype of ``content``,  |
|                    |``None``   |if any.  (:ref:`See           |
|                    |           |Below<modelcontent>`)         |
+--------------------+-----------+------------------------------+
|**format**          |unicode or |The format of ``content``,    |
|                    |``None``   |if any. (:ref:`See            |
|                    |           |Below<modelcontent>`)         |
+--------------------+-----------+------------------------------+

.. _modelcontent:

Certain model fields vary in structure depending on the ``type`` field of the
model. There are three model types: **notebook**, **file**, and **directory**.

- ``notebook`` models
    - The ``format`` field is always ``"json"``.
    - The ``mimetype`` field is always ``None``.
    - The ``content`` field contains a
      :class:`nbformat.notebooknode.NotebookNode` representing the .ipynb file
      represented by the model.  See the `NBFormat`_ documentation for a full
      description.

- ``file`` models
    - The ``format`` field is either ``"text"`` or ``"base64"``.
    - The ``mimetype`` field can be any mimetype string, but defaults to 
      ``text/plain`` for text-format models and
      ``application/octet-stream`` for base64-format models. For files with
      unknown mime types (e.g. unknown file extensions), this field may be
      `None`.
    - The ``content`` field is always of type ``unicode``.  For text-format
      file models, ``content`` simply contains the file's bytes after decoding
      as UTF-8.  Non-text (``base64``) files are read as bytes, base64 encoded,
      and then decoded as UTF-8.

- ``directory`` models
    - The ``format`` field is always ``"json"``.
    - The ``mimetype`` field is always ``None``.
    - The ``content`` field contains a list of :ref:`content-free<contentfree>`
      models representing the entities in the directory.

.. note::

   .. _contentfree:

   In certain circumstances, we don't need the full content of an entity to
   complete a Contents API request. In such cases, we omit the ``content``, and
   ``format`` keys from the model. The default values for the ``mimetype``
   field will might also not be evaluated, in which case it will be set as `None`.
   This reduced reply most commonly occurs when listing a directory, in
   which circumstance we represent files within the directory as content-less
   models to avoid having to recursively traverse and serialize the entire
   filesystem.

**Sample Models**

.. code-block:: python

    # Notebook Model with Content
    {
        'content': {
            'metadata': {},
            'nbformat': 4,
            'nbformat_minor': 0,
            'cells': [
                {
                    'cell_type': 'markdown',
                    'metadata': {},
                    'source': 'Some **Markdown**',
                },
            ],
        },
        'created': datetime(2015, 7, 25, 19, 50, 19, 19865),
        'format': 'json',
        'last_modified': datetime(2015, 7, 25, 19, 50, 19, 19865),
        'mimetype': None,
        'name': 'a.ipynb',
        'path': 'foo/a.ipynb',
        'type': 'notebook',
        'writable': True,
    }

    # Notebook Model without Content
    {
        'content': None,
        'created': datetime.datetime(2015, 7, 25, 20, 17, 33, 271931),
        'format': None,
        'last_modified': datetime.datetime(2015, 7, 25, 20, 17, 33, 271931),
        'mimetype': None,
        'name': 'a.ipynb',
        'path': 'foo/a.ipynb',
        'type': 'notebook',
        'writable': True
    }


API Paths
~~~~~~~~~
.. _apipaths:

ContentsManager methods represent the locations of filesystem resources as
**API-style paths**.  Such paths are interpreted as relative to the root
directory of the notebook server.  For compatibility across systems, the
following guarantees are made:

* Paths are always ``unicode``, not ``bytes``.
* Paths are not URL-escaped.
* Paths are always forward-slash (/) delimited, even on Windows.
* Leading and trailing slashes are stripped.  For example, ``/foo/bar/buzz/``
  becomes ``foo/bar/buzz``.
* The empty string (``""``) represents the root directory.


Writing a Custom ContentsManager
--------------------------------

The default ContentsManager is designed for users running the notebook as an
application on a personal computer.  It stores notebooks as .ipynb files on the
local filesystem, and it maps files and directories in the Notebook UI to files
and directories on disk.  It is possible to override how notebooks are stored
by implementing your own custom subclass of ``ContentsManager``. For example,
if you deploy the notebook in a context where you don't trust or don't have
access to the filesystem of the notebook server, it's possible to write your
own ContentsManager that stores notebooks and files in a database.


Required Methods
~~~~~~~~~~~~~~~~

A minimal complete implementation of a custom
:class:`~manager.ContentsManager` must implement the following
methods:

.. autosummary::
   ContentsManager.get
   ContentsManager.save
   ContentsManager.delete_file
   ContentsManager.rename_file
   ContentsManager.file_exists
   ContentsManager.dir_exists
   ContentsManager.is_hidden

You may be required to specify a Checkpoints object, as the default one,
``FileCheckpoints``, could be incompatible with your custom 
ContentsManager.


Chunked Saving
~~~~~~~~~~~~~~

The contents API allows for "chunked" saving of files, i.e.
saving/transmitting in partial pieces:

* This can only be used when the ``type`` of the model is ``file``.
* The model should be as otherwise expected for
  :meth:`~manager.ContentsManager.save`, with an added field ``chunk``.
* The value of ``chunk`` should be an integer starting at ``1``, and incrementing
  for each subsequent chunk, except for the final chunk, which should be
  indicated with a value of ``-1``.
* The model returned from using :meth:`~manager.ContentsManager.save` with
  ``chunk`` should be treated as unreliable for all chunks except the final one.
* Any interaction with a file being saved in a chunked manner is unreliable
  until the final chunk has been saved. This includes directory listings.


Customizing Checkpoints
-----------------------
.. currentmodule:: notebook.services.contents.checkpoints

Customized Checkpoint definitions allows behavior to be 
altered and extended.

The ``Checkpoints`` and ``GenericCheckpointsMixin`` classes
(from :mod:`notebook.services.contents.checkpoints`)
have reusable code and are intended to be used together, 
but require the following methods to be implemented.

.. autosummary::
   Checkpoints.rename_checkpoint
   Checkpoints.list_checkpoints
   Checkpoints.delete_checkpoint
   GenericCheckpointsMixin.create_file_checkpoint
   GenericCheckpointsMixin.create_notebook_checkpoint
   GenericCheckpointsMixin.get_file_checkpoint
   GenericCheckpointsMixin.get_notebook_checkpoint

No-op example
~~~~~~~~~~~~~

Here is an example of a no-op checkpoints object - note the mixin
comes first. The docstrings indicate what each method should do or 
return for a more complete implementation.

.. code-block:: python

    class NoOpCheckpoints(GenericCheckpointsMixin, Checkpoints):
        """requires the following methods:"""
        def create_file_checkpoint(self, content, format, path):
            """ -> checkpoint model"""
        def create_notebook_checkpoint(self, nb, path):
            """ -> checkpoint model"""
        def get_file_checkpoint(self, checkpoint_id, path):
            """ -> {'type': 'file', 'content': <str>, 'format': {'text', 'base64'}}"""
        def get_notebook_checkpoint(self, checkpoint_id, path):
            """ -> {'type': 'notebook', 'content': <output of nbformat.read>}"""
        def delete_checkpoint(self, checkpoint_id, path):
            """deletes a checkpoint for a file"""
        def list_checkpoints(self, path):
            """returns a list of checkpoint models for a given file, 
            default just does one per file
            """
            return []
        def rename_checkpoint(self, checkpoint_id, old_path, new_path):
            """renames checkpoint from old path to new path"""

See ``GenericFileCheckpoints`` in :mod:`notebook.services.contents.filecheckpoints`
for a more complete example.

Testing
-------
.. currentmodule:: notebook.services.contents.tests

:mod:`notebook.services.contents.tests` includes several test suites written
against the abstract Contents API.  This means that an excellent way to test a
new ContentsManager subclass is to subclass our tests to make them use your
ContentsManager.

.. note::

   PGContents_ is an example of a complete implementation of a custom
   ``ContentsManager``.  It stores notebooks and files in PostgreSQL_ and encodes
   directories as SQL relations.  PGContents also provides an example of how to
   re-use the notebook's tests.

.. _NBFormat: https://nbformat.readthedocs.io/en/latest/index.html
.. _PGContents: https://github.com/quantopian/pgcontents
.. _PostgreSQL: https://www.postgresql.org/
Custom bundler extensions
=========================

The notebook server supports the writing of *bundler extensions* that
transform, package, and download/deploy notebook files. As a developer, you
need only write a single Python function to implement a bundler. The notebook
server automatically generates a *File -> Download as* or *File -> Deploy as*
menu item in the notebook front-end to trigger your bundler.

Here are some examples of what you can implement using bundler extensions:

* Convert a notebook file to a HTML document and publish it as a post on a
  blog site
* Create a snapshot of the current notebook environment and bundle that
  definition plus notebook into a zip download
* Deploy a notebook as a standalone, interactive `dashboard <https://github.com/jupyter-incubator/dashboards_bundlers>`_

To implement a bundler extension, you must do all of the following:

* Declare bundler extension metadata in your Python package
* Write a `bundle` function that responds to bundle requests
* Instruct your users on how to enable/disable your bundler extension

The following sections describe these steps in detail.

Declaring bundler metadata
--------------------------

You must provide information about the bundler extension(s) your package
provides by implementing a `_jupyter_bundlerextensions_paths` function. This
function can reside anywhere in your package so long as it can be imported
when enabling the bundler extension. (See :ref:`enabling-bundlers`.)

.. code:: python

    # in mypackage.hello_bundler

    def _jupyter_bundlerextension_paths():
        """Example "hello world" bundler extension"""
        return [{
            'name': 'hello_bundler',                    # unique bundler name
            'label': 'Hello Bundler',                   # human-readable menu item label
            'module_name': 'mypackage.hello_bundler',   # module containing bundle()
            'group': 'deploy'                           # group under 'deploy' or 'download' menu
        }]

Note that the return value is a list. By returning multiple dictionaries in
the list, you allow users to enable/disable sets of bundlers all at once.

Writing the `bundle` function
-----------------------------

At runtime, a menu item with the given label appears either in the
*File ->  Deploy as* or *File -> Download as* menu depending on the `group`
value in your metadata. When a user clicks the menu item, a new browser tab
opens and notebook server invokes a `bundle` function in the `module_name`
specified in the metadata.

You must implement a `bundle` function that matches the signature of the
following example:

.. code:: python

    # in mypackage.hello_bundler

    def bundle(handler, model):
        """Transform, convert, bundle, etc. the notebook referenced by the given
        model.

        Then issue a Tornado web response using the `handler` to redirect
        the user's browser, download a file, show a HTML page, etc. This function
        must finish the handler response before returning either explicitly or by
        raising an exception.

        Parameters
        ----------
        handler : tornado.web.RequestHandler
            Handler that serviced the bundle request
        model : dict
            Notebook model from the configured ContentManager
        """
        handler.finish('I bundled {}!'.format(model['path']))

Your `bundle` function is free to do whatever it wants with the request and
respond in any manner. For example, it may read additional query parameters
from the request, issue a redirect to another site, run a local process (e.g.,
`nbconvert`), make a HTTP request to another service, etc.

The caller of the `bundle` function is `@tornado.gen.coroutine` decorated and
wraps its call with `torando.gen.maybe_future`. This behavior means you may
handle the web request synchronously, as in the example above, or
asynchronously using `@tornado.gen.coroutine` and `yield`, as in the example
below.

.. code:: python

    from tornado import gen

    @gen.coroutine
    def bundle(handler, model):
      # simulate a long running IO op (e.g., deploying to a remote host)
      yield gen.sleep(10)

      # now respond
      handler.finish('I spent 10 seconds bundling {}!'.format(model['path']))

You should prefer the second, asynchronous approach when your bundle operation
is long-running and would otherwise block the notebook server main loop if
handled synchronously.

For more details about the data flow from menu item click to bundle function
invocation, see :ref:`bundler-details`.

.. _enabling-bundlers:

Enabling/disabling bundler extensions
-------------------------------------

The notebook server includes a command line interface (CLI) for enabling and
disabling bundler extensions.

You should document the basic commands for enabling and disabling your
bundler. One possible command for enabling the `hello_bundler` example is the
following:

.. code:: bash

    jupyter bundlerextension enable --py mypackage.hello_bundler --sys-prefix

The above updates the notebook configuration file in the current
conda/virtualenv environment (`--sys-prefix`) with the metadata returned by
the `mypackage.hellow_bundler._jupyter_bundlerextension_paths` function.

The corresponding command to later disable the bundler extension is the
following:

.. code:: bash

    jupyter bundlerextension disable --py mypackage.hello_bundler --sys-prefix

For more help using the `bundlerextension` subcommand, run the following.

.. code:: bash

    jupyter bundlerextension --help

The output describes options for listing enabled bundlers, configuring
bundlers for single users, configuring bundlers system-wide, etc.

Example: IPython Notebook bundle (.zip)
---------------------------------------

The `hello_bundler` example in this documentation is simplistic in the name
of brevity. For more meaningful examples, see
`notebook/bundler/zip_bundler.py` and `notebook/bundler/tarball_bundler.py`.
You can enable them to try them like so:

.. code:: bash

    jupyter bundlerextension enable --py notebook.bundler.zip_bundler --sys-prefix
    jupyter bundlerextension enable --py notebook.bundler.tarball_bundler --sys-prefix

.. _bundler-details:

Bundler invocation details
--------------------------

Support for bundler extensions comes from Python modules in `notebook/bundler`
and JavaScript in `notebook/static/notebook/js/menubar.js`. The flow of data
between the various components proceeds roughly as follows:

1. User opens a notebook document
2. Notebook front-end JavaScript loads notebook configuration
3. Bundler front-end JS creates menu items for all bundler extensions in the
   config
4. User clicks a bundler menu item
5. JS click handler opens a new browser window/tab to
   `<notebook base_url>/bundle/<path/to/notebook>?bundler=<name>` (i.e., a
   HTTP GET request)
6. Bundle handler validates the notebook path and bundler `name`
7. Bundle handler delegates the request to the `bundle` function in the
   bundler's `module_name`
8. `bundle` function finishes the HTTP request
======================
Extending the Notebook
======================

Certain subsystems of the notebook server are designed to be extended or
overridden by users. These documents explain these systems, and show how to
override the notebook's defaults with your own custom behavior.

.. toctree::
    :maxdepth: 2

    contents
    savehooks
    handlers
    frontend_extensions
    keymaps
    bundler_extensions
Customize keymaps
=================

.. note::

    Declarative Custom Keymaps is a provisional feature with unstable API
    which is not guaranteed to be kept in future versions of the notebook,
    and can be removed or changed without warnings.

The notebook shortcuts that are defined by jupyter both in edit mode and 
command mode are configurable in the frontend configuration file
``~/.jupyter/nbconfig/notebook.json``. The modification of keyboard 
shortcuts suffers from several limitations, mainly that your Browser and OS 
might prevent certain shortcuts from working correctly. If this is the case,
there is unfortunately not much that can be done. The second issue can arise
with keyboards that have a layout different than US English. Again, even if 
we are aware of the issue, there is not much that can be done.

Shortcuts are also limited by the underlying library that handles code and 
text editing: CodeMirror. If some keyboard shortcuts are conflicting, the 
method described below might not work to create new keyboard shortcuts, 
especially in the ``edit`` mode of the notebook.


The 4 sections of interest in ``~/.jupyter/nbconfig/notebook.json`` are the
following:

  - ``keys.command.unbind``
  - ``keys.edit.unbind``
  - ``keys.command.bind``
  - ``keys.edit.bind``

The first two sections describe which default keyboard shortcuts not to 
register at notebook startup time. These are mostly useful if you need to 
``unbind`` a default keyboard shortcut before binding it to a new 
``command``.

The first two sections apply respectively to the ``command`` and ``edit``
mode of the notebook. They take a list of shortcuts to ``unbind``.

For example, to unbind the shortcut to split a cell at the position of the
cursor (``Ctrl-Shift-Minus``) use the following:

.. code:: javascript

    // file ~/.jupyter/nbconfig/notebook.json

    {
      "keys": {
        "edit": {
          "unbind": [
            "Ctrl-Shift-Minus"
          ]
        },
      },
    }




The last two sections describe which new keyboard shortcuts to register
at notebook startup time and which actions they trigger.

The last two sections apply respectively to the ``command`` and ``edit``
mode of the notebook. They take a dictionary with shortcuts as ``keys`` and
``commands`` name as value.

For example, to bind the shortcut ``G,G,G`` (Press G three time in a row) in
command mode to the command that restarts the kernel and runs all cells, use
the following:


.. code:: javascript

    // file ~/.jupyter/nbconfig/notebook.json

    {
      "keys": {
        "command": {
            "bind": {
                "G,G,G":"jupyter-notebook:restart-kernel-and-run-all-cells"
            }
        }
      },
    }




The name of the available ``commands`` can be find by hovering over the 
right end of a row in the command palette.
=================
Notebook Examples
=================

The pages in this section are all converted notebook files. You can also
`view these notebooks on nbviewer`__.

__ https://nbviewer.jupyter.org/github/jupyter/notebook/blob/master/
   docs/source/examples/Notebook/

.. toctree::
   :maxdepth: 2

   What is the Jupyter Notebook
   Notebook Basics
   Running Code
   Working With Markdown Cells
   Custom Keyboard Shortcuts
   JavaScript Notebook Extensions
   Importing Notebooks
   Connecting with the Qt Console
   Typesetting Equations
