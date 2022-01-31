# Changelog

## [6.15.0](https://www.github.com/snakemake/snakemake/compare/v6.14.0...v6.15.0) (2022-01-29)


### Features

* adding default_target directive for declaring default target rules that are not the first rule in the workflow. ([#1358](https://www.github.com/snakemake/snakemake/issues/1358)) ([638ec1a](https://www.github.com/snakemake/snakemake/commit/638ec1a983741cd7ba8faaf1a9dc76ae43d012e5))


### Bug Fixes

* Draft notebook filename with wildcards and params. ([#1352](https://www.github.com/snakemake/snakemake/issues/1352)) ([11d4dc8](https://www.github.com/snakemake/snakemake/commit/11d4dc88598ffb901450bd4e076b91f4e27d37b0))
* proper error message when defining cache eligibility for rules with multiple output files and no multiext declaration. ([#1357](https://www.github.com/snakemake/snakemake/issues/1357)) ([47b5096](https://www.github.com/snakemake/snakemake/commit/47b5096ebbdd3d94a9c99b443064b1b0de389c64))


### Documentation

* Command line arguments for configuration files ([#1343](https://www.github.com/snakemake/snakemake/issues/1343)) ([ad8aaa4](https://www.github.com/snakemake/snakemake/commit/ad8aaa4853a150211513baecc474956575d326eb))
* fix broken link in executor_tutorial/tutorial.rst ([#1360](https://www.github.com/snakemake/snakemake/issues/1360)) ([c9be764](https://www.github.com/snakemake/snakemake/commit/c9be76482d05577c4b1528b0e52ba15fc17a1dd5))

## [6.14.0](https://www.github.com/snakemake/snakemake/compare/v6.13.1...v6.14.0) (2022-01-26)


### Features

* Added timestamp to each log message ([#1304](https://www.github.com/snakemake/snakemake/issues/1304)) ([a5769f0](https://www.github.com/snakemake/snakemake/commit/a5769f0baeaa829b7813dee8c78902edbb42cf4b))
* implement support for removing GFAL remote files ([#1103](https://www.github.com/snakemake/snakemake/issues/1103)) ([25943e5](https://www.github.com/snakemake/snakemake/commit/25943e5630ff6d83afa5cba28edf473ce2ca87da))
* specify conda environments via their name ([#1340](https://www.github.com/snakemake/snakemake/issues/1340)) ([735ab23](https://www.github.com/snakemake/snakemake/commit/735ab2301d0905ea054ad6efa3150acb296d0e78))
* support for post deploy scripts ([#1325](https://www.github.com/snakemake/snakemake/issues/1325)) ([e5dac4f](https://www.github.com/snakemake/snakemake/commit/e5dac4ff297b7aeeb1e1a0bbdd03cb967cee3011))


### Documentation

* link to list of dependencies from installation ([#1336](https://www.github.com/snakemake/snakemake/issues/1336)) ([99d7bfe](https://www.github.com/snakemake/snakemake/commit/99d7bfef1285f131d0e60331511bc4833e7e414a))
* update URL to emacs snakemake-mode ([#1339](https://www.github.com/snakemake/snakemake/issues/1339)) ([dae7b8f](https://www.github.com/snakemake/snakemake/commit/dae7b8fb0e580a1878d36881cfb5ffc8adeaeb9f))

### [6.13.1](https://www.github.com/snakemake/snakemake/compare/v6.13.0...v6.13.1) (2022-01-11)


### Bug Fixes

* --conda-frontend value not passed on to cluster jobs ([#1317](https://www.github.com/snakemake/snakemake/issues/1317)) ([df46ddb](https://www.github.com/snakemake/snakemake/commit/df46ddb37022b291a4feca22fd0fbcf8773e7d03))
* atomic job error display ([#1326](https://www.github.com/snakemake/snakemake/issues/1326)) ([aa2c265](https://www.github.com/snakemake/snakemake/commit/aa2c2652608d3e95ad7fb568df09ef1ae09e1def))
* fix source cache handling for remote source files retrieved via github() or gitlab() tags. ([#1322](https://www.github.com/snakemake/snakemake/issues/1322)) ([6e2ecd2](https://www.github.com/snakemake/snakemake/commit/6e2ecd26e48eb64fa04c9c38dde591857e03c722))
* typos in code examples ([#1324](https://www.github.com/snakemake/snakemake/issues/1324)) ([60010e4](https://www.github.com/snakemake/snakemake/commit/60010e4ef07b7ba9b89aa5f48ee90ff3cec85b75))

## [6.13.0](https://www.github.com/snakemake/snakemake/compare/v6.12.3...v6.13.0) (2021-12-21)


### Features

* allow prefix definition in module statements ([#1310](https://www.github.com/snakemake/snakemake/issues/1310)) ([29e6540](https://www.github.com/snakemake/snakemake/commit/29e6540aac95b08b5e386a8478bd2013334e5954))

### [6.12.3](https://www.github.com/snakemake/snakemake/compare/v6.12.2...v6.12.3) (2021-12-09)


### Bug Fixes

* fixed display of any exceptions and errors from within a workflow definition ([23d40d9](https://www.github.com/snakemake/snakemake/commit/23d40d99614a88fd3c596d05e6915509ae43d4ce))

### [6.12.2](https://www.github.com/snakemake/snakemake/compare/v6.12.1...v6.12.2) (2021-12-07)


### Bug Fixes

* rule inheritance within modules (did previously lead to key errors) ([#1292](https://www.github.com/snakemake/snakemake/issues/1292)) ([603e0a8](https://www.github.com/snakemake/snakemake/commit/603e0a87d2c7af57a8f1d397605bc501c50934e0))


### Documentation

* Fix typo in rules.rst (—draft-notebook) ([#1290](https://www.github.com/snakemake/snakemake/issues/1290)) ([f5c42cf](https://www.github.com/snakemake/snakemake/commit/f5c42cfdc68f1516cec71b8ead8d78225ae915e5))

### [6.12.1](https://www.github.com/snakemake/snakemake/compare/v6.12.0...v6.12.1) (2021-11-29)


### Bug Fixes

* set default number of nodes to 1 in test cases ([#1288](https://www.github.com/snakemake/snakemake/issues/1288)) ([f6e12b4](https://www.github.com/snakemake/snakemake/commit/f6e12b4798485be3a1bb240b4af44d57dd5c84b2))

## [6.12.0](https://www.github.com/snakemake/snakemake/compare/v6.11.1...v6.12.0) (2021-11-29)


### Features

* add flag --draft-notebook for generating a skeleton notebook for manual editing (e.g. in VSCode). ([#1284](https://www.github.com/snakemake/snakemake/issues/1284)) ([d279322](https://www.github.com/snakemake/snakemake/commit/d2793223f914790c07b25363cb9b314ef166cb3e))


### Bug Fixes

* issue [#1257](https://www.github.com/snakemake/snakemake/issues/1257) (missing logfile failure when using shadow directory) ([#1258](https://www.github.com/snakemake/snakemake/issues/1258)) ([426d92f](https://www.github.com/snakemake/snakemake/commit/426d92fd9610b61b414b7f0152d777c463c939a2))
* keep empty output and input dirs of --draft-notebook job ([f1181bd](https://www.github.com/snakemake/snakemake/commit/f1181bd41ea8b20fafd3975c2733ca1d439381dc))
* SameFileError [#1153](https://www.github.com/snakemake/snakemake/issues/1153) ([#1220](https://www.github.com/snakemake/snakemake/issues/1220)) ([ede313d](https://www.github.com/snakemake/snakemake/commit/ede313dcd31ea5f136b3b8f743e2265331475342))
* snakemake API using only 1 job as default ([#1283](https://www.github.com/snakemake/snakemake/issues/1283)) ([e92ad48](https://www.github.com/snakemake/snakemake/commit/e92ad4867feb456ce8ef3dc57fd8528affa64ae9))


### Documentation

* short tutorial updates ([#1286](https://www.github.com/snakemake/snakemake/issues/1286)) ([b653a44](https://www.github.com/snakemake/snakemake/commit/b653a44d105e4b3799425a695d75a08239dc0d6b))

### [6.11.1](https://www.github.com/snakemake/snakemake/compare/v6.11.0...v6.11.1) (2021-11-26)


### Bug Fixes

* provide temporary IPYTHONDIR for notebook execution in order to avoid race conditions in https://github.com/ipython/ipython/blob/master/IPython/paths.py#L20 upon execution of multiple notebooks at the same time. ([#1280](https://www.github.com/snakemake/snakemake/issues/1280)) ([4d70da1](https://www.github.com/snakemake/snakemake/commit/4d70da11f810224ddce192ae1472a6380898865f))


### Documentation

* move psutil import into benchmark methods to avoid needing it as a dependency for doc building ([6ffe38d](https://www.github.com/snakemake/snakemake/commit/6ffe38d1740294a7170765ab875b363f4ae82cd4))
* require sphinx>=3 ([1773875](https://www.github.com/snakemake/snakemake/commit/1773875fc8f2fddb09362410afb7c49c4406bfa3))
* skip lazy property ([2883718](https://www.github.com/snakemake/snakemake/commit/28837183fa55a6764621580983b3d724f3881a6a))

## [6.11.0](https://www.github.com/snakemake/snakemake/compare/v6.10.0...v6.11.0) (2021-11-25)


### Features

* fail with an error if snakemake cannot write job metadata. ([#1273](https://www.github.com/snakemake/snakemake/issues/1273)) ([cd968cd](https://www.github.com/snakemake/snakemake/commit/cd968cd03437ad6db1d791f5d7ae5295b9754137))


### Bug Fixes

* Adds fixes for the first two MREs in [#823](https://www.github.com/snakemake/snakemake/issues/823) ([#1215](https://www.github.com/snakemake/snakemake/issues/1215)) ([cfd2f89](https://www.github.com/snakemake/snakemake/commit/cfd2f890a0af57628f7b9278d8d43f59b7006825))
* env file usage after changes to source file handling (inspired by [#1233](https://www.github.com/snakemake/snakemake/issues/1233) and [#1211](https://www.github.com/snakemake/snakemake/issues/1211)). ([#1236](https://www.github.com/snakemake/snakemake/issues/1236)) ([3ac8e85](https://www.github.com/snakemake/snakemake/commit/3ac8e858a7b908326922c8f68cae512b1250e906))
* fixed code change detection when using modules ([#1264](https://www.github.com/snakemake/snakemake/issues/1264)) ([b571e09](https://www.github.com/snakemake/snakemake/commit/b571e09ce452f6a1a95395e1c3c8b9e3f83867ad))
* handle config file extension/overwriting more explicitly ([#1251](https://www.github.com/snakemake/snakemake/issues/1251)) ([d0a7bf2](https://www.github.com/snakemake/snakemake/commit/d0a7bf243c5df204136fa1f14706aab793793c68))
* Issue [#1253](https://www.github.com/snakemake/snakemake/issues/1253) (problems editing Jupyter Notebooks) ([#1255](https://www.github.com/snakemake/snakemake/issues/1255)) ([3398ddf](https://www.github.com/snakemake/snakemake/commit/3398ddffd1f68182af768ef4ea519e9a9ad4efaf))
* more informative nothing to be done message ([#1234](https://www.github.com/snakemake/snakemake/issues/1234)) ([368d265](https://www.github.com/snakemake/snakemake/commit/368d265ff3da984bd3a53b319dcb882d6916975b))
* only consider context of shell command for technical switches if called from snakemake rules. ([#1213](https://www.github.com/snakemake/snakemake/issues/1213)) ([4816a58](https://www.github.com/snakemake/snakemake/commit/4816a58653e466ca94b1482a1d947a856f5381b3))
* R encoding of pathlib.Path objects ([#1201](https://www.github.com/snakemake/snakemake/issues/1201)) ([bd516e9](https://www.github.com/snakemake/snakemake/commit/bd516e958af22e57c18cacf0cb22552c2a237bd8))
* Use 'snakemake.utils.update_config' instead of 'dict.update' ([#1126](https://www.github.com/snakemake/snakemake/issues/1126)) ([2658027](https://www.github.com/snakemake/snakemake/commit/2658027458dde4c10b3d6e1af7671564d175f9cb))

## [6.10.0](https://www.github.com/snakemake/snakemake/compare/v6.9.1...v6.10.0) (2021-10-21)


### Features

* Add more informative errors when evaluation of `--default-resources` fails ([#1192](https://www.github.com/snakemake/snakemake/issues/1192)) ([b3c4e68](https://www.github.com/snakemake/snakemake/commit/b3c4e687c87c75075393cef842b129dcec70e7f6))


### Bug Fixes

* add quotes to each item of the wait_for_files list ([#1160](https://www.github.com/snakemake/snakemake/issues/1160)) ([72856ed](https://www.github.com/snakemake/snakemake/commit/72856edd12fbe29d723731c6f596f05cd2b59c0e))
* caching process ([#1225](https://www.github.com/snakemake/snakemake/issues/1225)) ([0825a29](https://www.github.com/snakemake/snakemake/commit/0825a29e46c08b200efe6bd0c66acf1e6828eed8))
* enable usage of job grouping in GLS ([#1054](https://www.github.com/snakemake/snakemake/issues/1054)) ([d243c22](https://www.github.com/snakemake/snakemake/commit/d243c22ff494b63bd5e07b7c5bf1f6ff32539cde))
* Only --bind Snakemake when we're working with a Python script ([#1206](https://www.github.com/snakemake/snakemake/issues/1206)) ([1d79f62](https://www.github.com/snakemake/snakemake/commit/1d79f625b7262d66def71c779f2a2c091bc418d8))
* run dependencies with non-existent ancient files before the consuming job ([#1202](https://www.github.com/snakemake/snakemake/issues/1202)) ([84d1f64](https://www.github.com/snakemake/snakemake/commit/84d1f6451b12352eba5a8bfefcfcce8b2d98c5aa)), closes [#946](https://www.github.com/snakemake/snakemake/issues/946)
* status cmd repeats until killed by 11 *different* signals ([#1207](https://www.github.com/snakemake/snakemake/issues/1207)) ([8b28b57](https://www.github.com/snakemake/snakemake/commit/8b28b5740c34149c9b5df56dbbfa034219eb1574))
* typo in sourcecache use ([#1229](https://www.github.com/snakemake/snakemake/issues/1229)) ([8b54bc5](https://www.github.com/snakemake/snakemake/commit/8b54bc5db9d8e5c0bcb8f2c2ff141dc075e3e659))
* wms monitor arg parsing now accepts any kind of value ([#1181](https://www.github.com/snakemake/snakemake/issues/1181)) ([313de93](https://www.github.com/snakemake/snakemake/commit/313de932e2e2a4f2c530df18c1abb15d37eb3217))


### Documentation

* Clarification of --cluster-stats docs  &  elaborating on the situation where job ids are not passed to the status script ([#1221](https://www.github.com/snakemake/snakemake/issues/1221)) ([ed0e4a2](https://www.github.com/snakemake/snakemake/commit/ed0e4a27a2167a69a4fe1bcdf237dd27bb3732ca))
* Combine CHANGELOG.rst with CHANGELOG.md ([#1228](https://www.github.com/snakemake/snakemake/issues/1228)) ([19f5a43](https://www.github.com/snakemake/snakemake/commit/19f5a43261bd6ba548d6f01080640f0d4119871e))
* Mention required openssl dep for rust-script ([#1216](https://www.github.com/snakemake/snakemake/issues/1216)) ([fc8c5f6](https://www.github.com/snakemake/snakemake/commit/fc8c5f62c397a0239ef213ab45a26a1def50f9eb))
* Unpin docutils version ([#1230](https://www.github.com/snakemake/snakemake/issues/1230)) ([15a82bf](https://www.github.com/snakemake/snakemake/commit/15a82bfe402b3577bf19e6d2eca3b2fb86109628))

### [6.9.1](https://www.github.com/snakemake/snakemake/compare/v6.9.0...v6.9.1) (2021-09-30)


### Bug Fixes

* fix function call when creating report and hashes for between workflow caching ([#1198](https://www.github.com/snakemake/snakemake/issues/1198)) ([a4f6836](https://www.github.com/snakemake/snakemake/commit/a4f68365125c357f30510d0e61036f98b9d3aa69))

## [6.9.0](https://www.github.com/snakemake/snakemake/compare/v6.8.2...v6.9.0) (2021-09-29)


### Features

* autoconvert Path objects to str when passing to R or Julia scripts ([80ec513](https://www.github.com/snakemake/snakemake/commit/80ec51322f8134180c52c20b0a9dc6980df6c1bc))


### Bug Fixes

* fix source retrieval during between workflow caching and report generation ([2394ca4](https://www.github.com/snakemake/snakemake/commit/2394ca4a23a6b2792397bc9efc09945f01d1963b))

### [6.8.2](https://www.github.com/snakemake/snakemake/compare/v6.8.1...v6.8.2) (2021-09-29)


### Bug Fixes

* fix path returned by get_source() ([ee05315](https://www.github.com/snakemake/snakemake/commit/ee053153d2f44156171c127307cb110791b7624a))

### [6.8.1](https://www.github.com/snakemake/snakemake/compare/v6.8.0...v6.8.1) (2021-09-24)


### Bug Fixes

* async_run to allow nested event loops. ([#1170](https://www.github.com/snakemake/snakemake/issues/1170)) ([5dc6bbd](https://www.github.com/snakemake/snakemake/commit/5dc6bbd440ac46e81a926b6749969b98b7e33a9f))
* merging of pipe groups when multiple rules are chained together via pipes ([#1173](https://www.github.com/snakemake/snakemake/issues/1173)) ([de91d2c](https://www.github.com/snakemake/snakemake/commit/de91d2ccf53bd844b4dbf4f64dd087f4ee935be5))
* potential memory corruption caused by Google storage objects accessed from different threads ([#1174](https://www.github.com/snakemake/snakemake/issues/1174)) ([41a5071](https://www.github.com/snakemake/snakemake/commit/41a5071b750dca5d7fceec324d81d9a93c86bdb6))


### Performance Improvements

* more extensive caching of source files, including wrappers. ([#1182](https://www.github.com/snakemake/snakemake/issues/1182)) ([bdb75f8](https://www.github.com/snakemake/snakemake/commit/bdb75f828a3ae27ba97ea6cd5e71a34ac7b27eea))


### Documentation

* move note ([75a544b](https://www.github.com/snakemake/snakemake/commit/75a544ba528b30b43b861abc0ad464db4d6ae16f))
* polish ([47a7b62](https://www.github.com/snakemake/snakemake/commit/47a7b628686258a28dd870f20bf1f121b3a881c3))
* tutorial formatting ([594f5fb](https://www.github.com/snakemake/snakemake/commit/594f5fbb342e0722318641dea07d7da4c5eb8116))

## [6.8.0](https://www.github.com/snakemake/snakemake/compare/v6.7.0...v6.8.0) (2021-09-06)


### Features

* Add `shadow: "copy-minimal"` directive ([#1155](https://www.github.com/snakemake/snakemake/issues/1155)) ([1803f0b](https://www.github.com/snakemake/snakemake/commit/1803f0b4090d812df0c164653b26502fd130d326))
* support XRootD as a default remote provider ([#1017](https://www.github.com/snakemake/snakemake/issues/1017)) ([fe03157](https://www.github.com/snakemake/snakemake/commit/fe03157c31210984fce53c35d5fb87b20d278fe7))


### Bug Fixes

* AmbiguousRuleException bug caused by weak ordering of rules ([#1124](https://www.github.com/snakemake/snakemake/issues/1124)) ([7f54c39](https://www.github.com/snakemake/snakemake/commit/7f54c391f2821655ed168bcdafad6d07b96fcec7))
* Bugfix tes add files ([#1133](https://www.github.com/snakemake/snakemake/issues/1133)) ([8892bf2](https://www.github.com/snakemake/snakemake/commit/8892bf25d9d981a4032d5a1b525960ba3bdd1aec))
* Disable Persistence cache for snakemake jobs ([#1159](https://www.github.com/snakemake/snakemake/issues/1159)) ([7110f9d](https://www.github.com/snakemake/snakemake/commit/7110f9d2e7ee3f350bd1da3c5b4aab98c06725a1))
* efficient job status checking when using DRMAA API (this should yield much better parallelization and performance when using --drmaa) ([#1156](https://www.github.com/snakemake/snakemake/issues/1156)) ([ac004cb](https://www.github.com/snakemake/snakemake/commit/ac004cb19cebd4efb5e38f6039861a2810c702ff))
* improved error handling for cluster status scripts and smarter job selector choice in case of cluster submission (use greedy for single jobs). ([#1142](https://www.github.com/snakemake/snakemake/issues/1142)) ([48d2dd9](https://www.github.com/snakemake/snakemake/commit/48d2dd99a745fd54b74b1435cbb7e41e0ee1b4ac))
* Initialize assignments dictionary when setting rule-based resources ([#1154](https://www.github.com/snakemake/snakemake/issues/1154)) ([68c13fd](https://www.github.com/snakemake/snakemake/commit/68c13fd6fb2ad458e79bafe146499b601bf4bd0e))
* key error when handling FileNotFoundError in input functions. ([#1138](https://www.github.com/snakemake/snakemake/issues/1138)) ([d25f04d](https://www.github.com/snakemake/snakemake/commit/d25f04db820c9651835b7323baef5931d4f8dc0a))
* linting of remote snakefiles ([#1131](https://www.github.com/snakemake/snakemake/issues/1131)) ([2104e10](https://www.github.com/snakemake/snakemake/commit/2104e10d1d2c5e0f368e9c0fe95cc50f9d4847f1))


### Performance Improvements

* improve job selection performance in case of potential ambiguity that is resolved by comprehensive ruleorder statements. ([#1147](https://www.github.com/snakemake/snakemake/issues/1147)) ([921f4f7](https://www.github.com/snakemake/snakemake/commit/921f4f715e3814fc2e22a4f6527ff62e066cc5da))

## [6.7.0](https://www.github.com/snakemake/snakemake/compare/v6.6.1...v6.7.0) (2021-08-12)


### Features

* Add support for rust scripts (enabling directly integrated ad-hoc robust high performance scripting) ([#1053](https://www.github.com/snakemake/snakemake/issues/1053)) ([f0e8fa2](https://www.github.com/snakemake/snakemake/commit/f0e8fa285437a02ca7edcf87334bf00cb347064a))


### Bug Fixes

* Ga4gh tes bugfixes ([#1127](https://www.github.com/snakemake/snakemake/issues/1127)) ([af21d6c](https://www.github.com/snakemake/snakemake/commit/af21d6c2b125c22ef3dbc36a0a6a67a1874549c7))
* improved display of percentage of done jobs ([1fee8c0](https://www.github.com/snakemake/snakemake/commit/1fee8c06d6ed229d7e3757de3c693e755d01d1bb))
* improved error message in case of target rule misspecification ([83b1f5b](https://www.github.com/snakemake/snakemake/commit/83b1f5bbde437e13641be2160f4855f54043c046))


### Documentation

* fix contributing executors link ([#1112](https://www.github.com/snakemake/snakemake/issues/1112)) ([4bb58d1](https://www.github.com/snakemake/snakemake/commit/4bb58d12a44f77f79d47c5443f927cb6061677f5))
* Fix typo in file path in remote files documentation ([#1110](https://www.github.com/snakemake/snakemake/issues/1110)) ([9ce294f](https://www.github.com/snakemake/snakemake/commit/9ce294f6d5bdf72055a824ab610488a7f832a4d3))

### [6.6.1](https://www.github.com/snakemake/snakemake/compare/v6.6.0...v6.6.1) (2021-07-19)


### Bug Fixes

* avoid superfluous calls of conda info that have slowed down Snakemake since 6.4.1. ([#1099](https://www.github.com/snakemake/snakemake/issues/1099)) ([e990927](https://www.github.com/snakemake/snakemake/commit/e9909273c22a316dbd7301a243498e3c2a372642))

## [6.6.0](https://www.github.com/snakemake/snakemake/compare/v6.5.5...v6.6.0) (2021-07-16)


### Features

* Allow to mark all output files as temp with --all-temp ([#1097](https://www.github.com/snakemake/snakemake/issues/1097)) ([0ac3b38](https://www.github.com/snakemake/snakemake/commit/0ac3b3806c065d0ec3a551a5992faf30ddcf0576))

### [6.5.5](https://www.github.com/snakemake/snakemake/compare/v6.5.4...v6.5.5) (2021-07-16)


### Bug Fixes

* dummy release ([e4dca50](https://www.github.com/snakemake/snakemake/commit/e4dca508f6cbd3427d8580ef61f274f909ec8bab))

### [6.5.4](https://www.github.com/snakemake/snakemake/compare/v6.5.3...v6.5.4) (2021-07-16)


### Fixes

* Fixed --touch in combination with temp files (issue #1028) (@johanneskoester, @iromeo).

### Documentation

* Fix syntax error in docs/conf.py and update sphinx.ext.napoleon import ([#1084](https://www.github.com/snakemake/snakemake/issues/1084)) ([3e3fac2](https://www.github.com/snakemake/snakemake/commit/3e3fac2dbd5a8abad67e252f6181ad14bcfcb711))
* Improved pepfile (pepschema) documentation (@stolarczyk).
### \[6.5.3\] - 2021-07-06

-   Fixed a bug occuring when using --resources in the command line
    interface (@johanneskoester).
-   Minor improvements in the docs (@johanneskoester).

### \[6.5.2\] - 2021-07-02

-   Create directory pointed to by tmpdir resource if it does not yet
    exist (@johanneskoester).
-   Use a single core again in dryrun if --cores is not specified
    (@johanneskoester).
-   Bugfix for FTP remote provider (@jmeppley).
-   Improved documentation (@corneliusroemer).

### \[6.5.1\] - 2021-06-24

-   Extended best practices document (@johanneskoester)
-   Restore `-j all` behavior for local execution as a (deprecated) way
    of running Snakemake on all cores. Recommended now: `--cores all`
    (@johanneskoester).
-   Improved handling and better error messages for checkpoints
    (@johanneskoester).

## \[6.5.0\] - 2021-06-22

-   Allow to set the default profile via the environment variable
    $SNAKEMAKE_PROFILE.
-   There is a new default resource tmpdir (by default reflects the
    system setting), which is automatically used for temporary files by
    shell commands and scripts which properly consider the usual
    environment variables like $TMP, $TEMP, $TMPDIR (@johanneskoester).
-   The CLI flags --jobs and --cores are now separated, with --cores
    being responsible for local cores and global cores in the cluster
    case, and --jobs being responsible for number of jobs. Still -j and
    --jobs works as a fallback for local execution (@johanneskoester).
-   Added the ability to overwrite resources via --set-resources
    (@johanneskoester).
-   Various fixes for Windows execution (@melund).
-   Fixed a bug with fractional resources (@johanneskoester).
-   Fixed timeouts and other issues in google life science backend
    (@johanneskoester).
-   Fixed a bug with missing conda frontend definitions in subworkflows
    (@johanneskoester).
-   Skip envvar checking during linting (@johanneskoester).
-   Fixed a bug causing container images in modules to be ignored
    (@johanneskoester).

### \[6.4.1\] - 2021-05-27

-   Fixed bug in `workflow.source_path()` that occurred with modules
    included from remote locations (@johanneskoester).
-   Inform cluster jobs about conda/mamba/activate path such that they
    don't need to determine this themselves (@johanneskoester).

## \[6.4.0\] - 2021-05-20

-   Improvements in the docs (resource usage, best practices, remote
    files) (@johanneskoester, @admorris).
-   functions given to params can now safely open input files generated
    by previous rules. If they are not present, TBD will be displayed
    and function will be reevaluated immediately before the job is
    executed (i.e. when files are present) (@ASLeonard).
-   Connection pool for SFTP and FTP remote files, increasing download
    performance (@jmeppley).
-   Require correct minimum version of smart_open
    (@Redmar-van-den-Berg).
-   Added workflow.source_path(path), allowing to get the correct path
    relative to the current Snakefile, even when Snakefile is included
    via URL (@johanneskoester).
-   Fixed bugs in module system (@johanneskoester, @dlaehnemann).
-   Fixed issue with checkpoints and ruleorder where phantom
    dependencies are not properly removed from the DAG (@jmeppley,
    @johanneskoester).
-   Disable tibanna behavior that opens a browser window for each job
    (@nigiord).
-   Allow `Paramspace(..., filename_params="*")`, meaning that all
    columns of the paramspace will be encoded into the filename (@kpj).
-   Avoid PATH modification in cluster jobs (@johanneskoester).
-   For large sets of input files, pass files to wait for (FS latency)
    as a file instead of command line args (@kpj, @epruesse).

## \[6.3.0\] - 2021-04-29

-   Changed behavior of `workflow.snakefile` to always point to the
    current file instead of the main Snakefile (also in case of includes
    and modules) (@johanneskoester).
-   Fixed a typo in an error message (@nikostr).

## \[6.2.0\] - 2021-04-22

-   Support for integration of foreign workflow management systems by
    introducing a `handover` directive that passes on all resources to a
    particular rule (which can then invoke another workflow management
    system). See the docs ("Integrating foreign workflow management
    systems") (@johanneskoester).
-   Behavior improvement for temp handling of checkpoint rules
    (@epruesse).
-   Several improvements in the docs (@johanneskoester).

### \[6.2.1\] - 2021-04-20

-   Fixed a minor bug in the linter.

## \[6.2.0\] - 2021-04-20

-   Fixed several glitches in paramspace implementation (handling of
    bools, returning scalar values) (@kpj).
-   Fixed bugs in module implementation (@dlaehnemann,
    @johanneskoester).
-   Fall back to greedy scheduling solver if ILP solver needs more than
    10 sec (@johanneskoester).

### \[6.1.1\] - 2021-04-07

-   Fixed several small bugs of the new module system (@johanneskoester,
    @dlaehnemann).
-   Fixed archive based conda deployment (@johanneskoester).
-   Better handling of download and target attributed in the interactive
    report (@johanneskoester).

## \[6.1.0\] - 2021-04-01

-   Snakemake now uses **mamba** as the default conda frontend (which
    can be overwritten by specifying to use conda via the
    --conda-frontend flag) (@johanneskoester).
-   Profiles using --cluster option can now handle relative submit
    script paths in combination with arguments (@kdm9).
-   New AutoRemoteProvider, which infers the type of remote file
    protocol from the given URL (@kpj).
-   When using global container directive, container usage can be
    deactivated on a per rule base (@bilke).
-   Bugfixes for checkpoint handling (@johanneskoester).
-   Bugfixes for the module system (@johanneskoester, @dlaehnemann).
-   Various improvements for the tutorial.

### \[6.0.5\] - 2021-03-11

-   Fix bug (introduced with 6.0) when handling of HTML directories in
    report (@johanneskoester).

### \[6.0.4\] - 2021-03-11

-   Various textual improvements in the tutorial (@dlaehnemann).

### \[6.0.3\] - 2021-03-08

-   No longer use a shortened hash for naming conda environments in
    .snakemake/conda (@johanneskoester).
-   Various little updates to the docs (@johanneskoester).

### \[6.0.2\] - 2021-03-03

-   Fix race condition in conda checking code (@johanneskoester).

### \[6.0.1\] - 2021-03-03

-   Restored Python 3.5 compatibility by removing f-strings (@mbhall88)
-   Fix rendering issue in the docs.
-   Add gitpod dev environment and gitpod environment for the tutorial.

## \[6.0.0\] - 2021-02-26

-   Introduced a new module system, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules>
    (@johanneskoester).
-   Introduced a rule inheritance mechanism, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#rule-inheritance>
    (@johanneskoester).
-   Automatically containerize a conda-based pipeline with
    `--containerize`, see
    <https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html#containerization-of-conda-based-workflows>
    (@johanneskoester).
-   Use temporary files for long shell commands (@epruesse).
-   Various fixes in the documentation (@ctb, @SilasK, @EthanHolleman).
-   Fixed a bug in job grouping that led to non-deterministic behavior
    (@johanneskoester).

### \[5.32.2\] - 2021-02-11

### Changed

-   Fixed infinite loading of results in Snakemake
    reports (@FelixMoelder)

### \[5.32.1\] - 2021-02-08

### Changed

-   Improved warning on wildcard constraints (@jheuel)
-   Improved logging from the new scheduler
    implementation (@johanneskoester)
-   Restored Python 3.5 compatibility by removing f-strings (@mbhall88)
-   Snakemake now automatically adds a global wildcard constraint for
    {scatteritem}, when scatter/gather support is used.
-   The zip variant of Snakemake reports is now compressed
    (@FelixMoelder).
-   Improved docs (@ctb).
-   Make output file removal in cluster mode more robust (@sebschmi).

## \[5.32.0\] - 2021-01-15

#### Changed

-   Handle accidental use of GLS backend with singularity (@vsoch).
-   Improved and extended WMS-monitor implementation (@vsoch).
-   Display index and total count in `{scatteritem}` when using the
    scatter-gather helper (@johanneskoester).
-   Fixed problems with jobid display when handling checkpoint updates
    (@johanneskoester, @jmeppley).
-   Fixed bug when checking for directory containment of output files
    (@jmeppley).
-   Implement --no-subworkflows treatment in combination with --cluster
    (@goi42).

### \[5.31.1\] - 2020-12-21

#### Changed

-   added wget again to the container image

## \[5.31.0\] - 2020-12-21

#### Added

-   The `Paramspace` helper for automatically exploring parameter spaces
    given as Pandas dataframes.
-   A new directive `name:` for setting rule names from variables.

#### Changed

-   Various small bug fixes for scheduling and checkpoint handling.
-   Automatically block R_LIBS, PYTHONPATH, PERL5LIB, and PERLLIB when
    using conda with --use-conda. This behavior can be deactivated with
    --conda-not-block-envvars.
-   Update container image to latest singularity.

### \[5.30.2\] - 2020-12-16

#### Changed

-   Fix permission issues with jobscripts on some systems (@Phhere).
-   Added notes on WSL to the tutorial (@RomainFeron).
-   Scheduler fixes (@johanneskoester).
-   Fixed a bug in checkpoint handling that led to hanging workflow
    execution (@jmeppley).
-   Pass cluster nodes to subworkflows (@votti).
-   Fix start time recording in metadata (@lparsons).
-   Fix time retrieval in reports (@johanneskoester).
-   Fix error when returning a Path from an input function (@sappjw).
-   Extending monitoring docs with some notes about future api changes
    (@vsoch).

## \[5.30.0\] - 2020-11-23

#### Added

-   Benchmarks now also report CPU time (@natir).

#### Changed

-   Fixed a reauthentication bug in Kubernetes support (@haizi-zh).

## \[5.29.0\] - 2020-11-19

#### Changed

-   Fixed several bugs in reports and scheduler.
-   Remove automatic (but buggy) encoding of csv/tsv files into HTML
    tables in the report (we will soon have a better alternative).
-   Fixed bug in kubernetes executor occurring with large source files.

## \[5.28.0\] - 2020-11-12

#### Added

-   Execution backend for GA4GH TES (task execution scheduler) an
    abstraction layer for various cluster and cloud queuing systems
    (@svedziok, @uniqueg).
-   script, notebook, wrapper and cwl directives now permit to use
    wildcards and params for composing paths (@johanneskoester).

#### Changed

-   Restored compatibility with Python 3.5 and 3.6 (@cclienti).
-   Various usability bug fixes (@goi43, @johanneskoester, @dcroote).
-   Better and more secure parsing of values when using --config
    (@bingxiao).

### \[5.27.4\] - 2020-11-03

#### Changed

-   Further speed improvements for DAG computation.
-   Fixed metadata migration errors occuring with long output file
    paths.
-   Add WorkflowHub specifications to the docs.
-   Fix group assignments.

### \[5.27.3\] - 2020-10-30

#### Changed

-   Added missing files to source distribution.

### \[5.27.2\] - 2020-10-30

#### Changed

-   DAG computation runtime has been improved by orders of magnitude, it
    is linear in the number of jobs now (@mhulsmann, @johanneskoester).
-   Stat calls have been dramatically reduced and are now performed in
    parallel (@johanneskoester).
-   Scheduler fixes (@FelixMoelder).
-   Directory support and other fixes for Google Life Sciences backend
    (@vsoch, @millerdz).
-   Support for panoptes monitor server (@fgypas).
-   Extended pathlib support (@mbhall88).
-   Vim plugin improvements (@troycomi).
-   Prevent jobs being rerun when input files are marked as ancient and
    another job in the DAG creates them.
-   Fixed --list-code-changes for included rules (@jbloom).

#### Added

-   Syntax highlighting for nano (@baileythegreen).

### \[5.26.1\] - 2020-10-01

#### Changed

-   Use coin ILP solver for scheduling by default (GLPK has bugs that
    can cause it to fail in certain situations).
-   If coin is not available, fall back to greedy scheduler.

## \[5.26.0\] - 2020-09-30

#### Added

-   Flag --max-inventory-time for setting maximum time spend on creating
    file inventory.
-   Flag --scheduler-ilp-solver for defining which solver to use for the
    ILP scheduler.

#### Changed

-   Fixed various bugs with the new scheduler (@FelixMoelder).
-   Fixed bug causing certain parameters not to be passed to the cluster
    (--set-scatter, --scheduler, --set-threads).
-   Updated docs and fixed of google backend (@vsoch).
-   Display jupyter notebook code in reports.
-   Improved scheduler behavior in order to directly remove temporary
    files if possible.

## \[5.25.0\] - 2020-09-18

#### Added

-   Simplified and more configurable support for scatter-gather
    processes (see docs).
-   Fully configurable DAG partitioning by grouping jobs at the command
    line. This should provide a vast additional improvement to
    scalability in cluster and cloud settings.

#### Changed

-   Depend on latest pulp, thereby enable Python >=3.8 compatibility
    again.
-   Fixes for snakefile handling in google life sciences backend
    (@vsoch).

### \[5.24.2\] - 2020-09-15

#### Changed

-   Fixed a bug in the linter that caused a false warning when using
    resources in shell commands.

### \[5.24.1\] - 2020-09-13

#### Changed

-   Depend on pulp \< 2.0, which includes the default coin cbc solver
    for all platforms.

## \[5.24.0\] - 2020-09-09

#### Added

-   Preemtion support for google cloud backend (@vsoch).

#### Changed

-   Fixed compatibility issues in new scheduler code (@dtrodrigues and
    @johanneskoester).
-   Improved error messages (@Sam-Tygier, @terrycojones)
-   Various small bug fixes.
-   Improved profile documentation (@johanneskoester).

## \[5.23.0\] - 2020-08-24

#### Added

-   Support for workflow configuration via portable encapsulated
    projects (PEPs, <https://pep.databio.org>).
-   A new ILP based default scheduler now ensures that temporary files
    are deleted as fast as possible (@FelixMoelder, @johanneskoester).

#### Changed

-   Fixed bug in modification date comparison for files in google
    storage (@vsoch).
-   Various small documentation improvements (@dcroote, @erjel,
    @dlaehnemann, @goi42).

### \[5.22.1\] - 2020-08-14

#### Changed

-   Fixed a missing dependency for google storage in cloud execution.

## \[5.22.0\] - 2020-08-13

#### Added

-   Added short option `-T` for CLI parameter `--restart-times`
    (@mbhall88).

#### Changed

-   Various small fixes for google storage and life sciences backends
    (@vsoch).

## \[5.21.0\] - 2020-08-11

#### Changed

-   Added default-remote-provider support for Azure storage
    (@andreas-wilm).
-   Various small bug fixes and documentation improvements.

### \[5.20.1\] - 2020-07-08

#### Changed

-   Fixed a bug that caused singularity args to be not passed on
    correctly when using script or conda.

## \[5.20.0\] - 2020-07-08

#### Changed

-   Exceptions in input functions are now handled in a smarter way, by
    choosing alternative paths in the DAG if available.
-   Debugging dag creation (--debug-dag) now gives more hints if
    alternative DAG paths are chosen.
-   Fixes for XRootD remote file implementation.
-   Improved CLI documentation.
-   Improved docs.
-   Various minor bug fixes.
-   Restored Python 3.5 compatibility.
-   Speed improvements for workdir cleanup.
-   Allow Path objects to be passed to expand.

### \[5.19.3\] - 2020-06-16

#### Changed

-   Performance improvements for DAG generation (up to 7x in the google
    cloud, anything from a little to massive in a cluster, depending on
    the overall filesystem performance).
-   Made harcoded bucket in google cloud executor configurable.
-   Improved speed of --unlock command.

### \[5.19.2\] - 2020-06-04

#### Changed

-   Fixed a bug in script and wrapper directives. Tried to decode a str.

### \[5.19.1\] - 2020-06-03

#### Changed

-   Fixed an issue with the parameter linting code, that could cause an
    index out of bounds exception.

## \[5.19.0\] - 2020-06-02

#### Added

-   The multiext function now allows arbitrary file extensions (no
    longer required to start with a "." (thanks to @jafors)
-   The include directive can now also take a Pathlib Path object
    (thanks to @mbhall88).

#### Changed

-   Jupyter notebook integration no longer automatically starts a
    browser.
-   Empty directories are cleaned up after workflow execution.
-   Fixed directory handling: no longer fail if the same job writes both
    a dir and a contained file.
-   Linter now recommends using spaces only for indentation.
-   Persistence dir "aux" has been renamed to "auxilliary" in order to
    make windows happy.
-   Linter now distinguishes awk syntax from regular variable usage.
-   Various bug fixes for Windows (thanks to @melund).

## \[5.18.0\] - 2020-05-21

#### Added

-   Native Google Cloud support via the (despite the name generic)
    lifesciences API.
-   Ability to optionally exchange the conda frontend to mamba (faster
    and sometimes more correct) instead of conda.

#### Changed

- Improved notebook integration experience, with various removed bugs and
  pitfalls.
- Auto-retry google storage API calls on transient or checksum errors.

## \[5.17.0\] - 2020-05-07

#### Added

- --envvars flag for passing secrets to cloud executors

#### Changed

- Wider thumbnail dialogs in report.
- Updated installation instructions.
- Various small kubernetes bug fixes.
- Bug fix for iRods remote files.

## \[5.16.0\] - 2020-04-29

#### Added

- Interactive jupyter notebook editing. Notebooks defined by rules can
  be interactively drafted and updated using snakemake --edit-notebook
  (see docs).

#### Changed

- Fixed group resource usage to occupy one cluster/cloud node.
- Minor bug fixes.

## \[5.15.0\] - 2020-04-21

#### Changed

-   The resource directive can now take strings, e.g. for defining a GPU
    model (see docs). This will e.g. be used for upcoming updates to
    cloud executors.
-   More extensive conda cleanup with --conda-cleanup-packages, meant
    for CI usage.
-   Further polish for reports.

## \[5.14.0\] - 2020-04-08

#### Changed

-   Redesigned HTML reports, with improved interface and performance.
-   For big data, HTML reports can now be stored as ZIP, where files are
    not anymore embedded but rather are stored in an auxilliary folder,
    such that they don't have to be in memory during report rendering.
-   Added subcategories to report (see docs).
-   Fixed a bug linter, leading to only one rule or snakefile to be
    linted.
-   Breaking change in CLI: added flags --conda-cleanup-envs and
    --conda-cleanup-pkgs, removed flag --cleanup-conda.
-   Fixed scheduling of pipe jobs, they are now always scheduled, fixing
    a hangup.
-   Corrected quoting of shell command for cluster submission.

## \[5.13.0\] - 2020-03-27

#### Added

- Allow to flag directories for inclusion in the report.

#### Changed

- Fixed hash computation for --cache in case of positional params
  arguments.
- Automatically restrict thread usage of linear algebra libraries to whatever
  is specified in the rule/job.

### \[5.12.3\] - 2020-03-24

#### Changed

-   Various minor bug fixes.

### \[5.12.2\] - 2020-03-24

#### Changed

-   Further improved linter output.

### \[5.12.1\] - 2020-03-24

#### Changed

-   Linter fixes

## \[5.12.0\] - 2020-03-24

#### Changed

-   Fixed the ability to supply functions for the thread directive.
-   Improved error messages for caching.

#### Added

-   A new "cache: true" directive that allows to annotate between
    workflow caching eligibility for rules in the workflow.

### \[5.11.2\] - 2020-03-19

#### Changed

-   Fixed a spurious error message complaining about missing singularity
    image if --use-singularity is not activated.

### \[5.11.1\] - 2020-03-16

#### Changed

-   Fixed a KeyError bug when executing a workflow that defines
    containers without --use-singularity.

## \[5.11.0\] - 2020-03-16

#### Changed

-   Fixes for environment modules and tibanna-based AWS execution.
-   Fixes for --default-resources defaults.
-   --cores is now a mandatory argument!
-   Automatic checksum validation for google storage.

#### Added

-   Azure storage authentication via SAS
-   A generic container directive that will in the future allow for
    other backends than just singularity. This deprecates the
    singularity directive, which will however stay functional at least
    until the next major release.
-   envvars directive for asserting environment variable existence. See
    docs.
-   support for AWS spot instances via --tibanna-config spot=true.
-   Automatic code quality linting via --lint.

## \[5.10.0\] - 2020-01-20

#### Added

-   Jupyter notebook integration, see docs. This enables interactive
    development of certain data analysis parts (e.g. for plotting).
-   Ability to overwrite thread definitions at the command line
    (`--threads rulename=3`), thereby improving scalability.
-   Requester pays configuration for google storage remote files.
-   Add keyword `allow_missing` to expand function, thereby allowing
    partical expansion by skipping wildcards for which no keywords are
    defined.

#### Changed

-   Various bug fixes, e.g. for between workflow caching and script
    execution.

### \[5.9.1\] - 2019-12-20

#### Changed

-   Added a missing module.

## \[5.9.0\] - 2019-12-20

#### Added

-   Support for per-rule environment module definitions to enable HPC
    specific software deployment (see docs).
-   Allow custom log handler defitions via --log-handler-script (e.g.
    post errors and progress to a slack channel or send emails).

-   Allow setting threads as a function of the given cores (see docs).

#### Changed

- Various minor fixes.

### \[5.8.2\] - 2019-12-16

#### Added

- Implemented a `multiext` helper, allowing to define a set of output
files that just differ by extension.

#### Changed

- Fixed a failure when caching jobs with conda environments.
- Fixed various minor bugs.
- Caching now allows to cache the output of rules using `multiext`.

### \[5.8.1\] - 2019-11-15

#### Changed

-   Fixed a bug by adding a missing module.

## \[5.8.0\] - 2019-11-15

#### Added

-   Blockchain based caching between workflows (in collaboration with
    Sven Nahnsen from QBiC), see [the
    docs](https://snakemake.readthedocs.io/en/v5.8.0/executing/caching.html).
-   New flag --skip-cleanup-scripts, that leads to temporary scripts
    (coming from script or wrapper directive) are not deleted (by Vanessa
    Sochat).

#### Changed

- Various bug fixes.

### \[5.7.4\] - 2019-10-23

#### Changed

-   Various fixes and adaptations in the docker container image and the
    test suite.

### \[5.7.1\] - 2019-10-16

#### Added

- Ability to print log files of failed jobs with --show-failed-logs.

#### Changed

- Fixed bugs in tibanna executor.
- Fixed handling of symbolic links.
- Fixed typos in help texts.
- Fixed handling of default resources.
- Fixed bugs in azure storage backend.

## \[5.7.0\] - 2019-10-07

#### Changed

-   Fixed various corner case bugs. Many thanks to the community for
    pull requests and reporting!
-   Container execution adapted to latest singularity.

#### Added

-   First class support for Amazon cloud execution via a new
    [Tibanna backend](https://snakemake.readthedocs.io/en/v5.7.0/executable.html#executing-a-snakemake-workflow-via-tibanna-on-amazon-web-services).
    Thanks to Soo Lee from Harvard Biomedical Informatics!
-   Allow multiple config files to be passed via the command line.
-   A new, more detailed way to visualize the DAG (--filegraph). Thanks
    to Henning Timm!
-   Pathlib compatibility added. Input and output files can now also be
    Path objects. Thanks to Frederik Boulund!
-   New azure storage remote provider. Transparently access input and
    output files on Microsoft Azure. Thanks to Sebastian Kurscheid!

## \[5.6.0\] - 2019-09-06

#### Changed

-   Fix compatibility with latest singularity versions.
-   Various bug fixes (e.g. in cluster error handling, remote providers,
    kubernetes backend).

#### Added

- Add --default-resources flag, that
  allows to define default resources for jobs (e.g. mem_mb, disk_mb), see
  [docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#resources).
- Accept `--dry-run` as a synonym of `--dryrun`. Other Snakemake options
  are similarly hyphenated, so other documentation now refers to
  `--dry-run` but both (and also `-n`) will always be accepted
  equivalently.

### \[5.5.4\] - 2019-07-21

#### Changed

-   Reports now automatically include workflow code and configuration
    for improved transparency.

### \[5.5.3\] - 2019-07-11

#### Changed

-   Various bug fixes.
-   Polished reports.

### \[5.5.2\] - 2019-06-25

#### Changed

-   Various minor bug fixes in reports.
-   Speed improvements when using checkpoints.

### \[5.5.1\] - 2019-06-18

#### Changed

-   Improved report interface. In particular for large files.
-   Small TSV tables are automatically rendered as HTML with datatables.
-   Be more permissive with Snakefile choices: allow "Snakefile",
    "snakefile", "workflow/Snakefile", "workflow/snakefile".

## \[5.5.0\] - 2019-05-31

#### Added

- Script directives now also support Julia.

#### Changed

- Various small bug fixes.

### \[5.4.5\] - 2019-04-12

#### Changed

-   Fixed a bug with pipe output.
-   Cleaned up error output.

### \[5.4.4\] - 2019-03-22

#### Changed

-   Vastly improved performance of HTML reports generated with --report,
    via a more efficient encoding of dara-uri based download links.
-   Tighter layout, plus thumbnails and a lightbox for graphical results
    in HTML reports.
-   Bug fix for pipe groups.
-   Updated docs.
-   Better error handling in DRMAA executor.

### \[5.4.3\] - 2019-03-11

#### Changed

-   More robust handling of conda environment activation that should
    work with all setups where the conda is available when starting
    snakemake.
-   Fixed bugs on windows.

### \[5.4.2\] - 2019-02-15

#### Changed

-   Fixed a bug where git module cannot be imported from wrapper.

### \[5.4.1\] - 2019-02-14

#### Added

-   Warning when R script is used in combination with conda and R_LIBS
    environment variable is set. This can cause unexpected results and
    should be avoided.

#### Changed

-   Improved quoting of paths in conda commands.
-   Fixed various issues with checkpoints.
-   Improved error messages when combining groups with cluster config.
-   Fixed bugs in group implementation.
-   Fixed singularity in combination with shadow.

## \[5.4.0\] - 2018-12-18

#### Added

-   Snakemake now allows for data-dependent conditional re-evaluation of
    the job DAG via checkpoints. This feature also deprecates the
    `dynamic` flag. See [the
    docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#data-dependent-conditional-execution).

### \[5.3.1\] - 2018-12-06

#### Changed

-   Various fixed bugs and papercuts, e.g., in group handling,
    kubernetes execution, singularity support, wrapper and script usage,
    benchmarking, schema validation.

## \[5.3.0\] - 2018-09-18

#### Added

-   Snakemake workflows can now be exported to CWL via the flag
    --export-cwl, see [the
    docs](https://snakemake.readthedocs.io/en/stable/executing/interoperability.html).

#### Changed

-   Fixed bug in script and wrapper execution when using
    `--use-singularity --use-conda`.
-   Add host argument to S3RemoteProvider.
-   Various minor bug fixes.

### \[5.2.4\] - 2018-09-10

#### Added

-   New command line flag --shadow-prefix

#### Changed

-   Fixed permission issue when using the script directive. This is a
    breaking change for scripts referring to files relative to the
    script directory (see the
    [docs](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#external-scripts)).
-   Fixed various minor bugs and papercuts.
-   Allow URL to local git repo with wrapper directive
    (`git+file:///path/to/your/repo/path_to_file@@version`)

### \[5.2.2\] - 2018-08-01

#### Changed

-   Always print timestamps, removed the --timestamps CLI option.
-   more robust detection of conda command
-   Fixed bug in RMarkdown script execution.
-   Fixed a bug in detection of group jobs.

## \[5.2.0\] - 2018-06-28

#### Changed

-   Directory outputs have to marked with `directory`. This ensures
    proper handling of timestamps and cleanup. This is a breaking
    change. Implemented by Rasmus Ågren.
-   Fixed kubernetes tests, fixed kubernetes volume handling.
    Implemented by Andrew Schriefer.
-   jinja2 and networkx are not optional dependencies when installing
    via pip.
-   When conda or singularity directives are used and the corresponding
    CLI flags are not specified, the user is notified at the beginning
    of the log output.
-   Fixed numerous small bugs and papercuts and extended documentation.

### \[5.1.5\] - 2018-06-24

#### Changed

-   fixed missing version info in docker image.
-   several minor fixes to EGA support.

### \[5.1.4\] - 2018-05-28

#### Added

-   Allow `category` to be set.

#### Changed

-   Various cosmetic changes to reports.
-   Fixed encoding issues in reports.

### \[5.1.3\] - 2018-05-22

#### Changed

-   Fixed various bugs in job groups, shadow directive, singularity
    directive, and more.

### \[5.1.2\] - 2018-05-18

#### Changed

-   Fixed a bug in the report stylesheet.

## \[5.1.0\] - 2018-05-17

#### Added

-   A new framework for self-contained HTML reports, including results,
    statistics and topology information. In future releases this will be
    further extended.
-   A new utility snakemake.utils.validate() which allows to validate
    config and pandas data frames using JSON schemas.
-   Two new flags --cleanup-shadow and --cleanup-conda to clean up old
    unused conda and shadow data.

#### Changed

-   Benchmark repeats are now specified inside the workflow via a new
    flag repeat().
-   Command line interface help has been refactored into groups for
    better readability.

## \[5.0.0\] - 2018-05-11

#### Added

-   Group jobs for reduced queuing and network overhead, in particular
    with short running jobs.
-   Output files can be marked as pipes, such that producing and
    consuming job are executed simultaneously and interfomation is
    transferred directly without using disk.
-   Command line flags to clean output files.
-   Command line flag to list files in working directory that are not
    tracked by Snakemake.

#### Changed

-   Fix of --default-remote-prefix in case of input functions returning
    lists or dicts.
-   Scheduler no longer prefers jobs with many downstream jobs.

### \[4.8.1\] - 2018-04-25

#### Added

-   Allow URLs for the conda directive. # Changed
-   Various minor updates in the docs.
-   Several bug fixes with remote file handling.
-   Fix ImportError occuring with script directive.
-   Use latest singularity.
-   Improved caching for file existence checks. We first check existence
    of parent directories and cache these results. By this, large parts
    of the generated FS tree can be pruned if files are not yet present.
    If files are present, the overhead is minimal, since the checks for
    the parents are cached.
-   Various minor bug fixes.

## \[4.8.0\] - 2018-03-13

#### Added

-   Integration with CWL: the `cwl` directive allows to use CWL tool
    definitions in addition to shell commands or Snakemake wrappers.
-   A global `singularity` directive allows to define a global
    singularity container to be used for all rules that don't specify
    their own.
-   Singularity and Conda can now be combined. This can be used to
    specify the operating system (via singularity), and the software
    stack (via conda), without the overhead of creating specialized
    container images for workflows or tasks.

## \[4.7.0\] - 2018-02-19

#### Changed

-   Speedups when calculating dry-runs.
-   Speedups for workflows with many rules when calculating the DAG.
-   Accept SIGTERM to gracefully finish all running jobs and exit.
-   Various minor bug fixes.

## \[4.6.0\] - 2018-02-06

#### Changed

-   Log files can now be used as input files for other rules.
-   Adapted to changes in Kubernetes client API.
-   Fixed minor issues in --archive option.
-   Search path order in scripts was changed to fix a bug with leaked
    packages from root env when using script directive together with
    conda.

### \[4.5.1\] - 2018-02-01

#### Added

-   Input and output files can now tag pathlib objects. # ## Changed
-   Various minor bug fixes.

## \[4.5.0\] - 2018-01-18

#### Added

-   iRODS remote provider # ## Changed
-   Bug fix in shell usage of scripts and wrappers.
-   Bug fixes for cluster execution, --immediate-submit and
    subworkflows.

## \[4.4.0\] - 2017-12-21

#### Added

-   A new shadow mode (minimal) that only symlinks input files has been
    added.

#### Changed

-   The default shell is now bash on linux and macOS. If bash is not
    installed, we fall back to sh. Previously, Snakemake used the
    default shell of the user, which defeats the purpose of portability.
    If the developer decides so, the shell can be always overwritten
    using shell.executable().
-   Snakemake now requires Singularity 2.4.1 at least (only when running
    with --use-singularity).
-   HTTP remote provider no longer automatically unpacks gzipped files.
-   Fixed various smaller bugs.

### \[4.3.1\] - 2017-11-16

#### Added

-   List all conda environments with their location on disk via
    --list-conda-envs.

#### Changed

-   Do not clean up shadow on dry-run.
-   Allow R wrappers.

## \[4.3.0\] - 2017-10-27

#### Added

-   GridFTP remote provider. This is a specialization of the GFAL remote
    provider that uses globus-url-copy to download or upload files. # ##
    Changed
-   Scheduling and execution mechanisms have undergone a major revision
    that removes several potential (but rare) deadlocks.
-   Several bugs and corner cases of the singularity support have been
    fixed.
-   Snakemake now requires singularity 2.4 at least.

## \[4.2.0\] - 2017-10-10

#### Added

-   Support for executing jobs in per-rule singularity images. This is
    meant as an alternative to the conda directive (see docs), providing
    even more guarantees for reproducibility.

#### Changed

-   In cluster mode, jobs that are still running after Snakemake has
    been killed are automatically resumed.
-   Various fixes to GFAL remote provider.
-   Fixed --summary and --list-code-changes.
-   Many other small bug fixes.

## \[4.1.0\] - 2017-09-26

#### Added

-   Support for configuration profiles. Profiles allow to specify
    default options, e.g., a cluster submission command. They can be
    used via 'snakemake --profile myprofile'. See the docs for details.
-   GFAL remote provider. This allows to use GridFTP, SRM and any other
    protocol supported by GFAL for remote input and output files.
-   Added --cluster-status flag that allows to specify a command that
    returns jobs status. # ## Changed
-   The scheduler now tries to get rid of the largest temp files first.
-   The Docker image used for kubernetes support can now be configured
    at the command line.
-   Rate-limiting for cluster interaction has been unified.
-   S3 remote provider uses boto3.
-   Resource functions can now use an additional `attempt` parameter,
    that contains the number of times this job has already been tried.
-   Various minor fixes.

## \[4.0.0\] - 2017-07-24

#### Added

-   Cloud computing support via Kubernetes. Snakemake workflows can be
    executed transparently in the cloud, while storing input and output
    files within the cloud storage (e.g. S3 or Google Storage). I.e.,
    this feature does not need a shared filesystem between the cloud
    notes, and thereby makes the setup really simple.
-   WebDAV remote file support: Snakemake can now read and write from
    WebDAV. Hence, it can now, e.g., interact with Nextcloud or
    Owncloud.
-   Support for default remote providers: define a remote provider to
    implicitly use for all input and output files.
-   Added an option to only create conda environments instead of
    executing the workflow. # ## Changed
-   The number of files used for the metadata tracking of Snakemake
    (e.g., code, params, input changes) in the .snakemake directory has
    been reduced by a factor of 10, which should help with NFS and IO
    bottlenecks. This is a breaking change in the sense that Snakemake
    4.x won't see the metadata of workflows executed with Snakemake 3.x.
    However, old metadata won't be overwritten, so that you can always
    go back and check things by installing an older version of Snakemake
    again.
-   The google storage (GS) remote provider has been changed to use the
    google SDK. This is a breaking change, since the remote provider
    invocation has been simplified (see docs).
-   Due to WebDAV support (which uses asyncio), Snakemake now requires
    Python 3.5 at least.
-   Various minor bug fixes (e.g. for dynamic output files).

### \[3.13.3\] - 2017-06-23

#### Changed

-   Fix a followup bug in Namedlist where a single item was not returned
    as string.

### \[3.13.2\] - 2017-06-20

#### Changed

-   The --wrapper-prefix flag now also affects where the corresponding
    environment definition is fetched from.
-   Fix bug where empty output file list was recognized as containing
    duplicates (issue #574).

### \[3.13.1\] - 2017-06-20

#### Changed

-   Fix --conda-prefix to be passed to all jobs.
-   Fix cleanup issue with scripts that fail to download.

## \[3.13.0\] - 2017-06-12

#### Added

-   An NCBI remote provider. By this, you can seamlessly integrate any
    NCBI resouce (reference genome, gene/protein sequences, ...) as
    input file. # ## Changed
-   Snakemake now detects if automatically generated conda environments
    have to be recreated because the workflow has been moved to a new
    path.
-   Remote functionality has been made more robust, in particular to
    avoid race conditions.
-   `--config` parameter evaluation has been fixed for non-string types.
-   The Snakemake docker container is now based on the official debian
    image.

## \[3.12.0\] - 2017-05-09

#### Added

-   Support for RMarkdown (.Rmd) in script directives.
-   New option --debug-dag that prints all decisions while building the
    DAG of jobs. This helps to debug problems like cycles or unexpected
    MissingInputExceptions.
-   New option --conda-prefix to specify the place where conda
    environments are stored.

#### Changed

-   Benchmark files now also include the maximal RSS and VMS size of the
    Snakemake process and all sub processes.
-   Speedup conda environment creation.
-   Allow specification of DRMAA log dir.
-   Pass cluster config to subworkflow.

### \[3.11.2\] - 2017-03-15

#### Changed

-   Fixed fix handling of local URIs with the wrapper directive.

### \[3.11.1\] - 2017-03-14

#### Changed

-   --touch ignores missing files
-   Fixed handling of local URIs with the wrapper directive.

## \[3.11.0\] - 2017-03-08

#### Added

-   Param functions can now also refer to threads. # ## Changed
-   Improved tutorial and docs.
-   Made conda integration more robust.
-   None is converted to NULL in R scripts.

### \[3.10.2\] - 2017-02-28

#### Changed

-   Improved config file handling and merging.
-   Output files can be referred in params functions (i.e. lambda
    wildcards, output: ...)
-   Improved conda-environment creation.
-   Jobs are cached, leading to reduced memory footprint.
-   Fixed subworkflow handling in input functions.

## \[3.10.0\] - 2017-01-18

#### Added

-   Workflows can now be archived to a tarball with
    `snakemake --archive my-workflow.tar.gz`. The archive contains all
    input files, source code versioned with git and all software
    packages that are defined via conda environments. Hence, the archive
    allows to fully reproduce a workflow on a different machine. Such an
    archive can be uploaded to Zenodo, such that your workflow is
    secured in a self-contained, executable way for the future. # ##
    Changed
-   Improved logging.
-   Reduced memory footprint.
-   Added a flag to automatically unpack the output of input functions.
-   Improved handling of HTTP redirects with remote files.
-   Improved exception handling with DRMAA.
-   Scripts referred by the script directive can now use locally defined
    external python modules.

### \[3.9.1\] - 2016-12-23

#### Added

-   Jobs can be restarted upon failure (--restart-times). # ## Changed
-   The docs have been restructured and improved. Now available under
    snakemake.readthedocs.org.
-   Changes in scripts show up with --list-code-changes.
-   Duplicate output files now cause an error.
-   Various bug fixes.

## \[3.9.0\] - 2016-11-15

#### Added

-   Ability to define isolated conda software environments (YAML) per
    rule. Environments will be deployed by Snakemake upon workflow
    execution.
-   Command line argument --wrapper-prefix in order to overwrite the
    default URL for looking up wrapper scripts. # ## Changed
-   --summary now displays the log files correspoding to each output
    file.
-   Fixed hangups when using run directive and a large number of jobs
-   Fixed pickling errors with anonymous rules and run directive.
-   Various small bug fixes

### \[3.8.2\] - 2016-09-23

#### Changed

-   Add missing import in rules.py.
-   Use threading only in cluster jobs.

### \[3.8.1\] - 2016-09-14

#### Changed

-   Snakemake now warns when using relative paths starting with "./".
-   The option -R now also accepts an empty list of arguments.
-   Bug fix when handling benchmark directive.
-   Jobscripts exit with code 1 in case of failure. This should improve
    the error messages of cluster system.
-   Fixed a bug in SFTP remote provider.

## \[3.8.0\] - 2016-08-26

#### Added

-   Wildcards can now be constrained by rule and globally via the new
    `wildcard_constraints` directive (see the
    [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wildcards)).
-   Subworkflows now allow to overwrite their config file via the
    configfile directive in the calling Snakefile.
-   A method `log_fmt_shell` in the snakemake proxy object that is
    available in scripts and wrappers allows to obtain a formatted
    string to redirect logging output from STDOUT or STDERR.
-   Functions given to resources can now optionally contain an
    additional argument `input` that refers to the input files.
-   Functions given to params can now optionally contain additional
    arguments `input` (see above) and `resources`. The latter refers to
    the resources.
-   It is now possible to let items in shell commands be automatically
    quoted (see the
    [docs](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-rules)).
    This is usefull when dealing with filenames that contain
    whitespaces.

#### Changed

-   Snakemake now deletes output files before job exection. Further, it
    touches output files after job execution. This solves various
    problems with slow NFS filesystems.
-   A bug was fixed that caused dynamic output rules to be executed
    multiple times when forcing their execution with -R.
-   A bug causing double uploads with remote files was fixed. Various
    additional bug fixes related to remote files.
-   Various minor bug fixes.

### \[3.7.1\] - 2016-05-16

#### Changed

-   Fixed a missing import of the multiprocessing module.

## \[3.7.0\] - 2016-05-05

#### Added

-   The entries in `resources` and the `threads` job attribute can now
    be callables that must return `int` values.
-   Multiple `--cluster-config` arguments can be given to the Snakemake
    command line. Later one override earlier ones.
-   In the API, multiple `cluster_config` paths can be given as a list,
    alternatively to the previous behaviour of expecting one string for
    this parameter.
-   When submitting cluster jobs (either through `--cluster` or
    `--drmaa`), you can now use `--max-jobs-per-second` to limit the
    number of jobs being submitted (also available through Snakemake
    API). Some cluster installations have problems with too many jobs
    per second.
-   Wildcard values are now printed upon job execution in addition to
    input and output files. # ## Changed
-   Fixed a bug with HTTP remote providers.

### \[3.6.1\] - 2016-04-08

#### Changed

-   Work around missing RecursionError in Python \< 3.5
-   Improved conversion of numpy and pandas data structures to R
    scripts.
-   Fixed locking of working directory.

## \[3.6.0\] - 2016-03-10

#### Added

-   onstart handler, that allows to add code that shall be only executed
    before the actual workflow execution (not on dryrun).
-   Parameters defined in the cluster config file are now accessible in
    the job properties under the key "cluster".
-   The wrapper directive can be considered stable. # ## Changed
-   Allow to use rule/job parameters with braces notation in cluster
    config.
-   Show a proper error message in case of recursion errors.
-   Remove non-empty temp dirs.
-   Don't set the process group of Snakemake in order to allow kill
    signals from parent processes to be propagated.
-   Fixed various corner case bugs.
-   The params directive no longer converts a list `l` implicitly to
    `" ".join(l)`.

### \[3.5.5\] - 2016-01-23

#### Added

-   New experimental wrapper directive, which allows to refer to
    re-usable [wrapper
    scripts](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-wrappers).
    Wrappers are provided in the [Snakemake Wrapper
    Repository](https://bitbucket.org/snakemake/snakemake-wrappers).
-   David Koppstein implemented two new command line options to
    constrain the execution of the DAG of job to sub-DAGs (--until and
    --omit-from). # ## Changed
-   Fixed various bugs, e.g. with shadow jobs and --latency-wait.

### \[3.5.4\] - 2015-12-04

#### Changed

-   The params directive now fully supports non-string parameters.
    Several bugs in the remote support were fixed.

### \[3.5.3\] - 2015-11-24

#### Changed

-   The missing remote module was added to the package.

### \[3.5.2\] - 2015-11-24

#### Added

-   Support for easy integration of external R and Python scripts via
    the new [script
    directive](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-external-scripts).
-   Chris Tomkins-Tinch has implemented support for remote files:
    Snakemake can now handle input and output files from Amazon S3,
    Google Storage, FTP, SFTP, HTTP and Dropbox.
-   Simon Ye has implemented support for sandboxing jobs with [shadow
    rules](https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-shadow-rules).

#### Changed

-   Manuel Holtgrewe has fixed dynamic output files in combination with
    multiple wildcards.
-   It is now possible to add suffixes to all shell commands with
    shell.suffix("mysuffix").
-   Job execution has been refactored to spawn processes only when
    necessary, resolving several problems in combination with huge
    workflows consisting of thousands of jobs and reducing the memory
    footprint.
-   In order to reflect the new collaborative development model,
    Snakemake has moved from my personal bitbucket account to
    <http://snakemake.bitbucket.org>.

### \[3.4.2\] - 2015-09-12

#### Changed

-   Willem Ligtenberg has reduced the memory usage of Snakemake.
-   Per Unneberg has improved config file handling to provide a more
    intuitive overwrite behavior.
-   Simon Ye has improved the test suite of Snakemake and helped with
    setting up continuous integration via Codeship.
-   The cluster implementation has been rewritten to use only a single
    thread to wait for jobs. This avoids failures with large numbers of
    jobs.
-   Benchmarks are now writing tab-delimited text files instead of JSON.
-   Snakemake now always requires to set the number of jobs with -j when
    in cluster mode. Set this to a high value if your cluster does not
    have restrictions.
-   The Snakemake Conda package has been moved to the bioconda channel.
-   The handling of Symlinks was improved, which made a switch to Python
    3.3 as the minimum required Python version necessary.

### \[3.4.1\] - 2015-08-05

#### Changed

-   This release fixes a bug that caused named input or output files to
    always be returned as lists instead of single files.

## \[3.4\] - 2015-07-18

#### Added

-   This release adds support for executing jobs on clusters in
    synchronous mode (e.g. qsub -sync). Thanks to David Alexander for
    implementing this.
-   There is now vim syntax highlighting support (thanks to Jay
    Hesselberth).
-   Snakemake is now available as Conda package.

#### Changed

-   Lots of bugs have been fixed. Thanks go to e.g. David Koppstein,
    Marcel Martin, John Huddleston and Tao Wen for helping with useful
    reports and debugging.

See [here](https://bitbucket.org/snakemake/snakemake/wiki/News-Archive)
for older changes.
[![Gitpod Ready-to-Code](https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod)](https://gitpod.io/#https://github.com/snakemake/snakemake)
[![test status](https://github.com/snakemake/snakemake/workflows/CI/badge.svg?branch=main)](https://github.com/snakemake/snakemake/actions?query=branch%3Amain+workflow%3ACI)
[![Sonarcloud Status](https://sonarcloud.io/api/project_badges/measure?project=snakemake_snakemake&metric=alert_status)](https://sonarcloud.io/dashboard?id=snakemake_snakemake)
[![Bioconda](https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda)](https://bioconda.github.io/recipes/snakemake/README.html)
[![Pypi](https://img.shields.io/pypi/pyversions/snakemake.svg)](https://pypi.org/project/snakemake)
[![docker container status](https://img.shields.io/github/workflow/status/snakemake/snakemake/Publish%20to%20Docker%20Hub?color=blue&label=docker%20container)](https://hub.docker.com/r/snakemake/snakemake)
[![Stack Overflow](https://img.shields.io/badge/stack-overflow-orange.svg)](https://stackoverflow.com/questions/tagged/snakemake)
[![Twitter](https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow)](https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester)
[![Discord](https://img.shields.io/discord/753690260830945390?label=discord%20chat)](https://discord.gg/NUdMtmr)
[![Github stars](https://img.shields.io/github/stars/snakemake/snakemake?style=social)](https://github.com/snakemake/snakemake/stargazers)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) 

# Snakemake

The Snakemake workflow management system is a tool to create **reproducible and scalable** data analyses.
Snakemake is highly popular, with on average more than 6 new citations per week, and over 200k downloads.
Workflows are described via a human readable, Python based language.
They can be seamlessly scaled to server, cluster, grid and cloud environments without the need to modify the workflow definition.
Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.

**Homepage: https://snakemake.github.io**

Copyright (c) 2012-2021 Johannes Köster <johannes.koester@uni-due.com> (see LICENSE)

# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

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
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
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
reported via email to johannes.koester@uni-due.de.
All complaints will be reviewed and investigated promptly and fairly.

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

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

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
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.

# Test Remote iRODS

These are the instructions for testing Snakemake remote iRODS support locally
on your computer.

## Prerequisites

1. The latest Docker installation (a.k.a. not the one from the Ubuntu 16.04
repository). Follow the instructions on
https://docs.docker.com/engine/installation/linux/docker-ce/ubuntu/
2. The irods icommand tools as described in https://packages.irods.org/ to
setup the repository and `apt-get install irods-icommands` for installation.
The password file has to be generated with `iinit`. A valid environment file
is located in `setup-data/irods_environment.json` (in this case the
authentication file is expected in `~/.irods/.irodsA`). This is necessary,
because the obfuscation of the password uses the uid, so the `.irodsA` file
can't be shipped.
3. Snakemake

## Build and run the Docker container

```
make run
```

## Stop and delete the Docker container

```
make stop
```

## Run the test

```
snakemake -s Snakefile.local
```

## Touch the input file (for a new test)

```
make touch
```

## Show iRODS content

```
make ls
```

## Example

```
make run
make ls
snakemake -s Snakefile.local
make ls
snakemake -s Snakefile.local
# nothing to do here
make touch
snakemake -s Snakefile.local
```
```
# the following has been performed by aws admin already
tibanna deploy_unicorn --usergroup=johannes --buckets=snakemake-tibanna-test,snakemake-tibanna-test2
```
```
# run the following to do a test run
pip install -U tibanna
python cleanup.py  # first delete output files that are already on s3
export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=tibanna_unicorn_johannes
snakemake --tibanna --use-conda --configfile=config.json --default-remote-prefix=snakemake-tibanna-test/1
```

# Executing this test case

To run this test, you need a running kubernetes setup.
For google cloud, see [here](https://snakemake.readthedocs.io/en/stable/executing/cloud.html#setup-kubernetes-on-google-cloud-engine).
With this, you can execute in case of google cloud:

    snakemake --kubernetes --use-conda --default-remote-provider GS --default-remote-prefix my-bucket

while replacing ``my-bucket`` with your storage bucket. The same test should also work on amazon (given that kubernetes is setup):

    snakemake --kubernetes --use-conda --default-remote-provider S3 --default-remote-prefix my-bucket
# Instruction for testing of Azure Storage integration

In order to perform this test, you need an Azure Storage account
with read/write access.
Both the storage account and associated key or SAS token (without
leading questionmark) need to be
passed to snakemake at runtime, by exporting
environment variables `AZURE_ACCOUNT` and either `AZURE_KEY` or
`SAS_TOKEN`.

Furthermore, in the storage account, a container "snakemake-test"
needs to be created prior to running the test.

And lastly, a local file called `test.txt.gz` needs to be created.

### Description

<!--Add a description of your PR here-->

### QC
<!-- Make sure that you can tick the boxes below. -->

* [ ] The PR contains a test case for the changes or the changes are already covered by an existing test case.
* [ ] The documentation (`docs/`) is updated to reflect the changes or this is not necessary (e.g. if the change does neither modify the language nor the behavior or functionalities of Snakemake).
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

<!-- Please do not post usage questions here. Ask them on Stack Overflow: https://stackoverflow.com/questions/tagged/snakemake -->

**Snakemake version**
<!--Note the Snakemake version for which you experience the bug.
Please only report bugs of the **latest stable release of Snakemake**.
If possible please check whether the bug has been already fixed in the main branch.-->

**Describe the bug**
<!--A clear and concise description of what the bug is.-->

**Logs**
<!--If applicable, any terminal output to help explain your problem.-->

**Minimal example**
<!--Add a minimal example for reproducing the bug.-->

**Additional context**
<!--Add any other context about the problem here.-->
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
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
(project_info-history)=

(_changelog)=

```{include} ../../CHANGELOG.md
```
A vim syntax highlighting definition for Snakemake.

To install via Vundle use:

    Plugin 'https://github.com/snakemake/snakemake.git', {'rtp': 'misc/vim/'}

To install via [vim-plug]( https://github.com/junegunn/vim-plug):

    Plug 'snakemake/snakemake', {'rtp': 'misc/vim'}

To manually install, copy `syntax/snakemake.vim` file to `$HOME/.vim/syntax`
directory and `ftdetect/snakemake.vim` file to `$HOME/.vim/ftdetect`.

Highlighting can be forced in a vim session with `:set syntax=snakemake`.
A nano syntax highlighting definition for Snakemake.

To use this file in nano, copy the `syntax/snakemake.nanorc` file
to your `$HOME` directory and add a line to your ~/.nanorc saying:

    include $HOME/snakemake.nanorc


NB. Line 12 of the syntax file contains a regular expression for
identifying a shebang (#!) header line. This command is not
supported in some versions of nano, so is commented out. If you
wish to enable it, simply remove the # at the beginning of the
line. If you are uncertain if your version of nano supports this
command, you may remove it and try. You will see an error if it
does not. Leaving the line commented out will result in the header
being treated as a regular comment for highlighting purposes. 
---
title: "Test Report"
author:
    - "Your Name"
date: "`r format(Sys.time(), '%d %B, %Y')`"
params:
   rmd: "report.Rmd"
output:
  html_document:
  highlight: tango
  number_sections: no
  theme: default
  toc: yes
  toc_depth: 3
  toc_float:
    collapsed: no
    smooth_scroll: yes
---

## R Markdown

This is an R Markdown document.

Test include from snakemake ``r snakemake@input``.


```{r}
print(getwd())
```

```{r}
# This fails
data <- read.table(snakemake@input[[1]])
data
```

## Source
<a download="report.Rmd" href="`r base64enc::dataURI(file = params$rmd, mime = 'text/rmd', encoding = 'base64')`">R Markdown source file (to produce this document)</a>
---
title: "Test Report"
author: "Mattias"
date: "March 22, 2017"
output:
  html_document:
    highlight: null
    number_sections: no
    theme: null
    mathjax: null
---


## R Markdown

This is an R Markdown document.

Test include from snakemake `r snakemake@params[["test"]]`.
Files obtained from a directory. This file starts with {{ snakemake.wildcards.name }}. This value has been dynamically inferred from the given pattern.This is the workflow description. Test reference fig1.svg_.
This is a caption with a `link <https://www.google.com>`_.
Some math :math:`\sum_i i^2`.
An example table.
the caption
Files obtained from a directory. This file starts with {{ snakemake.wildcards.name }}. This value has been dynamically inferred from the given pattern.This is the workflow description. Test reference fig1.svg_.
This is a caption with a `link <https://www.google.com>`_.
Some math :math:`\sum_i i^2`.
An example table.
.. _manual-main:

=========
Snakemake
=========

.. image:: https://img.shields.io/badge/Gitpod-ready--to--code-blue?logo=gitpod
    :target: https://gitpod.io/#https://github.com/snakemake/snakemake

.. image:: https://img.shields.io/conda/dn/bioconda/snakemake.svg?label=Bioconda
    :target: https://bioconda.github.io/recipes/snakemake/README.html

.. image:: https://img.shields.io/pypi/pyversions/snakemake.svg
    :target: https://www.python.org

.. image:: https://img.shields.io/pypi/v/snakemake.svg
    :target: https://pypi.python.org/pypi/snakemake

.. image:: https://img.shields.io/github/workflow/status/snakemake/snakemake/Publish%20to%20Docker%20Hub?color=blue&label=docker%20container&branch=main
    :target: https://hub.docker.com/r/snakemake/snakemake

.. image:: https://github.com/snakemake/snakemake/workflows/CI/badge.svg?branch=main&label=tests
    :target: https://github.com/snakemake/snakemake/actions?query=branch%3Amain+workflow%3ACI

.. image:: https://img.shields.io/badge/stack-overflow-orange.svg
    :target: https://stackoverflow.com/questions/tagged/snakemake

.. image:: https://img.shields.io/twitter/follow/johanneskoester.svg?style=social&label=Follow
    :target: https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester

.. image:: https://img.shields.io/discord/753690260830945390?label=discord%20chat   
    :alt: Discord
    :target: https://discord.gg/NUdMtmr

.. image:: https://img.shields.io/github/stars/snakemake/snakemake?style=social
    :alt: GitHub stars
    :target: https://github.com/snakemake/snakemake/stargazers

.. .. raw:: html
          <span class="__dimensions_badge_embed__" data-doi="https://doi.org/10.1093/bioinformatics/bts480" data-legend="always" data-style="large_rectangle"></span><script async src="https://badge.dimensions.ai/badge.js" charset="utf-8"></script>

The Snakemake workflow management system is a tool to create **reproducible and scalable** data analyses.
Workflows are described via a human readable, Python based language.
They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition.
Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment.

Snakemake is **highly popular**, with `>5 new citations per week <https://badge.dimensions.ai/details/id/pub.1018944052>`_.
For an introduction, please visit https://snakemake.github.io.


.. _main-getting-started:

---------------
Getting started
---------------

* To get a first impression, please visit https://snakemake.github.io.
* To properly understand what Snakemake can do for you please read our `"rolling" paper <https://doi.org/10.12688/f1000research.29032.1>`_.
* News about Snakemake are published via `Twitter <https://twitter.com/search?l=&q=%23snakemake%20from%3Ajohanneskoester>`_.
* To learn Snakemake, please do the :ref:`tutorial`, and see the :ref:`FAQ <project_info-faq>`.
* **Best practices** for writing Snakemake workflows can be found :ref:`here <snakefiles-best_practices>`.
* For more advanced usage on various platforms, see the :ref:`executor_tutorial`.

.. _main-support:

-------
Support
-------

* For releases, see :ref:`Changelog <changelog>`.
* Check :ref:`frequently asked questions (FAQ) <project_info-faq>`.
* In case of **questions**, please post on `stack overflow <https://stackoverflow.com/questions/tagged/snakemake>`_.
* To **discuss** with other Snakemake users, use the `discord server <https://discord.gg/kHvtG6N>`_. **Please do not post questions there. Use stack overflow for questions.**
* For **bugs and feature requests**, please use the `issue tracker <https://github.com/snakemake/snakemake/issues>`_.
* For **contributions**, visit Snakemake on `Github <https://github.com/snakemake/snakemake>`_ and read the :ref:`guidelines <project_info-contributing>`.

--------
Citation
--------
When using Snakemake, please cite our "rolling" paper

`Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33. <https://doi.org/10.12688/f1000research.29032.1>`_

This paper will also be regularly updated when Snakemake receives new features.
See :doc:`Citations <project_info/citations>` for more information.

---------
Resources
---------

`Snakemake Wrappers Repository <https://snakemake-wrappers.readthedocs.org>`_
    The Snakemake Wrapper Repository is a collection of reusable wrappers that allow to quickly use popular tools from Snakemake rules and workflows.

`Snakemake Workflow Catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_
    An automatically scraped catalog of publicly available Snakemake workflows for any kind of data analysis.

`Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_
    This project provides a collection of high quality modularized and re-usable workflows.
    The provided code should also serve as a best-practices of how to build production ready workflows with Snakemake.
    Everybody is invited to contribute.

`Snakemake Profiles Project <https://github.com/snakemake-profiles/doc>`_
    This project provides Snakemake configuration profiles for various execution environments.
    Please consider contributing your own if it is still missing.

`Conda-Forge <https://conda-forge.org>`_
    Conda-Forge is a community driven distribution of Conda packages that can be used from Snakemake for creating completely reproducible workflows by defining the used software versions and providing binaries.

`Bioconda <https://bioconda.github.io/>`_
    Bioconda, a partner project of conda-forge, is a community driven distribution of bioinformatics-related Conda packages that can be used from Snakemake for creating completely reproducible workflows by defining the used software versions and providing binaries.


.. toctree::
   :caption: Getting started
   :name: getting_started
   :hidden:
   :maxdepth: 1

   getting_started/installation
   tutorial/tutorial
   tutorial/short
   executor_tutorial/tutorial
   snakefiles/best_practices

.. toctree::
  :caption: Executing workflows
  :name: execution
  :hidden:
  :maxdepth: 1

  executing/cli
  executing/cluster
  executing/cloud
  executing/grouping
  executing/caching
  executing/interoperability
  executing/monitoring

.. toctree::
    :caption: Defining workflows
    :name: snakefiles
    :hidden:
    :maxdepth: 1

    snakefiles/writing_snakefiles
    snakefiles/rules
    snakefiles/configuration
    snakefiles/modularization
    snakefiles/remote_files
    snakefiles/utils
    snakefiles/deployment
    snakefiles/reporting
    snakefiles/testing
    snakefiles/foreign_wms


.. toctree::
    :caption: API Reference
    :name: api-reference
    :hidden:
    :maxdepth: 1

    api_reference/snakemake
    api_reference/snakemake_utils
    api_reference/internal/modules


.. toctree::
    :caption: Project Info
    :name: project-info
    :hidden:
    :maxdepth: 1

    project_info/citations
    project_info/more_resources
    project_info/faq
    project_info/contributing
    project_info/authors
    project_info/history
    project_info/license
.. tutorial-advanced:

Advanced: Decorating the example workflow
-----------------------------------------

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Conda: https://conda.pydata.org
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: https://www.graphviz.org
.. _RestructuredText: https://docutils.sourceforge.io/docs/user/rst/quickstart.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: https://json.org
.. _YAML: https://yaml.org
.. _DRMAA: https://www.drmaa.org
.. _rpy2: https://rpy2.github.io
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: https://slides.com/johanneskoester/deck-1

Now that the basic concepts of Snakemake have been illustrated, we can introduce some advanced functionality.

Step 1: Specifying the number of used threads
:::::::::::::::::::::::::::::::::::::::::::::

For some tools, it is advisable to use more than one thread in order to speed up the computation.
**Snakemake can be made aware of the threads a rule needs** with the ``threads`` directive.
In our example workflow, it makes sense to use multiple threads for the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

The number of threads can be propagated to the shell command with the familiar braces notation (i.e. ``{threads}``).
If no ``threads`` directive is given, a rule is assumed to need 1 thread.

When a workflow is executed, **the number of threads the jobs need is considered by the Snakemake scheduler**.
In particular, the scheduler ensures that the sum of the threads of all jobs running at the same time does not exceed a given number of available CPU cores.
This number is given with the ``--cores`` command line argument, which is mandatory for ``snakemake`` calls that actually run the workflow.
For example

.. code:: console

    $ snakemake --cores 10


.. sidebar:: Note

  Apart from the very common thread resource, Snakemake provides a ``resources`` directive that can be used to **specify arbitrary resources**, e.g., memory usage or auxiliary computing devices like GPUs.
  Similar to threads, these can be considered by the scheduler when an available amount of that resource is given with the command line argument ``--resources`` (see :ref:`snakefiles-resources`).

would execute the workflow with 10 cores.
Since the rule ``bwa_map`` needs 8 threads, only one job of the rule can run at a time, and the Snakemake scheduler will try to saturate the remaining cores with other jobs like, e.g., ``samtools_sort``.
The threads directive in a rule is interpreted as a maximum: when **less cores than threads** are provided, the number of threads a rule uses will be **reduced to the number of given cores**.

If ``--cores`` is given without a number, all available cores are used.

Exercise
........

* With the flag ``--forceall`` you can enforce a complete re-execution of the workflow. Combine this flag with different values for ``--cores`` and examine how the scheduler selects jobs to run in parallel.

Step 2: Config files
::::::::::::::::::::

So far, we specified which samples to consider by providing a Python list in the Snakefile.
However, often you want your workflow to be customizable, so that it can easily be adapted to new data.
For this purpose, Snakemake provides a `config file mechanism <https://snakemake.readthedocs.io/en/latest/snakefiles/configuration.html>`_.
Config files can be written in JSON_ or YAML_, and are used with the ``configfile`` directive.
In our example workflow, we add the line

.. code:: python

    configfile: "config.yaml"

to the top of the Snakefile.
Snakemake will load the config file and store its contents into a globally available `dictionary <https://docs.python.org/3/tutorial/datastructures.html#dictionaries>`_ named ``config``.
In our case, it makes sense to specify the samples in ``config.yaml`` as

.. code:: yaml

    samples:
        A: data/samples/A.fastq
        B: data/samples/B.fastq

Now, we can remove the statement defining ``SAMPLES`` from the Snakefile and change the rule ``bcftools_call`` to

.. code:: python

    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
        output:
            "calls/all.vcf"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

.. _tutorial-input_functions:

Step 3: Input functions
:::::::::::::::::::::::

Since we have stored the path to the FASTQ files in the config file, we can also generalize the rule ``bwa_map`` to use these paths.
This case is different to the rule ``bcftools_call`` we modified above.
To understand this, it is important to know that Snakemake workflows are executed in three phases.

1. In the **initialization** phase, the files defining the workflow are parsed and all rules are instantiated.
2. In the **DAG** phase, the directed acyclic dependency graph of all jobs is built by filling wildcards and matching input files to output files.
3. In the **scheduling** phase, the DAG of jobs is executed, with jobs started according to the available resources.

The expand functions in the list of input files of the rule ``bcftools_call`` are executed during the initialization phase.
In this phase, we don't know about jobs, wildcard values and rule dependencies.
Hence, we cannot determine the FASTQ paths for rule ``bwa_map`` from the config file in this phase, because we don't even know which jobs will be generated from that rule.
Instead, we need to defer the determination of input files to the DAG phase.
This can be achieved by specifying an **input function** instead of a string as inside of the input directive.
For the rule ``bwa_map`` this works as follows:

.. code:: python

    def get_bwa_map_input_fastqs(wildcards):
        return config["samples"][wildcards.sample]
    
    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        threads: 8
        shell:
            "bwa mem -t {threads} {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

  Snakemake does not automatically rerun jobs when new input files are added as
  in the excercise below. However, you can get a list of output files that
  are affected by such changes with ``snakemake --list-input-changes``.
  To trigger a rerun, this bit of bash magic helps:

  .. code:: console

    snakemake -n --forcerun $(snakemake --list-input-changes)

Any normal function would work as well.
Input functions take as **single argument** a ``wildcards`` object, that allows to access the wildcards values via attributes (here ``wildcards.sample``).
They have to **return a string or a list of strings**, that are interpreted as paths to input files (here, we return the path that is stored for the sample in the config file).
Input functions are evaluated once the wildcard values of a job are determined.


Exercise
........

* In the ``data/samples`` folder, there is an additional sample ``C.fastq``. Add that sample to the config file and see how Snakemake wants to recompute the part of the workflow belonging to the new sample, when invoking with ``snakemake -n --reason --forcerun bcftools_call``.

Step 4: Rule parameters
:::::::::::::::::::::::

Sometimes, shell commands are not only composed of input and output files and some static flags.
In particular, it can happen that additional parameters need to be set depending on the wildcard values of the job.
For this, Snakemake allows to **define arbitrary parameters** for rules with the ``params`` directive.
In our workflow, it is reasonable to annotate aligned reads with so-called read groups, that contain metadata like the sample name.
We modify the rule ``bwa_map`` accordingly:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        threads: 8
        shell:
            "bwa mem -R '{params.rg}' -t {threads} {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

  The ``params`` directive can also take functions like in Step 3 to defer
  initialization to the DAG phase. In contrast to input functions, these can
  optionally take additional arguments ``input``, ``output``, ``threads``, and ``resources``.

Similar to input and output files, ``params`` can be accessed from the shell command, the Python based ``run`` block, or the script directive (see :ref:`tutorial-script`).

Exercise
........

* Variant calling can consider a lot of parameters. A particularly important one is the prior mutation rate (1e-3 per default). It is set via the flag ``-P`` of the ``bcftools call`` command. Consider making this flag configurable via adding a new key to the config file and using the ``params`` directive in the rule ``bcftools_call`` to propagate it to the shell command.

Step 5: Logging
:::::::::::::::

When executing a large workflow, it is usually desirable to store the logging output of each job into a separate file, instead of just printing all logging output to the terminal---when multiple jobs are run in parallel, this would result in chaotic output.
For this purpose, Snakemake allows to **specify log files** for rules.
Log files are defined via the ``log`` directive and handled similarly to output files, but they are not subject of rule matching and are not cleaned up when a job fails.
We modify our rule ``bwa_map`` as follows:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            "mapped_reads/{sample}.bam"
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

.. sidebar:: Note

  It is best practice to store all log files in a subdirectory ``logs/``, prefixed by the rule or tool name.

The shell command is modified to `collect STDERR output <https://tldp.org/LDP/abs/html/io-redirection.html>`_ of both ``bwa`` and ``samtools`` and pipe it into the file referred to by ``{log}``.
Log files must contain exactly the same wildcards as the output files to avoid file name clashes between different jobs of the same rule.

Exercise
........

* Add a log directive to the ``bcftools_call`` rule as well.
* Time to re-run the whole workflow (remember the command line flags to force re-execution). See how log files are created for variant calling and read mapping.
* The ability to track the provenance of each generated result is an important step towards reproducible analyses. Apart from the ``report`` functionality discussed before, Snakemake can summarize various provenance information for all output files of the workflow. The flag ``--summary`` prints a table associating each output file with the rule used to generate it, the creation date and optionally the version of the tool used for creation is provided. Further, the table informs about updated input files and changes to the source code of the rule after creation of the output file. Invoke Snakemake with ``--summary`` to examine the information for our example.

.. _tutorial_temp-and-protected-files:

Step 6: Temporary and protected files
:::::::::::::::::::::::::::::::::::::

In our workflow, we create two BAM files for each sample, namely
the output of the rules ``bwa_map`` and ``samtools_sort``.
When not dealing with examples, the underlying data is usually huge.
Hence, the resulting BAM files need a lot of disk space and their creation takes some time.
To save disk space, you can **mark output files as temporary**.
Snakemake will delete the marked files for you, once all the consuming jobs (that need it as input) have been executed.
We use this mechanism for the output file of the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

This results in the deletion of the BAM file once the corresponding ``samtools_sort`` job has been executed.
Since the creation of BAM files via read mapping and sorting is computationally expensive, it is reasonable to **protect** the final BAM file **from accidental deletion or modification**.
We modify the rule ``samtools_sort`` to mark its output file as ``protected``:

.. code:: python

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            protected("sorted_reads/{sample}.bam")
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

After successful execution of the job, Snakemake will write-protect the output file in the filesystem, so that it can't be overwritten or deleted by accident.

Exercise
........

* Re-execute the whole workflow and observe how Snakemake handles the temporary and protected files.
* Run Snakemake with the target ``mapped_reads/A.bam``. Although the file is marked as temporary, you will see that Snakemake does not delete it because it is specified as a target file.
* Try to re-execute the whole workflow again with the dry-run option. You will see that it fails (as intended) because Snakemake cannot overwrite the protected output files.

Summary
:::::::

For this advanced part of the tutorial, we have now created a ``config.yaml`` configuration file:

.. code:: yaml

    samples:
        A: data/samples/A.fastq
        B: data/samples/B.fastq
    
    prior_mutation_rate: 0.001


With this, the final version of our workflow in the ``Snakefile`` looks like this:

.. code:: python

    configfile: "config.yaml"


    rule all:
        input:
            "plots/quals.svg"


    def get_bwa_map_input_fastqs(wildcards):
        return config["samples"][wildcards.sample]
    

    rule bwa_map:
        input:
            "data/genome.fa",
            get_bwa_map_input_fastqs
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"


    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            protected("sorted_reads/{sample}.bam")
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"


    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"


    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=config["samples"]),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=config["samples"])
        output:
            "calls/all.vcf"
        params:
            rate=config["prior_mutation_rate"]
        log:
            "logs/bcftools_call/all.log"
        shell:
            "(samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv -P {params.rate} - > {output}) 2> {log}"


    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        script:
            "scripts/plot-quals.py"

.. _tutorial-setup:

Setup
-----

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Graphviz: https://www.graphviz.org
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Jinja2: https://jinja.palletsprojects.com
.. _NetworkX: https://networkx.github.io
.. _Matplotlib: https://matplotlib.org
.. _Pysam: https://pysam.readthedocs.io
.. _Bioconda: https://bioconda.github.io
.. _WSL: https://docs.microsoft.com/en-us/windows/wsl/about
.. _WSL Documentation: https://docs.microsoft.com/en-us/windows/wsl/install-win10
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html

Requirements
::::::::::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ ≥5.24.1
* BWA_ 0.7
* SAMtools_ 1.9
* Pysam_ 0.15
* BCFtools_ 1.9
* Graphviz_ 2.42
* Jinja2_ 2.11
* NetworkX_ 2.5
* Matplotlib_ 3.3

However, don't install any of these this manually now, we guide you through better ways below.

Run tutorial for free in the cloud via Gitpod
:::::::::::::::::::::::::::::::::::::::::::::

.. sidebar:: Note

    A common thing to happen while using the development environment in GitPod is to hit ``Ctrl-s`` while in the terminal window, because you wanted to save a file in the editor window.
    This will freeze up you terminal.
    To get it back, make sure you selected the terminal window by clicking on it and then hit ``Ctrl-q``.

The easiest way to run this tutorial is to use Gitpod, which enables performing the excercises via your browser---including all required software, for free and in the cloud.
In order to do this, simply open the predefined `snakemake-tutorial GitPod workspace <https://gitpod.io/#https://github.com/snakemake/snakemake-tutorial-data>`_ in your browser.
GitPod provides you with a `Theia development environment <https://theia-ide.org/docs>`_, which you can learn about in the linked documentation.
Once you have a basic understanding of this environment, you can go on directy with :ref:`tutorial-basics`.

Running the tutorial on your local machine
::::::::::::::::::::::::::::::::::::::::::

If you prefer to run the tutorial on your local machine, please follow the steps below.

The easiest way to set these prerequisites up, is to use the Mambaforge_ Python 3 distribution 
(Mambaforge_ is a Conda based distribution like Miniconda_, which however uses Mamba_ a fast and more robust replacement for the Conda_ package manager).
The tutorial assumes that you are using either Linux or MacOS X.
Both Snakemake and Mambaforge_ work also under Windows, but the Windows shell is too different to be able to provide generic examples.

Setup on Windows
::::::::::::::::

If you already use Linux or MacOS X, go on with **Step 1**.

Windows Subsystem for Linux
"""""""""""""""""""""""""""

If you use Windows 10, you can set up the Windows Subsystem for Linux (`WSL`_) to natively run linux applications.
Install the WSL following the instructions in the `WSL Documentation`_. You can chose any Linux distribution available for the WSL, but the most popular and accessible one is Ubuntu.
Start the WSL and set up your account; now, you can follow the steps of our tutorial from within your Linux environment in the WSL.

Vagrant virtual machine
"""""""""""""""""""""""

If you are using a version of Windows older than 10 or if you do not wish to install the WSL, you can instead setup a Linux virtual machine (VM) with Vagrant_.
First, install Vagrant following the installation instructions in the `Vagrant Documentation`_.
Then, create a new directory you want to share with your Linux VM, for example, create a folder named ``vagrant-linux`` somewhere.
Open a command line prompt, and change into that directory.
Here, you create a 64-bit Ubuntu Linux environment with

.. code:: console

    > vagrant init hashicorp/precise64
    > vagrant up

If you decide to use a 32-bit image, you will need to download the 32-bit version of Miniconda in the next step.
The contents of the ``vagrant-linux`` folder will be shared with the virtual machine that is set up by vagrant.
You can log into the virtual machine via

.. code:: console

    > vagrant ssh

If this command tells you to install an SSH client, you can follow the instructions in this Blogpost_.
Now, you can follow the steps of our tutorial from within your Linux VM.


Step 1: Installing Mambaforge
:::::::::::::::::::::::::::::

First, please **open a terminal** or make sure you are logged into your Vagrant Linux VM.
Assuming that you have a 64-bit system, on Linux, download and install Miniconda 3 with

.. code:: console

    $ wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh
    $ bash Mambaforge-Linux-x86_64.sh

On MacOS with x86_64 architecture, download and install with

.. code:: console

    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-x86_64.sh -o Mambaforge-MacOSX-x86_64.sh
    $ bash Mambaforge-MacOSX-x86_64.sh

On MacOS with ARM/M1 architecture, download and install with

.. code:: console

    $ curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-MacOSX-arm64.sh -o Mambaforge-MacOSX-arm64.sh
    $ bash Mambaforge-MacOSX-arm64.sh

When you are asked the question

.. code::

    Do you wish the installer to prepend the install location to PATH ...? [yes|no]

answer with **yes**.
Along with a minimal Python 3 environment, Mambaforge contains the package manager Mamba_.
After closing your current terminal and opening a **new terminal**, you can use the new ``conda`` command to install software packages and create isolated environments to, for example, use different versions of the same package.
We will later use Conda_ to create an isolated environment with all the required software for this tutorial.

Step 2: Preparing a working directory
:::::::::::::::::::::::::::::::::::::

First, **create a new directory** ``snakemake-tutorial`` at a **place you can easily remember** and change into that directory in your terminal:

.. code:: console

    $ mkdir snakemake-tutorial
    $ cd snakemake-tutorial

If you use a Vagrant Linux VM from Windows as described above, create that directory under ``/vagrant/``, so that the contents are shared with your host system (you can then edit all files from within Windows with an editor that supports Unix line breaks).
Then, **change to the newly created directory**.
In this directory, we will later create an example workflow that illustrates the Snakemake syntax and execution environment.
First, we download some example data on which the workflow shall be executed:

.. code:: console

    $ wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.24.1.tar.gz
    $ tar --wildcards -xf v5.24.1.tar.gz --strip 1 "*/data" "*/environment.yaml"

This will create a folder ``data`` and a file ``environment.yaml`` in the working directory.
If your tar command does not provide a ``--wildcards`` flag, you can also just unpack the file without it (which will just leave some more unneeded files in the working directory).

Step 3: Creating an environment with the required software
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

First, make sure to activate the conda base environment with

.. code:: console

    $ conda activate base

The ``environment.yaml`` file that you have obtained with the previous step (Step 2) can be used to install all required software into an isolated Conda environment with the name ``snakemake-tutorial`` via

.. code:: console

    $ mamba env create --name snakemake-tutorial --file environment.yaml

If you don't have the Mamba_ command because you used a different conda distribution than Mambaforge_, you can also first install Mamba_ 
(which is a faster and more robust replacement for Conda_) in your base environment with

.. code:: console

    $ conda install -n base -c conda-forge mamba

and then run the `mamba env create` command shown above.

Step 4: Activating the environment
::::::::::::::::::::::::::::::::::

To activate the ``snakemake-tutorial`` environment, execute

.. code:: console

    $ conda activate snakemake-tutorial

Now you can use the installed tools.
Execute

.. code:: console

    $ snakemake --help

to test this and get information about the command-line interface of Snakemake.
To exit the environment, you can execute

.. code:: console

    $ conda deactivate

but **don't do that now**, since we finally want to start working with Snakemake :-).
Short tutorial
==============

Here we provide a short tutorial that guides you through the main features of Snakemake.
Note that this is not suited to learn Snakemake from scratch, rather to give a first impression.
To really learn Snakemake (starting from something simple, and extending towards advanced features), use the main :ref:`tutorial`.

This document shows all steps performed in the official `Snakemake live demo <https://youtu.be/hPrXcUUp70Y>`_,
such that it becomes possible to follow them at your own pace.
Solutions to each step can be found at the bottom of this document.

The examples presented in this tutorial come from Bioinformatics.
However, Snakemake is a general-purpose workflow management system for any discipline.
For an explanation of the steps you will perform here, have a look at :ref:`tutorial-background`.
More thorough explanations are provided in the full :ref:`tutorial`.


Prerequisites
-------------

First, install Snakemake via Conda, as outlined in :ref:`conda-install`.
The minimal version of Snakemake is sufficient for this demo.

Second, download and unpack the test data needed for this example from
`here <https://github.com/snakemake/snakemake-tutorial-data>`_,
e.g., via

::

   mkdir snakemake-demo
   cd snakemake-demo
   wget https://github.com/snakemake/snakemake-tutorial-data/archive/v5.4.5.tar.gz
   tar --wildcards -xf v5.4.5.tar.gz --strip 1 "*/data"

Step 1
------

First, create an empty workflow in the current directory with:

::

   mkdir workflow
   touch workflow/Snakefile

Once a Snakefile is present, you can perform a dry run of Snakemake
with:

::

   snakemake -n

Since the Snakefile is empty, it will report that nothing has to be
done. In the next steps, we will gradually fill the Snakefile with an
example analysis workflow.
 
Step 2
------

The data folder in your working directory looks as follows:

::

   data
   ├── genome.fa
   ├── genome.fa.amb
   ├── genome.fa.ann
   ├── genome.fa.bwt
   ├── genome.fa.fai
   ├── genome.fa.pac
   ├── genome.fa.sa
   └── samples
       ├── A.fastq
       ├── B.fastq
       └── C.fastq

You will create a workflow that maps the sequencing samples in the
``data/samples`` folder to the reference genome ``data/genome.fa``.
Then, you will call genomic variants over the mapped samples, and create
an example plot.

First, create a rule called ``map_reads``, with input files

-  ``data/genome.fa``
-  ``data/samples/A.fastq``

and output file

-  ``results/mapped/A.bam``

To generate output from input, use the shell command

.. code:: python

       "bwa mem {input} | samtools view -Sb - > {output}"

Providing a shell command is not enough to run your workflow on an
unprepared system. For reproducibility, you also have to provide the
required software stack and define the desired version. This can be done
with the `Conda package manager <https://conda.io>`__, which is directly
integrated with Snakemake: add a directive
``conda: "envs/mapping.yaml"`` that points to a `Conda environment
definition <https://conda.io/docs/user-guide/tasks/manage-environments.html?highlight=environment#creating-an-environment-file-manually>`__,
with the following content

.. code:: yaml

       channels:
         - bioconda
         - conda-forge
       dependencies:
         - bwa =0.7.17
         - samtools =1.9

Upon execution, Snakemake will automatically create that environment,
and execute the shell command within.

Now, test your workflow by simulating the creation of the file
``results/mapped/A.bam`` via

::

   snakemake --use-conda -n results/mapped/A.bam

to perform a dry-run and

::

   snakemake --use-conda results/mapped/A.bam --cores 1

to perform the actual execution.
 
Step 3
------

Now, generalize the rule ``map_reads`` by replacing the concrete sample name
``A`` with a wildcard ``{sample}`` in input and output file the rule
``map_reads``. This way, Snakemake can apply the rule to map any of the three
available samples to the reference genome.

Test this by creating the file ``results/mapped/B.bam``.

Step 4
------

Next, create a rule ``sort_alignments`` that sorts the obtained ``.bam`` file by
genomic coordinate. The rule should have the input file

-  ``results/mapped/{sample}.bam``

and the output file

-  ``results/mapped/{sample}.sorted.bam``

and uses the shell command

::

   samtools sort -o {output} {input}

to perform the sorting. Moreover, use the same ``conda:`` directive as
for the previous rule.

Test your workflow with

::

   snakemake --use-conda -n results/mapped/A.sorted.bam

and

::

   snakemake --use-conda results/mapped/A.sorted.bam --cores 1

Step 5
------

Now, we aggregate over all samples to perform a joint calling of genomic
variants. First, we define a variable

.. code:: python

       samples = ["A", "B", "C"]

at the top of the ``Snakefile``. This serves as a definition of the
samples over which we would want to aggregate. In real life, you would
want to use an external sample sheet or a `config
file <https://snakemake.readthedocs.io/en/stable/tutorial/advanced.html#step-2-config-files>`__
for things like this.

For aggregation over many files, Snakemake provides the helper function
``expand`` (see `the
docs <https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-5-calling-genomic-variants>`__).
Create a rule ``call`` with input files

-  ``fa="data/genome.fa"``
-  ``bam=expand("results/mapped/{sample}.sorted.bam", sample=samples)``

output file

-  ``"results/calls/all.vcf"``

and shell command

::

   bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}

Further, define a new conda environment file with the following content:

.. code:: yaml

       channels:
         - bioconda
         - conda-forge
       dependencies:
         - bcftools =1.9

Step 6
------

Finally, we strive to calculate some exemplary statistics. This time, we
don’t use a shell command, but rather employ Snakemake’s ability to
integrate with scripting languages like R and Python, and Jupyter notebooks.

First, we create a rule ``plot_quals`` with input file

-  ``"results/calls/all.vcf"``

and output file

-  ``"results/plots/quals.svg"``.

Instead of a shell command, we use Snakemake's Jupyter notebook integration by specifying

.. code:: python

       notebook:
           "notebooks/plot-quals.py"

instead of using the ``shell`` directive as before.

Next, we have to define a conda environment for the rule, say
``workflow/envs/stats.yaml``, that provides the required Python packages to
execute the script:

.. code:: yaml

    channels:
      - bioconda
      - conda-forge
    dependencies:
      - pysam =0.17
      - altair =4.1
      - altair_saver =0.5
      - pandas =1.3
      - jupyter =1.0

Then, we let Snakemake generate a skeleton notebook for us with

.. code:: console

    snakemake --draft-notebook results/plots/quals.svg --cores 1 --use-conda

Snakemake will print instructions on how to open, edit and execute the notebook.

We open the notebook in the editor and add the following content

.. code:: python

    import pandas as pd
    import altair as alt
    from pysam import VariantFile

    quals = pd.DataFrame({"qual": [record.qual for record in VariantFile(snakemake.input[0])]})

    chart = alt.Chart(quals).mark_bar().encode(
        alt.X("qual", bin=True),
        alt.Y("count()")
    )

    chart.save(snakemake.output[0])

As you can see, instead of writing a command line parser for passing
parameters like input and output files, you have direct access to the
properties of the rule via a magic ``snakemake`` object, that Snakemake
automatically inserts into the notebook before executing the rule.

Make sure to test your workflow with

::

   snakemake --use-conda --force results/plots/quals.svg --cores 1

Here, the force ensures that the readily drafted notebook is re-executed even if you had already generated the output plot in the interactive mode.
 
Step 7
------

So far, we have always specified a target file at the command line when
invoking Snakemake. When no target file is specified, Snakemake tries to
execute the first rule in the ``Snakefile``. We can use this property to
define default target files.

At the top of your ``Snakefile`` define a rule ``all``, with input files

-  ``"results/calls/all.vcf"``
-  ``"results/plots/quals.svg"``

and neither a shell command nor output files. This rule simply serves as
an indicator of what shall be collected as results.

Step 8
------

As a last step, we strive to annotate our workflow with some additional
information.

Automatic reports
~~~~~~~~~~~~~~~~~

Snakemake can automatically create HTML reports with

::

   snakemake --report report.html

Such a report contains runtime statistics, a visualization of the
workflow topology, used software and data provenance information.

In addition, you can mark any output file generated in your workflow for
inclusion into the report. It will be encoded directly into the report,
such that it can be, e.g., emailed as a self-contained document. The
reader (e.g., a collaborator of yours) can at any time download the
enclosed results from the report for further use, e.g., in a manuscript
you write together. In this example, please mark the output file
``"results/plots/quals.svg"`` for inclusion by replacing it with
``report("results/plots/quals.svg", caption="report/calling.rst")`` and adding a
file ``report/calling.rst``, containing some description of the output
file. This description will be presented as caption in the resulting
report.

Threads
~~~~~~~

The first rule ``map_reads`` can in theory use multiple threads. You can make
Snakemake aware of this, such that the information can be used for
scheduling. Add a directive ``threads: 8`` to the rule and alter the
shell command to

::

   bwa mem -t {threads} {input} | samtools view -Sb - > {output}

This passes the threads defined in the rule as a command line argument
to the ``bwa`` process.

Temporary files
~~~~~~~~~~~~~~~

The output of the ``map_reads`` rule becomes superfluous once the sorted
version of the ``.bam`` file is generated by the rule ``sort``.
Snakemake can automatically delete the superfluous output once it is not
needed anymore. For this, mark the output as temporary by replacing
``"results/mapped/{sample}.bam"`` in the rule ``bwa`` with
``temp("results/mapped/{sample}.bam")``.

Solutions
---------

Only read this if you have a problem with one of the steps.

.. _step-2-1:

Step 2
~~~~~~

The rule should look like this:

.. code:: python

    rule map_reads:
        input:
            "data/genome.fa",
            "data/samples/A.fastq"
        output:
            "results/mapped/A.bam"
        conda:
            "envs/mapping.yaml"
        shell:
            "bwa mem {input} | samtools view -b - > {output}"

.. _step-3-1:

Step 3
~~~~~~

The rule should look like this:

.. code:: python

    rule map_reads:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "results/mapped/{sample}.bam"
        conda:
            "envs/mapping.yaml"
        shell:
            "bwa mem {input} | samtools view -b - > {output}"

.. _step-4-1:

Step 4
~~~~~~

The rule should look like this:

.. code:: python

    rule sort_alignments:
        input:
            "results/mapped/{sample}.bam"
        output:
            "results/mapped/{sample}.sorted.bam"
        conda:
            "envs/mapping.yaml"
        shell:
            "samtools sort -o {output} {input}"

.. _step-5-1:

Step 5
~~~~~~

The rule should look like this:

.. code:: python

    samples = ["A", "B", "C"]

    rule call_variants:
        input:
            fa="data/genome.fa",
            bam=expand("results/mapped/{sample}.sorted.bam", sample=SAMPLES)
        output:
            "results/calls/all.vcf"
        conda:
            "envs/calling.yaml"
        shell:
            "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}"

.. _step-6-1:

Step 6
~~~~~~

The rule should look like this:

.. code:: python

    rule plot_quals:
        input:
            "results/calls/all.vcf"
        output:
            "results/plots/quals.svg"
        conda:
            "envs/stats.yaml"
        notebook:
            "notebooks/plot-quals.py.ipynb"

.. _step-7-1:

Step 7
~~~~~~

The rule should look like this:

.. code:: python

    rule all:
        input:
            "results/calls/all.vcf",
            "results/plots/quals.svg"

It has to appear as first rule in the ``Snakefile``.

.. _step-8-1:

Step 8
~~~~~~

The complete workflow should look like this:

.. code:: python

    SAMPLES = ["A", "B", "C"]

    rule all:
        input:
            "results/calls/all.vcf",
            "results/plots/quals.svg"

    rule map_reads:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "results/mapped/{sample}.bam"
        conda:
            "envs/mapping.yaml"
        shell:
            "bwa mem {input} | samtools view -b - > {output}"


    rule sort_alignments:
        input:
            "results/mapped/{sample}.bam"
        output:
            "results/mapped/{sample}.sorted.bam"
        conda:
            "envs/mapping.yaml"
        shell:
            "samtools sort -o {output} {input}"


    rule call_variants:
        input:
            fa="data/genome.fa",
            bam=expand("results/mapped/{sample}.sorted.bam", sample=SAMPLES)
        output:
            "results/calls/all.vcf"
        conda:
            "envs/calling.yaml"
        shell:
            "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}"


    rule plot_quals:
        input:
            "results/calls/all.vcf"
        output:
            "results/plots/quals.svg"
        conda:
            "envs/stats.yaml"
        notebook:
            "notebooks/plot-quals.py.ipynb"
.. _tutorial-basics:

Basics: An example workflow
---------------------------

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Conda: https://conda.pydata.org
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: https://www.graphviz.org
.. _RestructuredText: https://docutils.sourceforge.io/docs/user/rst/quickstart.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: https://json.org
.. _YAML: https://yaml.org
.. _DRMAA: https://www.drmaa.org
.. _rpy2: https://rpy2.github.io
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: https://slides.com/johanneskoester/deck-1

Please make sure that you have **activated** the environment we created before, and that you have an open terminal in the working directory you have created.

**A Snakemake workflow is defined by specifying rules in a Snakefile**.
**Rules decompose the workflow into small steps** (for example, the application of a single tool) by specifying how to create sets of **output files** from sets of **input files**.
Snakemake automatically **determines the dependencies** between the rules by matching file names.

The Snakemake language extends the Python language, adding syntactic structures for rule definition and additional controls.
All added syntactic structures begin with a keyword followed by a code block that is either in the same line or indented and consisting of multiple lines.
The resulting syntax resembles that of original Python constructs.

In the following, we will introduce the Snakemake syntax by creating an example workflow.
The workflow comes from the domain of genome analysis.
It maps sequencing reads to a reference genome and calls variants on the mapped reads.
The tutorial does not require you to know what this is about.
Nevertheless, we provide some background in the following paragraph.

.. _tutorial-background:

Background
::::::::::

The genome of a living organism encodes its hereditary information.
It serves as a blueprint for proteins, which form living cells, carry information
and drive chemical reactions.
Differences between species, populations or individuals can be reflected by differences in the genome.
Certain variants can cause syndromes or predisposition for certain diseases, or cause cancerous growth in the case of tumour cells that have accumulated changes with respect to healthy cells.
This makes the genome a major target of biological and medical research.
Today, it is often analyzed with DNA sequencing, producing gigabytes of data from
a single biological sample (for example a biopsy of some tissue).
For technical reasons, DNA sequencing cuts the DNA of a sample into millions
of small pieces, called **reads**.
In order to recover the genome of the sample, one has to map these reads against
a known **reference genome** (for example, the human one obtained during the famous
`human genome project <https://en.wikipedia.org/wiki/Human_Genome_Project>`_).
This task is called **read mapping**.
Often, it is of interest where an individual genome is different from the species-wide consensus
represented with the reference genome.
Such differences are called **variants**. They are responsible for harmless individual
differences (like eye color), but can also cause diseases like cancer.
By investigating the differences between the mapped reads
and the reference sequence at a particular genome position, variants can be detected.
This is a statistical challenge, because they have
to be distinguished from artifacts generated by the sequencing process.

Step 1: Mapping reads
:::::::::::::::::::::

Our first Snakemake rule maps reads of a given sample to a given reference genome (see :ref:`tutorial-background`).
For this, we will use the tool bwa_, specifically the subcommand ``bwa mem``.
In the working directory, **create a new file** called ``Snakefile`` with an editor of your choice.
We propose to use the Atom_ editor, since it provides out-of-the-box syntax highlighting for Snakemake.
In the Snakefile, define the following rule:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/A.fastq"
        output:
            "mapped_reads/A.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

    A common error is to forget the comma between the input or output items.
    Since Python concatenates subsequent strings, this can lead to unexpected behavior.

A Snakemake rule has a name (here ``bwa_map``) and a number of directives, here ``input``, ``output`` and ``shell``.
The ``input`` and ``output`` directives are followed by lists of files that are expected to be used or created by the rule.
In the simplest case, these are just explicit Python strings.
The ``shell`` directive is followed by a Python string containing the shell command to execute.
In the shell command string, we can refer to elements of the rule via braces notation (similar to the Python format function).
Here, we refer to the output file by specifying ``{output}`` and to the input files by specifying ``{input}``.
Since the rule has multiple input files, Snakemake will concatenate them, separated by a whitespace.
In other words, Snakemake will replace ``{input}`` with ``data/genome.fa data/samples/A.fastq`` before executing the command.
The shell command invokes ``bwa mem`` with reference genome and reads, and pipes the output into ``samtools`` which creates a compressed `BAM <https://en.wikipedia.org/wiki/Binary_Alignment_Map>`_ file containing the alignments.
The output of ``samtools`` is redirected into the output file defined by the rule with ``>``.

.. sidebar:: Note

  It is best practice to have subsequent steps of a workflow in separate, unique, output folders. This keeps the working directory structured. Further, such unique prefixes allow Snakemake to quickly discard most rules in its search for rules that can provide the requested input. This accelerates the resolution of the rule dependencies in a workflow.

When a workflow is executed, Snakemake tries to generate given **target** files.
Target files can be specified via the command line.
By executing

.. code:: console

    $ snakemake -np mapped_reads/A.bam

in the working directory containing the Snakefile, we tell Snakemake to generate the target file ``mapped_reads/A.bam``.
Since we used the ``-n`` (or ``--dry-run``) flag, Snakemake will only show the execution plan instead of actually performing the steps.
The ``-p`` flag instructs Snakemake to also print the resulting shell command for illustration.
To generate the target files, **Snakemake applies the rules given in the Snakefile in a top-down way**.
The application of a rule to generate a set of output files is called **job**.
For each input file of a job, Snakemake again (i.e. recursively) determines rules that can be applied to generate it.
This yields a `directed acyclic graph (DAG) <https://en.wikipedia.org/wiki/Directed_acyclic_graph>`_ of jobs where the edges represent dependencies.
So far, we only have a single rule, and the DAG of jobs consists of a single node.
Nevertheless, we can **execute our workflow** with

.. code:: console

    $ snakemake --cores 1 mapped_reads/A.bam

Whenever executing a workflow, you need to specify the number of cores to use.
For this tutorial, we will use a single core for now. 
Later you will see how parallelization works.
Note that, after completion of above command, Snakemake will not try to create ``mapped_reads/A.bam`` again, because it is already present in the file system.
Snakemake **only re-runs jobs if one of the input files is newer than one of the output files or one of the input files will be updated by another job**.

Step 2: Generalizing the read mapping rule
::::::::::::::::::::::::::::::::::::::::::

Obviously, the rule will only work for a single sample with reads in the file ``data/samples/A.fastq``.
However, Snakemake allows **generalizing rules by using named wildcards**.
Simply replace the ``A`` in the second input file and in the output file with the wildcard ``{sample}``, leading to

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"

.. sidebar:: Note

  Note that if a rule has multiple output files, Snakemake requires them to all
  have exactly the same wildcards. Otherwise, it could happen that two jobs
  running the same rule in parallel want to write to the same file.

When Snakemake determines that this rule can be applied to generate a target file by replacing the wildcard ``{sample}`` in the output file with an appropriate value, it will propagate that value to all occurrences of ``{sample}`` in the input files and thereby determine the necessary input for the resulting job.
Note that you can have multiple wildcards in your file paths, however, to avoid conflicts with other jobs of the same rule, **all output files** of a rule have to **contain exactly the same wildcards**.

When executing

.. code:: console

    $ snakemake -np mapped_reads/B.bam

Snakemake will determine that the rule ``bwa_map`` can be applied to generate the target file by replacing the wildcard ``{sample}`` with the value ``B``.
In the output of the dry-run, you will see how the wildcard value is propagated to the input files and all filenames in the shell command.
You can also **specify multiple targets**, for example:

.. code:: console

    $ snakemake -np mapped_reads/A.bam mapped_reads/B.bam

Some Bash_ magic can make this particularly handy. For example, you can alternatively compose our multiple targets in a single pass via

.. code:: console

    $ snakemake -np mapped_reads/{A,B}.bam

Note that this is not a special Snakemake syntax.
Bash_ is just applying its `brace expansion <https://tldp.org/LDP/Bash-Beginners-Guide/html/sect_03_04.html>`_ to the set ``{A,B}``, creating the given path for each element and separating the resulting paths by a whitespace.

In both cases, you will see that Snakemake only proposes to create the output file ``mapped_reads/B.bam``.
This is because you already executed the workflow before (see the previous step) and no input file is newer than the output file ``mapped_reads/A.bam``.
You can update the file modification date of the input file
``data/samples/A.fastq`` via

.. code:: console

    $ touch data/samples/A.fastq

and see how Snakemake wants to re-run the job to create the file ``mapped_reads/A.bam`` by executing

.. code:: console

    $ snakemake -np mapped_reads/A.bam mapped_reads/B.bam


Step 3: Sorting read alignments
:::::::::::::::::::::::::::::::

For later steps, we need the read alignments in the BAM files to be sorted.
This can be achieved with the samtools_ ``sort`` command.
We add the following rule beneath the ``bwa_map`` rule:

.. code:: python

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

.. sidebar:: Note

  In the shell command above we split the string into two lines, which are however automatically concatenated into one by Python.
  This is a handy pattern to avoid too long shell command lines. When using this, make sure to have a trailing whitespace in each line but the last, 
  in order to avoid arguments to become not properly separated.

This rule will take the input file from the ``mapped_reads`` directory and store a sorted version in the ``sorted_reads`` directory.
Note that Snakemake **automatically creates missing directories** before jobs are executed.
For sorting, ``samtools`` requires a prefix specified with the flag ``-T``.
Here, we need the value of the wildcard ``sample``.
Snakemake allows to access wildcards in the shell command via the ``wildcards`` object that has an attribute with the value for each wildcard.

When issuing

.. code:: console

    $ snakemake -np sorted_reads/B.bam

you will see how Snakemake wants to run first the rule ``bwa_map`` and then the rule ``samtools_sort`` to create the desired target file:
as mentioned before, the dependencies are resolved automatically by matching file names.

Step 4: Indexing read alignments and visualizing the DAG of jobs
::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

.. sidebar:: Note

  Snakemake uses the `Python format mini language <https://docs.python.org/3/library/string.html#formatexamples>`_ to format shell commands.
  Sometimes you have to use braces (``{}``) for something else in a shell command.
  In that case, you have to escape them by doubling, for example when relying on the bash brace expansion we mentioned above:
  ``ls {{A,B}}.txt``.

Next, we need to use samtools_ again to index the sorted read alignments so that we can quickly access reads by the genomic location they were mapped to.
This can be done with the following rule:

.. code:: python

    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"

Having three steps already, it is a good time to take a closer look at the resulting directed acyclic graph (DAG) of jobs.
By executing

.. code:: console

    $ snakemake --dag sorted_reads/{A,B}.bam.bai | dot -Tsvg > dag.svg


.. sidebar:: Note

  If you went with: `Run tutorial for free in the cloud via Gitpod`_, you can easily view the resulting ``dag.svg`` by right-clicking on the file in the explorer panel on the left and selecting ``Open With -> Preview``.


we create a **visualization of the DAG** using the ``dot`` command provided by Graphviz_.
For the given target files, Snakemake specifies the DAG in the dot language and pipes it into the ``dot`` command, which renders the definition into `SVG format <https://en.wikipedia.org/wiki/Scalable_Vector_Graphics>`_.
The rendered DAG is piped into the file ``dag.svg`` and will look similar to this:

.. image:: workflow/dag_index.png
   :align: center

The DAG contains a node for each job with the edges connecting them representing the dependencies.
The frames of jobs that don't need to be run (because their output is up-to-date) are dashed.
For rules with wildcards, the value of the wildcard for the particular job is displayed in the job node.

Exercise
........

* Run parts of the workflow using different targets. Recreate the DAG and see how different rules' frames become dashed because their output is present and up-to-date.

Step 5: Calling genomic variants
::::::::::::::::::::::::::::::::

The next step in our workflow will aggregate the mapped reads from all samples and jointly call genomic variants on them (see :ref:`tutorial-background`).
For the variant calling, we will combine the two utilities samtools_ and bcftools_.
Snakemake provides a **helper function for collecting input files** that helps us to describe the aggregation in this step.
With

.. code:: python

    expand("sorted_reads/{sample}.bam", sample=SAMPLES)

we obtain a list of files where the given pattern ``"sorted_reads/{sample}.bam"`` was formatted with the values in a given list of samples ``SAMPLES``, i.e.

.. code:: python

    ["sorted_reads/A.bam", "sorted_reads/B.bam"]

The function is particularly useful when the pattern contains multiple wildcards.
For example,

.. code:: python

    expand("sorted_reads/{sample}.{replicate}.bam", sample=SAMPLES, replicate=[0, 1])

would create the product of all elements of ``SAMPLES`` and the list ``[0, 1]``, yielding

.. code:: python

    ["sorted_reads/A.0.bam", "sorted_reads/A.1.bam", "sorted_reads/B.0.bam", "sorted_reads/B.1.bam"]

Here, we use only the simple case of ``expand``.
We first let Snakemake know which samples we want to consider.
Remember that Snakemake works backwards from requested output, and not from available input.
Thus, it does not automatically infer all possible output from, for example, the fastq files in the data folder.
Also remember that Snakefiles are in principle Python code enhanced by some declarative statements to define workflows.
Hence, we can define the list of samples ad-hoc in plain Python at the top of the Snakefile:

.. code:: python

    SAMPLES = ["A", "B"]


.. sidebar:: Note

  If you name input or output files like above, their order won't be preserved when referring to them as ``{input}``.
  Further, note that named and unnamed (i.e., positional) input and output files can be combined, but the positional ones must come first, equivalent to Python functions with keyword arguments.

Later, we will learn about more sophisticated ways like **config files**.
But for now, this is enough so that we can add the following rule to our Snakefile:

.. code:: python

    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

With multiple input or output files, it is sometimes handy to refer to them separately in the shell command.
This can be done by **specifying names for input or output files**, for example with ``fa=...``.
The files can then be referred to in the shell command by name, for example with ``{input.fa}``.
For **long shell commands** like this one, it is advisable to **split the string over multiple indented lines**.
Python will automatically merge it into one.
Further, you will notice that the **input or output file lists can contain arbitrary Python statements**, as long as it returns a string, or a list of strings.
Here, we invoke our ``expand`` function to aggregate over the aligned reads of all samples.


Exercise
........

* obtain the updated DAG of jobs for the target file ``calls/all.vcf``, it should look like this:

.. image:: workflow/dag_call.png
   :align: center


.. _tutorial-script:

Step 6: Using custom scripts
::::::::::::::::::::::::::::

Usually, a workflow not only consists of invoking various tools, but also contains custom code to for example calculate summary statistics or create plots.
While Snakemake also allows you to directly :ref:`write Python code inside a rule <.. _snakefiles-rules>`, it is usually reasonable to move such logic into separate scripts.
For this purpose, Snakemake offers the ``script`` directive.
Add the following rule to your Snakefile:

.. code:: python

    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        script:
            "scripts/plot-quals.py"


.. sidebar:: Note

  ``snakemake.input`` and ``snakemake.output`` always contain a list of file names, even if the lists each contain only one file name.
  Therefore, to refer to a particular file name, you have to index into that list.
  ``snakemake.output[0]`` will give you the first element of the output file name list, something that always has to be there.

With this rule, we will eventually generate a histogram of the quality scores that have been assigned to the variant calls in the file ``calls/all.vcf``.
The actual Python code to generate the plot is hidden in the script ``scripts/plot-quals.py``.
Script paths are always relative to the referring Snakefile.
In the script, all properties of the rule like ``input``, ``output``, ``wildcards``, etc. are available as attributes of a global ``snakemake`` object.
Create the file ``scripts/plot-quals.py``, with the following content:

.. code:: python

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pysam import VariantFile

    quals = [record.qual for record in VariantFile(snakemake.input[0])]
    plt.hist(quals)

    plt.savefig(snakemake.output[0])


.. sidebar:: Note

  It is best practice to use the script directive whenever an inline code block would have
  more than a few lines of code.

Although there are other strategies to invoke separate scripts from your workflow
(for example, invoking them via shell commands), the benefit of this is obvious:
the script logic is separated from the workflow logic (and can even be shared between workflows),
but **boilerplate code like the parsing of command line arguments is unnecessary**.

Apart from Python scripts, it is also possible to use R scripts. In R scripts,
an S4 object named ``snakemake`` analogous to the Python case above is available and
allows access to input and output files and other parameters. Here, the syntax
follows that of S4 classes with attributes that are R lists, for example we can access
the first input file with ``snakemake@input[[1]]`` (note that the first file does
not have index 0 here, because R starts counting from 1). Named input and output
files can be accessed in the same way, by just providing the name instead of an
index, for example ``snakemake@input[["myfile"]]``.

For details and examples, see the :ref:`snakefiles-external_scripts` section in the Documentation.


Step 7: Adding a target rule
::::::::::::::::::::::::::::

So far, we always executed the workflow by specifying a target file at the command line.
Apart from filenames, Snakemake **also accepts rule names as targets** if the requested rule does not have wildcards.
Hence, it is possible to write target rules collecting particular subsets of the desired results or all results.
Moreover, if no target is given at the command line, Snakemake will define the **first rule** of the Snakefile as the target.
Hence, it is best practice to have a rule ``all`` at the top of the workflow which has all typically desired target files as input files.

Here, this means that we add a rule

.. code:: python

    rule all:
        input:
            "plots/quals.svg"

to the top of our workflow.
When executing Snakemake with

.. code:: console

    $ snakemake -n

.. sidebar:: Note

   In case you have mutliple reasonable sets of target files,
   you can add multiple target rules at the top of the Snakefile. While
   Snakemake will execute the first per default, you can target any of them via
   the command line (for example, ``snakemake -n mytarget``).

the execution plan for creating the file ``plots/quals.svg``, which contains and summarizes all our results, will be shown.
Note that, apart from Snakemake considering the first rule of the workflow as the default target, **the order of rules in the Snakefile is arbitrary and does not influence the DAG of jobs**.

Exercise
........

* Create the DAG of jobs for the complete workflow.
* Execute the complete workflow and have a look at the resulting ``plots/quals.svg``.
* Snakemake provides handy flags for forcing re-execution of parts of the workflow. Have a look at the command line help with ``snakemake --help`` and search for the flag ``--forcerun``. Then, use this flag to re-execute the rule ``samtools_sort`` and see what happens.
* With ``--reason`` it is possible to display the execution reason for each job. Try this flag together with a dry-run and the ``--forcerun`` flag to understand the decisions of Snakemake.

Summary
:::::::

In total, the resulting workflow looks like this:

.. code:: console

    SAMPLES = ["A", "B"]


    rule all:
        input:
            "plots/quals.svg"


    rule bwa_map:
        input:
            "data/genome.fa",
            "data/samples/{sample}.fastq"
        output:
            "mapped_reads/{sample}.bam"
        shell:
            "bwa mem {input} | samtools view -Sb - > {output}"


    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"


    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        shell:
            "samtools index {input}"


    rule bcftools_call:
        input:
            fa="data/genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"


    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        script:
            "scripts/plot-quals.py"
.. _tutorial:

==================
Snakemake Tutorial
==================

.. _Snakemake: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _slides: https://slides.com/johanneskoester/snakemake-tutorial
.. _Conda: https://conda.io
.. _Singularity: https://www.sylabs.io

This tutorial introduces the text-based workflow system Snakemake_.
Snakemake follows the `GNU Make`_ paradigm: workflows are defined in terms of rules that define how to create output files from input files.
Dependencies between the rules are determined automatically, creating a DAG (directed acyclic graph) of jobs that can be automatically parallelized.

Snakemake sets itself apart from other text-based workflow systems in the following way.
Hooking into the Python interpreter, Snakemake offers a definition language that is an extension of Python_ with syntax to define rules and workflow specific properties.
This allows to combine the flexibility of a plain scripting language with a pythonic workflow definition.
The Python language is known to be concise yet readable and can appear almost like pseudo-code.
The syntactic extensions provided by Snakemake maintain this property for the definition of the workflow.
Further, Snakemake's scheduling algorithm can be constrained by priorities, provided cores and customizable resources and it provides a generic support for distributed computing (e.g., cluster or batch systems).
Hence, a Snakemake workflow scales without modification from single core workstations and multi-core servers to cluster or batch systems.
Finally, Snakemake integrates with the package manager Conda_ and the container engine Singularity_ such that defining the software stack becomes part of the workflow itself.

The examples presented in this tutorial come from Bioinformatics.
However, Snakemake is a general-purpose workflow management system for any discipline.
We ensured that no bioinformatics knowledge is needed to understand the tutorial.

Also have a look at the corresponding slides_.


.. toctree::
   :maxdepth: 2

   setup
   basics
   advanced
   additional_features
.. tutorial-additional_features:

Additional features
-------------------

.. _Snakemake: https://snakemake.readthedocs.io
.. _Snakemake homepage: https://snakemake.readthedocs.io
.. _GNU Make: https://www.gnu.org/software/make
.. _Python: https://www.python.org
.. _BWA: http://bio-bwa.sourceforge.net
.. _SAMtools: https://www.htslib.org
.. _BCFtools: https://www.htslib.org
.. _Pandas: https://pandas.pydata.org
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Conda: https://conda.pydata.org
.. _Bash: https://www.tldp.org/LDP/Bash-Beginners-Guide/html
.. _Atom: https://atom.io
.. _Anaconda: https://anaconda.org
.. _Graphviz: https://www.graphviz.org
.. _RestructuredText: https://docutils.sourceforge.io/docs/user/rst/quickstart.html
.. _data URI: https://developer.mozilla.org/en-US/docs/Web/HTTP/data_URIs
.. _JSON: https://json.org
.. _YAML: https://yaml.org
.. _DRMAA: https://www.drmaa.org
.. _rpy2: https://rpy2.github.io
.. _R: https://www.r-project.org
.. _Rscript: https://stat.ethz.ch/R-manual/R-devel/library/utils/html/Rscript.html
.. _PyYAML: https://pyyaml.org
.. _Docutils: https://docutils.sourceforge.io
.. _Bioconda: https://bioconda.github.io
.. _Vagrant: https://www.vagrantup.com
.. _Vagrant Documentation: https://docs.vagrantup.com
.. _Blogpost: https://blog.osteel.me/posts/2015/01/25/how-to-use-vagrant-on-windows.html
.. _slides: https://slides.com/johanneskoester/deck-1

In the following, we introduce some features that are beyond the scope of above example workflow.
For details and even more features, see :ref:`user_manual-writing_snakefiles`, :ref:`project_info-faq` and the command line help (``snakemake --help``).


Benchmarking
::::::::::::

With the ``benchmark`` directive, Snakemake can be instructed to **measure the wall clock time of a job**.
We activate benchmarking for the rule ``bwa_map``:

.. code:: python

    rule bwa_map:
        input:
            "data/genome.fa",
            lambda wildcards: config["samples"][wildcards.sample]
        output:
            temp("mapped_reads/{sample}.bam")
        params:
            rg="@RG\tID:{sample}\tSM:{sample}"
        log:
            "logs/bwa_mem/{sample}.log"
        benchmark:
            "benchmarks/{sample}.bwa.benchmark.txt"
        threads: 8
        shell:
            "(bwa mem -R '{params.rg}' -t {threads} {input} | "
            "samtools view -Sb - > {output}) 2> {log}"

The ``benchmark`` directive takes a string that points to the file where benchmarking results shall be stored.
Similar to output files, the path can contain wildcards (it must be the same wildcards as in the output files).
When a job derived from the rule is executed, Snakemake will measure the wall clock time and memory usage (in MiB) and store it in the file in tab-delimited format.
It is possible to repeat a benchmark multiple times in order to get a sense for the variability of the measurements.
This can be done by annotating the benchmark file, e.g., with ``repeat("benchmarks/{sample}.bwa.benchmark.txt", 3)`` Snakemake can be told to run the job three times.
The repeated measurements occur as subsequent lines in the tab-delimited benchmark file.

Modularization
::::::::::::::

In order to re-use building blocks or simply to structure large workflows, it is sometimes reasonable to **split a workflow into modules**.
For this, Snakemake provides the ``include`` directive to include another Snakefile into the current one, e.g.:

.. code:: python

    include: "path/to/other.snakefile"

Alternatively, Snakemake allows to **define sub-workflows**.
A sub-workflow refers to a working directory with a complete Snakemake workflow.
Output files of that sub-workflow can be used in the current Snakefile.
When executing, Snakemake ensures that the output files of the sub-workflow are up-to-date before executing the current workflow.
This mechanism is particularly useful when you want to extend a previous analysis without modifying it.
For details about sub-workflows, see the :ref:`documentation <snakefiles-sub_workflows>`.


Exercise
........

* Put the read mapping related rules into a separate Snakefile and use the ``include`` directive to make them available in our example workflow again.


.. _tutorial-conda:

Automatic deployment of software dependencies
:::::::::::::::::::::::::::::::::::::::::::::

In order to get a fully reproducible data analysis, it is not sufficient to
be able to execute each step and document all used parameters.
The used software tools and libraries have to be documented as well.
In this tutorial, you have already seen how Conda_ can be used to specify an
isolated software environment for a whole workflow. With Snakemake, you can
go one step further and specify Conda environments per rule.
This way, you can even make use of conflicting software versions (e.g. combine
Python 2 with Python 3).

In our example, instead of using an external environment we can specify
environments per rule, e.g.:

.. code:: python

  rule samtools_index:
    input:
        "sorted_reads/{sample}.bam"
    output:
        "sorted_reads/{sample}.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

with ``envs/samtools.yaml`` defined as

.. code:: yaml

  channels:
    - bioconda
    - conda-forge
  dependencies:
    - samtools =1.9

.. sidebar:: Note

  The conda directive does not work in combination with ``run`` blocks, because
  they have to share their Python environment with the surrounding snakefile.

When Snakemake is executed with

.. code:: console

  snakemake --use-conda --cores 1

it will automatically create required environments and
activate them before a job is executed.
It is best practice to specify at least the `major and minor version <https://semver.org/>`_ of any packages
in the environment definition. Specifying environments per rule in this way has two
advantages.
First, the workflow definition also documents all used software versions.
Second, a workflow can be re-executed (without admin rights)
on a vanilla system, without installing any
prerequisites apart from Snakemake and Miniconda_.


Tool wrappers
:::::::::::::

In order to simplify the utilization of popular tools, Snakemake provides a
repository of so-called wrappers
(the `Snakemake wrapper repository <https://snakemake-wrappers.readthedocs.io>`_).
A wrapper is a short script that wraps (typically)
a command line application and makes it directly addressable from within Snakemake.
For this, Snakemake provides the ``wrapper`` directive that can be used instead of
``shell``, ``script``, or ``run``.
For example, the rule ``bwa_map`` could alternatively look like this:

.. code:: python

  rule bwa_mem:
    input:
        ref="data/genome.fa",
        sample=lambda wildcards: config["samples"][wildcards.sample]
    output:
        temp("mapped_reads/{sample}.bam")
    log:
        "logs/bwa_mem/{sample}.log"
    params:
        "-R '@RG\tID:{sample}\tSM:{sample}'"
    threads: 8
    wrapper:
        "0.15.3/bio/bwa/mem"

.. sidebar:: Note

  Updates to the Snakemake wrapper repository are automatically tested via
  `continuous integration <https://en.wikipedia.org/wiki/Continuous_integration>`_.

The wrapper directive expects a (partial) URL that points to a wrapper in the repository.
These can be looked up in the corresponding `database <https://snakemake-wrappers.readthedocs.io>`_.
The first part of the URL is a Git version tag. Upon invocation, Snakemake
will automatically download the requested version of the wrapper.
Furthermore, in combination with ``--use-conda`` (see :ref:`tutorial-conda`),
the required software will be automatically deployed before execution.

Cluster execution
:::::::::::::::::

By default, Snakemake executes jobs on the local machine it is invoked on.
Alternatively, it can execute jobs in **distributed environments, e.g., compute clusters or batch systems**.
If the nodes share a common file system, Snakemake supports three alternative execution modes.

In cluster environments, compute jobs are usually submitted as shell scripts via commands like ``qsub``.
Snakemake provides a **generic mode** to execute on such clusters.
By invoking Snakemake with

.. code:: console

    $ snakemake --cluster qsub --jobs 100

each job will be compiled into a shell script that is submitted with the given command (here ``qsub``).
The ``--jobs`` flag limits the number of concurrently submitted jobs to 100.
This basic mode assumes that the submission command returns immediately after submitting the job.
Some clusters allow to run the submission command in **synchronous mode**, such that it waits until the job has been executed.
In such cases, we can invoke e.g.

.. code:: console

    $ snakemake --cluster-sync "qsub -sync yes" --jobs 100

The specified submission command can also be **decorated with additional parameters taken from the submitted job**.
For example, the number of used threads can be accessed in braces similarly to the formatting of shell commands, e.g.

.. code:: console

    $ snakemake --cluster "qsub -pe threaded {threads}" --jobs 100

Alternatively, Snakemake can use the Distributed Resource Management Application API (DRMAA_).
This API provides a common interface to control various resource management systems.
The **DRMAA support** can be activated by invoking Snakemake as follows:

.. code:: console

    $ snakemake --drmaa --jobs 100

If available, **DRMAA is preferable over the generic cluster modes** because it provides better control and error handling.
To support additional cluster specific parametrization, a Snakefile can be complemented by a :ref:`snakefiles-cluster_configuration` file.

Using --cluster-status
::::::::::::::::::::::

Sometimes you need specific detection to determine if a cluster job completed successfully, failed or is still running.
Error detection with ``--cluster`` can be improved for edge cases such as timeouts and jobs exceeding memory that are silently terminated by
the queueing system.
This can be achieved with the ``--cluster-status`` option. The value of this option should be a executable script which takes a job id as the first argument and prints to stdout only one of [running|success|failed]. Importantly, the job id snakemake passes on is captured from the stdout of the cluster submit tool. This string will often include more than the job id, but snakemake does not modify this string and will pass this string to the status script unchanged. In the situation where snakemake has received more than the job id these are 3 potential solutions to consider: parse the string received by the script and extract the job id within the script, wrap the submission tool to intercept its stdout and return just the job code, or ideally, the cluster may offer an option to only return the job id upon submission and you can instruct snakemake to use that option. For sge this would look like  ``snakemake --cluster "qsub -terse"``.

The following (simplified) script detects the job status on a given SLURM cluster (>= 14.03.0rc1 is required for ``--parsable``).

.. code:: python

    #!/usr/bin/env python
    import subprocess
    import sys

    jobid = sys.argv[1]

    output = str(subprocess.check_output("sacct -j %s --format State --noheader | head -1 | awk '{print $1}'" % jobid, shell=True).strip())

    running_status=["PENDING", "CONFIGURING", "COMPLETING", "RUNNING", "SUSPENDED"]
    if "COMPLETED" in output:
      print("success")
    elif any(r in output for r in running_status):
      print("running")
    else:
      print("failed")

To use this script call snakemake similar to below, where ``status.py`` is the script above.

.. code:: console

    $ snakemake all --jobs 100 --cluster "sbatch --cpus-per-task=1 --parsable" --cluster-status ./status.py


Constraining wildcards
::::::::::::::::::::::

Snakemake uses regular expressions to match output files to input files and determine dependencies between the jobs.
Sometimes it is useful to constrain the values a wildcard can have.
This can be achieved by adding a regular expression that describes the set of allowed wildcard values.
For example, the wildcard ``sample`` in the output file ``"sorted_reads/{sample}.bam"`` can be constrained to only allow alphanumeric sample names as ``"sorted_reads/{sample,[A-Za-z0-9]+}.bam"``.
Constraints may be defined per rule or globally using the ``wildcard_constraints`` keyword, as demonstrated in :ref:`snakefiles-wildcards`.
This mechanism helps to solve two kinds of ambiguity.

* It can help to avoid ambiguous rules, i.e. two or more rules that can be applied to generate the same output file. Other ways of handling ambiguous rules are described in the Section :ref:`snakefiles-ambiguous-rules`.
* It can help to guide the regular expression based matching so that wildcards are assigned to the right parts of a file name. Consider the output file ``{sample}.{group}.txt`` and assume that the target file is ``A.1.normal.txt``. It is not clear whether ``dataset="A.1"`` and ``group="normal"`` or ``dataset="A"`` and ``group="1.normal"`` is the right assignment. Here, constraining the dataset wildcard by ``{sample,[A-Z]+}.{group}`` solves the problem.

When dealing with ambiguous rules, it is best practice to first try to solve the ambiguity by using a proper file structure, for example, by separating the output files of different steps in different directories.

.. _snakefiles-foreign-wms:

===============================================
Integrating foreign workflow management systems
===============================================

Snakemake 6.2 and later allows to hand over execution steps to other workflow management systems.
By this, it is possible to make use of workflows written for other systems, while performing any further pre- or postprocessing within Snakemake.
Such a handover is indicated with the ``handover`` directive.
Consider the following example:

.. code-block:: python

    rule chipseq_pipeline:
        input:
            input="design.csv",
            fasta="data/genome.fasta",
            gtf="data/genome.gtf",
        output:
            "multiqc/broadPeaks/multiqc_report.html",
        params:
            pipeline="nf-core/chipseq",
            revision="1.2.1",
            profile=["conda"],
        handover: True
        wrapper:
            "0.74.0/utils/nextflow"

Here, the workflow is executed as usual until this rule is reached.
Then, Snakemake passes all resources to the nextflow workflow management system, which generates certain files.
The rule is executed as a :ref:`local rule <snakefiles-local-rule>`, meaning that it would not be submitted to a cluster or cloud system by Snakemake.
Instead, the invoked other workflow management system is responsible for that.
E.g., in case of `Nextflow <https://nextflow.io>`_, submission behavior can be configured via a ``nextflow.conf`` file or environment variables.
After the step is done, Snakemake continues execution with the output files produced by the foreign workflow... _snakefiles-utils:

=====
Utils
=====

The module ``snakemake.utils`` provides a collection of helper functions for common tasks in Snakemake workflows. Details can be found in :ref:`utils-api`.
.. _user_manual-writing_snakefiles:

=================
Writing Workflows
=================

In Snakemake, workflows are specified as Snakefiles.
Inspired by GNU Make, a Snakefile contains rules that denote how to create output files from input files.
Dependencies between rules are handled implicitly, by matching filenames of input files against output files.
Thereby wildcards can be used to write general rules.

.. _snakefiles-grammar:

-------
Grammar
-------

The Snakefile syntax obeys the following grammar, given in extended Backus-Naur form (EBNF)

.. code-block:: text

    snakemake    = statement | rule | include | workdir | module | configfile | container
    rule         = "rule" (identifier | "") ":" ruleparams
    include      = "include:" stringliteral
    workdir      = "workdir:" stringliteral
    module       = "module" identifier ":" moduleparams
    configfile   = "configfile" ":" stringliteral
    userule      = "use" "rule" (identifier | "*") "from" identifier ["as" identifier] ["with" ":" norunparams]
    ni           = NEWLINE INDENT
    norunparams  = [ni input] [ni output] [ni params] [ni message] [ni threads] [ni resources] [ni log] [ni conda] [ni container] [ni benchmark] [ni cache]
    ruleparams   = norunparams [ni (run | shell | script | notebook)] NEWLINE snakemake
    input        = "input" ":" parameter_list
    output       = "output" ":" parameter_list
    params       = "params" ":" parameter_list
    log          = "log" ":" parameter_list
    benchmark    = "benchmark" ":" statement
    cache        = "cache" ":" bool
    message      = "message" ":" stringliteral
    threads      = "threads" ":" integer
    resources    = "resources" ":" parameter_list
    version      = "version" ":" statement
    conda        = "conda" ":" stringliteral
    container    = "container" ":" stringliteral
    run          = "run" ":" ni statement
    shell        = "shell" ":" stringliteral
    script       = "script" ":" stringliteral
    notebook     = "notebook" ":" stringliteral
    moduleparams = [ni snakefile] [ni metawrapper] [ni config] [ni skipval]
    snakefile    = "snakefile" ":" stringliteral
    metawrapper  = "meta_wrapper" ":" stringliteral
    config       = "config" ":" stringliteral
    skipval      = "skip_validation" ":" stringliteral
    

while all not defined non-terminals map to their Python equivalents.

.. _snakefiles-depend_version:

Depend on a Minimum Snakemake Version
-------------------------------------

From Snakemake 3.2 on, if your workflow depends on a minimum Snakemake version, you can easily ensure that at least this version is installed via

.. code-block:: python

    from snakemake.utils import min_version

    min_version("3.2")

given that your minimum required version of Snakemake is 3.2. The statement will raise a WorkflowError (and therefore abort the workflow execution) if the version is not met... _snakefiles_configuration:

=============
Configuration
=============

Snakemake allows you to use configuration files for making your workflows more flexible and also for abstracting away direct dependencies to a fixed HPC cluster scheduler.


.. _snakefiles_standard_configuration:

----------------------
Standard Configuration
----------------------

Snakemake directly supports the configuration of your workflow.
A configuration is provided as a JSON or YAML file and can be loaded with:

.. code-block:: python

    configfile: "path/to/config.yaml"

The config file can be used to define a dictionary of configuration parameters and their values.
In the workflow, the configuration is accessible via the global variable `config`, e.g.

.. code-block:: python

    rule all:
        input:
            expand("{sample}.{param}.output.pdf", sample=config["samples"], param=config["yourparam"])

If the `configfile` statement is not used, the config variable provides an empty array.
In addition to the `configfile` statement, config values can be overwritten via the command line or the :ref:`api_reference_snakemake`, e.g.:

.. code-block:: console

    $ snakemake --config yourparam=1.5

Further, you can manually alter the config dictionary using any Python code **outside** of your rules. Changes made from within a rule won't be seen from other rules.
Finally, you can use the ``--configfile`` command line argument to overwrite values from the `configfile` statement.
Note that any values parsed into the ``config`` dictionary with any of above mechanisms are merged, i.e., all keys defined via a ``configfile``
statement, or the ``--configfile`` and ``--config`` command line arguments will end up in the final `config` dictionary, but if two methods define the same key, command line
overwrites the ``configfile`` statement.

For adding config placeholders into a shell command, Python string formatting syntax requires you to leave out the quotes around the key name, like so:

.. code-block:: python

    shell:
        "mycommand {config[foo]} ..."

.. _snakefiles_tabular_configuration

---------------------
Tabular configuration
---------------------

It is usually advisable to complement YAML based configuration (see above) by a sheet based approach for meta-data that is of tabular form. For example, such
a sheet can contain per-sample information.
With the `Pandas library <https://pandas.pydata.org/>`_ such data can be read and used with minimal overhead, e.g.,

.. code-block:: python

    import pandas as pd

    samples = pd.read_table("samples.tsv").set_index("samples", drop=False)

reads in a table ``samples.tsv`` in TSV format and makes every record accessible by the sample name.
For details, see the `Pandas documentation <https://pandas.pydata.org/pandas-docs/stable/generated/pandas.read_table.html?highlight=read_table#pandas-read-table>`_.
A fully working real-world example containing both types of configuration can be found `here <https://github.com/snakemake-workflows/rna-seq-star-deseq2>`_.

---------------------
Environment variables
---------------------

Sometimes, it is not desirable to put configuration information into text files.
For example, this holds for secrets like access tokens or passwords.
Here, `environment variables <https://en.wikipedia.org/wiki/Environment_variable>`_ are the method of choice.
Snakemake allows to assert the existence of environment variables by adding a statement like:

.. code-block:: python

    envvars:
        "SOME_VARIABLE",
        "SOME_OTHER_VARIABLE"

When executing, Snakemake will fail with a reasonable error message if the variables ``SOME_VARIABLE`` and ``SOME_OTHER_VARIABLE`` are undefined.
Otherwise, it will take care of passing them to cluster and cloud environments. However, note that this does **not** mean that Snakemake makes them available e.g. in the jobs shell command.
Instead, for data provenance and reproducibility reasons, you are required to pass them explicitly to your job via the params directive, e.g. like this:

.. code-block:: python

    envvars:
        "SOME_VARIABLE"

    rule do_something:
        output:
             "test.txt"
        params:
            x=os.environ["SOME_VARIABLE"]
        shell:
            "echo {params.x} > {output}"


.. _snakefiles_config_validation:

----------
Validation
----------

With Snakemake 5.1, it is possible to validate both types of configuration via `JSON schemas <https://json-schema.org>`_.
The function ``snakemake.utils.validate`` takes a loaded configuration (a config dictionary or a Pandas data frame) and validates it with a given JSON schema.
Thereby, the schema can be provided in JSON or YAML format. Also, by using the defaults property it is possible to populate entries with default values. See `jsonschema FAQ on setting default values <https://python-jsonschema.readthedocs.io/en/latest/faq/>`_ for details.
In case of the data frame, the schema should model the record that is expected in each row of the data frame.
In the following example,

.. code-block:: python

  import pandas as pd
  from snakemake.utils import validate

  configfile: "config.yaml"
  validate(config, "config.schema.yaml")

  samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
  validate(samples, "samples.schema.yaml")


  rule all:
      input:
          expand("test.{sample}.txt", sample=samples.index)


  rule a:
      output:
          "test.{sample}.txt"
      shell:
          "touch {output}"

the schema for validating the samples data frame looks like this:

.. code-block:: yaml

  $schema: "https://json-schema.org/draft-06/schema#"
  description: an entry in the sample sheet
  properties:
    sample:
      type: string
      description: sample name/identifier
    condition:
      type: string
      description: sample condition that will be compared during differential expression analysis (e.g. a treatment, a tissue time, a disease)
    case:
      type: boolean
      default: true
      description: boolean that indicates if sample is case or control

  required:
    - sample
    - condition

Here, in case the case column is missing, the validate function will
populate it with True for all entries.

.. _snakefiles-peps:

-------------------------------------------
Configuring scientific experiments via PEPs
-------------------------------------------

Often scientific experiments consist of a set of samples (with optional subsamples), for which raw data and metainformation is known.
Instead of writing custom sample sheets as shown above, Snakemake allows to use `portable encapsulated project (PEP) <http://pep.databio.org>`_ definitions to configure a workflow.
This is done via a special directive `pepfile`, that can optionally complemented by a schema for validation (which is recommended for production workflows):

.. code-block:: python

    pepfile: "pep/config.yaml"
    pepschema: "schemas/pep.yaml"

    rule all:
        input:
            expand("{sample}.txt", sample=pep.sample_table["sample_name"])

    rule a:
        output:
            "{sample}.txt"
        shell:
            "touch {output}"

Using the ``pepfile`` directive leads to parsing of the provided PEP with `peppy <http://peppy.databio.org>`_.
The resulting project object is made globally available under the name ``pep``.
Here, we use it to aggregate over the set of sample names that is defined in the corresponding PEP.

**Importantly**, note that PEPs are meant to contain sample metadata and any global information about a project or experiment. 
They should **not** be used to encode workflow specific configuration options.
For those, one should always complement the pepfile with an ordinary :ref:`config file <snakefiles_standard_configuration>`.
The rationale is that PEPs should be portable between different data analysis workflows (that could be applied to the same data) and even between workflow management systems.
In other words, a PEP should describe everything needed about the data, while a workflow and its configuration should describe everything needed about the analysis that is applied to it.

^^^^^^^^^^^^^^^
Validating PEPs
^^^^^^^^^^^^^^^

Using the ``pepschema`` directive leads to an automatic parsing of the provided schema *and* PEP validation with the PEP validation tool -- `eido <http://eido.databio.org>`_. Eido schemas extend `JSON Schema <https://json-schema.org>`_ vocabulary to accommodate the powerful PEP features. Follow the `How to write a PEP schema <http://eido.databio.org/en/latest/writing-a-schema>`_ guide to learn more.

.. _snakefiles-cluster_configuration:

----------------------------------
Cluster Configuration (deprecated)
----------------------------------

While still being possible, **cluster configuration has been deprecated** by the introduction of :ref:`profiles`.

Snakemake supports a separate configuration file for execution on a cluster.
A cluster config file allows you to specify cluster submission parameters outside the Snakefile.
The cluster config is a JSON- or YAML-formatted file that contains objects that match names of rules in the Snakefile.
The parameters in the cluster config are then accessed by the ``cluster.*`` wildcard when you are submitting jobs.
Note that a workflow shall never depend on a cluster configuration, because this would limit its portability.
Therefore, it is also not intended to access the cluster configuration from **within** the workflow.

For example, say that you have the following Snakefile:

.. code-block:: python

    rule all:
        input: "input1.txt", "input2.txt"

    rule compute1:
        output: "input1.txt"
        shell: "touch input1.txt"

    rule compute2:
        output: "input2.txt"
        shell: "touch input2.txt"

This Snakefile can then be configured by a corresponding cluster config, say "cluster.json":


.. code-block:: json

    {
        "__default__" :
        {
            "account" : "my account",
            "time" : "00:15:00",
            "n" : 1,
            "partition" : "core"
        },
        "compute1" :
        {
            "time" : "00:20:00"
        }
    }

Any string in the cluster configuration can be formatted in the same way as shell commands, e.g. ``{rule}.{wildcards.sample}`` is formatted to ``a.xy`` if the rulename is ``a`` and the wildcard value is ``xy``.
Here ``__default__`` is a special object that specifies default parameters, these will be inherited by the other configuration objects. The ``compute1`` object here changes the ``time`` parameter, but keeps the other parameters from ``__default__``. The rule ``compute2`` does not have any configuration, and will therefore use the default configuration. You can then run the Snakefile with the following command on a SLURM system.

.. code-block:: console

    $ snakemake -j 999 --cluster-config cluster.json --cluster "sbatch -A {cluster.account} -p {cluster.partition} -n {cluster.n}  -t {cluster.time}"


For cluster systems using LSF/BSUB, a cluster config may look like this:

.. code-block:: json

    {
        "__default__" :
        {
            "queue"     : "medium_priority",
            "nCPUs"     : "16",
            "memory"    : 20000,
            "resources" : "\"select[mem>20000] rusage[mem=20000] span[hosts=1]\"",
            "name"      : "JOBNAME.{rule}.{wildcards}",
            "output"    : "logs/cluster/{rule}.{wildcards}.out",
            "error"     : "logs/cluster/{rule}.{wildcards}.err"
        },


        "trimming_PE" :
        {
            "memory"    : 30000,
            "resources" : "\"select[mem>30000] rusage[mem=30000] span[hosts=1]\"",
        }
    }

The advantage of this setup is that it is already pretty general by exploiting the wildcard possibilities that Snakemake provides via ``{rule}`` and ``{wildcards}``. So job names, output and error files all have reasonable and trackable default names, only the directies (``logs/cluster``) and job names (``JOBNAME``) have to adjusted accordingly.
If a rule named ``bamCoverage`` is executed with the wildcard ``basename = sample1``, for example, the output and error files will be ``bamCoverage.basename=sample1.out`` and ``bamCoverage.basename=sample1.err``, respectively.


---------------------------
Configure Working Directory
---------------------------

All paths in the snakefile are interpreted relative to the directory snakemake is executed in. This behaviour can be overridden by specifying a workdir in the snakefile:

.. code-block:: python

    workdir: "path/to/workdir"

Usually, it is preferred to only set the working directory via the command line, because above directive limits the portability of Snakemake workflows.
.. snakefiles-modularization:

.. _Snakemake Wrapper Repository: https://snakemake-wrappers.readthedocs.io

==============
Modularization
==============

Modularization in Snakemake comes at four different levels.

1. The most fine-grained level are wrappers. They are available and can be published at the `Snakemake Wrapper Repository`_. These wrappers can then be composed and customized according to your needs, by copying skeleton rules into your workflow. In combination with conda integration, wrappers also automatically deploy the needed software dependencies into isolated environments.
2. For larger, reusable parts that shall be integrated into a common workflow, it is recommended to write small Snakefiles and include them into a main Snakefile via the include statement. In such a setup, all rules share a common config file.
3. The third level is provided via the :ref:`module statement <snakefiles-modules>`, which enables arbitrary combination and reuse of rules.
4. Finally, Snakemake provides a syntax for defining :ref:`subworkflows <snakefiles-sub_workflows>`, which is however deprecated in favor of the module statement.


.. _snakefiles-wrappers:

--------
Wrappers
--------

The wrapper directive allows to have re-usable wrapper scripts around e.g. command line tools.
In contrast to modularization strategies like ``include`` or subworkflows, the wrapper directive allows to re-wire the DAG of jobs.
For example

.. code-block:: python

    rule samtools_sort:
        input:
            "mapped/{sample}.bam"
        output:
            "mapped/{sample}.sorted.bam"
        params:
            "-m 4G"
        threads: 8
        wrapper:
            "0.0.8/bio/samtools/sort"

.. note::

    It is possible to refer to wildcards and params in the wrapper identifier, e.g. by specifying ``"0.0.8/bio/{params.wrapper}"`` or ``"0.0.8/bio/{wildcards.wrapper}"``.

Refers to the wrapper ``"0.0.8/bio/samtools/sort"`` to create the output from the input.
Snakemake will automatically download the wrapper from the `Snakemake Wrapper Repository`_.
Thereby, ``0.0.8`` can be replaced with the git `version tag <https://github.com/snakemake/snakemake-wrappers/releases>`_ you want to use, or a `commit id <https://github.com/snakemake/snakemake-wrappers/commits>`_.
This ensures reproducibility since changes in the wrapper implementation will only be propagated to your workflow once you update the version tag.
Examples for each wrapper can be found in the READMEs located in the wrapper subdirectories at the `Snakemake Wrapper Repository`_.

Alternatively, for example during development, the wrapper directive can also point to full URLs, including URLs to local files with absolute paths ``file://`` or relative paths ``file:``.
Such a URL will have to point to the folder containing the ``wrapper.*`` and ``environment.yaml`` files.
In the above example, the full GitHub URL could for example be provided with ``wrapper: https://github.com/snakemake/snakemake-wrappers/raw/0.0.8/bio/samtools/sort``.
Note that it needs to point to the ``/raw/`` version of the folder, not the rendered HTML version.

In addition, the `Snakemake Wrapper Repository`_ offers so-called meta-wrappers, which can be used as modules, see :ref:`snakefiles-meta-wrappers`.

The `Snakemake Wrapper Repository`_ is meant as a collaborative project and pull requests are very welcome.


.. _cwl:

--------------------------------------
Common-Workflow-Language (CWL) support
--------------------------------------

With Snakemake 4.8.0, it is possible to refer to `CWL <https://www.commonwl.org/>`_ tool definitions in rules instead of specifying a wrapper or a plain shell command.
A CWL tool definition can be used as follows.

.. code-block:: python

    rule samtools_sort:
        input:
            input="mapped/{sample}.bam"
        output:
            output_name="mapped/{sample}.sorted.bam"
        params:
            threads=lambda wildcards, threads: threads,
            memory="4G"
        threads: 8
        cwl:
            "https://github.com/common-workflow-language/workflows/blob/"
            "fb406c95/tools/samtools-sort.cwl"

.. note::

    It is possible to refer to wildcards and params in the tool definition URL, e.g. by specifying something like ``"https://.../tools/{params.tool}.cwl"`` or ``"https://.../tools/{wildcards.tool}.cwl"``.

It is advisable to use a github URL that includes the commit as above instead of a branch name, in order to ensure reproducible results.
Snakemake will execute the rule by invoking `cwltool`, which has to be available via your `$PATH` variable, and can be, e.g., installed via `conda` or `pip`.
When using in combination with :ref:`--use-singularity <singularity>`, Snakemake will instruct `cwltool` to execute the command via Singularity in user space.
Otherwise, `cwltool` will in most cases use a Docker container, which requires Docker to be set up properly.

The advantage is that predefined tools available via any `repository of CWL tool definitions <https://www.commonwl.org/#Repositories_of_CWL_Tools_and_Workflows>`_ can be used in any supporting workflow management system.
In contrast to a :ref:`Snakemake wrapper <snakefiles-wrappers>`, CWL tool definitions are in general not suited to alter the behavior of a tool, e.g., by normalizing output names or special input handling.
As you can see in comparison to the analog :ref:`wrapper declaration <snakefiles-wrappers>` above, the rule becomes slightly more verbose, because input, output, and params have to be dispatched to the specific expectations of the CWL tool definition.

.. _snakefiles-includes:

--------
Includes
--------

Another Snakefile with all its rules can be included into the current:

.. code-block:: python

    include: "path/to/other/snakefile"

The default target rule (often called the ``all``-rule), won't be affected by the include.
I.e. it will always be the first rule in your Snakefile, no matter how many includes you have above your first rule.
Includes are relative to the directory of the Snakefile in which they occur.
For example, if above Snakefile resides in the directory ``my/dir``, then Snakemake will search for the include at ``my/dir/path/to/other/snakefile``, regardless of the working directory.


.. _snakefiles-modules:

-------
Modules
-------

With Snakemake 6.0 and later, it is possible to define external workflows as modules, from which rules can be used by explicitly "importing" them.

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        snakefile:
            # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
            "other_workflow/Snakefile"
    
    use rule * from other_workflow as other_*

The ``module other_workflow:`` statement registers the external workflow as a module, by defining the path to the main snakefile of ``other_workflow``.
Here, plain paths, HTTP/HTTPS URLs and special markers for code hosting providers like Github or Gitlab are possible (see :ref:`snakefile-code-hosting-providers`).
The second statement, ``use rule * from other_workflow as other_*``, declares all rules of that module to be used in the current one.
Thereby, the ``as other_*`` at the end renames all those rules with a common prefix.
This can be handy to avoid rule name conflicts (note that rules from modules can otherwise overwrite rules from your current workflow or other modules).
The module is evaluated in a separate namespace, and only the selected rules are added to the current workflow.
Non-rule Python statements inside the module are also evaluated in that separate namespace.
They are available in the module-defining workflow under the name of the module (e.g. here ``other_workflow.myfunction()`` would call the function ``myfunction`` that has been defined in the model, e.g. in ``other_workflow/Snakefile``).

It is possible to overwrite the global config dictionary for the module, which is usually filled by the ``configfile`` statement (see :ref:`snakefiles_standard_configuration`):

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]
    
    use rule * from other_workflow as other_*

In this case, any ``configfile`` statements inside the module are ignored.
In addition, it is possible to skip any :ref:`validation <snakefiles_config_validation>` statements in the module, by specifying ``skip_validation: True`` in the module statment.
Moreover, one can automatically move all relative input and output files of a module into a dedicated folder: by specifying ``prefix: "foo"`` in the module definition, e.g. any output file ``path/to/output.txt`` in the module would be stored under ``foo/path/to/output.txt`` instead.
This becomes particularly usefull when combining multiple modules, see :ref:`use_with_modules`.

Instead of using all rules, it is possible to import specific rules.
Specific rules may even be modified before using them, via a final ``with:`` followed by a block that lists items to overwrite.
This modification can be performed after a general import, and will overwrite any unmodified import of the same rule.

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    module other_workflow:
        # here, plain paths, URLs and the special markers for code hosting providers (see below) are possible.
        snakefile: "other_workflow/Snakefile"
        config: config["other-workflow"]

    use rule * from other_workflow as other_*

    use rule some_task from other_workflow as other_some_task with:
        output:
            "results/some-result.txt"

By such a modifying use statement, any properties of the rule (``input``, ``output``, ``log``, ``params``, ``benchmark``, ``threads``, ``resources``, etc.) can be overwritten, except the actual execution step (``shell``, ``notebook``, ``script``, ``cwl``, or ``run``).

Note that the second use statement has to use the original rule name, not the one that has been prefixed with ``other_`` via the first use statement (there is no rule ``other_some_task`` in the module ``other_workflow``).
In order to overwrite the rule ``some_task`` that has been imported with the first ``use rule`` statement, it is crucial to ensure that the rule is used with the same name in the second statement, by adding an equivalent ``as`` clause (here ``other_some_task``).
Otherwise, you will have two versions of the same rule, which might be unintended (a common symptom of such unintended repeated uses would be ambiguous rule exceptions thrown by Snakemake).

Of course, it is possible to combine the use of rules from multiple modules (see :ref:`use_with_modules`), and via modifying statements they can be rewired and reconfigured in an arbitrary way.

..  _snakefiles-meta-wrappers:

~~~~~~~~~~~~~
Meta-Wrappers
~~~~~~~~~~~~~

Snakemake wrappers offer a simple way to include commonly used tools in Snakemake workflows.
In addition the `Snakemake Wrapper Repository`_ offers so-called meta-wrappers, which are combinations of wrappers, meant to perform common tasks.
Both wrappers and meta-wrappers are continously tested.
The module statement also allows to easily use meta-wrappers, for example:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config.yaml"


    module bwa_mapping:
        meta_wrapper: "0.72.0/meta/bio/bwa_mapping"


    use rule * from bwa_mapping


    def get_input(wildcards):
        return config["samples"][wildcards.sample]


    use rule bwa_mem from bwa_mapping with:
        input:
            get_input


First, we define the meta-wrapper as a module.
Next, we declare all rules from the module to be used.
And finally, we overwrite the input directive of the rule ``bwa_mem`` such that the raw data is taken from the place where our workflow configures it via it's config file.

.. _snakefiles-sub_workflows:

-------------
Sub-Workflows
-------------

In addition to including rules of another workflow, Snakemake allows to depend on the output of other workflows as sub-workflows.
A sub-workflow is executed independently before the current workflow is executed.
Thereby, Snakemake ensures that all files the current workflow depends on are created or updated if necessary.
This allows to create links between otherwise separate data analyses.

.. code-block:: python

    subworkflow otherworkflow:
        workdir:
            "../path/to/otherworkflow"
        snakefile:
            "../path/to/otherworkflow/Snakefile"
        configfile:
            "path/to/custom_configfile.yaml"

    rule a:
        input:
            otherworkflow("test.txt")
        output: ...
        shell:  ...

Here, the subworkflow is named "otherworkflow" and it is located in the working directory ``../path/to/otherworkflow``.
The snakefile is in the same directory and called ``Snakefile``.
If ``snakefile`` is not defined for the subworkflow, it is assumed be located in the workdir location and called ``Snakefile``, hence, above we could have left the ``snakefile`` keyword out as well.
If ``workdir`` is not specified, it is assumed to be the same as the current one.
The (optional) definition of a ``configfile`` allows to parameterize the subworkflow as needed.
Files that are output from the subworkflow that we depend on are marked with the ``otherworkflow`` function (see the input of rule a).
This function automatically determines the absolute path to the file (here ``../path/to/otherworkflow/test.txt``).

When executing, snakemake first tries to create (or update, if necessary) ``test.txt`` (and all other possibly mentioned dependencies) by executing the subworkflow.
Then the current workflow is executed.
This can also happen recursively, since the subworkflow may have its own subworkflows as well.


.. _snakefile-code-hosting-providers:

----------------------
Code hosting providers
----------------------

To obtain the correct URL to an external source code resource (e.g. a snakefile, see :ref:`snakefiles-modules`), Snakemake provides markers for code hosting providers.
Currently, Github 

.. code-block:: python

    github("owner/repo", path="workflow/Snakefile", tag="v1.0.0")


and Gitlab are supported:

.. code-block:: python

    gitlab("owner/repo", path="workflow/Snakefile", tag="v1.0.0")

For the latter, it is also possible to specify an alternative host, e.g.

.. code-block:: python

    gitlab("owner/repo", path="workflow/Snakefile", tag="v1.0.0", host="somecustomgitlab.org")


While specifying a tag is highly encouraged, it is alternatively possible to specify a `commit` or a `branch` via respective keyword arguments.
Note that only when specifying a tag or a commit, Snakemake is able to persistently cache the source, thereby avoiding to repeatedly query it in case of multiple executions.
.. _snakefiles-remote_files:

============
Remote files
============

In versions ``snakemake>=3.5``.

The ``Snakefile`` supports a wrapper function, ``remote()``, indicating a file is on a remote storage provider (this is similar to ``temp()`` or ``protected()``). In order to use all types of remote files, the Python packages ``boto``, ``moto``, ``filechunkio``, ``pysftp``, ``dropbox``, ``requests``, ``ftputil``, ``XRootD``, and ``biopython`` must be installed.

During rule execution, a remote file (or object) specified is downloaded to the local ``cwd``, within a sub-directory bearing the same name as the remote provider. This sub-directory naming lets you have multiple remote origins with reduced likelihood of name collisions, and allows Snakemake to easily translate remote objects to local file paths. You can think of each local remote sub-directory as a local mirror of the remote system. The ``remote()`` wrapper is mutually-exclusive with the ``temp()`` and ``protected()`` wrappers.

Snakemake includes the following remote providers, supported by the corresponding classes:

* Amazon Simple Storage Service (AWS S3): ``snakemake.remote.S3``
* Google Cloud Storage (GS): ``snakemake.remote.GS``
* Microsoft Azure Blob Storage: ``snakemake.remote.AzBlob``
* File transfer over SSH (SFTP): ``snakemake.remote.SFTP``
* Read-only web (HTTP[S]): ``snakemake.remote.HTTP``
* File transfer protocol (FTP): ``snakemake.remote.FTP``
* Dropbox: ``snakemake.remote.dropbox``
* XRootD: ``snakemake.remote.XRootD``
* GenBank / NCBI Entrez: ``snakemake.remote.NCBI``
* WebDAV: ``snakemake.remote.webdav``
* GFAL: ``snakemake.remote.gfal``
* GridFTP: ``snakemake.remote.gridftp``
* iRODS: ``snakemake.remote.iRODS``
* EGA: ``snakemake.remote.EGA``
* AUTO: an automated remote selector

Amazon Simple Storage Service (S3)
==================================

This section describes usage of the S3 RemoteProvider, and also provides an intro to remote files and their usage.

It is important to note that you must have credentials (``access_key_id`` and ``secret_access_key``) which permit read/write access. If a file only serves as input to a Snakemake rule, read access is sufficient. You may specify credentials as environment variables or in the file ``~/.aws/credentials``, prefixed with ``AWS_*``, as with a standard `boto config <https://boto.readthedocs.org/en/latest/boto_config_tut.html>`_. Credentials may also be explicitly listed in the ``Snakefile``, as shown below:

For the Amazon S3 and Google Cloud Storage providers, the sub-directory used must be the bucket name.

Using remote files is easy (AWS S3 shown):

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input:
            S3.remote("bucket-name/file.txt")

Expand still works as expected, just wrap the expansion:


.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider()

    rule all:
        input:
            S3.remote(expand("bucket-name/{letter}-2.txt", letter=["A", "B", "C"]))

Only remote files needed to satisfy the DAG build are downloaded for the workflow. By default, remote files are downloaded prior to rule execution and are removed locally as soon as no rules depend on them. Remote files can be explicitly kept by setting the ``keep_local=True`` keyword argument:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input: S3.remote('bucket-name/prefix{split_id}.txt', keep_local=True)

If you wish to have a rule to simply download a file to a local copy, you can do so by declaring the same file path locally as is used by the remote file:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")

    rule all:
        input:
            S3.remote("bucket-name/out.txt")
        output:
            "bucket-name/out.txt"
        run:
            shell("cp {output[0]} ./")

In some cases the rule can use the data directly on the remote provider, in these cases ``stay_on_remote=True`` can be set to avoid downloading/uploading data unnecessarily. Additionally, if the backend supports it, any potentially corrupt output files will be removed from the remote. The default for ``stay_on_remote`` and ``keep_local`` can be configured by setting these properties on the remote provider object:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET", keep_local=True, stay_on_remote=True)

The remote provider also supports a new ``glob_wildcards()`` (see :ref:`glob-wildcards`) which acts the same as the local version of ``glob_wildcards()``, but for remote files:

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEY", secret_access_key="MYSECRET")
    S3.glob_wildcards("bucket-name/{file_prefix}.txt")

    # (result looks just like as if the local glob_wildcards() function were used on a locally with a folder called "bucket-name")

If the AWS CLI is installed it is possible to configure your keys globally. This removes the necessity of hardcoding the keys in the Snakefile. The interactive AWS credentials setup can be done using the following command:

.. code-block:: python

    aws configure

S3 then can be used without the keys.

.. code-block:: python

    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider
    S3 = S3RemoteProvider()

Finally, it is also possible to overwrite the S3 host via adding a ``host`` argument (taking a URL string) to ``S3RemoteProvider``.

Google Cloud Storage (GS)
=========================

Usage of the GS provider is the same as the S3 provider.
For authentication, one simply needs to login via the ``gcloud`` tool before
executing Snakemake, i.e.:

.. code-block:: console

    $ gcloud auth application-default login

In the Snakefile, no additional authentication information has to be provided:

.. code-block:: python

    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    GS = GSRemoteProvider()

    rule all:
        input:
            GS.remote("bucket-name/file.txt")


Microsoft Azure Blob Storage
=============================

Usage of the Azure Blob Storage provider is similar to the S3 provider. For
authentication, an account name and shared access signature (SAS) or key can be used. If these
variables are not passed directly to AzureRemoteProvider (see
[BlobServiceClient
class](https://docs.microsoft.com/en-us/python/api/azure-storage-blob/azure.storage.blob.blobserviceclient?view=azure-python)
for naming), they will be read from environment variables, named
`AZ_BLOB_ACCOUNT_URL` and `AZ_BLOB_CREDENTIAL`. `AZ_BLOB_ACCOUNT_URL` takes the form
`https://<accountname>.blob.core.windows.net` and may also contain a SAS. If
a SAS is not part of the URL, `AZ_BLOB_CREDENTIAL` has to be set to the SAS or alternatively to
the storage account key.

When using AzBlob as default remote provider you will almost always want to
pass these environment variables on to the remote execution environment (e.g.
Kubernetes) with `--envvars`, e.g
`--envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL`.

.. code-block:: python

    from snakemake.remote.AzBlob import RemoteProvider as AzureRemoteProvider
    AS = AzureRemoteProvider()# assumes env vars AZ_BLOB_ACCOUNT_URL and possibly AZ_BLOB_CREDENTIAL are set

    rule a:
        input:
            AS.remote("path/to/file.txt")




File transfer over SSH (SFTP)
=============================

Snakemake can use files on remove servers accessible via SFTP (i.e. most \*nix servers).
It uses `pysftp <https://pysftp.readthedocs.org/en/release_0.2.8/pysftp.html#pysftp.Connection>`_ for the underlying support of SFTP, so the same connection options exist.
Assuming you have SSH keys already set up for the server you are using in the ``Snakefile``, usage is simple:


.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

If you need to create the output directories in the remote server, you can specify ``mkdir_remote=True``  in the ``RemoteProvider`` constructor.

.. code-block:: python

   from snakemake.remote.SFTP import RemoteProvider
   SFTP = RemoteProvider(mkdir_remote=True)

   rule all:
       input:
           "/home/foo/bar.txt"
       output:
           SFTP.remote('example.com/home/foo/create/dir/bar.txt')
       shell:
           "cp {input} {output}"

The remote file addresses used must be specified with the host (domain or IP address) and the absolute path to the file on the remote server. A port may be specified if the SSH daemon on the server is listening on a port other than 22, in either the ``RemoteProvider`` or in each instance of ``remote()``:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(port=4040)

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

.. code-block:: python


    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()

    rule all:
        input:
            SFTP.remote("example.com:4040/path/to/file.bam")

The standard keyword arguments used by `pysftp <https://pysftp.readthedocs.org/en/release_0.2.8/pysftp.html#pysftp.Connection>`_ may be provided to the RemoteProvider to specify credentials (either password or private key):

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", private_key="/Users/myusername/.ssh/particular_id_rsa")

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            SFTP.remote("example.com/path/to/file.bam")

If you share credentials between servers but connect to one on a different port, the alternate port may be specified in the ``remote()`` wrapper:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            SFTP.remote("some-example-server-1.com/path/to/file.bam"),
            SFTP.remote("some-example-server-2.com:2222/path/to/file.bam")

There is a ``glob_wildcards()`` function:

.. code-block:: python

    from snakemake.remote.SFTP import RemoteProvider
    SFTP = RemoteProvider()
    SFTP.glob_wildcards("example.com/path/to/{sample}.bam")

Read-only web (HTTP[s])
=======================

Snakemake can access web resources via a read-only HTTP(S) provider.
This provider can be helpful for including public web data in a workflow.

Web addresses must be specified without protocol, so if your URI looks like this:

.. code-block:: text

    https://server3.example.com/path/to/myfile.tar.gz

The URI used in the ``Snakefile`` must look like this:

.. code-block:: text

    server3.example.com/path/to/myfile.tar.gz

It is straightforward to use the HTTP provider to download a file to the `cwd`:

.. code-block:: python

    import os
    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider

    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com/path/to/document.pdf", keep_local=True)
        run:
            outputName = os.path.basename(input[0])
            shell("mv {input} {outputName}")

To connect on a different port, specify the port as part of the URI string:

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com:8080/path/to/document.pdf", keep_local=True)

By default, the HTTP provider always uses HTTPS (TLS). If you need to connect to a resource with regular HTTP (no TLS), you must explicitly include ``insecure`` as a ``kwarg`` to ``remote()``:

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("www.example.com/path/to/document.pdf", insecure=True, keep_local=True)

If the URI used includes characters not permitted in a local file path, you may include them as part of the ``additional_request_string`` in the ``kwargs`` for ``remote()``. This may also be useful for including additional parameters you don not want to be part of the local filename (since the URI string becomes the local file name).

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/query.php", additional_request_string="?range=2;3")

If the file requires authentication, you can specify a username and password for HTTP Basic Auth with the Remote Provider, or with each instance of `remote()`.
For different types of authentication, you can pass in a Python ```requests.auth`` object (see `here <https://requests.readthedocs.io/en/master/api/#authentication>`_) the `auth` ``kwarg``.

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", keep_local=True)

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", username="myusername", password="mypassword", keep_local=True)

.. code-block:: python

    from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
    HTTP = HTTPRemoteProvider()

    rule all:
        input:
            HTTP.remote("example.com/interactive.php", auth=requests.auth.HTTPDigestAuth("myusername", "mypassword"), keep_local=True)

Since remote servers do not present directory contents uniformly, ``glob_wildcards()`` is __not__ supported by the HTTP provider.

File Transfer Protocol (FTP)
============================

Snakemake can work with files stored on regular FTP.
Currently supported are authenticated FTP and anonymous FTP, excluding FTP via TLS.

Usage is similar to the SFTP provider, however the paths specified are relative to the FTP home directory (since this is typically a chroot):

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            FTP.remote("example.com/rel/path/to/file.tar.gz")

The port may be specified in either the provider, or in each instance of `remote()`:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword", port=2121)

    rule all:
        input:
            FTP.remote("example.com/rel/path/to/file.tar.gz")

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider

    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    rule all:
        input:
            FTP.remote("example.com:2121/rel/path/to/file.tar.gz")

Anonymous download of FTP resources is possible:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()

    rule all:
        input:
            # only keeping the file so we can move it out to the cwd
            FTP.remote("example.com/rel/path/to/file.tar.gz", keep_local=True)
        run:
            shell("mv {input} ./")

``glob_wildcards()``:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    print(FTP.glob_wildcards("example.com/somedir/{file}.txt"))

Setting `immediate_close=True` allows the use of a large number of remote FTP input files in a job where the endpoint server limits the number of concurrent connections. When `immediate_close=True`, Snakemake will terminate FTP connections after each remote file action (`exists()`, `size()`, `download()`, `mtime()`, etc.). This is in contrast to the default behavior which caches FTP details and leaves the connection open across actions to improve performance (closing the connection upon job termination).  :

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider()

    rule all:
        input:
            # only keep the file so we can move it out to the cwd
            # This server limits the number of concurrent connections so we need to have Snakemake close each after each FTP action.
            FTP.remote(expand("ftp.example.com/rel/path/to/{file}", file=large_list), keep_local=True, immediate_close=True)
        run:
            shell("mv {input} ./")

``glob_wildcards()``:

.. code-block:: python

    from snakemake.remote.FTP import RemoteProvider as FTPRemoteProvider
    FTP = FTPRemoteProvider(username="myusername", password="mypassword")

    print(FTP.glob_wildcards("example.com/somedir/{file}.txt"))

Dropbox
=======

The Dropbox remote provider allows you to upload and download from your `Dropbox <https://www.dropbox.com>`_ account without having the client installed on your machine. In order to use the provider you  first need to register an "app" on the `Dropbox developer website <https://www.dropbox.com/developers/apps/create>`_, with access to the Full Dropbox. After registering, generate an OAuth2 access token. You will need the token to use the Snakemake Dropbox remote provider.

Using the Dropbox provider is straightforward:

.. code-block:: python

    from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
    DBox = DropboxRemoteProvider(oauth2_access_token="mytoken")

    rule all:
        input:
            DBox.remote("path/to/input.txt")

``glob_wildcards()`` is supported:

.. code-block:: python

    from snakemake.remote.dropbox import RemoteProvider as DropboxRemoteProvider
    DBox = DropboxRemoteProvider(oauth2_access_token="mytoken")

    DBox.glob_wildcards("path/to/{title}.txt")

Note that Dropbox paths are case-insensitive.

XRootD
=======

Snakemake can be used with `XRootD <https://xrootd.slac.stanford.edu/>`_ backed storage provided the python bindings are installed.
This is typically most useful when combined with the ``stay_on_remote`` flag to minimise local storage requirements.
This flag can be overridden on a file by file basis as described in the S3 remote. Additionally ``glob_wildcards()`` is supported:

.. code-block:: python

    from snakemake.remote.XRootD import RemoteProvider as XRootDRemoteProvider

    XRootD = XRootDRemoteProvider(stay_on_remote=True)
    file_numbers = XRootD.glob_wildcards("root://eospublic.cern.ch//eos/opendata/lhcb/MasterclassDatasets/D0lifetime/2014/mclasseventv2_D0_{n}.root").n

    rule all:
        input:
            expand("local_data/mclasseventv2_D0_{n}.root", n=file_numbers)

    rule make_data:
        input:
            XRootD.remote("root://eospublic.cern.ch//eos/opendata/lhcb/MasterclassDatasets/D0lifetime/2014/mclasseventv2_D0_{n}.root")
        output:
            'local_data/mclasseventv2_D0_{n}.root'
        shell:
            'xrdcp {input[0]} {output[0]}'

GenBank / NCBI Entrez
=====================

Snakemake can directly source input files from `GenBank <https://www.ncbi.nlm.nih.gov/genbank/>`_ and other `NCBI Entrez databases <https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly>`_ if the Biopython library is installed.

.. code-block:: python

    from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
    NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse

    rule all:
        input:
            "size.txt"

    rule download_and_count:
        input:
            NCBI.remote("KY785484.1.fasta", db="nuccore")
        output:
            "size.txt"
        run:
            shell("wc -c {input} > {output}")

The output format and source database of a record retrieved from GenBank is inferred from the file extension specified. For example, ``NCBI.RemoteProvider().remote("KY785484.1.fasta", db="nuccore")`` will download a FASTA file while ``NCBI.RemoteProvider().remote("KY785484.1.gb", db="nuccore")`` will download a GenBank-format file. If the options are ambiguous, Snakemake will raise an exception and inform the user of possible format choices. To see available formats, consult the `Entrez EFetch documentation <https://www.ncbi.nlm.nih.gov/books/NBK25499/table/chapter4.T._valid_values_of__retmode_and/?report=objectonly>`_. To view the valid file extensions for these formats, access ``NCBI.RemoteProvider()._gb.valid_extensions``, or instantiate an ``NCBI.NCBIHelper()`` and access ``NCBI.NCBIHelper().valid_extensions`` (this is a property).

When used in conjunction with ``NCBI.RemoteProvider().search()``, Snakemake and ``NCBI.RemoteProvider().remote()`` can be used to find accessions by query and download them:

.. code-block:: python

    from snakemake.remote.NCBI import RemoteProvider as NCBIRemoteProvider
    NCBI = NCBIRemoteProvider(email="someone@example.com") # email required by NCBI to prevent abuse

    # get accessions for the first 3 results in a search for full-length Zika virus genomes
    # the query parameter accepts standard GenBank search syntax
    query = '"Zika virus"[Organism] AND (("9000"[SLEN] : "20000"[SLEN]) AND ("2017/03/20"[PDAT] : "2017/03/24"[PDAT])) '
    accessions = NCBI.search(query, retmax=3)

    # give the accessions a file extension to help the RemoteProvider determine the
    # proper output type.
    input_files = expand("{acc}.fasta", acc=accessions)

    rule all:
        input:
            "sizes.txt"

    rule download_and_count:
        input:
            # Since *.fasta files could come from several different databases, specify the database here.
            # if the input files are ambiguous, the provider will alert the user with possible options
            # standard options like "seq_start" are supported
            NCBI.remote(input_files, db="nuccore", seq_start=5000)

        output:
            "sizes.txt"
        run:
            shell("wc -c {input} > sizes.txt")

Normally, all accessions for a query are returned from ``NCBI.RemoteProvider.search()``. To truncate the results, specify ``retmax=<desired_number>``. Standard Entrez `fetch query options <https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch>`_ are supported as kwargs, and may be passed in to ``NCBI.RemoteProvider.remote()`` and ``NCBI.RemoteProvider.search()``.

WebDAV
======

WebDAV support is currently ``experimental`` and available in Snakemake 4.0 and later.

Snakemake supports reading and writing WebDAV remote files. The protocol defaults to ``https://``, but insecure connections
can be used by specifying ``protocol=="http://"``. Similarly, the port defaults to 443, and can be overridden by specifying ``port=##`` or by including the port as part of the file address.

.. code-block:: python

    from snakemake.remote import webdav

    webdav = webdav.RemoteProvider(username="test", password="test", protocol="http://")

    rule a:
        input:
            webdav.remote("example.com:8888/path/to/input_file.csv"),
        shell:
            # do something


GFAL
====

GFAL support is available in Snakemake 4.1 and later.

Snakemake supports reading and writing remote files via the `GFAL <https://dmc.web.cern.ch/projects/gfal-2/home>`_ command line client (gfal-* commands).
By this, it supports various grid storage protocols like `GridFTP <https://en.wikipedia.org/wiki/GridFTP>`_.
In general, if you are able to use the `gfal-*` commands directly, Snakemake support for GFAL will work as well.

.. code-block:: python

    from snakemake.remote import gfal

    gfal = gfal.RemoteProvider(retry=5)

    rule a:
        input:
            gfal.remote("gridftp.grid.sara.nl:2811/path/to/infile.txt")
        output:
            gfal.remote("gridftp.grid.sara.nl:2811/path/to/outfile.txt")
        shell:
            # do something

Authentication has to be setup in the system, e.g. via certificates in the ``.globus`` directory.
Usually, this is already the case and no action has to be taken.
The keyword argument to the remote provider allows to set the number of retries (10 per default) in case of failed commands (the GRID is usually relatively unreliable).
The latter may be unsupported depending on the system configuration.

Note that GFAL support used together with the flags ``--no-shared-fs`` and ``--default-remote-provider`` enables you
to transparently use Snakemake in a grid computing environment without a shared network filesystem.
For an example see the `surfsara-grid configuration profile <https://github.com/Snakemake-Profiles/surfsara-grid>`_.

GridFTP
=======

GridFTP support is available in Snakemake 4.3.0 and later.

As a more specialized alternative to the GFAL remote provider, Snakemake provides a `GridFTP <https://en.wikipedia.org/wiki/GridFTP>`_ remote provider.
This provider only supports the GridFTP protocol. Internally, it uses the `globus-url-copy <http://toolkit.globus.org/toolkit/docs/latest-stable/gridftp/user/#globus-url-copy>`_ command for downloads and uploads, while all other tasks are delegated to the GFAL remote provider.

.. code-block:: python

    from snakemake.remote import gridftp

    gridftp = gridftp.RemoteProvider(retry=5)

    rule a:
        input:
            gridftp.remote("gridftp.grid.sara.nl:2811/path/to/infile.txt")
        output:
            gridftp.remote("gridftp.grid.sara.nl:2811/path/to/outfile.txt")
        shell:
            # do something

Authentication has to be setup in the system, e.g. via certificates in the ``.globus`` directory.
Usually, this is already the case and no action has to be taken.
The keyword argument to the remote provider allows to set the number of retries (10 per default) in case of failed commands (the GRID is usually relatively unreliable).
The latter may be unsupported depending on the system configuration.

Note that GridFTP support used together with the flags ``--no-shared-fs`` and ``--default-remote-provider`` enables you
to transparently use Snakemake in a grid computing environment without a shared network filesystem.
For an example see the `surfsara-grid configuration profile <https://github.com/Snakemake-Profiles/surfsara-grid>`_.


Remote cross-provider transfers
===============================

It is possible to use Snakemake to transfer files between remote providers (using the local machine as an intermediary), as long as the sub-directory (bucket) names differ:

.. code-block:: python

    from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
    from snakemake.remote.S3 import RemoteProvider as S3RemoteProvider

    GS = GSRemoteProvider(access_key_id="MYACCESSKEYID", secret_access_key="MYSECRETACCESSKEY")
    S3 = S3RemoteProvider(access_key_id="MYACCESSKEYID", secret_access_key="MYSECRETACCESSKEY")

    fileList, = S3.glob_wildcards("source-bucket/{file}.bam")
    rule all:
        input:
            GS.remote( expand("destination-bucket/{file}.bam", file=fileList) )
    rule transfer_S3_to_GS:
        input:
            S3.remote( expand("source-bucket/{file}.bam", file=fileList) )
        output:
            GS.remote( expand("destination-bucket/{file}.bam", file=fileList) )
        run:
            shell("cp {input} {output}")


iRODS
=====

You can access an iRODS server to retrieve data from and upload data to it.
If your iRODS server is not set to a certain timezone, it is using UTC. It is
advised to shift the modification time provided by iRODS (``modify_time``)
then to your timezone by providing the ``timezone`` parameter such that
timestamps coming from iRODS are converted to the correct time.

iRODS actually does not save the timestamp from your original file but creates
its own timestamp of the upload time. When iRODS downloads the file for
processing, it does not take the timestamp from the remote file. Instead,
the file will have the timestamp when it was downloaded. To get around this,
we create a metadata entry to store the original file stamp from your system
and alter the timestamp of the downloaded file accordingly. While uploading,
the metadata entries ``atime``, ``ctime`` and ``mtime`` are added. When this
entry does not exist (because this module didn't upload the file), we fall back
to the timestamp provided by iRODS with the above mentioned strategy.

To access the iRODS server you need to have an iRODS environment configuration
file available and in this file the authentication needs to be configured.
The iRODS configuration file can be created by following the `official
instructions
<https://docs.irods.org/master/system_overview/configuration/#irodsirods_environmentjson>`_).

The default location for the configuration file is
``~/.irods/irods_environment.json``.  The ``RemoteProvider()`` class accepts
the parameter ``irods_env_file`` where an alternative path to the
``irods_environment.json`` file can be specified.  Another way is to export the
environment variable ``IRODS_ENVIRONMENT_FILE`` in your shell to specify the
location.

There are several ways to configure the authentication against the iRODS
server, depending on what your iRODS server offers. If you are using the
authentication via password, the default location of the authentication file is
``~/.irods/.irodsA``. Usually this file is generated with the ``iinit`` command
from the ``iCommands`` program suite. Inside the ``irods_environment.json``
file, the parameter ``"irods_authentication_file"`` can be set to specifiy an
alternative location for the ``.irodsA`` file. Another possibility to change
the location is to export the environment variable
``IRODS_AUTHENTICATION_FILE``.

The ``glob_wildcards()`` function is supported.

.. code-block:: python

    from snakemake.remote.iRODS import RemoteProvider

    irods = RemoteProvider(irods_env_file='setup-data/irods_environment.json',
                           timezone="Europe/Berlin") # all parameters are optional

    # please note the comma after the variable name!
    # access: irods.remote(expand('home/rods/{f}), f=files))
    files, = irods.glob_wildcards('home/rods/{files})

    rule all:
        input:
            irods.remote('home/rods/testfile.out'),

    rule gen:
        input:
            irods.remote('home/rods/testfile.in')
        output:
            irods.remote('home/rods/testfile.out')
        shell:
            r"""
            touch {output}
            """

An example for the iRODS configuration file (``irods_environment.json``):

.. code-block:: json

    {
        "irods_host": "localhost",
        "irods_port": 1247,
        "irods_user_name": "rods",
        "irods_zone_name": "tempZone",
        "irods_authentication_file": "setup-data/.irodsA"
    }


Please note that the ``zone`` folder is not included in the path as it will be
taken from the configuration file. The path also must not start with a ``/``.

By default, temporarily stored local files are removed. You can specify anyway
the parameter ``overwrite`` to tell iRODS to overwrite existing files that are
downloaded, because iRODS complains if a local file already exists when a
download attempt is issued (uploading is not a problem, though).

In the Snakemake source directory in ``snakemake/tests/test_remote_irods`` you
can find a working example.


EGA
===

The European Genome-phenome Archive (EGA) is a service for permanent archiving
and sharing of all types of personally identifiable genetic and phenotypic data
resulting from biomedical research projects.

From version 5.2 on, Snakemake provides experimental support to use EGA as a remote provider, such that
EGA hosted files can be transparently used as input.
For this to work, you need to define your username and password as environment
variables ``EGA_USERNAME`` and ``EGA_PASSWORD``.

Files in a dataset are addressed via the pattern ``ega/<dataset_id>/<filename>``.
Note that the filename should not include the ``.cip`` ending that is sometimes displayed in EGA listings:

.. code-block:: python

    import snakemake.remote.EGA as EGA

    ega = EGA.RemoteProvider()


    rule a:
        input:
            ega.remote("ega/EGAD00001002142/COLO_829_EPleasance_TGENPipe.bam.bai")
        output:
            "data/COLO_829BL_BCGSC_IlluminaPipe.bam.bai"
        shell:
            "cp {input} {output}"

Upon download, Snakemake will automatically decrypt the file and check the MD5 hash.


AUTO
====

A wrapper which automatically selects an appropriate remote provider based on the url's scheme.
It removes some of the boilerplate code required to download remote files from various providers:

.. code-block:: python

    from snakemake.remote import AUTO


    rule all:
        input:
            'foo'


    rule download:
        input:
            ftp_file_list=AUTO.remote([
                'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxcat.tar.gz',
                'ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
            ], keep_local=True),
            http_file=AUTO.remote(
                'https://github.com/hetio/hetionet/raw/master/hetnet/tsv/hetionet-v1.0-nodes.tsv'
            )
        output:
            touch('foo')
        shell:
            """
            head {input.http_file}
            """
.. _distribution_and_reproducibility:

================================
Distribution and Reproducibility
================================

It is recommended to store each workflow in a dedicated git repository of the
following structure:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    |   └── Snakefile
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources

In other words, the workflow code goes into a subfolder ``workflow``, while the configuration is stored in a subfolder ``config``. 
Inside of the ``workflow`` subfolder, the central ``Snakefile`` marks the entrypoint of the workflow (it will be automatically discovered when running snakemake from the root of above structure. 
In addition to the central ``Snakefile``, rules can be stored in a modular way, using the optional subfolder ``workflow/rules``. Such modules should end with ``.smk`` the recommended file extension of Snakemake.
Further, :ref:`scripts <snakefiles-external_scripts>` should be stored in a subfolder ``workflow/scripts`` and notebooks in a subfolder ``workflow/notebooks``.
Conda environments (see :ref:`integrated_package_management`) should be stored in a subfolder ``workflow/envs`` (make sure to keep them as finegrained as possible to improve transparency and maintainability).
Finally, :ref:`report caption files <snakefiles-reports>` should be stored in ``workflow/report``.
All output files generated in the workflow should be stored under ``results``, unless they are rather retrieved resources, in which case they should be stored under ``resources``. The latter subfolder may also contain small resources that shall be delivered along with the workflow via git (although it might be tempting, please refrain from trying to generate output file paths with string concatenation of a central ``outdir`` variable or so, as this hampers readability).

Workflows setup in above structure can be easily used and combined via :ref:`the Snakemake module system <use_with_modules>`.
Such deployment can even be automated via  `Snakedeploy <https://snakedeploy.readthedocs.io>`_.
Moreover, by publishing a workflow on `Github <https://github.com>`_ and following a set of additional `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_ the workflow will be automatically included in the `Snakemake workflow catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_, thereby easing discovery and even automating its usage documentation.
For an example of such automated documentation, see `here <https://snakemake.github.io/snakemake-workflow-catalog/?usage=snakemake-workflows%2Fdna-seq-varlociraptor>`_.

Visit the `Snakemake Workflows Project <https://github.com/snakemake-workflows/docs>`_ for more best-practice workflows.

.. _use_with_modules:

-----------------------------------------
Using and combining pre-exising workflows
-----------------------------------------

Via the :ref:`module/use <snakefiles-modules>` system introduced with Snakemake 6.0, it is very easy to deploy existing workflows for new projects.
This ranges from the simple application to new data to the complex combination of several complementary workflows in order to perfom an integrated analysis over multiple data types.

Consider the following example:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            # here, it is also possible to provide a plain raw URL like "https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile"
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config:
            config

    use rule * from dna_seq

First, we load a local configuration file.
Next, we define the module ``dna_seq`` to be loaded from the URL ``https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile``, while using the contents of the local configuration file.
Note that it is possible to either specify the full URL pointing to the raw Snakefile as a string or to use the github marker as done here.
With the latter, Snakemake can however cache the used source files persistently (if a tag is given), such that they don't have to be downloaded on each invocation.
Finally we declare all rules of the dna_seq module to be used.

This kind of deployment is equivalent to just cloning the original repository and modifying the configuration in it.
However, the advantage here is that we are (a) able to easily extend of modify the workflow, while making the changes transparent, and (b) we can store this workflow in a separate (e.g. private) git repository, along with for example configuration and meta data, without the need to duplicate the workflow code.
Finally, we are always able to later combine another module into the current workflow, e.g. when further kinds of analyses are needed.
The ability to modify rules upon using them (see :ref:`snakefiles-modules`) allows for arbitrary rewiring and configuration of the combined modules.

For example, we can easily add another rule to extend the given workflow:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            # here, it is also possible to provide a plain raw URL like "https://github.com/snakemake-workflows/dna-seq-gatk-variant-calling/raw/v2.0.1/workflow/Snakefile"
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config: config

    use rule * from dna_seq as dna_seq_*

    # easily extend the workflow
    rule plot_vafs:
        input:
            "filtered/all.vcf.gz"
        output:
            "results/plots/vafs.svg"
        notebook:
            "notebooks/plot-vafs.py.ipynb"

    # Define a new default target that collects both the targets from the dna_seq module as well as
    # the new plot.
    rule all:
        input:
            rules.dna_seq_all.input,
            "results/plots/vafs.svg",
        default_target: True

Above, we have added a prefix to all rule names of the dna_seq module, such that there is no name clash with the added rules (``as dna_seq_*`` in the ``use rule`` statement).
In addition, we have added a new rule ``all``, defining the default target in case the workflow is executed (as usually) without any specific target files or rule.
The new target rule collects both all input files of the rule ``all`` from the dna_seq workflow, as well as additionally collecting the new plot.

It is possible to further extend the workflow with other modules, thereby generating an integrative analysis.
Here, let us assume that we want to conduct another kind of analysis, say RNA-seq, using a different external workflow.
We can extend above example in the following way:

.. code-block:: python

    from snakemake.utils import min_version
    min_version("6.0")

    configfile: "config/config.yaml"

    module dna_seq:
        snakefile:
            github("snakemake-workflows/dna-seq-gatk-variant-calling", path="workflow/Snakefile", tag="v2.0.1")
        config: config["dna-seq"]
        prefix: "dna-seq"

    use rule * from dna_seq as dna_seq_*

    rule plot_vafs:
        input:
            "filtered/all.vcf.gz"
        output:
            "results/plots/vafs.svg"
        notebook:
            "notebooks/plot-vafs.py.ipynb"

    module rna_seq:
        snakefile:
            github("snakemake-workflows/rna-seq-kallisto-sleuth", path="workflow/Snakefile", tag="v2.0.1")
        config: config["rna-seq"]
        prefix: "rna-seq"

    use rule * from rna_seq as rna_seq_*


    # Define a new default target that collects all the targets from the dna_seq and rna_seq module.
    rule all:
        input:
            rules.dna_seq_all.input,
            rules.rna_seq_all.input,
        default_target: True

Above, several things have changed. 

* First, we have added another module ``rna_seq``.
* Second, we have added a prefix to all non-absolute input and output file names of both modules (``prefix: "dna-seq"`` and ``prefix: "rna-seq"``) in order to avoid file name clashes.
* Third, we have added a default target rule that collects both the default targets from the module ``dna_seq`` as well as the module ``rna_seq``.
* Finally, we provide the config of the two modules via two separate sections in the common config file (``config["dna-seq"]`` and ``config["rna-seq"]``).

----------------------------------
Uploading workflows to WorkflowHub
----------------------------------

In order to share a workflow with the scientific community it is advised to upload the repository to `WorkflowHub <https://workflowhub.eu/>`_, where each submission will be automatically parsed and encapsulated into a `Research Object Crate <https://w3id.org/ro/crate>`_. That way a *snakemake* workflow is annotated with proper metatada and thus complies with the `FAIR <https://en.wikipedia.org/wiki/FAIR_data>`_ principles of scientific data.

To adhere to the high WorkflowHub standards of scientific workflows the recommended *snakemake* repository structure presented above needs to be extended by the following elements:

- Code of Conduct
- Contribution instructions
- Workflow rule graph
- Workflow documentation
- Test directory

A code of conduct for the repository developers as well as instruction on how to contribute to the project should be placed in the top-level files: ``CODE_OF_CONDUCT.md`` and ``CONTRIBUTING.md``, respectively. Each *snakemake* workflow repository needs to contain an SVG-formatted rule graph placed in a subdirectory ``images/rulegraph.svg``. Additionally, the workflow should be annotated with a technical documentation of all of its subsequent steps, described in ``workflow/documentation.md``. Finally, the repository should contain a ``.tests`` directory with two subdirectories: ``.tests/integration`` and ``.tests/unit``. The former has to contain all the input data, configuration specifications and shell commands required to run an integration test of the whole workflow. The latter shall contain subdirectories dedicated to testing each of the separate workflow steps independently. To simplify the testing procedure *snakemake* can automatically generate unit tests from a successful workflow execution (see :ref:`snakefiles-testing`).

Therefore, the repository structure should comply with:

.. code-block:: none

    ├── .gitignore
    ├── README.md
    ├── LICENSE.md
    ├── CODE_OF_CONDUCT.md
    ├── CONTRIBUTING.md
    ├── .tests
    │   ├── integration
    │   └── unit
    ├── images
    │   └── rulegraph.svg
    ├── workflow
    │   ├── rules
    |   │   ├── module1.smk
    |   │   └── module2.smk
    │   ├── envs
    |   │   ├── tool1.yaml
    |   │   └── tool2.yaml
    │   ├── scripts
    |   │   ├── script1.py
    |   │   └── script2.R
    │   ├── notebooks
    |   │   ├── notebook1.py.ipynb
    |   │   └── notebook2.r.ipynb
    │   ├── report
    |   │   ├── plot1.rst
    |   │   └── plot2.rst
    │   ├── Snakefile
    |   └── documentation.md
    ├── config
    │   ├── config.yaml
    │   └── some-sheet.tsv
    ├── results
    └── resources


.. _integrated_package_management:

-----------------------------
Integrated Package Management
-----------------------------

With Snakemake 3.9.0 it is possible to define isolated software environments per rule.
Upon execution of a workflow, the `Conda package manager <https://conda.pydata.org>`_ is used to obtain and deploy the defined software packages in the specified versions. Packages will be installed into your working directory, without requiring any admin/root priviledges.
Given that conda is available on your system (see `Miniconda <https://conda.pydata.org/miniconda.html>`_), to use the Conda integration, add the ``--use-conda`` flag to your workflow execution command, e.g. ``snakemake --cores 8 --use-conda``.
When ``--use-conda`` is activated, Snakemake will automatically create software environments for any used wrapper (see :ref:`snakefiles-wrappers`).
Further, you can manually define environments via the ``conda`` directive, e.g.:

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

with the following `environment definition <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually>`_:


.. code-block:: yaml

    channels:
     - r
    dependencies:
     - r=3.3.1
     - r-ggplot2=2.1.0

The path to the environment definition is interpreted as **relative to the Snakefile that contains the rule** (unless it is an absolute path, which is discouraged).

Instead of using a concrete path, it is also possible to provide a path containing wildcards (which must also occur in the output files of the rule), analogous to the specification of input files.

.. sidebar:: Note

   Note that conda environments are only used with ``shell``, ``script`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.
   
   Further, note that search path modifying environment variables like ``R_LIBS`` and ``PYTHONPATH`` can interfere with your conda environments. 
   Therefore, Snakemake automatically deactivates them for a job when a conda environment definition is used.
   If you know what you are doing, in order to deactivate this behavior, you can use the flag ``--conda-not-block-search-path-envvars``.

Snakemake will store the environment persistently in ``.snakemake/conda/$hash`` with ``$hash`` being the MD5 hash of the environment definition file content. This way, updates to the environment definition are automatically detected.
Note that you need to clean up environments manually for now. However, in many cases they are lightweight and consist of symlinks to your central conda installation. 

Conda deployment also works well for offline or air-gapped environments. Running ``snakemake --use-conda --conda-create-envs-only`` will only install the required conda environments without running the full workflow. Subsequent runs with ``--use-conda`` will make use of the local environments without requiring internet access.


.. _conda_named_env:

Using already existing named conda environments
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Sometimes it can be handy to refer to an already existing named conda environment from a rule, instead of defining a new one from scratch.
Importantly, one should be aware that this can **hamper reproducibility**, because the workflow then relies on this environment to be present
**in exactly the same way** on any new system where the workflow is executed. Essentially, you will have to take care of this manually in such a case.
Therefore, the approach using environment definition files described above is highly recommended and preferred.

Nevertheless, in case you are still sure that you want to use an existing named environment, it can simply be put into the conda directive, e.g.

.. code-block:: python
    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "some-env-name"
        script:
            "scripts/plot-stuff.R"

For such a rule, Snakemake will just activate the given environment, instead of automatically deploying anything.
Instead of using a concrete name, it is also possible to provide a name containing wildcards (which must also occur in the output files of the rule), analogous to the specification of input files.

Note that Snakemake distinguishes file based environments from named ones as follows: 
if the given specification ends on ``.yaml`` or ``.yml``, Snakemake assumes it to be a path to an environment definition file; otherwise, it assumes the given specification
to be the name of an existing environment.

.. _singularity:


-------------------------
Providing post-deployment scripts
-------------------------

From Snakemake 6.14 onwards post-deployment shell-scripts can be provided to perform additional adjustments of a conda environment.
This might be helpful in case a conda package is missing components or requires further configuration for execution.
Post-deployment scripts must be placed next to their corresponding environment-file and require the suffix ``.post-deploy.sh``, e.g.:

.. code-block:: python

    rule NAME:
        input:
            "seqs.fastq"
        output:
            "results.tsv"
        conda:
            "envs/interproscan.yaml"
        shell:
            "interproscan.sh -i {input} -f tsv -o {output}"

.. code-block:: none

    ├── Snakefile
    └── envs
        ├── interproscan.yaml
        └── interproscan.post-deploy.sh

The path of the conda environment can be accessed within the script via ``$CONDA_PREFIX``.

--------------------------
Running jobs in containers
--------------------------

As an alternative to using Conda (see above), it is possible to define, for each rule, a (docker) container to use, e.g.,

.. code-block:: python

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        container:
            "docker://joseespinosa/docker-r-ggplot2"
        script:
            "scripts/plot-stuff.R"

When executing Snakemake with

.. code-block:: bash

    snakemake --use-singularity

it will execute the job within a container that is spawned from the given image.
Allowed image urls entail everything supported by singularity (e.g., ``shub://`` and ``docker://``).
However, ``docker://`` is preferred, as other container runtimes will be supported in the future (e.g. podman).

.. sidebar:: Note

   Note that singularity integration is only used with ``shell``, ``script`` and the ``wrapper`` directive, not the ``run`` directive.
   The reason is that the ``run`` directive has access to the rest of the Snakefile (e.g. globally defined variables) and therefore must be executed in the same process as Snakemake itself.


When ``--use-singularity`` is combined with ``--kubernetes`` (see :ref:`kubernetes`), cloud jobs will be automatically configured to run in priviledged mode, because this is a current requirement of the singularity executable.
Importantly, those privileges won't be shared by the actual code that is executed in the singularity container though.

A global definition of a container image can be given:

.. code-block:: python

    container: "docker://joseespinosa/docker-r-ggplot2"

    rule NAME:
        ...

In this case all jobs will be executed in a container. You can disable execution in container
by setting the container directive of the rule to ``None``.

.. code-block:: python

    container: "docker://joseespinosa/docker-r-ggplot2"

    rule NAME:
        container: None

-----------------------------------------
Containerization of Conda based workflows
-----------------------------------------
While :ref:`integrated_package_management` provides control over the used software in exactly
the desired versions, it does not control the underlying operating system.
However, given a workflow with conda environments for each rule, Snakemake can automatically
generate a container image specification (in the form of a ``Dockerfile``) that contains
all required environments via the flag --containerize:

.. code-block:: bash

    snakemake --containerize > Dockerfile

The container image specification generated by Snakemake aims to be transparent and readable, e.g. by displaying each contained environment in a human readable way.
Via the special directive ``containerized`` this container image can be used in the workflow (both globally or per rule) such that no further conda package downloads are necessary, for example:

.. code-block:: python

    containerized: "docker://username/myworkflow:1.0.0"

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

Using the containerization of Snakemake has three advantages over manually crafting a container image for a workflow:

1. A workflow with conda environment definitions is much more transparent to the reader than a black box container image, as each rule directly shows which software stack is used. Containerization just persistently projects those environments into a container image.
2. It remains possible to run the workflow without containers, just via the conda environments.
3. During development, testing can first happen without the container and just on the conda environments. When releasing a production version of the workflow the image can be uploaded just once and for future stable releases, thereby limiting the overhead created in container registries.

--------------------------------------------------------------
Ad-hoc combination of Conda package management with containers
--------------------------------------------------------------

While :ref:`integrated_package_management` provides control over the used software in exactly
the desired versions, it does not control the underlying operating system.
Here, it becomes handy that Snakemake >=4.8.0 allows to combine Conda-based package management
with :ref:`singularity`.
For example, you can write

.. code-block:: python

    container: "docker://continuumio/miniconda3:4.4.10"

    rule NAME:
        input:
            "table.txt"
        output:
            "plots/myplot.pdf"
        conda:
            "envs/ggplot.yaml"
        script:
            "scripts/plot-stuff.R"

in other words, a global definition of a container image can be combined with a
per-rule conda directive.
Then, upon invocation with

.. code-block:: bash

    snakemake --use-conda --use-singularity

Snakemake will first pull the defined container image, and then create the requested conda environment from within the container.
The conda environments will still be stored in your working environment, such that they don't have to be recreated unless they have changed.
The hash under which the environments are stored includes the used container image url, such that changes to the container image also lead to new environments to be created.
When a job is executed, Snakemake will first enter the container and then activate the conda environment.

By this, both packages and OS can be easily controlled without the overhead of creating and distributing specialized container images.
Of course, it is also possible (though less common) to define a container image per rule in this scenario.

The user can, upon execution, freely choose the desired level of reproducibility:

* no package management (use whatever is on the system)
* Conda based package management (use versions defined by the workflow developer)
* Conda based package management in containerized OS (use versions and OS defined by the workflow developer)

-------------------------
Using environment modules
-------------------------

In high performace cluster systems (HPC), it can be preferable to use environment modules for deployment of optimized versions of certain standard tools.
Snakemake allows to define environment modules per rule:

.. code-block:: python

    rule bwa:
        input:
            "genome.fa"
            "reads.fq"
        output:
            "mapped.bam"
        conda:
            "envs/bwa.yaml"
        envmodules:
            "bio/bwa/0.7.9",
            "bio/samtools/1.9"
        shell:
            "bwa mem {input} | samtools view -Sbh - > {output}"

Here, when Snakemake is executed with ``snakemake --use-envmodules``, it will load the defined modules in the given order, instead of using the also defined conda environment.
Note that although not mandatory, one should always provide either a conda environment or a container (see above), along with environment module definitions.
The reason is that environment modules are often highly platform specific, and cannot be assumed to be available somewhere else, thereby limiting reproducibility.
By definition an equivalent conda environment or container as a fallback, people outside of the HPC system where the workflow has been designed can still execute it, e.g. by running ``snakemake --use-conda`` instead of ``snakemake --use-envmodules``.

--------------------------------------
Sustainable and reproducible archiving
--------------------------------------

With Snakemake 3.10.0 it is possible to archive a workflow into a
`tarball <https://en.wikipedia.org/wiki/Tar_(computing)>`_
(`.tar`, `.tar.gz`, `.tar.bz2`, `.tar.xz`), via

.. code-block:: bash

    snakemake --archive my-workflow.tar.gz

If above layout is followed, this will archive any code and config files that
is under git version control. Further, all input files will be included into the
archive. Finally, the software packages of each defined conda environment are included.
This results in a self-contained workflow archive that can be re-executed on a
vanilla machine that only has Conda and Snakemake installed via

.. code-block:: bash

    tar -xf my-workflow.tar.gz
    snakemake -n

Note that the archive is platform specific. For example, if created on Linux, it will
run on any Linux newer than the minimum version that has been supported by the used
Conda packages at the time of archiving (e.g. CentOS 6).

A useful pattern when publishing data analyses is to create such an archive,
upload it to `Zenodo <https://zenodo.org/>`_ and thereby obtain a
`DOI <https://en.wikipedia.org/wiki/Digital_object_identifier>`_.
Then, the DOI can be cited in manuscripts, and readers are able to download
and reproduce the data analysis at any time in the future.
.. _snakefiles-best_practices:

==============
Best practices
==============

* Snakemake (>=5.11) comes with a code quality checker (a so called linter), that analyzes your workflow and highlights issues that should be solved in order to follow best practices, achieve maximum readability, and reproducibility.
  The linter can be invoked with 

  .. code-block:: bash

      snakemake --lint

  given that a ``Snakefile`` or ``workflow/Snakefile`` is accessible from your working directory.
  It is **highly recommended** to run the linter before publishing any workflow, asking questions on Stack Overflow or filing issues on Github.
* There is an automatic formatter for Snakemake workflows, called `Snakefmt <https://github.com/snakemake/snakefmt>`_, which should be applied to any Snakemake workflow before publishing it.
* When publishing your workflow in a `Github <https://github.com>`_ repository, it is a good idea to add some minimal test data and configure `Github Actions <https://github.com/features/actions>`_ for continuously testing the workflow on each new commit.
  For this purpose, we provide predefined Github actions for both running tests and linting `here <https://github.com/snakemake/snakemake-github-action>`_, as well as formatting `here <https://github.com/snakemake/snakefmt#github-actions>`_.
* For publishing and distributing a Snakemake workflow, it is a good idea to stick to a :ref:`standardized structure <distribution_and_reproducibility>` that is expected by frequent users of Snakemake.
  The `Snakemake workflow catalog <https://snakemake.github.io/snakemake-workflow-catalog>`_ automatically lists Snakemake workflows hosted on `Github <https://github.com>`_ if they follow certain `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_.
  By complying to these `rules <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_ you can make your workflow more discoverable and even automate its usage documentation (see `"Standardized usage" <https://snakemake.github.io/snakemake-workflow-catalog/?rules=true>`_).
* Configuration of a workflow should be handled via :ref:`config files <snakefiles_standard_configuration>` and, if needed, tabular configuration like sample sheets (either via :ref:`Pandas <snakefiles_tabular_configuration>` or :ref:`PEPs <snakefiles_peps>`).
  Use such configuration for metadata and experiement information, **not for runtime specific configuration** like threads, resources and output folders.
  For those, just rely on Snakemake's CLI arguments like ``--set-threads``, ``--set-resources``, ``--set-default-resources``, and ``--directory``. 
  This makes workflows more readable, scalable, and portable.
* Try to keep filenames short (thus easier on the eye), but informative. Avoid mixing of too many special characters (e.g. decide whether to use ``_`` or ``-`` as a separator and do that consistently throughout the workflow).
* Try to keep Python code like helper functions separate from rules (e.g. in a ``workflow/rules/common.smk`` file). This way, you help non-experts to read the workflow without needing to parse internals that are irrelevant for them. The helper function names should be chosen in a way that makes them sufficiently informative without looking at their content. Also avoid ``lambda`` expressions inside of rules.
* Make use of `Snakemake wrappers <https://snakemake-wrappers.readthedocs.io>`_ whenever possible. Consider contributing to the wrapper repo whenever you have a rule that reoccurs in at least two of your workflows... _snakefiles-rules:

====================
Snakefiles and Rules
====================

A Snakemake workflow defines a data analysis in terms of rules that are specified in the Snakefile.
Most commonly, rules consist of a name, input files, output files, and a shell command to generate the output from the input:

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        shell: "somecommand {input} {output}"

The name is optional and can be left out, creating an anonymous rule. It can also be overridden by setting a rule's ``name`` attribute.

.. sidebar:: Note

    Note that any placeholders in the shell command (like ``{input}``) are always evaluated and replaced
    when the corresponding job is executed, even if they are occuring inside a comment.
    To avoid evaluation and replacement, you have to mask the braces by doubling them,
    i.e. ``{{input}}``.

Inside the shell command, all local and global variables, especially input and output files can be accessed via their names in the `python format minilanguage <https://docs.python.org/py3k/library/string.html#formatspec>`_. 
Here, input and output (and in general any list or tuple) automatically evaluate to a space-separated list of files (i.e. ``path/to/inputfile path/to/other/inputfile``).
From Snakemake 3.8.0 on, adding the special formatting instruction ``:q`` (e.g. ``"somecommand {input:q} {output:q}")``) will let Snakemake quote each of the list or tuple elements that contains whitespace.


Instead of a shell command, a rule can run some python code to generate the output:

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", somename = "path/to/another/outputfile"
        run:
            for f in input:
                ...
                with open(output[0], "w") as out:
                    out.write(...)
            with open(output.somename, "w") as out:
                out.write(...)

As can be seen, instead of accessing input and output as a whole, we can also access by index (``output[0]``) or by keyword (``output.somename``).
Note that, when adding keywords or names for input or output files, their order won't be preserved when accessing them as a whole via e.g. ``{output}`` in a shell command.

Shell commands like above can also be invoked inside a python based rule, via the function ``shell`` that takes a string with the command and allows the same formatting like in the rule above, e.g.:

.. code-block:: python

    shell("somecommand {output.somename}")

Further, this combination of python and shell commands allows us to iterate over the output of the shell command, e.g.:

.. code-block:: python

    for line in shell("somecommand {output.somename}", iterable=True):
        ... # do something in python

Note that shell commands in Snakemake use the bash shell in `strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_ by default.

.. _snakefiles-wildcards:

Wildcards
---------

Usually, it is useful to generalize a rule to be applicable to a number of e.g. datasets. For this purpose, wildcards can be used.
Automatically resolved multiple named wildcards are a key feature and strength of Snakemake in comparison to other systems.
Consider the following example.

.. code-block:: python

    rule complex_conversion:
        input:
            "{dataset}/inputfile"
        output:
            "{dataset}/file.{group}.txt"
        shell:
            "somecommand --group {wildcards.group} < {input} > {output}"

Here, we define two wildcards, ``dataset`` and ``group``. By this, the rule can produce all files that follow the regular expression pattern ``.+/file\..+\.txt``, i.e. the wildcards are replaced by the regular expression ``.+``. If the rule's output matches a requested file, the substrings matched by the wildcards are propagated to the input files and to the variable wildcards, that is here also used in the shell command. The wildcards object can be accessed in the same way as input and output, which is described above.

For example, if another rule in the workflow requires the file ``101/file.A.txt``, Snakemake recognizes that this rule is able to produce it by setting ``dataset=101`` and ``group=A``.
Thus, it requests file ``101/inputfile`` as input and executes the command ``somecommand --group A  < 101/inputfile  > 101/file.A.txt``.
Of course, the input file might have to be generated by another rule with different wildcards.

Importantly, the wildcard names in input and output must be named identically. Most typically, the same wildcard is present in both input and output, but it is of course also possible to have wildcards only in the output but not the input section.


Multiple wildcards in one filename can cause ambiguity.
Consider the pattern ``{dataset}.{group}.txt`` and assume that a file ``101.B.normal.txt`` is available.
It is not clear whether ``dataset=101.B`` and ``group=normal`` or ``dataset=101`` and ``group=B.normal`` in this case.

Hence wildcards can be constrained to given regular expressions.
Here we could restrict the wildcard ``dataset`` to consist of digits only using ``\d+`` as the corresponding regular expression.
With Snakemake 3.8.0, there are three ways to constrain wildcards.
First, a wildcard can be constrained within the file pattern, by appending a regular expression separated by a comma:

.. code-block:: python

    output: "{dataset,\d+}.{group}.txt"

Second, a wildcard can be constrained within the rule via the keyword ``wildcard_constraints``:

.. code-block:: python

    rule complex_conversion:
        input:
            "{dataset}/inputfile"
        output:
            "{dataset}/file.{group}.txt"
        wildcard_constraints:
            dataset="\d+"
        shell:
            "somecommand --group {wildcards.group}  < {input}  > {output}"

Finally, you can also define global wildcard constraints that apply for all rules:

.. code-block:: python

    wildcard_constraints:
        dataset="\d+"

    rule a:
        ...

    rule b:
        ...

See the `Python documentation on regular expressions <https://docs.python.org/py3k/library/re.html>`_ for detailed information on regular expression syntax.

.. _snakefiles_aggregation:

Aggregation
-----------

Input files can be Python lists, allowing to easily aggregate over parameters or samples:

.. code-block:: python

    rule aggregate:
        input: 
            ["{dataset}/a.txt".format(dataset=dataset) for dataset in DATASETS]
        output:
            "aggregated.txt"
        shell:
            ...

The above expression can be simplified in two ways.

.. _snakefiles_expand:

The expand function
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    rule aggregate:
        input: 
            expand("{dataset}/a.txt", dataset=DATASETS)
        output:
            "aggregated.txt"
        shell:
            ...


Note that *dataset* is NOT a wildcard here because it is resolved by Snakemake due to the ``expand`` statement.
The ``expand`` function also allows us to combine different variables, e.g.

.. code-block:: python

    rule aggregate:
        input: 
            expand("{dataset}/a.{ext}", dataset=DATASETS, ext=FORMATS)
        output:
            "aggregated.txt"
        shell:
            ...

If ``FORMATS=["txt", "csv"]`` contains a list of desired output formats then expand will automatically combine any dataset with any of these extensions.

Furthermore, the first argument can also be a list of strings. In that case, the transformation is applied to all elements of the list. E.g.

.. code-block:: python

    expand(["{dataset}/a.{ext}", "{dataset}/b.{ext}"], dataset=DATASETS, ext=FORMATS)

leads to

.. code-block:: python

    ["ds1/a.txt", "ds1/b.txt", "ds2/a.txt", "ds2/b.txt", "ds1/a.csv", "ds1/b.csv", "ds2/a.csv", "ds2/b.csv"]

Per default, ``expand`` uses the python itertools function ``product`` that yields all combinations of the provided wildcard values. However by inserting a second positional argument this can be replaced by any combinatoric function, e.g. ``zip``:

.. code-block:: python

    expand(["{dataset}/a.{ext}", "{dataset}/b.{ext}"], zip, dataset=DATASETS, ext=FORMATS)

leads to

.. code-block:: python

    ["ds1/a.txt", "ds1/b.txt", "ds2/a.csv", "ds2/b.csv"]

You can also mask a wildcard expression in ``expand`` such that it will be kept, e.g.

.. code-block:: python

    expand("{{dataset}}/a.{ext}", ext=FORMATS)

will create strings with all values for ext but starting with the wildcard ``"{dataset}"``.


.. _snakefiles-multiext:

The multiext function
~~~~~~~~~~~~~~~~~~~~~

``multiext`` provides a simplified variant of ``expand`` that allows us to define a set of output or input files that just differ by their extension:


.. code-block:: python

    rule plot:
        input: 
            ...
        output:
            multiext("some/plot", ".pdf", ".svg", ".png")
        shell:
            ...

The effect is the same as if you would write ``expand("some/plot.{ext}", ext=[".pdf", ".svg", ".png"])``, however, using a simpler syntax.
Moreover, defining output with ``multiext`` is the only way to use :ref:`between workflow caching <caching>` for rules with multiple output files.


.. _snakefiles-targets:

Targets and aggregation
-----------------------

By default snakemake executes the first rule in the snakefile. This gives rise to pseudo-rules at the beginning of the file that can be used to define build-targets similar to GNU Make:

.. code-block:: python

    rule all:
      input:
        expand("{dataset}/file.A.txt", dataset=DATASETS)


Here, for each dataset in a python list ``DATASETS`` defined before, the file ``{dataset}/file.A.txt`` is requested.
In this example, Snakemake recognizes automatically that these can be created by multiple applications of the rule ``complex_conversion`` shown above.

It is possible to overwrite this behavior to use the first rule as a default target, by explicitly marking a rule as being the default target via the ``default_target`` directive:

.. code-block:: python

    rule xy:
        input:
            expand("{dataset}/file.A.txt", dataset=DATASETS)
        default_target: True

Regardless of where this rule appears in the Snakefile, it will be the default target.
Usually, it is still recommended to keep the default target rule (and in fact all other rules that could act as optional targets) at the top of the file, such that it can be easily found.
The ``default_target`` directive becomes particularly useful when :ref:`combining several pre-existing workflows <use_with_modules>`.

.. _snakefiles-threads:

Threads
-------

Further, a rule can be given a number of threads to use, i.e.

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: 8
        shell: "somecommand --threads {threads} {input} {output}"

.. sidebar:: Note

    On a cluster node, Snakemake uses as many cores as available on that node.
    Hence, the number of threads used by a rule never exceeds the number of physically available cores on the node. 
    Note: This behavior is not affected by ``--local-cores``, which only applies to jobs running on the main node.

Snakemake can alter the number of cores available based on command line options. Therefore it is useful to propagate it via the built in variable ``threads`` rather than hardcoding it into the shell command.
In particular, it should be noted that the specified threads have to be seen as a maximum. When Snakemake is executed with fewer cores, the number of threads will be adjusted, i.e. ``threads = min(threads, cores)`` with ``cores`` being the number of cores specified at the command line (option ``--cores``). 

Hardcoding a particular maximum number of threads like above is useful when a certain tool has a natural maximum beyond which parallelization won't help to further speed it up.
This is often the case, and should be evaluated carefully for production workflows.
Also, setting a ``threads:`` maximum is required to achieve parallelism in tools that (often implicitly and without the user knowing) rely on an environment variable for the maximum of cores to use.
For example, this is the case for many linear algebra libraries and for OpenMP.
Snakemake limits the respective environment variables to one core by default, to avoid unexpected and unlimited core-grabbing, but will override this with the ``threads:`` you specify in a rule (the parameters set to ``threads:``, or defaulting to ``1``, are: ``OMP_NUM_THREADS``, ``GOTO_NUM_THREADS``, ``OPENBLAS_NUM_THREADS``, ``MKL_NUM_THREADS``, ``VECLIB_MAXIMUM_THREADS``, ``NUMEXPR_NUM_THREADS``).

If it is certain that no maximum for efficient parallelism exists for a tool, one can instead define threads as a function of the number of cores given to Snakemake:

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: workflow.cores * 0.75
        shell: "somecommand --threads {threads} {input} {output}"

The number of given cores is globally available in the Snakefile as an attribute of the workflow object: ``workflow.cores``.
Any arithmetic operation can be performed to derive a number of threads from this. E.g., in the above example, we reserve 75% of the given cores for the rule.
Snakemake will always round the calculated value down (while enforcing a minimum of 1 thread).

Starting from version 3.7, threads can also be a callable that returns an ``int`` value. The signature of the callable should be ``callable(wildcards[, input])`` (input is an optional parameter).  It is also possible to refer to a predefined variable (e.g, ``threads: threads_max``) so that the number of cores for a set of rules can be changed with one change only by altering the value of the variable ``threads_max``.


.. _snakefiles-resources:

Resources
---------

In addition to threads, a rule can use arbitrary user-defined resources by specifying them with the resources-keyword:

.. code-block:: python

    rule a:
        input:     ...
        output:    ...
        resources:
            mem_mb=100
        shell:
            "..."

If limits for the resources are given via the command line, e.g.

.. code-block:: console

    $ snakemake --resources mem_mb=100


the scheduler will ensure that the given resources are not exceeded by running jobs.
Resources are always meant to be specified as total per job, not by thread (i.e. above ``mem_mb=100`` in rule ``a`` means that any job from rule ``a`` will require ``100`` megabytes of memory in total, and not per thread).

In general, resources are just names to the Snakemake scheduler, i.e., Snakemake does not check whether a job exceeds a certain resource.
However, resources are used to determine which jobs can be executed at a time while not exceeding the given limits at the command line.
If no limits are given, the resources are ignored in local execution.
In cluster or cloud execution, resources are always passed to the backend, even if ``--resources`` is not specified.
Apart from making Snakemake aware of hybrid-computing architectures (e.g. with a limited number of additional devices like GPUs) this allows us to control scheduling in various ways, e.g. to limit IO-heavy jobs by assigning an artificial IO-resource to them and limiting it via the ``--resources`` flag.
Resources must be ``int`` or ``str`` values. Note that you are free to choose any names for the given resources.


Resources can also be callables that return ``int`` or ``str`` values.
The signature of the callable has to be ``callable(wildcards [, input] [, threads] [, attempt])`` (``input``, ``threads``, and ``attempt`` are optional parameters).

The parameter ``attempt`` allows us to adjust resources based on how often the job has been restarted (see :ref:`all_options`, option ``--restart-times``).
This is handy when executing a Snakemake workflow in a cluster environment, where jobs can e.g. fail because of too limited resources.
When Snakemake is executed with ``--restart-times 3``, it will try to restart a failed job 3 times before it gives up.
Thereby, the parameter ``attempt`` will contain the current attempt number (starting from ``1``).
This can be used to adjust the required memory as follows

.. code-block:: python

    def get_mem_mb(wildcards, attempt):
        return attempt * 100

    rule:
        input:    ...
        output:   ...
        resources:
            mem_mb=get_mem_mb
        shell:
            "..."

Here, the first attempt will require 100 MB memory, the second attempt will require 200 MB memory and so on.
When passing memory requirements to the cluster engine, you can by this automatically try out larger nodes if it turns out to be necessary.

Another application of callables as resources is when memory usage depends on the number of threads:

.. code-block:: python

    def get_mem_mb(wildcards, threads):
        return threads * 150

    rule b:
        input:     ...
        output:    ...
        threads: 8
        resources:
            mem_mb=get_mem_mb
        shell:
            "..."

Here, the value the function ``get_mem_mb`` returns grows linearly with the number of threads.
Of course, any other arithmetic could be performed in that function.

Both threads and resources can be overwritten upon invocation via `--set-threads` and `--set-resources`, see :ref:`user_manual-snakemake_options`.

Standard Resources
~~~~~~~~~~~~~~~~~~

There are three **standard resources**, for total memory, disk usage and the temporary directory of a job: ``mem_mb`` and ``disk_mb`` and ``tmpdir``.
The ``tmpdir`` resource automatically leads to setting the TMPDIR variable for shell commands, scripts, wrappers and notebooks.
When defining memory constraints, it is advised to use ``mem_mb``, because some execution modes make direct use of this information (e.g., when using :ref:`Kubernetes <kubernetes>`).

Since it would be cumbersome to define such standard resources them for every rule, you can set default values at 
the terminal or in a :ref:`profile <profiles>`.
This works via the command line flag ``--default-resources``, see ``snakemake --help`` for more information.
If those resource definitions are mandatory for a certain execution mode, Snakemake will fail with a hint if they are missing.
Any resource definitions inside a rule override what has been defined with ``--default-resources``.
If ``--default-resources`` are not specified, Snakemake uses ``'mem_mb=max(2*input.size_mb, 1000)'``, 
``'disk_mb=max(2*input.size_mb, 1000)'``, and ``'tmpdir=system_tmpdir'``.
The latter points to whatever is the default of the operating system or specified by any of the environment variables ``$TMPDIR``, ``$TEMP``, or ``$TMP`` as outlined `here <https://docs.python.org/3/library/tempfile.html#tempfile.gettempdir>`_.


Preemptible Jobs
~~~~~~~~~~~~~~~~


You can specify parameters ``preemptible-rules`` and ``preemption-default`` to request a `Google Cloud preemptible virtual machine <https://cloud.google.com/life-sciences/docs/reference/gcloud-examples#using_preemptible_vms>`_ for use with the `Google Life Sciences Executor <https://snakemake.readthedocs.io/en/stable/executing/cloud.html#executing-a-snakemake-workflow-via-google-cloud-life-sciences>`_. There are
several ways to go about doing this. This first example will use preemptible instances for all rules, with 10 repeats (restarts
of the instance if it stops unexpectedly).

.. code-block:: console

    snakemake --preemption-default 10


If your preference is to set a default but then overwrite some rules with a custom value, this is where you can use ``--preemtible-rules``:

.. code-block:: console

    snakemake --preemption-default 10 --preemptible-rules map_reads=3 call_variants=0


The above statement says that we want to use preemtible instances for all steps, defaulting to 10 retries,
but for the steps "map_reads" and "call_variants" we want to apply 3 and 0 retries, respectively. The final
option is to not use preemptible instances by default, but only for a particular rule:


.. code-block:: console

    snakemake --preemptible-rules map_reads=10


Note that this is currently implemented for the Google Life Sciences API.


GPU Resources
~~~~~~~~~~~~~

The Google Life Sciences API currently has support for 
`NVIDIA GPUs <https://cloud.google.com/compute/docs/gpus#restrictions>`_, meaning that you can request a number of NVIDIA GPUs explicitly by adding ``nvidia_gpu`` or ``gpu`` to your Snakefile resources for a step:


.. code-block:: python

    rule a:
        output:
            "test.txt"
        resources:
            nvidia_gpu=1
        shell:
            "somecommand ..."


A specific `gpu model <https://cloud.google.com/compute/docs/gpus#introduction>`_ can be requested using ``gpu_model`` and lowercase identifiers like ``nvidia-tesla-p100`` or ``nvidia-tesla-p4``, for example: ``gpu_model="nvidia-tesla-p100"``. If you don't specify ``gpu`` or ``nvidia_gpu`` with a count, but you do specify a ``gpu_model``, the count will default to 1.



Messages
--------

When executing snakemake, a short summary for each running rule is given to the console. This can be overridden by specifying a message for a rule:


.. code-block:: python

    rule NAME:
        input: "path/to/inputfile", "path/to/other/inputfile"
        output: "path/to/outputfile", "path/to/another/outputfile"
        threads: 8
        message: "Executing somecommand with {threads} threads on the following files {input}."
        shell: "somecommand --threads {threads} {input} {output}"

Note that access to wildcards is also possible via the variable ``wildcards`` (e.g, ``{wildcards.sample}``), which is the same as with shell commands. It is important to have a namespace around wildcards in order to avoid clashes with other variable names.

Priorities
----------

Snakemake allows for rules that specify numeric priorities:


.. code-block:: python

    rule:
      input: ...
      output: ...
      priority: 50
      shell: ...

Per default, each rule has a priority of 0. Any rule that specifies a higher priority, will be preferred by the scheduler over all rules that are ready to execute at the same time without having at least the same priority.

Furthermore, the ``--prioritize`` or ``-P`` command line flag allows to specify files (or rules) that shall be created with highest priority during the workflow execution. This means that the scheduler will assign the specified target and all its dependencies highest priority, such that the target is finished as soon as possible.
The ``--dry-run`` (equivalently ``--dryrun``) or ``-n`` option allows you to see the scheduling plan including the assigned priorities.



Log-Files
---------

Each rule can specify a log file where information about the execution is written to:

.. code-block:: python

    rule abc:
        input: "input.txt"
        output: "output.txt"
        log: "logs/abc.log"
        shell: "somecommand --log {log} {input} {output}"

Log files can be used as input for other rules, just like any other output file.
However, unlike output files, log files are not deleted upon error.
This is obviously necessary in order to discover causes of errors which might become visible in the log file.

The variable ``log`` can be used inside a shell command to tell the used tool to which file to write the logging information.
The log file has to use the same wildcards as output files, e.g.

.. code-block:: python

    log: "logs/abc.{dataset}.log"


For programs that do not have an explicit ``log`` parameter, you may always use ``2> {log}`` to redirect standard output to a file (here, the ``log`` file) in Linux-based systems.
Note that it is also supported to have multiple (named) log files being specified:

.. code-block:: python

    rule abc:
        input: "input.txt"
        output: "output.txt"
        log: log1="logs/abc.log", log2="logs/xyz.log"
        shell: "somecommand --log {log.log1} METRICS_FILE={log.log2} {input} {output}"

Non-file parameters for rules
-----------------------------

Sometimes you may want to define certain parameters separately from the rule body. Snakemake provides the ``params`` keyword for this purpose:


.. code-block:: python

    rule:
        input:
            ...
        params:
            prefix="somedir/{sample}"
        output:
            "somedir/{sample}.csv"
        shell:
            "somecommand -o {params.prefix}"

The ``params`` keyword allows you to specify additional parameters depending on the wildcards values. This allows you to circumvent the need to use ``run:`` and python code for non-standard commands like in the above case.
Here, the command ``somecommand`` expects the prefix of the output file instead of the actual one. The ``params`` keyword helps here since you cannot simply add the prefix as an output file (as the file won't be created, Snakemake would throw an error after execution of the rule).

Furthermore, for enhanced readability and clarity, the ``params`` section is also an excellent place to name and assign parameters and variables for your subsequent command.

Similar to ``input``, ``params`` can take functions as well (see :ref:`snakefiles-input_functions`), e.g. you can write

.. code-block:: python

    rule:
        input:
            ...
        params:
            prefix=lambda wildcards, output: output[0][:-4]
        output:
            "somedir/{sample}.csv"
        shell:
            "somecommand -o {params.prefix}"

.. sidebar:: Note

    When accessing auxiliary source files (i.e. files that are located relative to the current Snakefile, e.g. some additional configuration)
    it is crucial to not manually build their path but rather rely on Snakemake's special registration for these files, see :ref:`snakefiles-aux_source_files`.

to get the same effect as above. Note that in contrast to the ``input`` directive, the
``params`` directive can optionally take more arguments than only ``wildcards``, namely ``input``, ``output``, ``threads``, and ``resources``.
From the Python perspective, they can be seen as optional keyword arguments without a default value.
Their order does not matter, apart from the fact that ``wildcards`` has to be the first argument.
In the example above, this allows you to derive the prefix name from the output file.

.. _snakefiles-external_scripts:

External scripts
----------------

A rule can also point to an external script instead of a shell command or inline Python code, e.g.

Python
~~~~~~

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "scripts/script.py"

.. sidebar:: Note

    It is possible to refer to wildcards and params in the script path, e.g. by specifying ``"scripts/{params.scriptname}.py"`` or ``"scripts/{wildcards.scriptname}.py"``.

The script path is always relative to the Snakefile containing the directive (in contrast to the input and output file paths, which are relative to the working directory).
It is recommended to put all scripts into a subfolder ``scripts`` as above.
Inside the script, you have access to an object ``snakemake`` that provides access to the same objects that are available in the ``run`` and ``shell`` directives (input, output, params, wildcards, log, threads, resources, config), e.g. you can use ``snakemake.input[0]`` to access the first input file of above rule.

An example external Python script could look like this:

.. code-block:: python

    def do_something(data_path, out_path, threads, myparam):
        # python code

    do_something(snakemake.input[0], snakemake.output[0], snakemake.threads, snakemake.config["myparam"])

You can use the Python debugger from within the script if you invoke Snakemake with ``--debug``.

R and R Markdown
~~~~~~~~~~~~~~~~

Apart from Python scripts, this mechanism also allows you to integrate R_ and R Markdown_ scripts with Snakemake, e.g.

.. _R: https://www.r-project.org
.. _Markdown: https://rmarkdown.rstudio.com

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "scripts/script.R"

In the R script, an S4 object named ``snakemake`` analogous to the Python case above is available and allows access to input and output files and other parameters. Here the syntax follows that of S4 classes with attributes that are R lists, e.g. we can access the first input file with ``snakemake@input[[1]]`` (note that the first file does not have index ``0`` here, because R starts counting from ``1``). Named input and output files can be accessed in the same way, by just providing the name instead of an index, e.g. ``snakemake@input[["myfile"]]``.

An equivalent script (:ref:`to the Python one above <Python>`) written in R would look like this:

.. code-block:: r

    do_something <- function(data_path, out_path, threads, myparam) {
        # R code
    }

    do_something(snakemake@input[[1]], snakemake@output[[1]], snakemake@threads, snakemake@config[["myparam"]])


To debug R scripts, you can save the workspace with ``save.image()``, and invoke R after Snakemake has terminated. Then you can use the usual R debugging facilities while having access to the ``snakemake`` variable.
It is best practice to wrap the actual code into a separate function. This increases the portability if the code shall be invoked outside of Snakemake or from a different rule.
A convenience method, ``snakemake@source()``, acts as a wrapper for the normal R ``source()`` function, and can be used to source files relative to the original script directory.

An R Markdown file can be integrated in the same way as R and Python scripts, but only a single output (html) file can be used:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/report.html",
        script:
            "path/to/report.Rmd"

In the R Markdown file you can insert output from a R command, and access variables stored in the S4 object named ``snakemake``

.. code-block:: R

    ---
    title: "Test Report"
    author:
        - "Your Name"
    date: "`r format(Sys.time(), '%d %B, %Y')`"
    params:
       rmd: "report.Rmd"
    output:
      html_document:
      highlight: tango
      number_sections: no
      theme: default
      toc: yes
      toc_depth: 3
      toc_float:
        collapsed: no
        smooth_scroll: yes
    ---

    ## R Markdown

    This is an R Markdown document.

    Test include from snakemake `r snakemake@input`.

    ## Source
    <a download="report.Rmd" href="`r base64enc::dataURI(file = params$rmd, mime = 'text/rmd', encoding = 'base64')`">R Markdown source file (to produce this document)</a>

A link to the R Markdown document with the snakemake object can be inserted. Therefore a variable called ``rmd`` needs to be added to the ``params`` section in the header of the ``report.Rmd`` file. The generated R Markdown file with snakemake object will be saved in the file specified in this ``rmd`` variable. This file can be embedded into the HTML document using base64 encoding and a link can be inserted as shown in the example above.
Also other input and output files can be embedded in this way to make a portable report. Note that the above method with a data URI only works for small files. An experimental technology to embed larger files is using Javascript Blob `object <https://developer.mozilla.org/en-US/docs/Web/API/Blob>`_.

Julia_
~~~~~~

.. _Julia: https://julialang.org

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile"
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        script:
            "path/to/script.jl"

In the Julia_ script, a ``snakemake`` object is available, which can be accessed similar to the :ref:`Python case <Python>`, with the only difference that you have to index from 1 instead of 0.

Rust_
~~~~~

.. _Rust: https://www.rust-lang.org/

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile",
            "path/to/other/inputfile",
            named_input="path/to/named/inputfile",
        output:
            "path/to/outputfile",
            "path/to/another/outputfile"
        params:
            seed=4
        conda:
            "rust.yaml"
        log:
            stdout="path/to/stdout.log",
            stderr="path/to/stderr.log",
        script:
            "path/to/script.rs"

The ability to execute Rust scripts is facilitated by |rust-script|_.
As such, the script must be a valid ``rust-script`` script and ``rust-script``
(plus OpenSSL and a C compiler toolchain, provided by Conda packages ``openssl``, ``c-compiler``, ``pkg-config``)
must be available in the environment the rule is run in.
The minimum required ``rust-script`` version is 1.15.0, so in the example above, the contents of ``rust.yaml`` might look like this:

.. code block:: yaml

    channels:
      - conda-forge
      - bioconda
    dependencies:
      - rust-script>=0.15.0
      - openssl
      - c-compiler
      - pkg-config



Some example scripts can be found in the
`tests directory <https://github.com/snakemake/snakemake/tree/main/tests/test_script/scripts>`_.

In the Rust script, a ``snakemake`` instance is available, which is automatically generated from the python snakemake object using |json_typegen|_.
It usually looks like this:

.. code-block:: rust

    pub struct Snakemake {
        input: Input,
        output: Ouput,
        params: Params,
        wildcards: Wildcards,
        threads: u64,
        log: Log,
        resources: Resources,
        config: Config,
        rulename: String,
        bench_iteration: Option<usize>,
        scriptdir: String,
    }

Any named parameter is translated to a corresponding ``field_name: Type``, such that ``params.seed`` from the example above can be accessed just like in python, i.e.:

.. code-block:: rust

    let seed = snakemake.params.seed;
    assert_eq!(seed, 4);

Positional arguments for ``input``, ``output``, ``log`` and ``wildcards`` can be accessed by index and iterated over:

.. code-block:: rust

    let input = &snakemake.input;

    // Input implements Index<usize>
    let inputfile = input[0];
    assert_eq!(inputfile, "path/to/inputfile");

    // Input implements IntoIterator
    //
    // prints
    // > 'path/to/inputfile'
    // > 'path/to/other/inputfile'
    for f in input {
        println!("> '{}'", &f);
    }


It is also possible to redirect ``stdout`` and ``stderr``:

.. code-block:: rust

    println!("This will NOT be written to path/to/stdout.log");
    // redirect stdout to "path/to/stdout.log"
    let _stdout_redirect = snakemake.redirect_stdout(snakemake.log.stdout)?;
    println!("This will be written to path/to/stdout.log");

    // redirect stderr to "path/to/stderr.log"
    let _stderr_redirect = snakemake.redirect_stderr(snakemake.log.stderr)?;
    eprintln!("This will be written to path/to/stderr.log");
    drop(_stderr_redirect);
    eprintln!("This will NOT be written to path/to/stderr.log");

Redirection of stdout/stderr is only "active" as long as the returned ``Redirect`` instance is alive; in order to stop redirecting, drop the respective instance.

In order to work, rust-script support for snakemake has some dependencies enabled by default:

#. ``anyhow=1``, for its ``Result`` type
#. ``gag=1``, to enable stdout/stderr redirects
#. ``json_typegen=0.6``, for generating rust structs from a json representation of the snakemake object
#. ``lazy_static=1.4``, to make a ``snakemake`` instance easily accessible
#. ``serde=1``, explicit dependency of ``json_typegen``
#. ``serde_derive=1``, explicit dependency of ``json_typegen``
#. ``serde_json=1``, explicit dependency of ``json_typegen``

If your script uses any of these packages, you do not need to ``use`` them in your script. Trying to ``use`` them will cause a compilation error.

.. |rust-script| replace:: ``rust-script``
.. _rust-script: https://rust-script.org/
.. |json_typegen| replace:: ``json_typegen``
.. _json_typegen: https://github.com/evestera/json_typegen

----

For technical reasons, scripts are executed in ``.snakemake/scripts``. The original script directory is available as ``scriptdir`` in the ``snakemake`` object.

.. _snakefiles_notebook-integration:

Jupyter notebook integration
----------------------------

Instead of plain scripts (see above), one can integrate Jupyter_ Notebooks.
This enables the interactive development of data analysis components (e.g. for plotting).
Integration works as follows (note the use of `notebook:` instead of `script:`):

.. _Jupyter: https://jupyter.org/

.. code-block:: python

    rule hello:
        output:
            "test.txt"
        log:
            # optional path to the processed notebook
            notebook="logs/notebooks/processed_notebook.ipynb"
        notebook:
            "notebooks/hello.py.ipynb"

.. sidebar:: Note

    Consider Jupyter notebook integration as a way to get the best of both worlds.
    A modular, readable workflow definition with Snakemake, and the ability to quickly explore and plot data with Jupyter.
    The benefit will be maximal when integrating many small notebooks that each do a particular job, hence allowing to get away from large monolithic, and therefore unreadable notebooks.

It is recommended to prefix the ``.ipynb`` suffix with either ``.py`` or ``.r`` to indicate the notebook language.
In the notebook, a snakemake object is available, which can be accessed in the same way as the with :ref:`script integration <snakefiles_external-scripts>`.
In other words, you have access to input files via ``snakemake.input`` (in the Python case) and ``snakemake@input`` (in the R case) etc..
Optionally it is possible to automatically store the processed notebook.
This can be achieved by adding a named logfile ``notebook=...`` to the ``log`` directive.

.. sidebar:: Note

    It is possible to refer to wildcards and params in the notebook path, e.g. by specifying ``"notebook/{params.name}.py"`` or ``"notebook/{wildcards.name}.py"``.

In order to simplify the coding of notebooks given the automatically inserted ``snakemake`` object, Snakemake provides an interactive edit mode for notebook rules.
Let us assume you have written above rule, but the notebook does not yet exist.
By running

.. code-block:: console

    snakemake --cores 1 --edit-notebook test.txt

you instruct Snakemake to allow interactive editing of the notebook needed to create the file ``test.txt``.
Snakemake will run all dependencies of the notebook rule, such that all input files are present.
Then, it will start a jupyter notebook server with an empty draft of the notebook, in which you can interactively program everything needed for this particular step.
Once done, you should save the notebook from the jupyter web interface, go to the jupyter dashboard and hit the ``Quit`` button on the top right in order to shut down the jupyter server.
Snakemake will detect that the server is closed and automatically store the drafted notebook into the path given in the rule (here ``hello.py.ipynb``).
If the notebook already exists, above procedure can be used to easily modify it.
Note that Snakemake requires local execution for the notebook edit mode.
On a cluster or the cloud, you can generate all dependencies of the notebook rule via

.. code-block:: console

    snakemake --cluster ... --jobs 100 --until test.txt

Then, the notebook rule can easily be executed locally.
An demo of the entire interactive editing process can be found by clicking below:

.. image:: images/snakemake-notebook-demo.gif
    :scale: 20%
    :alt: Notebook integration demo
    :align: center

Finally, it is advisable to combine the ``notebook`` directive with the ``conda`` directive (see :ref:`integrated_package_management`) in order to define a software stack to use.
At least, this software stack should contain jupyter and the language to use (e.g. Python or R).
For the above case, this means

.. code-block:: python

    rule hello:
        output:
            "test.txt"
        conda:
            "envs/hello.yaml"
        notebook:
            "notebooks/hello.py.ipynb"

with

.. code-block:: yaml

    channels:
      - conda-forge
    dependencies:
      - python =3.8
      - jupyter =1.0
      - jupyterlab_code_formatter =1.4

The last dependency is advisable in order to enable autoformatting of notebook cells when editing.
When using other languages than Python in the notebook, one needs to additionally add the respective kernel, e.g. ``r-irkernel`` for R support.

When using an IDE with built-in Jupyter support, an alternative to ``--edit-notebook`` is ``--draft-notebook``.
Instead of firing up a notebook server, ``--draft-notebook`` just creates a skeleton notebook for editing within the IDE.
In addition, it prints instructions for configuring the IDE's notebook environment to use the interpreter from the 
Conda environment defined in the corresponding rule.
For example, running

.. code-block:: console

    snakemake --cores 1 --draft-notebook test.txt --use-conda

will generate skeleton code in ``notebooks/hello.py.ipynb`` and additionally print instructions on how to open and execute the notebook in VSCode.


Protected and Temporary Files
-----------------------------

A particular output file may require a huge amount of computation time. Hence one might want to protect it against accidental deletion or overwriting. Snakemake allows this by marking such a file as ``protected``:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            protected("path/to/outputfile")
        shell:
            "somecommand {input} {output}"

A protected file will be write-protected after the rule that produces it is completed.

Further, an output file marked as ``temp`` is deleted after all rules that use it as an input are completed:

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            temp("path/to/outputfile")
        shell:
            "somecommand {input} {output}"

Directories as outputs
----------------------

Sometimes it can be convenient to have directories, rather than files, as outputs of a rule. As of version 5.2.0, directories as outputs have to be explicitly marked with ``directory``. This is primarily for safety reasons; since all outputs are deleted before a job is executed, we don't want to risk deleting important directories if the user makes some mistake. Marking the output as ``directory`` makes the intent clear, and the output can be safely removed. Another reason comes down to how modification time for directories work. The modification time on a directory changes when a file or a subdirectory is added, removed or renamed. This can easily happen in not-quite-intended ways, such as when Apple macOS or MS Windows add ``.DS_Store`` or ``thumbs.db`` files to store parameters for how the directory contents should be displayed. When the ``directory`` flag is used a hidden file called ``.snakemake_timestamp`` is created in the output directory, and the modification time of that file is used when determining whether the rule output is up to date or if it needs to be rerun. Always consider if you can't formulate your workflow using normal files before resorting to using ``directory()``.

.. code-block:: python

    rule NAME:
        input:
            "path/to/inputfile"
        output:
            directory("path/to/outputdir")
        shell:
            "somecommand {input} {output}"

Ignoring timestamps
-------------------

For determining whether output files have to be re-created, Snakemake checks whether the file modification date (i.e. the timestamp) of any input file of the same job is newer than the timestamp of the output file.
This behavior can be overridden by marking an input file as ``ancient``.
The timestamp of such files is ignored and always assumed to be older than any of the output files:

.. code-block:: python

    rule NAME:
        input:
            ancient("path/to/inputfile")
        output:
            "path/to/outputfile"
        shell:
            "somecommand {input} {output}"

Here, this means that the file ``path/to/outputfile`` will not be triggered for re-creation after it has been generated once, even when the input file is modified in the future.
Note that any flag that forces re-creation of files still also applies to files marked as ``ancient``.

Shadow rules
------------

Shadow rules result in each execution of the rule to be run in isolated temporary directories.
This "shadow" directory contains symlinks to files and directories in the current workdir.
This is useful for running programs that generate lots of unused files which you don't want to manually cleanup in your snakemake workflow.
It can also be useful if you want to keep your workdir clean while the program executes,
or simplify your workflow by not having to worry about unique filenames for all outputs of all rules.

By setting ``shadow: "shallow"``, the top level files and directories are symlinked,
so that any relative paths in a subdirectory will be real paths in the filesystem.
The setting ``shadow: "full"`` fully shadows the entire subdirectory structure of the current workdir.
The setting ``shadow: "minimal"`` only symlinks the inputs to the rule,
and ``shadow: "copy-minimal"`` copies the inputs instead of just creating symlinks.
Once the rule successfully executes, the output file will be moved if necessary to the real path as indicated by ``output``.

Typically, you will not need to modify your rule for compatibility with ``shadow``,
unless you reference parent directories relative to your workdir in a rule.

.. code-block:: python

    rule NAME:
        input: "path/to/inputfile"
        output: "path/to/outputfile"
        shadow: "shallow"
        shell: "somecommand --other_outputs other.txt {input} {output}"

Shadow directories are stored one per rule execution in ``.snakemake/shadow/``,
and are cleared on successful execution.
Consider running with the ``--cleanup-shadow`` argument every now and then
to remove any remaining shadow directories from aborted jobs.
The base shadow directory can be changed with the ``--shadow-prefix`` command line argument.

Flag files
----------

Sometimes it is necessary to enforce some rule execution order without real file dependencies. This can be achieved by "touching" empty files that denote that a certain task was completed. Snakemake supports this via the `touch` flag:

.. code-block:: python

    rule all:
        input: "mytask.done"

    rule mytask:
        output: touch("mytask.done")
        shell: "mycommand ..."

With the ``touch`` flag, Snakemake touches (i.e. creates or updates) the file ``mytask.done`` after ``mycommand`` has finished successfully.


.. _snakefiles-job_properties:

Job Properties
--------------

When executing a workflow on a cluster using the ``--cluster`` parameter (see below), Snakemake creates a job script for each job to execute.
This script is then invoked using the provided cluster submission command (e.g. ``qsub``).
Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters.
As this might be based on properties of the job, Snakemake stores the job properties (e.g. rule name, threads, input files, params etc.) as JSON inside the job script.
For convenience, there exists a parser function ``snakemake.utils.read_job_properties`` that can be used to access the properties.
The following shows an example job submission wrapper:

.. code-block:: python

    #!/usr/bin/env python3
    import os
    import sys

    from snakemake.utils import read_job_properties

    jobscript = sys.argv[1]
    job_properties = read_job_properties(jobscript)

    # do something useful with the threads
    threads = job_properties[threads]

    # access property defined in the cluster configuration file (Snakemake >=3.6.0)
    job_properties["cluster"]["time"]

    os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))

.. _snakefiles-input_functions:

Functions as Input Files
------------------------

Instead of specifying strings or lists of strings as input files, snakemake can also make use of functions that return single **or** lists of input files:

.. code-block:: python

    def myfunc(wildcards):
        return [... a list of input files depending on given wildcards ...]

    rule:
        input: myfunc
        output: "someoutput.{somewildcard}.txt"
        shell: "..."

The function has to accept a single argument that will be the wildcards object generated from the application of the rule to create some requested output files.
Note that you can also use `lambda expressions <https://docs.python.org/3/tutorial/controlflow.html#lambda-expressions>`_ instead of full function definitions.
By this, rules can have entirely different input files (both in form and number) depending on the inferred wildcards. E.g. you can assign input files that appear in entirely different parts of your filesystem based on some wildcard value and a dictionary that maps the wildcard value to file paths.

Note that the function will be executed when the rule is evaluated and before the workflow actually starts to execute. Further note that using a function as input overrides the default mechanism of replacing wildcards with their values inferred from the output files. You have to take care of that yourself with the given wildcards object.

Finally, when implementing the input function, it is best practice to make sure that it can properly handle all possible wildcard values your rule can have.
In particular, input files should not be combined with very general rules that can be applied to create almost any file: Snakemake will try to apply the rule, and will report the exceptions of your input function as errors.

For a practical example, see the :ref:`tutorial` (:ref:`tutorial-input_functions`).

.. _snakefiles-unpack:

Input Functions and ``unpack()``
--------------------------------

In some cases, you might want to have your input functions return named input files.
This can be done by having them return ``dict()`` objects with the names as the dict keys and the file names as the dict values and using the ``unpack()`` keyword.

.. code-block:: python

    def myfunc(wildcards):
        return {'foo': '{wildcards.token}.txt'.format(wildcards=wildcards)}

    rule:
        input: unpack(myfunc)
        output: "someoutput.{token}.txt"
        shell: "..."

Note that ``unpack()`` is only necessary for input functions returning ``dict``.
While it also works for ``list``, remember that lists (and nested lists) of strings are automatically flattened.

Also note that if you do not pass in a *function* into the input list but you directly *call a function* then you shouldn't use ``unpack()``.
Here, you can simply use Python's double-star (``**``) operator for unpacking the parameters.

Note that as Snakefiles are translated into Python for execution, the same rules as for using the `star and double-star unpacking Python operators <https://docs.python.org/3/tutorial/controlflow.html#unpacking-argument-lists>`_ apply.
These restrictions do not apply when using ``unpack()``.

.. code-block:: python

    def myfunc1():
        return ['foo.txt']

    def myfunc2():
        return {'foo': 'nowildcards.txt'}

    rule:
        input:
            *myfunc1(),
            **myfunc2(),
        output: "..."
        shell: "..."

.. _snakefiles-version_tracking:

Version Tracking
----------------

Rules can specify a version that is tracked by Snakemake together with the output files. When the version changes snakemake informs you when using the flag ``--summary`` or ``--list-version-changes``.
The version can be specified by the version directive, which takes a string:

.. code-block:: python

    rule:
        input:   ...
        output:  ...
        version: "1.0"
        shell:   ...

The version can of course also be filled with the output of a shell command, e.g.:

.. code-block:: python

    SOMECOMMAND_VERSION = subprocess.check_output("somecommand --version", shell=True)

    rule:
        version: SOMECOMMAND_VERSION

Alternatively, you might want to use file modification times in case of local scripts:

.. code-block:: python

    SOMECOMMAND_VERSION = str(os.path.getmtime("path/to/somescript"))

    rule:
        version: SOMECOMMAND_VERSION

A re-run can be automated by invoking Snakemake as follows:

.. code-block:: console

    $ snakemake -R `snakemake --list-version-changes`

With the availability of the ``conda`` directive (see :ref:`integrated_package_management`)
the ``version`` directive has become **obsolete** in favor of defining isolated
software environments that can be automatically deployed via the conda package
manager.


.. _snakefiles-code_tracking:

Code Tracking
-------------

Snakemake tracks the code that was used to create your files.
In combination with ``--summary`` or ``--list-code-changes`` this can be used to see what files may need a re-run because the implementation changed.
Re-run can be automated by invoking Snakemake as follows:

.. code-block:: console

    $ snakemake -R `snakemake --list-code-changes`


.. _snakefiles-job_lifetime_handlers:

Onstart, onsuccess and onerror handlers
---------------------------------------

Sometimes, it is necessary to specify code that shall be executed when the workflow execution is finished (e.g. cleanup, or notification of the user).
With Snakemake 3.2.1, this is possible via the ``onsuccess`` and ``onerror`` keywords:

.. code-block:: python

    onsuccess:
        print("Workflow finished, no error")

    onerror:
        print("An error occurred")
        shell("mail -s "an error occurred" youremail@provider.com < {log}")

The ``onsuccess`` handler is executed if the workflow finished without error. Otherwise, the ``onerror`` handler is executed.
In both handlers, you have access to the variable ``log``, which contains the path to a logfile with the complete Snakemake output.
Snakemake 3.6.0 adds an ``onstart`` handler, that will be executed before the workflow starts.
Note that dry-runs do not trigger any of the handlers.


Rule dependencies
-----------------

From version 2.4.8 on, rules can also refer to the output of other rules in the Snakefile, e.g.:

.. code-block:: python

    rule a:
        input:  "path/to/input"
        output: "path/to/output"
        shell:  ...

    rule b:
        input:  rules.a.output
        output: "path/to/output/of/b"
        shell:  ...

Importantly, be aware that referring to rule ``a`` here requires that rule ``a`` was defined above rule ``b`` in the file, since the object has to be known already.
This feature also allows us to resolve dependencies that are ambiguous when using filenames.

Note that when the rule you refer to defines multiple output files but you want to require only a subset of those as input for another rule, you should name the output files and refer to them specifically:

.. code-block:: python

    rule a:
        input:  "path/to/input"
        output: a = "path/to/output", b = "path/to/output2"
        shell:  ...

    rule b:
        input:  rules.a.output.a
        output: "path/to/output/of/b"
        shell:  ...


.. _snakefiles-ambiguous-rules:

Handling Ambiguous Rules
------------------------

When two rules can produce the same output file, snakemake cannot decide which one to use without additional guidance. Hence an ``AmbiguousRuleException`` is thrown.
Note: ruleorder is not intended to bring rules in the correct execution order (this is solely guided by the names of input and output files you use), it only helps snakemake to decide which rule to use when multiple ones can create the same output file!
To deal with such ambiguity, provide a ``ruleorder`` for the conflicting rules, e.g.

.. code-block:: python

    ruleorder: rule1 > rule2 > rule3

Here, ``rule1`` is preferred over ``rule2`` and ``rule3``, and ``rule2`` is preferred over ``rule3``.
Only if rule1 and rule2 cannot be applied (e.g. due to missing input files), rule3 is used to produce the desired output file.

Alternatively, rule dependencies (see above) can also resolve ambiguities.

Another (quick and dirty) possiblity is to tell snakemake to allow ambiguity via a command line option

.. code-block:: console

    $ snakemake --allow-ambiguity

such that similar to GNU Make always the first matching rule is used. Here, a warning that summarizes the decision of snakemake is provided at the terminal.

.. _snakefiles-local-rules:

Local Rules
-----------

When working in a cluster environment, not all rules need to become a job that has to be submitted (e.g. downloading some file, or a target-rule like `all`, see :ref:`snakefiles-targets`).
The keyword `localrules` allows to mark a rule as local, so that it is not submitted to the cluster and instead executed on the host node:

.. code-block:: python

    localrules: all, foo

    rule all:
        input: ...

    rule foo:
        ...

    rule bar:
        ...

Here, only jobs from the rule ``bar`` will be submitted to the cluster, whereas all and foo will be run locally.
Note that you can use the localrules directive **multiple times**. The result will be the union of all declarations.

Benchmark Rules
---------------

Since version 3.1, Snakemake provides support for benchmarking the run times of rules.
This can be used to create complex performance analysis pipelines.
With the `benchmark` keyword, a rule can be declared to store a benchmark of its code into the specified location. E.g. the rule

.. code-block:: python

    rule benchmark_command:
        input:
            "path/to/input.{sample}.txt"
        output:
            "path/to/output.{sample}.txt"
        benchmark:
            "benchmarks/somecommand/{sample}.tsv"
        shell:
            "somecommand {input} {output}"

benchmarks the CPU and wall clock time of the command ``somecommand`` for the given output and input files.
For this, the shell or run body of the rule is executed on that data, and all run times are stored into the given benchmark tsv file (which will contain a tab-separated table of run times and memory usage in MiB).
Per default, Snakemake executes the job once, generating one run time.
However, the benchmark file can be annotated with the desired number of repeats, e.g.,

.. code-block:: python

    rule benchmark_command:
        input:
            "path/to/input.{sample}.txt"
        output:
            "path/to/output.{sample}.txt"
        benchmark:
            repeat("benchmarks/somecommand/{sample}.tsv", 3)
        shell:
            "somecommand {input} {output}"

will instruct Snakemake to run each job of this rule three times and store all measurements in the benchmark file.
The resulting tsv file can be used as input for other rules, just like any other output file.

.. sidebar:: Note

    Note that benchmarking is only possible in a reliable fashion for subprocesses (thus for tasks run through the ``shell``, ``script``, and ``wrapper`` directive).
    In the ``run`` block, the variable ``bench_record`` is available that you can pass to ``shell()`` as ``bench_record=bench_record``.
    When using ``shell(..., bench_record=bench_record)``, the maximum of all measurements of all ``shell()`` calls will be used but the running time of the rule execution including any Python code.


.. _snakefiles-scattergather:

Defining scatter-gather processes
---------------------------------

Via Snakemake's powerful and abitrary Python based aggregation abilities (via the ``expand`` function and arbitrary Python code, see :ref:`here <snakefiles_aggregation>`), scatter-gather workflows well supported.
Nevertheless, it can sometimes be handy to use Snakemake's specific scatter-gather support, which allows to avoid boilerplate and offers additional configuration options.
Scatter-gather processes can be defined via a global ``scattergather`` directive:

.. code-block:: python

    scattergather:
        split=8

Each process thereby defines a name (here e.g. ``split``) and a default number of scatter items.
Then, scattering and gathering can be implemented by using globally available ``scatter`` and ``gather`` objects:

.. code-block:: python


    rule all:
        input:
            "gathered/all.txt"


    rule split:
        output:
            scatter.split("splitted/{scatteritem}.txt")
        shell:
            "touch {output}"


    rule intermediate:
        input:
            "splitted/{scatteritem}.txt"
        output:
            "splitted/{scatteritem}.post.txt"
        shell:
            "cp {input} {output}"


    rule gather:
        input:
            gather.split("splitted/{scatteritem}.post.txt")
        output:
            "gathered/all.txt"
        shell:
            "cat {input} > {output}"

Thereby, ``scatter.split("splitted/{scatteritem}.txt")`` yields a list of paths ``"splitted/1-of-n.txt"``, ``"splitted/2-of-n.txt"``, ..., depending on the number ``n`` of scatter items defined.
Analogously, ``gather.split("splitted/{scatteritem}.post.txt")``, yields a list of paths ``"splitted/0.post.txt"``, ``"splitted/1.pos.txt"``, ..., which request the application of the rule ``intermediate`` to each scatter item.

The default number of scatter items can be overwritten via the command line interface.
For example

.. code-block:: bash

    snakemake --set-scatter split=2

would set the number of scatter items for the split process defined above to 2 instead of 8. 
This allows to adapt parallelization according to the needs of the underlying computing platform and the analysis at hand.

.. _snakefiles-grouping:

Defining groups for execution
-----------------------------

From Snakemake 5.0 on, it is possible to assign rules to groups.
Such groups will be executed together in **cluster** or **cloud mode**, as a so-called **group job**, i.e., all jobs of a particular group will be submitted at once, to the same computing node.
When executing locally, group definitions are ignored.

Groups can be defined via the ``group`` keyword.
This way, queueing and execution time can be saved, in particular if one or several short-running rules are involved.

.. code-block:: python

  samples = [1,2,3,4,5]


  rule all:
      input:
          "test.out"


  rule a:
      output:
          "a/{sample}.out"
      group: "mygroup"
      shell:
          "touch {output}"


  rule b:
      input:
          "a/{sample}.out"
      output:
          "b/{sample}.out"
      group: "mygroup"
      shell:
          "touch {output}"


  rule c:
      input:
          expand("b/{sample}.out", sample=samples)
      output:
          "test.out"
      shell:
          "touch {output}"

Here, jobs from rule ``a`` and ``b`` end up in one group ``mygroup``, whereas jobs from rule ``c`` are executed separately.
Note that Snakemake always determines a **connected subgraph** with the same group id to be a **group job**.
Here, this means that, e.g., the jobs creating ``a/1.out`` and ``b/1.out`` will be in one group, and the jobs creating ``a/2.out`` and ``b/2.out`` will be in a separate group.
However, if we would add ``group: "mygroup"`` to rule ``c``, all jobs would end up in a single group, including the one spawned from rule ``c``, because ``c`` connects all the other jobs.

Alternatively, groups can be defined via the command line interface.
This enables to almost arbitrarily partition the DAG, e.g. in order to safe network traffic, see :ref:`here <job_grouping>`.

For execution on the cloud using Google Life Science API and preemptible instances, we expect all rules in the group to be homogenously set as preemptible instances (e.g., with command-line option ``--preemptible-rules``), such that a preemptible VM is requested for the execution of the group job.

Piped output
------------

From Snakemake 5.0 on, it is possible to mark output files as pipes, via the ``pipe`` flag, e.g.:

.. code-block:: python

  rule all:
      input:
          expand("test.{i}.out", i=range(2))


  rule a:
      output:
          pipe("test.{i}.txt")
      shell:
          "for i in {{0..2}}; do echo {wildcards.i} >> {output}; done"


  rule b:
      input:
          "test.{i}.txt"
      output:
          "test.{i}.out"
      shell:
          "grep {wildcards.i} < {input} > {output}"

If an output file is marked to be a pipe, then Snakemake will first create a `named pipe <https://en.wikipedia.org/wiki/Named_pipe>`_ with the given name and then execute the creating job simultaneously with the consuming job, inside a **group job** (see above).
This works in all execution modes, local, cluster, and cloud.
Naturally, a pipe output may only have a single consumer.
It is possible to combine explicit group definition as above with pipe outputs.
Thereby, pipe jobs can live within, or (automatically) extend existing groups.
However, the two jobs connected by a pipe may not exist in conflicting groups.

.. _snakefiles-paramspace:

Parameter space exploration
---------------------------

The basic Snakemake functionality already provides everything to handle parameter spaces in any way (sub-spacing for certain rules and even depending on wildcard values, the ability to read or generate spaces on the fly or from files via pandas, etc.).
However, it usually would require some boilerplate code for translating a parameter space into wildcard patterns, and translate it back into concrete parameters for scripts and commands. 
From Snakemake 5.31 on (inspired by `JUDI <https://pyjudi.readthedocs.io>`_), this is solved via the Paramspace helper, which can be used as follows:

.. code-block:: python

    from snakemake.utils import Paramspace
    import pandas as pd

    # declare a dataframe to be a paramspace
    paramspace = Paramspace(pd.read_csv("params.tsv", sep="\t"))


    rule all:
        input:
            # Aggregate over entire parameter space (or a subset thereof if needed)
            # of course, something like this can happen anywhere in the workflow (not 
            # only at the end).
            expand("results/plots/{params}.pdf", params=paramspace.instance_patterns)


    rule simulate:
        output:
            # format a wildcard pattern like "alpha~{alpha}/beta~{beta}/gamma~{gamma}" 
            # into a file path, with alpha, beta, gamma being the columns of the data frame
            f"results/simulations/{paramspace.wildcard_pattern}.tsv"
        params:
            # automatically translate the wildcard values into an instance of the param space
            # in the form of a dict (here: {"alpha": ..., "beta": ..., "gamma": ...})
            simulation=paramspace.instance
        script:
            "scripts/simulate.py"


    rule plot:
        input:
            f"results/simulations/{paramspace.wildcard_pattern}.tsv"
        output:
            f"results/plots/{paramspace.wildcard_pattern}.pdf"
        shell:
            "touch {output}"


In above example, **please note** the Python ``f``-string formatting (the ``f`` before the initial quotes) applied to the input and output file strings that contain ``paramspace.wildcard_pattern``.
This means that the file that is registered as input or output file by Snakemake does not contain a wildcard ``{paramspace.wildcard_pattern}``, but instead this item is replaced by a pattern of multiple wildcards derived from the columns of the paramter space dataframe.
This is done by the Python ``f``-string formatting before the string is registered in the rule.
Given that `params.tsv` contains:

.. code-block:: none

    alpha	beta	gamma
    1.0	0.1	0.99
    2.0	0.0	3.9


This workflow will run as follows:

.. code-block:: none

    [Fri Nov 27 20:57:27 2020]
    rule simulate:
        output: results/simulations/alpha~2.0/beta~0.0/gamma~3.9.tsv                                                                                                                           
        jobid: 4                                                                                                                                                                               
        wildcards: alpha=2.0, beta=0.0, gamma=3.9                                                                                                                                              

    [Fri Nov 27 20:57:27 2020]
    rule simulate:
        output: results/simulations/alpha~1.0/beta~0.1/gamma~0.99.tsv                                                                                                                          
        jobid: 2                                                                                                                                                                               
        wildcards: alpha=1.0, beta=0.1, gamma=0.99                                                                                                                                             

    [Fri Nov 27 20:57:27 2020]
    rule plot:
        input: results/simulations/alpha~2.0/beta~0.0/gamma~3.9.tsv                                                                                                                            
        output: results/plots/alpha~2.0/beta~0.0/gamma~3.9.pdf                                                                                                                                 
        jobid: 3                                                                                                                                                                               
        wildcards: alpha=2.0, beta=0.0, gamma=3.9                                                                                                                                              


    [Fri Nov 27 20:57:27 2020]
    rule plot:
        input: results/simulations/alpha~1.0/beta~0.1/gamma~0.99.tsv                                                                                                                           
        output: results/plots/alpha~1.0/beta~0.1/gamma~0.99.pdf                                                                                                                                
        jobid: 1                                                                                                                                                                               
        wildcards: alpha=1.0, beta=0.1, gamma=0.99                                                                                                                                             


    [Fri Nov 27 20:57:27 2020]
    localrule all:
        input: results/plots/alpha~1.0/beta~0.1/gamma~0.99.pdf, results/plots/alpha~2.0/beta~0.0/gamma~3.9.pdf                                                                                 
        jobid: 0


Naturally, it is possible to create sub-spaces from ``Paramspace`` objects, simply by applying all the usual methods and attributes that Pandas data frames provide (e.g. ``.loc[...]``, ``.filter()`` etc.).
Further, the form of the created ``wildcard_pattern`` can be controlled via additional arguments of the ``Paramspace`` constructor (see :ref:`utils-api`).

.. _snakefiles-checkpoints:

Data-dependent conditional execution
------------------------------------

From Snakemake 5.4 on, conditional reevaluation of the DAG of jobs based on the content outputs is possible.
The key idea is that rules can be declared as checkpoints, e.g.,

.. code-block:: python

  checkpoint somestep:
      input:
          "samples/{sample}.txt"
      output:
          "somestep/{sample}.txt"
      shell:
          "somecommand {input} > {output}"

Snakemake allows to re-evaluate the DAG after the successful execution of every job spawned from a checkpoint.
For this, every checkpoint is registered by its name in a globally available ``checkpoints`` object.
The ``checkpoints`` object can be accessed by :ref:`input functions <snakefiles-input_functions>`.
Assuming that the checkpoint is named ``somestep`` as above, the output files for a particular job can be retrieved with

.. code-block:: python

  checkpoints.somestep.get(sample="a").output

.. sidebar:: Note

    Note that output files of checkpoints that are accessed via this mechanism should not be marked as temporary.
    Otherwise, they would require to trigger reruns of the checkpoint whenever the DAG shall be reevaluated (because they are already missing at that point).

Thereby, the ``get`` method throws ``snakemake.exceptions.IncompleteCheckpointException`` if the checkpoint has not yet been executed for these particular wildcard value(s).
Inside an input function, the exception will be automatically handled by Snakemake, and leads to a re-evaluation after the checkpoint has been successfully passed.

To illustrate the possibilities of this mechanism, consider the following complete example:

.. code-block:: python

  # a target rule to define the desired final output
  rule all:
      input:
          "aggregated/a.txt",
          "aggregated/b.txt"


  # the checkpoint that shall trigger re-evaluation of the DAG
  checkpoint somestep:
      input:
          "samples/{sample}.txt"
      output:
          "somestep/{sample}.txt"
      shell:
          # simulate some output value
          "echo {wildcards.sample} > somestep/{wildcards.sample}.txt"


  # intermediate rule
  rule intermediate:
      input:
          "somestep/{sample}.txt"
      output:
          "post/{sample}.txt"
      shell:
          "touch {output}"


  # alternative intermediate rule
  rule alt_intermediate:
      input:
          "somestep/{sample}.txt"
      output:
          "alt/{sample}.txt"
      shell:
          "touch {output}"


  # input function for the rule aggregate
  def aggregate_input(wildcards):
      # decision based on content of output file
      # Important: use the method open() of the returned file!
      # This way, Snakemake is able to automatically download the file if it is generated in
      # a cloud environment without a shared filesystem.
      with checkpoints.somestep.get(sample=wildcards.sample).output[0].open() as f:
          if f.read().strip() == "a":
              return "post/{sample}.txt"
          else:
              return "alt/{sample}.txt"


  rule aggregate:
      input:
          aggregate_input
      output:
          "aggregated/{sample}.txt"
      shell:
          "touch {output}"

As can be seen, the rule aggregate uses an input function.
Inside the function, we first retrieve the output files of the checkpoint ``somestep`` with the wildcards, passing through the value of the wildcard sample.
Upon execution, if the checkpoint is not yet complete, Snakemake will record ``somestep`` as a direct dependency of the rule ``aggregate``.
Once ``somestep`` has finished for a given sample, the input function will automatically be re-evaluated and the method ``get`` will no longer return an exception.
Instead, the output file will be opened, and depending on its contents either ``"post/{sample}.txt"`` or ``"alt/{sample}.txt"`` will be returned by the input function.
This way, the DAG becomes conditional on some produced data.

It is also possible to use checkpoints for cases where the output files are unknown before execution.
A typical example is a clustering process with an unknown number of clusters, where each cluster shall be saved into a separate file.
Consider the following example:

.. code-block:: python

  # a target rule to define the desired final output
  rule all:
      input:
          "aggregated/a.txt",
          "aggregated/b.txt"


  # the checkpoint that shall trigger re-evaluation of the DAG
  checkpoint clustering:
      input:
          "samples/{sample}.txt"
      output:
          clusters=directory("clustering/{sample}")
      shell:
          "mkdir clustering/{wildcards.sample}; "
          "for i in 1 2 3; do echo $i > clustering/{wildcards.sample}/$i.txt; done"


  # an intermediate rule
  rule intermediate:
      input:
          "clustering/{sample}/{i}.txt"
      output:
          "post/{sample}/{i}.txt"
      shell:
          "cp {input} {output}"


  def aggregate_input(wildcards):
      checkpoint_output = checkpoints.clustering.get(**wildcards).output[0]
      return expand("post/{sample}/{i}.txt",
             sample=wildcards.sample,
             i=glob_wildcards(os.path.join(checkpoint_output, "{i}.txt")).i)


  # an aggregation over all produced clusters
  rule aggregate:
      input:
          aggregate_input
      output:
          "aggregated/{sample}.txt"
      shell:
          "cat {input} > {output}"

Here, our checkpoint simulates a clustering.
We pretend that the number of clusters is unknown beforehand.
Hence, the checkpoint only defines an output ``directory``.
The rule ``aggregate`` again uses the ``checkpoints`` object to retrieve the output of the checkpoint.
This time, instead of explicitly writing

.. code-block:: python

  checkpoints.clustering.get(sample=wildcards.sample).output[0]

we use the shorthand

.. code-block:: python

  checkpoints.clustering.get(**wildcards).output[0]

which automatically unpacks the wildcards as keyword arguments (this is standard python argument unpacking).
If the checkpoint has not yet been executed, accessing ``checkpoints.clustering.get(**wildcards)`` ensure that Snakemake records the checkpoint as a direct dependency of the rule ``aggregate``.
Upon completion of the checkpoint, the input function is re-evaluated, and the code beyond its first line is executed.
Here, we retrieve the values of the wildcard ``i`` based on all files named ``{i}.txt`` in the output directory of the checkpoint.
These values are then used to expand the pattern ``"post/{sample}/{i}.txt"``, such that the rule ``intermediate`` is executed for each of the determined clusters.


.. _snakefiles-rule-inheritance:

Rule inheritance
----------------

With Snakemake 6.0 and later, it is possible to inherit from previously defined rules, or in other words, reuse an existing rule in a modified way.
This works via the ``use rule`` statement that also allows to declare the usage of rules from external modules (see :ref:`snakefiles-modules`).
Consider the following example:

.. code-block:: python

    rule a:
        output:
            "test.out"
        shell:
            "echo test > {output}"


    use rule a as b with:
        output:
            "test2.out"


As can be seen, we first declare a rule a, and then we reuse the rule a as rule b, while changing only the output file and keeping everything else the same.
In reality, one will often change more.
Analogously to the ``use rule`` from external modules, any properties of the rule (``input``, ``output``, ``log``, ``params``, ``benchmark``, ``threads``, ``resources``, etc.) can be modified, except the actual execution step (``shell``, ``notebook``, ``script``, ``cwl``, or ``run``).
All unmodified properties are inherited from the parent rule.

.. _snakefiles-aux_source_files:

Accessing auxiliary source files
--------------------------------

Snakemake workflows can refer to various other source files via paths relative to the current Snakefile.
This happens for example with the :ref:`script directive <snakefiles-external_scripts>` or the :ref:`conda directive <integrated_package_management>`.
Sometimes, it is necessary to access further source files that are in a directory relative to the current Snakefile.
Since workflows can be imported from remote locations (e.g. when using :ref:`modules <snakefiles-modules>`), it is important to not do this manually, so that Snakemake has the chance to cache these files locally before they are accessed.
This can be achieved by accessing their path via the ``workflow.get_source``, which (a) computes the correct path relative to the current Snakefile such that the file can be accessed from any working directory, and (b) downloads remote files to a local cache:

.. code-block:: python

    rule a:
        output:
            "test.out"
        params:
            json=workflow.source_path("../resources/test.json")
        shell:
            "somecommand {params.json} > {output}"
.. _snakefiles-reports:

-------
Reports
-------

From Snakemake 5.1 on, it is possible to automatically generate detailed self-contained HTML reports that encompass runtime statistics, provenance information, workflow topology and results.
**A realistic example report from a real workflow can be found** `here <https://koesterlab.github.io/resources/report.html>`_.

For including results into the report, the Snakefile has to be annotated with additional information.
Each output file that shall be part of the report has to be marked with the ``report`` flag, which optionally points to a caption in `restructured text format <https://docutils.sourceforge.io/docs/user/rst/quickstart.html>`_ and allows to define a ``category`` for grouping purposes.
Moreover, a global workflow description can be defined via the ``report`` directive.
Consider the following example:

.. code-block:: python

  report: "report/workflow.rst"


  rule all:
      input:
          ["fig1.svg", "fig2.png", "testdir"]


  rule c:
      output:
          "test.{i}.out"
      singularity:
          "docker://continuumio/miniconda3:4.4.10"
      conda:
          "envs/test.yaml"
      shell:
          "sleep `shuf -i 1-3 -n 1`; touch {output}"


  rule a:
      input:
          expand("test.{i}.out", i=range(10))
      output:
          report("fig1.svg", caption="report/fig1.rst", category="Step 1")
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig1.svg {output}"


  rule b:
      input:
          expand("{model}.{i}.out", i=range(10))
      output:
          report("fig2.png", caption="report/fig2.rst", category="Step 2", subcategory="{model}")
      shell:
          "sleep `shuf -i 1-3 -n 1`; cp data/fig2.png {output}"

  rule d:
      output:
          report(directory("testdir"), patterns=["{name}.txt"], caption="report/somedata.rst", category="Step 3")
      shell:
          "mkdir {output}; for i in 1 2 3; do echo $i > {output}/$i.txt; done"

As can be seen, we define a global description which is contained in the file ``report/workflow.rst``.
In addition, we mark ``fig1.svg`` and ``fig2.png`` for inclusion into the report, while in both cases specifying a caption text via again referring to a restructured text file.
Note the paths to the ``.rst``-files are interpreted relative to the current Snakefile.

Inside the ``.rst``-files you can use `Jinja2 <https://jinja.palletsprojects.com>`_ templating to access context information.
In case of the global description, you can access the config dictionary via ``{{ snakemake.config }}``, (e.g., use ``{{ snakemake.config["mykey"] }}`` to access the key ``mykey``).
In case of output files, you can access the same values as available with the :ref:`script directive <snakefiles-external_scripts>` (e.g., ``snakemake.wildcards``).

When marking files for inclusion in the report, a ``category`` and a ``subcategory`` can be given, allowing to group results in of the report.
For both, wildcards (like ``{model}`` see rule b in the example), are automatically replaced with the respective values from the corresponding job.

The last rule ``d`` creates a directory with several files, here mimicing the case that it is impossible to specify exactly which files will be created while writing the workflow (e.g. it might depend on the data).
Nevertheless, it is still possible to include those files one by one into the report by defining inclusion patterns (here ``patterns=["{name}.txt"]``) along with the report flag.
When creating the report, Snakemake will scan the directory for files matching the given patterns and include all of them in the report.
Wildcards in those patterns are made available in the jinja-templated caption document along with the rules wildcards in the ``snakemake.wildcards`` object.

If the output of a rule is a directory with an HTML file hierarchy, it is also possible to specify an entry-point HTML file for inclusion into the report, instead of the ``patterns`` approach from above.
This works as follows:

.. code-block:: python

    rule generate_html_hierarchy:
        output:
            report(directory("test"), caption="report/caption.rst", htmlindex="test.html")
        shell:
            """
            # mimic writing of an HTML hierarchy
            mkdir test
            cp template.html test/test.html
            mkdir test/js
            echo \"alert('test')\" > test/js/test.js
            """


Moreover, in every ``.rst`` document, you can link to

* the **Workflow** panel (with ``Rules_``),
* the **Statistics** panel (with ``Statistics_``),
* any **category** panel (with ``Mycategory_``, while ``Mycategory`` is the name given for the category argument of the report flag). E.g., with above example, you could write ``see `Step 2`_`` in order to link to the section with the results that have been assigned to the category ``Step 2``.
* any **file** marked with the report flag (with ``myfile.txt_``, while ``myfile.txt`` is the basename of the file, without any leading directories). E.g., with above example, you could write ``see fig2.png_`` in order to link to the result in the report document.

For details about the hyperlink mechanism of restructured text see `here <https://docutils.sourceforge.io/docs/user/rst/quickref.html#hyperlink-targets>`_.

To create the report simply run

.. code-block:: bash

    snakemake --report report.html

after your workflow has finished.
All other information contained in the report (e.g. runtime statistics) is automatically collected during creation.
These statistics are obtained from the metadata that is stored in the ``.snakemake`` directory inside your working directory.


You can define an institute specific stylesheet with:

.. code-block:: bash

    snakemake --report report.html --report-stylesheet custom-stylesheet.css

In particular, this allows you to e.g. set a logo at the top (by using CSS to inject a background for the placeholder ``<div id="brand">``, or overwrite colors.
For an example custom stylesheet defining the logo, see :download:`here <../../tests/test_report/custom-stylesheet.css>`.
The report for above example can be found :download:`here <../../tests/test_report/expected-results/report.html>` (with a custom branding for the University of Duisburg-Essen).
The full example source code can be found `here <https://github.com/snakemake/snakemake/tree/main/tests/test_report/>`_.

Note that the report can be restricted to particular jobs and results by specifying targets at the command line, analog to normal Snakemake execution.
For example, with

.. code-block:: bash

    snakemake fig1.svg --report report-short.html

the report contains only ``fig1.svg``.
.. _snakefiles-testing:

===================================
Automatically generating unit tests
===================================

Snakemake can automatically generate unit tests from a workflow that has already been successfully executed.
By running

.. code-block:: bash

    snakemake --generate-unit-tests

Snakemake is instructed to take one representative job for each rule and copy its input files to a hidden folder ``.tests/unit``,
along with generating test cases for Pytest_.

Importantly, note that such unit tests shall not be generated from big data, as they should usually be finished in a few seconds.
Further, it makes sense to store the generated unit tests in version control (e.g. git), such that huge files are not recommended.
Instead, we suggest to first execute the workflow that shall be tested with some kind of small dummy datasets, and then use the results thereof to generate the unit tests.
The small dummy datasets can in addition be used to generate an integration test, that could e.g. be stored under ``.tests/integration``, next to the unit tests.

Each auto-generated unit test is stored in a file ``.tests/unit/test_<rulename>.py``, and executes just the one representative job of the respective rule.
After successfull execution of the job, it will compare the obtained results with those that have been present when running ``snakemake --generate-unit-tests``.
By default, the comparison happens byte by byte (using ``cmp``). This behavior can be overwritten by modifying the test file.

.. _Pytest: https://pytest.org.. _api_reference_snakemake:

The Snakemake API
=================

.. autofunction:: snakemake.snakemake
.. _utils-api:

Additional utils
================

.. automodule:: snakemake.utils
    :members:
snakemake.report package
========================

Module contents
---------------

.. automodule:: snakemake.report
    :members:
    :undoc-members:
    :show-inheritance:
snakemake package
=================

Subpackages
-----------

.. toctree::

    snakemake.remote
    snakemake.report

Submodules
----------

snakemake.benchmark module
--------------------------

.. automodule:: snakemake.benchmark
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.checkpoints module
----------------------------

.. automodule:: snakemake.checkpoints
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.common module
-----------------------

.. automodule:: snakemake.common
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.conda module
----------------------

.. automodule:: snakemake.conda
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.cwl module
--------------------

.. automodule:: snakemake.cwl
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.dag module
--------------------

.. automodule:: snakemake.dag
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.decorators module
---------------------------

.. automodule:: snakemake.decorators
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.exceptions module
---------------------------

.. automodule:: snakemake.exceptions
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.executors module
--------------------------

.. automodule:: snakemake.executors
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.gui module
--------------------

.. automodule:: snakemake.gui
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.io module
-------------------

.. automodule:: snakemake.io
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.jobs module
---------------------

.. automodule:: snakemake.jobs
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.logging module
------------------------

.. automodule:: snakemake.logging
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.output\_index module
------------------------------

.. automodule:: snakemake.output_index
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.parser module
-----------------------

.. automodule:: snakemake.parser
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.persistence module
----------------------------

.. automodule:: snakemake.persistence
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.rules module
----------------------

.. automodule:: snakemake.rules
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.scheduler module
--------------------------

.. automodule:: snakemake.scheduler
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.script module
-----------------------

.. automodule:: snakemake.script
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.shell module
----------------------

.. automodule:: snakemake.shell
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.singularity module
----------------------------

.. automodule:: snakemake.singularity
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.stats module
----------------------

.. automodule:: snakemake.stats
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.utils module
----------------------

.. automodule:: snakemake.utils
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.workflow module
-------------------------

.. automodule:: snakemake.workflow
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.wrapper module
------------------------

.. automodule:: snakemake.wrapper
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: snakemake
    :members:
    :undoc-members:
    :show-inheritance:
Internal API
============

These pages document the entire internal API of Snakemake.

.. toctree::
   :maxdepth: 4

   snakemake
snakemake.remote package
========================

Submodules
----------

snakemake.remote.EGA module
---------------------------

.. automodule:: snakemake.remote.EGA
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.FTP module
---------------------------

.. automodule:: snakemake.remote.FTP
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.GS module
--------------------------

.. automodule:: snakemake.remote.GS
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.HTTP module
----------------------------

.. automodule:: snakemake.remote.HTTP
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.NCBI module
----------------------------

.. automodule:: snakemake.remote.NCBI
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.S3 module
--------------------------

.. automodule:: snakemake.remote.S3
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.S3Mocked module
--------------------------------

.. automodule:: snakemake.remote.S3Mocked
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.SFTP module
----------------------------

.. automodule:: snakemake.remote.SFTP
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.XRootD module
------------------------------

.. automodule:: snakemake.remote.XRootD
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.dropbox module
-------------------------------

.. automodule:: snakemake.remote.dropbox
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.gfal module
----------------------------

.. automodule:: snakemake.remote.gfal
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.gridftp module
-------------------------------

.. automodule:: snakemake.remote.gridftp
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.iRODS module
-----------------------------

.. automodule:: snakemake.remote.iRODS
    :members:
    :undoc-members:
    :show-inheritance:

snakemake.remote.webdav module
------------------------------

.. automodule:: snakemake.remote.webdav
    :members:
    :undoc-members:
    :show-inheritance:


Module contents
---------------

.. automodule:: snakemake.remote
    :members:
    :undoc-members:
    :show-inheritance:
.. _Miniconda: https://conda.pydata.org/miniconda.html
.. _Mambaforge: https://github.com/conda-forge/miniforge#mambaforge
.. _Mamba: https://github.com/mamba-org/mamba
.. _Conda: https://conda.pydata.org


.. _getting_started-installation:

============
Installation
============

Snakemake is available on PyPi as well as through Bioconda and also from source code.
You can use one of the following ways for installing Snakemake.

.. _conda-install:

Installation via Conda/Mamba
============================

This is the **recommended** way to install Snakemake,
because it also enables Snakemake to :ref:`handle software dependencies of your
workflow <integrated_package_management>`.

First, you have install a Conda-based Python3 distribution.
The recommended choice is Mambaforge_ which not only provides the required Python and Conda commands, 
but also includes Mamba_ an extremely fast and robust replacement for the Conda_ package manager which is highly recommended.
The default conda solver is a bit slow and sometimes has issues with `selecting the latest package releases <https://github.com/conda/conda/issues/9905>`_. 
Therefore, we recommend to in any case use Mamba_.

In case you don't use Mambaforge_ you can always install Mamba_ into any other Conda-based Python distribution with

.. code-block:: console

    $ conda install -n base -c conda-forge mamba

Full installation
-----------------

Snakemake can be installed with all goodies needed to run in any environment and for creating interactive reports via

.. code-block:: console

    $ conda activate base
    $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help

Installing into isolated environments is best practice in order to avoid side effects with other packages.

Note that full installation is not possible from **Windows**, because some of the dependencies are Unix (Linux/MacOS) only.
For Windows, please use the minimal installation below.

Minimal installation
--------------------

A minimal version of Snakemake which only depends on the bare necessities can be installed with

.. code-block:: console

    $ conda activate base
    $ mamba create -c bioconda -c conda-forge -n snakemake snakemake-minimal

In contrast to the full installation, which depends on some Unix (Linux/MacOS) only packages, this also works on Windows.

Notes on Bioconda as a package source
-------------------------------------

Note that Snakemake is available via Bioconda for historical, reproducibility, and continuity reasons (although it is not limited to biology applications at all).
However, it is easy to combine Snakemake installation with other channels, e.g., by prefixing the package name with ``::bioconda``, i.e.,

.. code-block:: console

    $ conda activate base
    $ mamba create -n some-env -c conda-forge bioconda::snakemake bioconda::snakemake-minimal ...

Installation via pip
====================

Instead of conda, snakemake can be installed with pip.
However, note that snakemake has non-python dependencies, such that the pip based installation has a limited functionality if those dependencies are not manually installed in addition.

A list of Snakemake's dependencies can be found within its `meta.yaml conda recipe <https://bioconda.github.io/recipes/snakemake/README.html>`_.
.. _project_info-contributing:

============
Contributing
============

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

You can contribute in many ways:


----------------------
Types of Contributions
----------------------


Report Bugs
===========

Report bugs at https://github.com/snakemake/snakemake/issues

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.


Fix Bugs
========

Look through the Github issues for bugs.
If you want to start working on a bug then please write short message on the issue tracker to prevent duplicate work.


Implement Features
==================

Look through the Github issues for features.
If you want to start working on an issue then please write short message on the issue tracker to prevent duplicate work.

Contributing a new cluster or cloud execution backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Execution backends are added by implementing a so-called ``Executor``.
All executors are located in `snakemake/executors/ <https://github.com/snakemake/snakemake/tree/main/snakemake/executors>`_.
In order to implement a new executor, you have to inherit from the class ``ClusterExecutor``.
Below you find a skeleton

.. code-block:: python

    class SkeletonExecutor(ClusterExecutor):
        def __init__(self, workflow, dag, cores,
                 jobname="snakejob.{name}.{jobid}.sh",
                 printreason=False,
                 quiet=False,
                 printshellcmds=False,
                 latency_wait=3,
                 cluster_config=None,
                 local_input=None,
                 restart_times=None,
                 exec_job=None,
                 assume_shared_fs=True,
                 max_status_checks_per_second=1):

            # overwrite the command to execute a single snakemake job if necessary
            # exec_job = "..."

            super().__init__(workflow, dag, None,
                             jobname=jobname,
                             printreason=printreason,
                             quiet=quiet,
                             printshellcmds=printshellcmds,
                             latency_wait=latency_wait,
                             cluster_config=cluster_config,
                             local_input=local_input,
                             restart_times=restart_times,
                             exec_job=exec_job,
                             assume_shared_fs=False,
                             max_status_checks_per_second=10)

            # add additional attributes

        def shutdown(self):
            # perform additional steps on shutdown if necessary
            super().shutdown()

        def cancel(self):
            for job in self.active_jobs:
                # cancel active jobs here
            self.shutdown()
        
        def run_jobs(self, jobs, callback=None, submit_callback=None, error_callback=None):
            """Run a list of jobs that is ready at a given point in time.

            By default, this method just runs each job individually.
            This behavior is inherited and therefore this method can be removed from the skeleton if the
            default behavior is intended.
            This method can be overwritten to submit many jobs in a more efficient way than one-by-one.

            Note that in any case, for each job, the callback functions have to be called individually!
            """
            for job in jobs:
                self.run(
                    job,
                    callback=callback,
                    submit_callback=submit_callback,
                    error_callback=error_callback,
                )

        def run(self, job,
                callback=None,
                submit_callback=None,
                error_callback=None):
            """Run an individual job or a job group.
            """

            super()._run(job)
            # obtain job execution command
            exec_job = self.format_job(
                self.exec_job, job, _quote_all=True,
                use_threads="--force-use-threads" if not job.is_group() else "")

            # submit job here, and obtain job ids from the backend

            # register job as active, using your own namedtuple.
            # The namedtuple must at least contain the attributes
            # job, jobid, callback, error_callback.
            self.active_jobs.append(MyJob(
                job, jobid, callback, error_callback))

        def _wait_for_jobs(self):
            # busy wait on job completion
            # This is only needed if your backend does not allow to use callbacks
            # for obtaining job status.
            while True:
                # always use self.lock to avoid race conditions
                with self.lock:
                    if not self.wait:
                        return
                    active_jobs = self.active_jobs
                    self.active_jobs = list()
                    still_running = list()
                for j in active_jobs:
                    # use self.status_rate_limiter to avoid too many API calls.
                    with self.status_rate_limiter:

                        # Retrieve status of job j from your backend via j.jobid
                        # Handle completion and errors, calling either j.callback(j.job)
                        # or j.error_callback(j.job)
                        # In case of error, add job j to still_running.
                with self.lock:
                    self.active_jobs.extend(still_running)
                sleep()


Write Documentation
===================

Snakemake could always use more documentation, whether as part of the official vcfpy docs, in docstrings, or even on the web in blog posts, articles, and such.

Snakemake uses `Sphinx <https://sphinx-doc.org>`_ for the user manual (that you are currently reading).
See `project_info-doc_guidelines` on how the documentation reStructuredText is used.


Submit Feedback
===============

The best way to send feedback is to file an issue at https://github.com/snakemake/snakemake/issues

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions are welcome :)

-----------------------
Pull Request Guidelines
-----------------------

To update the documentation, fix bugs or add new features you need to create a Pull Request
. A PR is a change you make to your local copy of the code for us to review and potentially integrate into the code base.

To create a Pull Request you need to do these steps:

1. Create a Github account.
2. Fork the repository.
3. Clone your fork locally.
4. Go to the created snakemake folder with :code:`cd snakemake`.
5. Create a new branch with :code:`git checkout -b <descriptive_branch_name>`.
6. Make your changes to the code or documentation.
7. Run :code:`git add .` to add all the changed files to the commit (to see what files will be added you can run :code:`git add . --dry-run`).
8. To commit the added files use :code:`git commit`. (This will open a command line editor to write a commit message. These should have a descriptive 80 line header, followed by an empty line, and then a description of what you did and why. To use your command line text editor of choice use (for example) :code:`export GIT_EDITOR=vim` before running :code:`git commit`).
9. Now you can push your changes to your Github copy of Snakemake by running :code:`git push origin <descriptive_branch_name>`.
10. If you now go to the webpage for your Github copy of Snakemake you should see a link in the sidebar called "Create Pull Request".
11. Now you need to choose your PR from the menu and click the "Create pull request" button. Be sure to change the pull request target branch to <descriptive_branch_name>!

If you want to create more pull requests, first run :code:`git checkout main` and then start at step 5. with a new branch name.

Feel free to ask questions about this if you want to contribute to Snakemake :)

------------------
Testing Guidelines
------------------

To ensure that you do not introduce bugs into Snakemake, you should test your code thouroughly.

To have integration tests run automatically when commiting code changes to Github, you need to sign up on wercker.com and register a user.

The easiest way to run your development version of Snakemake is perhaps to go to the folder containing your local copy of Snakemake and call:

.. code-block:: console

    $ conda env create -f test-environment.yml -n snakemake-testing
    $ conda activate snakemake-testing
    $ pip install -e .

This will make your development version of Snakemake the one called when running snakemake. You do not need to run this command after each time you make code changes.

From the base snakemake folder you call :code:`nosetests` to run all the tests, or choose one specific test. For this to work, Nose (the testing framework we use) can be installed to the conda environment using pip:

.. code-block:: console

   $ pip install nose
   $ nosetests
   $ nosetests tests.tests:test_log_input

If you introduce a new feature you should add a new test to the tests directory. See the folder for examples.

.. project_info-doc_guidelines:

------------------------
Documentation Guidelines
------------------------

For the documentation, please adhere to the following guidelines:

- Put each sentence on its own line, this makes tracking changes through Git SCM easier.
- Provide hyperlink targets, at least for the first two section levels.
  For this, use the format ``<document_part>-<section_name>``, e.g., ``project_info-doc_guidelines``.
- Use the section structure from below.

::

    .. document_part-heading_1:

    =========
    Heading 1
    =========


    .. document_part-heading_2:

    ---------
    Heading 2
    ---------


    .. document_part-heading_3:

    Heading 3
    =========


    .. document_part-heading_4:

    Heading 4
    ---------


    .. document_part-heading_5:

    Heading 5
    ~~~~~~~~~


    .. document_part-heading_6:

    Heading 6
    :::::::::

.. _doc_setup:

-------------------
Documentation Setup
-------------------

For building the documentation, you have to install the Sphinx.
If you have already installed Conda, all you need to do is to create a
Snakemake development environment via

.. code-block:: console

    $ git clone git@github.com:snakemake/snakemake.git
    $ cd snakemake
    $ conda env create -f doc-environment.yml -n snakemake

Then, the docs can be built with

.. code-block:: console

    $ conda activate snakemake
    $ cd docs
    $ make html
    $ make clean && make html  # force rebuild

Alternatively, you can use virtualenv.
The following assumes you have a working Python 3 setup.

.. code-block:: console

    $ git clone git@github.org:snakemake/snakemake.git
    $ cd snakemake/docs
    $ virtualenv -p python3 .venv
    $ source .venv/bin/activate
    $ pip install --upgrade -r requirements.txt

Afterwards, the docs can be built with

.. code-block:: console

    $ source .venv/bin/activate
    $ make html  # rebuild for changed files only
    $ make clean && make html  # force rebuild
.. _project_info-citations:

====================
Citing and Citations
====================

This section gives instructions on how to cite Snakemake and lists citing articles.


.. project_info-citing_snakemake:

----------------
Citing Snakemake
----------------

When using Snakemake for a publication, **please cite the following article** in you paper:

`Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33. <https://doi.org/10.12688/f1000research.29032.1>`_

This "rolling" paper will be regularly updated when Snakemake receives new features.

More References
===============

The initial Snakemake publication was:

`Köster, Johannes and Rahmann, Sven. "Snakemake - A scalable bioinformatics workflow engine". Bioinformatics 2012. <https://bioinformatics.oxfordjournals.org/content/28/19/2520>`_

Another publication describing more of Snakemake internals:

`Köster, Johannes and Rahmann, Sven. "Building and Documenting Bioinformatics Workflows with Python-based Snakemake". Proceedings of the GCB 2012. <https://drops.dagstuhl.de/opus/volltexte/oasics-complete/oasics-vol26-gcb2012-complete.pdf>`_

And my PhD thesis which describes all algorithmic details as of 2015:

`Johannes Köster, "Parallelization, Scalability, and Reproducibility in Next-Generation Sequencing Analysis", TU Dortmund 2014 <https://hdl.handle.net/2003/33940>`_

The most comprehensive publication is our "rolling" paper (see above):

`Mölder, F., Jablonski, K.P., Letcher, B., Hall, M.B., Tomkins-Tinch, C.H., Sochat, V., Forster, J., Lee, S., Twardziok, S.O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., Nahnsen, S., Köster, J., 2021. Sustainable data analysis with Snakemake. F1000Res 10, 33. <https://doi.org/10.12688/f1000research.29032.1>`_


Project Pages
=============

If you publish a Snakemake workflow, consider to add this badge to your project page:

.. image:: https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat
   :target: https://snakemake.readthedocs.io

The markdown syntax is

.. sourcecode:: text

    [![Snakemake](https://img.shields.io/badge/snakemake-≥5.6.0-brightgreen.svg?style=flat)](https://snakemake.readthedocs.io)

Replace the ``5.6.0`` with the minimum required Snakemake version.
You can also `change the style <https://shields.io/#styles>`_.
.. _project_info-more_resources:

==============
More Resources
==============

.. _project_info-talks_and_posters:

-----------------
Talks and Posters
-----------------

* `Poster at ECCB 2016, The Hague, Netherlands. <https://johanneskoester.bitbucket.io/posters/snakemake+bioconda-2016.pdf>`_
* `Invited talk by Johannes Köster at the Broad Institute, Boston 2015. <https://slides.com/johanneskoester/snakemake-broad-2015>`_
* `Introduction to Snakemake. Tutorial Slides presented by Johannes Köster at the GCB 2015, Dortmund, Germany. <https://slides.com/johanneskoester/deck-1>`_
* `Invited talk by Johannes Köster at the DTL Focus Meeting: "NGS Production Pipelines", Dutch Techcentre for Life Sciences, Utrecht 2014. <https://speakerdeck.com/johanneskoester/workflow-management-with-snakemake>`_
* `Taming Snakemake by Jeremy Leipzig, Bioinformatics software developer at Children's Hospital of Philadelphia, 2014. <https://de.slideshare.net/jermdemo/taming-snakemake>`_
* `"Snakemake makes ... snakes?" - An Introduction by Marcel Martin from SciLifeLab, Stockholm 2015 <https://marcelm.net/talks/2015/snakemake>`_
* `"Workflow Management with Snakemake" by Johannes Köster, 2015. Held at the Department of Biostatistics and Computational Biology, Dana-Farber Cancer Institute <https://speakerdeck.com/johanneskoester/workflow-management-with-snakemake-1>`_


.. _project_info-external_resources:

------------------
External Resources
------------------

These resources are not part of the official documentation.

* `A number of tutorials on the subject "Tools for reproducible research" <https://nbis-reproducible-research.readthedocs.io>`_
* `Snakemake workflow used for the Kallisto paper <https://github.com/pachterlab/kallisto_paper_analysis>`_
* `An alternative tutorial for Snakemake <https://slowkow.com/notes/snakemake-tutorial/>`_
* `An Emacs mode for Snakemake <https://melpa.org/#/snakemake-mode>`_
* `Flexible bioinformatics pipelines with Snakemake <http://watson.nci.nih.gov/~sdavis/blog/flexible_bioinformatics_pipelines_with_snakemake/>`_
* `Sandwiches with Snakemake <https://github.com/leipzig/SandwichesWithSnakemake>`_
* `A visualization of the past years of Snakemake development <https://youtu.be/bq3vXrWw1yk>`_
* `Japanese version of the Snakemake tutorial <https://github.com/joemphilips/Translate_Snakemake_Tutorial>`_
* `Basic <https://bioinfo-fr.net/snakemake-pour-les-nuls>`_ and `advanced <https://bioinfo-fr.net/snakemake-aller-plus-loin-avec-la-parallelisation>`_ french Snakemake tutorial.
* `Mini tutorial on Snakemake and Bioconda <https://github.com/dlaehnemann/TutMinicondaSnakemake>`_
* `Snakeparse: a utility to expose Snakemake workflow configuation via a command line interface <https://github.com/nh13/snakeparse>`_
.. project_info-authors:

=======
Credits
=======


Development Lead
----------------

- Johannes Köster

Development Team
----------------

- Christopher Tomkins-Tinch
- David Koppstein
- Tim Booth
- Manuel Holtgrewe
- Christian Arnold
- Wibowo Arindrarto
- Rasmus Ågren
- Soohyun Lee
- Vanessa Sochat

Contributors
------------

In alphabetical order

- Andreas Wilm
- Anthony Underwood
- Ryan Dale
- David Alexander
- Elias Kuthe
- Elmar Pruesse
- Hyeshik Chang
- Jay Hesselberth
- Jesper Foldager
- John Huddleston
- Joona Lehtomäki
- Karel Brinda
- Karl Gutwin
- Kemal Eren
- Kostis Anagnostopoulos
- Kyle A. Beauchamp
- Kyle Meyer
- Lance Parsons
- Manuel Holtgrewe
- Marcel Martin
- Matthew Shirley
- Mattias Franberg
- Matt Shirley
- Paul Moore
- Per Unneberg
- Ryan C. Thompson
- Ryan Dale
- Sean Davis
- Simon Ye
- Tobias Marschall
- Vanessa Sochat
- Willem Ligtenberg
.. _project_info-faq:

==========================
Frequently Asked Questions
==========================

.. contents::

What is the key idea of Snakemake workflows?
--------------------------------------------

The key idea is very similar to GNU Make. The workflow is determined automatically from top (the files you want) to bottom (the files you have), by applying very general rules with wildcards you give to Snakemake:

.. image:: img/idea.png
    :alt: Snakemake idea

When you start using Snakemake, please make sure to walk through the :ref:`official tutorial <tutorial>`.
It is crucial to understand how to properly use the system.

Snakemake does not connect my rules as I have expected, what can I do to debug my dependency structure?
-------------------------------------------------------------------------------------------------------

Since dependencies are inferred implicitly, results can sometimes be suprising when little errors are made in filenames or when input functions raise unexpected errors.
For debugging such cases, Snakemake provides the command line flag ``--debug-dag`` that leads to printing details each decision that is taken while determining the dependencies.

In addition, it is advisable to check whether certain intermediate files would be created by targetting them individually via the command line.

Finally, it is possible to constrain the rules that are considered for DAG creating via ``--allowed-rules``. 
This way, you can easily check rule by rule if it does what you expect.
However, note that ``--allowed-rules`` is only meant for debugging.
A workflow should always work fine without it.

My shell command fails with with errors about an "unbound variable", what's wrong?
----------------------------------------------------------------------------------

This happens often when calling virtual environments from within Snakemake. Snakemake is using `bash strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_, to ensure e.g. proper error behavior of shell scripts.
Unfortunately, virtualenv and some other tools violate bash strict mode.
The quick fix for virtualenv is to temporarily deactivate the check for unbound variables

.. code-block:: bash

    set +u; source /path/to/venv/bin/activate; set -u

For more details on bash strict mode, see the `here <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_.


My shell command fails with exit code != 0 from within a pipe, what's wrong?
----------------------------------------------------------------------------

Snakemake is using `bash strict mode <http://redsymbol.net/articles/unofficial-bash-strict-mode/>`_ to ensure best practice error reporting in shell commands.
This entails the pipefail option, which reports errors from within a pipe to outside. If you don't want this, e.g., to handle empty output in the pipe, you can disable pipefail via prepending

.. code-block:: bash

    set +o pipefail;

to your shell command in the problematic rule.


I don't want Snakemake to detect an error if my shell command exits with an exitcode > 1. What can I do?
---------------------------------------------------------------------------------------------------------

Sometimes, tools encode information in exit codes bigger than 1. Snakemake by default treats anything > 0 as an error. Special cases have to be added by yourself. For example, you can write

.. code-block:: python

    shell:
        """
        set +e
        somecommand ...
        exitcode=$?
        if [ $exitcode -eq 1 ]
        then
            exit 1
        else
            exit 0
        fi
        """

This way, Snakemake only treats exit code 1 as an error, and thinks that everything else is fine.
Note that such tools are an excellent use case for contributing a `wrapper <https://snakemake-wrappers.readthedocs.io>`_.


.. _glob-wildcards:

How do I run my rule on all files of a certain directory?
---------------------------------------------------------

In Snakemake, similar to GNU Make, the workflow is determined from the top, i.e. from the target files. Imagine you have a directory with files ``1.fastq, 2.fastq, 3.fastq, ...``, and you want to produce files ``1.bam, 2.bam, 3.bam, ...`` you should specify these as target files, using the ids ``1,2,3,...``. You could end up with at least two rules like this (or any number of intermediate steps):


.. code-block:: python

    IDS = "1 2 3 ...".split() # the list of desired ids

    # a pseudo-rule that collects the target files
    rule all:
        input:  expand("otherdir/{id}.bam", id=IDS)

    # a general rule using wildcards that does the work
    rule:
        input:  "thedir/{id}.fastq"
        output: "otherdir/{id}.bam"
        shell:  "..."

Snakemake will then go down the line and determine which files it needs from your initial directory.

In order to infer the IDs from present files, Snakemake provides the ``glob_wildcards`` function, e.g.

.. code-block:: python

    IDS, = glob_wildcards("thedir/{id}.fastq")

The function matches the given pattern against the files present in the filesystem and thereby infers the values for all wildcards in the pattern. A named tuple that contains a list of values for each wildcard is returned. Here, this named tuple has only one item, that is the list of values for the wildcard ``{id}``.

I don't want expand to use the product of every wildcard, what can I do?
------------------------------------------------------------------------

By default the expand function uses ``itertools.product`` to create every combination of the supplied wildcards.
Expand takes an optional, second positional argument which can customize how wildcards are combined.
To create the list ``["a_1.txt", "b_2.txt", "c_3.txt"]``, invoke expand as:
``expand("{sample}_{id}.txt", zip, sample=["a", "b", "c"], id=["1", "2", "3"])``

I don't want expand to use every wildcard, what can I do?
---------------------------------------------------------

Sometimes partially expanding wildcards is useful to define inputs which still depend on some wildcards.
Expand takes an optional keyword argument, allow_missing=True, that will format only wildcards which are supplied, leaving others as is.
To create the list ``["{sample}_1.txt", "{sample}_2.txt"]``, invoke expand as:
``expand("{sample}_{id}.txt", id=["1", "2"], allow_missing=True)``
If the filename contains the wildcard ``allow_missing``, it will be formatted normally:
``expand("{allow_missing}.txt", allow_missing=True)`` returns ``["True.txt"]``.


Snakemake complains about a cyclic dependency or a PeriodicWildcardError. What can I do?
----------------------------------------------------------------------------------------

One limitation of Snakemake is that graphs of jobs have to be acyclic (similar to GNU Make). This means, that no path in the graph may be a cycle. Although you might have considered this when designing your workflow, Snakemake sometimes runs into situations where a cyclic dependency cannot be avoided without further information, although the solution seems obvious for the developer. Consider the following example:

.. code-block:: text

    rule all:
        input:
            "a"

    rule unzip:
        input:
            "{sample}.tar.gz"
        output:
            "{sample}"
        shell:
            "tar -xf {input}"

If this workflow is executed with

.. code-block:: console

    snakemake -n

two things may happen.

1. If the file ``a.tar.gz`` is present in the filesystem, Snakemake will propose the following (expected and correct) plan:

    .. code-block:: text

        rule a:
	        input: a.tar.gz
    	    output: a
    	    wildcards: sample=a
        localrule all:
	        input: a
        Job counts:
	        count	jobs
	        1	a
	        1	all
	        2

2. If the file ``a.tar.gz`` is not present and cannot be created by any other rule than rule ``a``, Snakemake will try to run rule ``a`` again, with ``{sample}=a.tar.gz``. This would infinitely go on recursively. Snakemake detects this case and produces a ``PeriodicWildcardError``.

In summary, ``PeriodicWildcardErrors`` hint to a problem where a rule or a set of rules can be applied to create its own input. If you are lucky, Snakemake can be smart and avoid the error by stopping the recursion if a file exists in the filesystem. Importantly, however, bugs upstream of that rule can manifest as ``PeriodicWildcardError``, although in reality just a file is missing or named differently.
In such cases, it is best to restrict the wildcard of the output file(s), or follow the general rule of putting output files of different rules into unique subfolders of your working directory. This way, you can discover the true source of your error.


Is it possible to pass variable values to the workflow via the command line?
----------------------------------------------------------------------------

Yes, this is possible. Have a look at :ref:`snakefiles_configuration`.
Previously it was necessary to use environment variables like so:
E.g. write

.. code-block:: bash

    $ SAMPLES="1 2 3 4 5" snakemake

and have in the Snakefile some Python code that reads this environment variable, i.e.

.. code-block:: python

    SAMPLES = os.environ.get("SAMPLES", "10 20").split()

I get a NameError with my shell command. Are braces unsupported?
----------------------------------------------------------------

You can use the entire Python `format minilanguage <https://docs.python.org/3/library/string.html#formatspec>`_ in shell commands. Braces in shell commands that are not intended to insert variable values thus have to be escaped by doubling them:

This:

.. code-block:: python

    ...
    shell: "awk '{print $1}' {input}"

becomes:

.. code-block:: python

    ...
    shell: "awk '{{print $1}}' {input}"

Here the double braces are escapes, i.e. there will remain single braces in the final command. In contrast, ``{input}`` is replaced with an input filename.

In addition, if your shell command has literal backslashes, ``\\``, you must escape them with a backslash, ``\\\\``. For example:

This:

.. code-block:: python

    shell: """printf \">%s\"" {{input}}"""

becomes:

.. code-block:: python

    shell: """printf \\">%s\\"" {{input}}"""

How do I incorporate files that do not follow a consistent naming scheme?
-------------------------------------------------------------------------

The best solution is to have a dictionary that translates a sample id to the inconsistently named files and use a function (see :ref:`snakefiles-input_functions`) to provide an input file like this:

.. code-block:: python

    FILENAME = dict(...)  # map sample ids to the irregular filenames here

    rule:
        # use a function as input to delegate to the correct filename
        input: lambda wildcards: FILENAME[wildcards.sample]
        output: "somefolder/{sample}.csv"
        shell: ...

How do I force Snakemake to rerun all jobs from the rule I just edited?
-----------------------------------------------------------------------

This can be done by invoking Snakemake with the ``--forcerules`` or ``-R`` flag, followed by the rules that should be re-executed:

.. code-block:: console

    $ snakemake -R somerule

This will cause Snakemake to re-run all jobs of that rule and everything downstream (i.e. directly or indirectly depending on the rules output).

How should Snakefiles be formatted?
--------------------------------------

To ensure readability and consistency, you can format Snakefiles with our tool `snakefmt <https://github.com/snakemake/snakefmt>`_. 

Python code gets formatted with `black <https://github.com/psf/black>`_ and Snakemake-specific blocks are formatted using similar principles (such as `PEP8 <https://www.python.org/dev/peps/pep-0008/>`_).

How do I enable syntax highlighting in Vim for Snakefiles?
----------------------------------------------------------

Instructions for doing this are located `here
<https://github.com/snakemake/snakemake/tree/main/misc/vim>`_.

Note that you can also format Snakefiles in Vim using :ref:`snakefmt
<How should Snakefiles be formatted?>`, with instructions located `here
<https://github.com/snakemake/snakefmt/blob/master/docs/editor_integration.md#vim>`_!

I want to import some helper functions from another python file. Is that possible?
----------------------------------------------------------------------------------

Yes, from version 2.4.8 on, Snakemake allows to import python modules (and also simple python files) from the same directory where the Snakefile resides.

How can I run Snakemake on a cluster where its main process is not allowed to run on the head node?
---------------------------------------------------------------------------------------------------

This can be achived by submitting the main Snakemake invocation as a job to the cluster. If it is not allowed to submit a job from a non-head cluster node, you can provide a submit command that goes back to the head node before submitting:

.. code-block:: bash

    qsub -N PIPE -cwd -j yes python snakemake --cluster "ssh user@headnode_address 'qsub -N pipe_task -j yes -cwd -S /bin/sh ' " -j

This hint was provided by Inti Pedroso.

Can the output of a rule be a symlink?
--------------------------------------

Yes. As of Snakemake 3.8, output files are removed before running a rule and then touched after the rule completes to ensure they are newer than the input.  Symlinks are treated just the same as normal files in this regard, and Snakemake ensures that it only modifies the link and not the target when doing this.

Here is an example where you want to merge N files together, but if N == 1 a symlink will do.  This is easier than attempting to implement workflow logic that skips the step entirely.  Note the **-r** flag, supported by modern versions of ln, is useful to achieve correct linking between files in subdirectories.

.. code-block:: python

    rule merge_files:
        output: "{foo}/all_merged.txt"
        input: my_input_func  # some function that yields 1 or more files to merge
        run:
            if len(input) > 1:
                shell("cat {input} | sort > {output}")
            else:
                shell("ln -sr {input} {output}")

Do be careful with symlinks in combination with :ref:`tutorial_temp-and-protected-files`.
When the original file is deleted, this can cause various errors once the symlink does not point to a valid file any more.

If you get a message like ``Unable to set utime on symlink .... Your Python build does not support it.`` this means that Snakemake is unable to properly adjust the modification time of the symlink.
In this case, a workaround is to add the shell command `touch -h {output}` to the end of the rule.

Can the input of a rule be a symlink?
-------------------------------------

Yes.  In this case, since Snakemake 3.8, one extra consideration is applied.  If *either* the link itself or the target of the link is newer than the output files for the rule then it will trigger the rule to be re-run.

I would like to receive a mail upon snakemake exit. How can this be achieved?
-----------------------------------------------------------------------------

On unix, you can make use of the commonly pre-installed `mail` command:

.. code-block:: bash

    snakemake 2> snakemake.log
    mail -s "snakemake finished" youremail@provider.com < snakemake.log

In case your administrator does not provide you with a proper configuration of the sendmail framework, you can configure `mail` to work e.g. via Gmail (see `here <https://www.cyberciti.biz/tips/linux-use-gmail-as-a-smarthost.html>`_).

I want to pass variables between rules. Is that possible?
---------------------------------------------------------

Because of the cluster support and the ability to resume a workflow where you stopped last time, Snakemake in general should be used in a way that information is stored in the output files of your jobs. Sometimes it might though be handy to have a kind of persistent storage for simple values between jobs and rules. Using plain python objects like a global dict for this will not work as each job is run in a separate process by snakemake. What helps here is the `PersistentDict` from the `pytools <https://github.com/inducer/pytools>`_ package. Here is an example of a Snakemake workflow using this facility:

.. code-block:: python

    from pytools.persistent_dict import PersistentDict

    storage = PersistentDict("mystorage")

    rule a:
        input: "test.in"
        output: "test.out"
        run:
            myvar = storage.fetch("myvar")
            # do stuff

    rule b:
        output: temp("test.in")
        run:
            storage.store("myvar", 3.14)

Here, the output rule b has to be temp in order to ensure that ``myvar`` is stored in each run of the workflow as rule a relies on it. In other words, the PersistentDict is persistent between the job processes, but not between different runs of this workflow. If you need to conserve information between different runs, use output files for them.

Why do my global variables behave strangely when I run my job on a cluster?
---------------------------------------------------------------------------

This is closely related to the question above.  Any Python code you put outside of a rule definition is normally run once before Snakemake starts to process rules, but on a cluster it is re-run again for each submitted job, because Snakemake implements jobs by re-running itself.

Consider the following...

.. code-block:: python

    from mydatabase import get_connection

    dbh = get_connection()
    latest_parameters = dbh.get_params().latest()

    rule a:
        input: "{foo}.in"
        output: "{foo}.out"
        shell: "do_op -params {latest_parameters}  {input} {output}"


When run a single machine, you will see a single connection to your database and get a single value for *latest_parameters* for the duration of the run.  On a cluster you will see a connection attempt from the cluster node for each job submitted, regardless of whether it happens to involve rule a or not, and the parameters will be recalculated for each job.

I want to configure the behavior of my shell for all rules. How can that be achieved with Snakemake?
----------------------------------------------------------------------------------------------------

You can set a prefix that will prepended to all shell commands by adding e.g.

.. code-block:: python

    shell.prefix("set -o pipefail; ")

to the top of your Snakefile. Make sure that the prefix ends with a semicolon, such that it will not interfere with the subsequent commands.
To simulate a bash login shell, you can do the following:

.. code-block:: python

    shell.executable("/bin/bash")
    shell.prefix("source ~/.bashrc; ")

Some command line arguments like --config cannot be followed by rule or file targets. Is that intended behavior?
----------------------------------------------------------------------------------------------------------------

This is a limitation of the argparse module, which cannot distinguish between the perhaps next arg of ``--config`` and a target.
As a solution, you can put the `--config` at the end of your invocation, or prepend the target with a single ``--``, i.e.


.. code-block:: console

    $ snakemake --config foo=bar -- mytarget
    $ snakemake mytarget --config foo=bar


How do I enforce config values given at the command line to be interpreted as strings?
--------------------------------------------------------------------------------------

When passing config values like this

.. code-block:: console

    $ snakemake --config version=2018_1

Snakemake will first try to interpret the given value as number.
Only if that fails, it will interpret the value as string.
Here, it does not fail, because the underscore `_` is interpreted as thousand separator.
In order to ensure that the value is interpreted as string, you have to pass it in quotes.
Since bash otherwise automatically removes quotes, you have to also wrap the entire entry into quotes, e.g.:

.. code-block:: console

    $ snakemake --config 'version="2018_1"'


How do I make my rule fail if an output file is empty?
------------------------------------------------------

Snakemake expects shell commands to behave properly, meaning that failures should cause an exit status other than zero. If a command does not exit with a status other than zero, Snakemake assumes everything worked fine, even if output files are empty. This is because empty output files are also a reasonable tool to indicate progress where no real output was produced. However, sometimes you will have to deal with tools that do not properly report their failure with an exit status. Here, the recommended way is to use bash to check for non-empty output files, e.g.:

.. code-block:: python

    rule:
        input:  ...
        output: "my/output/file.txt"
        shell:  "somecommand {input} {output} && [[ -s {output} ]]"


How does Snakemake lock the working directory?
----------------------------------------------

Per default, Snakemake will lock a working directory by output and input files. Two Snakemake instances that want to create the same output file are not possible. Two instances creating disjoint sets of output files are possible.
With the command line option ``--nolock``, you can disable this mechanism on your own risk. With ``--unlock``, you can be remove a stale lock. Stale locks can appear if your machine is powered off with a running Snakemake instance.


Snakemake does not trigger re-runs if I add additional input files. What can I do?
----------------------------------------------------------------------------------

Snakemake has a kind of "lazy" policy about added input files if their modification date is older than that of the output files. One reason is that information what to do cannot be inferred just from the input and output files. You need additional information about the last run to be stored. Since behaviour would be inconsistent between cases where that information is available and where it is not, this functionality has been encoded as an extra switch. To trigger updates for jobs with changed input files, you can use the command line argument --list-input-changes in the following way:

.. code-block:: console

    $ snakemake -n -R `snakemake --list-input-changes`

Here, ``snakemake --list-input-changes`` returns the list of output files with changed input files, which is fed into ``-R`` to trigger a re-run.


How do I trigger re-runs for rules with updated code or parameters?
-------------------------------------------------------------------

Similar to the solution above, you can use

.. code-block:: console

    $ snakemake -n -R `snakemake --list-params-changes`

and

.. code-block:: console


    $ snakemake -n -R `snakemake --list-code-changes`

Again, the list commands in backticks return the list of output files with changes, which are fed into ``-R`` to trigger a re-run.


How do I remove all files created by snakemake, i.e. like ``make clean``
------------------------------------------------------------------------

To remove all files created by snakemake as output files to start from scratch, you can use

.. code-block:: console

    $ snakemake some_target --delete-all-output

Only files that are output of snakemake rules will be removed, not those that serve as primary inputs to the workflow.
Note that this will only affect the files involved in reaching the specified target(s).
It is strongly advised to first run together with ``--dry-run`` to list the files that would be removed without actually deleting anything.
The flag ``--delete-temp-output`` can be used in a similar manner to only delete files flagged as temporary.


Why can't I use the conda directive with a run block?
-----------------------------------------------------

The run block of a rule (see :ref:`snakefiles-rules`) has access to anything defined in the Snakefile, outside of the rule.
Hence, it has to share the conda environment with the main Snakemake process.
To avoid confusion we therefore disallow the conda directive together with the run block.
It is recommended to use the script directive instead (see :ref:`snakefiles-external_scripts`).


My workflow is very large, how do I stop Snakemake from printing all this rule/job information in a dry-run?
------------------------------------------------------------------------------------------------------------

Indeed, the information for each individual job can slow down a dry-run if there are tens of thousands of jobs.
If you are just interested in the final summary, you can use the ``--quiet`` flag to suppress this.

.. code-block:: console

    $ snakemake -n --quiet

Git is messing up the modification times of my input files, what can I do?
--------------------------------------------------------------------------

When you checkout a git repository, the modification times of updated files are set to the time of the checkout. If you rely on these files as input **and** output files in your workflow, this can cause trouble. For example, Snakemake could think that a certain (git-tracked) output has to be re-executed, just because its input has been checked out a bit later. In such cases, it is advisable to set the file modification dates to the last commit date after an update has been pulled. One solution is to add the following lines to your ``.bashrc`` (or similar):

.. code-block:: bash

    gitmtim(){
        local f
        for f; do
            touch -d @0`git log --pretty=%at -n1 -- "$f"` "$f"
        done
    }
    gitmodtimes(){
        for f in $(git ls-tree -r $(git rev-parse --abbrev-ref HEAD) --name-only); do
            gitmtim $f
        done
    }

(inspired by the answer `here <https://stackoverflow.com/questions/2458042/restore-files-modification-time-in-git/22638823#22638823>`_). You can then run ``gitmodtimes`` to update the modification times of all tracked files on the current branch to their last commit time in git; BE CAREFUL--this does not account for local changes that have not been commited.

How do I exit a running Snakemake workflow?
-------------------------------------------

There are two ways to exit a currently running workflow.

1. If you want to kill all running jobs, hit Ctrl+C. Note that when using ``--cluster``, this will only cancel the main Snakemake process.
2. If you want to stop the scheduling of new jobs and wait for all running jobs to be finished, you can send a TERM signal, e.g., via

   .. code-block:: bash

       killall -TERM snakemake

How can I make use of node-local storage when running cluster jobs?
-------------------------------------------------------------------
When running jobs on a cluster you might want to make use of a node-local scratch
directory in order to reduce cluster network traffic and/or get more efficient disk
storage for temporary files. There is currently no way of doing this in Snakemake,
but a possible workaround involves the ``shadow`` directive and setting the
``--shadow-prefix`` flag to e.g. ``/scratch``.

.. code-block:: python

  rule:
      output:
          "some_summary_statistics.txt"
      shadow: "minimal"
      shell:
          """
          generate huge_file.csv
          summarize huge_file.csv > {output}
          """

The following would then lead to the job being executed in ``/scratch/shadow/some_unique_hash/``, and the
temporary file ``huge_file.csv`` could be kept at the compute node.

.. code-block:: console

   $ snakemake --shadow-prefix /scratch some_summary_statistics.txt --cluster ...

If you want the input files of your rule to be copied to the node-local scratch directory
instead of just using symbolic links, you can use ``copy-minimal`` in the ``shadow`` directive.
This is useful for example for benchmarking tools as a black-box.

.. code-block:: python

  rule:
      input:
          "input_file.txt"
      output:
          file = "output_file.txt",
          benchmark = "benchmark_results.txt",
      shadow: "copy-minimal"
      shell:
          """
          /usr/bin/time -v command "{input}" "{output.file}" > "{output.benchmark}"
          """

Executing snakemake as above then leads to the shell script accessing only node-local storage.

How do I access elements of input or output by a variable index?
----------------------------------------------------------------

Assuming you have something like the following rule

   .. code-block:: python

      rule a:
          output:
              expand("test.{i}.out", i=range(20))
          run:
              for i in range(20):
                  shell("echo test > {output[i]}")

Snakemake will fail upon execution with the error ``'OutputFiles' object has no attribute 'i'``. The reason is that the shell command is using the `Python format mini language <https://docs.python.org/3/library/string.html#formatspec>`_, which does only allow indexing via constants, e.g., ``output[1]``, but not via variables. Variables are treated as attribute names instead. The solution is to write

   .. code-block:: python

      rule a:
          output:
              expand("test.{i}.out", i=range(20))
          run:
              for i in range(20):
                  f = output[i]
                  shell("echo test > {f}")

or, more concise in this special case:

   .. code-block:: python

      rule a:
          output:
              expand("test.{i}.out", i=range(20))
          run:
              for f in output:
                  shell("echo test > {f}")

There is a compiler error when installing Snakemake with pip or easy_install, what shall I do?
----------------------------------------------------------------------------------------------

Snakemake itself is plain Python, hence the compiler error must come from one of the dependencies, like e.g., datrie.
You should have a look if maybe you are missing some library or a certain compiler package.
If everything seems fine, please report to the upstream developers of the failing dependency.

Note that in general it is recommended to install Snakemake via `Conda <https://conda.io>`_ which gives you precompiled packages and the additional benefit of having :ref:`automatic software deployment <integrated_package_management>` integrated into your workflow execution.

How to enable autocompletion for the zsh shell?
-----------------------------------------------

For users of the `Z shell <https://www.zsh.org/>`_ (zsh), just run the following (assuming an activated zsh) to activate autocompletion for snakemake:

.. code-block:: console

    compdef _gnu_generic snakemake

Example:
Say you have forgotten how to use the various options starting ``force``, just type the partial match i.e. ``--force`` which results in a list of all potential hits along with a description:


.. code-block:: console

    $snakemake --force**pressing tab**

    --force              -- Force the execution of the selected target or the
    --force-use-threads  -- Force threads rather than processes. Helpful if shared
    --forceall           -- Force the execution of the selected (or the first)
    --forcerun           -- (TARGET (TARGET ...)), -R (TARGET (TARGET ...))

To activate this autocompletion permanently, put this line in ``~/.zshrc``.

`Here <https://github.com/zsh-users/zsh-completions/blob/master/zsh-completions-howto.org>`_ is some further reading.
.. _tutorial-azure-aks:

Auto-scaling Azure Kubernetes cluster without shared filesystem
---------------------------------------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Python: https://www.python.org/

In this tutorial we will show how to execute a Snakemake workflow
on an auto-scaling Azure Kubernetes cluster without a shared file-system.
While Kubernetes is mainly known as microservice orchestration system with
self-healing properties, we will use it here simply as auto-scaling
compute orchestrator. One could use `persistent volumes in
Kubernetes <https://docs.microsoft.com/en-us/azure/aks/azure-files-dynamic-pv>`__
as shared file system, but this adds an unnecessary level of complexity
and most importantly costs. Instead we use cheap Azure Blob storage,
which is used by Snakemake to automatically stage data in and out for
every job.

Following the steps below you will

#. set up Azure Blob storage, download the Snakemake tutorial data and upload to Azure
#. then create an Azure Kubernetes (AKS) cluster
#. and finally run the analysis with Snakemake on the cluster 


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ ≥5.17

You should install conda as outlined in the :ref:`tutorial <tutorial-setup>`,
and then install full snakemake with:

.. code:: console

    conda create -c bioconda -c conda-forge -n snakemake snakemake

Make sure that the ``kubernetes`` and ``azure-storage-blob`` modules are installed
in this environment. Should they be missing install with:

.. code:: console

   pip install kubernetes
   pip install azure-storage-blob

In addition you will need the
`Azure CLI command <https://docs.microsoft.com/en-us/cli/azure/install-azure-cli?view=azure-cli-latest>`__ 
installed.

Create an Azure storage account and upload data
:::::::::::::::::::::::::::::::::::::::::::::::

We will be starting from scratch, i.e. we will 
create a new resource group and storage account. You can obviously reuse 
existing resources instead.

.. code:: console

   # change the following names as required
   # azure region where to run:
   region=southeastasia
   # name of the resource group to create:
   resgroup=snakemaks-rg
   # name of storage account to create (all lowercase, no hyphens etc.):
   stgacct=snakemaksstg

   # create a resource group with name and in region as defined above
   az group create --name $resgroup --location $region
   # create a general purpose storage account with cheapest SKU
   az storage account create -n $stgacct -g $resgroup --sku Standard_LRS -l $region

Get a key for that account and save it as ``stgkey`` for later use:

.. code:: console

   stgkey=$(az storage account keys list -g $resgroup -n $storageacct | head -n1 | cut -f 3)

Next, you will create a storage container (think: bucket) to upload the Snakemake tutorial data to:

.. code:: console

   az storage container create --resource-group $resgroup --account-name $stgacct \
       --account-key $stgkey --name snakemake-tutorial
   cd /tmp
   git clone https://github.com/snakemake/snakemake-tutorial-data.git
   cd snakemake-tutorial-data
   az storage blob upload-batch -d snakemake-tutorial --account-name $stgacct \
       --account-key $stgkey -s data/ --destination-path data

We are using `az storage blob` for uploading, because that `az` is already installed.
A  more efficient way of uploading would be to use
`azcopy <https://docs.microsoft.com/en-us/azure/storage/common/storage-use-azcopy-v10>`__.

Create an auto-scaling Kubernetes cluster
:::::::::::::::::::::::::::::::::::::::::

.. code:: console

   # change the cluster name as you like
   clustername=snakemaks-aks
   az aks create --resource-group $resgroup --name $clustername \
       --vm-set-type VirtualMachineScaleSets --load-balancer-sku standard --enable-cluster-autoscaler \
       --node-count 1 --min-count 1 --max-count 3 --node-vm-size Standard_D3_v2

There is a lot going on here, so let’s unpack it: this creates an
`auto-scaling Kubernetes
cluster <https://docs.microsoft.com/en-us/azure/aks/cluster-autoscaler>`__
(``--enable-cluster-autoscaler``) called ``$clustername`` (i.e. ``snakemaks-aks``), which starts
out with one node (``--node-count 1``) and has a maximum of three nodes
(``--min-count 1 --max-count 3``). For real world applications you will
want to increase the maximum count and also increase the VM size. You
could for example choose a large instance from the DSv2 series and add a
larger disk with (``--node-osdisk-size``) if needed. See `here for more
info on Linux VM
sizes <https://docs.microsoft.com/en-us/azure/virtual-machines/linux/sizes>`__.

Note, if you are creating the cluster in the Azure portal, click on the
ellipsis under node-pools to find the auto-scaling option.

Next, let’s fetch the credentials for this cluster, so that we can
actually interact with it.

.. code:: console

   az aks get-credentials --resource-group $resgroup --name $clustername
   # print basic cluster info
   kubectl cluster-info



Running the workflow
::::::::::::::::::::

Below we will task Snakemake to install software on the fly with conda.
For this we need a Snakefile with corresponding conda environment
yaml files. You can download the package containing all those files `here <workflow/snakedir.zip>`__.
After downloading, unzip it and cd into the newly created directory.

.. code:: console

   $ cd /tmp
   $ unzip ~/Downloads/snakedir.zip
   $ cd snakedir
   $ find .
   .
   ./Snakefile
   ./envs
   ./envs/calling.yaml
   ./envs/mapping.yaml


Now, we will need to setup the credentials that allow the Kubernetes nodes to
read and write from blob storage. For the AzBlob storage provider in
Snakemake this is done through the environment variables
``AZ_BLOB_ACCOUNT_URL`` and optionally ``AZ_BLOB_CREDENTIAL``. See the
`documentation <snakefiles/remote_files.html#microsoft-azure-storage>`__ for more info.
``AZ_BLOB_ACCOUNT_URL`` takes the form
``https://<accountname>.blob.core.windows.net`` and may also contain a
shared access signature (SAS), which is a powerful way to define fine grained
and even time controlled access to storage on Azure. The SAS can be part of the
URL, but if it’s missing, then you can set it with
``AZ_BLOB_CREDENTIAL`` or alternatively use the storage account key. To
keep things simple we’ll use the storage key here, since we already have it available,
but a SAS is generally more powerful. We’ll pass those variables on to the Kubernetes
with ``--envvars`` (see below).

Now you are ready to run the analysis:

.. code:: console

   export AZ_BLOB_ACCOUNT_URL="https://${stgacct}.blob.core.windows.net"
   export AZ_BLOB_CREDENTIAL="$stgkey"
   snakemake --kubernetes \
       --default-remote-prefix snakemake-tutorial --default-remote-provider AzBlob \
       --envvars AZ_BLOB_ACCOUNT_URL AZ_BLOB_CREDENTIAL --use-conda --jobs 3

This will use the default Snakemake image from Dockerhub. If you would like to use your
own, make sure that the image contains the same Snakemake version as installed locally
and also supports Azure Blob storage. If you plan to use your own image hosted on
 Azure Container Registries (ACR), make sure to attach the ACR to your Kubernetes 
 cluster. See `here <https://docs.microsoft.com/en-us/azure/aks/cluster-container-registry-integration>`__ for more info.


While Snakemake is running the workflow, it prints handy debug
statements per job, e.g.:

.. code:: console

   kubectl describe pod snakejob-c4d9bf9e-9076-576b-a1f9-736ec82afc64
   kubectl logs snakejob-c4d9bf9e-9076-576b-a1f9-736ec82afc64

With these you can also follow the scale-up of the cluster:

.. code:: console

   Events:
   Type     Reason             Age                From                Message
   ----     ------             ----               ----                -------
   Warning  FailedScheduling   60s (x3 over 62s)  default-scheduler   0/1 nodes are available: 1 Insufficient cpu.
   Normal   TriggeredScaleUp   50s                cluster-autoscaler  pod triggered scale-up: [{aks-nodepool1-17839284-vmss 1->3 (max: 3)}]

After a while you will see three nodes (each running one BWA job), which
was defined as the maximum above while creating your Kubernetes cluster:

.. code:: console

   $ kubectl get nodes
   NAME                                STATUS   ROLES   AGE   VERSION
   aks-nodepool1-17839284-vmss000000   Ready    agent   74m   v1.15.11
   aks-nodepool1-17839284-vmss000001   Ready    agent   11s   v1.15.11
   aks-nodepool1-17839284-vmss000002   Ready    agent   62s   v1.15.11

To get detailed information including historical data about used
resources, check Insights in the Azure portal under your AKS cluster
Monitoring/Insights. The alternative is an instant snapshot on the
command line:

::

   $ kubectl top node
   NAME                                CPU(cores)   CPU%   MEMORY(bytes)   MEMORY%
   aks-nodepool1-17839284-vmss000000   217m         5%     1796Mi          16%
   aks-nodepool1-17839284-vmss000001   1973m        51%    529Mi           4%
   aks-nodepool1-17839284-vmss000002   698m         18%    1485Mi          13%

After completion all results including
logs can be found in the blob container. You will also find results
listed in the first Snakefile target downloaded to the working directoy.

::

   $ find snakemake-tutorial/
   snakemake-tutorial/
   snakemake-tutorial/calls
   snakemake-tutorial/calls/all.vcf


   $ az storage blob list  --container-name snakemake-tutorial --account-name $stgacct --account-key $stgkey -o table
   Name                     Blob Type    Blob Tier    Length    Content Type                       Last Modified              Snapshot
   -----------------------  -----------  -----------  --------  ---------------------------------  -------------------------  ----------
   calls/all.vcf            BlockBlob    Hot          90986     application/octet-stream           2020-06-08T05:11:31+00:00
   data/genome.fa           BlockBlob    Hot          234112    application/octet-stream           2020-06-08T03:26:54+00:00
   # etc.
   logs/mapped_reads/A.log  BlockBlob    Hot          346       application/octet-stream           2020-06-08T04:59:50+00:00
   mapped_reads/A.bam       BlockBlob    Hot          2258058   application/octet-stream           2020-06-08T04:59:50+00:00
   sorted_reads/A.bam       BlockBlob    Hot          2244660   application/octet-stream           2020-06-08T05:03:41+00:00
   sorted_reads/A.bam.bai   BlockBlob    Hot          344       application/octet-stream           2020-06-08T05:06:25+00:00
   # same for samples B and C

Now that the execution is complete, the AKS cluster will scale down
automatically. If you are not planning to run anything else, it makes
sense to shut down it down entirely:

::

   az aks delete --name akscluster --resource-group $resgroup


.. _tutorial-google-lifesciences:

Google Life Sciences Tutorial
------------------------------

.. _Snakemake: http://snakemake.readthedocs.io
.. _Snakemake Remotes: https://snakemake.readthedocs.io/en/stable/snakefiles/remote_files.html
.. _Python: https://www.python.org/


Setup
:::::

To go through this tutorial, you need the following software installed:

* Python_ ≥3.5
* Snakemake_ ≥5.16
* git

First, you have to install the Miniconda Python3 distribution.
See `here <https://conda.io/en/latest/miniconda.html>`_ for installation instructions.
Make sure to ...

* Install the **Python 3** version of Miniconda.
* Answer yes to the question whether conda shall be put into your PATH.

The default conda solver is a bit slow and sometimes has issues with `selecting the latest package releases <https://github.com/conda/conda/issues/9905>`_. Therefore, we recommend to install `Mamba <https://github.com/QuantStack/mamba>`_ as a drop-in replacement via

.. code-block:: console

    $ conda install -c conda-forge mamba

Then, you can install Snakemake with

.. code-block:: console

    $ mamba create -c conda-forge -c bioconda -n snakemake snakemake

from the `Bioconda <https://bioconda.github.io>`_ channel.
This will install snakemake into an isolated software environment, that has to be activated with

.. code-block:: console

    $ conda activate snakemake
    $ snakemake --help



Credentials
:::::::::::

Using the Google Life Sciences executor with Snakemake requires the environment 
variable `GOOGLE_APPLICATION_CREDENTIALS` exported, which should point to
the full path of the file on your local machine. To generate this file, you
can refer to the page under iam-admin to `download your service account <https://console.cloud.google.com/iam-admin/iam>`_ key and export it to the environment.

.. code:: console

    export GOOGLE_APPLICATION_CREDENTIALS="/home/[username]/credentials.json"


The suggested, minimal permissions required for this role include the following:

 - Compute Storage Admin(Can potentially be restricted further)
 - Compute Viewer
 - Service Account User
 - Cloud Life Sciences Workflows Runner
 - Service Usage Consumer


Step 1: Upload Your Data
::::::::::::::::::::::::

We will be obtaining inputs from Google Cloud Storage, as well as saving
outputs there. You should first clone the repository with the Snakemake tutorial data:


.. code:: console

    git clone https://github.com/snakemake/snakemake-lsh-tutorial-data
    cd snakemake-lsh-tutorial-data


And then either manually create a bucket and upload data files there, or
use the `provided script and instructions <https://github.com/snakemake/snakemake-lsh-tutorial-data#google-cloud-storage>`_
to do it programatically from the command line. The script generally works like:

.. code:: console

    python upload_google_storage.py <bucket>/<subpath>   <folder>

And you aren't required to provide a subfolder path if you want to upload
to the root of a bucket. As an example, for this tutorial we upload the contents of
"data" to the root of the bucket `snakemake-testing-data`

.. code:: console

    export GOOGLE_APPLICATION_CREDENTIALS="/path/to/credentials.json"
    python upload_google_storage.py snakemake-testing-data data/

If you wanted to upload to a "subfolder" path in a bucket, you would do that as follows:

.. code:: console

    export GOOGLE_APPLICATION_CREDENTIALS="/path/to/credentials.json"
    python upload_google_storage.py snakemake-testing-data/subfolder data/

Your bucket (and the folder prefix) will be referred to as the
`--default-remote-prefix` when you run snakemake. You can visually
browse your data in the `storage browser <https://console.cloud.google.com/storage/>_`.


.. image:: workflow/upload-google-storage.png


Step 2: Write your Snakefile, Environment File, and Scripts
:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Now that we've exported our credentials and have all dependencies installed, let's
get our workflow! This is the exact same workflow from the :ref:`basic tutorial<tutorial-basics>`,
so if you need a refresher on the design or basics, please see those pages.
You can find the Snakefile, supporting scripts for plotting and environment in the `snakemake-lsh-tutorial-data <https://github.com/snakemake/snakemake-lsh-tutorial-data>`_ repository.

First, how does a working directory work for this executor? The present
working directory, as identified by Snakemake that has the Snakefile, and where
a more advanced setup might have a folder of environment specifications (env) a folder of scripts 
(scripts), and rules (rules), is considered within the context of the build.
When the Google Life Sciences executor is used, it generates a build package of all
of the files here (within a reasonable size) and uploads those to storage. This
package includes the .snakemake folder that would have been generated locally.
The build package is then downloaded and extracted by each cloud executor, which
is a Google Compute instance.

We next need an `environment.yaml` file that will define the dependencies
that we want installed with conda for our job. If you cloned the "snakemake-lsh-tutorial-data"
repository you will already have this, and you are good to go. If not, save this to `environment.yaml`
in your working directory:

.. code:: yaml

    channels:
      - conda-forge
      - bioconda
    dependencies:
      - python =3.6
      - jinja2 =2.10
      - networkx =2.1
      - matplotlib =2.2.3
      - graphviz =2.38.0
      - bcftools =1.9
      - samtools =1.9
      - bwa =0.7.17
      - pysam =0.15.0
    

Notice that we reference this `environment.yaml` file in the Snakefile below.
Importantly, if you were optimizing a pipeline, you would likely have a folder
"envs" with more than one environment specification, one for each step.
This workflow uses the same environment (with many dependencies) instead of
this strategy to minimize the number of files for you.

The Snakefile (also included in the repository) then has the following content. It's important to note
that we have not customized this file from the basic tutorial to hard code 
any storage. We will be telling snakemake to use the remote bucket as 
storage instead of the local filesystem.

.. code:: python

    SAMPLES = ["A", "B"]

    rule all:
        input:
            "plots/quals.svg"

    rule bwa_map:
        input:
            fastq="samples/{sample}.fastq",
            idx=multiext("genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
        conda:
            "environment.yaml"
        output:
            "mapped_reads/{sample}.bam"
        params:
            idx=lambda w, input: os.path.splitext(input.idx[0])[0]
        shell:
            "bwa mem {params.idx} {input.fastq} | samtools view -Sb - > {output}"

    rule samtools_sort:
        input:
            "mapped_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam"
        conda:
            "environment.yaml"
        shell:
            "samtools sort -T sorted_reads/{wildcards.sample} "
            "-O bam {input} > {output}"

    rule samtools_index:
        input:
            "sorted_reads/{sample}.bam"
        output:
            "sorted_reads/{sample}.bam.bai"
        conda:
            "environment.yaml"
        shell:
            "samtools index {input}"

    rule bcftools_call:
        input:
            fa="genome.fa",
            bam=expand("sorted_reads/{sample}.bam", sample=SAMPLES),
            bai=expand("sorted_reads/{sample}.bam.bai", sample=SAMPLES)
        output:
            "calls/all.vcf"
        conda:
            "environment.yaml"
        shell:
            "samtools mpileup -g -f {input.fa} {input.bam} | "
            "bcftools call -mv - > {output}"

    rule plot_quals:
        input:
            "calls/all.vcf"
        output:
            "plots/quals.svg"
        conda:
            "environment.yaml"
        script:
            "plot-quals.py"



And make sure you also have the script `plot-quals.py` in your present working directory for the last step.
This script will help us to do the plotting, and is also included in the `snakemake-lsh-tutorial-data <https://github.com/snakemake/snakemake-lsh-tutorial-data>`_ repository.


.. code:: python

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from pysam import VariantFile

    quals = [record.qual for record in VariantFile(snakemake.input[0])]
    plt.hist(quals)

    plt.savefig(snakemake.output[0])


Step 3: Run Snakemake
:::::::::::::::::::::

Now let's run Snakemake with the Google Life Sciences Executor.


.. code:: console

    snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1


The flags above refer to:

 - `--google-lifesciences`: to indicate that we want to use the Google Life Sciences API
 - `--default-remote-prefix`: refers to the Google Storage bucket. The bucket name is "snakemake-testing-data" and the "subfolder" (or path) (not defined above) would be a subfolder, if needed.
 - `--google-lifesciences-region`: the region that you want the instances to deploy to. Your storage bucket should be accessible from here, and your selection can have a small influence on the machine type selected.


Once you submit the job, you'll immediately see the familiar Snakemake console output,
but with additional lines for inspecting google compute instances with gcloud:

.. code:: console

    Building DAG of jobs...
    Unable to retrieve additional files from git. This is not a git repository.
    Using shell: /bin/bash
    Rules claiming more threads will be scaled down.
    Job counts:
    	count	jobs
    	1	all
    	1	bcftools_call
    	2	bwa_map
	1	plot_quals
	2	samtools_index
	2	samtools_sort
	9

    [Thu Apr 16 19:16:24 2020]
    rule bwa_map:
        input: snakemake-testing-data/genome.fa, snakemake-testing-data/samples/B.fastq
        output: snakemake-testing-data/mapped_reads/B.bam
        jobid: 8
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe 13586583122112209762
    gcloud beta lifesciences operations list


Take note of those last three lines to describe and list operations - this is how you
get complete error and output logs for the run, which we will demonstrate using later.


And you'll see a block like that for each rule. Here is what the entire workflow looks
like after completion:

.. code:: console

    Building DAG of jobs...
    Unable to retrieve additional files from git. This is not a git repository.
    Using shell: /bin/bash
    Rules claiming more threads will be scaled down.
    Job counts:
    	count	jobs
   	1	all
	1	bcftools_call
	2	bwa_map
	1	plot_quals
	2	samtools_index
	2	samtools_sort
	9

    [Fri Apr 17 20:27:51 2020]
    rule bwa_map:
        input: snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa
        output: snakemake-testing-data/mapped_reads/B.bam
        jobid: 8
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/16135317625786219242
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:31:16 2020]
    Finished job 8.
    1 of 9 steps (11%) done

    [Fri Apr 17 20:31:16 2020]
    rule bwa_map:
        input: snakemake-testing-data/samples/A.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa
        output: snakemake-testing-data/mapped_reads/A.bam
        jobid: 7
        wildcards: sample=A
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/5458247376121133509
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:34:30 2020]
    Finished job 7.
    2 of 9 steps (22%) done

    [Fri Apr 17 20:34:30 2020]
    rule samtools_sort:
        input: snakemake-testing-data/mapped_reads/B.bam
        output: snakemake-testing-data/sorted_reads/B.bam
        jobid: 4
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/13750029425473765929
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:37:34 2020]
    Finished job 4.
    3 of 9 steps (33%) done

    [Fri Apr 17 20:37:35 2020]
    rule samtools_sort:
        input: snakemake-testing-data/mapped_reads/A.bam
        output: snakemake-testing-data/sorted_reads/A.bam
        jobid: 3
        wildcards: sample=A
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/15643873965497084056
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:40:37 2020]
    Finished job 3.
    4 of 9 steps (44%) done

    [Fri Apr 17 20:40:38 2020]
    rule samtools_index:
        input: snakemake-testing-data/sorted_reads/B.bam
        output: snakemake-testing-data/sorted_reads/B.bam.bai
        jobid: 6
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/6525320566174651173
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:43:41 2020]
    Finished job 6.
    5 of 9 steps (56%) done

    [Fri Apr 17 20:43:41 2020]
    rule samtools_index:
        input: snakemake-testing-data/sorted_reads/A.bam
        output: snakemake-testing-data/sorted_reads/A.bam.bai
        jobid: 5
        wildcards: sample=A
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/9175497885319251567
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:46:44 2020]
    Finished job 5.
    6 of 9 steps (67%) done

    [Fri Apr 17 20:46:44 2020]
    rule bcftools_call:
        input: snakemake-testing-data/genome.fa, snakemake-testing-data/sorted_reads/A.bam, snakemake-testing-data/sorted_reads/B.bam, snakemake-testing-data/sorted_reads/A.bam.bai, snakemake-testing-data/sorted_reads/B.bam.bai
        output: snakemake-testing-data/calls/all.vcf
        jobid: 2
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/622600526583374352
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:49:57 2020]
    Finished job 2.
    7 of 9 steps (78%) done

    [Fri Apr 17 20:49:57 2020]
    rule plot_quals:
        input: snakemake-testing-data/calls/all.vcf
        output: snakemake-testing-data/plots/quals.svg
        jobid: 1
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us-west2/operations/9350722561866518561
    gcloud beta lifesciences operations list
    [Fri Apr 17 20:53:10 2020]
    Finished job 1.
    8 of 9 steps (89%) done

    [Fri Apr 17 20:53:10 2020]
    localrule all:
        input: snakemake-testing-data/plots/quals.svg
        jobid: 0
        resources: mem_mb=15360, disk_mb=128000

    Downloading from remote: snakemake-testing-data/plots/quals.svg
    Finished download.
    [Fri Apr 17 20:53:10 2020]
    Finished job 0.
    9 of 9 steps (100%) done
    Complete log: /home/vanessa/snakemake-work/tutorial/.snakemake/log/2020-04-17T202749.218820.snakemake.log


We've finished the run, great! Let's inspect our results.

Step 4: View Results
::::::::::::::::::::

The entirety of the log that was printed to the terminal will be available
on your local machine where you submit the job in the hidden `.snakemake`
folder under "log" and timestamped accordingly. If you look at the last line
in the output above, you'll see the full path to this file.

You also might notice a line about downloading results:

.. code:: console

    Downloading from remote: snakemake-testing-data/plots/quals.svg


Since we defined this to be the target of our run

.. code:: console


    rule all:
        input:
            "plots/quals.svg"


this fill is downloaded to our host too. Actually, you'll notice
that paths in storage are mirrored on your filesystem (this is what the workers
do too):


.. code:: console

    $ tree snakemake-testing-data/
    snakemake-testing-data/
    └── plots
        └── quals.svg


We can see the result of our run, quals.svg, below:

.. image:: workflow/quals.svg


And if we look at the remote storage, we see that the result file (under plots) and intermediate
results (under sorted_reads and calls) are saved there too!

.. image:: workflow/results-google-storage.png

The source folder contains a cache folder with archives that contain your working directories
that are extracted on the worker instances. You can safely delete this folder, or keep it if you want to reproduce
the run in the future.


Step 5: Debugging
:::::::::::::::::

Let's introduce an error (purposefully) into our Snakefile to practice debugging.
Let's remove the conda environment.yaml file for the first rule, so we would
expect that Snakemake won't be able to find the executables for bwa and samtools.
In your Snakefile, change this:

.. code:: python

    rule bwa_map:
        input:
            fastq="samples/{sample}.fastq",
            idx=multiext("genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
        conda:
            "environment.yaml"
        output:
            "mapped_reads/{sample}.bam"
        params:
            idx=lambda w, input: os.path.splitext(input.idx[0])[0]
        shell:
            "bwa mem {params.idx} {input.fastq} | samtools view -Sb - > {output}"


to this:

.. code:: python

    rule bwa_map:
        input:
            fastq="samples/{sample}.fastq",
            idx=multiext("genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
        output:
            "mapped_reads/{sample}.bam"
        params:
            idx=lambda w, input: os.path.splitext(input.idx[0])[0]
        shell:
            "bwa mem {params.idx} {input.fastq} | samtools view -Sb - > {output}"


And then for the same command to run everything again, you would need to remove the 
plots, mapped_reads, and calls folders. Instead, we can make this request more easily
by adding the argument `--forceall`:

.. code:: console

    snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1 --forceall

Everything will start out okay as it did before, and it will pause on the first 
step when it's deploying the first container image. The last part of the 
log will look somethig like this:


.. code:: console

    [Fri Apr 17 22:01:38 2020]
    rule bwa_map:
        input: snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa
        output: snakemake-testing-data/mapped_reads/B.bam
        jobid: 8
        wildcards: sample=B
        resources: mem_mb=15360, disk_mb=128000

    Get status with:
    gcloud config set project snakemake-testing
    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us/operations/11698975339184312706
    gcloud beta lifesciences operations list


Since we removed an important dependency to install libraries with conda, 
we are definitely going to hit an error! That looks like this:

.. code:: console

    [Fri Apr 17 22:03:08 2020]
    Error in rule bwa_map:
        jobid: 8
        output: snakemake-testing-data/mapped_reads/B.bam
        shell:
            bwa mem snakemake-testing-data/genome.fa snakemake-testing-data/samples/B.fastq | samtools view -Sb - > snakemake-testing-data/mapped_reads/B.bam
            (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)
        jobid: 11698975339184312706

    Shutting down, this might take some time.


Oh no! How do we debug it? The error above just indicates that "one of the commands
exised with a non-zero exit code," and that isn't really enough to know what happened,
and how to fix it. Debugging is actually quite simple, we can copy paste the gcloud
command to describe our operation into the console. This will spit out an entire structure
that shows every step of the rule running, from pulling a container, to downloading
the working directory, to running the step.

.. code:: console

    gcloud beta lifesciences operations describe projects/snakemake-testing/locations/us/operations/11698975339184312706
    done: true
    error:
      code: 9
      message: 'Execution failed: generic::failed_precondition: while running "snakejob-bwa_map-8":
        unexpected exit status 1 was not ignored'
    metadata:
      '@type': type.googleapis.com/google.cloud.lifesciences.v2beta.Metadata
      createTime: '2020-04-17T22:01:39.642966Z'
      endTime: '2020-04-17T22:02:59.149914114Z'
      events:
      - description: Worker released
        timestamp: '2020-04-17T22:02:59.149914114Z'
        workerReleased:
          instance: google-pipelines-worker-b1cdd36c743c3b477af8114d2511e37e
          zone: us-west1-c
      - description: 'Execution failed: generic::failed_precondition: while running "snakejob-bwa_map-8":
          unexpected exit status 1 was not ignored'
        failed:
          cause: 'Execution failed: generic::failed_precondition: while running "snakejob-bwa_map-8":
            unexpected exit status 1 was not ignored'
          code: FAILED_PRECONDITION
        timestamp: '2020-04-17T22:02:57.950752682Z'
      - description: Unexpected exit status 1 while running "snakejob-bwa_map-8"
        timestamp: '2020-04-17T22:02:57.842529458Z'
        unexpectedExitStatus:
          actionId: 1
          exitStatus: 1
      - containerStopped:
          actionId: 1
          exitStatus: 1
          stderr: |
            me.fa.bwt
            Finished download.
            /bin/bash: bwa: command not found
            /bin/bash: samtools: command not found
            [Fri Apr 17 22:02:57 2020]
            Error in rule bwa_map:
                jobid: 0
                output: snakemake-testing-data/mapped_reads/B.bam
                shell:
                    bwa mem snakemake-testing-data/genome.fa snakemake-testing-data/samples/B.fastq | samtools view -Sb - > snakemake-testing-data/mapped_reads/B.bam
                    (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

            Removing output files of failed job bwa_map since they might be corrupted:
            snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa, snakemake-testing-data/mapped_reads/B.bam
            Shutting down, this might take some time.
            Exiting because a job execution failed. Look above for error message
            Complete log: /workdir/.snakemake/log/2020-04-17T220254.129519.snakemake.log
        description: |-
          Stopped running "snakejob-bwa_map-8": exit status 1: me.fa.bwt
          Finished download.
          /bin/bash: bwa: command not found
          /bin/bash: samtools: command not found
          [Fri Apr 17 22:02:57 2020]
          Error in rule bwa_map:
              jobid: 0
              output: snakemake-testing-data/mapped_reads/B.bam
              shell:
                  bwa mem snakemake-testing-data/genome.fa snakemake-testing-data/samples/B.fastq | samtools view -Sb - > snakemake-testing-data/mapped_reads/B.bam
                  (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

          Removing output files of failed job bwa_map since they might be corrupted:
          snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa, snakemake-testing-data/mapped_reads/B.bam
          Shutting down, this might take some time.
          Exiting because a job execution failed. Look above for error message
          Complete log: /workdir/.snakemake/log/2020-04-17T220254.129519.snakemake.log
        timestamp: '2020-04-17T22:02:57.842442588Z'
      - containerStarted:
          actionId: 1
        description: Started running "snakejob-bwa_map-8"
        timestamp: '2020-04-17T22:02:51.724433437Z'
      - description: Stopped pulling "snakemake/snakemake:v5.10.0"
        pullStopped:
          imageUri: snakemake/snakemake:v5.10.0
        timestamp: '2020-04-17T22:02:43.696978950Z'
      - description: Started pulling "snakemake/snakemake:v5.10.0"
        pullStarted:
          imageUri: snakemake/snakemake:v5.10.0
        timestamp: '2020-04-17T22:02:10.339950219Z'
      - description: Worker "google-pipelines-worker-b1cdd36c743c3b477af8114d2511e37e"
          assigned in "us-west1-c"
        timestamp: '2020-04-17T22:01:43.232858222Z'
        workerAssigned:
          instance: google-pipelines-worker-b1cdd36c743c3b477af8114d2511e37e
          machineType: n2-highmem-2
          zone: us-west1-c
      labels:
        app: snakemake
        name: snakejob-b346c449-9fd6-4f1e-8043-17c300cc9c0d-bwa_map-8
      pipeline:
        actions:
        - commands:
          - /bin/bash
          - -c
          - 'mkdir -p /workdir && cd /workdir && wget -O /download.py https://gist.githubusercontent.com/vsoch/84886ef6469bedeeb9a79a4eb7aec0d1/raw/181499f8f17163dcb2f89822079938cbfbd258cc/download.py
            && chmod +x /download.py && source activate snakemake || true && pip install
            crc32c && python /download.py download snakemake-testing-data source/cache/snakeworkdir-5f4f325b9ddb188d5da8bfab49d915f023509c0b1986eb72cb4a2540d7991c12.tar.gz
            /tmp/workdir.tar.gz && tar -xzvf /tmp/workdir.tar.gz && snakemake snakemake-testing-data/mapped_reads/B.bam
            --snakefile Snakefile --force -j --keep-target-files --keep-remote --latency-wait
            0 --attempt 1 --force-use-threads  --allowed-rules bwa_map --nocolor --notemp
            --no-hooks --nolock  --use-conda  --default-remote-provider GS --default-remote-prefix
            snakemake-testing-data  --default-resources "mem_mb=15360" "disk_mb=128000" '
          containerName: snakejob-bwa_map-8
          imageUri: snakemake/snakemake:v5.10.0
          labels:
            app: snakemake
            name: snakejob-b346c449-9fd6-4f1e-8043-17c300cc9c0d-bwa_map-8
        resources:
          regions:
          - us-west1
          virtualMachine:
            bootDiskSizeGb: 135
            bootImage: projects/cos-cloud/global/images/family/cos-stable
            labels:
              app: snakemake
              goog-pipelines-worker: 'true'
            machineType: n2-highmem-2
            serviceAccount:
              email: default
              scopes:
              - https://www.googleapis.com/auth/cloud-platform
        timeout: 604800s
      startTime: '2020-04-17T22:01:43.232858222Z'
    name: projects/411393320858/locations/us/operations/11698975339184312706


The log is hefty, so let's break it into pieces to talk about. Firstly, it's
intended to be read from the bottom up if you want to see things in chronological order.
The very bottom line is the unique id of the operation, and this is what you used 
(with the project identifier string, the number after projects, replaced with your project
name) to query for the log. Let's look at the next section, `pipeline`. This was
the specification built up by Snakemake and sent to the Google Life Sciences API
as a request:

.. code:: console

      pipeline:
        actions:
        - commands:
          - /bin/bash
          - -c
          - 'mkdir -p /workdir && cd /workdir && wget -O /download.py https://gist.githubusercontent.com/vsoch/84886ef6469bedeeb9a79a4eb7aec0d1/raw/181499f8f17163dcb2f89822079938cbfbd258cc/download.py
            && chmod +x /download.py && source activate snakemake || true && pip install
            crc32c && python /download.py download snakemake-testing-data source/cache/snakeworkdir-5f4f325b9ddb188d5da8bfab49d915f023509c0b1986eb72cb4a2540d7991c12.tar.gz
            /tmp/workdir.tar.gz && tar -xzvf /tmp/workdir.tar.gz && snakemake snakemake-testing-data/mapped_reads/B.bam
            --snakefile Snakefile --force -j --keep-target-files --keep-remote --latency-wait
            0 --attempt 1 --force-use-threads  --allowed-rules bwa_map --nocolor --notemp
            --no-hooks --nolock  --use-conda  --default-remote-provider GS --default-remote-prefix
            snakemake-testing-data  --default-resources "mem_mb=15360" "disk_mb=128000" '
          containerName: snakejob-bwa_map-8
          imageUri: snakemake/snakemake:v5.10.0
          labels:
            app: snakemake
            name: snakejob-b346c449-9fd6-4f1e-8043-17c300cc9c0d-bwa_map-8
        resources:
          regions:
          - us-west1
          virtualMachine:
            bootDiskSizeGb: 135
            bootImage: projects/cos-cloud/global/images/family/cos-stable
            labels:
              app: snakemake
              goog-pipelines-worker: 'true'
            machineType: n2-highmem-2
            serviceAccount:
              email: default
              scopes:
              - https://www.googleapis.com/auth/cloud-platform
        timeout: 604800s
      startTime: '2020-04-17T22:01:43.232858222Z'


There is a lot of useful information here. Under *resources*:

- **virtualMachine** shows the **machineType** that should correspond to the instance type. You can specify a full name or prefix with `--machine-type-prefix` or "machine_type" defined under resources for a step. Since we didn't set any requirements, it chose a reasonable choice for us. This section also shows the size of the boot disk (in GB) and if you added hardware accelerators (GPU) they should show up here too.
- **regions** is the region that the instance was deployed in, which is important to know if you need to specify to run from a particular region. This parameter defalts to regions in the US, and can be modified with the `--google-lifesciences-regions` parameter.

Under *actions* you'll find a few important fields:

- **imageUri** is important to know to see the version of Snakemake (or another container base) that was used. You can customize this with `--container-image`, and it will default to the latest snakemake.
- **commands** are the commands run to execute the container (also known as the entrypoint). For example, if you wanted to bring up your own instance manually and pull the container defined by `imageUri`, you could execute the commands to the container (or shell inside and run them interactively) to interactively debug. Notice that the commands ends with a call to snakemake, and shows the arguments that are used. Make sure that this matches your expectation.

The next set of steps pertain to assigning the worker, pulling the container, and starting it. 
That looks something like this, and it's fairly straight forward. You can again see
that earlier timestamps are on the bottom.

.. code:: console

      - containerStarted:
          actionId: 1
        description: Started running "snakejob-bwa_map-8"
        timestamp: '2020-04-17T22:02:51.724433437Z'
      - description: Stopped pulling "snakemake/snakemake:v5.10.0"
        pullStopped:
          imageUri: snakemake/snakemake:v5.10.0
        timestamp: '2020-04-17T22:02:43.696978950Z'
      - description: Started pulling "snakemake/snakemake:v5.10.0"
        pullStarted:
          imageUri: snakemake/snakemake:v5.10.0
        timestamp: '2020-04-17T22:02:10.339950219Z'
      - description: Worker "google-pipelines-worker-b1cdd36c743c3b477af8114d2511e37e"
          assigned in "us-west1-c"
        timestamp: '2020-04-17T22:01:43.232858222Z'
        workerAssigned:
          instance: google-pipelines-worker-b1cdd36c743c3b477af8114d2511e37e
          machineType: n2-highmem-2
          zone: us-west1-c


The next section, when the container is stopped, have the meat of the information
that we need to debug! This is the step where there was a non-zero exit code.

.. code:: console

      - containerStopped:
          actionId: 1
          exitStatus: 1
          stderr: |
            me.fa.bwt
            Finished download.
            /bin/bash: bwa: command not found
            /bin/bash: samtools: command not found
            [Fri Apr 17 22:02:57 2020]
            Error in rule bwa_map:
                jobid: 0
                output: snakemake-testing-data/mapped_reads/B.bam
                shell:
                    bwa mem snakemake-testing-data/genome.fa snakemake-testing-data/samples/B.fastq | samtools view -Sb - > snakemake-testing-data/mapped_reads/B.bam
                    (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

            Removing output files of failed job bwa_map since they might be corrupted:
            snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa, snakemake-testing-data/mapped_reads/B.bam
            Shutting down, this might take some time.
            Exiting because a job execution failed. Look above for error message
            Complete log: /workdir/.snakemake/log/2020-04-17T220254.129519.snakemake.log
        description: |-
          Stopped running "snakejob-bwa_map-8": exit status 1: me.fa.bwt
          Finished download.
          /bin/bash: bwa: command not found
          /bin/bash: samtools: command not found
          [Fri Apr 17 22:02:57 2020]
          Error in rule bwa_map:
              jobid: 0
              output: snakemake-testing-data/mapped_reads/B.bam
              shell:
                  bwa mem snakemake-testing-data/genome.fa snakemake-testing-data/samples/B.fastq | samtools view -Sb - > snakemake-testing-data/mapped_reads/B.bam
                  (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

          Removing output files of failed job bwa_map since they might be corrupted:
          snakemake-testing-data/samples/B.fastq, snakemake-testing-data/genome.fa.amb, snakemake-testing-data/genome.fa.ann, snakemake-testing-data/genome.fa.bwt, snakemake-testing-data/genome.fa.pac, snakemake-testing-data/genome.fa.sa, snakemake-testing-data/mapped_reads/B.bam
          Shutting down, this might take some time.
          Exiting because a job execution failed. Look above for error message
          Complete log: /workdir/.snakemake/log/2020-04-17T220254.129519.snakemake.log
        timestamp: '2020-04-17T22:02:57.842442588Z'


Along with seeing the error in `stderr`, the description key holds the same error. We see
what we would have seen if we were running the bwa mem command on our own command line,
that the executables weren't found:

.. code:: console

      stderr: |
        me.fa.bwt
        Finished download.
        /bin/bash: bwa: command not found
        /bin/bash: samtools: command not found


But we shouldn't be surprised, we on purpose removed the environment file to install
them! This is where you would read the error, and respond by updating your Snakefile with
a fix. 


Step 6: Adding a Log File
:::::::::::::::::::::::::

How might we do better at debugging in the future? The answer is to 
add a log file for each step, which is where any stderr will be written 
in the case of failure. For the same step above, we would update the rule
to look like this:


.. code:: python

    rule bwa_map:
        input:
            fastq="samples/{sample}.fastq",
            idx=multiext("genome.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
        output:
            "mapped_reads/{sample}.bam"
        params:
            idx=lambda w, input: os.path.splitext(input.idx[0])[0]
        shell:
            "bwa mem {params.idx} {input.fastq} | samtools view -Sb - > {output}"
        log:
            "logs/bwa_map/{sample}.log" 


In the above, we would write a log file to storage in a "subfolder" of the
snakemake prefix located at "logs/bwa_map." The log file will be named according
to the sample. You could also imagine a flatted structure with a path like
`logs/bwa_map-{sample}.log`. It's up to you how you want to organize your output.
This means that when you see the error appear in your terminal, you can quickly
look at this log file instead of resorting to using the gcloud tool. It's generally
good to remember when debugging that:

 - You should not make assumptions about anything's existence. Use print statements to verify.
 - The biggest errors tend to be syntax and/or path errors
 - If you want to test a different snakemake container, you can use the `--container` flag.
 - If the error is especially challenging, set up a small toy example that implements the most basic functionality that you want to achieve.
 - If you need help, reach out to ask for it! If there is an issue with the Google Life Sciences workflow executor, please `open an issue <https://github.com/snakemake/snakemake/issues>`_.
 - It also sometimes helps to take a break from working on something, and coming back with fresh eyes.
.. _executor_tutorial:

============================
Snakemake Executor Tutorials
============================

.. _cloud executors: https://snakemake.readthedocs.io/en/stable/executing/cloud.html
.. _tutorial: https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html

This set of tutorials are intended to introduce you to executing `cloud executors`_.
We start with the original Snakemake `tutorial`_ and expand upon it to be run
in different cloud environments. For each run, we show you how to:

 - authenticate with credentials, if required
 - prepare your workspace
 - submit a basic job
 - generate an error and debug

The examples presented in these tutorials come from Bioinformatics.
However, Snakemake is a general-purpose workflow management system for any discipline.
We ensured that no bioinformatics knowledge is needed to understand the tutorial.

.. toctree::
   :maxdepth: 2

   google_lifesciences

   azure_aks

   
.. _job_grouping:

============
Job Grouping
============

The graph of jobs that Snakemake determines before execution can be partitioned into groups.
Such groups will be executed together in **cluster** or **cloud mode**, as a so-called **group job**, i.e., all jobs of a particular group will be submitted at once, to the same computing node.
When executing locally, group definitions are ignored.

Groups can be defined along with the workflow definition via the ``group`` keyword, see :ref:`snakefiles-grouping`.
This way, queueing and execution time can be saved, in particular by attaching short-running downstream jobs to long running upstream jobs.

However, often the benefit grouping should be heavily dependent on the specifics of the underlying computing platform.
Hence, it is possible to assign groups via the command line.
For example, with

.. code-block:: bash

    snakemake --groups somerule=group0 someotherrule=group0

we assign the two rules ``somerule`` and ``someotherrule`` to the same group ``group0``.
by default, groups do not span disconnected parts of the DAG.
Here, this means that by default only jobs of ``somerule`` and ``someotherrule`` end in the same group that are directly connected.
It is however possible to configure the number of connected DAG components that are spanned by a group via the flag ``--group-components``.
This way, it is e.g. possible to define batches of jobs of the same kind that shall be executed within one group, e.g.


.. code-block:: bash

    snakemake --groups somerule=group0 --group-components group0=5

means that given that there exist ``n`` jobs spawned from rule ``somerule``, Snakemake will create ``n / 5`` groups which each execute 5 jobs of ``somerule`` together.
For example, with 10 jobs from ``somerule`` you would end up with 2 groups of 5 jobs that are submitted as one piece each.

================
Interoperability
================

.. _cwl_export:

----------
CWL export
----------

Snakemake workflows can be exported to `CWL <https://www.commonwl.org/>`_, such that they can be executed in any `CWL-enabled workflow engine <https://www.commonwl.org/#Implementations>`_.
Since, CWL is less powerful for expressing workflows than Snakemake (most importantly Snakemake offers more flexible scatter-gather patterns, since full Python can be used), export works such that every Snakemake job is encoded into a single step in the CWL workflow.
Moreover, every step of that workflow calls Snakemake again to execute the job. The latter enables advanced Snakemake features like scripts, benchmarks and remote files to work inside CWL.
So, when exporting keep in mind that the resulting CWL file can become huge, depending on the number of jobs in your workflow.
To export a Snakemake workflow to CWL, simply run

.. code-block:: console

    $ snakemake --export-cwl workflow.cwl

The resulting workflow will by default use the `Snakemake docker image <https://hub.docker.com/r/snakemake/snakemake>`_ for every step, but this behavior can be overwritten via the CWL execution environment.
Then, the workflow can be executed in the same working directory with, e.g.,

.. code-block:: console

    $ cwltool workflow.cwl

Note that due to limitations in CWL, it seems currently impossible to avoid that all target files (output files of target jobs), are written directly to the workdir, regardless of their relative paths in the Snakefile.

Note that export is impossible in case the workflow contains :ref:`checkpoints <snakefiles-checkpoints>` or output files with absolute paths... _executable:

======================
Command line interface
======================

This part of the documentation describes the ``snakemake`` executable.  Snakemake
is primarily a command-line tool, so the ``snakemake`` executable is the primary way
to execute, debug, and visualize workflows.

.. user_manual-snakemake_options:

-----------------------------
Useful Command Line Arguments
-----------------------------

If called with the number of cores to use, i.e.

.. code-block:: console

    $ snakemake --cores 1

Snakemake tries to execute the workflow specified in a file called ``Snakefile`` in the same directory (the Snakefile can be given via the parameter ``-s``).

By issuing

.. code-block:: console

    $ snakemake -n

a dry-run can be performed.
This is useful to test if the workflow is defined properly and to estimate the amount of needed computation.
Further, the reason for each rule execution can be printed via


.. code-block:: console

    $ snakemake -n -r

Importantly, Snakemake can automatically determine which parts of the workflow can be run in parallel.
By specifying more than one available core, i.e.

.. code-block:: console

    $ snakemake --cores 4

one can tell Snakemake to use up to 4 cores and solve a binary knapsack problem to optimize the scheduling of jobs.
If the number is omitted (i.e., only ``--cores`` is given), the number of used cores is determined as the number of available CPU cores in the machine.

Snakemake workflows usually define the number of used threads of certain rules. Sometimes, it makes sense to overwrite the defaults given in the workflow definition.
This can be done by using the ``--set-threads`` argument, e.g.,

.. code-block:: console

    $ snakemake --cores 4 --set-threads myrule=2

would overwrite whatever number of threads has been defined for the rule ``myrule`` and use ``2`` instead.
Similarly, it is possible to overwrite other resource definitions in rules, via

.. code-block:: console

    $ snakemake --cores 4 --set-resources myrule:partition="foo"

Both mechanisms can be particularly handy when used in combination with :ref:`cluster execution <cluster>`.

Dealing with very large workflows
---------------------------------

If your workflow has a lot of jobs, Snakemake might need some time to infer the dependencies (the job DAG) and which jobs are actually required to run.
The major bottleneck involved is the filesystem, which has to be queried for existence and modification dates of files.
To overcome this issue, Snakemake allows to run large workflows in batches.
This way, fewer files have to be evaluated at once, and therefore the job DAG can be inferred faster.
By running

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=1/3

you instruct to only compute the first of three batches of the inputs of the rule `myrule`.
To generate the second batch, run

.. code-block:: console

    $ snakemake --cores 4 --batch myrule=2/3

Finally, when running


.. code-block:: console

    $ snakemake --cores 4 --batch myrule=3/3

Snakemake will process beyond the rule `myrule`, because all of its input files have been generated, and complete the workflow.
Obviously, a good choice of the rule to perform the batching is a rule that has a lot of input files and upstream jobs, for example a central aggregation step within your workflow.
We advice all workflow developers to inform potential users of the best suited batching rule.

.. _profiles:

--------
Profiles
--------

Adapting Snakemake to a particular environment can entail many flags and options.
Therefore, since Snakemake 4.1, it is possible to specify a configuration profile
to be used to obtain default options:

.. code-block:: console

   $ snakemake --profile myprofile

Here, a folder ``myprofile`` is searched in per-user and global configuration directories (on Linux, this will be ``$HOME/.config/snakemake`` and ``/etc/xdg/snakemake``, you can find the answer for your system via ``snakemake --help``).
Alternatively, an absolute or relative path to the folder can be given.
The profile folder is expected to contain a file ``config.yaml`` that defines default values for the Snakemake command line arguments.
For example, the file

.. code-block:: yaml

    cluster: qsub
    jobs: 100

would setup Snakemake to always submit to the cluster via the ``qsub`` command, and never use more than 100 parallel jobs in total.
The profile can be used to set a default for each option of the Snakemake command line interface.
For this, option ``--someoption`` becomes ``someoption:`` in the profile.
If options accept multiple arguments these must be given as YAML list in the profile.
Under https://github.com/snakemake-profiles/doc, you can find publicly available profiles.
Feel free to contribute your own.

The profile folder can additionally contain auxilliary files, e.g., jobscripts, or any kind of wrappers.
See https://github.com/snakemake-profiles/doc for examples.


.. _getting_started-visualization:

-------------
Visualization
-------------

To visualize the workflow, one can use the option ``--dag``.
This creates a representation of the DAG in the graphviz dot language which has to be postprocessed by the graphviz tool ``dot``.
E.g. to visualize the DAG that would be executed, you can issue:

.. code-block:: console

    $ snakemake --dag | dot | display

For saving this to a file, you can specify the desired format:

.. code-block:: console

    $ snakemake --dag | dot -Tpdf > dag.pdf

To visualize the whole DAG regardless of the eventual presence of files, the ``forceall`` option can be used:

.. code-block:: console

    $ snakemake --forceall --dag | dot -Tpdf > dag.pdf

Of course the visual appearance can be modified by providing further command line arguments to ``dot``.

**Note:** The DAG is printed in DOT format straight to the standard output, along with other ``print`` statements you may have in your Snakefile. Make sure to comment these other ``print`` statements so that ``dot`` can build a visual representation of your DAG.


.. _all_options:

-----------
All Options
-----------

.. argparse::
   :module: snakemake
   :func: get_argument_parser
   :prog: snakemake

   All command line options can be printed by calling ``snakemake -h``.

.. _getting_started-bash_completion:

---------------
Bash Completion
---------------

Snakemake supports bash completion for filenames, rulenames and arguments.
To enable it globally, just append

.. code-block:: bash

    `snakemake --bash-completion`

including the backticks to your ``.bashrc``.
This only works if the ``snakemake`` command is in your path.
.. _caching:

========================
Between workflow caching
========================

Within certain data analysis fields, there are certain intermediate results that reoccur in exactly the same way in many analysis.
For example, in bioinformatics, reference genomes or annotations are downloaded, and read mapping indexes are built.
Since such steps are independent of the actual data or measurements that are analyzed, but still computationally or timely expensive to conduct, it has been common practice to externalize their computation and assume the presence of the resulting files before execution of a workflow.

From version 5.8.0 on, Snakemake offers a way to keep those steps inside the actual analysis without requiring from redundant computations.
By hashing all steps, parameters, software stacks (in terms of conda environments or containers), and raw input required up to a certain step in a `blockchain <https://en.wikipedia.org/wiki/Blockchain>`_, Snakemake is able to recognize **before** the computation whether a certain result is already available in a central cache on the same system.
**Note that this is explicitly intended for caching results between workflows! There is no need to use this feature to avoid redundant computations within a workflow. Snakemake does this already out of the box.**

Such caching has to be explicitly activated per rule, which can be done via the command line interface.
For example,

.. code-block:: console

    $ export SNAKEMAKE_OUTPUT_CACHE=/mnt/snakemake-cache/
    $ snakemake --cache download_data create_index

would instruct Snakemake to cache and reuse the results of the rules ``download_data`` and ``create_index``.
The environment variable definition that happens in the first line (defining the location of the cache) should of course be done only once and system wide in reality.
When Snakemake is executed without a shared filesystem (e.g., in the cloud, see :ref:`cloud`), the environment variable has to point to a location compatible with the given remote provider (e.g. an S3 or Google Storage bucket).
In any case, the provided location should be shared between all workflows of your group, institute or computing environment, in order to benefit from the reuse of previously obtained intermediate results.

Alternatively, rules can be marked as eligible for caching via the ``cache`` directive:

.. code-block:: python

    rule download_data:
        output:
            "results/data/worldcitiespop.csv"
        cache: True
        shell:
            "curl -L https://burntsushi.net/stuff/worldcitiespop.csv > {output}"

For workflows defining cache rules like this, it is enough to invoke Snakemake with

.. code-block:: console

    $ snakemake --cache

without explicit rulenames listed.

Note that only rules with just a single output file (or directory) or with :ref:`multiext output files <snakefiles-multiext>` are eligible for caching.
The reason is that for other rules it would be impossible to unambiguously assign the output files to cache entries while being agnostic of the actual file names.
Also note that the rules need to retrieve all their parameters via the ``params`` directive (except input files).
It is not allowed to directly use ``wildcards``, ``config`` or any global variable in the shell command or script, because these are not captured in the hash (otherwise, reuse would be unnecessarily limited).

Also note that Snakemake will store everything in the cache as readable and writeable for **all users** on the system (except in the remote case, where permissions are not enforced and depend on your storage configuration).
Hence, caching is not intended for private data, just for steps that deal with publicly available resources.

Finally, be aware that the implementation should be considered **experimental** until this note is removed.
.. _monitoring:

==========
Monitoring
==========

Snakemake supports `panoptes <https://github.com/panoptes-organization/panoptes>`_ a server (under development) that lets you monitor the execution of snakemake workflows.
Snakemake communicates with panoptes via the :code:`--wms-monitor` flag. The flag specifies the ip and port where panoptes is running (e.g. :code:`--wms-monitor http://127.0.0.1:5000`).

For panoptes versions 0.1.1 and lower, Snakemake sends the following requests to wms monitor:

.. csv-table::
   :header: "API", "Method", "Data", "Description"
   :widths: 40, 20, 20, 60

   ":code:`/api/service-info`", "GET", "json", "Snakemake gets the status of panoptes. Snakemake continues to run if the status (:code:`json['status']`) is :code:`'running'`. In all other cases snakemake exits with an error message."
   ":code:`/create_workflow`", "GET", "json", "Snakemake gets a unique id/name :code:`str(uuid.uuid4())` for each workflow triggered."
   ":code:`/update_workflow_status`", "POST", "dictionary", "Snakemake posts updates for workflows/jobs. The dictionary sent contains the log message dictionary , the current timestamp and the unique id/name of the workflow.
   
    .. code:: python

        {
            'msg': repr(msg), 
            'timestamp': time.asctime(), 
            'id': id
        }"


For future versions, Panoptes will implement`a more structured schema <https://github.com/panoptes-organization/monitor-schema>`_
to interact with the server. This means that for Snakemake 3.30.1 and lower, you should use Panoptes 0.1.1 and lower.
The documentation here will be updated when a new version of Panoptes with the Monitor Schema is released.
.. _cloud:

===========================
Cloud execution
===========================



------------------------------------
Generic cloud support via Kubernetes
------------------------------------

Snakemake 4.0 and later supports execution in the cloud via Kubernetes.
This is independent of the cloud provider, but we provide the setup steps for GCE below.

Setup Kubernetes on Google cloud engine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, install the `Google Cloud SDK <https://cloud.google.com/sdk/docs/quickstarts>`_.
Then, run

.. code-block:: console

    $ gcloud init

to setup your access.
Then, you can create a new kubernetes cluster via

.. code-block:: console

    $ gcloud container clusters create $CLUSTER_NAME --num-nodes=$NODES --scopes storage-rw

with ``$CLUSTER_NAME`` being the cluster name and ``$NODES`` being the number of cluster
nodes. If you intend to use google storage, make sure that ``--scopes storage-rw`` is set.
This enables Snakemake to write to the google storage from within the cloud nodes.
Next, you configure Kubernetes to use the new cluster via

.. code-block:: console

    $ gcloud container clusters get-credentials $CLUSTER_NAME


If you are having issues with authentication, please refer to the help text:

.. code-block:: console

    $ gcloud container clusters get-credentials --help

You likely also want to use google storage for reading and writing files.
For this, you will additionally need to authenticate with your google cloud account via

.. code-block:: console

    $ gcloud auth application-default login

This enables Snakemake to access google storage in order to check existence and modification dates of files.
Now, Snakemake is ready to use your cluster.

**Important:** After finishing your work, do not forget to delete the cluster with

.. code-block:: console

    $ gcloud container clusters delete $CLUSTER_NAME

in order to avoid unnecessary charges.


.. _kubernetes:


Executing a Snakemake workflow via kubernetes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming that kubernetes has been properly configured (see above), you can
execute a workflow via:

.. code-block:: console

    snakemake --kubernetes --use-conda --default-remote-provider $REMOTE --default-remote-prefix $PREFIX

In this mode, Snakemake will assume all input and output files to be stored in a given
remote location, configured by setting ``$REMOTE`` to your provider of choice
(e.g. ``GS`` for Google cloud storage or ``S3`` for Amazon S3) and ``$PREFIX``
to a bucket name or subfolder within that remote storage.
After successful execution, you find your results in the specified remote storage.
Of course, if any input or output already defines a different remote location, the latter will be used instead.
Importantly, this means that Snakemake does **not** require a shared network
filesystem to work in the cloud.


.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.

Currently, this mode requires that the Snakemake workflow is stored in a git repository.
Snakemake uses git to query necessary source files (the Snakefile, scripts, config, ...)
for workflow execution and encodes them into the kubernetes job.
Importantly, this also means that you should not put large non-source files into the git repo, since Snakemake will try to upload them to kubernetes with every job.
With large files in the git repo, this can lead to performance issues or even random SSL errors from kubernetes.

It is further possible to forward arbitrary environment variables to the kubernetes
jobs via the flag ``--envvars`` (see ``snakemake --help``) or the ``envvars`` directive in the Snakefile.
The former should be used e.g. for platform specific variables (e.g. secrets that are only needed for your kubernetes setup), whereas the latter should be used for variables that are needed for the workflow itself, regardless of whether it is executed on kubernetes or with a different backend.

When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as ``mem_mb`` to kubernetes. Further, it will propagate the number of threads
a job intends to use, such that kubernetes can allocate it to the correct cloud
computing node.


-------------------------------------------------------------
Executing a Snakemake workflow via Google Cloud Life Sciences
-------------------------------------------------------------

The `Google Cloud Life Sciences <https://cloud.google.com/life-sciences/docs/>`_
provides a rich application programming interface to design pipelines.
You'll first need to `follow instructions here <https://cloud.google.com/life-sciences/docs/quickstart>`_  to
create a Google Cloud Project and enable Life Sciences, Storage, and Compute Engine APIs,
and continue with the prompts to create credentials. You'll want to create
a service account for your host (it's easiest to give project Owner permissions), 
and save the json credentials. You'll want to export the full path to this file to ``GOOGLE_APPLICATION_CREDENTIALS`` :

.. code-block:: console

      $ export GOOGLE_APPLICATION_CREDENTIALS=$HOME/path/snakemake-credentials.json

If you lose the link to the credentials interface, you can `find it here <https://console.cloud.google.com/apis/credentials>`_.

Optionally, you can export ``GOOGLE_CLOUD_PROJECT`` as the name of your Google Cloud Project. By default, the project associated with your application credentials will be used.

.. code-block:: console

      $ export GOOGLE_CLOUD_PROJECT=my-project-name


Data in Google Storage
~~~~~~~~~~~~~~~~~~~~~~

Using this executor typically requires you to start with large data files
already in Google Storage, and then interact with them via the Google Storage
remote executor. An easy way to do this is to use the
`gsutil <https://cloud.google.com/storage/docs/uploading-objects>`_
command line client. For example, here is how we might upload a file
to storage using it:

.. code-block:: console

    $ gsutil -m cp mydata.txt gs://snakemake-bucket/1/mydata.txt

The ``-m`` parameter enables multipart uploads for large files, so you
can remove it if you are uploading one or more smaller files.
And note that you'll need to modify the file and bucket names.
Note that you can also easily use the Google Cloud Console interface, if
a graphical interface is preferable to you.

Environment Variables
~~~~~~~~~~~~~~~~~~~~~

**Important:** Google Cloud Life Sciences uses Google Compute, and does
**not** encrypt environment variables. If you specify environment
variables with the envvars directive or ``--envvars`` they will **not** be secrets.

Container Bases
~~~~~~~~~~~~~~~

By default, Google Life Sciences uses the latest stable version of
`snakemake/snakemake <https://hub.docker.com/r/snakemake/snakemake/tags>`_
on Docker Hub. You can choose to modify the container base with
the ``--container-image`` (or ``container_image`` from within Python),
however if you do so, your container must meet the following requirements:

 - have an entrypoint that can execute a ``/bin/bash`` command
 - have snakemake installed, either via ``conda activate snakemake`` or already on the path
 - also include snakemake Python dependencies for google.cloud

If you use any Snakemake container as a base, you should be good to go. If you'd
like to get a reference for requirements, it's helpful to look at the
`Dockerfile <https://github.com/snakemake/snakemake/blob/main/Dockerfile>`_
for Snakemake.

Requesting GPUs
~~~~~~~~~~~~~~~

The Google Life Sciences API currently has support for 
`NVIDIA GPUs <https://cloud.google.com/compute/docs/gpus#restrictions>`_, meaning that you can request a number of NVIDIA GPUs explicitly by adding ``nvidia_gpu`` or ``gpu`` to your Snakefile resources for a step:

.. code-block:: python

    rule a:
        output:
            "test.txt"
        resources:
            nvidia_gpu=1
        shell:
            "somecommand ..."

A specific `gpu model <https://cloud.google.com/compute/docs/gpus#introduction>`_ can be requested using ``gpu_model`` and lowercase identifiers like ``nvidia-tesla-p100`` or ``nvidia-tesla-p4``, for example: ``gpu_model="nvidia-tesla-p100"``. If you don't specify ``gpu`` or ``nvidia_gpu`` with a count, but you do specify a ``gpu_model``, the count will default to 1.

In addition to GPU for the Google Lifesciences Executor, you can request a `Google Cloud preemptible virtual machine <https://cloud.google.com/life-sciences/docs/reference/gcloud-examples#using_preemptible_vms>`_ for one or more steps. See the `rules documentation <https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#preemptible-virtual-machine>`_ for how to add one or more preemptible arguments.


Machine Types
~~~~~~~~~~~~~

To specify an exact `machine type <https://cloud.google.com/compute/docs/machine-types>`_
or a prefix to filter down to and then select based on other resource needs, 
you can set a default resource on the command line, either for a prefix or 
a full machine type:

.. code-block:: console

    --default-resources "machine_type=n1-standard"


If you want to specify the machine type as a resource, you can do that too:

.. code-block:: python

    rule a:
        output:
            "test.txt"
        resources:
            machine_type="n1-standard-8"
        shell:
            "somecommand ..."


If you request a gpu, this requires the "n1" prefix and your preference from
the file or command line will be overridden. Note that the default resources
for Google Life Sciences (memory and disk) are the same as for Tibanna.

Running the Life Sciences Executor
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When your Snakefile is ready, you can run snakemake to specify the life
sciences executor. Notice that we are also providing a remote prefix for our storage path,
along with a region.

.. code-block:: console

    $ snakemake --google-lifesciences --default-remote-prefix snakemake-testing-data --use-conda --google-lifesciences-region us-west1


For more details and examples, we recommend you reference the 
`Google Life Sciences Executor Tutorial <https://snakemake.readthedocs.io/en/stable/executor_tutorial/google_lifesciences.html>`_.


-----------------------------------------------------------------
Executing a Snakemake workflow via Tibanna on Amazon Web Services
-----------------------------------------------------------------

First, install `Tibanna <https://tibanna.readthedocs.io/en/latest/>`_.

.. code-block:: console

    $ pip install -U tibanna


Set up aws configuration either by creating files ``~/.aws/credentials`` and ``~/.aws/config`` 
or by setting up environment variables as below (see Tibanna or AWS documentation for more details):

.. code-block:: console

    $ export AWS_ACCESS_KEY_ID=<AWS_ACCESS_KEY>
    $ export AWS_SECRET_ACCESS_KEY=<AWS_SECRET_ACCESS_KEY>
    $ export AWS_DEFAULT_REGION=<AWS_DEFAULT_REGION>


As an AWS admin, deploy Tibanna Unicorn to Cloud with permissions to a specific S3 bucket.
Name the Unicorn / Unicorn usergroup with the ``--usergroup`` option.
Unicorn is a serverless scheduler, and keeping unicorn on the cloud does not incur extra cost. 
One may have many different unicorns with different names and different bucket permissions.
Then, add other (IAM) users to the user group that has permission to use this unicorn / buckets.

.. code-block:: console

    $ tibanna deploy_unicorn -g <name> -b <bucket>
    $ tibanna add_user -u <username> -g <name>


As a user that has been added to the group (or as an admin), set up the default unicorn.

.. code-block:: console

    $ export TIBANNA_DEFAULT_STEP_FUNCTION_NAME=tibanna_unicorn_<name>


Then, you can run as many snakemake runs as you wish as below, inside a directory that contains
Snakefile and other necessary components (e.g. ``env.yml``, ``config.json``, ...).

.. code-block:: console

    $ snakemake --tibanna --default-remote-prefix=<bucketname>/<subdir> [<other options>]


In this mode, Snakemake will assume all input and output files to be stored in the specified remote location
(a subdirectory inside a given S3 bucket.)
After successful execution, you find your results in the specified remote storage.
Of course, if any input or output already defines a different remote location, the latter will be used instead.
In that case, Tibanna Unicorn must be deployed with all the relevant buckets (``-b bucket1,bucket2,bucket3,...``)
to allow access to the Unicorn serverless components.
Snakemake will assign 3x of the total input size as the allocated space for each execution. The execution may fail
if the total input + output + temp file sizes exceed this estimate.

In addition to regular snakemake options, ``--precommand=<command>`` option allows sending a command to execute before
executing each step on an isolated environment. This kind of command could involve downloading or installing
necessary files that cannot be handled using conda (e.g. the command may begin with ``wget``, ``git clone``, etc.) 


To check Tibanna execution logs, first use ``tibanna stat`` to see the list of all the individual runs.

.. code-block:: console

    $ tibanna stat -n <number_of_executions_to_view> -l


Then, check the detailed log for each job using the Tibanna job id that can be obtained from the first column
of the output of ``tibanna stat``.


.. code-block:: console

    $ tibanna log -j <jobid>


.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.


When executing, Snakemake will make use of the defined resources and threads
to schedule jobs to the correct nodes. In particular, it will forward memory requirements
defined as `mem_mb` to Tibanna. Further, it will propagate the number of threads
a job intends to use, such that Tibanna can allocate it to the most cost-effective
cloud compute instance available.

-----------------------------------------------------------------
Executing a Snakemake workflow via GA4GH TES
-----------------------------------------------------------------

The task execution service (`TES <https://github.com/ga4gh/task-execution-schemas>`_) is an application programming interface developed by the Global Alliance for Genomics and Health (`GA4GH <https://www.ga4gh.org/>`_).
It is used to process workflow tasks in a cloud environment.
A TES server can be easily implemented in a public cloud or at a commercial cloud provider.
Here, the TES standard provides an additional abstraction layer between the execution of a workflow (e.g. on your local machine) and technologies for execution of single tasks (e.g. based Kubernetes or HPC).
We recommend using either `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ or `TESK <https://github.com/EMBL-EBI-TSI/TESK/>`_  to install a TES server.
The guide here is based on Funnel (0.10.0).
To install and configure Funnel follow its official `documentation <https://ohsu-comp-bio.github.io/funnel/docs/>`_.

Configuration
~~~~~~~~~~~~~

Three steps are required to make a Snakemake workflow TES ready:

**Attach conda to rules:**
Execution of Snakemake tasks via TES means, Snakemake is running in a container in the cloud and it executes a specific rule (or a group of rules) with defined input/output data.
By default, the TES module uses the latest Snakemake container.
Running Snakemake within a container requires having all external tools installed within this container.
This can be done by providing a custom container image having installed Snakemake and other all required tools (e.g. BWA).
Or it can be done by attaching a conda environment to each rule, such that those tools will be installed within the running container.
For simplicity, this guide recommends to attach a specific conda environment to each rule, although it is more efficient in the long term to provide custom container images.

**Use remote files:**
The TES module requires using a remote file storage system for input/output files such that all files are available on the cloud machines and within their running container.
There are several options available in Snakemake to use remote files.
This guide recommends to use S3 (or SWIFT) object storage.

**Install py-tes module:**
TES backend requires py-tes to be installed. Please install py-tes, e.g. via Conda or Pip.

.. code-block:: console

    $ pip install py-tes 

Execution
~~~~~~~~~

Funnel starts container in read only mode, which is good practice.
Anyhow, using the default Snakemake container image will likely require installing additional software within the running container.
Therefore, we need to set two conda specific variables such that new environments will be installed at `/tmp` which will be mounted as a writable volume in the container.

.. code-block:: console

    $ export CONDA_PKGS_DIRS=/tmp/conda
    $ export CONDA_ENVS_PATH=/tmp/conda

Next, using S3 or SWIFT storage, we also need to set credentials.

.. code-block:: console

    $ export AWS_ACCESS_KEY_ID=YOUR_ACCESS_KEY
    $ export AWS_SECRET_ACCESS_KEY=YOUR_SECRET_ACCESS_KEY

Now we can run Snakemake using:

.. code-block:: console

    $ snakemake \
        --tes $TES_URL \
        --use-conda \
        --envvars CONDA_PKGS_DIRS CONDA_ENVS_PATH AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY \
        --conda-prefix $CONDA_ENVS_PATH \
        all
.. _cluster:

=================
Cluster Execution
=================


Snakemake can make use of cluster engines that support shell scripts and have access to a common filesystem, (e.g. the Sun Grid Engine).
In this case, Snakemake simply needs to be given a submit command that accepts a shell script as first positional argument:

.. code-block:: console

    $ snakemake --cluster qsub --jobs 32


Here, ``--jobs`` denotes the number of jobs submitted to the cluster at the same time (here 32).
The cluster command can be decorated with job specific information, e.g.

.. sidebar:: Note

  Consider to :ref:`group jobs <snakefiles-grouping>` in order to minimize overhead, in particular for short-running jobs.


.. code-block:: console

    $ snakemake --cluster "qsub {threads}"

Thereby, all keywords of a rule are allowed (e.g. rulename, params, input, output, threads, priority, resources, ...).
For example, you could encode the expected running time in minutes into a :ref:`resource <snakefiles-resources>` ``runtime_min``:

.. code-block:: python

    rule:
        input:  
            ...
        output:
            ...
        resources: 
            runtime_min=240
        shell:
            ...

and forward it to the cluster scheduler:

.. code-block:: console

    $ snakemake --cluster "qsub --runtime {resources.runtime}"

In order to avoid specifying ``runtime_min`` for each rule, you can make use of the ``--default-resources`` flag, see ``snakemake --help``.

If your cluster system supports `DRMAA <https://www.drmaa.org/>`_, Snakemake can make use of that to increase the control over jobs.
E.g. jobs can be cancelled upon pressing ``Ctrl+C``, which is not possible with the generic ``--cluster`` support.
With DRMAA, no ``qsub`` command needs to be provided, but system specific arguments can still be given as a string, e.g.

.. code-block:: console

    $ snakemake --drmaa " -q username" -j 32

Note that the string has to contain a leading whitespace.
Else, the arguments will be interpreted as part of the normal Snakemake arguments, and execution will fail.

Adapting to a specific cluster can involve quite a lot of options. It is therefore a good idea to setup a :ref:`a profile <profiles>`.

--------------
Job Properties
--------------

When executing a workflow on a cluster using the ``--cluster`` parameter (see below), Snakemake creates a job script for each job to execute. This script is then invoked using the provided cluster submission command (e.g. ``qsub``). Sometimes you want to provide a custom wrapper for the cluster submission command that decides about additional parameters. As this might be based on properties of the job, Snakemake stores the job properties (e.g. name, rulename, threads, input, output, params etc.) as JSON inside the job script (for group jobs, the rulename will be "GROUP", otherwise it will be the same as the job name). For convenience, there exists a parser function `snakemake.utils.read_job_properties` that can be used to access the properties. The following shows an example job submission wrapper:

.. code-block:: python

    #!python

    #!/usr/bin/env python3
    import os
    import sys

    from snakemake.utils import read_job_properties

    jobscript = sys.argv[1]
    job_properties = read_job_properties(jobscript)

    # do something useful with the threads
    threads = job_properties[threads]

    # access property defined in the cluster configuration file (Snakemake >=3.6.0)
    job_properties["cluster"]["time"]

    os.system("qsub -t {threads} {script}".format(threads=threads, script=jobscript))
