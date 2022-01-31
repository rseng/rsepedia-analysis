# 0.15.4 (Thu Dec 16 2021)

#### 🐛 Bug Fix

- BF: autorc - replace incorrect releaseTypes with "none" [#6320](https://github.com/datalad/datalad/pull/6320) ([@yarikoptic](https://github.com/yarikoptic))
- Minor enhancement to CONTRIBUTING.md [#6309](https://github.com/datalad/datalad/pull/6309) ([@bpoldrack](https://github.com/bpoldrack))
- UX: If a clean repo is dirty after a failed run, give clean-up hints [#6112](https://github.com/datalad/datalad/pull/6112) ([@adswa](https://github.com/adswa))
- Stop using distutils [#6113](https://github.com/datalad/datalad/pull/6113) ([@jwodder](https://github.com/jwodder))
- BF: RIARemote - set UI backend to annex to make it interactive [#6287](https://github.com/datalad/datalad/pull/6287) ([@yarikoptic](https://github.com/yarikoptic) [@bpoldrack](https://github.com/bpoldrack))
- Fix invalid escape sequences [#6293](https://github.com/datalad/datalad/pull/6293) ([@jwodder](https://github.com/jwodder))
- CI: Update environment for windows CI builds [#6292](https://github.com/datalad/datalad/pull/6292) ([@bpoldrack](https://github.com/bpoldrack))
- bump the python version used for mac os tests [#6288](https://github.com/datalad/datalad/pull/6288) ([@christian-monch](https://github.com/christian-monch) [@bpoldrack](https://github.com/bpoldrack))
- ENH(UX): log a hint to use ulimit command in case of "Too long" exception [#6173](https://github.com/datalad/datalad/pull/6173) ([@yarikoptic](https://github.com/yarikoptic))
- Report correct HTTP URL for RIA store content [#6091](https://github.com/datalad/datalad/pull/6091) ([@mih](https://github.com/mih))
- BF: Don't overwrite subdataset source candidates [#6168](https://github.com/datalad/datalad/pull/6168) ([@bpoldrack](https://github.com/bpoldrack))
- Bump sphinx requirement to bypass readthedocs defaults [#6189](https://github.com/datalad/datalad/pull/6189) ([@mih](https://github.com/mih))
- infra: Provide custom prefix to auto-related labels [#6151](https://github.com/datalad/datalad/pull/6151) ([@adswa](https://github.com/adswa))
- Remove all usage of exc_str() [#6142](https://github.com/datalad/datalad/pull/6142) ([@mih](https://github.com/mih))
- BF: obtain information about annex special remotes also from annex journal [#6135](https://github.com/datalad/datalad/pull/6135) ([@yarikoptic](https://github.com/yarikoptic) [@mih](https://github.com/mih))
- BF: clone tried to save new subdataset despite failing to clone [#6140](https://github.com/datalad/datalad/pull/6140) ([@bpoldrack](https://github.com/bpoldrack))

#### 🧪 Tests

- RF+BF: use skip_if_no_module helper instead of try/except for libxmp and boto [#6148](https://github.com/datalad/datalad/pull/6148) ([@yarikoptic](https://github.com/yarikoptic))
- git://github.com -> https://github.com [#6134](https://github.com/datalad/datalad/pull/6134) ([@mih](https://github.com/mih))

#### Authors: 6

- Adina Wagner ([@adswa](https://github.com/adswa))
- Benjamin Poldrack ([@bpoldrack](https://github.com/bpoldrack))
- Christian Mönch ([@christian-monch](https://github.com/christian-monch))
- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Michael Hanke ([@mih](https://github.com/mih))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.15.3 (Sat Oct 30 2021)

#### 🐛 Bug Fix

- BF: Don't make create-sibling recursive by default [#6116](https://github.com/datalad/datalad/pull/6116) ([@adswa](https://github.com/adswa))
- BF: Add dashes to 'force' option in non-empty directory error message [#6078](https://github.com/datalad/datalad/pull/6078) ([@DisasterMo](https://github.com/DisasterMo))
- DOC: Add supported URL types to download-url's docstring [#6098](https://github.com/datalad/datalad/pull/6098) ([@adswa](https://github.com/adswa))
- BF: Retain git-annex error messages & don't show them if operation successful [#6070](https://github.com/datalad/datalad/pull/6070) ([@DisasterMo](https://github.com/DisasterMo))
- Remove uses of `__full_version__` and `datalad.version` [#6073](https://github.com/datalad/datalad/pull/6073) ([@jwodder](https://github.com/jwodder))
- BF: ORA shouldn't crash while handling a failure [#6063](https://github.com/datalad/datalad/pull/6063) ([@bpoldrack](https://github.com/bpoldrack))
- DOC: Refine --reckless docstring on usage and wording [#6043](https://github.com/datalad/datalad/pull/6043) ([@adswa](https://github.com/adswa))
- BF: archives upon strip - use rmtree which retries etc instead of rmdir [#6064](https://github.com/datalad/datalad/pull/6064) ([@yarikoptic](https://github.com/yarikoptic))
- BF: do not leave test in a tmp dir destined for removal [#6059](https://github.com/datalad/datalad/pull/6059) ([@yarikoptic](https://github.com/yarikoptic))
- Next wave of exc_str() removals [#6022](https://github.com/datalad/datalad/pull/6022) ([@mih](https://github.com/mih))

#### ⚠️ Pushed to `maint`

- CI: Enable new codecov uploader in Appveyor CI ([@adswa](https://github.com/adswa))

#### 🏠 Internal

- UX: Log clone-candidate number and URLs [#6092](https://github.com/datalad/datalad/pull/6092) ([@adswa](https://github.com/adswa))
- UX/ENH: Disable reporting, and don't do superfluous internal subdatasets calls [#6094](https://github.com/datalad/datalad/pull/6094) ([@adswa](https://github.com/adswa))
- Update codecov action to v2 [#6072](https://github.com/datalad/datalad/pull/6072) ([@jwodder](https://github.com/jwodder))

#### 📝 Documentation

- Design document on URL substitution feature [#6065](https://github.com/datalad/datalad/pull/6065) ([@mih](https://github.com/mih))

#### 🧪 Tests

- BF(TST): remove reuse of the same tape across unrelated tests [#6127](https://github.com/datalad/datalad/pull/6127) ([@yarikoptic](https://github.com/yarikoptic))
- Fail Travis tests on deprecation warnings [#6074](https://github.com/datalad/datalad/pull/6074) ([@jwodder](https://github.com/jwodder))
- Ux get result handling broken [#6052](https://github.com/datalad/datalad/pull/6052) ([@christian-monch](https://github.com/christian-monch))
- enable metalad tests again [#6060](https://github.com/datalad/datalad/pull/6060) ([@christian-monch](https://github.com/christian-monch))

#### Authors: 7

- Adina Wagner ([@adswa](https://github.com/adswa))
- Benjamin Poldrack ([@bpoldrack](https://github.com/bpoldrack))
- Christian Mönch ([@christian-monch](https://github.com/christian-monch))
- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Michael Burgardt ([@DisasterMo](https://github.com/DisasterMo))
- Michael Hanke ([@mih](https://github.com/mih))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.15.2 (Wed Oct 06 2021)

#### 🐛 Bug Fix

- BF: Don't suppress datalad subdatasets output [#6035](https://github.com/datalad/datalad/pull/6035) ([@DisasterMo](https://github.com/DisasterMo) [@mih](https://github.com/mih))
- Honor datalad.runtime.use-patool if set regardless of OS (was Windows only) [#6033](https://github.com/datalad/datalad/pull/6033) ([@mih](https://github.com/mih))
- Discontinue usage of deprecated (public) helper [#6032](https://github.com/datalad/datalad/pull/6032) ([@mih](https://github.com/mih))
- BF: ProgressHandler - close the other handler if was specified [#6020](https://github.com/datalad/datalad/pull/6020) ([@yarikoptic](https://github.com/yarikoptic))
- UX: Report GitLab weburl of freshly created projects in the result [#6017](https://github.com/datalad/datalad/pull/6017) ([@adswa](https://github.com/adswa))
- Ensure there's a blank line between the class `__doc__` and "Parameters" in `build_doc` docstrings [#6004](https://github.com/datalad/datalad/pull/6004) ([@jwodder](https://github.com/jwodder))
- Large code-reorganization of everything runner-related [#6008](https://github.com/datalad/datalad/pull/6008) ([@mih](https://github.com/mih))
- Discontinue exc_str() in all modern parts of the code base [#6007](https://github.com/datalad/datalad/pull/6007) ([@mih](https://github.com/mih))

#### 🧪 Tests

- TST: Add test to ensure functionality with subdatasets starting with a hyphen (-) [#6042](https://github.com/datalad/datalad/pull/6042) ([@DisasterMo](https://github.com/DisasterMo))
- BF(TST): filter away warning from coverage from analysis of stderr of --help [#6028](https://github.com/datalad/datalad/pull/6028) ([@yarikoptic](https://github.com/yarikoptic))
- BF: disable outdated SSL root certificate breaking chain on older/buggy clients [#6027](https://github.com/datalad/datalad/pull/6027) ([@yarikoptic](https://github.com/yarikoptic))
- BF: start global test_http_server only if not running already [#6023](https://github.com/datalad/datalad/pull/6023) ([@yarikoptic](https://github.com/yarikoptic))

#### Authors: 5

- Adina Wagner ([@adswa](https://github.com/adswa))
- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Michael Burgardt ([@DisasterMo](https://github.com/DisasterMo))
- Michael Hanke ([@mih](https://github.com/mih))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.15.1 (Fri Sep 24 2021)

#### 🐛 Bug Fix

- BF: downloader - fail to download even on non-crippled FS if symlink exists [#5991](https://github.com/datalad/datalad/pull/5991) ([@yarikoptic](https://github.com/yarikoptic))
- ENH: import datalad.api to bind extensions methods for discovery of dataset methods [#5999](https://github.com/datalad/datalad/pull/5999) ([@yarikoptic](https://github.com/yarikoptic))
- Restructure cmdline API presentation [#5988](https://github.com/datalad/datalad/pull/5988) ([@mih](https://github.com/mih))
- Close file descriptors after process exit [#5983](https://github.com/datalad/datalad/pull/5983) ([@mih](https://github.com/mih))

#### ⚠️ Pushed to `maint`

- Discontinue testing of hirni extension ([@mih](https://github.com/mih))

#### 🏠 Internal

- Add debugging information to release step [#5980](https://github.com/datalad/datalad/pull/5980) ([@jwodder](https://github.com/jwodder))

#### 📝 Documentation

- Coarse description of the credential subsystem's functionality [#5998](https://github.com/datalad/datalad/pull/5998) ([@mih](https://github.com/mih))

#### 🧪 Tests

- BF(TST): use sys.executable, mark test_ria_basics.test_url_keys as requiring network [#5986](https://github.com/datalad/datalad/pull/5986) ([@yarikoptic](https://github.com/yarikoptic))

#### Authors: 3

- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Michael Hanke ([@mih](https://github.com/mih))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.15.0 (Tue Sep 14 2021) --  We miss you Kyle!

#### Enhancements and new features

- Command execution is now performed by a new `Runner` implementation that is
  no longer based on the `asyncio` framework, which was found to exhibit
  fragile performance in interaction with other `asyncio`-using code, such as
  Jupyter notebooks. The new implementation is based on threads. It also supports
  the specification of "protocols" that were introduced with the switch to the
  `asyncio` implementation in 0.14.0. ([#5667][]) 

- `clone` now supports arbitrary URL transformations based on regular
  expressions. One or more transformation steps can be defined via
  `datalad.clone.url-substitute.<label>` configuration settings. The feature can
  be (and is now) used to support convenience mappings, such as
  `https://osf.io/q8xnk/` (displayed in a browser window) to `osf://q8xnk`
  (clonable via the `datalad-osf` extension. ([#5749][])

- Homogenize SSH use and configurability between DataLad and git-annex, by
  instructing git-annex to use DataLad's `sshrun` for SSH calls (instead of SSH
  directly). ([#5389][])

- The ORA special remote has received several new features:

  - It now support a `push-url` setting as an alternative to `url` for write
    access. An analog parameter was also added to `create-sibling-ria`.
    ([#5420][], [#5428][])

  - Access of RIA stores now performs homogeneous availability checks,
    regardless of access protocol. Before, broken HTTP-based access due to
    misspecified URLs could have gone unnoticed. ([#5459][], [#5672][])

  - Error reporting was introduce to inform about undesirable conditions in
    remote RIA stores. ([#5683][])

- `create-sibling-ria` now supports `--alias` for the specification of a
  convenience dataset alias name in a RIA store. ([#5592][])

- Analog to `git commit`, `save` now features an `--amend` mode to support
  incremental updates of a dataset state. ([#5430][])

- `run` now supports a dry-run mode that can be used to inspect the result of
  parameter expansion on the effective command to ease the composition of more
  complicated command lines. ([#5539][])

- `run` now supports a `--assume-ready` switch to avoid the (possibly
  expensive) preparation of inputs and outputs with large datasets that have
  already been readied through other means. ([#5431][])

- `update` now features `--how` and `--how-subds` parameters to configure how
  an update shall be performed. Supported modes are `fetch` (unchanged
  default), and `merge` (previously also possible via `--merge`), but also new
  strategies like `reset` or `checkout`. ([#5534][])

- `update` has a new `--follow=parentds-lazy` mode that only performs a fetch
  operation in subdatasets when the desired commit is not yet present. During
  recursive updates involving many subdatasets this can substantially speed up
  performance. ([#5474][])

- DataLad's command line API can now report the version for individual commands
  via `datalad <cmd> --version`. The output has been homogenized to
  `<providing package> <version>`. ([#5543][])

- `create-sibling` now logs information on an auto-generated sibling name, in
  the case that no `--name/-s` was provided. ([#5550][])

- `create-sibling-github` has been updated to emit result records like any
  standard DataLad command. Previously it was implemented as a "plugin", which
  did not support all standard API parameters. ([#5551][])

- `copy-file` now also works with content-less files in datasets on crippled
  filesystems (adjusted mode), when a recent enough git-annex (8.20210428 or
  later) is available. ([#5630][])

- `addurls` can now be instructed how to behave in the event of file name
  collision via a new parameter `--on-collision`. ([#5675][])

- `addurls` reporting now informs which particular subdatasets were created.
  ([#5689][])

- Credentials can now be provided or overwritten via all means supported by
  `ConfigManager`. Importantly, `datalad.credential.<name>.<field>`
  configuration settings and analog specification via environment variables are
  now supported (rather than custom environment variables only). Previous
  specification methods are still supported too. ([#5680][])

- A new `datalad.credentials.force-ask` configuration flag can now be used to
  force re-entry of already known credentials. This simplifies credential
  updates without having to use an approach native to individual credential
  stores. ([#5777][])

- Suppression of rendering repeated similar results is now configurable via the
  configuration switches `datalad.ui.suppress-similar-results` (bool), and
  `datalad.ui.suppress-similar-results-threshold` (int). ([#5681][])

- The performance of `status` and similar functionality when determining local
  file availability has been improved. ([#5692][])

- `push` now renders a result summary on completion. ([#5696][])

- A dedicated info log message indicates when dataset repositories are
  subjected to an annex version upgrade. ([#5698][])

- Error reporting improvements:

  - The `NoDatasetFound` exception now provides information for which purpose a
    dataset is required. ([#5708][])

  - Wording of the `MissingExternalDependeny` error was rephrased to account
    for cases of non-functional installations. ([#5803][])

  - `push` reports when a `--to` parameter specification was (likely)
    forgotten. ([#5726][])

  - Detailed information is now given when DataLad fails to obtain a lock for
    credential entry in a timely fashion. Previously only a generic debug log
    message was emitted. ([#5884][])

  - Clarified error message when `create-sibling-gitlab` was called without
    `--project`. ([#5907][])

- `add-readme` now provides a README template with more information on the
  nature and use of DataLad datasets. A README file is no longer annex'ed by
  default, but can be using the new `--annex` switch. ([#5723][], [#5725][])

- `clean` now supports a `--dry-run` mode to inform about cleanable content.
  ([#5738][])

- A new configuration setting `datalad.locations.locks` can be used to control
  the placement of lock files. ([#5740][])

- `wtf` now also reports branch names and states. ([#5804][])

- `AnnexRepo.whereis()` now supports batch mode. ([#5533][])

### Deprecations and removals

- The minimum supported git-annex version is now 8.20200309. ([#5512][])

- ORA special remote configuration items `ssh-host`, and `base-path` are
  deprecated. They are completely replaced by `ria+<protocol>://` URL
  specifications. ([#5425][])

- The deprecated `no_annex` parameter of `create()` was removed from the Python
  API. ([#5441][])

- The unused `GitRepo.pull()` method has been removed. ([#5558][])

- Residual support for "plugins" (a mechanism used before DataLad supported
  extensions) was removed. This includes the configuration switches
  `datalad.locations.{system,user}-plugins`. ([#5554][], [#5564][])

- Several features and comments have been moved to the `datalad-deprecated`
  package. This package must now be installed to be able to use keep using this
  functionality.

  - The `publish` command. Use `push` instead. ([#5837][])

  - The `ls` command. ([#5569][])

  - The web UI that is deployable via `datalad create-sibling --ui`. ([#5555][])

  - The "automagic IO" feature. ([#5577][])

- `AnnexRepo.copy_to()` has been deprecated. The `push` command should be used
  instead. ([#5560][])

- `AnnexRepo.sync()` has been deprecated. `AnnexRepo.call_annex(['sync', ...])`
  should be used instead. ([#5461][])

- All `GitRepo.*_submodule()` methods have been deprecated and will be removed
  in a future release. ([#5559][])

- `create-sibling-github`'s `--dryrun` switch was deprecated, use `--dry-run` instead.
  ([#5551][])

- The `datalad --pbs-runner` option has been deprecated, use `condor_run`
  (or similar) instead. ([#5956][])

#### 🐛 Fixes

- Prevent invalid declaration of a publication dependencies for 'origin' on any
  auto-detected ORA special remotes, when cloing from a RIA store. An ORA
  remote is now checked whether it actually points to the RIA store the clone was
  made from. ([#5415][])

- The ORA special remote implementation has received several fixes:

  - It can now handle HTTP redirects. ([#5792][])

  - Prevents failure when URL-type annex keys contain the '/' character.
    ([#5823][])

  - Properly support the specification of usernames, passwords and ports in
    `ria+<protocol>://` URLs. ([#5902][])

- It is now possible to specifically select the default (or generic) result
  renderer via `datalad -f default` and with that override a `tailored` result
  renderer that may be preconfigured for a particular command. ([#5476][])

- Starting with 0.14.0, original URLs given to `clone` were recorded in a
  subdataset record. This was initially done in a second commit, leading to
  inflation of commits and slowdown in superdatasets with many subdatasets. Such
  subdataset record annotation is now collapsed into a single commits.
  ([#5480][]) 

- `run` now longer removes leading empty directories as part of the output
  preparation. This was surprising behavior for commands that do not ensure on
  their own that output directories exist. ([#5492][])

- A potentially existing `message` property is no longer removed when using the
  `json` or `json_pp` result renderer to avoid undesired withholding of
  relevant information. ([#5536][])

- `subdatasets` now reports `state=present`, rather than `state=clean`, for
  installed subdatasets to complement `state=absent` reports for uninstalled
  dataset. ([#5655][])

- `create-sibling-ria` now executes commands with a consistent environment
  setup that matches all other command execution in other DataLad commands.
  ([#5682][])

- `save` no longer saves unspecified subdatasets when called with an explicit
  path (list). The fix required a behavior change of
  `GitRepo.get_content_info()` in its interpretation of `None` vs. `[]` path
  argument values that now aligns the behavior of `GitRepo.diff|status()` with
  their respective documentation. ([#5693][])

- `get` now prefers the location of a subdatasets that is recorded in a
  superdataset's `.gitmodules` record. Previously, DataLad tried to obtain a
  subdataset from an assumed checkout of the superdataset's origin. This new
  default order is (re-)configurable via the
  `datalad.get.subdataset-source-candidate-<priority-label>` configuration
  mechanism. ([#5760][])

- `create-sibling-gitlab` no longer skips the root dataset when `.` is given as
  a path. ([#5789][])

- `siblings` now rejects a value given to `--as-common-datasrc` that clashes
  with the respective Git remote. ([#5805][])

- The usage synopsis reported by `siblings` now lists all supported actions.
  ([#5913][])

- `siblings` now renders non-ok results to avoid silent failure. ([#5915][])

- `.gitattribute` file manipulations no longer leave the file without a
  trailing newline. ([#5847][])

- Prevent crash when trying to delete a non-existing keyring credential field.
  ([#5892][])

- git-annex is no longer called with an unconditional `annex.retry=3`
  configuration. Instead, this parameterization is now limited to `annex get`
  and `annex copy` calls. ([#5904][])

#### 🧪 Tests

- `file://` URLs are no longer the predominant test case for `AnnexRepo`
  functionality. A built-in HTTP server now used in most cases. ([#5332][])

---

# 0.14.8 (Sun Sep 12 2021)

#### 🐛 Bug Fix

- BF: add-archive-content on .xz and other non-.gz stream compressed files [#5930](https://github.com/datalad/datalad/pull/5930) ([@yarikoptic](https://github.com/yarikoptic))
- BF(UX): do not keep logging ERROR possibly present in progress records [#5936](https://github.com/datalad/datalad/pull/5936) ([@yarikoptic](https://github.com/yarikoptic))
- Annotate datalad_core as not needing actual data -- just uses annex whereis [#5971](https://github.com/datalad/datalad/pull/5971) ([@yarikoptic](https://github.com/yarikoptic))
- BF: limit CMD_MAX_ARG if obnoxious value is encountered. [#5945](https://github.com/datalad/datalad/pull/5945) ([@yarikoptic](https://github.com/yarikoptic))
- Download session/credentials locking -- inform user if locking is "failing" to be obtained, fail upon ~5min timeout [#5884](https://github.com/datalad/datalad/pull/5884) ([@yarikoptic](https://github.com/yarikoptic))
- Render siblings()'s non-ok results with the default renderer [#5915](https://github.com/datalad/datalad/pull/5915) ([@mih](https://github.com/mih))
- BF: do not crash, just skip whenever trying to delete non existing field in the underlying keyring [#5892](https://github.com/datalad/datalad/pull/5892) ([@yarikoptic](https://github.com/yarikoptic))
- Fix argument-spec for `siblings` and improve usage synopsis [#5913](https://github.com/datalad/datalad/pull/5913) ([@mih](https://github.com/mih))
- Clarify error message re unspecified gitlab project [#5907](https://github.com/datalad/datalad/pull/5907) ([@mih](https://github.com/mih))
- Support username, password and port specification in RIA URLs [#5902](https://github.com/datalad/datalad/pull/5902) ([@mih](https://github.com/mih))
- BF: take path from SSHRI, test URLs not only on Windows [#5881](https://github.com/datalad/datalad/pull/5881) ([@yarikoptic](https://github.com/yarikoptic))
- ENH(UX): warn user if keyring returned a "null" keyring [#5875](https://github.com/datalad/datalad/pull/5875) ([@yarikoptic](https://github.com/yarikoptic))
- ENH(UX): state original purpose in NoDatasetFound exception + detail it for get [#5708](https://github.com/datalad/datalad/pull/5708) ([@yarikoptic](https://github.com/yarikoptic))

#### ⚠️ Pushed to `maint`

- Merge branch 'bf-http-headers-agent' into maint ([@yarikoptic](https://github.com/yarikoptic))
- RF(BF?)+DOC: provide User-Agent to entire session headers + use those if provided ([@yarikoptic](https://github.com/yarikoptic))

#### 🏠 Internal

- Pass `--no-changelog` to `auto shipit` if changelog already has entry [#5952](https://github.com/datalad/datalad/pull/5952) ([@jwodder](https://github.com/jwodder))
- Add isort config to match current convention + run isort via pre-commit (if configured) [#5923](https://github.com/datalad/datalad/pull/5923) ([@jwodder](https://github.com/jwodder))
- .travis.yml: use python -m {nose,coverage} invocations, and always show combined report [#5888](https://github.com/datalad/datalad/pull/5888) ([@yarikoptic](https://github.com/yarikoptic))
- Add project URLs into the package metadata for convenience links on Pypi [#5866](https://github.com/datalad/datalad/pull/5866) ([@adswa](https://github.com/adswa) [@yarikoptic](https://github.com/yarikoptic))

#### 🧪 Tests

- BF: do use OBSCURE_FILENAME instead of hardcoded unicode [#5944](https://github.com/datalad/datalad/pull/5944) ([@yarikoptic](https://github.com/yarikoptic))
- BF(TST): Skip testing for having PID listed if no psutil [#5920](https://github.com/datalad/datalad/pull/5920) ([@yarikoptic](https://github.com/yarikoptic))
- BF(TST): Boost version of git-annex to 8.20201129 to test an error message [#5894](https://github.com/datalad/datalad/pull/5894) ([@yarikoptic](https://github.com/yarikoptic))

#### Authors: 4

- Adina Wagner ([@adswa](https://github.com/adswa))
- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Michael Hanke ([@mih](https://github.com/mih))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.14.7 (Tue Aug 03 2021)

#### 🐛 Bug Fix

- UX: When two or more clone URL templates are found, error out more gracefully [#5839](https://github.com/datalad/datalad/pull/5839) ([@adswa](https://github.com/adswa))
- BF: http_auth - follow redirect (just 1) to re-authenticate after initial attempt [#5852](https://github.com/datalad/datalad/pull/5852) ([@yarikoptic](https://github.com/yarikoptic))
- addurls Formatter - provide value repr in exception [#5850](https://github.com/datalad/datalad/pull/5850) ([@yarikoptic](https://github.com/yarikoptic))
- ENH: allow for "patch" level semver for "master" branch [#5839](https://github.com/datalad/datalad/pull/5839) ([@yarikoptic](https://github.com/yarikoptic))
- BF: Report info from annex JSON error message in CommandError [#5809](https://github.com/datalad/datalad/pull/5809) ([@mih](https://github.com/mih))
- RF(TST): do not test for no EASY and pkg_resources in shims [#5817](https://github.com/datalad/datalad/pull/5817) ([@yarikoptic](https://github.com/yarikoptic))
- http downloaders: Provide custom informative User-Agent, do not claim to be "Authenticated access" [#5802](https://github.com/datalad/datalad/pull/5802) ([@yarikoptic](https://github.com/yarikoptic))
- ENH(UX,DX): inform user with a warning if version is 0+unknown [#5787](https://github.com/datalad/datalad/pull/5787) ([@yarikoptic](https://github.com/yarikoptic))
- shell-completion: add argcomplete to 'misc' extra_depends, log an ERROR if argcomplete fails to import [#5781](https://github.com/datalad/datalad/pull/5781) ([@yarikoptic](https://github.com/yarikoptic))
- ENH (UX): add python-gitlab dependency [#5776](https://github.com/datalad/datalad/pull/5776) (s.heunis@fz-juelich.de)

#### 🏠 Internal

- BF: Fix reported paths in ORA remote [#5821](https://github.com/datalad/datalad/pull/5821) ([@adswa](https://github.com/adswa))
- BF: import importlib.metadata not importlib_metadata whenever available [#5818](https://github.com/datalad/datalad/pull/5818) ([@yarikoptic](https://github.com/yarikoptic))

#### 🧪 Tests

- TST: set --allow-unrelated-histories in the mk_push_target setup for Windows [#5855](https://github.com/datalad/datalad/pull/5855) ([@adswa](https://github.com/adswa))
- Tests: Allow for version to contain + as a separator and provide more information for version related comparisons [#5786](https://github.com/datalad/datalad/pull/5786) ([@yarikoptic](https://github.com/yarikoptic))

#### Authors: 4

- Adina Wagner ([@adswa](https://github.com/adswa))
- Michael Hanke ([@mih](https://github.com/mih))
- Stephan Heunis ([@jsheunis](https://github.com/jsheunis))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.14.6 (Sun Jun 27 2021)

#### 🏠 Internal

- BF: update changelog conversion from .md to .rst (for sphinx) [#5757](https://github.com/datalad/datalad/pull/5757) ([@yarikoptic](https://github.com/yarikoptic) [@jwodder](https://github.com/jwodder))

#### Authors: 2

- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.14.5 (Mon Jun 21 2021)

#### 🐛 Bug Fix

- BF(TST): parallel - take longer for producer to produce [#5747](https://github.com/datalad/datalad/pull/5747) ([@yarikoptic](https://github.com/yarikoptic))
- add --on-failure default value and document it [#5690](https://github.com/datalad/datalad/pull/5690) ([@christian-monch](https://github.com/christian-monch) [@yarikoptic](https://github.com/yarikoptic))
- ENH: harmonize "purpose" statements to imperative form [#5733](https://github.com/datalad/datalad/pull/5733) ([@yarikoptic](https://github.com/yarikoptic))
- ENH(TST): populate heavy tree with 100 unique keys (not just 1) among 10,000 [#5734](https://github.com/datalad/datalad/pull/5734) ([@yarikoptic](https://github.com/yarikoptic))
- BF: do not use .acquired - just get state from acquire() [#5718](https://github.com/datalad/datalad/pull/5718) ([@yarikoptic](https://github.com/yarikoptic))
- BF: account for annex now "scanning for annexed" instead of "unlocked" files [#5705](https://github.com/datalad/datalad/pull/5705) ([@yarikoptic](https://github.com/yarikoptic))
- interface: Don't repeat custom summary for non-generator results [#5688](https://github.com/datalad/datalad/pull/5688) ([@kyleam](https://github.com/kyleam))
- RF: just pip install datalad-installer [#5676](https://github.com/datalad/datalad/pull/5676) ([@yarikoptic](https://github.com/yarikoptic))
- DOC: addurls.extract: Drop mention of removed 'stream' parameter [#5690](https://github.com/datalad/datalad/pull/5690) ([@kyleam](https://github.com/kyleam))
- Merge pull request #5674 from kyleam/test-addurls-copy-fix [#5674](https://github.com/datalad/datalad/pull/5674) ([@kyleam](https://github.com/kyleam))
- Merge pull request #5663 from kyleam/status-ds-equal-path [#5663](https://github.com/datalad/datalad/pull/5663) ([@kyleam](https://github.com/kyleam))
- Merge pull request #5671 from kyleam/update-fetch-fail [#5671](https://github.com/datalad/datalad/pull/5671) ([@kyleam](https://github.com/kyleam))
- BF: update: Honor --on-failure if fetch fails [#5671](https://github.com/datalad/datalad/pull/5671) ([@kyleam](https://github.com/kyleam))
- RF: update: Avoid fetch's deprecated kwargs [#5671](https://github.com/datalad/datalad/pull/5671) ([@kyleam](https://github.com/kyleam))
- CLN: update: Drop an unused import [#5671](https://github.com/datalad/datalad/pull/5671) ([@kyleam](https://github.com/kyleam))
- Merge pull request #5664 from kyleam/addurls-better-url-parts-error [#5664](https://github.com/datalad/datalad/pull/5664) ([@kyleam](https://github.com/kyleam))
- Merge pull request #5661 from kyleam/sphinx-fix-plugin-refs [#5661](https://github.com/datalad/datalad/pull/5661) ([@kyleam](https://github.com/kyleam))
- BF: status: Provide special treatment of "this dataset" path [#5663](https://github.com/datalad/datalad/pull/5663) ([@kyleam](https://github.com/kyleam))
- BF: addurls: Provide better placeholder error for special keys [#5664](https://github.com/datalad/datalad/pull/5664) ([@kyleam](https://github.com/kyleam))
- RF: addurls: Simply construction of placeholder exception message [#5664](https://github.com/datalad/datalad/pull/5664) ([@kyleam](https://github.com/kyleam))
- RF: addurls._get_placeholder_exception: Rename a parameter [#5664](https://github.com/datalad/datalad/pull/5664) ([@kyleam](https://github.com/kyleam))
- RF: status: Avoid repeated Dataset.path access [#5663](https://github.com/datalad/datalad/pull/5663) ([@kyleam](https://github.com/kyleam))
- DOC: Reference plugins via datalad.api [#5661](https://github.com/datalad/datalad/pull/5661) ([@kyleam](https://github.com/kyleam))
- download-url: Set up datalad special remote if needed [#5648](https://github.com/datalad/datalad/pull/5648) ([@kyleam](https://github.com/kyleam) [@yarikoptic](https://github.com/yarikoptic))

#### ⚠️ Pushed to `maint`

- MNT: Post-release dance ([@kyleam](https://github.com/kyleam))

#### 🏠 Internal

- Switch to versioneer and auto [#5669](https://github.com/datalad/datalad/pull/5669) ([@jwodder](https://github.com/jwodder) [@yarikoptic](https://github.com/yarikoptic))
- MNT: setup.py: Temporarily avoid Sphinx 4 [#5649](https://github.com/datalad/datalad/pull/5649) ([@kyleam](https://github.com/kyleam))

#### 🧪 Tests

- BF(TST): skip testing for showing "Scanning for ..." since not shown if too quick [#5727](https://github.com/datalad/datalad/pull/5727) ([@yarikoptic](https://github.com/yarikoptic))
- Revert "TST: test_partial_unlocked: Document and avoid recent git-annex failure" [#5651](https://github.com/datalad/datalad/pull/5651) ([@kyleam](https://github.com/kyleam))

#### Authors: 4

- Christian Mönch ([@christian-monch](https://github.com/christian-monch))
- John T. Wodder II ([@jwodder](https://github.com/jwodder))
- Kyle Meyer ([@kyleam](https://github.com/kyleam))
- Yaroslav Halchenko ([@yarikoptic](https://github.com/yarikoptic))

---

# 0.14.4 (May 10, 2021) -- .

## Fixes

- Following an internal call to `git-clone`, [clone][] assumed that
  the remote name was "origin", but this may not be the case if
  `clone.defaultRemoteName` is configured (available as of Git 2.30).
  ([#5572][])

- Several test fixes, including updates for changes in git-annex.
  ([#5612][]) ([#5632][]) ([#5639][])


# 0.14.3 (April 28, 2021) -- .

## Fixes

- For outputs that include a glob, [run][] didn't re-glob after
  executing the command, which is necessary to catch changes if
  `--explicit` or `--expand={outputs,both}` is specified.  ([#5594][])

- [run][] now gives an error result rather than a warning when an
  input glob doesn't match.  ([#5594][])

- The procedure for creating a RIA store checks for an existing
  ria-layout-version file and makes sure its version matches the
  desired version.  This check wasn't done correctly for SSH hosts.
  ([#5607][])

- A helper for transforming git-annex JSON records into DataLad
  results didn't account for the unusual case where the git-annex
  record doesn't have a "file" key.  ([#5580][])

- The test suite required updates for recent changes in PyGithub and
  git-annex.  ([#5603][]) ([#5609][])

## Enhancements and new features

- The DataLad source repository has long had a
  tools/cmdline-completion helper.  This functionality is now exposed
  as a command, `datalad shell-completion`.  ([#5544][])


# 0.14.2 (April 14, 2021) -- .

## Fixes

- [push][] now works bottom-up, pushing submodules first so that hooks
  on the remote can aggregate updated subdataset information. ([#5416][])

- [run-procedure][] didn't ensure that the configuration of
  subdatasets was reloaded.  ([#5552][])


# 0.14.1 (April 01, 2021) -- .

## Fixes

- The recent default branch changes on GitHub's side can lead to
  "git-annex" being selected over "master" as the default branch on
  GitHub when setting up a sibling with [create-sibling-github][].  To
  work around this, the current branch is now pushed first.
  ([#5010][])

- The logic for reading in a JSON line from git-annex failed if the
  response exceeded the buffer size (256 KB on *nix systems).

- Calling [unlock][] with a path of "." from within an untracked
  subdataset incorrectly aborted, complaining that the "dataset
  containing given paths is not underneath the reference dataset".
  ([#5458][])

- [clone][] didn't account for the possibility of multiple accessible
  ORA remotes or the fact that none of them may be associated with the
  RIA store being cloned.  ([#5488][])

- [create-sibling-ria][] didn't call `git update-server-info` after
  setting up the remote repository and, as a result, the repository
  couldn't be fetched until something else (e.g., a push) triggered a
  call to `git update-server-info`.  ([#5531][])

- The parser for git-config output didn't properly handle multi-line
  values and got thrown off by unexpected and unrelated lines.  ([#5509][])

- The 0.14 release introduced regressions in the handling of progress
  bars for git-annex actions, including collapsing progress bars for
  concurrent operations.  ([#5421][]) ([#5438][])

- [save][] failed if the user configured Git's `diff.ignoreSubmodules`
  to a non-default value.  ([#5453][])

- A interprocess lock is now used to prevent a race between checking
  for an SSH socket's existence and creating it.  ([#5466][])

- If a Python procedure script is executable, [run-procedure][]
  invokes it directly rather than passing it to `sys.executable`.  The
  non-executable Python procedures that ship with DataLad now include
  shebangs so that invoking them has a chance of working on file
  systems that present all files as executable.  ([#5436][])

- DataLad's wrapper around `argparse` failed if an underscore was used
  in a positional argument.  ([#5525][])

## Enhancements and new features

- DataLad's method for mapping environment variables to configuration
  options (e.g., `DATALAD_FOO_X__Y` to `datalad.foo.x-y`) doesn't work
  if the subsection name ("FOO") has an underscore.  This limitation
  can be sidestepped with the new `DATALAD_CONFIG_OVERRIDES_JSON`
  environment variable, which can be set to a JSON record of
  configuration values.  ([#5505][])


# 0.14.0 (February 02, 2021) -- .

## Major refactoring and deprecations

- Git versions below v2.19.1 are no longer supported.  ([#4650][])

- The minimum git-annex version is still 7.20190503, but, if you're on
  Windows (or use adjusted branches in general), please upgrade to at
  least 8.20200330 but ideally 8.20210127 to get subdataset-related
  fixes.  ([#4292][]) ([#5290][])

- The minimum supported version of Python is now 3.6.  ([#4879][])

- [publish][] is now deprecated in favor of [push][].  It will be
  removed in the 0.15.0 release at the earliest.

- A new command runner was added in v0.13.  Functionality related to
  the old runner has now been removed: `Runner`, `GitRunner`, and
  `run_gitcommand_on_file_list_chunks` from the `datalad.cmd` module
  along with the `datalad.tests.protocolremote`,
  `datalad.cmd.protocol`, and `datalad.cmd.protocol.prefix`
  configuration options.  ([#5229][])

- The `--no-storage-sibling` switch of `create-sibling-ria` is
  deprecated in favor of `--storage-sibling=off` and will be removed
  in a later release.  ([#5090][])

- The `get_git_dir` static method of `GitRepo` is deprecated and will
  be removed in a later release.  Use the `dot_git` attribute of an
  instance instead.  ([#4597][])

- The `ProcessAnnexProgressIndicators` helper from
  `datalad.support.annexrepo` has been removed.  ([#5259][])

- The `save` argument of [install][], a noop since v0.6.0, has been
  dropped.  ([#5278][])

- The `get_URLS` method of `AnnexCustomRemote` is deprecated and will
  be removed in a later release.  ([#4955][])

- `ConfigManager.get` now returns a single value rather than a tuple
  when there are multiple values for the same key, as very few callers
  correctly accounted for the possibility of a tuple return value.
  Callers can restore the old behavior by passing `get_all=True`.
  ([#4924][])

- In 0.12.0, all of the `assure_*` functions in `datalad.utils` were
  renamed as `ensure_*`, keeping the old names around as compatibility
  aliases.  The `assure_*` variants are now marked as deprecated and
  will be removed in a later release.  ([#4908][])

- The `datalad.inteface.run` module, which was deprecated in 0.12.0
  and kept as a compatibility shim for `datalad.core.local.run`, has
  been removed.  ([#4583][])

- The `saver` argument of `datalad.core.local.run.run_command`, marked
  as obsolete in 0.12.0, has been removed.  ([#4583][])

- The `dataset_only` argument of the `ConfigManager` class was
  deprecated in 0.12 and has now been removed.  ([#4828][])

- The `linux_distribution_name`, `linux_distribution_release`, and
  `on_debian_wheezy` attributes in `datalad.utils` are no longer set
  at import time and will be removed in a later release.  Use
  `datalad.utils.get_linux_distribution` instead.  ([#4696][])

- `datalad.distribution.clone`, which was marked as obsolete in v0.12
  in favor of `datalad.core.distributed.clone`, has been removed.
  ([#4904][])

- `datalad.support.annexrepo.N_AUTO_JOBS`, announced as deprecated in
  v0.12.6, has been removed.  ([#4904][])

- The `compat` parameter of `GitRepo.get_submodules`, added in v0.12
  as a temporary compatibility layer, has been removed.  ([#4904][])

- The long-deprecated (and non-functional) `url` parameter of
  `GitRepo.__init__` has been removed.  ([#5342][])

## Fixes

- Cloning onto a system that enters adjusted branches by default (as
  Windows does) did not properly record the clone URL.  ([#5128][])

- The RIA-specific handling after calling [clone][] was correctly
  triggered by `ria+http` URLs but not `ria+https` URLs.  ([#4977][])

- If the registered commit wasn't found when cloning a subdataset, the
  failed attempt was left around.  ([#5391][])

- The remote calls to `cp` and `chmod` in [create-sibling][] were not
  portable and failed on macOS.  ([#5108][])

- A more reliable check is now done to decide if configuration files
  need to be reloaded.  ([#5276][])

- The internal command runner's handling of the event loop has been
  improved to play nicer with outside applications and scripts that
  use asyncio.  ([#5350][]) ([#5367][])

## Enhancements and new features

- The subdataset handling for adjusted branches, which is particularly
  important on Windows where git-annex enters an adjusted branch by
  default, has been improved.  A core piece of the new approach is
  registering the commit of the primary branch, not its checked out
  adjusted branch, in the superdataset.  Note: This means that `git
  status` will always consider a subdataset on an adjusted branch as
  dirty while `datalad status` will look more closely and see if the
  tip of the primary branch matches the registered commit.
  ([#5241][])

- The performance of the [subdatasets][] command has been improved,
  with substantial speedups for recursive processing of many
  subdatasets.  ([#4868][]) ([#5076][])

- Adding new subdatasets via [save][] has been sped up.  ([#4793][])

- [get][], [save][], and [addurls][] gained support for parallel
  operations that can be enabled via the `--jobs` command-line option
  or the new `datalad.runtime.max-jobs` configuration option.  ([#5022][])

- [addurls][]
  - learned how to read data from standard input.  ([#4669][])
  - now supports tab-separated input.  ([#4845][])
  - now lets Python callers pass in a list of records rather than a
    file name.  ([#5285][])
  - gained a `--drop-after` switch that signals to drop a file's
    content after downloading and adding it to the annex.  ([#5081][])
  - is now able to construct a tree of files from known checksums
    without downloading content via its new `--key` option.  ([#5184][])
  - records the URL file in the commit message as provided by the
    caller rather than using the resolved absolute path. ([#5091][])
  - is now speedier.  ([#4867][]) ([#5022][])

- [create-sibling-github][] learned how to create private repositories
  (thanks to Nolan Nichols).  ([#4769][])

- [create-sibling-ria][] gained a `--storage-sibling` option.  When
  `--storage-sibling=only` is specified, the storage sibling is
  created without an accompanying Git sibling.  This enables using
  hosts without Git installed for storage.  ([#5090][])

- The download machinery (and thus the `datalad` special remote)
  gained support for a new scheme, `shub://`, which follows the same
  format used by `singularity run` and friends.  In contrast to the
  short-lived URLs obtained by querying Singularity Hub directly,
  `shub://` URLs are suitable for registering with git-annex.  ([#4816][])

- A provider is now included for https://registry-1.docker.io URLs.
  This is useful for storing an image's blobs in a dataset and
  registering the URLs with git-annex.  ([#5129][])

- The `add-readme` command now links to the [DataLad
  handbook][handbook] rather than <http://docs.datalad.org>.  ([#4991][])

- New option `datalad.locations.extra-procedures` specifies an
  additional location that should be searched for procedures.  ([#5156][])

- The class for handling configuration values, `ConfigManager`, now
  takes a lock before writes to allow for multiple processes to modify
  the configuration of a dataset.  ([#4829][])

- [clone][] now records the original, unresolved URL for a subdataset
  under `submodule.<name>.datalad-url` in the parent's .gitmodules,
  enabling later [get][] calls to use the original URL.  This is
  particularly useful for `ria+` URLs.  ([#5346][])

- Installing a subdataset now uses custom handling rather than calling
  `git submodule update --init`.  This avoids some locking issues when
  running [get][] in parallel and enables more accurate source URLs to
  be recorded.  ([#4853][])

- `GitRepo.get_content_info`, a helper that gets triggered by many
  commands, got faster by tweaking its `git ls-files` call.  ([#5067][])

- [wtf][] now includes credentials-related information (e.g. active
  backends) in the its output.  ([#4982][])

- The `call_git*` methods of `GitRepo` now have a `read_only`
  parameter.  Callers can set this to `True` to promise that the
  provided command does not write to the repository, bypassing the
  cost of some checks and locking.  ([#5070][])

- New `call_annex*` methods in the `AnnexRepo` class provide an
  interface for running git-annex commands similar to that of the
  `GitRepo.call_git*` methods.  ([#5163][])

- It's now possible to register a custom metadata indexer that is
  discovered by [search][] and used to generate an index.  ([#4963][])

- The `ConfigManager` methods `get`, `getbool`, `getfloat`, and
  `getint` now return a single value (with same precedence as `git
  config --get`) when there are multiple values for the same key (in
  the non-committed git configuration, if the key is present there, or
  in the dataset configuration).  For `get`, the old behavior can be
  restored by specifying `get_all=True`.  ([#4924][])

- Command-line scripts are now defined via the `entry_points` argument
  of `setuptools.setup` instead of the `scripts` argument.  ([#4695][])

- Interactive use of `--help` on the command-line now invokes a pager
  on more systems and installation setups.  ([#5344][])

- The `datalad` special remote now tries to eliminate some unnecessary
  interactions with git-annex by being smarter about how it queries
  for URLs associated with a key.  ([#4955][])

- The `GitRepo` class now does a better job of handling bare
  repositories, a step towards bare repositories support in DataLad.
  ([#4911][])

- More internal work to move the code base over to the new command
  runner.  ([#4699][]) ([#4855][]) ([#4900][]) ([#4996][]) ([#5002][])
  ([#5141][]) ([#5142][]) ([#5229][])


# 0.13.7 (January 04, 2021) -- .

## Fixes

- Cloning from a RIA store on the local file system initialized annex
  in the Git sibling of the RIA source, which is problematic because
  all annex-related functionality should go through the storage
  sibling.  [clone][] now sets `remote.origin.annex-ignore` to `true`
  after cloning from RIA stores to prevent this.  ([#5255][])

- [create-sibling][] invoked `cp` in a way that was not compatible
  with macOS.  ([#5269][])

- Due to a bug in older Git versions (before 2.25), calling [status][]
  with a file under .git/ (e.g., `datalad status .git/config`)
  incorrectly reported the file as untracked.  A workaround has been
  added.  ([#5258][])

- Update tests for compatibility with latest git-annex.  ([#5254][])

## Enhancements and new features

- [copy-file][] now aborts if .git/ is in the target directory, adding
  to its existing .git/ safety checks.  ([#5258][])


# 0.13.6 (December 14, 2020) -- .

## Fixes

- An assortment of fixes for Windows compatibility.  ([#5113][]) ([#5119][])
  ([#5125][]) ([#5127][]) ([#5136][]) ([#5201][]) ([#5200][]) ([#5214][])

- Adding a subdataset on a system that defaults to using an adjusted
  branch (i.e. doesn't support symlinks) didn't properly set up the
  submodule URL if the source dataset was not in an adjusted state.
  ([#5127][])

- [push][] failed to push to a remote that did not have an
  `annex-uuid` value in the local `.git/config`.  ([#5148][])

- The default renderer has been improved to avoid a spurious leading
  space, which led to the displayed path being incorrect in some
  cases.  ([#5121][])

- [siblings][] showed an uninformative error message when asked to
  configure an unknown remote.  ([#5146][])

- [drop][] confusingly relayed a suggestion from `git annex drop` to
  use `--force`, an option that does not exist in `datalad drop`.
  ([#5194][])

- [create-sibling-github][] no longer offers user/password
  authentication because it is no longer supported by GitHub.
  ([#5218][])

- The internal command runner's handling of the event loop has been
  tweaked to hopefully fix issues with runnning DataLad from IPython.
  ([#5106][])

- SSH cleanup wasn't reliably triggered by the ORA special remote on
  failure, leading to a stall with a particular version of git-annex,
  8.20201103.  (This is also resolved on git-annex's end as of
  8.20201127.)  ([#5151][])

## Enhancements and new features

- The credential helper no longer asks the user to repeat tokens or
  AWS keys.  ([#5219][])

- The new option `datalad.locations.sockets` controls where Datalad
  stores SSH sockets, allowing users to more easily work around file
  system and path length restrictions.  ([#5238][])

# 0.13.5 (October 30, 2020) -- .

## Fixes

- SSH connection handling has been reworked to fix cloning on Windows.
  A new configuration option, `datalad.ssh.multiplex-connections`,
  defaults to false on Windows.  ([#5042][])

- The ORA special remote and post-clone RIA configuration now provide
  authentication via DataLad's credential mechanism and better
  handling of HTTP status codes.  ([#5025][]) ([#5026][])

- By default, if a git executable is present in the same location as
  git-annex, DataLad modifies `PATH` when running git and git-annex so
  that the bundled git is used.  This logic has been tightened to
  avoid unnecessarily adjusting the path, reducing the cases where the
  adjustment interferes with the local environment, such as special
  remotes in a virtual environment being masked by the system-wide
  variants.  ([#5035][])

- git-annex is now consistently invoked as "git annex" rather than
  "git-annex" to work around failures on Windows.  ([#5001][])

- [push][] called `git annex sync ...` on plain git repositories.
  ([#5051][])

- [save][] in genernal doesn't support registering multiple levels of
  untracked subdatasets, but it can now properly register nested
  subdatasets when all of the subdataset paths are passed explicitly
  (e.g., `datalad save -d. sub-a sub-a/sub-b`).  ([#5049][])

- When called with `--sidecar` and `--explicit`, [run][] didn't save
  the sidecar.  ([#5017][])

- A couple of spots didn't properly quote format fields when combining
  substrings into a format string.  ([#4957][])

- The default credentials configured for `indi-s3` prevented anonymous
  access.  ([#5045][])

## Enhancements and new features

- Messages about suppressed similar results are now rate limited to
  improve performance when there are many similar results coming
  through quickly.  ([#5060][])

- [create-sibling-github][] can now be told to replace an existing
  sibling by passing `--existing=replace`.  ([#5008][])

- Progress bars now react to changes in the terminal's width (requires
  tqdm 2.1 or later).  ([#5057][])


# 0.13.4 (October 6, 2020) -- .

## Fixes

- Ephemeral clones mishandled bare repositories.  ([#4899][])

- The post-clone logic for configuring RIA stores didn't consider
  `https://` URLs.  ([#4977][])

- DataLad custom remotes didn't escape newlines in messages sent to
  git-annex.  ([#4926][])

- The datalad-archives special remote incorrectly treated file names
  as percent-encoded.  ([#4953][])

- The result handler didn't properly escape "%" when constructing its
  message template.  ([#4953][])

- In v0.13.0, the tailored rendering for specific subtypes of external
  command failures (e.g., "out of space" or "remote not available")
  was unintentionally switched to the default rendering.  ([#4966][])

- Various fixes and updates for the NDA authenticator.  ([#4824][])

- The helper for getting a versioned S3 URL did not support anonymous
  access or buckets with "." in their name.  ([#4985][])

- Several issues with the handling of S3 credentials and token
  expiration have been addressed.  ([#4927][]) ([#4931][]) ([#4952][])

## Enhancements and new features

- A warning is now given if the detected Git is below v2.13.0 to let
  users that run into problems know that their Git version is likely
  the culprit.  ([#4866][])

- A fix to [push][] in v0.13.2 introduced a regression that surfaces
  when `push.default` is configured to "matching" and prevents the
  git-annex branch from being pushed.  Note that, as part of the fix,
  the current branch is now always pushed even when it wouldn't be
  based on the configured refspec or `push.default` value. ([#4896][])

- [publish][]
  - now allows spelling the empty string value of `--since=` as `^`
    for consistency with [push][].  ([#4683][])
  - compares a revision given to `--since=` with `HEAD` rather than
    the working tree to speed up the operation.  ([#4448][])

- [rerun][]
  - emits more INFO-level log messages.  ([#4764][])
  - provides better handling of adjusted branches and aborts with a
    clear error for cases that are not supported.  ([#5328][])

- The archives are handled with p7zip, if available, since DataLad
  v0.12.0.  This implementation now supports .tgz and .tbz2 archives.
  ([#4877][])


# 0.13.3 (August 28, 2020) -- .

## Fixes

- Work around a Python bug that led to our asyncio-based command
  runner intermittently failing to capture the output of commands that
  exit very quickly.  ([#4835][])

- [push][] displayed an overestimate of the transfer size when
  multiple files pointed to the same key.  ([#4821][])

- When [download-url][] calls `git annex addurl`, it catches and
  reports any failures rather than crashing.  A change in v0.12.0
  broke this handling in a particular case.  ([#4817][])

## Enhancements and new features

- The wrapper functions returned by decorators are now given more
  meaningful names to hopefully make tracebacks easier to digest.
  ([#4834][])


# 0.13.2 (August 10, 2020) -- .

## Deprecations

- The `allow_quick` parameter of `AnnexRepo.file_has_content` and
  `AnnexRepo.is_under_annex` is now ignored and will be removed in a
  later release.  This parameter was only relevant for git-annex
  versions before 7.20190912.  ([#4736][])

## Fixes

- Updates for compatibility with recent git and git-annex releases.
  ([#4746][]) ([#4760][]) ([#4684][])

- [push][] didn't sync the git-annex branch when `--data=nothing` was
  specified.  ([#4786][])

- The `datalad.clone.reckless` configuration wasn't stored in
  non-annex datasets, preventing the values from being inherited by
  annex subdatasets.  ([#4749][])

- Running the post-update hook installed by `create-sibling --ui`
  could overwrite web log files from previous runs in the unlikely
  event that the hook was executed multiple times in the same second.
  ([#4745][])

- [clone][] inspected git's standard error in a way that could cause
  an attribute error.  ([#4775][])

- When cloning a repository whose `HEAD` points to a branch without
  commits, [clone][] tries to find a more useful branch to check out.
  It unwisely considered adjusted branches.  ([#4792][])

- Since v0.12.0, `SSHManager.close` hasn't closed connections when the
  `ctrl_path` argument was explicitly given.  ([#4757][])

- When working in a dataset in which `git annex init` had not yet been
  called, the `file_has_content` and `is_under_annex` methods of
  `AnnexRepo` incorrectly took the "allow quick" code path on file
  systems that did not support it ([#4736][])

## Enhancements

- [create][] now assigns version 4 (random) UUIDs instead of version 1
  UUIDs that encode the time and hardware address.  ([#4790][])

- The documentation for [create][] now does a better job of describing
  the interaction between `--dataset` and `PATH`.  ([#4763][])

- The `format_commit` and `get_hexsha` methods of `GitRepo` have been
  sped up.  ([#4807][]) ([#4806][])

- A better error message is now shown when the `^` or `^.` shortcuts
  for `--dataset` do not resolve to a dataset.  ([#4759][])

- A more helpful error message is now shown if a caller tries to
  download an `ftp://` link but does not have `request_ftp` installed.
  ([#4788][])

- [clone][] now tries harder to get up-to-date availability
  information after auto-enabling `type=git` special remotes.  ([#2897][])


# 0.13.1 (July 17, 2020) -- .

## Fixes

- Cloning a subdataset should inherit the parent's
  `datalad.clone.reckless` value, but that did not happen when cloning
  via `datalad get` rather than `datalad install` or `datalad clone`.
  ([#4657][])

- The default result renderer crashed when the result did not have a
  `path` key.  ([#4666][]) ([#4673][])

- `datalad push` didn't show information about `git push` errors when
  the output was not in the format that it expected.  ([#4674][])

- `datalad push` silently accepted an empty string for `--since` even
  though it is an invalid value.  ([#4682][])

- Our JavaScript testing setup on Travis grew stale and has now been
  updated.  (Thanks to Xiao Gui.)  ([#4687][])

- The new class for running Git commands (added in v0.13.0) ignored
  any changes to the process environment that occurred after
  instantiation.  ([#4703][])

## Enhancements and new features

- `datalad push` now avoids unnecessary `git push` dry runs and pushes
  all refspecs with a single `git push` call rather than invoking `git
  push` for each one.  ([#4692][]) ([#4675][])

- The readability of SSH error messages has been improved.  ([#4729][])

- `datalad.support.annexrepo` avoids calling
  `datalad.utils.get_linux_distribution` at import time and caches the
  result once it is called because, as of Python 3.8, the function
  uses `distro` underneath, adding noticeable overhead.  ([#4696][])

  Third-party code should be updated to use `get_linux_distribution`
  directly in the unlikely event that the code relied on the
  import-time call to `get_linux_distribution` setting the
  `linux_distribution_name`, `linux_distribution_release`, or
  `on_debian_wheezy` attributes in `datalad.utils.


# 0.13.0 (June 23, 2020) -- .

A handful of new commands, including `copy-file`, `push`, and
`create-sibling-ria`, along with various fixes and enhancements

## Major refactoring and deprecations

- The `no_annex` parameter of [create][], which is exposed in the
  Python API but not the command line, is deprecated and will be
  removed in a later release.  Use the new `annex` argument instead,
  flipping the value.  Command-line callers that use `--no-annex` are
  unaffected.  ([#4321][])

- `datalad add`, which was deprecated in 0.12.0, has been removed.
  ([#4158][]) ([#4319][])

- The following `GitRepo` and `AnnexRepo` methods have been removed:
  `get_changed_files`, `get_missing_files`, and `get_deleted_files`.
  ([#4169][]) ([#4158][])

- The `get_branch_commits` method of `GitRepo` and `AnnexRepo` has
  been renamed to `get_branch_commits_`.  ([#3834][])

- The custom `commit` method of `AnnexRepo` has been removed, and
  `AnnexRepo.commit` now resolves to the parent method,
  `GitRepo.commit`.  ([#4168][])

- GitPython's `git.repo.base.Repo` class is no longer available via
  the `.repo` attribute of `GitRepo` and `AnnexRepo`.  ([#4172][])

- `AnnexRepo.get_corresponding_branch` now returns `None` rather than
  the current branch name when a managed branch is not checked out.
  ([#4274][])

- The special UUID for git-annex web remotes is now available as
  `datalad.consts.WEB_SPECIAL_REMOTE_UUID`.  It remains accessible as
  `AnnexRepo.WEB_UUID` for compatibility, but new code should use
  `consts.WEB_SPECIAL_REMOTE_UUID` ([#4460][]).

## Fixes

- Widespread improvements in functionality and test coverage on
  Windows and crippled file systems in general.  ([#4057][])
  ([#4245][]) ([#4268][]) ([#4276][]) ([#4291][]) ([#4296][])
  ([#4301][]) ([#4303][]) ([#4304][]) ([#4305][]) ([#4306][])

- `AnnexRepo.get_size_from_key` incorrectly handled file chunks.
  ([#4081][])

- [create-sibling][] would too readily clobber existing paths when
  called with `--existing=replace`.  It now gets confirmation from the
  user before doing so if running interactively and unconditionally
  aborts when running non-interactively.  ([#4147][])

- [update][]  ([#4159][])
  - queried the incorrect branch configuration when updating non-annex
    repositories.
  - didn't account for the fact that the local repository can be
    configured as the upstream "remote" for a branch.

- When the caller included `--bare` as a `git init` option, [create][]
  crashed creating the bare repository, which is currently
  unsupported, rather than aborting with an informative error message.
  ([#4065][])

- The logic for automatically propagating the 'origin' remote when
  cloning a local source could unintentionally trigger a fetch of a
  non-local remote.  ([#4196][])

- All remaining `get_submodules()` call sites that relied on the
  temporary compatibility layer added in v0.12.0 have been updated.
  ([#4348][])

- The custom result summary renderer for [get][], which was visible
  with `--output-format=tailored`, displayed incorrect and confusing
  information in some cases.  The custom renderer has been removed
  entirely.  ([#4471][])

- The documentation for the Python interface of a command listed an
  incorrect default when the command overrode the value of command
  parameters such as `result_renderer`.  ([#4480][])

## Enhancements and new features

- The default result renderer learned to elide a chain of results
  after seeing ten consecutive results that it considers similar,
  which improves the display of actions that have many results (e.g.,
  saving hundreds of files).  ([#4337][])

- The default result renderer, in addition to "tailored" result
  renderer, now triggers the custom summary renderer, if any.  ([#4338][])

- The new command [create-sibling-ria][] provides support for creating
  a sibling in a [RIA store][handbook-scalable-datastore]. ([#4124][])

- DataLad ships with a new special remote, git-annex-remote-ora, for
  interacting with [RIA stores][handbook-scalable-datastore] and a new
  command [export-archive-ora][] for exporting an archive from a local
  annex object store.  ([#4260][]) ([#4203][])

- The new command [push][] provides an alternative interface to
  [publish][] for pushing a dataset hierarchy to a sibling.
  ([#4206][]) ([#4581][]) ([#4617][]) ([#4620][])

- The new command [copy-file][] copies files and associated
  availability information from one dataset to another.  ([#4430][])

- The command examples have been expanded and improved.  ([#4091][])
  ([#4314][]) ([#4464][])

- The tooling for linking to the [DataLad Handbook][handbook] from
  DataLad's documentation has been improved.  ([#4046][])

- The `--reckless` parameter of [clone][] and [install][] learned two
  new modes:
  - "ephemeral", where the .git/annex/ of the cloned repository is
    symlinked to the local source repository's.  ([#4099][])
  - "shared-{group|all|...}" that can be used to set up datasets for
    collaborative write access.  ([#4324][])

- [clone][]
  - learned to handle dataset aliases in RIA stores when given a URL
    of the form `ria+<protocol>://<storelocation>#~<aliasname>`.
    ([#4459][])
  - now checks `datalad.get.subdataset-source-candidate-NAME` to see
    if `NAME` starts with three digits, which is taken as a "cost".
    Sources with lower costs will be tried first.  ([#4619][])

- [update][]  ([#4167][])
  - learned to disallow non-fast-forward updates when `ff-only` is
    given to the `--merge` option.
  - gained a `--follow` option that controls how `--merge` behaves,
    adding support for merging in the revision that is registered in
    the parent dataset rather than merging in the configured branch
    from the sibling.
  - now provides a result record for merge events.

- [create-sibling][] now supports local paths as targets in addition
  to SSH URLs.  ([#4187][])

- [siblings][] now
  - shows a warning if the caller requests to delete a sibling that
    does not exist.  ([#4257][])
  - phrases its warning about non-annex repositories in a less
    alarming way.  ([#4323][])

- The rendering of command errors has been improved.  ([#4157][])

- [save][] now
  - displays a message to signal that the working tree is clean,
    making it more obvious that no results being rendered corresponds
    to a clean state.  ([#4106][])
  - provides a stronger warning against using `--to-git`.  ([#4290][])

- [diff][] and [save][] learned about scenarios where they could avoid
  unnecessary and expensive work.  ([#4526][]) ([#4544][]) ([#4549][])

- Calling [diff][] without `--recursive` but with a path constraint
  within a subdataset ("<subdataset>/<path>") now traverses into the
  subdataset, as "<subdataset>/" would, restricting its report to
  "<subdataset>/<path>".  ([#4235][])

- New option `datalad.annex.retry` controls how many times git-annex
  will retry on a failed transfer.  It defaults to 3 and can be set to
  0 to restore the previous behavior.  ([#4382][])

- [wtf][] now warns when the specified dataset does not exist.
  ([#4331][])

- The `repr` and `str` output of the dataset and repo classes got a
  facelift.  ([#4420][]) ([#4435][]) ([#4439][])

- The DataLad Singularity container now comes with p7zip-full.

- DataLad emits a log message when the current working directory is
  resolved to a different location due to a symlink.  This is now
  logged at the DEBUG rather than WARNING level, as it typically does
  not indicate a problem.  ([#4426][])

- DataLad now lets the caller know that `git annex init` is scanning
  for unlocked files, as this operation can be slow in some
  repositories.  ([#4316][])

- The `log_progress` helper learned how to set the starting point to a
  non-zero value and how to update the total of an existing progress
  bar, two features needed for planned improvements to how some
  commands display their progress.  ([#4438][])

- The `ExternalVersions` object, which is used to check versions of
  Python modules and external tools (e.g., git-annex), gained an `add`
  method that enables DataLad extensions and other third-party code to
  include other programs of interest.  ([#4441][])

- All of the remaining spots that use GitPython have been rewritten
  without it.  Most notably, this includes rewrites of the `clone`,
  `fetch`, and `push` methods of `GitRepo`.  ([#4080][]) ([#4087][])
  ([#4170][]) ([#4171][]) ([#4175][]) ([#4172][])

- When `GitRepo.commit` splits its operation across multiple calls to
  avoid exceeding the maximum command line length, it now amends to
  initial commit rather than creating multiple commits.  ([#4156][])

- `GitRepo` gained a `get_corresponding_branch` method (which always
   returns None), allowing a caller to invoke the method without
   needing to check if the underlying repo class is `GitRepo` or
   `AnnexRepo`.  ([#4274][])

- A new helper function `datalad.core.local.repo.repo_from_path`
  returns a repo class for a specified path.  ([#4273][])

- New `AnnexRepo` method `localsync` performs a `git annex sync` that
  disables external interaction and is particularly useful for
  propagating changes on an adjusted branch back to the main branch.
  ([#4243][])


# 0.12.7 (May 22, 2020) -- .

## Fixes

- Requesting tailored output (`--output=tailored`) from a command with
  a custom result summary renderer produced repeated output. ([#4463][])

- A longstanding regression in argcomplete-based command-line
  completion for Bash has been fixed.  You can enable completion by
  configuring a Bash startup file to run `eval
  "$(register-python-argcomplete datalad)"` or source DataLad's
  `tools/cmdline-completion`.  The latter should work for Zsh as well.
  ([#4477][])

- [publish][] didn't prevent `git-fetch` from recursing into
  submodules, leading to a failure when the registered submodule was
  not present locally and the submodule did not have a remote named
  'origin'.  ([#4560][])

- [addurls][] botched path handling when the file name format started
  with "./" and the call was made from a subdirectory of the dataset.
  ([#4504][])

- Double dash options in manpages were unintentionally escaped.
  ([#4332][])

- The check for HTTP authentication failures crashed in situations
  where content came in as bytes rather than unicode.  ([#4543][])

- A check in `AnnexRepo.whereis` could lead to a type error.  ([#4552][])

- When installing a dataset to obtain a subdataset, [get][]
  confusingly displayed a message that described the containing
  dataset as "underneath" the subdataset.  ([#4456][])

- A couple of Makefile rules didn't properly quote paths.  ([#4481][])

- With DueCredit support enabled (`DUECREDIT_ENABLE=1`), the query for
  metadata information could flood the output with warnings if
  datasets didn't have aggregated metadata.  The warnings are now
  silenced, with the overall failure of a [metadata][] call logged at
  the debug level.  ([#4568][])

## Enhancements and new features

- The resource identifier helper learned to recognize URLs with
  embedded Git transport information, such as
  gcrypt::https://example.com.  ([#4529][])

- When running non-interactively, a more informative error is now
  signaled when the UI backend, which cannot display a question, is
  asked to do so.  ([#4553][])


# 0.12.6 (April 23, 2020) -- .

## Major refactoring and deprecations

- The value of `datalad.support.annexrep.N_AUTO_JOBS` is no longer
  considered.  The variable will be removed in a later release.
  ([#4409][])

## Fixes

- Staring with v0.12.0, `datalad save` recorded the current branch of
  a parent dataset as the `branch` value in the .gitmodules entry for
  a subdataset.  This behavior is problematic for a few reasons and
  has been reverted.  ([#4375][])

- The default for the `--jobs` option, "auto", instructed DataLad to
  pass a value to git-annex's `--jobs` equal to `min(8, max(3, <number
  of CPUs>))`, which could lead to issues due to the large number of
  child processes spawned and file descriptors opened.  To avoid this
  behavior, `--jobs=auto` now results in git-annex being called with
  `--jobs=1` by default.  Configure the new option
  `datalad.runtime.max-annex-jobs` to control the maximum value that
  will be considered when `--jobs='auto'`.  ([#4409][])

- Various commands have been adjusted to better handle the case where
  a remote's HEAD ref points to an unborn branch.  ([#4370][])

- [search]
  - learned to use the query as a regular expression that restricts
    the keys that are shown for `--show-keys short`. ([#4354][])
  - gives a more helpful message when query is an invalid regular
    expression.  ([#4398][])

- The code for parsing Git configuration did not follow Git's behavior
  of accepting a key with no value as shorthand for key=true.  ([#4421][])

- `AnnexRepo.info` needed a compatibility update for a change in how
  git-annex reports file names.  ([#4431][])

- [create-sibling-github][] did not gracefully handle a token that did
  not have the necessary permissions.  ([#4400][])

## Enhancements and new features

- [search] learned to use the query as a regular expression that
  restricts the keys that are shown for `--show-keys short`. ([#4354][])

- `datalad <subcommand>` learned to point to the [datalad-container][]
  extension when a subcommand from that extension is given but the
  extension is not installed.  ([#4400][]) ([#4174][])


# 0.12.5 (Apr 02, 2020) -- a small step for datalad ...
￼
Fix some bugs and make the world an even better place.

## Fixes

- Our `log_progress` helper mishandled the initial display and step of
  the progress bar.  ([#4326][])

- `AnnexRepo.get_content_annexinfo` is designed to accept `init=None`,
  but passing that led to an error.  ([#4330][])

- Update a regular expression to handle an output change in Git
  v2.26.0.  ([#4328][])

- We now set `LC_MESSAGES` to 'C' while running git to avoid failures
  when parsing output that is marked for translation.  ([#4342][])

- The helper for decoding JSON streams loaded the last line of input
  without decoding it if the line didn't end with a new line, a
  regression introduced in the 0.12.0 release.  ([#4361][])

- The clone command failed to git-annex-init a fresh clone whenever
  it considered to add the origin of the origin as a remote.  ([#4367][])


# 0.12.4 (Mar 19, 2020) -- Windows?!
￼
The main purpose of this release is to have one on PyPi that has no
associated wheel to enable a working installation on Windows ([#4315][]).

## Fixes

- The description of the `log.outputs` config switch did not keep up
  with code changes and incorrectly stated that the output would be
  logged at the DEBUG level; logging actually happens at a lower
  level.  ([#4317][])

# 0.12.3 (March 16, 2020) -- .

Updates for compatibility with the latest git-annex, along with a few
miscellaneous fixes

## Major refactoring and deprecations

- All spots that raised a `NoDatasetArgumentFound` exception now raise
  a `NoDatasetFound` exception to better reflect the situation: it is
  the _dataset_ rather than the _argument_ that is not found.  For
  compatibility, the latter inherits from the former, but new code
  should prefer the latter.  ([#4285][])

## Fixes

- Updates for compatibility with git-annex version 8.20200226. ([#4214][])

- `datalad export-to-figshare` failed to export if the generated title
  was fewer than three characters.  It now queries the caller for the
  title and guards against titles that are too short.  ([#4140][])

- Authentication was requested multiple times when git-annex launched
  parallel downloads from the `datalad` special remote. ([#4308][])

- At verbose logging levels, DataLad requests that git-annex display
  debugging information too.  Work around a bug in git-annex that
  prevented that from happening.  ([#4212][])

- The internal command runner looked in the wrong place for some
  configuration variables, including `datalad.log.outputs`, resulting
  in the default value always being used.  ([#4194][])

- [publish][] failed when trying to publish to a git-lfs special
  remote for the first time.  ([#4200][])

- `AnnexRepo.set_remote_url` is supposed to establish shared SSH
  connections but failed to do so.  ([#4262][])

## Enhancements and new features

- The message provided when a command cannot determine what dataset to
  operate on has been improved.  ([#4285][])

- The "aws-s3" authentication type now allows specifying the host
  through "aws-s3_host", which was needed to work around an
  authorization error due to a longstanding upstream bug.  ([#4239][])

- The xmp metadata extractor now recognizes ".wav" files.


# 0.12.2 (Jan 28, 2020) -- Smoothen the ride

Mostly a bugfix release with various robustifications, but also makes
the first step towards versioned dataset installation requests.

## Major refactoring and deprecations

- The minimum required version for GitPython is now 2.1.12. ([#4070][])

## Fixes

- The class for handling configuration values, `ConfigManager`,
  inappropriately considered the current working directory's dataset,
  if any, for both reading and writing when instantiated with
  `dataset=None`.  This misbehavior is fairly inaccessible through
  typical use of DataLad.  It affects `datalad.cfg`, the top-level
  configuration instance that should not consider repository-specific
  values.  It also affects Python users that call `Dataset` with a
  path that does not yet exist and persists until that dataset is
  created. ([#4078][])

- [update][] saved the dataset when called with `--merge`, which is
  unnecessary and risks committing unrelated changes.  ([#3996][])

- Confusing and irrelevant information about Python defaults have been
  dropped from the command-line help.  ([#4002][])

- The logic for automatically propagating the 'origin' remote when
  cloning a local source didn't properly account for relative paths.
  ([#4045][])

- Various fixes to file name handling and quoting on Windows.
  ([#4049][]) ([#4050][])

- When cloning failed, error lines were not bubbled up to the user in
  some scenarios.  ([#4060][])

## Enhancements and new features

- [clone][] (and thus [install][])
  - now propagates the `reckless` mode from the superdataset when
    cloning a dataset into it.  ([#4037][])
  - gained support for `ria+<protocol>://` URLs that point to
    [RIA][handbook-scalable-datastore] stores.  ([#4022][])
  - learned to read "@version" from `ria+` URLs and install that
    version of a dataset ([#4036][]) and to apply URL rewrites
    configured through Git's `url.*.insteadOf` mechanism ([#4064][]).
  - now copies `datalad.get.subdataset-source-candidate-<name>`
    options configured within the superdataset into the subdataset.
    This is particularly useful for RIA data stores. ([#4073][])

- Archives are now (optionally) handled with 7-Zip instead of
  `patool`.  7-Zip will be used by default, but `patool` will be used
  on non-Windows systems if the `datalad.runtime.use-patool` option is
  set or the `7z` executable is not found.  ([#4041][])


# 0.12.1 (Jan 15, 2020) -- Small bump after big bang

Fix some fallout after major release.

## Fixes

- Revert incorrect relative path adjustment to URLs in [clone][]. ([#3538][])

- Various small fixes to internal helpers and test to run on Windows
  ([#2566][]) ([#2534][])

# 0.12.0 (Jan 11, 2020) -- Krakatoa

This release is the result of more than a year of development that includes
fixes for a large number of issues, yielding more robust behavior across a
wider range of use cases, and introduces major changes in API and behavior. It
is the first release for which extensive user documentation is available in a
dedicated [DataLad Handbook][handbook].  Python 3 (3.5 and later) is now the
only supported Python flavor.

## Major changes 0.12 vs 0.11

- [save][] fully replaces [add][] (which is obsolete now, and will be removed
  in a future release).

- A new Git-annex aware [status][] command enables detailed inspection of dataset
  hierarchies. The previously available [diff][] command has been adjusted to
  match [status][] in argument semantics and behavior.

- The ability to configure dataset procedures prior and after the execution of
  particular commands has been replaced by a flexible "hook" mechanism that is able
  to run arbitrary DataLad commands whenever command results are detected that match
  a specification.

- Support of the Windows platform has been improved substantially. While performance
  and feature coverage on Windows still falls behind Unix-like systems, typical data
  consumer use cases, and standard dataset operations, such as [create][] and [save][],
  are now working. Basic support for data provenance capture via [run][] is also
  functional.

- Support for Git-annex direct mode repositories has been removed, following the
  end of support in Git-annex itself.

- The semantics of relative paths in command line arguments have changed. Previously,
  a call `datalad save --dataset /tmp/myds some/relpath` would have been interpreted
  as saving a file at `/tmp/myds/some/relpath` into dataset `/tmp/myds`. This has
  changed to saving `$PWD/some/relpath` into dataset `/tmp/myds`. More generally,
  relative paths are now always treated as relative to the current working directory,
  except for path arguments of [Dataset][] class instance methods of the Python API.
  The resulting partial duplication of path specifications between path and dataset
  arguments is mitigated by the introduction of two special symbols that can be given
  as dataset argument: `^` and `^.`, which identify the topmost superdataset and the
  closest dataset that contains the working directory, respectively.

- The concept of a "core API" has been introduced. Commands situated in the module
  `datalad.core` (such as [create][], [save][], [run][], [status][], [diff][])
  receive additional scrutiny regarding API and implementation, and are
  meant to provide longer-term stability. Application developers are encouraged to
  preferentially build on these commands.

## Major refactoring and deprecations since 0.12.0rc6

- [clone][] has been incorporated into the growing core API. The public
  `--alternative-source` parameter has been removed, and a `clone_dataset`
  function with multi-source capabilities is provided instead. The
  `--reckless` parameter can now take literal mode labels instead of just
  beeing a binary flag, but backwards compatibility is maintained.

- The `get_file_content` method of `GitRepo` was no longer used
  internally or in any known DataLad extensions and has been removed.
  ([#3812][])

- The function `get_dataset_root` has been replaced by
  `rev_get_dataset_root`.  `rev_get_dataset_root` remains as a
  compatibility alias and will be removed in a later release.  ([#3815][])

- The `add_sibling` module, marked obsolete in v0.6.0, has been
  removed.  ([#3871][])

- `mock` is no longer declared as an external dependency because we
   can rely on it being in the standard library now that our minimum
   required Python version is 3.5. ([#3860][])

- [download-url][] now requires that directories be indicated with a
  trailing slash rather than interpreting a path as directory when it
  doesn't exist.  This avoids confusion that can result from typos and
  makes it possible to support directory targets that do not exist.
  ([#3854][])

- The `dataset_only` argument of the `ConfigManager` class is
  deprecated.  Use `source="dataset"` instead.  ([#3907][])

- The `--proc-pre` and `--proc-post` options have been removed, and
  configuration values for `datalad.COMMAND.proc-pre` and
  `datalad.COMMAND.proc-post` are no longer honored.  The new result
  hook mechanism provides an alternative for `proc-post`
  procedures. ([#3963][])

## Fixes since 0.12.0rc6

- [publish][] crashed when called with a detached HEAD.  It now aborts
  with an informative message.  ([#3804][])

- Since 0.12.0rc6 the call to [update][] in [siblings][] resulted in a
  spurious warning.  ([#3877][])

- [siblings][] crashed if it encountered an annex repository that was
  marked as dead.  ([#3892][])

- The update of [rerun][] in v0.12.0rc3 for the rewritten [diff][]
  command didn't account for a change in the output of `diff`, leading
  to `rerun --report` unintentionally including unchanged files in its
  diff values.  ([#3873][])

- In 0.12.0rc5 [download-url][] was updated to follow the new path
  handling logic, but its calls to AnnexRepo weren't properly
  adjusted, resulting in incorrect path handling when the called from
  a dataset subdirectory.  ([#3850][])

- [download-url][] called `git annex addurl` in a way that failed to
  register a URL when its header didn't report the content size.
  ([#3911][])

- With Git v2.24.0, saving new subdatasets failed due to a bug in that
  Git release.  ([#3904][])

- With DataLad configured to stop on failure (e.g., specifying
  `--on-failure=stop` from the command line), a failing result record
  was not rendered.  ([#3863][])

- Installing a subdataset yielded an "ok" status in cases where the
  repository was not yet in its final state, making it ineffective for
  a caller to operate on the repository in response to the result.
  ([#3906][])

- The internal helper for converting git-annex's JSON output did not
  relay information from the "error-messages" field.  ([#3931][])

- [run-procedure][] reported relative paths that were confusingly not
  relative to the current directory in some cases.  It now always
  reports absolute paths. ([#3959][])

- [diff][] inappropriately reported files as deleted in some cases
  when `to` was a value other than `None`.  ([#3999][])

- An assortment of fixes for Windows compatibility.  ([#3971][]) ([#3974][])
  ([#3975][]) ([#3976][]) ([#3979][])

- Subdatasets installed from a source given by relative path will now
  have this relative path used as 'url' in their .gitmodules record,
  instead of an absolute path generated by Git. ([#3538][])

- [clone][] will now correctly interpret '~/...' paths as absolute path
  specifications. ([#3958][])

- [run-procedure][] mistakenly reported a directory as a procedure.
  ([#3793][])

- The cleanup for batched git-annex processes has been improved.
  ([#3794][]) ([#3851][])

- The function for adding a version ID to an AWS S3 URL doesn't
  support URLs with an "s3://" scheme and raises a
  `NotImplementedError` exception when it encounters one.  The
  function learned to return a URL untouched if an "s3://" URL comes
  in with a version ID.  ([#3842][])

- A few spots needed to be adjusted for compatibility with git-annex's
  new `--sameas` [feature][gx-sameas], which allows special remotes to
  share a data store. ([#3856][])

- The `swallow_logs` utility failed to capture some log messages due
  to an incompatibility with Python 3.7.  ([#3935][])

- [siblings][]
  - crashed if `--inherit` was passed but the parent dataset did not
    have a remote with a matching name.  ([#3954][])
  - configured the wrong pushurl and annexurl values in some
    cases. ([#3955][])

## Enhancements and new features since 0.12.0rc6

- By default, datasets cloned from local source paths will now get a
  configured remote for any recursively discoverable 'origin' sibling that
  is also available from a local path in order to maximize automatic file
  availability across local annexes. ([#3926][])

- The new [result hooks mechanism][hooks] allows callers to specify,
  via local Git configuration values, DataLad command calls that will
  be triggered in response to matching result records (i.e., what you
  see when you call a command with `-f json_pp`).  ([#3903][])

- The command interface classes learned to use a new `_examples_`
  attribute to render documentation examples for both the Python and
  command-line API.  ([#3821][])

- Candidate URLs for cloning a submodule can now be generated based on
  configured templates that have access to various properties of the
  submodule, including its dataset ID.  ([#3828][])

- DataLad's check that the user's Git identity is configured has been
  sped up and now considers the appropriate environment variables as
  well.  ([#3807][])

- The `tag` method of `GitRepo` can now tag revisions other than
  `HEAD` and accepts a list of arbitrary `git tag` options.
  ([#3787][])

- When `get` clones a subdataset and the subdataset's HEAD differs
  from the commit that is registered in the parent, the active branch
  of the subdataset is moved to the registered commit if the
  registered commit is an ancestor of the subdataset's HEAD commit.
  This handling has been moved to a more central location within
  `GitRepo`, and now applies to any `update_submodule(..., init=True)`
  call.  ([#3831][])

- The output of `datalad -h` has been reformatted to improve
  readability.  ([#3862][])

- [unlock][] has been sped up.  ([#3880][])

- [run-procedure][] learned to provide and render more information
  about discovered procedures, including whether the procedure is
  overridden by another procedure with the same base name.  ([#3960][])

- [save][] now ([#3817][])
  - records the active branch in the superdataset when registering a
    new subdataset.
  - calls `git annex sync` when saving a dataset on an adjusted branch
    so that the changes are brought into the mainline branch.

- [subdatasets][] now aborts when its `dataset` argument points to a
  non-existent dataset.  ([#3940][])

- [wtf][] now
  - reports the dataset ID if the current working directory is
    visiting a dataset.  ([#3888][])
  - outputs entries deterministically.  ([#3927][])

- The `ConfigManager` class
  - learned to exclude ``.datalad/config`` as a source of
    configuration values, restricting the sources to standard Git
    configuration files, when called with `source="local"`.
    ([#3907][])
  - accepts a value of "override" for its `where` argument to allow
    Python callers to more convenient override configuration.
    ([#3970][])

- Commands now accept a `dataset` value of "^."  as shorthand for "the
  dataset to which the current directory belongs".  ([#3242][])

# 0.12.0rc6 (Oct 19, 2019) -- some releases are better than the others

bet we will fix some bugs and make a world even a better place.

## Major refactoring and deprecations

- DataLad no longer supports Python 2.  The minimum supported version
  of Python is now 3.5.  ([#3629][])

- Much of the user-focused content at http://docs.datalad.org has been
  removed in favor of more up to date and complete material available
  in the [DataLad Handbook][handbook].  Going forward, the plan is to
  restrict http://docs.datalad.org to technical documentation geared
  at developers.  ([#3678][])

- [update][] used to allow the caller to specify which dataset(s) to
  update as a `PATH` argument or via the the `--dataset` option; now
  only the latter is supported.  Path arguments only serve to restrict
  which subdataset are updated when operating recursively.
  ([#3700][])

- Result records from a [get][] call no longer have a "state" key.
  ([#3746][])

- [update][] and [get][] no longer support operating on independent
  hierarchies of datasets.  ([#3700][]) ([#3746][])

- The [run][] update in 0.12.0rc4 for the new path resolution logic
  broke the handling of inputs and outputs for calls from a
  subdirectory.  ([#3747][])

- The `is_submodule_modified` method of `GitRepo` as well as two
  helper functions in gitrepo.py, `kwargs_to_options` and
  `split_remote_branch`, were no longer used internally or in any
  known DataLad extensions and have been removed.  ([#3702][])
  ([#3704][])

- The `only_remote` option of `GitRepo.is_with_annex` was not used
  internally or in any known extensions and has been dropped.
  ([#3768][])

- The `get_tags` method of `GitRepo` used to sort tags by committer
  date.  It now sorts them by the tagger date for annotated tags and
  the committer date for lightweight tags.  ([#3715][])

- The `rev_resolve_path` substituted `resolve_path` helper. ([#3797][])


## Fixes

- Correctly handle relative paths in [publish][]. ([#3799][]) ([#3102][])

- Do not errorneously discover directory as a procedure. ([#3793][])

- Correctly extract version from manpage to trigger use of manpages for
  `--help`. ([#3798][])

- The `cfg_yoda` procedure saved all modifications in the repository
  rather than saving only the files it modified.  ([#3680][])

- Some spots in the documentation that were supposed appear as two
  hyphens were incorrectly rendered in the HTML output en-dashs.
  ([#3692][])

- [create][], [install][], and [clone][] treated paths as relative to
  the dataset even when the string form was given, violating the new
  path handling rules.  ([#3749][]) ([#3777][]) ([#3780][])

- Providing the "^" shortcut to `--dataset` didn't work properly when
  called from a subdirectory of a subdataset.  ([#3772][])

- We failed to propagate some errors from git-annex when working with
  its JSON output.  ([#3751][])

- With the Python API, callers are allowed to pass a string or list of
  strings as the `cfg_proc` argument to [create][], but the string
  form was mishandled.  ([#3761][])

- Incorrect command quoting for SSH calls on Windows that rendered
  basic SSH-related functionality (e.g., [sshrun][]) on Windows
  unusable.  ([#3688][])

- Annex JSON result handling assumed platform-specific paths on Windows
  instead of the POSIX-style that is happening across all platforms.
  ([#3719][])

- `path_is_under()` was incapable of comparing Windows paths with different
  drive letters.  ([#3728][])

## Enhancements and new features

- Provide a collection of "public" `call_git*` helpers within GitRepo
  and replace use of "private" and less specific `_git_custom_command`
  calls.  ([#3791][])

- [status][] gained a `--report-filetype`.  Setting it to "raw" can
  give a performance boost for the price of no longer distinguishing
  symlinks that point to annexed content from other symlinks.
  ([#3701][])

- [save][] disables file type reporting by [status][] to improve
  performance.  ([#3712][])

- [subdatasets][] ([#3743][])
  - now extends its result records with a `contains` field that lists
    which `contains` arguments matched a given subdataset.
  - yields an 'impossible' result record when a `contains` argument
    wasn't matched to any of the reported subdatasets.

- [install][] now shows more readable output when cloning fails.
  ([#3775][])

- `SSHConnection` now displays a more informative error message when
  it cannot start the `ControlMaster` process.  ([#3776][])

- If the new configuration option `datalad.log.result-level` is set to
  a single level, all result records will be logged at that level.  If
  you've been bothered by DataLad's double reporting of failures,
  consider setting this to "debug".  ([#3754][])

- Configuration values from `datalad -c OPTION=VALUE ...` are now
  validated to provide better errors.  ([#3695][])

- [rerun][] learned how to handle history with merges.  As was already
  the case when cherry picking non-run commits, re-creating merges may
  results in conflicts, and `rerun` does not yet provide an interface
  to let the user handle these.  ([#2754][])

- The `fsck` method of `AnnexRepo` has been enhanced to expose more
  features of the underlying `git fsck` command.  ([#3693][])

- `GitRepo` now has a `for_each_ref_` method that wraps `git
  for-each-ref`, which is used in various spots that used to rely on
  GitPython functionality.  ([#3705][])

- Do not pretend to be able to work in optimized (`python -O`) mode,
  crash early with an informative message. ([#3803][])

# 0.12.0rc5 (September 04, 2019) -- .

Various fixes and enhancements that bring the 0.12.0 release closer.

## Major refactoring and deprecations

- The two modules below have a new home.  The old locations still
  exist as compatibility shims and will be removed in a future
  release.
  - `datalad.distribution.subdatasets` has been moved to
    `datalad.local.subdatasets` ([#3429][])
  - `datalad.interface.run` has been moved to `datalad.core.local.run`
    ([#3444][])

- The `lock` method of `AnnexRepo` and the `options` parameter of
  `AnnexRepo.unlock` were unused internally and have been removed.
  ([#3459][])

- The `get_submodules` method of `GitRepo` has been rewritten without
  GitPython.  When the new `compat` flag is true (the current
  default), the method returns a value that is compatible with the old
  return value.  This backwards-compatible return value and the
  `compat` flag will be removed in a future release.  ([#3508][])

- The logic for resolving relative paths given to a command has
  changed ([#3435][]).  The new rule is that relative paths are taken
  as relative to the dataset only if a dataset _instance_ is passed by
  the caller.  In all other scenarios they're considered relative to
  the current directory.

  The main user-visible difference from the command line is that using
  the `--dataset` argument does _not_ result in relative paths being
  taken as relative to the specified dataset.  (The undocumented
  distinction between "rel/path" and "./rel/path" no longer exists.)

  All commands under `datalad.core` and `datalad.local`, as well as
  `unlock` and `addurls`, follow the new logic.  The goal is for all
  commands to eventually do so.

## Fixes

- The function for loading JSON streams wasn't clever enough to handle
  content that included a Unicode line separator like
  U2028. ([#3524][])

- When [unlock][] was called without an explicit target (i.e., a
  directory or no paths at all), the call failed if any of the files
  did not have content present.  ([#3459][])

- `AnnexRepo.get_content_info` failed in the rare case of a key
  without size information.  ([#3534][])

- [save][] ignored `--on-failure` in its underlying call to
  [status][].  ([#3470][])

- Calling [remove][] with a subdirectory displayed spurious warnings
  about the subdirectory files not existing.  ([#3586][])

- Our processing of `git-annex --json` output mishandled info messages
  from special remotes.  ([#3546][])

- [create][]
  - didn't bypass the "existing subdataset" check when called with
    `--force` as of 0.12.0rc3 ([#3552][])
  - failed to register the up-to-date revision of a subdataset when
    `--cfg-proc` was used with `--dataset` ([#3591][])

- The base downloader had some error handling that wasn't compatible
  with Python 3.  ([#3622][])

- Fixed a number of Unicode py2-compatibility issues. ([#3602][])

- `AnnexRepo.get_content_annexinfo` did not properly chunk file
  arguments to avoid exceeding the command-line character limit.
  ([#3587][])

## Enhancements and new features

- New command `create-sibling-gitlab` provides an interface for
  creating a publication target on a GitLab instance.  ([#3447][])

- [subdatasets][]  ([#3429][])
  - now supports path-constrained queries in the same manner as
    commands like `save` and `status`
  - gained a `--contains=PATH` option that can be used to restrict the
    output to datasets that include a specific path.
  - now narrows the listed subdatasets to those underneath the current
    directory when called with no arguments

- [status][] learned to accept a plain `--annex` (no value) as
  shorthand for `--annex basic`.  ([#3534][])

- The `.dirty` property of `GitRepo` and `AnnexRepo` has been sped up.
  ([#3460][])

- The `get_content_info` method of `GitRepo`, used by `status` and
  commands that depend on `status`, now restricts its git calls to a
  subset of files, if possible, for a performance gain in repositories
  with many files.  ([#3508][])

- Extensions that do not provide a command, such as those that provide
  only metadata extractors, are now supported.  ([#3531][])

- When calling git-annex with `--json`, we log standard error at the
  debug level rather than the warning level if a non-zero exit is
  expected behavior.  ([#3518][])

- [create][] no longer refuses to create a new dataset in the odd
  scenario of an empty .git/ directory upstairs.  ([#3475][])

- As of v2.22.0 Git treats a sub-repository on an unborn branch as a
  repository rather than as a directory.  Our documentation and tests
  have been updated appropriately.  ([#3476][])

- [addurls][] learned to accept a `--cfg-proc` value and pass it to
  its `create` calls.  ([#3562][])

# 0.12.0rc4 (May 15, 2019) -- the revolution is over

With the replacement of the `save` command implementation with `rev-save`
the revolution effort is now over, and the set of key commands for
local dataset operations (`create`, `run`, `save`, `status`, `diff`) is
 now complete. This new core API is available from `datalad.core.local`
(and also via `datalad.api`, as any other command).
￼
## Major refactoring and deprecations

- The `add` command is now deprecated. It will be removed in a future
  release.

## Fixes

- Remove hard-coded dependencies on POSIX path conventions in SSH support
  code ([#3400][])

- Emit an `add` result when adding a new subdataset during [save][] ([#3398][])

- SSH file transfer now actually opens a shared connection, if none exists
  yet ([#3403][])

## Enhancements and new features

- `SSHConnection` now offers methods for file upload and dowload (`get()`,
  `put()`. The previous `copy()` method only supported upload and was
  discontinued ([#3401][])


# 0.12.0rc3 (May 07, 2019) -- the revolution continues
￼
Continues API consolidation and replaces the `create` and `diff` command
with more performant implementations.

## Major refactoring and deprecations

- The previous `diff` command has been replaced by the diff variant
  from the [datalad-revolution][] extension.  ([#3366][])

- `rev-create` has been renamed to `create`, and the previous `create`
  has been removed.  ([#3383][])

- The procedure `setup_yoda_dataset` has been renamed to `cfg_yoda`
  ([#3353][]).

- The `--nosave` of `addurls` now affects only added content, not
  newly created subdatasets ([#3259][]).

- `Dataset.get_subdatasets` (deprecated since v0.9.0) has been
  removed.  ([#3336][])

- The `.is_dirty` method of `GitRepo` and `AnnexRepo` has been
  replaced by `.status` or, for a subset of cases, the `.dirty`
  property.  ([#3330][])

- `AnnexRepo.get_status` has been replaced by `AnnexRepo.status`.
  ([#3330][])

## Fixes

- [status][]
  - reported on directories that contained only ignored files ([#3238][])
  - gave a confusing failure when called from a subdataset with an
    explicitly specified dataset argument and "." as a path ([#3325][])
  - misleadingly claimed that the locally present content size was
    zero when `--annex basic` was specified ([#3378][])

- An informative error wasn't given when a download provider was
  invalid.  ([#3258][])

- Calling `rev-save PATH` saved unspecified untracked subdatasets.
  ([#3288][])

- The available choices for command-line options that take values are
  now displayed more consistently in the help output.  ([#3326][])

- The new pathlib-based code had various encoding issues on Python 2.
  ([#3332][])

## Enhancements and new features

- [wtf][] now includes information about the Python version.  ([#3255][])

- When operating in an annex repository, checking whether git-annex is
  available is now delayed until a call to git-annex is actually
  needed, allowing systems without git-annex to operate on annex
  repositories in a restricted fashion.  ([#3274][])

- The `load_stream` on helper now supports auto-detection of
  compressed files.  ([#3289][])

- `create` (formerly `rev-create`)
  - learned to be speedier by passing a path to `status` ([#3294][])
  - gained a `--cfg-proc` (or `-c`) convenience option for running
    configuration procedures (or more accurately any procedure that
    begins with "cfg_") in the newly created dataset ([#3353][])

- `AnnexRepo.set_metadata` now returns a list while
  `AnnexRepo.set_metadata_` returns a generator, a behavior which is
  consistent with the `add` and `add_` method pair.  ([#3298][])

- `AnnexRepo.get_metadata` now supports batch querying of known annex
   files.  Note, however, that callers should carefully validate the
   input paths because the batch call will silently hang if given
   non-annex files.  ([#3364][])

- [status][]
  - now reports a "bytesize" field for files tracked by Git ([#3299][])
  - gained a new option `eval_subdataset_state` that controls how the
    subdataset state is evaluated.  Depending on the information you
    need, you can select a less expensive mode to make `status`
    faster.  ([#3324][])
  - colors deleted files "red" ([#3334][])

- Querying repository content is faster due to batching of `git
  cat-file` calls.  ([#3301][])

- The dataset ID of a subdataset is now recorded in the superdataset.
  ([#3304][])

- `GitRepo.diffstatus`
  - now avoids subdataset recursion when the comparison is not with
    the working tree, which substantially improves performance when
    diffing large dataset hierarchies  ([#3314][])
  - got smarter and faster about labeling a subdataset as "modified"
    ([#3343][])

- `GitRepo.get_content_info` now supports disabling the file type
  evaluation, which gives a performance boost in cases where this
  information isn't needed.  ([#3362][])

- The XMP metadata extractor now filters based on file name to improve
  its performance.  ([#3329][])

# 0.12.0rc2 (Mar 18, 2019) -- revolution!

## Fixes

- `GitRepo.dirty` does not report on nested empty directories ([#3196][]).

- `GitRepo.save()` reports results on deleted files.

## Enhancements and new features

- Absorb a new set of core commands from the datalad-revolution extension:
  - `rev-status`: like `git status`, but simpler and working with dataset
     hierarchies
  - `rev-save`: a 2-in-1 replacement for save and add
  - `rev-create`: a ~30% faster create

- JSON support tools can now read and write compressed files.


# 0.12.0rc1 (Mar 03, 2019) -- to boldly go ...

## Major refactoring and deprecations

- Discontinued support for git-annex direct-mode (also no longer
  supported upstream).

## Enhancements and new features

- Dataset and Repo object instances are now hashable, and can be
  created based on pathlib Path object instances

- Imported various additional methods for the Repo classes to query
  information and save changes.


# 0.11.8 (Oct 11, 2019) -- annex-we-are-catching-up

## Fixes

- Our internal command runner failed to capture output in some cases.
  ([#3656][])
- Workaround in the tests around python in cPython >= 3.7.5 ';' in
  the filename confusing mimetypes ([#3769][]) ([#3770][])

## Enhancements and new features

- Prepared for upstream changes in git-annex, including support for
  the latest git-annex
  - 7.20190912 auto-upgrades v5 repositories to v7.  ([#3648][]) ([#3682][])
  - 7.20191009 fixed treatment of (larger/smaller)than in .gitattributes ([#3765][])

- The `cfg_text2git` procedure, as well the `--text-no-annex` option
  of [create][], now configure .gitattributes so that empty files are
  stored in git rather than annex.  ([#3667][])


# 0.11.7 (Sep 06, 2019) -- python2-we-still-love-you-but-...

Primarily bugfixes with some optimizations and refactorings.

## Fixes

- [addurls][]
  - now provides better handling when the URL file isn't in the
    expected format.  ([#3579][])
  - always considered a relative file for the URL file argument as
    relative to the current working directory, which goes against the
    convention used by other commands of taking relative paths as
    relative to the dataset argument.  ([#3582][])

- [run-procedure][]
  - hard coded "python" when formatting the command for non-executable
    procedures ending with ".py".  `sys.executable` is now used.
    ([#3624][])
  - failed if arguments needed more complicated quoting than simply
    surrounding the value with double quotes.  This has been resolved
    for systems that support `shlex.quote`, but note that on Windows
    values are left unquoted. ([#3626][])

- [siblings][] now displays an informative error message if a local
  path is given to `--url` but `--name` isn't specified.  ([#3555][])

- [sshrun][], the command DataLad uses for `GIT_SSH_COMMAND`, didn't
  support all the parameters that Git expects it to.  ([#3616][])

- Fixed a number of Unicode py2-compatibility issues. ([#3597][])

- [download-url][] now will create leading directories of the output path
  if they do not exist ([#3646][])

## Enhancements and new features

- The [annotate-paths][] helper now caches subdatasets it has seen to
  avoid unnecessary calls.  ([#3570][])

- A repeated configuration query has been dropped from the handling of
  `--proc-pre` and `--proc-post`.  ([#3576][])

- Calls to `git annex find` now use `--in=.` instead of the alias
  `--in=here` to take advantage of an optimization that git-annex (as
  of the current release, 7.20190730) applies only to the
  former. ([#3574][])

- [addurls][] now suggests close matches when the URL or file format
  contains an unknown field.  ([#3594][])

- Shared logic used in the setup.py files of Datalad and its
  extensions has been moved to modules in the _datalad_build_support/
  directory.  ([#3600][])

- Get ready for upcoming git-annex dropping support for direct mode
  ([#3631][])


# 0.11.6 (Jul 30, 2019) -- am I the last of 0.11.x?

Primarily bug fixes to achieve more robust performance

## Fixes

- Our tests needed various adjustments to keep up with upstream
  changes in Travis and Git. ([#3479][]) ([#3492][]) ([#3493][])

- `AnnexRepo.is_special_annex_remote` was too selective in what it
  considered to be a special remote.  ([#3499][])

- We now provide information about unexpected output when git-annex is
  called with `--json`.  ([#3516][])

- Exception logging in the `__del__` method of `GitRepo` and
  `AnnexRepo` no longer fails if the names it needs are no longer
  bound.  ([#3527][])

- [addurls][] botched the construction of subdataset paths that were
  more than two levels deep and failed to create datasets in a
  reliable, breadth-first order.  ([#3561][])

- Cloning a `type=git` special remote showed a spurious warning about
  the remote not being enabled.  ([#3547][])

## Enhancements and new features

- For calls to git and git-annex, we disable automatic garbage
  collection due to past issues with GitPython's state becoming stale,
  but doing so results in a larger .git/objects/ directory that isn't
  cleaned up until garbage collection is triggered outside of DataLad.
  Tests with the latest GitPython didn't reveal any state issues, so
  we've re-enabled automatic garbage collection.  ([#3458][])

- [rerun][] learned an `--explicit` flag, which it relays to its calls
  to [run][[]].  This makes it possible to call `rerun` in a dirty
  working tree ([#3498][]).

- The [metadata][] command aborts earlier if a metadata extractor is
  unavailable.  ([#3525][])

# 0.11.5 (May 23, 2019) -- stability is not overrated

Should be faster and less buggy, with a few enhancements.

## Fixes

- [create-sibling][]  ([#3318][])
  - Siblings are no longer configured with a post-update hook unless a
    web interface is requested with `--ui`.
  - `git submodule update --init` is no longer called from the
    post-update hook.
  - If `--inherit` is given for a dataset without a superdataset, a
    warning is now given instead of raising an error.
- The internal command runner failed on Python 2 when its `env`
  argument had unicode values.  ([#3332][])
- The safeguard that prevents creating a dataset in a subdirectory
  that already contains tracked files for another repository failed on
  Git versions before 2.14.  For older Git versions, we now warn the
  caller that the safeguard is not active.  ([#3347][])
- A regression introduced in v0.11.1 prevented [save][] from committing
  changes under a subdirectory when the subdirectory was specified as
  a path argument.  ([#3106][])
- A workaround introduced in v0.11.1 made it possible for [save][] to
  do a partial commit with an annex file that has gone below the
  `annex.largefiles` threshold.  The logic of this workaround was
  faulty, leading to files being displayed as typechanged in the index
  following the commit.  ([#3365][])
- The resolve_path() helper confused paths that had a semicolon for
  SSH RIs.  ([#3425][])
- The detection of SSH RIs has been improved.  ([#3425][])

## Enhancements and new features

- The internal command runner was too aggressive in its decision to
  sleep.  ([#3322][])
- The "INFO" label in log messages now retains the default text color
  for the terminal rather than using white, which only worked well for
  terminals with dark backgrounds.  ([#3334][])
- A short flag `-R` is now available for the `--recursion-limit` flag,
  a flag shared by several subcommands.  ([#3340][])
- The authentication logic for [create-sibling-github][] has been
  revamped and now supports 2FA.  ([#3180][])
- New configuration option `datalad.ui.progressbar` can be used to
  configure the default backend for progress reporting ("none", for
  example, results in no progress bars being shown).  ([#3396][])
- A new progress backend, available by setting datalad.ui.progressbar
  to "log", replaces progress bars with a log message upon completion
  of an action.  ([#3396][])
- DataLad learned to consult the [NO_COLOR][] environment variable and
  the new `datalad.ui.color` configuration option when deciding to
  color output.  The default value, "auto", retains the current
  behavior of coloring output if attached to a TTY ([#3407][]).
- [clean][] now removes annex transfer directories, which is useful
  for cleaning up failed downloads. ([#3374][])
- [clone][] no longer refuses to clone into a local path that looks
  like a URL, making its behavior consistent with `git clone`.
  ([#3425][])
- [wtf][]
  - Learned to fall back to the `dist` package if `platform.dist`,
    which has been removed in the yet-to-be-release Python 3.8, does
    not exist.  ([#3439][])
  - Gained a `--section` option for limiting the output to specific
    sections and a `--decor` option, which currently knows how to
    format the output as GitHub's `<details>` section.  ([#3440][])

# 0.11.4 (Mar 18, 2019) -- get-ready

Largely a bug fix release with a few enhancements

## Important

- 0.11.x series will be the last one with support for direct mode of [git-annex][]
  which is used on crippled (no symlinks and no locking) filesystems.
  v7 repositories should be used instead.

## Fixes

- Extraction of .gz files is broken without p7zip installed.  We now
  abort with an informative error in this situation.  ([#3176][])

- Committing failed in some cases because we didn't ensure that the
  path passed to `git read-tree --index-output=...` resided on the
  same filesystem as the repository.  ([#3181][])

- Some pointless warnings during metadata aggregation have been
  eliminated.  ([#3186][])

- With Python 3 the LORIS token authenticator did not properly decode
  a response ([#3205][]).

- With Python 3 downloaders unnecessarily decoded the response when
  getting the status, leading to an encoding error.  ([#3210][])

- In some cases, our internal command Runner did not adjust the
  environment's `PWD` to match the current working directory specified
  with the `cwd` parameter.  ([#3215][])

- The specification of the pyliblzma dependency was broken.  ([#3220][])

- [search] displayed an uninformative blank log message in some
  cases.  ([#3222][])

- The logic for finding the location of the aggregate metadata DB
  anchored the search path incorrectly, leading to a spurious warning.
  ([#3241][])

- Some progress bars were still displayed when stdout and stderr were
  not attached to a tty.  ([#3281][])

- Check for stdin/out/err to not be closed before checking for `.isatty`.
  ([#3268][])

## Enhancements and new features

- Creating a new repository now aborts if any of the files in the
  directory are tracked by a repository in a parent directory.
  ([#3211][])

- [run] learned to replace the `{tmpdir}` placeholder in commands with
  a temporary directory.  ([#3223][])
 
- [duecredit][] support has been added for citing DataLad itself as
  well as datasets that an analysis uses.  ([#3184][])

- The `eval_results` interface helper unintentionally modified one of
  its arguments.  ([#3249][])

- A few DataLad constants have been added, changed, or renamed ([#3250][]):
  - `HANDLE_META_DIR` is now `DATALAD_DOTDIR`.  The old name should be
     considered deprecated.
  - `METADATA_DIR` now refers to `DATALAD_DOTDIR/metadata` rather than
    `DATALAD_DOTDIR/meta` (which is still available as
    `OLDMETADATA_DIR`).
  - The new `DATASET_METADATA_FILE` refers to `METADATA_DIR/dataset.json`.
  - The new `DATASET_CONFIG_FILE` refers to `DATALAD_DOTDIR/config`.
  - `METADATA_FILENAME` has been renamed to `OLDMETADATA_FILENAME`.

# 0.11.3 (Feb 19, 2019) -- read-me-gently

Just a few of important fixes and minor enhancements.

## Fixes

- The logic for setting the maximum command line length now works
  around Python 3.4 returning an unreasonably high value for
  `SC_ARG_MAX` on Debian systems. ([#3165][])

- DataLad commands that are conceptually "read-only", such as
  `datalad ls -L`, can fail when the caller lacks write permissions
  because git-annex tries merging remote git-annex branches to update
  information about availability. DataLad now disables
  `annex.merge-annex-branches` in some common "read-only" scenarios to
  avoid these failures. ([#3164][])

## Enhancements and new features

- Accessing an "unbound" dataset method now automatically imports the
  necessary module rather than requiring an explicit import from the
  Python caller. For example, calling `Dataset.add` no longer needs to
  be preceded by `from datalad.distribution.add import Add` or an
  import of `datalad.api`. ([#3156][])

- Configuring the new variable `datalad.ssh.identityfile` instructs
  DataLad to pass a value to the `-i` option of `ssh`. ([#3149][])
  ([#3168][])

# 0.11.2 (Feb 07, 2019) -- live-long-and-prosper

A variety of bugfixes and enhancements

## Major refactoring and deprecations

- All extracted metadata is now placed under git-annex by default.
  Previously files smaller than 20 kb were stored in git. ([#3109][])
- The function `datalad.cmd.get_runner` has been removed. ([#3104][])

## Fixes

- Improved handling of long commands:
  - The code that inspected `SC_ARG_MAX` didn't check that the
    reported value was a sensible, positive number. ([#3025][])
  - More commands that invoke `git` and `git-annex` with file
    arguments learned to split up the command calls when it is likely
    that the command would fail due to exceeding the maximum supported
    length. ([#3138][])
- The `setup_yoda_dataset` procedure created a malformed
  .gitattributes line. ([#3057][])
- [download-url][] unnecessarily tried to infer the dataset when
  `--no-save` was given. ([#3029][])
- [rerun][] aborted too late and with a confusing message when a ref
  specified via `--onto` didn't exist. ([#3019][])
- [run][]:
  - `run` didn't preserve the current directory prefix ("./") on
     inputs and outputs, which is problematic if the caller relies on
     this representation when formatting the command. ([#3037][])
  - Fixed a number of unicode py2-compatibility issues. ([#3035][]) ([#3046][])
  - To proceed with a failed command, the user was confusingly
    instructed to use `save` instead of `add` even though `run` uses
    `add` underneath. ([#3080][])
- Fixed a case where the helper class for checking external modules
  incorrectly reported a module as unknown. ([#3051][])
- [add-archive-content][] mishandled the archive path when the leading
  path contained a symlink. ([#3058][])
- Following denied access, the credential code failed to consider a
  scenario, leading to a type error rather than an appropriate error
  message. ([#3091][])
- Some tests failed when executed from a `git worktree` checkout of the
  source repository. ([#3129][])
- During metadata extraction, batched annex processes weren't properly
  terminated, leading to issues on Windows. ([#3137][])
- [add][] incorrectly handled an "invalid repository" exception when
  trying to add a submodule. ([#3141][])
- Pass `GIT_SSH_VARIANT=ssh` to git processes to be able to specify
  alternative ports in SSH urls

## Enhancements and new features

- [search][] learned to suggest closely matching keys if there are no
  hits. ([#3089][])
- [create-sibling][]
  - gained a `--group` option so that the caller can specify the file
    system group for the repository. ([#3098][])
  - now understands SSH URLs that have a port in them (i.e. the
    "ssh://[user@]host.xz[:port]/path/to/repo.git/" syntax mentioned
    in `man git-fetch`). ([#3146][])
- Interface classes can now override the default renderer for
  summarizing results. ([#3061][])
- [run][]:
  - `--input` and `--output` can now be shortened to `-i` and `-o`.
    ([#3066][])
  - Placeholders such as "{inputs}" are now expanded in the command
    that is shown in the commit message subject. ([#3065][])
  - `interface.run.run_command` gained an `extra_inputs` argument so
    that wrappers like [datalad-container][] can specify additional inputs
    that aren't considered when formatting the command string. ([#3038][])
  - "--" can now be used to separate options for `run` and those for
    the command in ambiguous cases. ([#3119][])
- The utilities `create_tree` and `ok_file_has_content` now support
  ".gz" files. ([#3049][])
- The Singularity container for 0.11.1 now uses [nd_freeze][] to make
  its builds reproducible.
- A [publications][] page has been added to the documentation. ([#3099][])
- `GitRepo.set_gitattributes` now accepts a `mode` argument that
  controls whether the .gitattributes file is appended to (default) or
  overwritten. ([#3115][])
- `datalad --help` now avoids using `man` so that the list of
  subcommands is shown.  ([#3124][])

# 0.11.1 (Nov 26, 2018) -- v7-better-than-v6

Rushed out bugfix release to stay fully compatible with recent
[git-annex][] which introduced v7 to replace v6.

## Fixes

- [install][]: be able to install recursively into a dataset ([#2982][])
- [save][]: be able to commit/save changes whenever files potentially
  could have swapped their storage between git and annex
  ([#1651][]) ([#2752][]) ([#3009][])
- [aggregate-metadata][]:
  - dataset's itself is now not "aggregated" if specific paths are
    provided for aggregation ([#3002][]). That resolves the issue of
    `-r` invocation aggregating all subdatasets of the specified dataset
    as well
  - also compare/verify the actual content checksum of aggregated metadata
    while considering subdataset metadata for re-aggregation ([#3007][])
- `annex` commands are now chunked assuming 50% "safety margin" on the
  maximal command line length. Should resolve crashes while operating
  ot too many files at ones ([#3001][])
- `run` sidecar config processing ([#2991][])
- no double trailing period in docs ([#2984][])
- correct identification of the repository with symlinks in the paths
  in the tests ([#2972][])
- re-evaluation of dataset properties in case of dataset changes ([#2946][])
- [text2git][] procedure to use `ds.repo.set_gitattributes`
  ([#2974][]) ([#2954][])
- Switch to use plain `os.getcwd()` if inconsistency with env var
  `$PWD` is detected ([#2914][])
- Make sure that credential defined in env var takes precedence
  ([#2960][]) ([#2950][])

## Enhancements and new features

- [shub://datalad/datalad:git-annex-dev](https://singularity-hub.org/containers/5663/view)
  provides a Debian buster Singularity image with build environment for
  [git-annex][]. `tools/bisect-git-annex` provides a helper for running
  `git bisect` on git-annex using that Singularity container ([#2995][])
- Added `.zenodo.json` for better integration with Zenodo for citation
- [run-procedure][] now provides names and help messages with a custom
  renderer for ([#2993][])
- Documentation: point to [datalad-revolution][] extension (prototype of
  the greater DataLad future)
- [run][]
  - support injecting of a detached command ([#2937][])
- `annex` metadata extractor now extracts `annex.key` metadata record.
  Should allow now to identify uses of specific files etc ([#2952][])
- Test that we can install from http://datasets.datalad.org
- Proper rendering of `CommandError` (e.g. in case of "out of space"
  error) ([#2958][])


# 0.11.0 (Oct 23, 2018) -- Soon-to-be-perfect

[git-annex][] 6.20180913 (or later) is now required - provides a number of
fixes for v6 mode operations etc.

## Major refactoring and deprecations

- `datalad.consts.LOCAL_CENTRAL_PATH` constant was deprecated in favor
  of `datalad.locations.default-dataset` [configuration][config] variable
  ([#2835][])

## Minor refactoring

- `"notneeded"` messages are no longer reported by default results
  renderer
- [run][] no longer shows commit instructions upon command failure when
  `explicit` is true and no outputs are specified ([#2922][])
- `get_git_dir` moved into GitRepo ([#2886][])
- `_gitpy_custom_call` removed from GitRepo ([#2894][])
- `GitRepo.get_merge_base` argument is now called `commitishes` instead
  of `treeishes` ([#2903][])

## Fixes

- [update][] should not leave the dataset in non-clean state ([#2858][])
  and some other enhancements ([#2859][])
- Fixed chunking of the long command lines to account for decorators
  and other arguments ([#2864][])
- Progress bar should not crash the process on some missing progress
  information ([#2891][])
- Default value for `jobs` set to be `"auto"` (not `None`) to take
  advantage of possible parallel get if in `-g` mode ([#2861][])
- [wtf][] must not crash if `git-annex` is not installed etc ([#2865][]),
  ([#2865][]), ([#2918][]), ([#2917][])
- Fixed paths (with spaces etc) handling while reporting annex error
  output ([#2892][]), ([#2893][])
- `__del__` should not access `.repo` but `._repo` to avoid attempts
  for reinstantiation etc ([#2901][])
- Fix up submodule `.git` right in `GitRepo.add_submodule` to avoid
  added submodules being non git-annex friendly ([#2909][]), ([#2904][])
- [run-procedure][] ([#2905][])
  - now will provide dataset into the procedure if called within dataset
  - will not crash if procedure is an executable without `.py` or `.sh`
    suffixes
- Use centralized `.gitattributes` handling while setting annex backend
  ([#2912][])
- `GlobbedPaths.expand(..., full=True)` incorrectly returned relative
   paths when called more than once ([#2921][])

## Enhancements and new features

- Report progress on [clone][] when installing from "smart" git servers
  ([#2876][])
- Stale/unused `sth_like_file_has_content` was removed ([#2860][])
- Enhancements to [search][] to operate on "improved" metadata layouts
  ([#2878][])
- Output of `git annex init` operation is now logged ([#2881][])
- New
  - `GitRepo.cherry_pick` ([#2900][])
  - `GitRepo.format_commit` ([#2902][])
- [run-procedure][] ([#2905][])
  - procedures can now recursively be discovered in subdatasets as well.
    The uppermost has highest priority
  - Procedures in user and system locations now take precedence over
    those in datasets.

# 0.10.3.1 (Sep 13, 2018) -- Nothing-is-perfect

Emergency bugfix to address forgotten boost of version in
`datalad/version.py`.

# 0.10.3 (Sep 13, 2018) -- Almost-perfect

This is largely a bugfix release which addressed many (but not yet all)
issues of working with git-annex direct and version 6 modes, and operation
on Windows in general.  Among enhancements you will see the
support of public S3 buckets (even with periods in their names),
ability to configure new providers interactively, and improved `egrep`
search backend.

Although we do not require with this release, it is recommended to make
sure that you are using a recent `git-annex` since it also had a variety
of fixes and enhancements in the past months.

## Fixes

- Parsing of combined short options has been broken since DataLad
  v0.10.0. ([#2710][])
- The `datalad save` instructions shown by `datalad run` for a command
  with a non-zero exit were incorrectly formatted. ([#2692][])
- Decompression of zip files (e.g., through `datalad
  add-archive-content`) failed on Python 3.  ([#2702][])
- Windows:
  - colored log output was not being processed by colorama.  ([#2707][])
  - more codepaths now try multiple times when removing a file to deal
    with latency and locking issues on Windows.  ([#2795][])
- Internal git fetch calls have been updated to work around a
  GitPython `BadName` issue.  ([#2712][]), ([#2794][])
- The progess bar for annex file transferring was unable to handle an
  empty file.  ([#2717][])
- `datalad add-readme` halted when no aggregated metadata was found
  rather than displaying a warning.  ([#2731][])
- `datalad rerun` failed if `--onto` was specified and the history
  contained no run commits.  ([#2761][])
- Processing of a command's results failed on a result record with a
  missing value (e.g., absent field or subfield in metadata).  Now the
  missing value is rendered as "N/A".  ([#2725][]).
- A couple of documentation links in the "Delineation from related
  solutions" were misformatted.  ([#2773][])
- With the latest git-annex, several known V6 failures are no longer
  an issue.  ([#2777][])
- In direct mode, commit changes would often commit annexed content as
  regular Git files.  A new approach fixes this and resolves a good
  number of known failures.  ([#2770][])
- The reporting of command results failed if the current working
  directory was removed (e.g., after an unsuccessful `install`). ([#2788][])
- When installing into an existing empty directory, `datalad install`
  removed the directory after a failed clone.  ([#2788][])
- `datalad run` incorrectly handled inputs and outputs for paths with
  spaces and other characters that require shell escaping.  ([#2798][])
- Globbing inputs and outputs for `datalad run` didn't work correctly
  if a subdataset wasn't installed.  ([#2796][])
- Minor (in)compatibility with git 2.19 - (no) trailing period
  in an error message now. ([#2815][])

## Enhancements and new features

- Anonymous access is now supported for S3 and other downloaders.  ([#2708][])
- A new interface is available to ease setting up new providers.  ([#2708][])
- Metadata: changes to egrep mode search  ([#2735][])
  - Queries in egrep mode are now case-sensitive when the query
    contains any uppercase letters and are case-insensitive otherwise.
    The new mode egrepcs can be used to perform a case-sensitive query
    with all lower-case letters.
  - Search can now be limited to a specific key.
  - Multiple queries (list of expressions) are evaluated using AND to
    determine whether something is a hit.
  - A single multi-field query (e.g., `pa*:findme`) is a hit, when any
    matching field matches the query.
  - All matching key/value combinations across all (multi-field)
    queries are reported in the query_matched result field.
  - egrep mode now shows all hits rather than limiting the results to
    the top 20 hits.
- The documentation on how to format commands for `datalad run` has
  been improved.  ([#2703][])
- The method for determining the current working directory on Windows
  has been improved.  ([#2707][])
- `datalad --version` now simply shows the version without the
  license.  ([#2733][])
- `datalad export-archive` learned to export under an existing
  directory via its `--filename` option.  ([#2723][])
- `datalad export-to-figshare` now generates the zip archive in the
  root of the dataset unless `--filename` is specified.  ([#2723][])
- After importing `datalad.api`, `help(datalad.api)` (or
  `datalad.api?` in IPython) now shows a summary of the available
  DataLad commands.  ([#2728][])
- Support for using `datalad` from IPython has been improved.  ([#2722][])
- `datalad wtf` now returns structured data and reports the version of
  each extension.  ([#2741][])
- The internal handling of gitattributes information has been
  improved.  A user-visible consequence is that `datalad create
  --force` no longer duplicates existing attributes.  ([#2744][])
- The "annex" metadata extractor can now be used even when no content
  is present.  ([#2724][])
- The `add_url_to_file` method (called by commands like `datalad
  download-url` and `datalad add-archive-content`) learned how to
  display a progress bar.  ([#2738][])


# 0.10.2 (Jul 09, 2018) -- Thesecuriestever

Primarily a bugfix release to accommodate recent git-annex release
forbidding file:// and http://localhost/ URLs which might lead to
revealing private files if annex is publicly shared.

## Fixes

- fixed testing to be compatible with recent git-annex (6.20180626)
- [download-url][] will now download to current directory instead of the
  top of the dataset

## Enhancements and new features

- do not quote ~ in URLs to be consistent with quote implementation in
  Python 3.7 which now follows RFC 3986
- [run][] support for user-configured placeholder values
- documentation on native git-annex metadata support
- handle 401 errors from LORIS tokens
- `yoda` procedure will instantiate `README.md`
- `--discover` option added to [run-procedure][] to list available
  procedures

# 0.10.1 (Jun 17, 2018) -- OHBM polish

The is a minor bugfix release.

## Fixes

- Be able to use backports.lzma as a drop-in replacement for pyliblzma.
- Give help when not specifying a procedure name in `run-procedure`.
- Abort early when a downloader received no filename.
- Avoid `rerun` error when trying to unlock non-available files.

# 0.10.0 (Jun 09, 2018) -- The Release

This release is a major leap forward in metadata support.

## Major refactoring and deprecations

- Metadata
  - Prior metadata provided by datasets under `.datalad/meta` is no
    longer used or supported. Metadata must be reaggregated using 0.10
    version
  - Metadata extractor types are no longer auto-guessed and must be
    explicitly specified in `datalad.metadata.nativetype` config
    (could contain multiple values)
  - Metadata aggregation of a dataset hierarchy no longer updates all
    datasets in the tree with new metadata. Instead, only the target
    dataset is updated. This behavior can be changed via the --update-mode
    switch. The new default prevents needless modification of (3rd-party)
    subdatasets.
  - Neuroimaging metadata support has been moved into a dedicated extension:
    https://github.com/datalad/datalad-neuroimaging
- Crawler
  - moved into a dedicated extension:
    https://github.com/datalad/datalad-crawler
- `export_tarball` plugin has been generalized to `export_archive` and
  can now also generate ZIP archives.
- By default a dataset X is now only considered to be a super-dataset of
  another dataset Y, if Y is also a registered subdataset of X.

## Fixes

A number of fixes did not make it into the 0.9.x series:

- Dynamic configuration overrides via the `-c` option were not in effect.
- `save` is now more robust with respect to invocation in subdirectories
  of a dataset.
- `unlock` now reports correct paths when running in a dataset subdirectory.
- `get` is more robust to path that contain symbolic links.
- symlinks to subdatasets of a dataset are now correctly treated as a symlink,
  and not as a subdataset
- `add` now correctly saves staged subdataset additions.
- Running `datalad save` in a dataset no longer adds untracked content to the
  dataset. In order to add content a path has to be given, e.g. `datalad save .`
- `wtf` now works reliably with a DataLad that wasn't installed from Git (but,
  e.g., via pip)
- More robust URL handling in `simple_with_archives` crawler pipeline.

## Enhancements and new features

- Support for DataLad extension that can contribute API components from 3rd-party sources,
  incl. commands, metadata extractors, and test case implementations.
  See https://github.com/datalad/datalad-extension-template for a demo extension.
- Metadata (everything has changed!)
  - Metadata extraction and aggregation is now supported for datasets and individual
    files.
  - Metadata query via `search` can now discover individual files.
  - Extracted metadata can now be stored in XZ compressed files, is optionally
    annexed (when exceeding a configurable size threshold), and obtained on
    demand (new configuration option `datalad.metadata.create-aggregate-annex-limit`).
  - Status and availability of aggregated metadata can now be reported via
    `metadata --get-aggregates`
  - New configuration option `datalad.metadata.maxfieldsize` to exclude too large
    metadata fields from aggregation.
  - The type of metadata is no longer guessed during metadata extraction. A new
    configuration option `datalad.metadata.nativetype` was introduced to enable
    one or more particular metadata extractors for a dataset.
  - New configuration option `datalad.metadata.store-aggregate-content` to enable
    the storage of aggregated metadata for dataset content (i.e. file-based metadata)
    in contrast to just metadata describing a dataset as a whole.
- `search` was completely reimplemented. It offers three different modes now:
  - 'egrep' (default): expression matching in a plain string version of metadata
  - 'textblob': search a text version of all metadata using a fully featured
     query language (fast indexing, good for keyword search)
  - 'autofield': search an auto-generated index that preserves individual fields
     of metadata that can be represented in a tabular structure (substantial
     indexing cost, enables the most detailed queries of all modes)
- New extensions:
  - [addurls][], an extension for creating a dataset (and possibly subdatasets)
    from a list of URLs.
  - export_to_figshare
  - extract_metadata
- add_readme makes use of available metadata
- By default the wtf extension now hides sensitive information, which can be
  included in the output by passing `--senstive=some` or `--senstive=all`.
- Reduced startup latency by only importing commands necessary for a particular
  command line call.
- [create][]:
  - `-d <parent> --nosave` now registers subdatasets, when possible.
  - `--fake-dates` configures dataset to use fake-dates
- [run][] now provides a way for the caller to save the result when a
  command has a non-zero exit status.
- `datalad rerun` now has a `--script` option that can be used to extract
  previous commands into a file.
- A DataLad Singularity container is now available on
  [Singularity Hub](https://singularity-hub.org/collections/667).
- More casts have been embedded in the [use case section of the documentation](http://docs.datalad.org/en/docs/usecases/index.html).
- `datalad --report-status` has a new value 'all' that can be used to
  temporarily re-enable reporting that was disable by configuration settings.


# 0.9.3 (Mar 16, 2018) -- pi+0.02 release

Some important bug fixes which should improve usability

## Fixes

- `datalad-archives` special remote now will lock on acquiring or
  extracting an archive - this allows for it to be used with -J flag
  for parallel operation
- relax introduced in 0.9.2 demand on git being configured for datalad
  operation - now we will just issue a warning
- `datalad ls` should now list "authored date" and work also for datasets
  in detached HEAD mode
- `datalad save` will now save original file as well, if file was
  "git mv"ed, so you can now `datalad run git mv old new` and have
  changes recorded

## Enhancements and new features

- `--jobs` argument now could take `auto` value which would decide on
  # of jobs depending on the # of available CPUs.
  `git-annex` > 6.20180314 is recommended to avoid regression with -J.
- memoize calls to `RI` meta-constructor -- should speed up operation a
  bit
- `DATALAD_SEED` environment variable could be used to seed Python RNG
  and provide reproducible UUIDs etc (useful for testing and demos)


# 0.9.2 (Mar 04, 2018) -- it is (again) better than ever

Largely a bugfix release with a few enhancements.

## Fixes

- Execution of external commands (git) should not get stuck when
  lots of both stdout and stderr output, and should not loose remaining
  output in some cases
- Config overrides provided in the command line (-c) should now be
  handled correctly
- Consider more remotes (not just tracking one, which might be none)
  while installing subdatasets
- Compatibility with git 2.16 with some changed behaviors/annotations
  for submodules
- Fail `remove` if `annex drop` failed
- Do not fail operating on files which start with dash (-)
- URL unquote paths within S3, URLs and DataLad RIs (///)
- In non-interactive mode fail if authentication/access fails
- Web UI:
  - refactored a little to fix incorrect listing of submodules in
    subdirectories
  - now auto-focuses on search edit box upon entering the page
- Assure that extracted from tarballs directories have executable bit set

## Enhancements and new features

- A log message and progress bar will now inform if a tarball to be
  downloaded while getting specific files
  (requires git-annex > 6.20180206)
- A dedicated `datalad rerun` command capable of rerunning entire
  sequences of previously `run` commands.
  **Reproducibility through VCS. Use `run` even if not interested in `rerun`**
- Alert the user if `git` is not yet configured but git operations
  are requested
- Delay collection of previous ssh connections until it is actually
  needed.  Also do not require ':' while specifying ssh host
- AutomagicIO: Added proxying of isfile, lzma.LZMAFile and io.open
- Testing:
  - added DATALAD_DATASETS_TOPURL=http://datasets-tests.datalad.org to
    run tests against another website to not obscure access stats
  - tests run against temporary HOME to avoid side-effects
  - better unit-testing of interactions with special remotes
- CONTRIBUTING.md describes how to setup and use `git-hub` tool to
  "attach" commits to an issue making it into a PR
- DATALAD_USE_DEFAULT_GIT env variable could be used to cause DataLad
  to use default (not the one possibly bundled with git-annex) git
- Be more robust while handling not supported requests by annex in
  special remotes
- Use of `swallow_logs` in the code was refactored away -- less
  mysteries now, just increase logging level
- `wtf` plugin will report more information about environment, externals
  and the system


# 0.9.1 (Oct 01, 2017) -- "DATALAD!"(JBTM)

Minor bugfix release

## Fixes

- Should work correctly with subdatasets named as numbers of bool
  values (requires also GitPython >= 2.1.6)
- Custom special remotes should work without crashing with 
  git-annex >= 6.20170924


# 0.9.0 (Sep 19, 2017) -- isn't it a lucky day even though not a Friday?

## Major refactoring and deprecations

- the `files` argument of [save][] has been renamed to `path` to be uniform with
  any other command
- all major commands now implement more uniform API semantics and result reporting.
  Functionality for modification detection of dataset content has been completely replaced
  with a more efficient implementation
- [publish][] now features a `--transfer-data` switch that allows for a
  disambiguous specification of whether to publish data -- independent of
  the selection which datasets to publish (which is done via their paths).
  Moreover, [publish][] now transfers data before repository content is pushed.

## Fixes

- [drop][] no longer errors when some subdatasets are not installed
- [install][] will no longer report nothing when a Dataset instance was
  given as a source argument, but rather perform as expected
- [remove][] doesn't remove when some files of a dataset could not be dropped
- [publish][] 
  - no longer hides error during a repository push
  - publish behaves "correctly" for `--since=` in considering only the
    differences the last "pushed" state
  - data transfer handling while publishing with dependencies, to github
- improved robustness with broken Git configuration
- [search][] should search for unicode strings correctly and not crash
- robustify git-annex special remotes protocol handling to allow for spaces in
  the last argument
- UI credentials interface should now allow to Ctrl-C the entry
- should not fail while operating on submodules named with
  numerics only or by bool (true/false) names
- crawl templates should not now override settings for `largefiles` if 
  specified in `.gitattributes`


## Enhancements and new features

- **Exciting new feature** [run][] command to protocol execution of an external 
  command and rerun computation if desired. 
  See [screencast](http://datalad.org/features.html#reproducible-science)
- [save][] now uses Git for detecting with sundatasets need to be inspected for
  potential changes, instead of performing a complete traversal of a dataset tree
- [add][] looks for changes relative to the last commited state of a dataset
  to discover files to add more efficiently
- [diff][] can now report untracked files in addition to modified files
- [uninstall][] will check itself whether a subdataset is properly registered in a
  superdataset, even when no superdataset is given in a call
- [subdatasets][] can now configure subdatasets for exclusion from recursive
  installation (`datalad-recursiveinstall` submodule configuration property)
- precrafted pipelines of [crawl][] now will not override `annex.largefiles`
  setting if any was set within `.gitattribues` (e.g. by `datalad create --text-no-annex`)
- framework for screencasts: `tools/cast*` tools and sample cast scripts under
  `doc/casts` which are published at [datalad.org/features.html](http://datalad.org/features.html)
- new [project YouTube channel](https://www.youtube.com/channel/UCB8-Zf7D0DSzAsREoIt0Bvw) 
- tests failing in direct and/or v6 modes marked explicitly

# 0.8.1 (Aug 13, 2017) -- the best birthday gift

Bugfixes

## Fixes

- Do not attempt to [update][] a not installed sub-dataset
- In case of too many files to be specified for [get][] or [copy_to][], we
  will make multiple invocations of underlying git-annex command to not
  overfill command line
- More robust handling of unicode output in terminals which might not support it

## Enhancements and new features

- Ship a copy of numpy.testing to facilitate [test][] without requiring numpy
  as dependency. Also allow to pass to command which test(s) to run
- In [get][] and [copy_to][] provide actual original requested paths, not the
  ones we deduced need to be transferred, solely for knowing the total


# 0.8.0 (Jul 31, 2017) -- it is better than ever

A variety of fixes and enhancements

## Fixes

- [publish][] would now push merged `git-annex` branch even if no other changes
  were done
- [publish][] should be able to publish using relative path within SSH URI
  (git hook would use relative paths)
- [publish][] should better tollerate publishing to pure git and `git-annex` 
  special remotes 

## Enhancements and new features

- [plugin][] mechanism came to replace [export][]. See [export_tarball][] for the
  replacement of [export][].  Now it should be easy to extend datalad's interface
  with custom functionality to be invoked along with other commands.
- Minimalistic coloring of the results rendering
- [publish][]/`copy_to` got progress bar report now and support of `--jobs`
- minor fixes and enhancements to crawler (e.g. support of recursive removes)


# 0.7.0 (Jun 25, 2017) -- when it works - it is quite awesome!

New features, refactorings, and bug fixes.

## Major refactoring and deprecations

- [add-sibling][] has been fully replaced by the [siblings][] command
- [create-sibling][], and [unlock][] have been re-written to support the
  same common API as most other commands

## Enhancements and new features

- [siblings][] can now be used to query and configure a local repository by
  using the sibling name ``here``
- [siblings][] can now query and set annex preferred content configuration. This
  includes ``wanted`` (as previously supported in other commands), and now
  also ``required``
- New [metadata][] command to interface with datasets/files [meta-data][] 
- Documentation for all commands is now built in a uniform fashion
- Significant parts of the documentation of been updated
- Instantiate GitPython's Repo instances lazily

## Fixes

- API documentation is now rendered properly as HTML, and is easier to browse by
  having more compact pages
- Closed files left open on various occasions (Popen PIPEs, etc)
- Restored basic (consumer mode of operation) compatibility with Windows OS 


# 0.6.0 (Jun 14, 2017) -- German perfectionism

This release includes a **huge** refactoring to make code base and functionality
more robust and flexible

- outputs from API commands could now be highly customized.  See
  `--output-format`, `--report-status`, `--report-type`, and `--report-type`
  options for [datalad][] command.
- effort was made to refactor code base so that underlying functions behave as
  generators where possible
- input paths/arguments analysis was redone for majority of the commands to provide
  unified behavior

## Major refactoring and deprecations

- `add-sibling` and `rewrite-urls` were refactored in favor of new [siblings][]
  command which should be used for siblings manipulations
- 'datalad.api.alwaysrender' config setting/support is removed in favor of new
  outputs processing

## Fixes

- Do not flush manually git index in pre-commit to avoid "Death by the Lock" issue
- Deployed by [publish][] `post-update` hook script now should be more robust
  (tolerate directory names with spaces, etc.)
- A variety of fixes, see
  [list of pull requests and issues closed](https://github.com/datalad/datalad/milestone/41?closed=1)
  for more information

## Enhancements and new features

- new [annotate-paths][] plumbing command to inspect and annotate provided
  paths.  Use `--modified` to summarize changes between different points in
  the history
- new [clone][] plumbing command to provide a subset (install a single dataset
  from a URL) functionality of [install][]
- new [diff][] plumbing command
- new [siblings][] command to list or manipulate siblings
- new [subdatasets][] command to list subdatasets and their properties
- [drop][] and [remove][] commands were refactored
- `benchmarks/` collection of [Airspeed velocity](https://github.com/spacetelescope/asv/)
  benchmarks initiated.  See reports at http://datalad.github.io/datalad/
- crawler would try to download a new url multiple times increasing delay between
  attempts.  Helps to resolve problems with extended crawls of Amazon S3
- [CRCNS][] crawler pipeline now also fetches and aggregates meta-data for the
  datasets from datacite
- overall optimisations to benefit from the aforementioned refactoring and
  improve user-experience
- a few stub and not (yet) implemented commands (e.g. `move`) were removed from
  the interface
- Web frontend got proper coloring for the breadcrumbs and some additional
  caching to speed up interactions.  See http://datasets.datalad.org
- Small improvements to the online documentation.  See e.g.
  [summary of differences between git/git-annex/datalad](http://docs.datalad.org/en/latest/related.html#git-git-annex-datalad)

# 0.5.1 (Mar 25, 2017) -- cannot stop the progress

A bugfix release

## Fixes

- [add][] was forcing addition of files to annex regardless of settings
  in `.gitattributes`.  Now that decision is left to annex by default
- `tools/testing/run_doc_examples` used to run
  doc examples as tests, fixed up to provide status per each example
  and not fail at once
- `doc/examples`
  - [3rdparty_analysis_workflow.sh](http://docs.datalad.org/en/latest/generated/examples/3rdparty_analysis_workflow.html)
    was fixed up to reflect changes in the API of 0.5.0.
- progress bars
  - should no longer crash **datalad** and report correct sizes and speeds
  - should provide progress reports while using Python 3.x

## Enhancements and new features

- `doc/examples`
  - [nipype_workshop_dataset.sh](http://docs.datalad.org/en/latest/generated/examples/nipype_workshop_dataset.html)
    new example to demonstrate how new super- and sub- datasets were established
    as a part of our datasets collection


# 0.5.0 (Mar 20, 2017) -- it's huge

This release includes an avalanche of bug fixes, enhancements, and
additions which at large should stay consistent with previous behavior
but provide better functioning.  Lots of code was refactored to provide
more consistent code-base, and some API breakage has happened.  Further
work is ongoing to standardize output and results reporting
([#1350][])

## Most notable changes

- requires [git-annex][] >= 6.20161210 (or better even >= 6.20161210 for
  improved functionality)
- commands should now operate on paths specified (if any), without
  causing side-effects on other dirty/staged files
- [save][]
    - `-a` is deprecated in favor of `-u` or `--all-updates`
      so only changes known components get saved, and no new files
      automagically added
    - `-S` does no longer store the originating dataset in its commit
       message
- [add][]
    - can specify commit/save message with `-m`
- [add-sibling][] and [create-sibling][]
    - now take the name of the sibling (remote) as a `-s` (`--name`)
      option, not a positional argument
    - `--publish-depends` to setup publishing data and code to multiple
      repositories (e.g. github + webserve) should now be functional
      see [this comment](https://github.com/datalad/datalad/issues/335#issuecomment-277240733)
    - got `--publish-by-default` to specify what refs should be published
      by default
    - got `--annex-wanted`, `--annex-groupwanted` and `--annex-group`
      settings which would be used to instruct annex about preferred
      content. [publish][] then will publish data using those settings if
      `wanted` is set.
    - got `--inherit` option to automagically figure out url/wanted and
      other git/annex settings for new remote sub-dataset to be constructed
- [publish][]
    - got `--skip-failing` refactored into `--missing` option
      which could use new feature of [create-sibling][] `--inherit`

## Fixes

- More consistent interaction through ssh - all ssh connections go
  through [sshrun][] shim for a "single point of authentication", etc.
- More robust [ls][] operation outside of the datasets
- A number of fixes for direct and v6 mode of annex

## Enhancements and new features

- New [drop][] and [remove][] commands
- [clean][]
    - got `--what` to specify explicitly what cleaning steps to perform
      and now could be invoked with `-r`
- `datalad` and `git-annex-remote*` scripts now do not use setuptools
  entry points mechanism and rely on simple import to shorten start up time
- [Dataset][] is also now using [Flyweight pattern][], so the same instance is
  reused for the same dataset
- progressbars should not add more empty lines

## Internal refactoring

- Majority of the commands now go through `_prep` for arguments validation
  and pre-processing to avoid recursive invocations


# 0.4.1 (Nov 10, 2016) -- CA release

Requires now GitPython >= 2.1.0

## Fixes

- [save][]
     - to not save staged files if explicit paths were provided
- improved (but not yet complete) support for direct mode
- [update][] to not crash if some sub-datasets are not installed
- do not log calls to `git config` to avoid leakage of possibly 
  sensitive settings to the logs

## Enhancements and new features

- New [rfc822-compliant metadata][] format
- [save][]
    - -S to save the change also within all super-datasets
- [add][] now has progress-bar reporting
- [create-sibling-github][] to create a :term:`sibling` of a dataset on
  github
- [OpenfMRI][] crawler and datasets were enriched with URLs to separate
  files where also available from openfmri s3 bucket
  (if upgrading your datalad datasets, you might need to run
  `git annex enableremote datalad` to make them available)
- various enhancements to log messages
- web interface
    - populates "install" box first thus making UX better over slower
      connections


# 0.4 (Oct 22, 2016) -- Paris is waiting

Primarily it is a bugfix release but because of significant refactoring
of the [install][] and [get][] implementation, it gets a new minor release. 

## Fixes

- be able to [get][] or [install][] while providing paths while being 
  outside of a dataset
- remote annex datasets get properly initialized
- robust detection of outdated [git-annex][]

## Enhancements and new features

- interface changes
    - [get][] `--recursion-limit=existing` to not recurse into not-installed
       subdatasets
    - [get][] `-n` to possibly install sub-datasets without getting any data
    - [install][] `--jobs|-J` to specify number of parallel jobs for annex 
      [get][] call could use (ATM would not work when data comes from archives)
- more (unit-)testing
- documentation: see http://docs.datalad.org/en/latest/basics.html
  for basic principles and useful shortcuts in referring to datasets
- various webface improvements:  breadcrumb paths, instructions how
  to install dataset, show version from the tags, etc.

# 0.3.1 (Oct 1, 2016) -- what a wonderful week

Primarily bugfixes but also a number of enhancements and core
refactorings

## Fixes

- do not build manpages and examples during installation to avoid
  problems with possibly previously outdated dependencies
- [install][] can be called on already installed dataset (with `-r` or
  `-g`)

## Enhancements and new features

- complete overhaul of datalad configuration settings handling
  (see [Configuration documentation][]), so majority of the environment.
  Now uses git format and stores persistent configuration settings under
  `.datalad/config` and local within `.git/config`
  variables we have used were renamed to match configuration names
- [create-sibling][] does not now by default upload web front-end
- [export][] command with a plug-in interface and `tarball` plugin to export
  datasets
- in Python, `.api` functions with rendering of results in command line
  got a _-suffixed sibling, which would render results as well in Python
  as well (e.g., using `search_` instead of `search` would also render
  results, not only output them back as Python objects)
- [get][]
    - `--jobs` option (passed to `annex get`) for parallel downloads
    - total and per-download (with git-annex >= 6.20160923) progress bars
      (note that if content is to be obtained from an archive, no progress
      will be reported yet)
- [install][] `--reckless` mode option
- [search][]
    - highlights locations and fieldmaps for better readability
    - supports `-d^` or `-d///` to point to top-most or centrally
      installed meta-datasets
    - "complete" paths to the datasets are reported now
    - `-s` option to specify which fields (only) to search
- various enhancements and small fixes to [meta-data][] handling, [ls][],
  custom remotes, code-base formatting, downloaders, etc
- completely switched to `tqdm` library (`progressbar` is no longer
  used/supported)


# 0.3 (Sep 23, 2016) -- winter is coming

Lots of everything, including but not limited to

- enhanced index viewer, as the one on http://datasets.datalad.org
- initial new data providers support: [Kaggle][], [BALSA][], [NDA][], [NITRC][]
- initial [meta-data support and management][]
- new and/or improved crawler pipelines for [BALSA][], [CRCNS][], [OpenfMRI][]
- refactored [install][] command, now with separate [get][]
- some other commands renaming/refactoring (e.g., [create-sibling][])
- datalad [search][] would give you an option to install datalad's 
  super-dataset under ~/datalad if ran outside of a dataset

## 0.2.3 (Jun 28, 2016) -- busy OHBM

New features and bugfix release

- support of /// urls to point to http://datasets.datalad.org
- variety of fixes and enhancements throughout

## 0.2.2 (Jun 20, 2016) -- OHBM we are coming!

New feature and bugfix release

- greately improved documentation
- publish command API RFing allows for custom options to annex, and uses
  --to REMOTE for consistent with annex invocation
- variety of fixes and enhancements throughout

## 0.2.1 (Jun 10, 2016)

- variety of fixes and enhancements throughout

# 0.2 (May 20, 2016)

Major RFing to switch from relying on rdf to git native submodules etc

# 0.1 (Oct 14, 2015)

Release primarily focusing on interface functionality including initial
publishing

[git-annex]: http://git-annex.branchable.com/
[gx-sameas]: https://git-annex.branchable.com/tips/multiple_remotes_accessing_the_same_data_store/
[duecredit]: https://github.com/duecredit/duecredit

[Kaggle]: https://www.kaggle.com
[BALSA]: http://balsa.wustl.edu
[NDA]: http://data-archive.nimh.nih.gov
[NITRC]: https://www.nitrc.org
[CRCNS]: http://crcns.org
[FCON1000]: http://fcon_1000.projects.nitrc.org
[OpenfMRI]: http://openfmri.org

[Configuration documentation]: http://docs.datalad.org/config.html

[Dataset]: http://docs.datalad.org/en/latest/generated/datalad.api.Dataset.html
[Sibling]: http://docs.datalad.org/en/latest/glossary.html

[rfc822-compliant metadata]: http://docs.datalad.org/en/latest/metadata.html#rfc822-compliant-meta-data
[meta-data support and management]: http://docs.datalad.org/en/latest/cmdline.html#meta-data-handling
[meta-data]: http://docs.datalad.org/en/latest/cmdline.html#meta-data-handling

[add-archive-content]: https://datalad.readthedocs.io/en/latest/generated/man/datalad-add-archive-content.html
[add-sibling]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-add-sibling.html
[add]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html
[addurls]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html
[annotate-paths]: http://docs.datalad.org/en/latest/generated/man/datalad-annotate-paths.html
[clean]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-clean.html
[clone]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html
[config]: http://docs.datalad.org/en/latest/config.html
[configuration]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-configuration.html
[copy-file]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-copy-file.html
[copy_to]: http://docs.datalad.org/en/latest/_modules/datalad/support/annexrepo.html?highlight=%22copy_to%22
[create]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html
[create-sibling-github]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html
[create-sibling-ria]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-ria.html
[create-sibling]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html
[datalad]: http://docs.datalad.org/en/latest/generated/man/datalad.html
[datalad-container]: https://github.com/datalad/datalad-container
[datalad-revolution]: http://github.com/datalad/datalad-revolution
[download-url]: https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html
[diff]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html
[drop]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-drop.html
[export-archive-ora]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-export-archive-ora.html
[export]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-export.html
[export_tarball]: http://docs.datalad.org/en/latest/generated/datalad.plugin.export_tarball.html
[get]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html
[install]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html
[ls]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-ls.html
[metadata]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-metadata.html
[nd_freeze]: https://github.com/neurodebian/neurodebian/blob/master/tools/nd_freeze
[plugin]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-plugin.html
[publications]: https://datalad.readthedocs.io/en/latest/publications.html
[publish]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html
[push]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html
[remove]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-remove.html
[rerun]: https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html
[run]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html
[run-procedure]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html
[save]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html
[search]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html
[siblings]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html
[sshrun]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-sshrun.html
[status]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html
[subdatasets]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html
[unlock]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-unlock.html
[update]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html
[wtf]: http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html

[handbook]: http://handbook.datalad.org
[handbook-scalable-datastore]: http://handbook.datalad.org/en/latest/usecases/datastorage_for_institutions.html
[hooks]: http://handbook.datalad.org/en/latest/basics/101-145-hooks.html
[Flyweight pattern]: https://en.wikipedia.org/wiki/Flyweight_pattern
[NO_COLOR]: https://no-color.org/

[#5420]: https://github.com/datalad/datalad/issues/5420
[#5428]: https://github.com/datalad/datalad/issues/5428
[#5459]: https://github.com/datalad/datalad/issues/5459
[#5554]: https://github.com/datalad/datalad/issues/5554
[#5564]: https://github.com/datalad/datalad/issues/5564
[#5672]: https://github.com/datalad/datalad/issues/5672
[#1350]: https://github.com/datalad/datalad/issues/1350
[#1651]: https://github.com/datalad/datalad/issues/1651
[#2534]: https://github.com/datalad/datalad/issues/2534
[#2566]: https://github.com/datalad/datalad/issues/2566
[#2692]: https://github.com/datalad/datalad/issues/2692
[#2702]: https://github.com/datalad/datalad/issues/2702
[#2703]: https://github.com/datalad/datalad/issues/2703
[#2707]: https://github.com/datalad/datalad/issues/2707
[#2708]: https://github.com/datalad/datalad/issues/2708
[#2710]: https://github.com/datalad/datalad/issues/2710
[#2712]: https://github.com/datalad/datalad/issues/2712
[#2717]: https://github.com/datalad/datalad/issues/2717
[#2722]: https://github.com/datalad/datalad/issues/2722
[#2723]: https://github.com/datalad/datalad/issues/2723
[#2724]: https://github.com/datalad/datalad/issues/2724
[#2725]: https://github.com/datalad/datalad/issues/2725
[#2728]: https://github.com/datalad/datalad/issues/2728
[#2731]: https://github.com/datalad/datalad/issues/2731
[#2733]: https://github.com/datalad/datalad/issues/2733
[#2735]: https://github.com/datalad/datalad/issues/2735
[#2738]: https://github.com/datalad/datalad/issues/2738
[#2741]: https://github.com/datalad/datalad/issues/2741
[#2744]: https://github.com/datalad/datalad/issues/2744
[#2752]: https://github.com/datalad/datalad/issues/2752
[#2754]: https://github.com/datalad/datalad/issues/2754
[#2761]: https://github.com/datalad/datalad/issues/2761
[#2770]: https://github.com/datalad/datalad/issues/2770
[#2773]: https://github.com/datalad/datalad/issues/2773
[#2777]: https://github.com/datalad/datalad/issues/2777
[#2788]: https://github.com/datalad/datalad/issues/2788
[#2794]: https://github.com/datalad/datalad/issues/2794
[#2795]: https://github.com/datalad/datalad/issues/2795
[#2796]: https://github.com/datalad/datalad/issues/2796
[#2798]: https://github.com/datalad/datalad/issues/2798
[#2815]: https://github.com/datalad/datalad/issues/2815
[#2835]: https://github.com/datalad/datalad/issues/2835
[#2858]: https://github.com/datalad/datalad/issues/2858
[#2859]: https://github.com/datalad/datalad/issues/2859
[#2860]: https://github.com/datalad/datalad/issues/2860
[#2861]: https://github.com/datalad/datalad/issues/2861
[#2864]: https://github.com/datalad/datalad/issues/2864
[#2865]: https://github.com/datalad/datalad/issues/2865
[#2876]: https://github.com/datalad/datalad/issues/2876
[#2878]: https://github.com/datalad/datalad/issues/2878
[#2881]: https://github.com/datalad/datalad/issues/2881
[#2886]: https://github.com/datalad/datalad/issues/2886
[#2891]: https://github.com/datalad/datalad/issues/2891
[#2892]: https://github.com/datalad/datalad/issues/2892
[#2893]: https://github.com/datalad/datalad/issues/2893
[#2894]: https://github.com/datalad/datalad/issues/2894
[#2897]: https://github.com/datalad/datalad/issues/2897
[#2900]: https://github.com/datalad/datalad/issues/2900
[#2901]: https://github.com/datalad/datalad/issues/2901
[#2902]: https://github.com/datalad/datalad/issues/2902
[#2903]: https://github.com/datalad/datalad/issues/2903
[#2904]: https://github.com/datalad/datalad/issues/2904
[#2905]: https://github.com/datalad/datalad/issues/2905
[#2909]: https://github.com/datalad/datalad/issues/2909
[#2912]: https://github.com/datalad/datalad/issues/2912
[#2914]: https://github.com/datalad/datalad/issues/2914
[#2917]: https://github.com/datalad/datalad/issues/2917
[#2918]: https://github.com/datalad/datalad/issues/2918
[#2921]: https://github.com/datalad/datalad/issues/2921
[#2922]: https://github.com/datalad/datalad/issues/2922
[#2937]: https://github.com/datalad/datalad/issues/2937
[#2946]: https://github.com/datalad/datalad/issues/2946
[#2950]: https://github.com/datalad/datalad/issues/2950
[#2952]: https://github.com/datalad/datalad/issues/2952
[#2954]: https://github.com/datalad/datalad/issues/2954
[#2958]: https://github.com/datalad/datalad/issues/2958
[#2960]: https://github.com/datalad/datalad/issues/2960
[#2972]: https://github.com/datalad/datalad/issues/2972
[#2974]: https://github.com/datalad/datalad/issues/2974
[#2982]: https://github.com/datalad/datalad/issues/2982
[#2984]: https://github.com/datalad/datalad/issues/2984
[#2991]: https://github.com/datalad/datalad/issues/2991
[#2993]: https://github.com/datalad/datalad/issues/2993
[#2995]: https://github.com/datalad/datalad/issues/2995
[#3001]: https://github.com/datalad/datalad/issues/3001
[#3002]: https://github.com/datalad/datalad/issues/3002
[#3007]: https://github.com/datalad/datalad/issues/3007
[#3009]: https://github.com/datalad/datalad/issues/3009
[#3019]: https://github.com/datalad/datalad/issues/3019
[#3025]: https://github.com/datalad/datalad/issues/3025
[#3029]: https://github.com/datalad/datalad/issues/3029
[#3035]: https://github.com/datalad/datalad/issues/3035
[#3037]: https://github.com/datalad/datalad/issues/3037
[#3038]: https://github.com/datalad/datalad/issues/3038
[#3046]: https://github.com/datalad/datalad/issues/3046
[#3049]: https://github.com/datalad/datalad/issues/3049
[#3051]: https://github.com/datalad/datalad/issues/3051
[#3057]: https://github.com/datalad/datalad/issues/3057
[#3058]: https://github.com/datalad/datalad/issues/3058
[#3061]: https://github.com/datalad/datalad/issues/3061
[#3065]: https://github.com/datalad/datalad/issues/3065
[#3066]: https://github.com/datalad/datalad/issues/3066
[#3080]: https://github.com/datalad/datalad/issues/3080
[#3089]: https://github.com/datalad/datalad/issues/3089
[#3091]: https://github.com/datalad/datalad/issues/3091
[#3098]: https://github.com/datalad/datalad/issues/3098
[#3099]: https://github.com/datalad/datalad/issues/3099
[#3102]: https://github.com/datalad/datalad/issues/3102
[#3104]: https://github.com/datalad/datalad/issues/3104
[#3106]: https://github.com/datalad/datalad/issues/3106
[#3109]: https://github.com/datalad/datalad/issues/3109
[#3115]: https://github.com/datalad/datalad/issues/3115
[#3119]: https://github.com/datalad/datalad/issues/3119
[#3124]: https://github.com/datalad/datalad/issues/3124
[#3129]: https://github.com/datalad/datalad/issues/3129
[#3137]: https://github.com/datalad/datalad/issues/3137
[#3138]: https://github.com/datalad/datalad/issues/3138
[#3141]: https://github.com/datalad/datalad/issues/3141
[#3146]: https://github.com/datalad/datalad/issues/3146
[#3149]: https://github.com/datalad/datalad/issues/3149
[#3156]: https://github.com/datalad/datalad/issues/3156
[#3164]: https://github.com/datalad/datalad/issues/3164
[#3165]: https://github.com/datalad/datalad/issues/3165
[#3168]: https://github.com/datalad/datalad/issues/3168
[#3176]: https://github.com/datalad/datalad/issues/3176
[#3180]: https://github.com/datalad/datalad/issues/3180
[#3181]: https://github.com/datalad/datalad/issues/3181
[#3184]: https://github.com/datalad/datalad/issues/3184
[#3186]: https://github.com/datalad/datalad/issues/3186
[#3196]: https://github.com/datalad/datalad/issues/3196
[#3205]: https://github.com/datalad/datalad/issues/3205
[#3210]: https://github.com/datalad/datalad/issues/3210
[#3211]: https://github.com/datalad/datalad/issues/3211
[#3215]: https://github.com/datalad/datalad/issues/3215
[#3220]: https://github.com/datalad/datalad/issues/3220
[#3222]: https://github.com/datalad/datalad/issues/3222
[#3223]: https://github.com/datalad/datalad/issues/3223
[#3238]: https://github.com/datalad/datalad/issues/3238
[#3241]: https://github.com/datalad/datalad/issues/3241
[#3242]: https://github.com/datalad/datalad/issues/3242
[#3249]: https://github.com/datalad/datalad/issues/3249
[#3250]: https://github.com/datalad/datalad/issues/3250
[#3255]: https://github.com/datalad/datalad/issues/3255
[#3258]: https://github.com/datalad/datalad/issues/3258
[#3259]: https://github.com/datalad/datalad/issues/3259
[#3268]: https://github.com/datalad/datalad/issues/3268
[#3274]: https://github.com/datalad/datalad/issues/3274
[#3281]: https://github.com/datalad/datalad/issues/3281
[#3288]: https://github.com/datalad/datalad/issues/3288
[#3289]: https://github.com/datalad/datalad/issues/3289
[#3294]: https://github.com/datalad/datalad/issues/3294
[#3298]: https://github.com/datalad/datalad/issues/3298
[#3299]: https://github.com/datalad/datalad/issues/3299
[#3301]: https://github.com/datalad/datalad/issues/3301
[#3304]: https://github.com/datalad/datalad/issues/3304
[#3314]: https://github.com/datalad/datalad/issues/3314
[#3318]: https://github.com/datalad/datalad/issues/3318
[#3322]: https://github.com/datalad/datalad/issues/3322
[#3324]: https://github.com/datalad/datalad/issues/3324
[#3325]: https://github.com/datalad/datalad/issues/3325
[#3326]: https://github.com/datalad/datalad/issues/3326
[#3329]: https://github.com/datalad/datalad/issues/3329
[#3330]: https://github.com/datalad/datalad/issues/3330
[#3332]: https://github.com/datalad/datalad/issues/3332
[#3334]: https://github.com/datalad/datalad/issues/3334
[#3336]: https://github.com/datalad/datalad/issues/3336
[#3340]: https://github.com/datalad/datalad/issues/3340
[#3343]: https://github.com/datalad/datalad/issues/3343
[#3347]: https://github.com/datalad/datalad/issues/3347
[#3353]: https://github.com/datalad/datalad/issues/3353
[#3362]: https://github.com/datalad/datalad/issues/3362
[#3364]: https://github.com/datalad/datalad/issues/3364
[#3365]: https://github.com/datalad/datalad/issues/3365
[#3366]: https://github.com/datalad/datalad/issues/3366
[#3374]: https://github.com/datalad/datalad/issues/3374
[#3378]: https://github.com/datalad/datalad/issues/3378
[#3383]: https://github.com/datalad/datalad/issues/3383
[#3396]: https://github.com/datalad/datalad/issues/3396
[#3398]: https://github.com/datalad/datalad/issues/3398
[#3400]: https://github.com/datalad/datalad/issues/3400
[#3401]: https://github.com/datalad/datalad/issues/3401
[#3403]: https://github.com/datalad/datalad/issues/3403
[#3407]: https://github.com/datalad/datalad/issues/3407
[#3425]: https://github.com/datalad/datalad/issues/3425
[#3429]: https://github.com/datalad/datalad/issues/3429
[#3435]: https://github.com/datalad/datalad/issues/3435
[#3439]: https://github.com/datalad/datalad/issues/3439
[#3440]: https://github.com/datalad/datalad/issues/3440
[#3444]: https://github.com/datalad/datalad/issues/3444
[#3447]: https://github.com/datalad/datalad/issues/3447
[#3458]: https://github.com/datalad/datalad/issues/3458
[#3459]: https://github.com/datalad/datalad/issues/3459
[#3460]: https://github.com/datalad/datalad/issues/3460
[#3470]: https://github.com/datalad/datalad/issues/3470
[#3475]: https://github.com/datalad/datalad/issues/3475
[#3476]: https://github.com/datalad/datalad/issues/3476
[#3479]: https://github.com/datalad/datalad/issues/3479
[#3492]: https://github.com/datalad/datalad/issues/3492
[#3493]: https://github.com/datalad/datalad/issues/3493
[#3498]: https://github.com/datalad/datalad/issues/3498
[#3499]: https://github.com/datalad/datalad/issues/3499
[#3508]: https://github.com/datalad/datalad/issues/3508
[#3516]: https://github.com/datalad/datalad/issues/3516
[#3518]: https://github.com/datalad/datalad/issues/3518
[#3524]: https://github.com/datalad/datalad/issues/3524
[#3525]: https://github.com/datalad/datalad/issues/3525
[#3527]: https://github.com/datalad/datalad/issues/3527
[#3531]: https://github.com/datalad/datalad/issues/3531
[#3534]: https://github.com/datalad/datalad/issues/3534
[#3538]: https://github.com/datalad/datalad/issues/3538
[#3546]: https://github.com/datalad/datalad/issues/3546
[#3547]: https://github.com/datalad/datalad/issues/3547
[#3552]: https://github.com/datalad/datalad/issues/3552
[#3555]: https://github.com/datalad/datalad/issues/3555
[#3561]: https://github.com/datalad/datalad/issues/3561
[#3562]: https://github.com/datalad/datalad/issues/3562
[#3570]: https://github.com/datalad/datalad/issues/3570
[#3574]: https://github.com/datalad/datalad/issues/3574
[#3576]: https://github.com/datalad/datalad/issues/3576
[#3579]: https://github.com/datalad/datalad/issues/3579
[#3582]: https://github.com/datalad/datalad/issues/3582
[#3586]: https://github.com/datalad/datalad/issues/3586
[#3587]: https://github.com/datalad/datalad/issues/3587
[#3591]: https://github.com/datalad/datalad/issues/3591
[#3594]: https://github.com/datalad/datalad/issues/3594
[#3597]: https://github.com/datalad/datalad/issues/3597
[#3600]: https://github.com/datalad/datalad/issues/3600
[#3602]: https://github.com/datalad/datalad/issues/3602
[#3616]: https://github.com/datalad/datalad/issues/3616
[#3622]: https://github.com/datalad/datalad/issues/3622
[#3624]: https://github.com/datalad/datalad/issues/3624
[#3626]: https://github.com/datalad/datalad/issues/3626
[#3629]: https://github.com/datalad/datalad/issues/3629
[#3631]: https://github.com/datalad/datalad/issues/3631
[#3646]: https://github.com/datalad/datalad/issues/3646
[#3648]: https://github.com/datalad/datalad/issues/3648
[#3656]: https://github.com/datalad/datalad/issues/3656
[#3667]: https://github.com/datalad/datalad/issues/3667
[#3678]: https://github.com/datalad/datalad/issues/3678
[#3680]: https://github.com/datalad/datalad/issues/3680
[#3682]: https://github.com/datalad/datalad/issues/3682
[#3688]: https://github.com/datalad/datalad/issues/3688
[#3692]: https://github.com/datalad/datalad/issues/3692
[#3693]: https://github.com/datalad/datalad/issues/3693
[#3695]: https://github.com/datalad/datalad/issues/3695
[#3700]: https://github.com/datalad/datalad/issues/3700
[#3701]: https://github.com/datalad/datalad/issues/3701
[#3702]: https://github.com/datalad/datalad/issues/3702
[#3704]: https://github.com/datalad/datalad/issues/3704
[#3705]: https://github.com/datalad/datalad/issues/3705
[#3712]: https://github.com/datalad/datalad/issues/3712
[#3715]: https://github.com/datalad/datalad/issues/3715
[#3719]: https://github.com/datalad/datalad/issues/3719
[#3728]: https://github.com/datalad/datalad/issues/3728
[#3743]: https://github.com/datalad/datalad/issues/3743
[#3746]: https://github.com/datalad/datalad/issues/3746
[#3747]: https://github.com/datalad/datalad/issues/3747
[#3749]: https://github.com/datalad/datalad/issues/3749
[#3751]: https://github.com/datalad/datalad/issues/3751
[#3754]: https://github.com/datalad/datalad/issues/3754
[#3761]: https://github.com/datalad/datalad/issues/3761
[#3765]: https://github.com/datalad/datalad/issues/3765
[#3768]: https://github.com/datalad/datalad/issues/3768
[#3769]: https://github.com/datalad/datalad/issues/3769
[#3770]: https://github.com/datalad/datalad/issues/3770
[#3772]: https://github.com/datalad/datalad/issues/3772
[#3775]: https://github.com/datalad/datalad/issues/3775
[#3776]: https://github.com/datalad/datalad/issues/3776
[#3777]: https://github.com/datalad/datalad/issues/3777
[#3780]: https://github.com/datalad/datalad/issues/3780
[#3787]: https://github.com/datalad/datalad/issues/3787
[#3791]: https://github.com/datalad/datalad/issues/3791
[#3793]: https://github.com/datalad/datalad/issues/3793
[#3794]: https://github.com/datalad/datalad/issues/3794
[#3797]: https://github.com/datalad/datalad/issues/3797
[#3798]: https://github.com/datalad/datalad/issues/3798
[#3799]: https://github.com/datalad/datalad/issues/3799
[#3803]: https://github.com/datalad/datalad/issues/3803
[#3804]: https://github.com/datalad/datalad/issues/3804
[#3807]: https://github.com/datalad/datalad/issues/3807
[#3812]: https://github.com/datalad/datalad/issues/3812
[#3815]: https://github.com/datalad/datalad/issues/3815
[#3817]: https://github.com/datalad/datalad/issues/3817
[#3821]: https://github.com/datalad/datalad/issues/3821
[#3828]: https://github.com/datalad/datalad/issues/3828
[#3831]: https://github.com/datalad/datalad/issues/3831
[#3834]: https://github.com/datalad/datalad/issues/3834
[#3842]: https://github.com/datalad/datalad/issues/3842
[#3850]: https://github.com/datalad/datalad/issues/3850
[#3851]: https://github.com/datalad/datalad/issues/3851
[#3854]: https://github.com/datalad/datalad/issues/3854
[#3856]: https://github.com/datalad/datalad/issues/3856
[#3860]: https://github.com/datalad/datalad/issues/3860
[#3862]: https://github.com/datalad/datalad/issues/3862
[#3863]: https://github.com/datalad/datalad/issues/3863
[#3871]: https://github.com/datalad/datalad/issues/3871
[#3873]: https://github.com/datalad/datalad/issues/3873
[#3877]: https://github.com/datalad/datalad/issues/3877
[#3880]: https://github.com/datalad/datalad/issues/3880
[#3888]: https://github.com/datalad/datalad/issues/3888
[#3892]: https://github.com/datalad/datalad/issues/3892
[#3903]: https://github.com/datalad/datalad/issues/3903
[#3904]: https://github.com/datalad/datalad/issues/3904
[#3906]: https://github.com/datalad/datalad/issues/3906
[#3907]: https://github.com/datalad/datalad/issues/3907
[#3911]: https://github.com/datalad/datalad/issues/3911
[#3926]: https://github.com/datalad/datalad/issues/3926
[#3927]: https://github.com/datalad/datalad/issues/3927
[#3931]: https://github.com/datalad/datalad/issues/3931
[#3935]: https://github.com/datalad/datalad/issues/3935
[#3940]: https://github.com/datalad/datalad/issues/3940
[#3954]: https://github.com/datalad/datalad/issues/3954
[#3955]: https://github.com/datalad/datalad/issues/3955
[#3958]: https://github.com/datalad/datalad/issues/3958
[#3959]: https://github.com/datalad/datalad/issues/3959
[#3960]: https://github.com/datalad/datalad/issues/3960
[#3963]: https://github.com/datalad/datalad/issues/3963
[#3970]: https://github.com/datalad/datalad/issues/3970
[#3971]: https://github.com/datalad/datalad/issues/3971
[#3974]: https://github.com/datalad/datalad/issues/3974
[#3975]: https://github.com/datalad/datalad/issues/3975
[#3976]: https://github.com/datalad/datalad/issues/3976
[#3979]: https://github.com/datalad/datalad/issues/3979
[#3996]: https://github.com/datalad/datalad/issues/3996
[#3999]: https://github.com/datalad/datalad/issues/3999
[#4002]: https://github.com/datalad/datalad/issues/4002
[#4022]: https://github.com/datalad/datalad/issues/4022
[#4036]: https://github.com/datalad/datalad/issues/4036
[#4037]: https://github.com/datalad/datalad/issues/4037
[#4041]: https://github.com/datalad/datalad/issues/4041
[#4045]: https://github.com/datalad/datalad/issues/4045
[#4046]: https://github.com/datalad/datalad/issues/4046
[#4049]: https://github.com/datalad/datalad/issues/4049
[#4050]: https://github.com/datalad/datalad/issues/4050
[#4057]: https://github.com/datalad/datalad/issues/4057
[#4060]: https://github.com/datalad/datalad/issues/4060
[#4064]: https://github.com/datalad/datalad/issues/4064
[#4065]: https://github.com/datalad/datalad/issues/4065
[#4070]: https://github.com/datalad/datalad/issues/4070
[#4073]: https://github.com/datalad/datalad/issues/4073
[#4078]: https://github.com/datalad/datalad/issues/4078
[#4080]: https://github.com/datalad/datalad/issues/4080
[#4081]: https://github.com/datalad/datalad/issues/4081
[#4087]: https://github.com/datalad/datalad/issues/4087
[#4091]: https://github.com/datalad/datalad/issues/4091
[#4099]: https://github.com/datalad/datalad/issues/4099
[#4106]: https://github.com/datalad/datalad/issues/4106
[#4124]: https://github.com/datalad/datalad/issues/4124
[#4140]: https://github.com/datalad/datalad/issues/4140
[#4147]: https://github.com/datalad/datalad/issues/4147
[#4156]: https://github.com/datalad/datalad/issues/4156
[#4157]: https://github.com/datalad/datalad/issues/4157
[#4158]: https://github.com/datalad/datalad/issues/4158
[#4159]: https://github.com/datalad/datalad/issues/4159
[#4167]: https://github.com/datalad/datalad/issues/4167
[#4168]: https://github.com/datalad/datalad/issues/4168
[#4169]: https://github.com/datalad/datalad/issues/4169
[#4170]: https://github.com/datalad/datalad/issues/4170
[#4171]: https://github.com/datalad/datalad/issues/4171
[#4172]: https://github.com/datalad/datalad/issues/4172
[#4174]: https://github.com/datalad/datalad/issues/4174
[#4175]: https://github.com/datalad/datalad/issues/4175
[#4187]: https://github.com/datalad/datalad/issues/4187
[#4194]: https://github.com/datalad/datalad/issues/4194
[#4196]: https://github.com/datalad/datalad/issues/4196
[#4200]: https://github.com/datalad/datalad/issues/4200
[#4203]: https://github.com/datalad/datalad/issues/4203
[#4206]: https://github.com/datalad/datalad/issues/4206
[#4212]: https://github.com/datalad/datalad/issues/4212
[#4214]: https://github.com/datalad/datalad/issues/4214
[#4235]: https://github.com/datalad/datalad/issues/4235
[#4239]: https://github.com/datalad/datalad/issues/4239
[#4243]: https://github.com/datalad/datalad/issues/4243
[#4245]: https://github.com/datalad/datalad/issues/4245
[#4257]: https://github.com/datalad/datalad/issues/4257
[#4260]: https://github.com/datalad/datalad/issues/4260
[#4262]: https://github.com/datalad/datalad/issues/4262
[#4268]: https://github.com/datalad/datalad/issues/4268
[#4273]: https://github.com/datalad/datalad/issues/4273
[#4274]: https://github.com/datalad/datalad/issues/4274
[#4276]: https://github.com/datalad/datalad/issues/4276
[#4285]: https://github.com/datalad/datalad/issues/4285
[#4290]: https://github.com/datalad/datalad/issues/4290
[#4291]: https://github.com/datalad/datalad/issues/4291
[#4292]: https://github.com/datalad/datalad/issues/4292
[#4296]: https://github.com/datalad/datalad/issues/4296
[#4301]: https://github.com/datalad/datalad/issues/4301
[#4303]: https://github.com/datalad/datalad/issues/4303
[#4304]: https://github.com/datalad/datalad/issues/4304
[#4305]: https://github.com/datalad/datalad/issues/4305
[#4306]: https://github.com/datalad/datalad/issues/4306
[#4308]: https://github.com/datalad/datalad/issues/4308
[#4314]: https://github.com/datalad/datalad/issues/4314
[#4315]: https://github.com/datalad/datalad/issues/4315
[#4316]: https://github.com/datalad/datalad/issues/4316
[#4317]: https://github.com/datalad/datalad/issues/4317
[#4319]: https://github.com/datalad/datalad/issues/4319
[#4321]: https://github.com/datalad/datalad/issues/4321
[#4323]: https://github.com/datalad/datalad/issues/4323
[#4324]: https://github.com/datalad/datalad/issues/4324
[#4326]: https://github.com/datalad/datalad/issues/4326
[#4328]: https://github.com/datalad/datalad/issues/4328
[#4330]: https://github.com/datalad/datalad/issues/4330
[#4331]: https://github.com/datalad/datalad/issues/4331
[#4332]: https://github.com/datalad/datalad/issues/4332
[#4337]: https://github.com/datalad/datalad/issues/4337
[#4338]: https://github.com/datalad/datalad/issues/4338
[#4342]: https://github.com/datalad/datalad/issues/4342
[#4348]: https://github.com/datalad/datalad/issues/4348
[#4354]: https://github.com/datalad/datalad/issues/4354
[#4361]: https://github.com/datalad/datalad/issues/4361
[#4367]: https://github.com/datalad/datalad/issues/4367
[#4370]: https://github.com/datalad/datalad/issues/4370
[#4375]: https://github.com/datalad/datalad/issues/4375
[#4382]: https://github.com/datalad/datalad/issues/4382
[#4398]: https://github.com/datalad/datalad/issues/4398
[#4400]: https://github.com/datalad/datalad/issues/4400
[#4409]: https://github.com/datalad/datalad/issues/4409
[#4420]: https://github.com/datalad/datalad/issues/4420
[#4421]: https://github.com/datalad/datalad/issues/4421
[#4426]: https://github.com/datalad/datalad/issues/4426
[#4430]: https://github.com/datalad/datalad/issues/4430
[#4431]: https://github.com/datalad/datalad/issues/4431
[#4435]: https://github.com/datalad/datalad/issues/4435
[#4438]: https://github.com/datalad/datalad/issues/4438
[#4439]: https://github.com/datalad/datalad/issues/4439
[#4441]: https://github.com/datalad/datalad/issues/4441
[#4448]: https://github.com/datalad/datalad/issues/4448
[#4456]: https://github.com/datalad/datalad/issues/4456
[#4459]: https://github.com/datalad/datalad/issues/4459
[#4460]: https://github.com/datalad/datalad/issues/4460
[#4463]: https://github.com/datalad/datalad/issues/4463
[#4464]: https://github.com/datalad/datalad/issues/4464
[#4471]: https://github.com/datalad/datalad/issues/4471
[#4477]: https://github.com/datalad/datalad/issues/4477
[#4480]: https://github.com/datalad/datalad/issues/4480
[#4481]: https://github.com/datalad/datalad/issues/4481
[#4504]: https://github.com/datalad/datalad/issues/4504
[#4526]: https://github.com/datalad/datalad/issues/4526
[#4529]: https://github.com/datalad/datalad/issues/4529
[#4543]: https://github.com/datalad/datalad/issues/4543
[#4544]: https://github.com/datalad/datalad/issues/4544
[#4549]: https://github.com/datalad/datalad/issues/4549
[#4552]: https://github.com/datalad/datalad/issues/4552
[#4553]: https://github.com/datalad/datalad/issues/4553
[#4560]: https://github.com/datalad/datalad/issues/4560
[#4568]: https://github.com/datalad/datalad/issues/4568
[#4581]: https://github.com/datalad/datalad/issues/4581
[#4583]: https://github.com/datalad/datalad/issues/4583
[#4597]: https://github.com/datalad/datalad/issues/4597
[#4617]: https://github.com/datalad/datalad/issues/4617
[#4619]: https://github.com/datalad/datalad/issues/4619
[#4620]: https://github.com/datalad/datalad/issues/4620
[#4650]: https://github.com/datalad/datalad/issues/4650
[#4657]: https://github.com/datalad/datalad/issues/4657
[#4666]: https://github.com/datalad/datalad/issues/4666
[#4669]: https://github.com/datalad/datalad/issues/4669
[#4673]: https://github.com/datalad/datalad/issues/4673
[#4674]: https://github.com/datalad/datalad/issues/4674
[#4675]: https://github.com/datalad/datalad/issues/4675
[#4682]: https://github.com/datalad/datalad/issues/4682
[#4683]: https://github.com/datalad/datalad/issues/4683
[#4684]: https://github.com/datalad/datalad/issues/4684
[#4687]: https://github.com/datalad/datalad/issues/4687
[#4692]: https://github.com/datalad/datalad/issues/4692
[#4695]: https://github.com/datalad/datalad/issues/4695
[#4696]: https://github.com/datalad/datalad/issues/4696
[#4699]: https://github.com/datalad/datalad/issues/4699
[#4703]: https://github.com/datalad/datalad/issues/4703
[#4729]: https://github.com/datalad/datalad/issues/4729
[#4736]: https://github.com/datalad/datalad/issues/4736
[#4745]: https://github.com/datalad/datalad/issues/4745
[#4746]: https://github.com/datalad/datalad/issues/4746
[#4749]: https://github.com/datalad/datalad/issues/4749
[#4757]: https://github.com/datalad/datalad/issues/4757
[#4759]: https://github.com/datalad/datalad/issues/4759
[#4760]: https://github.com/datalad/datalad/issues/4760
[#4763]: https://github.com/datalad/datalad/issues/4763
[#4764]: https://github.com/datalad/datalad/issues/4764
[#4769]: https://github.com/datalad/datalad/issues/4769
[#4775]: https://github.com/datalad/datalad/issues/4775
[#4786]: https://github.com/datalad/datalad/issues/4786
[#4788]: https://github.com/datalad/datalad/issues/4788
[#4790]: https://github.com/datalad/datalad/issues/4790
[#4792]: https://github.com/datalad/datalad/issues/4792
[#4793]: https://github.com/datalad/datalad/issues/4793
[#4806]: https://github.com/datalad/datalad/issues/4806
[#4807]: https://github.com/datalad/datalad/issues/4807
[#4816]: https://github.com/datalad/datalad/issues/4816
[#4817]: https://github.com/datalad/datalad/issues/4817
[#4821]: https://github.com/datalad/datalad/issues/4821
[#4824]: https://github.com/datalad/datalad/issues/4824
[#4828]: https://github.com/datalad/datalad/issues/4828
[#4829]: https://github.com/datalad/datalad/issues/4829
[#4834]: https://github.com/datalad/datalad/issues/4834
[#4835]: https://github.com/datalad/datalad/issues/4835
[#4845]: https://github.com/datalad/datalad/issues/4845
[#4853]: https://github.com/datalad/datalad/issues/4853
[#4855]: https://github.com/datalad/datalad/issues/4855
[#4866]: https://github.com/datalad/datalad/issues/4866
[#4867]: https://github.com/datalad/datalad/issues/4867
[#4868]: https://github.com/datalad/datalad/issues/4868
[#4877]: https://github.com/datalad/datalad/issues/4877
[#4879]: https://github.com/datalad/datalad/issues/4879
[#4896]: https://github.com/datalad/datalad/issues/4896
[#4899]: https://github.com/datalad/datalad/issues/4899
[#4900]: https://github.com/datalad/datalad/issues/4900
[#4904]: https://github.com/datalad/datalad/issues/4904
[#4908]: https://github.com/datalad/datalad/issues/4908
[#4911]: https://github.com/datalad/datalad/issues/4911
[#4924]: https://github.com/datalad/datalad/issues/4924
[#4926]: https://github.com/datalad/datalad/issues/4926
[#4927]: https://github.com/datalad/datalad/issues/4927
[#4931]: https://github.com/datalad/datalad/issues/4931
[#4952]: https://github.com/datalad/datalad/issues/4952
[#4953]: https://github.com/datalad/datalad/issues/4953
[#4955]: https://github.com/datalad/datalad/issues/4955
[#4957]: https://github.com/datalad/datalad/issues/4957
[#4963]: https://github.com/datalad/datalad/issues/4963
[#4966]: https://github.com/datalad/datalad/issues/4966
[#4977]: https://github.com/datalad/datalad/issues/4977
[#4982]: https://github.com/datalad/datalad/issues/4982
[#4985]: https://github.com/datalad/datalad/issues/4985
[#4991]: https://github.com/datalad/datalad/issues/4991
[#4996]: https://github.com/datalad/datalad/issues/4996
[#5001]: https://github.com/datalad/datalad/issues/5001
[#5002]: https://github.com/datalad/datalad/issues/5002
[#5008]: https://github.com/datalad/datalad/issues/5008
[#5010]: https://github.com/datalad/datalad/issues/5010
[#5017]: https://github.com/datalad/datalad/issues/5017
[#5022]: https://github.com/datalad/datalad/issues/5022
[#5025]: https://github.com/datalad/datalad/issues/5025
[#5026]: https://github.com/datalad/datalad/issues/5026
[#5035]: https://github.com/datalad/datalad/issues/5035
[#5042]: https://github.com/datalad/datalad/issues/5042
[#5045]: https://github.com/datalad/datalad/issues/5045
[#5049]: https://github.com/datalad/datalad/issues/5049
[#5051]: https://github.com/datalad/datalad/issues/5051
[#5057]: https://github.com/datalad/datalad/issues/5057
[#5060]: https://github.com/datalad/datalad/issues/5060
[#5067]: https://github.com/datalad/datalad/issues/5067
[#5070]: https://github.com/datalad/datalad/issues/5070
[#5076]: https://github.com/datalad/datalad/issues/5076
[#5081]: https://github.com/datalad/datalad/issues/5081
[#5090]: https://github.com/datalad/datalad/issues/5090
[#5091]: https://github.com/datalad/datalad/issues/5091
[#5106]: https://github.com/datalad/datalad/issues/5106
[#5108]: https://github.com/datalad/datalad/issues/5108
[#5113]: https://github.com/datalad/datalad/issues/5113
[#5119]: https://github.com/datalad/datalad/issues/5119
[#5121]: https://github.com/datalad/datalad/issues/5121
[#5125]: https://github.com/datalad/datalad/issues/5125
[#5127]: https://github.com/datalad/datalad/issues/5127
[#5128]: https://github.com/datalad/datalad/issues/5128
[#5129]: https://github.com/datalad/datalad/issues/5129
[#5136]: https://github.com/datalad/datalad/issues/5136
[#5141]: https://github.com/datalad/datalad/issues/5141
[#5142]: https://github.com/datalad/datalad/issues/5142
[#5146]: https://github.com/datalad/datalad/issues/5146
[#5148]: https://github.com/datalad/datalad/issues/5148
[#5151]: https://github.com/datalad/datalad/issues/5151
[#5156]: https://github.com/datalad/datalad/issues/5156
[#5163]: https://github.com/datalad/datalad/issues/5163
[#5184]: https://github.com/datalad/datalad/issues/5184
[#5194]: https://github.com/datalad/datalad/issues/5194
[#5200]: https://github.com/datalad/datalad/issues/5200
[#5201]: https://github.com/datalad/datalad/issues/5201
[#5214]: https://github.com/datalad/datalad/issues/5214
[#5218]: https://github.com/datalad/datalad/issues/5218
[#5219]: https://github.com/datalad/datalad/issues/5219
[#5229]: https://github.com/datalad/datalad/issues/5229
[#5238]: https://github.com/datalad/datalad/issues/5238
[#5241]: https://github.com/datalad/datalad/issues/5241
[#5254]: https://github.com/datalad/datalad/issues/5254
[#5255]: https://github.com/datalad/datalad/issues/5255
[#5258]: https://github.com/datalad/datalad/issues/5258
[#5259]: https://github.com/datalad/datalad/issues/5259
[#5269]: https://github.com/datalad/datalad/issues/5269
[#5276]: https://github.com/datalad/datalad/issues/5276
[#5278]: https://github.com/datalad/datalad/issues/5278
[#5285]: https://github.com/datalad/datalad/issues/5285
[#5290]: https://github.com/datalad/datalad/issues/5290
[#5328]: https://github.com/datalad/datalad/issues/5328
[#5332]: https://github.com/datalad/datalad/issues/5332
[#5342]: https://github.com/datalad/datalad/issues/5342
[#5344]: https://github.com/datalad/datalad/issues/5344
[#5346]: https://github.com/datalad/datalad/issues/5346
[#5350]: https://github.com/datalad/datalad/issues/5350
[#5367]: https://github.com/datalad/datalad/issues/5367
[#5389]: https://github.com/datalad/datalad/issues/5389
[#5391]: https://github.com/datalad/datalad/issues/5391
[#5415]: https://github.com/datalad/datalad/issues/5415
[#5416]: https://github.com/datalad/datalad/issues/5416
[#5421]: https://github.com/datalad/datalad/issues/5421
[#5425]: https://github.com/datalad/datalad/issues/5425
[#5430]: https://github.com/datalad/datalad/issues/5430
[#5431]: https://github.com/datalad/datalad/issues/5431
[#5436]: https://github.com/datalad/datalad/issues/5436
[#5438]: https://github.com/datalad/datalad/issues/5438
[#5441]: https://github.com/datalad/datalad/issues/5441
[#5453]: https://github.com/datalad/datalad/issues/5453
[#5458]: https://github.com/datalad/datalad/issues/5458
[#5461]: https://github.com/datalad/datalad/issues/5461
[#5466]: https://github.com/datalad/datalad/issues/5466
[#5474]: https://github.com/datalad/datalad/issues/5474
[#5476]: https://github.com/datalad/datalad/issues/5476
[#5480]: https://github.com/datalad/datalad/issues/5480
[#5488]: https://github.com/datalad/datalad/issues/5488
[#5492]: https://github.com/datalad/datalad/issues/5492
[#5505]: https://github.com/datalad/datalad/issues/5505
[#5509]: https://github.com/datalad/datalad/issues/5509
[#5512]: https://github.com/datalad/datalad/issues/5512
[#5525]: https://github.com/datalad/datalad/issues/5525
[#5531]: https://github.com/datalad/datalad/issues/5531
[#5533]: https://github.com/datalad/datalad/issues/5533
[#5534]: https://github.com/datalad/datalad/issues/5534
[#5536]: https://github.com/datalad/datalad/issues/5536
[#5539]: https://github.com/datalad/datalad/issues/5539
[#5543]: https://github.com/datalad/datalad/issues/5543
[#5544]: https://github.com/datalad/datalad/issues/5544
[#5550]: https://github.com/datalad/datalad/issues/5550
[#5551]: https://github.com/datalad/datalad/issues/5551
[#5552]: https://github.com/datalad/datalad/issues/5552
[#5555]: https://github.com/datalad/datalad/issues/5555
[#5558]: https://github.com/datalad/datalad/issues/5558
[#5559]: https://github.com/datalad/datalad/issues/5559
[#5560]: https://github.com/datalad/datalad/issues/5560
[#5569]: https://github.com/datalad/datalad/issues/5569
[#5572]: https://github.com/datalad/datalad/issues/5572
[#5577]: https://github.com/datalad/datalad/issues/5577
[#5580]: https://github.com/datalad/datalad/issues/5580
[#5592]: https://github.com/datalad/datalad/issues/5592
[#5594]: https://github.com/datalad/datalad/issues/5594
[#5603]: https://github.com/datalad/datalad/issues/5603
[#5607]: https://github.com/datalad/datalad/issues/5607
[#5609]: https://github.com/datalad/datalad/issues/5609
[#5612]: https://github.com/datalad/datalad/issues/5612
[#5630]: https://github.com/datalad/datalad/issues/5630
[#5632]: https://github.com/datalad/datalad/issues/5632
[#5639]: https://github.com/datalad/datalad/issues/5639
[#5655]: https://github.com/datalad/datalad/issues/5655
[#5667]: https://github.com/datalad/datalad/issues/5667
[#5675]: https://github.com/datalad/datalad/issues/5675
[#5680]: https://github.com/datalad/datalad/issues/5680
[#5681]: https://github.com/datalad/datalad/issues/5681
[#5682]: https://github.com/datalad/datalad/issues/5682
[#5683]: https://github.com/datalad/datalad/issues/5683
[#5689]: https://github.com/datalad/datalad/issues/5689
[#5692]: https://github.com/datalad/datalad/issues/5692
[#5693]: https://github.com/datalad/datalad/issues/5693
[#5696]: https://github.com/datalad/datalad/issues/5696
[#5698]: https://github.com/datalad/datalad/issues/5698
[#5708]: https://github.com/datalad/datalad/issues/5708
[#5726]: https://github.com/datalad/datalad/issues/5726
[#5738]: https://github.com/datalad/datalad/issues/5738
[#5740]: https://github.com/datalad/datalad/issues/5740
[#5749]: https://github.com/datalad/datalad/issues/5749
[#5760]: https://github.com/datalad/datalad/issues/5760
[#5777]: https://github.com/datalad/datalad/issues/5777
[#5789]: https://github.com/datalad/datalad/issues/5789
[#5792]: https://github.com/datalad/datalad/issues/5792
[#5803]: https://github.com/datalad/datalad/issues/5803
[#5804]: https://github.com/datalad/datalad/issues/5804
[#5805]: https://github.com/datalad/datalad/issues/5805
[#5823]: https://github.com/datalad/datalad/issues/5823
[#5837]: https://github.com/datalad/datalad/issues/5837
[#5847]: https://github.com/datalad/datalad/issues/5847
[#5884]: https://github.com/datalad/datalad/issues/5884
[#5892]: https://github.com/datalad/datalad/issues/5892
[#5902]: https://github.com/datalad/datalad/issues/5902
[#5904]: https://github.com/datalad/datalad/issues/5904
[#5907]: https://github.com/datalad/datalad/issues/5907
[#5913]: https://github.com/datalad/datalad/issues/5913
[#5915]: https://github.com/datalad/datalad/issues/5915
[#5956]: https://github.com/datalad/datalad/issues/5956
     ____            _             _                   _ 
    |  _ \    __ _  | |_    __ _  | |       __ _    __| |
    | | | |  / _` | | __|  / _` | | |      / _` |  / _` |
    | |_| | | (_| | | |_  | (_| | | |___  | (_| | | (_| |
    |____/   \__,_|  \__|  \__,_| |_____|  \__,_|  \__,_|
                                                  Read me

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03262/status.svg)](https://doi.org/10.21105/joss.03262)
[![Travis tests status](https://app.travis-ci.com/datalad/datalad.svg?branch=master)](https://app.travis-ci.com/datalad/datalad)
[![Build status](https://ci.appveyor.com/api/projects/status/github/datalad/datalad?branch=master&svg=true)](https://ci.appveyor.com/project/mih/datalad/branch/master)
[![codecov.io](https://codecov.io/github/datalad/datalad/coverage.svg?branch=master)](https://codecov.io/github/datalad/datalad?branch=master)
[![Documentation](https://readthedocs.org/projects/datalad/badge/?version=latest)](http://datalad.rtfd.org)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/release/datalad/datalad.svg)](https://GitHub.com/datalad/datalad/releases/)
[![PyPI version fury.io](https://badge.fury.io/py/datalad.svg)](https://pypi.python.org/pypi/datalad/)
[![Supported Python versions](https://img.shields.io/pypi/pyversions/datalad)](https://pypi.org/project/datalad/)
[![Testimonials 4](https://img.shields.io/badge/testimonials-4-brightgreen.svg)](https://github.com/datalad/datalad/wiki/Testimonials)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/667)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://github.com/datalad/datalad/blob/master/CODE_OF_CONDUCT.md)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5624980.svg)](https://doi.org/10.5281/zenodo.5624980)
<!-- ALL-CONTRIBUTORS-BADGE:START - Do not remove or modify this section -->
[![All Contributors](https://img.shields.io/badge/all_contributors-38-orange.svg?style=flat-square)](#contributors-)
<!-- ALL-CONTRIBUTORS-BADGE:END -->


# 10000-ft. overview

DataLad makes data management and data distribution more accessible.
To do that, it stands on the shoulders of [Git] and [Git-annex] to deliver a
decentralized system for data exchange. This includes automated ingestion of
data from online portals and exposing it in readily usable form as Git(-annex)
repositories, so-called datasets. The actual data storage and permission
management, however, remains with the original data providers.

The full documentation is available at http://docs.datalad.org and
http://handbook.datalad.org provides a hands-on crash-course on DataLad.

# Extensions

A number of extensions are available that provide additional functionality for
DataLad. Extensions are separate packages that are to be installed in addition
to DataLad. In order to install DataLad customized for a particular domain, one
can simply install an extension directly, and DataLad itself will be
automatically installed with it. An [annotated list of
extensions](http://handbook.datalad.org/extension_pkgs.html) is available in
the [DataLad handbook](http://handbook.datalad.org).


# Support

The documentation for this project is found here:
http://docs.datalad.org

All bugs, concerns, and enhancement requests for this software can be submitted here:
https://github.com/datalad/datalad/issues

If you have a problem or would like to ask a question about how to use DataLad,
please [submit a question to
NeuroStars.org](https://neurostars.org/new-topic?body=-%20Please%20describe%20the%20problem.%0A-%20What%20steps%20will%20reproduce%20the%20problem%3F%0A-%20What%20version%20of%20DataLad%20are%20you%20using%20%28run%20%60datalad%20--version%60%29%3F%20On%20what%20operating%20system%20%28consider%20running%20%60datalad%20plugin%20wtf%60%29%3F%0A-%20Please%20provide%20any%20additional%20information%20below.%0A-%20Have%20you%20had%20any%20luck%20using%20DataLad%20before%3F%20%28Sometimes%20we%20get%20tired%20of%20reading%20bug%20reports%20all%20day%20and%20a%20lil'%20positive%20end%20note%20does%20wonders%29&tags=datalad)
with a `datalad` tag.  NeuroStars.org is a platform similar to StackOverflow
but dedicated to neuroinformatics.

All previous DataLad questions are available here:
http://neurostars.org/tags/datalad/


# Installation

## Debian-based systems

On Debian-based systems, we recommend enabling [NeuroDebian], via which we
provide recent releases of DataLad. Once enabled, just do:

    apt-get install datalad

## Gentoo-based systems

On Gentoo-based systems (i.e. all systems whose package manager can parse ebuilds as per the [Package Manager Specification]), we recommend [enabling the ::science overlay], via which we
provide recent releases of DataLad. Once enabled, just run:

    emerge datalad

## Other Linux'es via conda

    conda install -c conda-forge datalad

will install the most recently released version, and release candidates are
available via

    conda install -c conda-forge/label/rc datalad

## Other Linux'es, macOS via pip

Before you install this package, please make sure that you [install a recent
version of git-annex](https://git-annex.branchable.com/install).  Afterwards,
install the latest version of `datalad` from
[PyPI](https://pypi.org/project/datalad). It is recommended to use
a dedicated [virtualenv](https://virtualenv.pypa.io):

    # Create and enter a new virtual environment (optional)
    virtualenv --python=python3 ~/env/datalad
    . ~/env/datalad/bin/activate

    # Install from PyPI
    pip install datalad

By default, installation via pip installs the core functionality of DataLad,
allowing for managing datasets etc.  Additional installation schemes
are available, so you can request enhanced installation via
`pip install datalad[SCHEME]`, where `SCHEME` could be:

- `tests`
     to also install dependencies used by DataLad's battery of unit tests
- `full`
     to install all dependencies.

More details on installation and initial configuration can be found in the
[DataLad Handbook: Installation].

# License

MIT/Expat


# Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) if you are interested in internals or
contributing to the project. 

## Acknowledgements

DataLad development is supported by a US-German collaboration in computational
neuroscience (CRCNS) project "DataGit: converging catalogues, warehouses, and
deployment logistics into a federated 'data distribution'" (Halchenko/Hanke),
co-funded by the US National Science Foundation (NSF 1429999) and the German
Federal Ministry of Education and Research (BMBF 01GQ1411). Additional support
is provided by the German federal state of Saxony-Anhalt and the European
Regional Development Fund (ERDF), Project: Center for Behavioral Brain
Sciences, Imaging Platform.  This work is further facilitated by the ReproNim
project (NIH 1P41EB019936-01A1). Mac mini instance for development is provided
by [MacStadium](https://www.macstadium.com/).

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tr>
    <td align="center"><a href="https://github.com/glalteva"><img src="https://avatars2.githubusercontent.com/u/14296143?v=4?s=100" width="100px;" alt=""/><br /><sub><b>glalteva</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=glalteva" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/adswa"><img src="https://avatars1.githubusercontent.com/u/29738718?v=4?s=100" width="100px;" alt=""/><br /><sub><b>adswa</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=adswa" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/chrhaeusler"><img src="https://avatars0.githubusercontent.com/u/8115807?v=4?s=100" width="100px;" alt=""/><br /><sub><b>chrhaeusler</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=chrhaeusler" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/soichih"><img src="https://avatars3.githubusercontent.com/u/923896?v=4?s=100" width="100px;" alt=""/><br /><sub><b>soichih</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=soichih" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mvdoc"><img src="https://avatars1.githubusercontent.com/u/6150554?v=4?s=100" width="100px;" alt=""/><br /><sub><b>mvdoc</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=mvdoc" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mih"><img src="https://avatars1.githubusercontent.com/u/136479?v=4?s=100" width="100px;" alt=""/><br /><sub><b>mih</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=mih" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/yarikoptic"><img src="https://avatars3.githubusercontent.com/u/39889?v=4?s=100" width="100px;" alt=""/><br /><sub><b>yarikoptic</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=yarikoptic" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/loj"><img src="https://avatars2.githubusercontent.com/u/15157717?v=4?s=100" width="100px;" alt=""/><br /><sub><b>loj</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=loj" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/feilong"><img src="https://avatars2.githubusercontent.com/u/2242261?v=4?s=100" width="100px;" alt=""/><br /><sub><b>feilong</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=feilong" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/jhpoelen"><img src="https://avatars2.githubusercontent.com/u/1084872?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jhpoelen</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=jhpoelen" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/andycon"><img src="https://avatars1.githubusercontent.com/u/3965889?v=4?s=100" width="100px;" alt=""/><br /><sub><b>andycon</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=andycon" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/nicholsn"><img src="https://avatars3.githubusercontent.com/u/463344?v=4?s=100" width="100px;" alt=""/><br /><sub><b>nicholsn</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=nicholsn" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/adelavega"><img src="https://avatars0.githubusercontent.com/u/2774448?v=4?s=100" width="100px;" alt=""/><br /><sub><b>adelavega</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=adelavega" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/kskyten"><img src="https://avatars0.githubusercontent.com/u/4163878?v=4?s=100" width="100px;" alt=""/><br /><sub><b>kskyten</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=kskyten" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/TheChymera"><img src="https://avatars2.githubusercontent.com/u/950524?v=4?s=100" width="100px;" alt=""/><br /><sub><b>TheChymera</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=TheChymera" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/effigies"><img src="https://avatars0.githubusercontent.com/u/83442?v=4?s=100" width="100px;" alt=""/><br /><sub><b>effigies</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=effigies" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/jgors"><img src="https://avatars1.githubusercontent.com/u/386585?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jgors</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=jgors" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/debanjum"><img src="https://avatars1.githubusercontent.com/u/6413477?v=4?s=100" width="100px;" alt=""/><br /><sub><b>debanjum</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=debanjum" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/nellh"><img src="https://avatars3.githubusercontent.com/u/11369795?v=4?s=100" width="100px;" alt=""/><br /><sub><b>nellh</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=nellh" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/emdupre"><img src="https://avatars3.githubusercontent.com/u/15017191?v=4?s=100" width="100px;" alt=""/><br /><sub><b>emdupre</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=emdupre" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/aqw"><img src="https://avatars0.githubusercontent.com/u/765557?v=4?s=100" width="100px;" alt=""/><br /><sub><b>aqw</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=aqw" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/vsoch"><img src="https://avatars0.githubusercontent.com/u/814322?v=4?s=100" width="100px;" alt=""/><br /><sub><b>vsoch</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=vsoch" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/kyleam"><img src="https://avatars2.githubusercontent.com/u/1297788?v=4?s=100" width="100px;" alt=""/><br /><sub><b>kyleam</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=kyleam" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/driusan"><img src="https://avatars0.githubusercontent.com/u/498329?v=4?s=100" width="100px;" alt=""/><br /><sub><b>driusan</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=driusan" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/overlake333"><img src="https://avatars1.githubusercontent.com/u/28018084?v=4?s=100" width="100px;" alt=""/><br /><sub><b>overlake333</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=overlake333" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/akeshavan"><img src="https://avatars0.githubusercontent.com/u/972008?v=4?s=100" width="100px;" alt=""/><br /><sub><b>akeshavan</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=akeshavan" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/jwodder"><img src="https://avatars1.githubusercontent.com/u/98207?v=4?s=100" width="100px;" alt=""/><br /><sub><b>jwodder</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=jwodder" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/bpoldrack"><img src="https://avatars2.githubusercontent.com/u/10498301?v=4?s=100" width="100px;" alt=""/><br /><sub><b>bpoldrack</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=bpoldrack" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://github.com/yetanothertestuser"><img src="https://avatars0.githubusercontent.com/u/19335420?v=4?s=100" width="100px;" alt=""/><br /><sub><b>yetanothertestuser</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=yetanothertestuser" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/christian-monch"><img src="https://avatars.githubusercontent.com/u/17925232?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Christian Mönch</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=christian-monch" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mattcieslak"><img src="https://avatars.githubusercontent.com/u/170026?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Matt Cieslak</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=mattcieslak" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mikapfl"><img src="https://avatars.githubusercontent.com/u/7226087?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Mika Pflüger</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=mikapfl" title="Code">💻</a></td>
    <td align="center"><a href="https://me.ypid.de/"><img src="https://avatars.githubusercontent.com/u/1301158?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Robin Schneider</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=ypid" title="Code">💻</a></td>
    <td align="center"><a href="https://orcid.org/0000-0003-4652-3758"><img src="https://avatars.githubusercontent.com/u/7570456?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Sin Kim</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=AKSoo" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/DisasterMo"><img src="https://avatars.githubusercontent.com/u/49207524?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michael Burgardt</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=DisasterMo" title="Code">💻</a></td>
  </tr>
  <tr>
    <td align="center"><a href="https://remi-gau.github.io/"><img src="https://avatars.githubusercontent.com/u/6961185?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Remi Gau</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=Remi-Gau" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/mslw"><img src="https://avatars.githubusercontent.com/u/11985212?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Michał Szczepanik</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=mslw" title="Code">💻</a></td>
    <td align="center"><a href="https://github.com/bpinsard"><img src="https://avatars.githubusercontent.com/u/1155388?v=4?s=100" width="100px;" alt=""/><br /><sub><b>Basile</b></sub></a><br /><a href="https://github.com/datalad/datalad/commits?author=bpinsard" title="Code">💻</a></td>
  </tr>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

[![macstadium](https://uploads-ssl.webflow.com/5ac3c046c82724970fc60918/5c019d917bba312af7553b49_MacStadium-developerlogo.png)](https://www.macstadium.com/)

[Git]: https://git-scm.com
[Git-annex]: http://git-annex.branchable.com
[setup.py]: https://github.com/datalad/datalad/blob/master/setup.py
[NeuroDebian]: http://neuro.debian.net
[Package Manager Specification]: https://projects.gentoo.org/pms/latest/pms.html
[enabling the ::science overlay]: https://github.com/gentoo/sci#manual-install-

[DataLad Handbook: Installation]: http://handbook.datalad.org/en/latest/intro/installation.html

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
reported to the community leaders responsible for enforcement.
Please contact either Adina Wagner at adina.wagner@t-online.de or Michael Hanke
at michael.hanke@gmail.com.
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
[translations]: https://www.contributor-covenant.org/translations# Contributing to DataLad

[gh-datalad]: http://github.com/datalad/datalad

## Files organization

- [datalad/](./datalad) is the main Python module where major development is happening,
  with major submodules being:
    - `cmdline/` - helpers for accessing `interface/` functionality from
     command line
    - `customremotes/` - custom special remotes for annex provided by datalad
    - `downloaders/` - support for accessing data from various sources (e.g.
      http, S3, XNAT) via a unified interface.
        - `configs/` - specifications for known data providers and associated
          credentials
    - `interface/` - high level interface functions which get exposed via
      command line (`cmdline/`) or Python (`datalad.api`).
    - `tests/` - some unit- and regression- tests (more could be found under
      `tests/` of corresponding submodules. See [Tests](#tests))
        - [utils.py](./datalad/tests/utils.py) provides convenience helpers used by unit-tests such as
          `@with_tree`, `@serve_path_via_http` and other decorators
    - `ui/` - user-level interactions, such as messages about errors, warnings,
      progress reports, AND when supported by available frontend --
      interactive dialogs
    - `support/` - various support modules, e.g. for git/git-annex interfaces,
      constraints for the `interface/`, etc
- [benchmarks/](./benchmarks) - [asv] benchmarks suite (see [Benchmarking](#benchmarking))
- [docs/](./docs) - yet to be heavily populated documentation
    - `bash-completions` - bash and zsh completion setup for datalad (just
      `source` it)
- [fixtures/](./fixtures) currently not under git, contains generated by vcr fixtures
- [sandbox/](./sandbox) - various scripts and prototypes which are not part of
  the main/distributed with releases codebase
- [tools/](./tools) contains helper utilities used during development, testing, and
  benchmarking of DataLad.  Implemented in any most appropriate language
  (Python, bash, etc.)

Whenever a new top-level file or folder is added to the repository, it should
be listed in `MANIFEST.in` so that it will be either included in or excluded
from source distributions as appropriate.  [See
here](https://packaging.python.org/guides/using-manifest-in/) for information
about writing a `MANIFEST.in`.

## How to contribute

The preferred way to contribute to the DataLad code base is
to fork the [main repository][gh-datalad] on GitHub.  Here
we outline the workflow used by the developers:


0. Have a clone of our main [project repository][gh-datalad] as `origin`
   remote in your git:

          git clone git://github.com/datalad/datalad

1. Fork the [project repository][gh-datalad]: click on the 'Fork'
   button near the top of the page.  This creates a copy of the code
   base under your account on the GitHub server.

2. Add your forked clone as a remote to the local clone you already have on your
   local disk:

          git remote add gh-YourLogin git@github.com:YourLogin/datalad.git
          git fetch gh-YourLogin

    To ease addition of other github repositories as remotes, here is
    a little bash function/script to add to your `~/.bashrc`:

        ghremote () {
                url="$1"
                proj=${url##*/}
                url_=${url%/*}
                login=${url_##*/}
                git remote add gh-$login $url
                git fetch gh-$login
        }

    thus you could simply run:

         ghremote git@github.com:YourLogin/datalad.git

    to add the above `gh-YourLogin` remote.  Additional handy aliases
    such as `ghpr` (to fetch existing pr from someone's remote) and 
    `ghsendpr` could be found at [yarikoptic's bash config file](http://git.onerussian.com/?p=etc/bash.git;a=blob;f=.bash/bashrc/30_aliases_sh;hb=HEAD#l865)

3. Create a branch (generally off the `origin/master`) to hold your changes:

          git checkout -b nf-my-feature

    and start making changes. Ideally, use a prefix signaling the purpose of the
    branch
    - `nf-` for new features
    - `bf-` for bug fixes
    - `rf-` for refactoring
    - `doc-` for documentation contributions (including in the code docstrings).
    - `bm-` for changes to benchmarks
    We recommend to not work in the ``master`` branch!

4. Work on this copy on your computer using Git to do the version control. When
   you're done editing, do:

          git add modified_files
          git commit

   to record your changes in Git.  Ideally, prefix your commit messages with the
   `NF`, `BF`, `RF`, `DOC`, `BM` similar to the branch name prefixes, but you could
   also use `TST` for commits concerned solely with tests, and `BK` to signal
   that the commit causes a breakage (e.g. of tests) at that point.  Multiple
   entries could be listed joined with a `+` (e.g. `rf+doc-`).  See `git log` for
   examples.  If a commit closes an existing DataLad issue, then add to the end
   of the message `(Closes #ISSUE_NUMER)`

5. Push to GitHub with:

          git push -u gh-YourLogin nf-my-feature

   Finally, go to the web page of your fork of the DataLad repo, and click
   'Pull request' (PR) to send your changes to the maintainers for review. This
   will send an email to the committers.  You can commit new changes to this branch
   and keep pushing to your remote -- github automagically adds them to your
   previously opened PR.

(If any of the above seems like magic to you, then look up the
[Git documentation](http://git-scm.com/documentation) on the web.)

## Development environment

We support Python 3 only (>= 3.7).

See [README.md:Dependencies](README.md#Dependencies) for basic information
about installation of datalad itself.
On Debian-based systems we recommend to enable [NeuroDebian](http://neuro.debian.net)
since we use it to provide backports of recent fixed external modules we depend upon:

```sh
apt-get install -y -q git git-annex-standalone
apt-get install -y -q patool python3-scrapy python3-{argcomplete,git,humanize,keyring,lxml,msgpack,progressbar,requests,setuptools}
```

and additionally, for development we suggest to use tox and new
versions of dependencies from pypy:

```sh
apt-get install -y -q python3-{dev,httpretty,nose,pip,vcr,virtualenv} python3-tox
# Some libraries which might be needed for installing via pip
apt-get install -y -q lib{ffi,ssl,curl4-openssl,xml2,xslt1}-dev
```

some of which you could also install from PyPi using pip  (prior installation of those libraries listed above
might be necessary)

```sh
pip install -r requirements-devel.txt
```

and you will need to install recent git-annex using appropriate for your
OS means (for Debian/Ubuntu, once again, just use NeuroDebian).

Contributor Files History
-------------------------

The original repository provided a [.zenodo.json](.zenodo.json)
file, and we generate a [.contributors file](.all-contributorsrc) from that via:

```bash
pip install tributors
tributors --version
0.0.18
```

It helps to have a GitHub token to increase API limits:

```bash
export GITHUB_TOKEN=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
```

Instructions for these environment variables can be found [here](https://con.github.io/tributors/docs/getting-started#2-environment). 
Then update zenodo:

```bash
tributors update  zenodo
INFO:    zenodo:Updating .zenodo.json
INFO:    zenodo:Updating .tributors cache from .zenodo.json
WARNING:tributors:zenodo does not support updating from names.
```

In the case that there is more than one orcid found for a user, you will be given a list
to check. Others will be updated in the file. You can then curate the file as you see fit.
We next want to add the .allcontributors file:

```bash
$ tributors init allcontrib
INFO:allcontrib:Generating .all-contributorsrc for datalad/datalad
$ tributors update allcontrib
INFO:allcontrib:Updating .all-contributorsrc
INFO:allcontrib:Updating .tributors cache from .all-contributorsrc
INFO:allcontrib:⭐️ Found new contributor glalteva in .all-contributorsrc
INFO:allcontrib:⭐️ Found new contributor adswa in .all-contributorsrc
INFO:allcontrib:⭐️ Found new contributor chrhaeusler in .all-contributorsrc
...
INFO:allcontrib:⭐️ Found new contributor bpoldrack in .all-contributorsrc
INFO:allcontrib:⭐️ Found new contributor yetanothertestuser in .all-contributorsrc
WARNING:tributors:allcontrib does not support updating from orcids.
WARNING:tributors:allcontrib does not support updating from email.
```

We can then populate the shared .tributors file:

```bash
$ tributors update-lookup allcontrib
```

And then we can rely on the [GitHub action](.github/workflows/update-contributors.yml) to update contributors. The action is set to run on merges to master, meaning when the contributions are finalized. This means that we add new contributors, and we
look for new orcids as we did above.


## Documentation

### Docstrings

We use [NumPy standard] for the description of parameters docstrings.  If you are using
PyCharm, set your project settings (`Tools` -> `Python integrated tools` -> `Docstring format`).

[NumPy standard]: https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt#docstring-standard

In addition, we follow the guidelines of [Restructured Text] with the additional features and treatments
provided by [Sphinx].

[Restructured Text]: http://docutils.sourceforge.net/docs/user/rst/quickstart.html
[Sphinx]: http://www.sphinx-doc.org/en/stable/

## Additional Hints

### Merge commits

For merge commits to have more informative description, add to your
`.git/config` or `~/.gitconfig` following section:

    [merge]
    log = true

and if conflicts occur, provide short summary on how they were resolved
in "Conflicts" listing within the merge commit
(see [example](https://github.com/datalad/datalad/commit/eb062a8009d160ae51929998771964738636dcc2)).


## Quality Assurance

It is recommended to check that your contribution complies with the following
rules before submitting a pull request:

- All public methods should have informative docstrings with sample usage
  presented as doctests when appropriate.

- All other tests pass when everything is rebuilt from scratch.

- New code should be accompanied by tests.


### Tests

`datalad/tests` contains tests for the core portion of the project, and
more tests are provided under corresponding submodules in `tests/`
subdirectories to simplify re-running the tests concerning that portion
of the codebase.  To execute many tests, the codebase first needs to be
"installed" in order to generate scripts for the entry points.  For
that, the recommended course of action is to use `virtualenv`, e.g.

```sh
virtualenv --system-site-packages venv-tests
source venv-tests/bin/activate
pip install -r requirements.txt
python setup.py develop
```

and then use that virtual environment to run the tests, via

```sh
python -m nose -s -v datalad
```

or similarly,

```sh
nosetests -s -v datalad
```

then to later deactivate the virtualenv just simply enter

```sh
deactivate
```

**Note**: on Windows, please add `--traverse-namespace` option to the `nose`
call, or otherwise `nose` would not discover tests.

Alternatively, or complimentary to that, you can use `tox` -- there is a `tox.ini`
file which sets up a few virtual environments for testing locally, which you can
later reuse like any other regular virtualenv for troubleshooting.
Additionally, [tools/testing/test_README_in_docker](tools/testing/test_README_in_docker) script can
be used to establish a clean docker environment (based on any NtesteuroDebian-supported
release of Debian or Ubuntu) with all dependencies listed in README.md pre-installed.

#### Test attributes

[datalad/tests/utils.py]() defines many useful decorators. Some of those just to annotate tests
for various aspects to allow for easy sub-selection.

##### Speed

Please annotate with following decorators
- `@slow` if test runs over 10 seconds
- `@turtle` if test runs over 120 seconds (those would not typically be ran on CIs)

##### Purpose

As those tests also usually tend to be slower, use in conjunction with `@slow` or `@turtle` when slow
- `@integration` - tests verifying correct operation with external tools/services beyond git/git-annex
- `@usecase` - represents some (user) use-case, and not necessarily a "unit-test" of functionality

### CI setup

We are using Travis-CI and have [buildbot setup](https://github.com/datalad/buildbot) which also
exercises our tests battery for every PR and on the master.  Note that buildbot runs tests only submitted
by datalad developers, or if a PR acquires 'buildbot' label.

In case if you want to enter buildbot's environment

1. Login to our development server (`smaug`)

2. Find container ID associated with the environment you are interested in, e.g.

        docker ps | grep nd16.04

3. Enter that docker container environment using

        docker exec -it <CONTAINER ID> /bin/bash

4. Become buildbot user

        su - buildbot

5. Activate corresponding virtualenv using

        source <VENV/bin/activate>

   e.g. `source /home/buildbot/datalad-pr-docker-dl-nd15_04/build/venv-ci/bin/activate`

And now you should be in the same environment as the very last tested PR.
Note that the same path/venv is reused for all the PRs, so you might want
first to check using `git show` under the `build/` directory if it corresponds
to the commit you are interested to troubleshoot.

For developing on Windows you can use free [Windows VMs](https://developer.microsoft.com/en-us/microsoft-edge/tools/vms/).

### Coverage

You can also check for common programming errors with the following tools:

- Code with good unittest coverage (at least 80%), check with:

          pip install nose coverage
          nosetests --with-coverage path/to/tests_for_package

- We rely on https://codecov.io to provide convenient view of code coverage.
  Installation of the codecov extension for Firefox/Iceweasel or Chromium
  is strongly advised, since it provides coverage annotation of pull
  requests.

### Linting

We are not (yet) fully PEP8 compliant, so please use these tools as
guidelines for your contributions, but not to PEP8 entire code
base.

[beyond-pep8]: https://www.youtube.com/watch?v=wf-BqAjZb8M

*Sidenote*: watch [Raymond Hettinger - Beyond PEP 8][beyond-pep8]

- No pyflakes warnings, check with:

           pip install pyflakes
           pyflakes path/to/module.py

- No PEP8 warnings, check with:

           pip install pep8
           pep8 path/to/module.py

- AutoPEP8 can help you fix some of the easy redundant errors:

           pip install autopep8
           autopep8 path/to/pep8.py

Also, some team developers use
[PyCharm community edition](https://www.jetbrains.com/pycharm) which
provides built-in PEP8 checker and handy tools such as smart
splits/joins making it easier to maintain code following the PEP8
recommendations.  NeuroDebian provides `pycharm-community-sloppy`
package to ease pycharm installation even further.

### Benchmarking

We use [asv] to benchmark some core DataLad functionality.
The benchmarks suite is located under [benchmarks/](./benchmarks), and
periodically we publish results of running benchmarks on a dedicated host
to http://datalad.github.io/datalad/ .  Those results are collected
and available under the `.asv/` submodule of this repository, so to get started

- `git submodule update --init .asv`
- `pip install .[devel]` or just `pip install asv`
- `asv machine` - to configure asv for your host if you want to run
  benchmarks locally

And then you could use [asv] in multiple ways.

#### Quickly benchmark the working tree

- `asv run -E existing` - benchmark using the existing python environment
  and just print out results (not stored anywhere).  You can add `-q`
  to run each benchmark just once (thus less reliable estimates)
- `asv run -b api.supers.time_createadd_to_dataset -E existing`
  would run that specific benchmark using the existing python environment

Note: `--python=same` (`-E existing`) seems to have restricted
applicability, e.g. can't be used for a range of commits, so it can't
be used with `continuous`.

#### Compare results for two commits from recorded runs

Use [asv compare] to compare results from different runs, which should be
available under `.asv/results/<machine>`.  (Note that the example
below passes ref names instead of commit IDs, which requires asv v0.3
or later.)

```shell
> asv compare -m hopa maint master

All benchmarks:

       before           after         ratio
     [b619eca4]       [7635f467]
-           1.87s            1.54s     0.82  api.supers.time_createadd
-           1.85s            1.56s     0.84  api.supers.time_createadd_to_dataset
-           5.57s            4.40s     0.79  api.supers.time_installr
          145±6ms          145±6ms     1.00  api.supers.time_ls
-           4.59s            2.17s     0.47  api.supers.time_remove
          427±1ms          434±8ms     1.02  api.testds.time_create_test_dataset1
-           4.10s            3.37s     0.82  api.testds.time_create_test_dataset2x2
      1.81±0.07ms      1.73±0.04ms     0.96  core.runner.time_echo
       2.30±0.2ms      2.04±0.03ms    ~0.89  core.runner.time_echo_gitrunner
+        420±10ms          535±3ms     1.27  core.startup.time_help_np
          111±6ms          107±3ms     0.96  core.startup.time_import
+         334±6ms          466±4ms     1.39  core.startup.time_import_api
```


#### Run and compare results for two commits

[asv continuous] could be used to first run benchmarks for the to-be-tested
commits and then provide stats:

- `asv continuous maint master` - would run and compare `maint` and `master` branches
- `asv continuous HEAD` - would compare `HEAD` against `HEAD^`
- `asv continuous master HEAD` - would compare `HEAD` against state of master
- [TODO: contineous -E existing](https://github.com/airspeed-velocity/asv/issues/338#issuecomment-380520022)

Notes:
- only significant changes will be reported
- raw results from benchmarks are not stored (use `--record-samples` if
  desired)

#### Run and record benchmarks results (for later comparison etc)

- `asv run` would run all configured branches (see
  [asv.conf.json](./asv.conf.json))


#### Profile a benchmark and produce a nice graph visualization

Example (replace with the benchmark of interest)

    asv profile -v -o profile.gprof usecases.study_forrest.time_make_studyforrest_mockup
    gprof2dot -f pstats profile.gprof | dot -Tpng -o profile.png \
        && xdg-open profile.png

#### Common options

- `-E` to restrict to specific environment, e.g. `-E virtualenv:2.7`
- `-b` could be used to specify specific benchmark(s)
- `-q` to run benchmark just once for a quick assessment (results are
  not stored since too unreliable)


[asv compare]: http://asv.readthedocs.io/en/latest/commands.html#asv-compare
[asv continuous]: http://asv.readthedocs.io/en/latest/commands.html#asv-continuous

[asv]: http://asv.readthedocs.io


## Easy Issues

A great way to start contributing to DataLad is to pick an item from the list of
[Easy issues](https://github.com/datalad/datalad/labels/easy) in the issue
tracker.  Resolving these issues allows you to start contributing to the project
without much prior knowledge.  Your assistance in this area will be greatly
appreciated by the more experienced developers as it helps free up their time to
concentrate on other issues.

## Recognizing contributions

We welcome and recognize all contributions from documentation to testing to code development.

You can see a list of current contributors in our [zenodo file][link_zenodo].
If you are new to the project, don't forget to add your name and affiliation there!
We also have an .all-contributorsrc that is updated automatically on merges. Once it's
merged, if you helped in a non standard way (e.g., a contribution other than code)
you can open a pull request to add any [All Contributors Emoji][contrib_emoji] that
match your contribution types.

## Thank you!

You're awesome. :wave::smiley:



# Various hints for developers

## Useful tools

- While performing IO/net heavy operations use [dstat](http://dag.wieers.com/home-made/dstat)
  for quick logging of various health stats in a separate terminal window:
  
        dstat -c --top-cpu -d --top-bio --top-latency --net

- To monitor speed of any data pipelining [pv](http://www.ivarch.com/programs/pv.shtml) is really handy,
  just plug it in the middle of your pipe.

- For remote debugging epdb could be used (avail in pip) by using
  `import epdb; epdb.serve()` in Python code and then connecting to it with
  `python -c "import epdb; epdb.connect()".`

- We are using codecov which has extensions for the popular browsers
  (Firefox, Chrome) which annotates pull requests on github regarding changed coverage.

## Useful Environment Variables

Refer datalad/config.py for information on how to add these environment variables to the config file and their naming convention

- *DATALAD_DATASETS_TOPURL*:
  Used to point to an alternative location for `///` dataset. If running
  tests preferred to be set to https://datasets-tests.datalad.org
- *DATALAD_LOG_LEVEL*:
  Used for control the verbosity of logs printed to stdout while running datalad commands/debugging
- *DATALAD_LOG_NAME*:
  Whether to include logger name (e.g. `datalad.support.sshconnector`) in the log
- *DATALAD_LOG_OUTPUTS*:
  Used to control either both stdout and stderr of external commands execution are logged in detail (at DEBUG level)
- *DATALAD_LOG_PID*
  To instruct datalad to log PID of the process
- *DATALAD_LOG_TARGET*
  Where to log: `stderr` (default), `stdout`, or another filename
- *DATALAD_LOG_TIMESTAMP*:
  Used to add timestamp to datalad logs
- *DATALAD_LOG_TRACEBACK*:
  Runs TraceBack function with collide set to True, if this flag is set to 'collide'.
  This replaces any common prefix between current traceback log and previous invocation with "..."
- *DATALAD_LOG_VMEM*:
  Reports memory utilization (resident/virtual) at every log line, needs `psutil` module
- *DATALAD_EXC_STR_TBLIMIT*: 
  This flag is used by datalad to cap the number of traceback steps included in exception logging and result reporting to DATALAD_EXC_STR_TBLIMIT of pre-processed entries from traceback.
- *DATALAD_SEED*:
  To seed Python's `random` RNG, which will also be used for generation of dataset UUIDs to make
  those random values reproducible.  You might want also to set all the relevant git config variables
  like we do in one of the travis runs
- *DATALAD_TESTS_TEMP_KEEP*: 
  Function rmtemp will not remove temporary file/directory created for testing if this flag is set
- *DATALAD_TESTS_TEMP_DIR*: 
  Create a temporary directory at location specified by this flag.
  It is used by tests to create a temporary git directory while testing git annex archives etc
- *DATALAD_TESTS_NONETWORK*: 
  Skips network tests completely if this flag is set
  Examples include test for s3, git_repositories, openfmri etc
- *DATALAD_TESTS_SSH*: 
  Skips SSH tests if this flag is **not** set.  If you enable this,
  you need to set up a "datalad-test" and "datalad-test2" target in
  your SSH configuration.  The second target is used by only a couple
  of tests, so depending on the tests you're interested in, you can
  get by with only "datalad-test" configured.

  A Docker image that is used for DataLad's tests is available at
  <https://github.com/datalad-tester/docker-ssh-target>.  Note that
  the DataLad tests assume that target files exist in
  `DATALAD_TESTS_TEMP_DIR`, which restricts the "datalad-test" target
  to being either the localhost or a container that mounts
  `DATALAD_TESTS_TEMP_DIR`.
- *DATALAD_TESTS_NOTEARDOWN*: 
  Does not execute teardown_package which cleans up temp files and directories created by tests if this flag is set
- *DATALAD_TESTS_USECASSETTE*:
  Specifies the location of the file to record network transactions by the VCR module.
  Currently used by when testing custom special remotes
- *DATALAD_TESTS_OBSCURE_PREFIX*:
  A string to prefix the most obscure (but supported by the filesystem test filename
- *DATALAD_TESTS_PROTOCOLREMOTE*:
  Binary flag to specify whether to test protocol interactions of custom remote with annex
- *DATALAD_TESTS_RUNCMDLINE*:
  Binary flag to specify if shell testing using shunit2 to be carried out
- *DATALAD_TESTS_TEMP_FS*:
  Specify the temporary file system to use as loop device for testing DATALAD_TESTS_TEMP_DIR creation
- *DATALAD_TESTS_TEMP_FSSIZE*:
  Specify the size of temporary file system to use as loop device for testing DATALAD_TESTS_TEMP_DIR creation
- *DATALAD_TESTS_NONLO*:
  Specifies network interfaces to bring down/up for testing. Currently used by travis.
- *DATALAD_TESTS_KNOWNFAILURES_PROBE*:
  Binary flag to test whether "known failures" still actually are failures. That
  is - change behavior of tests, that decorated with any of the `known_failure`,
  to not skip, but executed and *fail* if they would pass (indicating that the
  decorator may be removed/reconsidered).
- *DATALAD_CMD_PROTOCOL*:
  Specifies the protocol number used by the Runner to note shell command or python function call times and allows for dry runs.
  'externals-time' for ExecutionTimeExternalsProtocol, 'time' for ExecutionTimeProtocol and 'null' for NullProtocol.
  Any new DATALAD_CMD_PROTOCOL has to implement datalad.support.protocol.ProtocolInterface
- *DATALAD_CMD_PROTOCOL_PREFIX*:
  Sets a prefix to add before the command call times are noted by DATALAD_CMD_PROTOCOL.
- *DATALAD_USE_DEFAULT_GIT*:
  Instructs to use `git` as available in current environment, and not the one which possibly comes with git-annex (default behavior).
- *DATALAD_ASSERT_NO_OPEN_FILES*:
  Instructs test helpers to check for open files at the end of a test. If set, remaining open files are logged at ERROR level. Alternative modes are: "assert" (raise AssertionError if any open file is found), "pdb"/"epdb" (drop into debugger when open files are found, info on files is provided in a "files" dictionary, mapping filenames to psutil process objects).
- *DATALAD_ALLOW_FAIL*:
  Instructs `@never_fail` decorator to allow to fail, e.g. to ease debugging.

# Release(s) workflow

## Branches

- `master`: changes toward the next `MAJOR.MINOR.0` release.
  Release candidates (tagged with an `rcX` suffix) are cut from this branch
- `maint`: bug fixes for the latest released `MAJOR.MINOR.PATCH`
- `maint-MAJOR.MINOR`: generally not used, unless some bug fix release with a critical bug fix is needed.

## Workflow

- upon release of `MAJOR.MINOR.0`, `maint` branch needs to be fast-forwarded to that release
- bug fixes to functionality released within the `maint` branch should be
  submitted against `maint` branch
- cherry-picking fixes from `master` into `maint` is allowed where needed
- `master` branch accepts PRs with new functionality
- `master` branch merges `maint` as frequently as needed

## Helpers

[Makefile](./Makefile) provides a number of useful `make` targets:

- `linkissues-changelog`: converts `(#ISSUE)` placeholders into proper markdown within [CHANGELOG.md]()
- `update-changelog`: uses above `linkissues-changelog` and updates .rst changelog
- `release-pypi`: ensures no `dist/` exists yet, creates a wheel and a source distribution and uploads to pypi.

## Releasing with GitHub Actions, auto, and pull requests

New releases of datalad are created via a GitHub Actions workflow built
around [`auto`](https://github.com/intuit/auto).  Whenever a pull request is
merged into `maint` that has the "`release`" label, `auto` updates the
changelog based on the pull requests since the last release, commits the
results, tags the new commit with the next version number, and creates a GitHub
release for the tag.  This in turn triggers a job for building an sdist & wheel
for the project and uploading them to PyPI.

### Labelling pull requests

The section that `auto` adds to the changelog on a new release consists of the
titles of all pull requests merged into master since the previous release,
organized by label.  `auto` recognizes the following PR labels:

- `semver-minor` — for changes corresponding to an increase in the minor version
  component
- `semver-patch` — for changes corresponding to an increase in the patch/micro version
  component; this is the default label for unlabelled PRs
- `semver-internal` — for changes only affecting the internal API
- `semver-documentation` — for changes only affecting the documentation
- `semver-tests` — for changes to tests
- `semver-dependencies` — for updates to dependency versions
- `semver-performance` — for performance improvements

[link_zenodo]: https://github.com/datalad/datalad/blob/master/.zenodo.json
[contrib_emoji]: https://allcontributors.org/docs/en/emoji-key
### Instructions

Please go through the following checklist before submitting the PR:

- [ ] Provide an overview of the changes you're making and explain why you're proposing them.

- [ ] Include `Fixes #NNN` or `Closes #NNN` somewhere in the description if this PR addresses an existing issue.

- [ ] If this PR is not complete, select "Create Draft Pull Request" in the pull request button's menu.

  Consider using a task list (e.g., `- [ ] add tests ...`) to indicate remaining to-do items.

- [ ] If you would like to list yourself as a DataLad contributor and your name is not mentioned please modify .zenodo.json file.

- [ ] **Delete these instructions**. :-)

Thanks for contributing!
Thoughts about redesign, well actually "design" since originally there
were none, of datalad crawl.

Global portion of the config
============================

::

  [dataset]
  path =
  description = 
  exec = 

Data providers
==============

`crawl` command collects data present possibly across
different remote data providers (regular HTTP websites, AWS S3
buckets, etc) and then consolidates access to them within a single
git-annex'ed repository.  `crawl` should also keep track of
status/versions of the files, so in case of the updates (changes,
removals, etc) on remote sites, git-annex repository could be
correspondingly updated.

Common config specs::

    [provider:XXX]
    type = (web|s3http|merge|git-annex) # default to web
    branch = master              # default to master
    commit_to_git =              # regexps of file names to commit directly to git
    ignore =                     # files to ignore entirely
    drop = False                 # either to drop the load upon 'completion'
    # some sanity checks
    check_entries_limit = -1     # no limit by default


(To be) Supported data providers
--------------------------------

Web
~~~

In many usecases data are hosted on a public portal, lab website,
personal page, etc.  Such data are often provided in tarballs, which
need to be downloaded and extracted later on.  Extraction will not be
a part of this provider -- only download from the web::

    [provider:incoming_http]
    type = web
    mode = (download|fast|relaxed)            # fast/relaxed/download
    filename = (url|request)                  # of cause also could be _e'valuated given the bs4 link get_video_filename(link, filename)
    recurse_(a|href) =                        # regexes to recurse
    # mimicing scrapy
    start_urls = http://...                   #
    certificates =                            # if there are https -- we need to allow specifying those
    allowed_domains = example.com/buga/duga   # to limit recursion
                      sample.com
    excluded_hrefs =                          # do not even search for "download" URLs on given pages.  Should also allow to be a function/callback to decide based on request?
    include_(a|href) =                        # what to download
    exclude_(a|href) =                        # and not (even if matches)
    ???generators = generate_readme              # Define some additional actions to be performed....

We need to separate options for crawling (recursion etc) and deciding
what to download/annex.

Q: should we just specify xpath's for information to get extracted
   from a response corresponding to a matching url?  just any crawled page?

Q: allow to use xpath syntax for better control of what to recurse/include?

Q: authentication -- we should here relate to the Hostings
!: scrapy's Spider provides start_requests() which could be used to
   initiate the connection, e.g. to authenticate and then use that connection.
   Authentication detail must not be a part of the configuration, BUT
   it must know HOW authentication should be achieved.  In many cases could
   be a regular netrc-style support (so username/password).

   Those authenticators should later be reused by "download clients"

Q: we might need to worry/test virtually about every possible associated
   to http downloads scenario, e.g. support proxy (with authentication).
   May be we could just switch to aria2 and allow to specify access options?

Q: may be (a new provider?) allow to use a scrapy spider's output to
   harvest the table of links which need to be fetched



Use cases to keep in mind
+++++++++++++++++++++++++

- versioning present in the file names
  ftp://ftp.ncbi.nlm.nih.gov/1000genomes/ftp/sequence_indices/

  - ha -- idea, all those should be referred in some other branch, like
    with archives, and then 'public' one would just take care about
    pointing to the "correct one" and serve a "compressed" view.
    Hence: monitor original, point "compressed" to a branch giving it
    a set of rules on how to determine version, i.e. on which files
    This way we could have both referenced in the same repository.


Amazon S3
~~~~~~~~~

Initial accent will be made on S3 buckets which have versioning
enabled, and which expose their content via regular http/https.

tricky points:
- versioning (must be enabled. If uploaded before enabled, version is Null)

- etags are md5s BUT only if upload was not multi-chunked, so
  it becomes difficult to identify files by md5sums (must be
  downloaded first then, or some meta-info of file should be modified so
  etag gets regenerated -- should result in file md5sum appearing as etag)

- S3 most probably would be just an additional -- not the primary provider


Generated
~~~~~~~~~

We should allow for files to be generated based on the content of the
repository and/or original information from the data providers,
e.g. content of the webpages containing the files to be
downloaded/referenced.  Originally envisioned as a separate branch,
where only archived content would be downloaded and later extracted
into corresponding locations of the "public" branch (e.g. master).

But may be it should be more similar to the stated above "versioning"
idea where it would simply be an alternative "view" of another branch,
where some content is simply extracted.  I.e. all those modifications
could be assembled as a set of "filters"::

    [generator:generate_readme]
    filename = README.txt
    content_e = generate_readme(link, filename)  # those should be obtained/provided while crawling

or

    [generator:fetch_license]
    filename = LICENSE.txt
    content_e = fetch_license(link, filename)  # those should be obtained/provided while crawling


Merge
~~~~~

Originally fetched Files might reside in e.g. 'incoming' branch while
'master' branch then could be 'assembled' from few other branches with
help of filtering::

    [provider:master]
    type = merge
    branch = master # here matches the name but see below if we need to repeat
    merge = incoming_data_http
            incoming_meta_website
    filters = extract_models
              extract_data
              generate_some_more_files_if_you_like


Q: should we may be 'git merge --no-commit' and then apply the
   filters???

   probably not since there could be conflicts if similarly named file
   is present in target branch (e.g. generated) and was present
   (moved/renamed via filters) in the original branch.

Q: but merging of branches is way too cool and better establishes the
   'timeline' and dependencies...
   So merge should be done "manually" by doing (there must be cleaner way)::

     git merge -s ours --no-commit
     git rm -r *
     # - collect and copy files for all the File's from branches to .
     # - stage all the files
     # - pipe those "File"s from all the branches through the filters
     #   (those should where necessary use git rm, mv, etc)
     # - add those File's to git/git-annex
     git commit

   but what if a filter (e.g. cmd) requires current state of files from
   different branches?...  all possible conflict problems could be
   mitigated by storing content in branches under some directories,
   then manipulating upon "merge" and renaming before actually 'git merging'


Q: what about filters per incoming branch???  we could options for
   filters specification
   (e.g. extract_models[branches=incoming_data_http]) or allow
   only regular 2-edge merge at a time but multiple times...


XNAT, COINS, ...
~~~~~~~~~~~~~~~~

Later ... but the idea should be the same I guess: they should expose
collections of File's with a set of URIs so they could be addurl'ed to
the files.  It is not clear yet either they would need to be crawled
or would provide some API similar to S3 to request all the necessary
information?


git/git-annex
~~~~~~~~~~~~~

If provider is already a Git(-annex) repository.  Usecase:
forrest_gump.  So it is pretty much a regular remote **but** it might
benefit from our filters etc.


torrent
~~~~~~~

I guess similar/identical to archives if torrent points to a single
file -- so just 'addurl'.  If torrent provides multiple files, would
need mapping of UUIDs I guess back to torrents/corresponding files.
So again -- similar to archives...?

aria2 seems to provide a single unified HTTP/HTTPS/FTP/BitTorrent
support, with fancy simultaneous fetching from multiple
remotes/feeding back to the torrent swarm (caution for non-free data).
It also has RPC support, which seems to be quite cool and might come
handy (e.g. to monitor progress etc)


Wild: Git repository for being rewritten
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

theoretically we could collect all the information to rewrite some
other Git repo but now injecting some files into git-annex (while
possibly even pointing for the load to e.g. original SVN repo).

Tricky:
- branches and merges -- would be really tricky and so far not
  envisioned how
- "updates" should correspond to commits in original repository
- all the commit information should be extracted/provided for the
  commit here


Filters
=======

Considering idea that all the modifications (archives extraction,
versioning etc) could be made through monitoring of another branch(es)
and applying a set of filters.

- files which aren't modified, should also propagate into target
  branch, along with all their urls

  file by file wouldn't work since filter might need to analyze the
  entire list of files...::

      def apply_filters(self):
       files_out = files_in
       for filter in self.filters:
        files_out = filter.apply(files_out)
       return files_out

  then each filter would decide on how to treat the list of files.
  May be some filters' subtyping would be desired
  (PerfileFilter/AllfilesFilter)

- filters should provide API to 'rerun' their action to obtain the
  same result.


Cross-branch
------------

Some filters to be applied on files from one branch to have results
placed into another:


Extract
~~~~~~~

Special kind of a beast: while keeping the original archive under
git-annex obtained from any other provider (e.g. 'Web'), we extract
the load (possibly with some filtering/selection):

  Q: how to deal with extract from archives -- extraction should
     better be queued to extract multiple files from the archive at
     once.  But ATM it would not happen since all those URIs will
     simply be requested by simple wget/curl calls by git-annex file
     at a time.
  A: upon such a first call, check if there is .../extracted_key/key/, if
     there is -- use.  If not -- extract and then use. use = hardlink
     into the target file.
     Upon completion of `datalad get` (or some other command) verify
     that all `/extracted/` are removed (and/or provide setting -- may
     be we could/should just keep those around)


Config Examples:
++++++++++++++++

::

    [filter:extract_models]
    filter = extract               # by default would be taken as the element after "filter:"
    input = *(\S+)_models\.tgz$    # and those files are not provided into output
    output_prefix = models/$1/     # somehow we should allow to reuse input regex's groups
    exclude =                      # regex for files to be excluded from extraction or straight for tar?
    strip_path = 1

Probably will just use patoolib (do not remember if has
strip_path... seems not:
https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=757483)

URI:  dl:extract:UID

and we keep information for what 'key' it came into what file (which
might later get renamed, so extraction from the archive shouldn't
later happen in-place, but rather outside and then moved accordingly)

Tricky point(s):

- may be by default should still extract all known archives types and
  just rely on the filename logic?
- the same file might be available from multiple archives.
  So we would need to keep track from previous updates, from which
  archive files could be fetched.
  - how to remove if archive is no longer avail?
    probably some fsck should take care about checking if archives
    are still avail, and if not -- remove the url

- keep track which files came from the archive, so we could later
  remove them happen if archive misses the file now.

Q: allow for 'relaxed' handling?
   If tarballs are not versioned at all, but we would like to create
   overall (? or just per files) 'relaxed' git-annex?

   Probably no complication if URIs will be based (natively) on the
   fast or relaxed keys.  Sure thing things would fail if archive was
   changed and lacks the file.

Q: hm -- what about MD5SUM checking? e.g. if archive was posted with
   the MD5SUMs file

   I guess some kind of additional filter which could be attached
   somehow?


Move/Rename/Delete
~~~~~~~~~~~~~~~~~~

Just move/rename/delete some files around e.g. for a custom view of
the dataset (e.g. to conform openfmri layout). Key would simply be
reused ;)

Q: should it be 'Within-branch' filter?


Command
~~~~~~~

A universal filter which would operate on some files and output
possibly in place or modified ones...

Then it would need to harvest and encode into file's URI the
provenance -- i.e. so it could later be recreated automagically.

For simple usecases (e.g. creation of lateralized atlas in HOX, some
data curation, etc)

URI:  dl:cmd:UID

while we keep a file providing the corresponding command for each UID,
where ARGUMENTS will would point to the original files keys in the git
annex.   Should it be kept in PROV format may be???

Config Examples::

    [filter:command_gunzip]
    in1 = *\.gz
    in2_e = in1.replace('.gz', '')
    #eval_order=in1 in2
    command = zcat {in1} > {in2}
    output_files = {in2}

Problems:

- might be tricky to provide generic enough interface?
- we need plentiful of use-cases to get it right, so this one is just
  to keep in mind for future -- might be quite cool after all.


Within-branch
-------------

Other "Filters" should operate within the branch, primarily simply for
checking the content


Checksum
~~~~~~~~

e.g. point to MD5SUMS file stored in the branch, provide how file
names must be augmented, run verification -- no files output, just the
status

Addurl
~~~~~~

If the repository is going/was published also online under some URL.
We might like to populate files with corresponding urls.

    [filter:addurl]
    prefix = http://psydata.ovgu.de/forrest_gump/.git/annex/
    check = (False|True)  # to verify presence or not ???

Usecase -- Michael's forrest_gump repository.  Now files are not
associated explicitly with that URL -- only via a regular git remote.
This cumbersomes work with clones which then all must have original
repository added as a remote.

`check = False` could be the one needed for a 'publish' operation
where this data present locally is not yet published anywhere.

Tagging
~~~~~~~

We might like to tag files... TODO: think what to provide/use to
develop nice tags.

Ideas:

- a tag given a set of filename regexps

      [tag:anatomicals]
      files = .*\_anat\.nii\.gz
      tag = modality=anatomy

   or just

      [tag:anatomicals]
      files = .*\_anat\.nii\.gz
      tag = anatomy

   if it is just a tag (anatomy) without a field

 - (full)filename regexp with groups defining possibly multiple
   tag/value pairs

      [tag:modality]
      files = .*\_(?P<modality>\S*)\.nii\.gz
      translate = anat: T1     #  might need some translation dictionary?
                  dwi: DTI


Notes:
- metadata cane be added only to files under git-annex control so those
  directly committed

Design thoughts
===============

Data providers should provide a unified interface

DataProvider
~~~~~~~~~~~~

Common Parameters
- add_to_git  - what files to commit to git directly (should we leverage
   git-annex largefiles option somehow?)
- ignore      - what files to ignore

- get_items(version=None) - return a list of Files
- get_item_by_name
- get_item_by_md5
  - should those be additional interfaces?
  - what if multiple items fulfill (content is the same, e.g. empty, names differ,
    we better get the most appropriate in the name or don't give a damn?)
  - what if a collision????
- get_item_by_sha256
  - e.g. natively provided by 'Branch' provider for annexed files
    (what to do about git committed ones -- compute/keep info?)
- get_versions(min_version=None)
  provider-wide version (i.e. not per file).  E.g. S3
  provider can have multiple versions of files.
  Might be that it needs to return a DAG of versions i.e. a
  (version, [prev_version1, prev_version2, ...]) to represent e.g.
  history of a Git repo.  In most of the cases would be degenerate to just
  one prev version, in which case could just be (version, ).
  We would need to store that meta-information for future updates at least
  for the last version so we could 'grow' next ones on top.
- ? get_release_versions() -- by default identical to above... but might
  differ (update was, but no new official release (yet), so no release
  tag)
- get_version_metainformation() -- primarily conceived when thinking
  about monitoring other VCS repos... so should be information to be
  used for a new Git commit into this new repository

.. note:

    Keep in mind
    - Web DataProvider must have an option to request the content filename
      (addressing use case with redirects etc)
    - Some providers might have multiple URIs (mirrors) so right away
      assign them per each file...  As such they might be from
      different Hostings!


File
~~~~

what would be saved as a file.  Should know about itself... and origins!

- filename
- URIs  - list containing origins (e.g. URLs) on where to fetch it from.
          First provided by the
          original DataProvider, but then might be expanded using
          other DataProviders
          Q: Those might need to be not just URIs but some classes associated
          with original Hosting's, e.g. for the cases of authentication etc?
          or we would associate with a Hosting based on the URI?
  # combination of known fields should be stored/used to detect changes
  # Different data providers might rely on a different subset of below
  # to see if there was a change.  We should probably assume some
  # "correspondence"
- key   # was thinking about Branch as DataProvider -- those must be reused
- md5
- sha256
- mtime
- size

It will be the job of a DataProvider to initiate File with the
appropriate filename.

URI
~~~

-> URL(URI):  will be our first and main "target" but it could
              also be direct S3, etc.

a URI should be associated with an "Hosting" (many-to-one), so we could
e.g. provide authentication information per actual "Hosting" as the
entity.  But now we are getting back to DataProvider, which is the
Hosting, or actually also a part of it (since Hosting could serve
multiple Providers, e.g. openfmri -> providers per each dataset?)
But also Provider might use/point to multiple Hostings (e.g. mirrors
listed on nitp-2013).

Hosting
~~~~~~~

Each DataProvider would be a factory of File's.


Ideas to not forget
~~~~~~~~~~~~~~~~~~~

- Before carrying out some operation, remember the state of all
  (involved) branches, so it would be very easy later on to "cancel"
  the entire transaction through a set of 'git reset --hard' or
  'update-ref's.

  Keep log of the above!

- multiple data providers could be specified but there should be
  'primary' and 'complimentary' ones:

  - primary provider(s) define the layout/content
  - complimentary providers just provide references to additional
    locations where that data (uniquely identified via checksums etc)
    could be obtained, so we could add more data providing urls
  - Q: should all DataProvider's be able to serve as primary and complimentary?
  - most probably we should allow for an option to 'fail' or issue a
    warning in some cases
    - secondary provider doesn't carry a requested load/file
    - secondary provider provides some files not provided by the primary
      data provider

- at the end of the crawl operation, verify that all the files have all
  and only urls from the provided data providers

- allow to add/specify conventional git/annex clones as additional,
  conventional (non special) remotes to be added.

- allow to prepopulate URLs given e.g. perspective hosting on HTTP.
  This way whenever content gets published there -- all files would
  have appropriate URLs associated and would 'transcend' through the
  clones without requiring adding original remote.

Updates
=======

- must track updates and removals of the files
- must verify presence (add, remove) of the urls associated with the
  files given a list of data providers


Meta information
================

Since a while `git annex` provides a neat feature allowing to assign
tags to the files and later use e.g. `git annex view` to quickly
generate customized views of the repository.


Some cool related tools
=======================

https://github.com/scrapy/scrapely
  Pure Python (no DOM, lxml, etc) scraping of pages, "training" the
  scraper given a sample.  May be could be handy???
https://github.com/scrapy/slybot
  Brings together scrapy + scrapely to provide json-specs for
  spiders/items/etc
  Might be worth at least adopting spiders specs...?
  Has a neat slybot/validation/schemas.json  which validates the schematic 
datalad metadata
================

**This was somehwat cool in the past -- ignore in the present. Kept for future references**

This is a documentation on datalad's approach to metadata. Especially on how
the metadata representation currently looks like.

DataLad uses RDF to represent metadata. However, this kind of representation is
required by datalad within collections only. A dataset may or may not contain
metadata, which is prepared that way. A collection's metadata about a dataset
can be imported from any location (within or not within the dataset itself) and
various metadata formats (rdf as well as non-rdf). This allows for different
collections containing the very same dataset but different metadata. It also
means, that any git-annex repository can be a dataset contained in a collection
without the need to be touched by datalad before.

dataset metadata
----------------

The metadata of a dataset has two levels. The first one contains the metadata
about the actual content of a dataset and is provided by whoever is creating or
maintaining the dataset. For this purpose, datalad is able to import metadata
from different metadata formats and represent this metadata as RDF statements.
There are a number of things datalad expects to be expressed by the use of
certain terms. In case of a non-rdf format datalad will generate statements,
that use these terms and in case of a rdf format already provided, datalad will
add statements using these terms while keeping the originally used ones in
order to provide the opportunity to use both in queries.

The second level is metadata about the dataset itself. This is generated by
datalad.

This sums up to a set of statements datalad expects to be present in the
metadata or has to generate respectively. I'll call this set of statements the
"datalad dataset descriptor".

Note: To be clear - "expects" means: If the information is available it is
provided by using this terms. It doesn't mean that certain information
necessarily is available, nor does it mean, that these information aren't
provided using other terms, too.

datalad dataset descriptor
~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the set of statements currently considered to be the datalad dataset
descriptor.
Note: There may be some minor changes or extensions soon.

Used prefixes::

rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
rdfs: <http://www.w3.org/2000/01/rdf-schema#>
xsd: <http://www.w3.org/2001/XMLSchema#>
prov: <http://www.w3.org/ns/prov#>
dcat: <http://www.w3.org/ns/dcat#>
dctypes: <http://purl.org/dc/dcmitype/>
dct: <http://purl.org/dc/terms/>
pav: <http://purl.org/pav/>
foaf: <http://xmlns.com/foaf/0.1/>
dlns: <http://datalad.org/terms/>
: <>

Using RDF we are talking about describing resources. The dataset is a
repository, so the resource is the path to (or the URL of) the repository.
Therefore, the first statement is to declare this resource is a datalad dataset::

<path/to/dataset> a dlns:Handle .

However, in case of a dataset descriptor, that is stored within the dataset itself,
this cannot be done this way for various reasons. In that case we use the
'special resource' dlns:this instead. When imported to a collection, it is
replaced by the URI, the collection uses to point to the dataset anyway.

To identify resources within the context of that dataset, like persons that play
a role described in the metadata, we use an 'empty prefix'::

@prefix : <>

Note: This will most likely slightly change, since it doesn't behave as
expected with rdflib.

So, we can now state who created the dataset. Note: We may call this the
"author" of a dataset, but this is not necessarily an author of the actual
content of the dataset::

:someone a prov:Person, foaf:Person ;
    foaf:name "someone"^^xsd:string;
.

:datalad a prov:SoftwareAgent ;
    rdfs:label "datalad"^^xsd:string;
    pav:version "1.0a"^^xsd:string;
.


<path/to/dataset> a dlns:Handle ;
    pav:createdBy :someone ;
    pav:createdWith :datalad ;

Additionally, a dataset has a description, a title and a license::

<path/to/dataset> a dlns:Handle ;
    pav:createdBy :someone ;
    pav:createdWith :datalad ;
    dct:title "the dataset's name"^^xsd:string;
    dct:description """This is a dataset and therefore it contains
    some kind of data. Probably the data is about a certain topic and was
    generated somehow.""" ;
    dct:license <uri/of/the/license> ;
.

Now, the dataset has some content, which can be described by different types of
data entities. There are a lot of terms that may be used to classify these
entities. That's not our concern. We just state, that the dataset contains these
entities and that these entities were authored by some people, so we can query
the metadata for that information. That very information may be stated using
other terms already (see 'content2' below). In that case we keep that statement,
but our own::

:content1 a dctypes:Dataset ;
    pav:authoredBy :someauthor ;
    pav:authoredBy :someotherauthor ;
.

:content2 a dcat:Distribution ;
    anotherNamespace:creator :someauthor ;
    pav:authoredBy :someauthor ;
.

<path/to/dataset> a dlns:Handle ;
    pav:createdBy :someone ;
    ... see above ...
    dct:hasPart :content1 ;
    dct:hasPart :content2 ;
.


In case the content's metadata doesn't provide data entities using certain terms
already, we create one data entity of type 'dctypes:Dataset' to describe the
content of the dataset.

# TODO reminders:

collection metadata
-------------------

(TODO)
(very similar)
dct:hasPart => dataset


datalad config data
-------------------

dlns:usesSource
Config file format
==================

Includes
--------

:mod:`~datalad.config` enhances class:`~ConfigParser.SafeConfigParser`
to also support includes.  It makes it possible to specify within the
`INCLUDES` section of the config file which other config file should
be considered before or after a currently read config file::

    [INCLUDES]
    before = defaults.cfg
    after = customizations.cfg

    [DEFAULTS]
    ....

Sections
--------

Download section
~~~~~~~~~~~~~~~~

It is the only type of a section at this point.  It specifies a single
resource which crawl/monitor and fetch specified content to be
deposited into the git/git-annex repository.  Following fields are
known and could either be specified in the specific section or in
`DEFAULT` section to be reused across different sections

mode
  Could be `download`, `fast` or `relaxed`. In `download` mode files
  are downloaded, and added to the git-annex, thus based on a checksum
  backend.  `fast` and `relaxed` modes correspond to the modes of `git
  annex addurl`
incoming
  Path to the `incoming` repository -- where everything gets initially
  imported, e.g. original archives.  If no archives to be extracted,
  it usually then matches with `public`.  Original idea for such a
  separation was to cover the cases where incoming materials
  (archives) might contain some non-distributable materials which
  should be stripped before being placed into `public` repository
public
  Path to the `public` repository which is the target repository to be
  shared
description
  Textual description to be placed under :file:`.git/description`
include_href
  Regular expression to specify which URLs, pointed to by HTML `<A>`
  should be considered to be added to the repository
include_href_a
  Regular expression to specify which links with matching text should
  be considered
exclude_href, exclude_href_a
  Similar setups to specify which ones to exclude (e.g. if `include_href=.*`)
recurse
  Regular expression to specify which URLs to consider for further
  traversal while crawling the website

TODO. Some additional documentation is currently within
:meth:`datalad.config.EnhancedConfigParser.get_default`
Delineation from related solutions
**********************************

To our knowledge, there is no other effort with a scope as broad as DataLad's.
DataLad aims to unify access to vast arrays of (scientific) data in a domain and
data modality agnostic fashion with as few and universally available software
dependencies as possible.

The most comparable project regarding the idea of federating access to various
data providers is the iRODS_-based `INCF Dataspace`_ project.  IRODS is a
powerful, NSF-supported framework, but it requires non-trivial deployment and
management procedures. As a representative of *data grid* technology, it is
more suitable for an institutional deployment, as data access, authentication,
permission management, and versioning are complex and not-feasible to be
performed directly by researchers. DataLad on the other hand federates
institutionally hosted data, but in addition enables individual researchers and
small labs to contribute datasets to the federation with minimal cost and
without the need for centralized coordination and permission management.

.. _IRODS: https://irods.org
.. _INCF Dataspace: http://www.incf.org/resources/data-space


Data catalogs
=============

Existing data-portals, such as DataDryad_, or domain-specific ones (e.g. `Human
Connectome`_, OpenfMRI_), concentrate on collecting, cataloging, and making
data available. They offer an abstraction from local data management
peculiarities (organization, updates, sharing).  Ad-hoc collections of pointers
to available data, such as `reddit datasets`_ and `Inside-R datasets`_, do not
provide any unified interface to assemble and manage such data.  Data portals
can be used as seed information and data providers for DataLad. These portals
could in turn adopt DataLad to expose readily usable data collections via a
federated infrastructure.

.. _Human Connectome: http://www.humanconnectomeproject.org
.. _OpenfMRI: http://openfmri.org
.. _DataDryad: http://datadryad.org
.. _reddit datasets: http://www.reddit.com/r/datasets
.. _Inside-R datasets: http://www.inside-r.org/howto/finding-data-internet


Data delivery/management middleware
===================================

Even though there are projects to manage data directly with dVCS (e.g. Git),
such as the `Rdatasets Git repository`_ this approach does not scale, for example
to the amount of data typically observed in a scientific context. DataLad
uses git-annex_ to support managing large amounts of data with Git, while
avoiding the scalability issues of putting data directly into Git repositories.

In scientific software development, frequently using Git for source code
management, many projects are also confronted with the problem of managing
large data arrays needed, for example, for software testing. An exemplar
project is `ITK Data`_ which is conceptually similar to git-annex: data content
is referenced by unique keys (checksums), which are made redundantly available
through multiple remote key-store farms and can be obtained using specialized
functionality in the CMake software build system.  However, the scope of this
project is limited to software QA, and only provides an ad-hoc collection of
guidelines and supporting scripts.

.. _Rdatasets Git repository: http://github.com/vincentarelbundock/Rdatasets
.. _ITK Data: http://www.itk.org/Wiki/ITK/Git/Develop/Data

The git-annex website provides a comparison_ of Git-annex to other available
distributed data management tools, such as git-media_, git-fat_, and others.
None of the alternative frameworks provides all of the features of git-annex,
such as integration with native Git workflows, distributed redundant storage,
and partial checkouts in one project.  Additional features of git-annex which
are not necessarily needed by DataLad (git-annex assistant, encryption support,
etc.) make it even more appealing for extended coverage of numerous scenarios.
Moreover, neither of the alternative solutions has already reached a maturity,
availability, and level of adoption that would be comparable to that of
git-annex.

.. _git-annex: http://git-annex.branchable.com
.. _comparison: http://git-annex.branchable.com/not
.. _git-media: https://github.com/schacon/git-media
.. _git-fat: https://github.com/jedbrown/git-fat

.. _chap-git-annex-datalad-comparison:

Git/Git-annex/DataLad
=====================

Although it is possible, and intended, to use DataLad without ever invoking git
or git-annex commands directly, it is useful to appreciate that DataLad is
build atop of very flexible and powerful tools.  Knowing basics of git and
git-annex in addition to DataLad helps to not only make better use of
DataLad but also to enable more advanced and more efficient data management
scenarios. DataLad makes use of lower-level configuration and data structures
as much as possible. Consequently, it is possible to manipulate DataLad
datasets with low-level tools if needed. Moreover, DataLad datasets are
compatible with tools and services designed to work with plain Git repositories,
such as the popular GitHub_ service.

.. _github: https://github.com

To better illustrate the different scopes, the following table provides an
overview of the features that are contributed by each software technology
layer.

================================================   =============  ===============   ==============
Feature                                             Git            Git-annex         DataLad
================================================   =============  ===============   ==============
Version control (text, code)                       |tup|          |tup| can mix     |tup| can mix
Version control (binary data)                      (not advised)  |tup|             |tup|
Auto-crawling available resources                                 |tup| RSS feeds   |tup| flexible
Unified dataset handling                                                            |tup|
- recursive operation on datasets                                                   |tup|
- seamless operation across datasets boundaries                                     |tup|
- meta-data support                                               |tup| per-file    |tup|
- meta-data aggregation                                                             |tup| flexible
Unified authentication interface                                                    |tup|
================================================   =============  ===============   ==============

.. |tup| unicode:: U+2713 .. check mark
   :trim:
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_cmdline:

**********************
Command line reference
**********************

Main command
============

.. toctree::
   :maxdepth: 1

   datalad: Main command entrypoint <generated/man/datalad>

Core commands
=============

A minimal set of commands that cover essential functionality. Core commands
receive special scrutiny with regard API composition and (breaking) changes.

Local operation
---------------

.. toctree::
   :maxdepth: 1

   datalad create: Create a new dataset <generated/man/datalad-create>
   datalad save: Save the state of a dataset <generated/man/datalad-save>
   datalad run: Run a shell command and record its impact on a dataset <generated/man/datalad-run>
   datalad status: Report on the state of dataset content <generated/man/datalad-status>
   datalad diff: Report differences between two states of a dataset <generated/man/datalad-diff>

Distributed operation
---------------------

.. toctree::
   :maxdepth: 1

   datalad clone: Obtain a dataset (sibling) from another location <generated/man/datalad-clone>
   datalad push: Push updates/data to a dataset sibling <generated/man/datalad-push>


Extended set of functionality
=============================

Dataset operations
------------------

.. toctree::
   :maxdepth: 1

   datalad add-readme: Add information on DataLad dataset to a README <generated/man/datalad-add-readme>
   datalad addurls: Update dataset content from a list of URLs <generated/man/datalad-addurls>
   datalad copy-file: Copy file identity and availability from one dataset to another <generated/man/datalad-copy-file>
   datalad drop: Drop datasets or dataset components <generated/man/datalad-drop>
   datalad get: Obtain any dataset content <generated/man/datalad-get>
   datalad install: Install a dataset from a (remote) source <generated/man/datalad-install>
   datalad no-annex: Configure a dataset to never put file content into an annex <generated/man/datalad-no-annex>
   datalad remove: Unlink components from a dataset <generated/man/datalad-remove>
   datalad subdatasets: Query and manipulate subdataset records of a dataset <generated/man/datalad-subdatasets>
   datalad unlock: Make dataset file content editable <generated/man/datalad-unlock>


Dataset siblings and 3rd-party platform support
-----------------------------------------------

.. toctree::
   :maxdepth: 1

   datalad siblings: Query and manipulate sibling configuration of a dataset <generated/man/datalad-siblings>
   datalad create-sibling: Create a sibling on an SSH-accessible machine <generated/man/datalad-create-sibling>
   datalad create-sibling-github: Create a sibling on GitHub <generated/man/datalad-create-sibling-github>
   datalad create-sibling-gitlab: Create a sibling on GitLab <generated/man/datalad-create-sibling-gitlab>
   datalad create-sibling-gogs: Create a sibling on GOGS <generated/man/datalad-create-sibling-gogs>
   datalad create-sibling-gitea: Create a sibling on Gitea <generated/man/datalad-create-sibling-gitea>
   datalad create-sibling-gin: Create a sibling on GIN (with content hosting) <generated/man/datalad-create-sibling-gin>
   datalad create-sibling-ria: Create a sibling in a RIA store <generated/man/datalad-create-sibling-ria>
   datalad export-archive: Export dataset content as a TAR/ZIP archive <generated/man/datalad-export-archive>
   datalad export-archive-ora: Export a local dataset annex for the ORA remote <generated/man/datalad-export-archive-ora>
   datalad export-to-figshare: Export dataset content as a ZIP archive to figshare <generated/man/datalad-export-to-figshare>
   datalad update: Obtain and incorporate updates from dataset siblings <generated/man/datalad-update>


Reproducible execution
----------------------

Extending the functionality of the core ``run`` command.

.. toctree::
   :maxdepth: 1

   datalad rerun: Re-execute previous datalad-run commands <generated/man/datalad-rerun>
   datalad run-procedure: Run prepared procedures (DataLad scripts) on a dataset <generated/man/datalad-run-procedure>


Metadata handling
-----------------

.. toctree::
   :maxdepth: 1

   datalad search: Query metadata of a dataset <generated/man/datalad-search>
   datalad metadata: Report known metadata on particular datasets or files <generated/man/datalad-metadata>
   datalad aggregate-metadata: Assemble metadata from datasets for later query <generated/man/datalad-aggregate-metadata>
   datalad extract-metadata: Run metadata extractor on a dataset or file <generated/man/datalad-extract-metadata>


Helpers and support utilities
-----------------------------

.. toctree::
   :maxdepth: 1

   datalad add-archive-content: Extract and add the content of an archive to a dataset <generated/man/datalad-add-archive-content>
   datalad clean: Remove temporary left-overs of DataLad operations <generated/man/datalad-clean>
   datalad check-dates: Scan a dataset for dates and timestamps <generated/man/datalad-check-dates>
   datalad create-test-dataset: Test helper <generated/man/datalad-create-test-dataset>
   datalad download-url: Download helper with support for DataLad's credential system <generated/man/datalad-download-url>
   datalad foreach-dataset: Run a command or Python code on the dataset and/or each of its sub-datasets <generated/man/datalad-foreach-dataset>
   datalad sshrun: Remote command execution using DataLad's connection management <generated/man/datalad-sshrun>
   datalad shell-completion: Helper to support command completion <generated/man/datalad-shell-completion>
   datalad test: Frontend for running DataLad's internal test battery <generated/man/datalad-test>
   datalad wtf: Report on a DataLad installation and its configuration <generated/man/datalad-wtf>


Deprecated commands
-------------------

.. toctree::
   :maxdepth: 1

   datalad uninstall: Drop subdatasets <generated/man/datalad-uninstall>
.. _configuration:

Configuration
*************

DataLad uses the same configuration mechanism and syntax as Git itself.
Consequently, datalad can be configured using the :command:`git config`
command. Both a *global* user configuration (typically at
:file:`~/.gitconfig`), and a *local* repository-specific configuration
(:file:`.git/config`) are inspected.

In addition, datalad supports a persistent dataset-specific configuration.
This configuration is stored at :file:`.datalad/config` in any dataset.  As it
is part of a dataset, settings stored there will also be in effect for any
consumer of such a dataset. Both *global* and *local* settings on a particular
machine always override configuration shipped with a dataset.

All datalad-specific configuration variables are prefixed with ``datalad.``.

It is possible to override or amend the configuration using environment
variables. Any variable with a name that starts with ``DATALAD_`` will
be available as the corresponding ``datalad.`` configuration variable,
replacing any ``__`` (two underscores) with a hyphen, then any ``_``
(single underscore) with a dot, and finally converting all letters to
lower case. Values from environment variables take precedence over
configuration file settings.

In addition, the ``DATALAD_CONFIG_OVERRIDES_JSON`` environment variable can
be set to a JSON record with configuration values.  This is
particularly useful for options that aren't accessible through the
naming scheme described above (e.g., an option name that includes an
underscore).

The following sections provide a (non-exhaustive) list of settings honored
by datalad. They are categorized according to the scope they are typically
associated with.


Global user configuration
=========================

.. include:: generated/cfginfo/global.rst.in

Local repository configuration
==============================

.. include:: generated/cfginfo/local.rst.in

Sticky dataset configuration
=============================

.. include:: generated/cfginfo/dataset.rst.in

Miscellaneous configuration
===========================

.. include:: generated/cfginfo/misc.rst.in
.. _chap_metadata:

Metadata
********

Overview
========

DataLad has built-in, modular, and extensible support for metadata in various
formats. Metadata is extracted from a dataset and its content by one or more
extractors that have to be enabled in a dataset's configuration. Extractors
yield metadata in a JSON-LD_-like structure that can be arbitrarily complex and
deeply nested. Metadata from each extractor is kept unmodified, unmangled, and
separate from metadata of other extractors. This design enables tailored
applications using particular metadata that can use Datalad as a
content-agnostic aggregation and transport layer without being limited or
impacted by other metadata sources and schemas.

Extracted metadata is stored in a dataset in (compressed) files using a JSON
stream format, separately for metadata describing a dataset as a whole, and
metadata describing individual files in a dataset. This limits the amount of
metadata that has to be obtained and processed for applications that do not
require all available metadata.

DataLad provides a content-agnostic metadata aggregation mechanism that
stores metadata of sub-datasets (with arbitrary nesting levels) in a
superdataset, where it can then be queried without having the subdatasets
locally present.

Lastly, DataLad comes with a `search` command that enable metadata queries
via a flexible query language. However, alternative applications for metadata
queries (e.g. graph-based queries) can be built on DataLad, by requesting
a complete or partial dump of aggregated metadata available in a dataset.

.. _JSON-LD: http://json-ld.org/
.. _linked data: https://en.wikipedia.org/wiki/Linked_data


Supported metadata sources
==========================

This following sections provide an overview of included metadata extractors for
particular types of data structures and file formats. Note that :ref:`DataLad
extension packages <chap_customization>`, such as the `neuroimaging extension
<https://github.com/datalad/datalad-neuroimaging>`_, can provide additional
extractors for particular domains and formats.

Only :ref:`annex <metadata-annex>` and :ref:`datalad_core <metadata-datalad_core>`
extractors are enabled by default.  Any additional metadata extractor should be
enabled by setting the :term:`datalad.metadata.nativetype` :ref:`configuration <configuration>` variable
via the ``git config`` command or by editing ``.datalad/config`` directly.
For example, ``git config -f .datalad/config --add datalad.metadata.nativetype audio``
would add :ref:`audio <metadata-audio>` metadata extractor to the list.


.. _metadata-annex:

Annex metadata (``annex``)
--------------------------

Content tracked by git-annex can have associated
`metadata records <http://git-annex.branchable.com/metadata/>`_.
From DataLad's perspective, git-annex metadata is just another source of
metadata that can be extracted and aggregated.

You can use the `git-annex metadata`_ command to assign git-annex
metadata.  And, if you have a table or records that contain data
sources and metadata, you can use :ref:`datalad addurls <man_datalad-addurls>`
to quickly populate a dataset with files and associated
git-annex metadata. (`///labs/openneurolab/metasearch
<https://datasets.datalad.org/?dir=/labs/openneurolab/metasearch>`_ is
an example of such a dataset.)


Pros of git-annex level metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- Many git-annex commands, such as `git-annex get`_ and `git-annex copy`_, can
  use metadata to decide which files (keys) to operate on, making it possible to
  automate file (re)distribution based on their metadata annotation
- Assigned metadata is available for use by git-annex right away without
  requiring any additional "aggregation" step
- `git-annex view`_ can be used to quickly generate completely new layouts
  of the repository solely based on the metadata fields associated with the files

.. _git-annex get: https://git-annex.branchable.com/git-annex-get/
.. _git-annex copy: https://git-annex.branchable.com/git-annex-copy/
.. _git-annex metadata: https://git-annex.branchable.com/git-annex-metadata/
.. _git-annex view: https://git-annex.branchable.com/git-annex-view/


Cons of git-annex level metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- Metadata fields are actually stored per git-annex key rather than per file.
  If multiple files contain the same content, metadata will be shared among them.
- Files whose content is tracked directly by git cannot have git-annex metadata assigned.
- No per repository/directory metadata, and no mechanism to use/aggregate
  metadata from sub-datasets
- Field names cannot contain some symbols, such as ':'
- Metadata is stored within the `git-annex` branch, so it is distributed
  across all clones of the dataset, making it hard to scale for large metadata
  sizes or to work with sensitive metadata (not intended to be redistributed)
- It is a generic storage with no prescribed vocabularly,
  making it very flexible but also requiring consistency and
  harmonization to make the stored metadata useful for search


Example uses of git-annex metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Annotating files for different purposes
#######################################

FreeSurfer project `uses <https://surfer.nmr.mgh.harvard.edu/fswiki/DevelopersGuide_git#GettheDataFiles>`_
`git-annex` for managing their source code+data base within a single
git/git-annex repository. Files necessary for different scenarios (deployment,
testing) are annotated and can be fetched selectively for the scenario at hand.

Automating "non-distribution" of sensitive files
################################################

In the `ReproIn <http://reproin.repronim.org>`_ framework for automated
conversion of BIDS dataset and in some manually prepared datasets
(such as
`///labs/gobbini/famface/data <https://datasets.datalad.org/?dir=/labs/gobbini/famface/data>`_
and
`///labs/haxby/raiders <https://datasets.datalad.org/?dir=/labs/haxby/raiders>`_),
we annotated materials that must not be publicly shared with a git-annex
metadata field ``distribution-restrictions``.  We used the following of values to
describe why any particular file (content) should not be redistributed:

- **sensitive** - files which potentially contain participant sensitive
  information, such as non-defaced anatomicals
- **proprietary** - files which contain proprietary data, which we have no
  permissions to share (e.g., movie video files)

Having annotated files this way, we could instruct git-annex
to publish all but those restricted files to our
server: ``git annex wanted datalad-public "not metadata=distribution-restrictions=*"``.

.. warning::
  The above setup depends on ``git annex copy --auto`` deciding to *not*
  copy the content.  To avoid inadvertently publishing sensitive data,
  make sure that public targets ("datalad-public" in the example
  above) do not want the content for another reason, in particular due
  to ``numcopies`` or required content configuration.  If ``numcopies``
  is set to a value greater than 1 (the default) and the requested
  number of copies cannot be verified, ``git annex copy --auto`` will
  transfer the data regardless of the preferred content expression set
  by the ``git annex wanted`` call above.


Flexible directory layout
#########################

If you are maintaining a collection of music files or PDFs for the lab, you
may want to display the files in an alternative or filtered hierarchy.
`git-annex view`_ could be of help. Example:

.. code-block:: sh

  datalad install ///labs/openneurolab/metasearch
  cd metasearch
  git annex view sex=* handedness=ambidextrous

would give you two directories (Male, Female) with only the files belonging to
ambidextrous subjects.


.. _metadata-audio:

Various audio file formats (``audio``)
--------------------------------------

This extractor uses the `mutagen <https://github.com/quodlibet/mutagen>`_
package to extract essential metadata from a range of audio file formats.  For
the most common metadata properties a constrained vocabulary, based on the
`Music Ontology <http://purl.org/ontology/mo/>`_ is employed.

datacite.org compliant datasets (``datacite``)
----------------------------------------------

This extractor can handle dataset-level metadata following the `datacite.org
<https://www.datacite.org>`_ specification. No constrained vocabulary is
identified at the moment.

.. _metadata-datalad_core:

Datalad's internal metadata storage (``datalad_core``)
------------------------------------------------------

This extractor can express Datalad's internal metadata representation, such
as the relationship of a super- and a subdataset. It uses DataLad's own
constrained vocabulary.

RFC822-compliant metadata (``datalad_rfc822``)
----------------------------------------------

This is a custom metadata format, inspired by the standard used for Debian
software packages that is particularly suited for manual entry. This format is
a good choice when metadata describing a dataset as a whole cannot be obtained
from some other structured format. The syntax is :rfc:`822`-compliant. In other
words: this is a text-based format that uses the syntax of email headers.
Metadata must be placed in ``DATASETROOT/.datalad/meta.rfc822`` for this format.

.. _RFC822: https://tools.ietf.org/html/rfc822

Here is an example:

.. code-block:: none

  Name: myamazingdataset
  Version: 1.0.0-rc3
  Description: Basic summary
   A text with arbitrary length and content that can span multiple
   .
   paragraphs (this is a new one)
  License: CC0-1.0
   The person who associated a work with this deed has dedicated the work to the
   public domain by waiving all of his or her rights to the work worldwide under
   copyright law, including all related and neighboring rights, to the extent
   allowed by law.
   .
   You can copy, modify, distribute and perform the work, even for commercial
   purposes, all without asking permission.
  Homepage: http://example.com
  Funding: Grandma's and Grandpa's support
  Issue-Tracker: https://github.com/datalad/datalad/issues
  Cite-As: Mike Author (2016). We made it. The breakthrough journal of unlikely
    events. 1, 23-453.
  DOI: 10.0000/nothere.48421

The following fields are supported:

``Audience``:
  A description of the target audience of the dataset.
``Author``:
  A comma-delimited list of authors of the dataset, preferably in the format.
  ``Firstname Lastname <Email Adress>``
``Cite-as``:
  Instructions on how to cite the dataset, or a structured citation.
``Description``:
  Description of the dataset as a whole. The first line should represent a
  compact short description with no more than 6-8 words.
``DOI``:
  A `digital object identifier <https://en.wikipedia.org/wiki/Digital_object_identifier>`_
  for the dataset.
``Funding``:
  Information on potential funding for the creation of the dataset and/or its
  content. This field can also be used to acknowledge non-monetary support.
``Homepage``:
  A URL to a project website for the dataset.
``Issue-tracker``:
  A URL to an issue tracker where known problems are documented and/or new
  reports can be submitted.
``License``:
  A description of the license or terms of use for the dataset. The first
  lines should be the SPDX License Identifier from the `SPDX License List <https://spdx.org/licenses/>`_
  (e.g. "CC0-1.0" or "PPDL-1.0"). More complex licensing situation can be expressed
  using
  `SPDX License Expressions <https://spdx.github.io/spdx-spec/appendix-IV-SPDX-license-expressions/>`_.
  Full license texts or term descriptions can be included.
``Maintainer``:
  Can be used in addition and analog to ``Author``, when authors (creators of
  the data) need to be distinguished from maintainers of the dataset.
``Name``:
  A short name for the dataset. It may be beneficial to avoid special
  characters, umlauts, spaces, etc. to enable widespread use of this name
  for URL, catalog keys, etc. in unmodified form.
``Version``:
  A version for the dataset. This should be in a format that is alphanumerically
  sortable and lead to a "greater" version for an update of a dataset.

Metadata keys used by this extractor are defined in DataLad's own constrained
vocabulary.

Friction-less data packages (``frictionless_datapackage``)
----------------------------------------------------------

DataLad has basic support for extraction of essential dataset-level metadata
from `friction-less data packages
<http://specs.frictionlessdata.io/data-packages>`_ (``datapackage.json``).
file. Metadata keys are constrained to DataLad's own vocabulary.

Exchangeable Image File Format (``exif``)
-----------------------------------------

The extractor yields EXIF metadata from any compatible file. It uses
the W3C EXIF vocabulary (http://www.w3.org/2003/12/exif/ns/).

Various image/photo formats (``image``)
---------------------------------------

Standard image metadata is extractor using the `Pillow package
<https://github.com/python-pillow/Pillow>`_. Core metadata is available using
an adhoc vocabulary defined by the extractor.

Extensible Metadata Platform (``xmp``)
--------------------------------------

This extractor yields any XMP-compliant metadata from any supported file (e.g.
PDFs, photos). XMP metadata uses fully qualified terms from standard
vocabularies that are simply passed through by the extractor. At the moment
metadata extraction from side-car files is not supported, but would be easy to
add.

Metadata aggregation and query
==============================

Metadata aggregation can be performed with the :ref:`aggregate-metadata
<man_datalad-aggregate-metadata>` command. Aggregation is done for two
interrelated but distinct reasons:

- Fast uniform metadata access, independent of local data availability
- Comprehensive data discovery without access to or knowledge of individual
  datasets

In an individual dataset, metadata aggregation engages any number of enabled
metadata extractors to build a JSON-LD based metadata representation that is
separate from the original data files. These metadata objects are added to the
dataset and are tracked with the same mechanisms that are used for any other
dataset content. Based on this metadata, DataLad can provide fast and uniform
access to metadata for any dataset component (individual files, subdatasets,
the whole dataset itself), based on the relative path of a component within a
dataset (available via the :ref:`metadata <man_datalad-metadata>` command).
This extracted metadata can be kept or made available locally for any such
query, even when it is impossible or undesirable to keep the associated data
files around (e.g. due to size constraints).

For any superdataset (a dataset that contains subdatasets as components),
aggregation can go one step further. In this case, aggregation imports
extracted metadata from subdatasets into the superdataset to offer the just
described query feature for any aggregated subdataset too. This works across
any number of levels of nesting. For example, a subdataset that contains the
aggregated metadata for eight other datasets (that might have never been
available locally) can be aggregated into a local superdataset with all its
metadata. In that superdataset, a DataLad user is then able to query
information on any content of any subdataset, regardless of their actual
availability. This principle also allows any user to install the superdataset
from https://datasets.datalad.org and perform *local and offline* queries about
any dataset available online from this server.

Besides full access to all aggregated metadata by path (via the :ref:`metadata
<man_datalad-metadata>` command), DataLad also comes with a :ref:`search
<man_datalad-search>` command that provides different search modes to query the
entirety of the locally available metadata. Its capabilities include simple
keyword searches as well as more complex queries using date ranges or logical
conjunctions.

Internal metadata representation
================================

.. warning::
  The information in this section is meant to provide insight into how
  DataLad structures extracted and aggregated metadata. However, this
  representation is not considered stable or part of the public API,
  hence these data should not be accessed directly. Instead, all
  metadata access should happen via the :command:`metadata` API command.

A dataset's metadata is stored in the `.datalad/metadata` directory. This
directory contains two main elements:

- a metadata inventory or catalog
- a store for metadata "objects"

The metadata inventory
----------------------

The inventory is kept in a JSON file, presently named ``aggregate_v1.json``.
It contains a single top-level dictionary/object. Each element in this
dictionary represents one subdataset from which metadata has been extracted
and aggregated into the dataset at hand. Keys in this dictionary are
paths to the respective (sub)datasets (relative to the root of the dataset).
If a dataset has no subdataset and metadata extraction was performed, the
dictionary will only have a single element under the key ``"."``.

Here is an excerpt of an inventory dictionary showing the record of the
root dataset itself.

.. code-block:: json

   {

      ".": {
         "content_info":
            "objects/0c/cn-b046b2c3a5e2b9c5599c980c7b5fab.xz",
         "datalad_version":
            "0.10.0.rc4.dev191",
         "dataset_info":
            "objects/0c/ds-b046b2c3a5e2b9c5599c980c7b5fab",
         "extractors": [
            "datalad_core",
            "annex",
            "bids",
            "nifti1"
         ],
         "id":
            "00ce405e-6589-11e8-b749-a0369fb55db0",
         "refcommit":
            "d170979ef33a82c67e6fefe3084b9fe7391b422b"
      },

   }

The record of each dataset contains the following elements:

``id``
  The DataLad dataset UUID of the dataset metadata was extracted and
  aggregated from.
``refcommit``
  The SHA sum of the last metadata-relevant commit in the history of
  the dataset metadata was extracted from. Metadata-relevant commits
  are any commits that modify dataset content that is not exclusively
  concerning DataLad's own internal status and configuration.
``datalad_version``
  The version string of the DataLad version that was used to perform
  the metadata extraction (not necessarily the metadata aggregation,
  as pre-extracted metadata can be aggregated from other superdatasets
  for a dataset that is itself not available locally).
``extractors``
  A list with the names of all enabled metadata extractors for this
  dataset. This list may include names for extractors that are provided
  by extensions, and may not be available for any given DataLad
  installation.
``content_info``, ``dataset_info``
  Path to the object files containing the actual metadata on the dataset
  as a whole, and on individual files in a dataset (content). Paths
  are to be interpreted relative to the inventory file, and point to
  the metadata object store.

Read-access to the metadata inventory is available via the ``metadata``
command and its ``--get-aggregates`` option.

The metadata object store
-------------------------

The object store holds the files containing dataset and content metadata for
each aggregated dataset. The object store is located in
`.datalad/metadata/objects`. However, this directory itself and the
subdirectory structure within it have no significance, they are completely
defined and exclusively discoverable via the ``content_info`` and
``dataset_info`` values in the metadata inventory records.

Metadata objects for datasets and content use a slightly different internal
format. Both files could be either compressed (XZ) or uncompressed. Current
practice uses compression for content metadata, but not for dataset metadata.
Any metadata object file could be directly committed to Git, or it could be
tracked via Git-annex. Reasons to choose one over the other could be file size,
or privacy concerns.

Read-access to the metadata objects of dataset and individual files is
available via the ``metadata`` command. Importantly, metadata can be requested


Metadata objects for datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These files have a single top-level JSON object/dictionary as content. A
JSON-LD ``@content`` field is used to assign a semantic markup to allow for
programmatic interpretation of metadata as linked data. Any other top-level key
identifies the name of a metadata extractor, and the value stored under this
key represents the output of the corresponding extractor.

Structure and content of an extractor's output are unconstrained and completely
up to the implementation of that particular extractor. Extractor can report
additional JSON-LD context information (but there is no requirement).

The output of one extractor does not interfere or collide with the output
of any other extractor.

Metadata objects for content/file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In contrast to metadata objects for entire datasets, these files use a JSON
stream format, i.e. one JSON object/dictionary per line (no surrounding list).
This makes it possible to process the content line-by-line instead of having
to load an entire files (with potentially millions of records).

The only other difference to dataset metadata objects is an additional top-level
key ``path`` that identifies the relative path (relative to the root of its parent
dataset) of the file the metadata record is associated with.

Otherwise, the extractor-specific metadata structure and content is unconstrained.

Content metadata objects tend to contain massively redundant information (e.g.
a dataset with a thousand 12 megapixel images will report the identical resolution
information a thousand times). Therefore, content metadata objects are by default
XZ compressed -- as this compressor is particularly capable discovering such
redundancy and yield a very compact file size.

The reason for gathering all metadata into a single file across all content files and
metadata extractors is to limit the impact on the performance of the underlying
Git repository. Large superdataset could otherwise quickly grow into dimensions
where tens of thousands of files would be required just to manage the metadata.
Such a configuration would also limit the compatibility of DataLad datasets with
constrained storage environments (think e.g. inode limits on super computers),
as these files are tracked in Git and would therefore be present in any copy,
regardless of whether metadata access is desired or not.


Vocabulary
==========

The following sections describe details and changes in the metadata
specifications implemented in datalad.

.. _2.0:

`v2.0 <http://docs.datalad.org/schema_v2.0.json>`_
--------------------------------------------------

* Current development version that will be released together with
  DataLad v0.10.

.. _1.0:

`v1.0 <http://docs.datalad.org/schema_v1.0.json>`_
--------------------------------------------------

* Original implementation that did not really see the light of the day.
.. -*- mode: rst; fill-column: 78; indent-tabs-mode: nil -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  #
  #   See COPYING file distributed along with the datalad package for the
  #   copyright and license terms.
  #
  ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

.. _chap_glossary:

********
Glossary
********

DataLad purposefully uses a terminology that is different from the one used by
its technological foundations Git_ and git-annex_. This glossary provides
definitions for terms used in the datalad documentation and API, and relates
them to the corresponding Git_/git-annex_ concepts.

.. glossary::
  :sorted:

  dataset
    A regular Git_ repository with an (optional) :term:`annex`.

  subdataset
    A :term:`dataset` that is part of another dataset, by means of being
    tracked as a Git_ submodule. As such, a subdataset is also a complete
    dataset and not different from a standalone dataset.

  superdataset
    A :term:`dataset` that contains at least one :term:`subdataset`.

  sibling
    A :term:`dataset` (location) that is related to a particular dataset,
    by sharing content and history. In Git_ terminology, this is a *clone*
    of a dataset that is configured as a *remote*.

  annex
    Extension to a Git_ repository, provided and managed by git-annex_ as
    means to track and distribute large (and small) files without having to
    inject them directly into a Git_ repository (which would slow Git
    operations significantly and impair handling of such repositories in
    general).

  CLI
    A `Command Line Interface`_. Could be used interactively by executing
    commands in a `shell`_, or as a programmable API for shell scripts.

  DataLad extension
    A Python package, developed outside of the core DataLad codebase, which
    (when installed) typically either provides additional top level `datalad`
    commands and/or additional metadata extractors.  Visit
    `Handbook, Ch.2. DataLad’s extensions <http://handbook.datalad.org/en/latest/basics/101-144-intro_extensions.html>`_
    for a representative list of extensions and instructions on how to install
    them.

.. _Git: https://git-scm.com
.. _Git-annex: http://git-annex.branchable.com
.. _`Command Line Interface`: https://en.wikipedia.org/wiki/Command-line_interface
.. _shell: https://en.wikipedia.org/wiki/Shell_(computing).. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_basic_principles:

****************
Basic principles
****************

DataLad is designed to be used both as a command-line tool, and as a Python
module. The sections :ref:`chap_cmdline` and :ref:`chap_modref` provide
detailed description of the commands and functions of the two interfaces.  This
section presents common concepts.  Although examples will frequently be
presented using command line interface commands, all functionality with
identically named functions and options are available through Python API as
well.

Datasets
========

A DataLad :term:`dataset` is a Git repository that may or may not have a data
:term:`annex` that is used to manage data referenced in a dataset. In practice,
most DataLad datasets will come with an annex.

Types of IDs used in datasets
-----------------------------

Four types of unique identifiers are used by DataLad to enable identification
of different aspects of datasets and their components.

Dataset ID
  A UUID that identifies a dataset as a whole across its entire history and
  flavors. This ID is stored in a dataset's own configuration file
  (``<dataset root>/.datalad/config``) under the configuration key
  ``datalad.dataset.id``.
  As this configuration is stored in a file that is part of the Git history of
  a dataset, this ID is identical for all "clones" of a dataset and across all
  its versions. If the purpose or scope of a dataset changes enough to warrant
  a new dataset ID, it can be changed by altering the dataset configuration
  setting.
Annex ID
  A UUID assigned to an annex of each individual clone of a dataset repository.
  Git-annex uses this UUID to track file content availability information. The
  UUID is available under the configuration key ``annex.uuid`` and is stored
  in the configuration file of a local clone (``<dataset root>/.git/config``).
  A single dataset instance (i.e. clone) can only have a single annex UUID,
  but a dataset with multiple clones will have multiple annex UUIDs.
Commit ID
  A Git hexsha or tag that identifies a version of a dataset. This ID uniquely
  identifies the content and history of a dataset up to its present state. As
  the dataset history also includes the dataset ID, a commit ID of a DataLad
  dataset is unique to a particular dataset.
Content ID
  Git-annex key (typically a checksum) assigned to the content of a file in
  a dataset's annex. The checksum reflects the content of a file, not its name.
  Hence the content of multiple identical files in a single (or across)
  dataset(s) will have the same checksum. Content IDs are managed by Git-annex
  in a dedicated ``annex`` branch of the dataset's Git repository.


Dataset nesting
---------------

Datasets can contain other datasets (:term:`subdataset`\s), which can in turn
contain subdatasets, and so on. There is no limit to the depth of nesting
datasets. Each dataset in such a hierarchy has its own annex and its own
history. The parent or :term:`superdataset` only tracks the specific state of a
subdataset, and information on where it can be obtained. This is a powerful yet
lightweight mechanism for combining multiple individual datasets for a specific
purpose, such as the combination of source code repositories with other
resources for a tailored application. In many cases DataLad can work with a
hierarchy of datasets just as if it were a single dataset. Here is a demo:

.. include:: basics_nesteddatasets.rst.in
   :start-after: Let's create a dataset
   :end-before:  ___________________________


Dataset collections
-------------------

A superdataset can also be seen as a curated collection of datasets, for example,
for a certain data modality, a field of science, a certain author, or from
one project (maybe the resource for a movie production). This lightweight
coupling between super and subdatasets enables scenarios where individual datasets
are maintained by a disjoint set of people, and the dataset collection itself can
be curated by a completely independent entity. Any individual dataset can be
part of any number of such collections.

Benefiting from Git's support for workflows based on decentralized "clones" of
a repository, DataLad's datasets can be (re-)published to a new location
without loosing the connection between the "original" and the new "copy". This
is extremely useful for collaborative work, but also in more mundane scenarios
such as data backup, or temporary deployment of a dataset on a compute cluster,
or in the cloud.  Using git-annex, data can also get synchronized across
different locations of a dataset (:term:`sibling`\s in DataLad terminology).
Using metadata tags, it is even possible to configure different levels of
desired data redundancy across the network of dataset, or to prevent
publication of sensitive data to publicly accessible repositories. Individual
datasets in a hierarchy of (sub)datasets need not be stored at the same location.
Continuing with an earlier example, it is possible to post a curated
collection of datasets, as a superdataset, on Github, while the actual datasets
live on different servers all around the world.

Basic command line usage
========================

.. include:: basics_cmdline.rst.in
   :end-before:  ___________________________


API principles
==============

You can use DataLad's ``install`` command to download datasets. The command accepts
URLs of different protocols (``http``, ``ssh``) as an argument. Nevertheless, the easiest way
to obtain a first dataset is downloading the default :term:`superdataset` from
https://datasets.datalad.org/ using a shortcut.

Downloading DataLad's default superdataset
--------------------------------------------

https://datasets.datalad.org provides a super-dataset consisting of datasets
from various portals and sites.  Many of them were crawled, and periodically
updated, using `datalad-crawler <https://github.com/datalad/datalad-crawler>`__
extension.  The argument ``///`` can be used
as a shortcut that points to the superdataset located at https://datasets.datalad.org/. 
Here are three common examples in command line notation:

``datalad install ///``
    installs this superdataset (metadata without subdatasets) in a
    `datasets.datalad.org/` subdirectory under the current directory
``datalad install -r ///openfmri``
    installs the openfmri superdataset into an `openfmri/` subdirectory.
    Additionally, the ``-r`` flag recursively downloads all metadata of datasets 
    available from http://openfmri.org as subdatasets into the `openfmri/` subdirectory
``datalad install -g -J3 -r ///labs/haxby``
    installs the superdataset of datasets released by the lab of Dr. James V. Haxby
    and all subdatasets' metadata. The ``-g`` flag indicates getting the actual data, too.
    It does so by using 3 parallel download processes (``-J3`` flag).

:ref:`datalad search <man_datalad-search>` command, if ran outside of any dataset,
will install this default superdataset under a path specified in
``datalad.locations.default-dataset`` :ref:`configuration <configuration>`
variable (by default ``$HOME/datalad``).

Downloading datasets via http
-----------------------------

In most places where DataLad accepts URLs as arguments these URLs can be
regular ``http`` or ``https`` protocol URLs. For example:

``datalad install https://github.com/psychoinformatics-de/studyforrest-data-phase2.git``

Downloading datasets via ssh
----------------------------
DataLad also supports SSH URLs, such as ``ssh://me@localhost/path``.

``datalad install ssh://me@localhost/path``

Finally, DataLad supports SSH login style resource identifiers, such as ``me@localhost:/path``.

``datalad install me@localhost:/path``


Commands `install` vs `get`
---------------------------

The ``install`` and ``get`` commands might seem confusingly similar at first.
Both of them could be used to install any number of subdatasets, and fetch
content of the data files.  Differences lie primarily in their default
behaviour and outputs, and thus intended use.  Both ``install`` and ``get``
take local paths as their arguments, but their default behavior and output
might differ;

- **install** primarily operates and reports at the level of **datasets**, and
  returns as a result dataset(s)
  which either were just installed, or were installed previously already under
  specified locations.   So result should be the same if the same ``install``
  command ran twice on the same datasets.  It **does not fetch** data files by
  default

- **get** primarily operates at the level of **paths** (datasets, directories, and/or
  files). As a result it returns only what was installed (datasets) or fetched
  (files).  So result of rerunning the same ``get`` command should report that
  nothing new was installed or fetched.  It **fetches** data files by default.

In how both commands operate on provided paths, it could be said that ``install
== get -n``, and ``install -g == get``.  But ``install`` also has ability to
install new datasets from remote locations given their URLs (e.g.,
``https://datasets.datalad.org/`` for our super-dataset) and SSH targets (e.g.,
``[login@]host:path``) if they are provided as the argument to its call or
explicitly as ``--source`` option.  If ``datalad install --source URL
DESTINATION`` (command line example) is used, then dataset from URL gets
installed under PATH. In case of ``datalad install URL`` invocation, PATH is
taken from the last name within URL similar to how ``git clone`` does it.  If
former specification allows to specify only a single URL and a PATH at a time,
later one can take multiple remote locations from which datasets could be
installed.

So, as a rule of thumb -- if you want to install from external URL or fetch a
sub-dataset without downloading data files stored under annex -- use ``install``.
In Python API ``install`` is also to be used when you want to receive in output the
corresponding Dataset object to operate on, and be able to use it even if you
rerun the script. In all other cases, use ``get``.
Publications
************

Further conceptual and technical information on DataLad, and applications built on DataLad,
are available from the publications listed below.

The best of both worlds: Using semantic web with JSOB-LD. An example with NIDM Results & DataLad [poster]
   - Camille Maumet, Satrajit Ghosh, Yaroslav O. Halchenko, Dorota Jarecka, Nolan Nichols, Jean-Baptist POline, Michael Hanke

One thing to bind them all: A complete raw data structure for auto-generation of BIDS datasets [poster]
   - Benjamin Poldrack, Kyle Meyer, Yaroslav O. Halchenko, Michael Hanke

Fantastic containers and how to tame them [poster]
   - Yaroslav O. Halchenko, Kyle Meyer, Matt Travers, Dorota Jarecka, Satrajit Ghosh, Jakub Kaczmarzyk, Michael Hanke

YODA: YODA's Organigram on Data Analysis [poster]
   - An outline of a simple approach to structuring and conducting data analyses that aims to
     tightly connect all their essential ingredients: data, code, and computational environments
     in a transparent, modular, accountable, and practical way.
   - Michael Hanke, Kyle A. Meyer, Matteo Visconti di Oleggio Castello, Benjamin Poldrack, Yaroslav O. Halchenko
   - F1000Research 2018, 7:1965 (https://doi.org/10.7490/f1000research.1116363.1)

Go FAIR with DataLad [talk]
   - On DataLad's capabilities to create and maintain Findable, Accessible, Interoperable, and Re-Usable (FAIR)
     resources.
   - Michael Hanke, Yaroslav O. Halchenko
   - Bernstein Conference 2018 workshop: Practical approaches to research data management and reproducibility
     (`slides <https://rawgit.com/psychoinformatics-de/talk-datalad-gofair/master/index.html>`__)
   - OpenNeuro kick-off meeting, 2018, Stanford (`slide sources <https://github.com/datalad/talk-openneuro-2018>`__)
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_customization:

********************************************
Customization and extension of functionality
********************************************

DataLad provides numerous commands that cover many use cases. However, there
will always be a demand for further customization or extensions of built-in
functionality at a particular site, or for an individual user. DataLad
addresses this need with a mechanism for extending particular Datalad
functionality, such as metadata extractor, or providing entire command suites
for a specialized purpose.

As the name suggests, a :term:`DataLad extension` package is a proper Python package.
Consequently, there is a significant amount of boilerplate code involved in the
creation of a new Datalad extension. However, this overhead enables a number of
useful features for extension developers:

- extensions can provide any number of additional commands that can be grouped into
  labeled command suites, and are automatically exposed via the standard DataLad commandline
  and Python API
- extensions can define `entry_points` for any number of additional metadata extractors
  that become automatically available to DataLad
- extensions can define `entry_points` for their test suites, such that the standard `datalad test`
  command will automatically run these tests in addition to the tests shipped with Datalad core
- extensions can ship additional dataset procedures by installing them into a
  directory ``resources/procedures`` underneath the extension module directory


Using an extension
==================

A :term:`DataLad extension` is a standard Python package. Beyond installation of the package there is
no additional setup required.


Writing your own extensions
===========================

A good starting point for implementing a new extension is the "helloworld" demo extension
available at https://github.com/datalad/datalad-extension-template. This repository can be cloned
and adjusted to suit one's needs. It includes:

- a basic Python package setup
- simple demo command implementation
- Travis test setup

A more complex extension setup can be seen in the DataLad Neuroimaging
extension: https://github.com/datalad/datalad-neuroimaging, including additional metadata extractors,
test suite registration, and a sphinx-based documentation setup for a DataLad extension.

As a DataLad extension is a standard Python package, an extension should declare
dependencies on an appropriate DataLad version, and possibly other extensions
via the standard mechanisms.
Acknowledgments
***************

DataLad development is being performed as part of a US-German collaboration in
computational neuroscience (CRCNS) project "DataGit: converging catalogues,
warehouses, and deployment logistics into a federated 'data distribution'"
(Halchenko_/Hanke_), co-funded by the US National Science Foundation (`NSF
1429999`_) and the German Federal Ministry of Education and Research (`BMBF
01GQ1411`_). Additional support is provided by the German federal state of
Saxony-Anhalt and the European Regional Development
Fund (ERDF), Project: `Center for Behavioral Brain Sciences`_, Imaging Platform

DataLad is built atop the git-annex_ software that is being developed and
maintained by `Joey Hess`_.

.. _Halchenko: http://haxbylab.dartmouth.edu/ppl/yarik.html
.. _Hanke: http://www.psychoinformatics.de
.. _NSF 1429999: http://www.nsf.gov/awardsearch/showAward?AWD_ID=1429999
.. _BMBF 01GQ1411: http://www.gesundheitsforschung-bmbf.de/de/2550.php
.. _Center for Behavioral Brain Sciences: http://cbbs.eu/en/
.. _git-annex: http://git-annex.branchable.com
.. _Joey Hess: https://joeyh.name
Background and motivation
*************************

Vision
======

Data is at the core of science, and unobstructed access promotes scientific
discovery through collaboration between data producers and consumers.  The last
years have seen dramatic improvements in availability of data resources for
collaborative research, and new data providers are becoming available all the
time.

However, despite the increased availability of data, their accessibility is far
from being optimal. Potential consumers of these public datasets have to
manually browse various disconnected warehouses with heterogeneous interfaces.
Once obtained, data is disconnected from its origin and data versioning is
often ad-hoc or completely absent. If data consumers can be reliably informed
about data updates at all, review of changes is difficult, and re-deployment is
tedious and error-prone. This leads to wasteful friction caused by outdated or
faulty data.

The vision for this project is to transform the state of data-sharing and
collaborative work by providing uniform access to available datasets --
independent of hosting solutions or authentication schemes -- with reliable
versioning and versatile deployment logistics. This is achieved by means of a
:term:`dataset` handle, a lightweight representation of a dataset
that is capable of tracking the identity and location of a dataset's content as
well as carry meta-data. Together with associated software tools, scientists
are able to obtain, use, extend, and share datasets (or parts thereof) in a
way that is traceable back to the original data producer and is therefore
capable of establishing a strong connection between data consumers and the
evolution of a dataset by future extension or error correction.

Moreover, DataLad aims to provide all tools necessary to create and publish
*data distributions* |---| an analog to software distributions or app-stores
that provide logistics middleware for software deployment. Scientific
communities can use these tools to gather, curate, and make publicly available
specialized collections of datasets for specific research topics or data
modalities. All of this is possible by leveraging existing data sharing
platforms and institutional resources without the need for funding extra
infrastructure of duplicate storage. Specifically, this project aims to provide
a comprehensive, extensible data distribution for neuroscientific datasets that
is kept up-to-date by an automated service.


Technological foundation: git-annex
===================================

The outlined task is not unique to the problem of data-sharing in science.
Logistical challenges such as delivering data, long-term storage and archiving,
identity tracking, and synchronization between multiple sites are rather
common. Consequently, solutions have been developed in other contexts that can
be adapted to benefit scientific data-sharing.

The closest match is the software tool git-annex_. It combines the features of
the distributed version control system (dVCS) Git_ |---| a technology that has
revolutionized collaborative software development -- with versatile data access
and delivery logistics. Git-annex was originally developed to address use cases
such as managing a collection of family pictures at home. With git-annex, any
family member can obtain an individual copy of such a picture library |---| the
:term:`annex`. The annex in this example is essentially an image repository
that presents individual pictures to users as files in a single directory
structure, even though the actual image file contents may be distributed across
multiple locations, including a home-server, cloud-storage, or even off-line
media such as external hard-drives.

Git-annex provides functionality to obtain file contents upon request and can
prompt users to make particular storage devices available when needed (e.g. a
backup hard-drive kept in a fire-proof compartment). Git-annex can also remove
files from a local copy of that image repository, for example to free up space
on a laptop, while ensuring a configurable level of data redundancy across all
known storage locations. Lastly, git-annex is able to synchronize the content
of multiple distributed copies of this image repository, for example in order
to incorporate images added with the git-annex on the laptop of another family
member. It is important to note that git-annex is agnostic of the actual file
types and is not limited to images.

We believe that the approach to data logistics taken by git-annex and the
functionality it is currently providing are an ideal middleware for scientific
data-sharing. Its data repository model :term:`annex` readily provides the
majority of principal features needed for a dataset handle such as history
recording, identity tracking, and item-based resource locators. Consequently,
instead of a from-scratch development, required features, such as dedicated
support for existing data-sharing portals and dataset meta-information, can be
added to a working solution that is already in production for several years.
As a result, DataLad focuses on the expansion of git-annex's functionality and
the development of tools that build atop Git and git-annex and enable the
creation, management, use, and publication of dataset handles and collections
thereof.

Objective
=========

Building atop git-annex, DataLad aims to provide a single, uniform interface to
access data from various data-sharing initiatives and data providers, and
functionality to create, deliver, update, and share datasets for individuals
and portal maintainers. As a command-line tool, it provides an abstraction
layer for the underlying Git-based middleware implementing the actual data
logistics, and serves as a foundation for other future user front-ends, such
as a web-interface.

.. |---| unicode:: U+02014 .. em dash

.. _Git: https://git-scm.com
.. _git-annex: http://git-annex.branchable.com
DataLad |---| data management and publication multitool
*******************************************************

Welcome to DataLad's **technical documentation**. Information here is targeting
software developers and is focused on the Python API and :term:`CLI`, as well
as software design, employed technologies, and key features.  Comprehensive
**user documentation** with information on installation, basic operation,
support, and (advanced) use case descriptions is available in the `DataLad
handbook <http://handbook.datalad.org>`_.

Content
^^^^^^^

.. toctree::
   :maxdepth: 1

   changelog
   acknowledgements
   publications

Concepts and technologies
=========================

.. toctree::
   :maxdepth: 2

   background
   related
   basics
   metadata
   customization
   design/index
   glossary

Commands and API
================

.. toctree::
   :maxdepth: 2

   cmdline
   modref
   config

Extension packages
==================

DataLad can be customized and additional functionality can be integrated via
extensions.  Each extension provides its own documentation:

- `Crawling web resources and automated data distributions <http://docs.datalad.org/projects/crawler>`_
- `Neuroimaging data and workflows <http://docs.datalad.org/projects/neuroimaging>`_
- `Containerized computational environments <http://docs.datalad.org/projects/container>`_
- `Advanced metadata tooling with JSON-LD reporting and additional metadata extractors <http://docs.datalad.org/projects/metalad>`_
- `Resources for working with the UKBiobank as a DataLad dataset <http://docs.datalad.org/projects/ukbiobank>`_
- `Deposit and retrieve DataLad datasets via the Open Science Framework <http://docs.datalad.org/projects/osf>`_
- `Functionality that has been phased out of the core package <http://docs.datalad.org/projects/deprecated>`_
- `Special interest functionality or drafts of future additions to DataLad proper <http://docs.datalad.org/projects/mihextras>`_

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

.. |---| unicode:: U+02014 .. em dash
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_modref:

***********************
Python module reference
***********************

This module reference extends the manual with a comprehensive overview of the
available functionality built into datalad.  Each module in the package is
documented by a general summary of its purpose and the list of classes and
functions it provides.


High-level user interface
=========================

Dataset operations
------------------

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   api.Dataset
   api.create
   api.create_sibling
   api.create_sibling_github
   api.create_sibling_gitlab
   api.create_sibling_gogs
   api.create_sibling_gitea
   api.create_sibling_gin
   api.drop
   api.get
   api.install
   api.push
   api.remove
   api.save
   api.update
   api.unlock

Metadata handling
-----------------

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   api.search
   api.metadata
   api.aggregate_metadata
   api.extract_metadata


Reproducible execution
----------------------

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   api.run
   api.rerun
   api.run_procedure


Plumbing commands
-----------------

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   api.clean
   api.clone
   api.copy_file
   api.create_test_dataset
   api.diff
   api.download_url
   api.sshrun
   api.siblings
   api.subdatasets

Miscellaneous commands
----------------------

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   api.add_archive_content
   api.test
   api.add_readme
   api.addurls
   api.check_dates
   api.export_archive
   api.export_to_figshare
   api.no_annex
   api.wtf

Support functionality
=====================

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   cmd
   consts
   log
   utils
   version
   support.gitrepo
   support.annexrepo
   support.archives
   customremotes.base
   customremotes.archives

Configuration management
========================

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   config

Test infrastructure
===================

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   tests.utils
   tests.utils_testrepos
   tests.heavyoutput

Command line interface infrastructure
=====================================

.. currentmodule:: datalad
.. autosummary::
   :toctree: generated

   cmdline.main
   cmdline.helpers
   cmdline.common_args
.. This file is auto-converted from CHANGELOG.md (make update-changelog) -- do not edit

Change log
**********
0.15.4 (Thu Dec 16 2021)
========================

Bug Fix
-------

-  BF: autorc - replace incorrect releaseTypes with “none”
   `#6320 <https://github.com/datalad/datalad/pull/6320>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Minor enhancement to CONTRIBUTING.md
   `#6309 <https://github.com/datalad/datalad/pull/6309>`__
   (`@bpoldrack <https://github.com/bpoldrack>`__)
-  UX: If a clean repo is dirty after a failed run, give clean-up hints
   `#6112 <https://github.com/datalad/datalad/pull/6112>`__
   (`@adswa <https://github.com/adswa>`__)
-  Stop using distutils
   `#6113 <https://github.com/datalad/datalad/pull/6113>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  BF: RIARemote - set UI backend to annex to make it interactive
   `#6287 <https://github.com/datalad/datalad/pull/6287>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__
   `@bpoldrack <https://github.com/bpoldrack>`__)
-  Fix invalid escape sequences
   `#6293 <https://github.com/datalad/datalad/pull/6293>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  CI: Update environment for windows CI builds
   `#6292 <https://github.com/datalad/datalad/pull/6292>`__
   (`@bpoldrack <https://github.com/bpoldrack>`__)
-  bump the python version used for mac os tests
   `#6288 <https://github.com/datalad/datalad/pull/6288>`__
   (`@christian-monch <https://github.com/christian-monch>`__
   `@bpoldrack <https://github.com/bpoldrack>`__)
-  ENH(UX): log a hint to use ulimit command in case of “Too long”
   exception `#6173 <https://github.com/datalad/datalad/pull/6173>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Report correct HTTP URL for RIA store content
   `#6091 <https://github.com/datalad/datalad/pull/6091>`__
   (`@mih <https://github.com/mih>`__)
-  BF: Don’t overwrite subdataset source candidates
   `#6168 <https://github.com/datalad/datalad/pull/6168>`__
   (`@bpoldrack <https://github.com/bpoldrack>`__)
-  Bump sphinx requirement to bypass readthedocs defaults
   `#6189 <https://github.com/datalad/datalad/pull/6189>`__
   (`@mih <https://github.com/mih>`__)
-  infra: Provide custom prefix to auto-related labels
   `#6151 <https://github.com/datalad/datalad/pull/6151>`__
   (`@adswa <https://github.com/adswa>`__)
-  Remove all usage of exc_str()
   `#6142 <https://github.com/datalad/datalad/pull/6142>`__
   (`@mih <https://github.com/mih>`__)
-  BF: obtain information about annex special remotes also from annex
   journal `#6135 <https://github.com/datalad/datalad/pull/6135>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__
   `@mih <https://github.com/mih>`__)
-  BF: clone tried to save new subdataset despite failing to clone
   `#6140 <https://github.com/datalad/datalad/pull/6140>`__
   (`@bpoldrack <https://github.com/bpoldrack>`__)

Tests
-----

-  RF+BF: use skip_if_no_module helper instead of try/except for libxmp
   and boto `#6148 <https://github.com/datalad/datalad/pull/6148>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  git://github.com -> https://github.com
   `#6134 <https://github.com/datalad/datalad/pull/6134>`__
   (`@mih <https://github.com/mih>`__)

Authors: 6
----------

-  Adina Wagner (`@adswa <https://github.com/adswa>`__)
-  Benjamin Poldrack (`@bpoldrack <https://github.com/bpoldrack>`__)
-  Christian Mnch
   (`@christian-monch <https://github.com/christian-monch>`__)
-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.15.3 (Sat Oct 30 2021)
========================

.. _bug-fix-1:

Bug Fix
-------

-  BF: Don’t make create-sibling recursive by default
   `#6116 <https://github.com/datalad/datalad/pull/6116>`__
   (`@adswa <https://github.com/adswa>`__)
-  BF: Add dashes to ‘force’ option in non-empty directory error message
   `#6078 <https://github.com/datalad/datalad/pull/6078>`__
   (`@DisasterMo <https://github.com/DisasterMo>`__)
-  DOC: Add supported URL types to download-url’s docstring
   `#6098 <https://github.com/datalad/datalad/pull/6098>`__
   (`@adswa <https://github.com/adswa>`__)
-  BF: Retain git-annex error messages & don’t show them if operation
   successful `#6070 <https://github.com/datalad/datalad/pull/6070>`__
   (`@DisasterMo <https://github.com/DisasterMo>`__)
-  Remove uses of ``__full_version__`` and ``datalad.version``
   `#6073 <https://github.com/datalad/datalad/pull/6073>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  BF: ORA shouldn’t crash while handling a failure
   `#6063 <https://github.com/datalad/datalad/pull/6063>`__
   (`@bpoldrack <https://github.com/bpoldrack>`__)
-  DOC: Refine –reckless docstring on usage and wording
   `#6043 <https://github.com/datalad/datalad/pull/6043>`__
   (`@adswa <https://github.com/adswa>`__)
-  BF: archives upon strip - use rmtree which retries etc instead of
   rmdir `#6064 <https://github.com/datalad/datalad/pull/6064>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: do not leave test in a tmp dir destined for removal
   `#6059 <https://github.com/datalad/datalad/pull/6059>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Next wave of exc_str() removals
   `#6022 <https://github.com/datalad/datalad/pull/6022>`__
   (`@mih <https://github.com/mih>`__)

Pushed to ``maint``
-------------------

-  CI: Enable new codecov uploader in Appveyor CI
   (`@adswa <https://github.com/adswa>`__)

Internal
--------

-  UX: Log clone-candidate number and URLs
   `#6092 <https://github.com/datalad/datalad/pull/6092>`__
   (`@adswa <https://github.com/adswa>`__)
-  UX/ENH: Disable reporting, and don’t do superfluous internal
   subdatasets calls
   `#6094 <https://github.com/datalad/datalad/pull/6094>`__
   (`@adswa <https://github.com/adswa>`__)
-  Update codecov action to v2
   `#6072 <https://github.com/datalad/datalad/pull/6072>`__
   (`@jwodder <https://github.com/jwodder>`__)

Documentation
-------------

-  Design document on URL substitution feature
   `#6065 <https://github.com/datalad/datalad/pull/6065>`__
   (`@mih <https://github.com/mih>`__)

.. _tests-1:

Tests
-----

-  BF(TST): remove reuse of the same tape across unrelated tests
   `#6127 <https://github.com/datalad/datalad/pull/6127>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Fail Travis tests on deprecation warnings
   `#6074 <https://github.com/datalad/datalad/pull/6074>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  Ux get result handling broken
   `#6052 <https://github.com/datalad/datalad/pull/6052>`__
   (`@christian-monch <https://github.com/christian-monch>`__)
-  enable metalad tests again
   `#6060 <https://github.com/datalad/datalad/pull/6060>`__
   (`@christian-monch <https://github.com/christian-monch>`__)

Authors: 7
----------

-  Adina Wagner (`@adswa <https://github.com/adswa>`__)
-  Benjamin Poldrack (`@bpoldrack <https://github.com/bpoldrack>`__)
-  Christian Mnch
   (`@christian-monch <https://github.com/christian-monch>`__)
-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Michael Burgardt (`@DisasterMo <https://github.com/DisasterMo>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.15.2 (Wed Oct 06 2021)
========================

.. _bug-fix-2:

Bug Fix
-------

-  BF: Don’t suppress datalad subdatasets output
   `#6035 <https://github.com/datalad/datalad/pull/6035>`__
   (`@DisasterMo <https://github.com/DisasterMo>`__
   `@mih <https://github.com/mih>`__)
-  Honor datalad.runtime.use-patool if set regardless of OS (was Windows
   only) `#6033 <https://github.com/datalad/datalad/pull/6033>`__
   (`@mih <https://github.com/mih>`__)
-  Discontinue usage of deprecated (public) helper
   `#6032 <https://github.com/datalad/datalad/pull/6032>`__
   (`@mih <https://github.com/mih>`__)
-  BF: ProgressHandler - close the other handler if was specified
   `#6020 <https://github.com/datalad/datalad/pull/6020>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  UX: Report GitLab weburl of freshly created projects in the result
   `#6017 <https://github.com/datalad/datalad/pull/6017>`__
   (`@adswa <https://github.com/adswa>`__)
-  Ensure there’s a blank line between the class ``__doc__`` and
   “Parameters” in ``build_doc`` docstrings
   `#6004 <https://github.com/datalad/datalad/pull/6004>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  Large code-reorganization of everything runner-related
   `#6008 <https://github.com/datalad/datalad/pull/6008>`__
   (`@mih <https://github.com/mih>`__)
-  Discontinue exc_str() in all modern parts of the code base
   `#6007 <https://github.com/datalad/datalad/pull/6007>`__
   (`@mih <https://github.com/mih>`__)

.. _tests-2:

Tests
-----

-  TST: Add test to ensure functionality with subdatasets starting with
   a hyphen (-) `#6042 <https://github.com/datalad/datalad/pull/6042>`__
   (`@DisasterMo <https://github.com/DisasterMo>`__)
-  BF(TST): filter away warning from coverage from analysis of stderr of
   –help `#6028 <https://github.com/datalad/datalad/pull/6028>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: disable outdated SSL root certificate breaking chain on
   older/buggy clients
   `#6027 <https://github.com/datalad/datalad/pull/6027>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: start global test_http_server only if not running already
   `#6023 <https://github.com/datalad/datalad/pull/6023>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

Authors: 5
----------

-  Adina Wagner (`@adswa <https://github.com/adswa>`__)
-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Michael Burgardt (`@DisasterMo <https://github.com/DisasterMo>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.15.1 (Fri Sep 24 2021)
========================

.. _bug-fix-3:

Bug Fix
-------

-  BF: downloader - fail to download even on non-crippled FS if symlink
   exists `#5991 <https://github.com/datalad/datalad/pull/5991>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH: import datalad.api to bind extensions methods for discovery of
   dataset methods
   `#5999 <https://github.com/datalad/datalad/pull/5999>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Restructure cmdline API presentation
   `#5988 <https://github.com/datalad/datalad/pull/5988>`__
   (`@mih <https://github.com/mih>`__)
-  Close file descriptors after process exit
   `#5983 <https://github.com/datalad/datalad/pull/5983>`__
   (`@mih <https://github.com/mih>`__)

.. _pushed-to-maint-1:

Pushed to ``maint``
-------------------

-  Discontinue testing of hirni extension
   (`@mih <https://github.com/mih>`__)

.. _internal-1:

Internal
--------

-  Add debugging information to release step
   `#5980 <https://github.com/datalad/datalad/pull/5980>`__
   (`@jwodder <https://github.com/jwodder>`__)

.. _documentation-1:

Documentation
-------------

-  Coarse description of the credential subsystem’s functionality
   `#5998 <https://github.com/datalad/datalad/pull/5998>`__
   (`@mih <https://github.com/mih>`__)

.. _tests-3:

Tests
-----

-  BF(TST): use sys.executable, mark test_ria_basics.test_url_keys as
   requiring network
   `#5986 <https://github.com/datalad/datalad/pull/5986>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

Authors: 3
----------

-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.15.0 (Tue Sep 14 2021) – We miss you Kyle!
============================================

Enhancements and new features
-----------------------------

-  Command execution is now performed by a new ``Runner`` implementation
   that is no longer based on the ``asyncio`` framework, which was found
   to exhibit fragile performance in interaction with other
   ``asyncio``-using code, such as Jupyter notebooks. The new
   implementation is based on threads. It also supports the
   specification of “protocols” that were introduced with the switch to
   the ``asyncio`` implementation in 0.14.0.
   (`#5667 <https://github.com/datalad/datalad/issues/5667>`__)

-  ``clone`` now supports arbitrary URL transformations based on regular
   expressions. One or more transformation steps can be defined via
   ``datalad.clone.url-substitute.<label>`` configuration settings. The
   feature can be (and is now) used to support convenience mappings,
   such as ``https://osf.io/q8xnk/`` (displayed in a browser window) to
   ``osf://q8xnk`` (clonable via the ``datalad-osf`` extension.
   (`#5749 <https://github.com/datalad/datalad/issues/5749>`__)

-  Homogenize SSH use and configurability between DataLad and git-annex,
   by instructing git-annex to use DataLad’s ``sshrun`` for SSH calls
   (instead of SSH directly).
   (`#5389 <https://github.com/datalad/datalad/issues/5389>`__)

-  The ORA special remote has received several new features:

   -  It now support a ``push-url`` setting as an alternative to ``url``
      for write access. An analog parameter was also added to
      ``create-sibling-ria``.
      (`#5420 <https://github.com/datalad/datalad/issues/5420>`__,
      `#5428 <https://github.com/datalad/datalad/issues/5428>`__)

   -  Access of RIA stores now performs homogeneous availability checks,
      regardless of access protocol. Before, broken HTTP-based access
      due to misspecified URLs could have gone unnoticed.
      (`#5459 <https://github.com/datalad/datalad/issues/5459>`__,
      `#5672 <https://github.com/datalad/datalad/issues/5672>`__)

   -  Error reporting was introduce to inform about undesirable
      conditions in remote RIA stores.
      (`#5683 <https://github.com/datalad/datalad/issues/5683>`__)

-  ``create-sibling-ria`` now supports ``--alias`` for the specification
   of a convenience dataset alias name in a RIA store.
   (`#5592 <https://github.com/datalad/datalad/issues/5592>`__)

-  Analog to ``git commit``, ``save`` now features an ``--amend`` mode
   to support incremental updates of a dataset state.
   (`#5430 <https://github.com/datalad/datalad/issues/5430>`__)

-  ``run`` now supports a dry-run mode that can be used to inspect the
   result of parameter expansion on the effective command to ease the
   composition of more complicated command lines.
   (`#5539 <https://github.com/datalad/datalad/issues/5539>`__)

-  ``run`` now supports a ``--assume-ready`` switch to avoid the
   (possibly expensive) preparation of inputs and outputs with large
   datasets that have already been readied through other means.
   (`#5431 <https://github.com/datalad/datalad/issues/5431>`__)

-  ``update`` now features ``--how`` and ``--how-subds`` parameters to
   configure how an update shall be performed. Supported modes are
   ``fetch`` (unchanged default), and ``merge`` (previously also
   possible via ``--merge``), but also new strategies like ``reset`` or
   ``checkout``.
   (`#5534 <https://github.com/datalad/datalad/issues/5534>`__)

-  ``update`` has a new ``--follow=parentds-lazy`` mode that only
   performs a fetch operation in subdatasets when the desired commit is
   not yet present. During recursive updates involving many subdatasets
   this can substantially speed up performance.
   (`#5474 <https://github.com/datalad/datalad/issues/5474>`__)

-  DataLad’s command line API can now report the version for individual
   commands via ``datalad <cmd> --version``. The output has been
   homogenized to ``<providing package> <version>``.
   (`#5543 <https://github.com/datalad/datalad/issues/5543>`__)

-  ``create-sibling`` now logs information on an auto-generated sibling
   name, in the case that no ``--name/-s`` was provided.
   (`#5550 <https://github.com/datalad/datalad/issues/5550>`__)

-  ``create-sibling-github`` has been updated to emit result records
   like any standard DataLad command. Previously it was implemented as a
   “plugin”, which did not support all standard API parameters.
   (`#5551 <https://github.com/datalad/datalad/issues/5551>`__)

-  ``copy-file`` now also works with content-less files in datasets on
   crippled filesystems (adjusted mode), when a recent enough git-annex
   (8.20210428 or later) is available.
   (`#5630 <https://github.com/datalad/datalad/issues/5630>`__)

-  ``addurls`` can now be instructed how to behave in the event of file
   name collision via a new parameter ``--on-collision``.
   (`#5675 <https://github.com/datalad/datalad/issues/5675>`__)

-  ``addurls`` reporting now informs which particular subdatasets were
   created. (`#5689 <https://github.com/datalad/datalad/issues/5689>`__)

-  Credentials can now be provided or overwritten via all means
   supported by ``ConfigManager``. Importantly,
   ``datalad.credential.<name>.<field>`` configuration settings and
   analog specification via environment variables are now supported
   (rather than custom environment variables only). Previous
   specification methods are still supported too.
   (`#5680 <https://github.com/datalad/datalad/issues/5680>`__)

-  A new ``datalad.credentials.force-ask`` configuration flag can now be
   used to force re-entry of already known credentials. This simplifies
   credential updates without having to use an approach native to
   individual credential stores.
   (`#5777 <https://github.com/datalad/datalad/issues/5777>`__)

-  Suppression of rendering repeated similar results is now configurable
   via the configuration switches
   ``datalad.ui.suppress-similar-results`` (bool), and
   ``datalad.ui.suppress-similar-results-threshold`` (int).
   (`#5681 <https://github.com/datalad/datalad/issues/5681>`__)

-  The performance of ``status`` and similar functionality when
   determining local file availability has been improved.
   (`#5692 <https://github.com/datalad/datalad/issues/5692>`__)

-  ``push`` now renders a result summary on completion.
   (`#5696 <https://github.com/datalad/datalad/issues/5696>`__)

-  A dedicated info log message indicates when dataset repositories are
   subjected to an annex version upgrade.
   (`#5698 <https://github.com/datalad/datalad/issues/5698>`__)

-  Error reporting improvements:

   -  The ``NoDatasetFound`` exception now provides information for
      which purpose a dataset is required.
      (`#5708 <https://github.com/datalad/datalad/issues/5708>`__)

   -  Wording of the ``MissingExternalDependeny`` error was rephrased to
      account for cases of non-functional installations.
      (`#5803 <https://github.com/datalad/datalad/issues/5803>`__)

   -  ``push`` reports when a ``--to`` parameter specification was
      (likely) forgotten.
      (`#5726 <https://github.com/datalad/datalad/issues/5726>`__)

   -  Detailed information is now given when DataLad fails to obtain a
      lock for credential entry in a timely fashion. Previously only a
      generic debug log message was emitted.
      (`#5884 <https://github.com/datalad/datalad/issues/5884>`__)

   -  Clarified error message when ``create-sibling-gitlab`` was called
      without ``--project``.
      (`#5907 <https://github.com/datalad/datalad/issues/5907>`__)

-  ``add-readme`` now provides a README template with more information
   on the nature and use of DataLad datasets. A README file is no longer
   annex’ed by default, but can be using the new ``--annex`` switch.
   ([#5723][], [#5725][])

-  ``clean`` now supports a ``--dry-run`` mode to inform about cleanable
   content. (`#5738 <https://github.com/datalad/datalad/issues/5738>`__)

-  A new configuration setting ``datalad.locations.locks`` can be used
   to control the placement of lock files.
   (`#5740 <https://github.com/datalad/datalad/issues/5740>`__)

-  ``wtf`` now also reports branch names and states.
   (`#5804 <https://github.com/datalad/datalad/issues/5804>`__)

-  ``AnnexRepo.whereis()`` now supports batch mode.
   (`#5533 <https://github.com/datalad/datalad/issues/5533>`__)

Deprecations and removals
~~~~~~~~~~~~~~~~~~~~~~~~~

-  The minimum supported git-annex version is now 8.20200309.
   (`#5512 <https://github.com/datalad/datalad/issues/5512>`__)

-  ORA special remote configuration items ``ssh-host``, and
   ``base-path`` are deprecated. They are completely replaced by
   ``ria+<protocol>://`` URL specifications.
   (`#5425 <https://github.com/datalad/datalad/issues/5425>`__)

-  The deprecated ``no_annex`` parameter of ``create()`` was removed
   from the Python API.
   (`#5441 <https://github.com/datalad/datalad/issues/5441>`__)

-  The unused ``GitRepo.pull()`` method has been removed.
   (`#5558 <https://github.com/datalad/datalad/issues/5558>`__)

-  Residual support for “plugins” (a mechanism used before DataLad
   supported extensions) was removed. This includes the configuration
   switches ``datalad.locations.{system,user}-plugins``.
   (`#5554 <https://github.com/datalad/datalad/issues/5554>`__,
   `#5564 <https://github.com/datalad/datalad/issues/5564>`__)

-  Several features and comments have been moved to the
   ``datalad-deprecated`` package. This package must now be installed to
   be able to use keep using this functionality.

   -  The ``publish`` command. Use ``push`` instead.
      (`#5837 <https://github.com/datalad/datalad/issues/5837>`__)

   -  The ``ls`` command.
      (`#5569 <https://github.com/datalad/datalad/issues/5569>`__)

   -  The web UI that is deployable via ``datalad create-sibling --ui``.
      (`#5555 <https://github.com/datalad/datalad/issues/5555>`__)

   -  The “automagic IO” feature.
      (`#5577 <https://github.com/datalad/datalad/issues/5577>`__)

-  ``AnnexRepo.copy_to()`` has been deprecated. The ``push`` command
   should be used instead.
   (`#5560 <https://github.com/datalad/datalad/issues/5560>`__)

-  ``AnnexRepo.sync()`` has been deprecated.
   ``AnnexRepo.call_annex(['sync', ...])`` should be used instead.
   (`#5461 <https://github.com/datalad/datalad/issues/5461>`__)

-  All ``GitRepo.*_submodule()`` methods have been deprecated and will
   be removed in a future release.
   (`#5559 <https://github.com/datalad/datalad/issues/5559>`__)

-  ``create-sibling-github``\ ’s ``--dryrun`` switch was deprecated, use
   ``--dry-run`` instead.
   (`#5551 <https://github.com/datalad/datalad/issues/5551>`__)

-  The ``datalad --pbs-runner`` option has been deprecated, use
   ``condor_run`` (or similar) instead.
   (`#5956 <https://github.com/datalad/datalad/issues/5956>`__)

Fixes
-----

-  Prevent invalid declaration of a publication dependencies for
   ‘origin’ on any auto-detected ORA special remotes, when cloing from a
   RIA store. An ORA remote is now checked whether it actually points to
   the RIA store the clone was made from.
   (`#5415 <https://github.com/datalad/datalad/issues/5415>`__)

-  The ORA special remote implementation has received several fixes:

   -  It can now handle HTTP redirects.
      (`#5792 <https://github.com/datalad/datalad/issues/5792>`__)

   -  Prevents failure when URL-type annex keys contain the ‘/’
      character.
      (`#5823 <https://github.com/datalad/datalad/issues/5823>`__)

   -  Properly support the specification of usernames, passwords and
      ports in ``ria+<protocol>://`` URLs.
      (`#5902 <https://github.com/datalad/datalad/issues/5902>`__)

-  It is now possible to specifically select the default (or generic)
   result renderer via ``datalad -f default`` and with that override a
   ``tailored`` result renderer that may be preconfigured for a
   particular command.
   (`#5476 <https://github.com/datalad/datalad/issues/5476>`__)

-  Starting with 0.14.0, original URLs given to ``clone`` were recorded
   in a subdataset record. This was initially done in a second commit,
   leading to inflation of commits and slowdown in superdatasets with
   many subdatasets. Such subdataset record annotation is now collapsed
   into a single commits.
   (`#5480 <https://github.com/datalad/datalad/issues/5480>`__)

-  ``run`` now longer removes leading empty directories as part of the
   output preparation. This was surprising behavior for commands that do
   not ensure on their own that output directories exist.
   (`#5492 <https://github.com/datalad/datalad/issues/5492>`__)

-  A potentially existing ``message`` property is no longer removed when
   using the ``json`` or ``json_pp`` result renderer to avoid undesired
   withholding of relevant information.
   (`#5536 <https://github.com/datalad/datalad/issues/5536>`__)

-  ``subdatasets`` now reports ``state=present``, rather than
   ``state=clean``, for installed subdatasets to complement
   ``state=absent`` reports for uninstalled dataset.
   (`#5655 <https://github.com/datalad/datalad/issues/5655>`__)

-  ``create-sibling-ria`` now executes commands with a consistent
   environment setup that matches all other command execution in other
   DataLad commands.
   (`#5682 <https://github.com/datalad/datalad/issues/5682>`__)

-  ``save`` no longer saves unspecified subdatasets when called with an
   explicit path (list). The fix required a behavior change of
   ``GitRepo.get_content_info()`` in its interpretation of ``None``
   vs. \ ``[]`` path argument values that now aligns the behavior of
   ``GitRepo.diff|status()`` with their respective documentation.
   (`#5693 <https://github.com/datalad/datalad/issues/5693>`__)

-  ``get`` now prefers the location of a subdatasets that is recorded in
   a superdataset’s ``.gitmodules`` record. Previously, DataLad tried to
   obtain a subdataset from an assumed checkout of the superdataset’s
   origin. This new default order is (re-)configurable via the
   ``datalad.get.subdataset-source-candidate-<priority-label>``
   configuration mechanism.
   (`#5760 <https://github.com/datalad/datalad/issues/5760>`__)

-  ``create-sibling-gitlab`` no longer skips the root dataset when ``.``
   is given as a path.
   (`#5789 <https://github.com/datalad/datalad/issues/5789>`__)

-  ``siblings`` now rejects a value given to ``--as-common-datasrc``
   that clashes with the respective Git remote.
   (`#5805 <https://github.com/datalad/datalad/issues/5805>`__)

-  The usage synopsis reported by ``siblings`` now lists all supported
   actions. (`#5913 <https://github.com/datalad/datalad/issues/5913>`__)

-  ``siblings`` now renders non-ok results to avoid silent failure.
   (`#5915 <https://github.com/datalad/datalad/issues/5915>`__)

-  ``.gitattribute`` file manipulations no longer leave the file without
   a trailing newline.
   (`#5847 <https://github.com/datalad/datalad/issues/5847>`__)

-  Prevent crash when trying to delete a non-existing keyring credential
   field. (`#5892 <https://github.com/datalad/datalad/issues/5892>`__)

-  git-annex is no longer called with an unconditional ``annex.retry=3``
   configuration. Instead, this parameterization is now limited to
   ``annex get`` and ``annex copy`` calls.
   (`#5904 <https://github.com/datalad/datalad/issues/5904>`__)

.. _tests-4:

Tests
-----

-  ``file://`` URLs are no longer the predominant test case for
   ``AnnexRepo`` functionality. A built-in HTTP server now used in most
   cases. (`#5332 <https://github.com/datalad/datalad/issues/5332>`__)

--------------

0.14.8 (Sun Sep 12 2021)
========================

.. _bug-fix-4:

Bug Fix
-------

-  BF: add-archive-content on .xz and other non-.gz stream compressed
   files `#5930 <https://github.com/datalad/datalad/pull/5930>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF(UX): do not keep logging ERROR possibly present in progress
   records `#5936 <https://github.com/datalad/datalad/pull/5936>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Annotate datalad_core as not needing actual data – just uses annex
   whereis `#5971 <https://github.com/datalad/datalad/pull/5971>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: limit CMD_MAX_ARG if obnoxious value is encountered.
   `#5945 <https://github.com/datalad/datalad/pull/5945>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Download session/credentials locking – inform user if locking is
   “failing” to be obtained, fail upon ~5min timeout
   `#5884 <https://github.com/datalad/datalad/pull/5884>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Render siblings()’s non-ok results with the default renderer
   `#5915 <https://github.com/datalad/datalad/pull/5915>`__
   (`@mih <https://github.com/mih>`__)
-  BF: do not crash, just skip whenever trying to delete non existing
   field in the underlying keyring
   `#5892 <https://github.com/datalad/datalad/pull/5892>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Fix argument-spec for ``siblings`` and improve usage synopsis
   `#5913 <https://github.com/datalad/datalad/pull/5913>`__
   (`@mih <https://github.com/mih>`__)
-  Clarify error message re unspecified gitlab project
   `#5907 <https://github.com/datalad/datalad/pull/5907>`__
   (`@mih <https://github.com/mih>`__)
-  Support username, password and port specification in RIA URLs
   `#5902 <https://github.com/datalad/datalad/pull/5902>`__
   (`@mih <https://github.com/mih>`__)
-  BF: take path from SSHRI, test URLs not only on Windows
   `#5881 <https://github.com/datalad/datalad/pull/5881>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH(UX): warn user if keyring returned a “null” keyring
   `#5875 <https://github.com/datalad/datalad/pull/5875>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH(UX): state original purpose in NoDatasetFound exception + detail
   it for get `#5708 <https://github.com/datalad/datalad/pull/5708>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

.. _pushed-to-maint-2:

Pushed to ``maint``
-------------------

-  Merge branch ‘bf-http-headers-agent’ into maint
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  RF(BF?)+DOC: provide User-Agent to entire session headers + use those
   if provided (`@yarikoptic <https://github.com/yarikoptic>`__)

.. _internal-2:

Internal
--------

-  Pass ``--no-changelog`` to ``auto shipit`` if changelog already has
   entry `#5952 <https://github.com/datalad/datalad/pull/5952>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  Add isort config to match current convention + run isort via
   pre-commit (if configured)
   `#5923 <https://github.com/datalad/datalad/pull/5923>`__
   (`@jwodder <https://github.com/jwodder>`__)
-  .travis.yml: use python -m {nose,coverage} invocations, and always
   show combined report
   `#5888 <https://github.com/datalad/datalad/pull/5888>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Add project URLs into the package metadata for convenience links on
   Pypi `#5866 <https://github.com/datalad/datalad/pull/5866>`__
   (`@adswa <https://github.com/adswa>`__
   `@yarikoptic <https://github.com/yarikoptic>`__)

.. _tests-5:

Tests
-----

-  BF: do use OBSCURE_FILENAME instead of hardcoded unicode
   `#5944 <https://github.com/datalad/datalad/pull/5944>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF(TST): Skip testing for having PID listed if no psutil
   `#5920 <https://github.com/datalad/datalad/pull/5920>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF(TST): Boost version of git-annex to 8.20201129 to test an error
   message `#5894 <https://github.com/datalad/datalad/pull/5894>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

Authors: 4
----------

-  Adina Wagner (`@adswa <https://github.com/adswa>`__)
-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.14.7 (Tue Aug 03 2021)
========================

.. _bug-fix-5:

Bug Fix
-------

-  UX: When two or more clone URL templates are found, error out more
   gracefully `#5839 <https://github.com/datalad/datalad/pull/5839>`__
   (`@adswa <https://github.com/adswa>`__)
-  BF: http_auth - follow redirect (just 1) to re-authenticate after
   initial attempt
   `#5852 <https://github.com/datalad/datalad/pull/5852>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  addurls Formatter - provide value repr in exception
   `#5850 <https://github.com/datalad/datalad/pull/5850>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH: allow for “patch” level semver for “master” branch
   `#5839 <https://github.com/datalad/datalad/pull/5839>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: Report info from annex JSON error message in CommandError
   `#5809 <https://github.com/datalad/datalad/pull/5809>`__
   (`@mih <https://github.com/mih>`__)
-  RF(TST): do not test for no EASY and pkg_resources in shims
   `#5817 <https://github.com/datalad/datalad/pull/5817>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  http downloaders: Provide custom informative User-Agent, do not claim
   to be “Authenticated access”
   `#5802 <https://github.com/datalad/datalad/pull/5802>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH(UX,DX): inform user with a warning if version is 0+unknown
   `#5787 <https://github.com/datalad/datalad/pull/5787>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  shell-completion: add argcomplete to ‘misc’ extra_depends, log an
   ERROR if argcomplete fails to import
   `#5781 <https://github.com/datalad/datalad/pull/5781>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH (UX): add python-gitlab dependency
   `#5776 <https://github.com/datalad/datalad/pull/5776>`__
   (s.heunis@fz-juelich.de)

.. _internal-3:

Internal
--------

-  BF: Fix reported paths in ORA remote
   `#5821 <https://github.com/datalad/datalad/pull/5821>`__
   (`@adswa <https://github.com/adswa>`__)
-  BF: import importlib.metadata not importlib_metadata whenever
   available `#5818 <https://github.com/datalad/datalad/pull/5818>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

.. _tests-6:

Tests
-----

-  TST: set –allow-unrelated-histories in the mk_push_target setup for
   Windows `#5855 <https://github.com/datalad/datalad/pull/5855>`__
   (`@adswa <https://github.com/adswa>`__)
-  Tests: Allow for version to contain + as a separator and provide more
   information for version related comparisons
   `#5786 <https://github.com/datalad/datalad/pull/5786>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)

.. _authors-4-1:

Authors: 4
----------

-  Adina Wagner (`@adswa <https://github.com/adswa>`__)
-  Michael Hanke (`@mih <https://github.com/mih>`__)
-  Stephan Heunis (`@jsheunis <https://github.com/jsheunis>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.14.6 (Sun Jun 27 2021)
========================

.. _internal-4:

Internal
--------

-  BF: update changelog conversion from .md to .rst (for sphinx)
   `#5757 <https://github.com/datalad/datalad/pull/5757>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__
   `@jwodder <https://github.com/jwodder>`__)

Authors: 2
----------

-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.14.5 (Mon Jun 21 2021)
========================

.. _bug-fix-6:

Bug Fix
-------

-  BF(TST): parallel - take longer for producer to produce
   `#5747 <https://github.com/datalad/datalad/pull/5747>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  add –on-failure default value and document it
   `#5690 <https://github.com/datalad/datalad/pull/5690>`__
   (`@christian-monch <https://github.com/christian-monch>`__
   `@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH: harmonize “purpose” statements to imperative form
   `#5733 <https://github.com/datalad/datalad/pull/5733>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  ENH(TST): populate heavy tree with 100 unique keys (not just 1) among
   10,000 `#5734 <https://github.com/datalad/datalad/pull/5734>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: do not use .acquired - just get state from acquire()
   `#5718 <https://github.com/datalad/datalad/pull/5718>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  BF: account for annex now “scanning for annexed” instead of
   “unlocked” files
   `#5705 <https://github.com/datalad/datalad/pull/5705>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  interface: Don’t repeat custom summary for non-generator results
   `#5688 <https://github.com/datalad/datalad/pull/5688>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  RF: just pip install datalad-installer
   `#5676 <https://github.com/datalad/datalad/pull/5676>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  DOC: addurls.extract: Drop mention of removed ‘stream’ parameter
   `#5690 <https://github.com/datalad/datalad/pull/5690>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  Merge pull request #5674 from kyleam/test-addurls-copy-fix
   `#5674 <https://github.com/datalad/datalad/pull/5674>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  Merge pull request #5663 from kyleam/status-ds-equal-path
   `#5663 <https://github.com/datalad/datalad/pull/5663>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  Merge pull request #5671 from kyleam/update-fetch-fail
   `#5671 <https://github.com/datalad/datalad/pull/5671>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  BF: update: Honor –on-failure if fetch fails
   `#5671 <https://github.com/datalad/datalad/pull/5671>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  RF: update: Avoid fetch’s deprecated kwargs
   `#5671 <https://github.com/datalad/datalad/pull/5671>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  CLN: update: Drop an unused import
   `#5671 <https://github.com/datalad/datalad/pull/5671>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  Merge pull request #5664 from kyleam/addurls-better-url-parts-error
   `#5664 <https://github.com/datalad/datalad/pull/5664>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  Merge pull request #5661 from kyleam/sphinx-fix-plugin-refs
   `#5661 <https://github.com/datalad/datalad/pull/5661>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  BF: status: Provide special treatment of “this dataset” path
   `#5663 <https://github.com/datalad/datalad/pull/5663>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  BF: addurls: Provide better placeholder error for special keys
   `#5664 <https://github.com/datalad/datalad/pull/5664>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  RF: addurls: Simply construction of placeholder exception message
   `#5664 <https://github.com/datalad/datalad/pull/5664>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  RF: addurls._get_placeholder_exception: Rename a parameter
   `#5664 <https://github.com/datalad/datalad/pull/5664>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  RF: status: Avoid repeated Dataset.path access
   `#5663 <https://github.com/datalad/datalad/pull/5663>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  DOC: Reference plugins via datalad.api
   `#5661 <https://github.com/datalad/datalad/pull/5661>`__
   (`@kyleam <https://github.com/kyleam>`__)
-  download-url: Set up datalad special remote if needed
   `#5648 <https://github.com/datalad/datalad/pull/5648>`__
   (`@kyleam <https://github.com/kyleam>`__
   `@yarikoptic <https://github.com/yarikoptic>`__)

.. _pushed-to-maint-3:

Pushed to ``maint``
-------------------

-  MNT: Post-release dance (`@kyleam <https://github.com/kyleam>`__)

.. _internal-5:

Internal
--------

-  Switch to versioneer and auto
   `#5669 <https://github.com/datalad/datalad/pull/5669>`__
   (`@jwodder <https://github.com/jwodder>`__
   `@yarikoptic <https://github.com/yarikoptic>`__)
-  MNT: setup.py: Temporarily avoid Sphinx 4
   `#5649 <https://github.com/datalad/datalad/pull/5649>`__
   (`@kyleam <https://github.com/kyleam>`__)

.. _tests-7:

Tests
-----

-  BF(TST): skip testing for showing “Scanning for …” since not shown if
   too quick `#5727 <https://github.com/datalad/datalad/pull/5727>`__
   (`@yarikoptic <https://github.com/yarikoptic>`__)
-  Revert “TST: test_partial_unlocked: Document and avoid recent
   git-annex failure”
   `#5651 <https://github.com/datalad/datalad/pull/5651>`__
   (`@kyleam <https://github.com/kyleam>`__)

.. _authors-4-2:

Authors: 4
----------

-  Christian Mnch
   (`@christian-monch <https://github.com/christian-monch>`__)
-  John T. Wodder II (`@jwodder <https://github.com/jwodder>`__)
-  Kyle Meyer (`@kyleam <https://github.com/kyleam>`__)
-  Yaroslav Halchenko (`@yarikoptic <https://github.com/yarikoptic>`__)

--------------

0.14.4 (May 10, 2021) – .
=========================

.. _fixes-1:

Fixes
-----

-  Following an internal call to ``git-clone``,
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   assumed that the remote name was “origin”, but this may not be the
   case if ``clone.defaultRemoteName`` is configured (available as of
   Git 2.30).
   (`#5572 <https://github.com/datalad/datalad/issues/5572>`__)

-  Several test fixes, including updates for changes in git-annex.
   (`#5612 <https://github.com/datalad/datalad/issues/5612>`__)
   (`#5632 <https://github.com/datalad/datalad/issues/5632>`__)
   (`#5639 <https://github.com/datalad/datalad/issues/5639>`__)

0.14.3 (April 28, 2021) – .
===========================

.. _fixes-2:

Fixes
-----

-  For outputs that include a glob,
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   didn’t re-glob after executing the command, which is necessary to
   catch changes if ``--explicit`` or ``--expand={outputs,both}`` is
   specified.
   (`#5594 <https://github.com/datalad/datalad/issues/5594>`__)

-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   now gives an error result rather than a warning when an input glob
   doesn’t match.
   (`#5594 <https://github.com/datalad/datalad/issues/5594>`__)

-  The procedure for creating a RIA store checks for an existing
   ria-layout-version file and makes sure its version matches the
   desired version. This check wasn’t done correctly for SSH hosts.
   (`#5607 <https://github.com/datalad/datalad/issues/5607>`__)

-  A helper for transforming git-annex JSON records into DataLad results
   didn’t account for the unusual case where the git-annex record
   doesn’t have a “file” key.
   (`#5580 <https://github.com/datalad/datalad/issues/5580>`__)

-  The test suite required updates for recent changes in PyGithub and
   git-annex.
   (`#5603 <https://github.com/datalad/datalad/issues/5603>`__)
   (`#5609 <https://github.com/datalad/datalad/issues/5609>`__)

.. _enhancements-and-new-features-1:

Enhancements and new features
-----------------------------

-  The DataLad source repository has long had a tools/cmdline-completion
   helper. This functionality is now exposed as a command,
   ``datalad shell-completion``.
   (`#5544 <https://github.com/datalad/datalad/issues/5544>`__)

0.14.2 (April 14, 2021) – .
===========================

.. _fixes-3:

Fixes
-----

-  `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   now works bottom-up, pushing submodules first so that hooks on the
   remote can aggregate updated subdataset information.
   (`#5416 <https://github.com/datalad/datalad/issues/5416>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   didn’t ensure that the configuration of subdatasets was reloaded.
   (`#5552 <https://github.com/datalad/datalad/issues/5552>`__)

0.14.1 (April 01, 2021) – .
===========================

.. _fixes-4:

Fixes
-----

-  The recent default branch changes on GitHub’s side can lead to
   “git-annex” being selected over “master” as the default branch on
   GitHub when setting up a sibling with
   `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__.
   To work around this, the current branch is now pushed first.
   (`#5010 <https://github.com/datalad/datalad/issues/5010>`__)

-  The logic for reading in a JSON line from git-annex failed if the
   response exceeded the buffer size (256 KB on \*nix systems).

-  Calling
   `unlock <http://datalad.readthedocs.io/en/latest/generated/man/datalad-unlock.html>`__
   with a path of “.” from within an untracked subdataset incorrectly
   aborted, complaining that the “dataset containing given paths is not
   underneath the reference dataset”.
   (`#5458 <https://github.com/datalad/datalad/issues/5458>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   didn’t account for the possibility of multiple accessible ORA remotes
   or the fact that none of them may be associated with the RIA store
   being cloned.
   (`#5488 <https://github.com/datalad/datalad/issues/5488>`__)

-  `create-sibling-ria <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-ria.html>`__
   didn’t call ``git update-server-info`` after setting up the remote
   repository and, as a result, the repository couldn’t be fetched until
   something else (e.g., a push) triggered a call to
   ``git update-server-info``.
   (`#5531 <https://github.com/datalad/datalad/issues/5531>`__)

-  The parser for git-config output didn’t properly handle multi-line
   values and got thrown off by unexpected and unrelated lines.
   (`#5509 <https://github.com/datalad/datalad/issues/5509>`__)

-  The 0.14 release introduced regressions in the handling of progress
   bars for git-annex actions, including collapsing progress bars for
   concurrent operations.
   (`#5421 <https://github.com/datalad/datalad/issues/5421>`__)
   (`#5438 <https://github.com/datalad/datalad/issues/5438>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   failed if the user configured Git’s ``diff.ignoreSubmodules`` to a
   non-default value.
   (`#5453 <https://github.com/datalad/datalad/issues/5453>`__)

-  A interprocess lock is now used to prevent a race between checking
   for an SSH socket’s existence and creating it.
   (`#5466 <https://github.com/datalad/datalad/issues/5466>`__)

-  If a Python procedure script is executable,
   `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   invokes it directly rather than passing it to ``sys.executable``. The
   non-executable Python procedures that ship with DataLad now include
   shebangs so that invoking them has a chance of working on file
   systems that present all files as executable.
   (`#5436 <https://github.com/datalad/datalad/issues/5436>`__)

-  DataLad’s wrapper around ``argparse`` failed if an underscore was
   used in a positional argument.
   (`#5525 <https://github.com/datalad/datalad/issues/5525>`__)

.. _enhancements-and-new-features-2:

Enhancements and new features
-----------------------------

-  DataLad’s method for mapping environment variables to configuration
   options (e.g., ``DATALAD_FOO_X__Y`` to ``datalad.foo.x-y``) doesn’t
   work if the subsection name (“FOO”) has an underscore. This
   limitation can be sidestepped with the new
   ``DATALAD_CONFIG_OVERRIDES_JSON`` environment variable, which can be
   set to a JSON record of configuration values.
   (`#5505 <https://github.com/datalad/datalad/issues/5505>`__)

0.14.0 (February 02, 2021) – .
==============================

Major refactoring and deprecations
----------------------------------

-  Git versions below v2.19.1 are no longer supported.
   (`#4650 <https://github.com/datalad/datalad/issues/4650>`__)

-  The minimum git-annex version is still 7.20190503, but, if you’re on
   Windows (or use adjusted branches in general), please upgrade to at
   least 8.20200330 but ideally 8.20210127 to get subdataset-related
   fixes. (`#4292 <https://github.com/datalad/datalad/issues/4292>`__)
   (`#5290 <https://github.com/datalad/datalad/issues/5290>`__)

-  The minimum supported version of Python is now 3.6.
   (`#4879 <https://github.com/datalad/datalad/issues/4879>`__)

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   is now deprecated in favor of
   `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__.
   It will be removed in the 0.15.0 release at the earliest.

-  A new command runner was added in v0.13. Functionality related to the
   old runner has now been removed: ``Runner``, ``GitRunner``, and
   ``run_gitcommand_on_file_list_chunks`` from the ``datalad.cmd``
   module along with the ``datalad.tests.protocolremote``,
   ``datalad.cmd.protocol``, and ``datalad.cmd.protocol.prefix``
   configuration options.
   (`#5229 <https://github.com/datalad/datalad/issues/5229>`__)

-  The ``--no-storage-sibling`` switch of ``create-sibling-ria`` is
   deprecated in favor of ``--storage-sibling=off`` and will be removed
   in a later release.
   (`#5090 <https://github.com/datalad/datalad/issues/5090>`__)

-  The ``get_git_dir`` static method of ``GitRepo`` is deprecated and
   will be removed in a later release. Use the ``dot_git`` attribute of
   an instance instead.
   (`#4597 <https://github.com/datalad/datalad/issues/4597>`__)

-  The ``ProcessAnnexProgressIndicators`` helper from
   ``datalad.support.annexrepo`` has been removed.
   (`#5259 <https://github.com/datalad/datalad/issues/5259>`__)

-  The ``save`` argument of
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__,
   a noop since v0.6.0, has been dropped.
   (`#5278 <https://github.com/datalad/datalad/issues/5278>`__)

-  The ``get_URLS`` method of ``AnnexCustomRemote`` is deprecated and
   will be removed in a later release.
   (`#4955 <https://github.com/datalad/datalad/issues/4955>`__)

-  ``ConfigManager.get`` now returns a single value rather than a tuple
   when there are multiple values for the same key, as very few callers
   correctly accounted for the possibility of a tuple return value.
   Callers can restore the old behavior by passing ``get_all=True``.
   (`#4924 <https://github.com/datalad/datalad/issues/4924>`__)

-  In 0.12.0, all of the ``assure_*`` functions in ``datalad.utils``
   were renamed as ``ensure_*``, keeping the old names around as
   compatibility aliases. The ``assure_*`` variants are now marked as
   deprecated and will be removed in a later release.
   (`#4908 <https://github.com/datalad/datalad/issues/4908>`__)

-  The ``datalad.inteface.run`` module, which was deprecated in 0.12.0
   and kept as a compatibility shim for ``datalad.core.local.run``, has
   been removed.
   (`#4583 <https://github.com/datalad/datalad/issues/4583>`__)

-  The ``saver`` argument of ``datalad.core.local.run.run_command``,
   marked as obsolete in 0.12.0, has been removed.
   (`#4583 <https://github.com/datalad/datalad/issues/4583>`__)

-  The ``dataset_only`` argument of the ``ConfigManager`` class was
   deprecated in 0.12 and has now been removed.
   (`#4828 <https://github.com/datalad/datalad/issues/4828>`__)

-  The ``linux_distribution_name``, ``linux_distribution_release``, and
   ``on_debian_wheezy`` attributes in ``datalad.utils`` are no longer
   set at import time and will be removed in a later release. Use
   ``datalad.utils.get_linux_distribution`` instead.
   (`#4696 <https://github.com/datalad/datalad/issues/4696>`__)

-  ``datalad.distribution.clone``, which was marked as obsolete in v0.12
   in favor of ``datalad.core.distributed.clone``, has been removed.
   (`#4904 <https://github.com/datalad/datalad/issues/4904>`__)

-  ``datalad.support.annexrepo.N_AUTO_JOBS``, announced as deprecated in
   v0.12.6, has been removed.
   (`#4904 <https://github.com/datalad/datalad/issues/4904>`__)

-  The ``compat`` parameter of ``GitRepo.get_submodules``, added in
   v0.12 as a temporary compatibility layer, has been removed.
   (`#4904 <https://github.com/datalad/datalad/issues/4904>`__)

-  The long-deprecated (and non-functional) ``url`` parameter of
   ``GitRepo.__init__`` has been removed.
   (`#5342 <https://github.com/datalad/datalad/issues/5342>`__)

.. _fixes-5:

Fixes
-----

-  Cloning onto a system that enters adjusted branches by default (as
   Windows does) did not properly record the clone URL.
   (`#5128 <https://github.com/datalad/datalad/issues/5128>`__)

-  The RIA-specific handling after calling
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   was correctly triggered by ``ria+http`` URLs but not ``ria+https``
   URLs. (`#4977 <https://github.com/datalad/datalad/issues/4977>`__)

-  If the registered commit wasn’t found when cloning a subdataset, the
   failed attempt was left around.
   (`#5391 <https://github.com/datalad/datalad/issues/5391>`__)

-  The remote calls to ``cp`` and ``chmod`` in
   `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   were not portable and failed on macOS.
   (`#5108 <https://github.com/datalad/datalad/issues/5108>`__)

-  A more reliable check is now done to decide if configuration files
   need to be reloaded.
   (`#5276 <https://github.com/datalad/datalad/issues/5276>`__)

-  The internal command runner’s handling of the event loop has been
   improved to play nicer with outside applications and scripts that use
   asyncio. (`#5350 <https://github.com/datalad/datalad/issues/5350>`__)
   (`#5367 <https://github.com/datalad/datalad/issues/5367>`__)

.. _enhancements-and-new-features-3:

Enhancements and new features
-----------------------------

-  The subdataset handling for adjusted branches, which is particularly
   important on Windows where git-annex enters an adjusted branch by
   default, has been improved. A core piece of the new approach is
   registering the commit of the primary branch, not its checked out
   adjusted branch, in the superdataset. Note: This means that
   ``git status`` will always consider a subdataset on an adjusted
   branch as dirty while ``datalad status`` will look more closely and
   see if the tip of the primary branch matches the registered commit.
   (`#5241 <https://github.com/datalad/datalad/issues/5241>`__)

-  The performance of the
   `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   command has been improved, with substantial speedups for recursive
   processing of many subdatasets.
   (`#4868 <https://github.com/datalad/datalad/issues/4868>`__)
   (`#5076 <https://github.com/datalad/datalad/issues/5076>`__)

-  Adding new subdatasets via
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   has been sped up.
   (`#4793 <https://github.com/datalad/datalad/issues/4793>`__)

-  `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__,
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__,
   and
   `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__
   gained support for parallel operations that can be enabled via the
   ``--jobs`` command-line option or the new
   ``datalad.runtime.max-jobs`` configuration option.
   (`#5022 <https://github.com/datalad/datalad/issues/5022>`__)

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__

   -  learned how to read data from standard input.
      (`#4669 <https://github.com/datalad/datalad/issues/4669>`__)
   -  now supports tab-separated input.
      (`#4845 <https://github.com/datalad/datalad/issues/4845>`__)
   -  now lets Python callers pass in a list of records rather than a
      file name.
      (`#5285 <https://github.com/datalad/datalad/issues/5285>`__)
   -  gained a ``--drop-after`` switch that signals to drop a file’s
      content after downloading and adding it to the annex.
      (`#5081 <https://github.com/datalad/datalad/issues/5081>`__)
   -  is now able to construct a tree of files from known checksums
      without downloading content via its new ``--key`` option.
      (`#5184 <https://github.com/datalad/datalad/issues/5184>`__)
   -  records the URL file in the commit message as provided by the
      caller rather than using the resolved absolute path.
      (`#5091 <https://github.com/datalad/datalad/issues/5091>`__)
   -  is now speedier.
      (`#4867 <https://github.com/datalad/datalad/issues/4867>`__)
      (`#5022 <https://github.com/datalad/datalad/issues/5022>`__)

-  `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   learned how to create private repositories (thanks to Nolan Nichols).
   (`#4769 <https://github.com/datalad/datalad/issues/4769>`__)

-  `create-sibling-ria <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-ria.html>`__
   gained a ``--storage-sibling`` option. When
   ``--storage-sibling=only`` is specified, the storage sibling is
   created without an accompanying Git sibling. This enables using hosts
   without Git installed for storage.
   (`#5090 <https://github.com/datalad/datalad/issues/5090>`__)

-  The download machinery (and thus the ``datalad`` special remote)
   gained support for a new scheme, ``shub://``, which follows the same
   format used by ``singularity run`` and friends. In contrast to the
   short-lived URLs obtained by querying Singularity Hub directly,
   ``shub://`` URLs are suitable for registering with git-annex.
   (`#4816 <https://github.com/datalad/datalad/issues/4816>`__)

-  A provider is now included for https://registry-1.docker.io URLs.
   This is useful for storing an image’s blobs in a dataset and
   registering the URLs with git-annex.
   (`#5129 <https://github.com/datalad/datalad/issues/5129>`__)

-  The ``add-readme`` command now links to the `DataLad
   handbook <http://handbook.datalad.org>`__ rather than
   http://docs.datalad.org.
   (`#4991 <https://github.com/datalad/datalad/issues/4991>`__)

-  New option ``datalad.locations.extra-procedures`` specifies an
   additional location that should be searched for procedures.
   (`#5156 <https://github.com/datalad/datalad/issues/5156>`__)

-  The class for handling configuration values, ``ConfigManager``, now
   takes a lock before writes to allow for multiple processes to modify
   the configuration of a dataset.
   (`#4829 <https://github.com/datalad/datalad/issues/4829>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   now records the original, unresolved URL for a subdataset under
   ``submodule.<name>.datalad-url`` in the parent’s .gitmodules,
   enabling later
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   calls to use the original URL. This is particularly useful for
   ``ria+`` URLs.
   (`#5346 <https://github.com/datalad/datalad/issues/5346>`__)

-  Installing a subdataset now uses custom handling rather than calling
   ``git submodule update --init``. This avoids some locking issues when
   running
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   in parallel and enables more accurate source URLs to be recorded.
   (`#4853 <https://github.com/datalad/datalad/issues/4853>`__)

-  ``GitRepo.get_content_info``, a helper that gets triggered by many
   commands, got faster by tweaking its ``git ls-files`` call.
   (`#5067 <https://github.com/datalad/datalad/issues/5067>`__)

-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__
   now includes credentials-related information (e.g. active backends)
   in the its output.
   (`#4982 <https://github.com/datalad/datalad/issues/4982>`__)

-  The ``call_git*`` methods of ``GitRepo`` now have a ``read_only``
   parameter. Callers can set this to ``True`` to promise that the
   provided command does not write to the repository, bypassing the cost
   of some checks and locking.
   (`#5070 <https://github.com/datalad/datalad/issues/5070>`__)

-  New ``call_annex*`` methods in the ``AnnexRepo`` class provide an
   interface for running git-annex commands similar to that of the
   ``GitRepo.call_git*`` methods.
   (`#5163 <https://github.com/datalad/datalad/issues/5163>`__)

-  It’s now possible to register a custom metadata indexer that is
   discovered by
   `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   and used to generate an index.
   (`#4963 <https://github.com/datalad/datalad/issues/4963>`__)

-  The ``ConfigManager`` methods ``get``, ``getbool``, ``getfloat``, and
   ``getint`` now return a single value (with same precedence as
   ``git config --get``) when there are multiple values for the same key
   (in the non-committed git configuration, if the key is present there,
   or in the dataset configuration). For ``get``, the old behavior can
   be restored by specifying ``get_all=True``.
   (`#4924 <https://github.com/datalad/datalad/issues/4924>`__)

-  Command-line scripts are now defined via the ``entry_points``
   argument of ``setuptools.setup`` instead of the ``scripts`` argument.
   (`#4695 <https://github.com/datalad/datalad/issues/4695>`__)

-  Interactive use of ``--help`` on the command-line now invokes a pager
   on more systems and installation setups.
   (`#5344 <https://github.com/datalad/datalad/issues/5344>`__)

-  The ``datalad`` special remote now tries to eliminate some
   unnecessary interactions with git-annex by being smarter about how it
   queries for URLs associated with a key.
   (`#4955 <https://github.com/datalad/datalad/issues/4955>`__)

-  The ``GitRepo`` class now does a better job of handling bare
   repositories, a step towards bare repositories support in DataLad.
   (`#4911 <https://github.com/datalad/datalad/issues/4911>`__)

-  More internal work to move the code base over to the new command
   runner. (`#4699 <https://github.com/datalad/datalad/issues/4699>`__)
   (`#4855 <https://github.com/datalad/datalad/issues/4855>`__)
   (`#4900 <https://github.com/datalad/datalad/issues/4900>`__)
   (`#4996 <https://github.com/datalad/datalad/issues/4996>`__)
   (`#5002 <https://github.com/datalad/datalad/issues/5002>`__)
   (`#5141 <https://github.com/datalad/datalad/issues/5141>`__)
   (`#5142 <https://github.com/datalad/datalad/issues/5142>`__)
   (`#5229 <https://github.com/datalad/datalad/issues/5229>`__)

0.13.7 (January 04, 2021) – .
=============================

.. _fixes-6:

Fixes
-----

-  Cloning from a RIA store on the local file system initialized annex
   in the Git sibling of the RIA source, which is problematic because
   all annex-related functionality should go through the storage
   sibling.
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   now sets ``remote.origin.annex-ignore`` to ``true`` after cloning
   from RIA stores to prevent this.
   (`#5255 <https://github.com/datalad/datalad/issues/5255>`__)

-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   invoked ``cp`` in a way that was not compatible with macOS.
   (`#5269 <https://github.com/datalad/datalad/issues/5269>`__)

-  Due to a bug in older Git versions (before 2.25), calling
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   with a file under .git/ (e.g., ``datalad status .git/config``)
   incorrectly reported the file as untracked. A workaround has been
   added. (`#5258 <https://github.com/datalad/datalad/issues/5258>`__)

-  Update tests for compatibility with latest git-annex.
   (`#5254 <https://github.com/datalad/datalad/issues/5254>`__)

.. _enhancements-and-new-features-4:

Enhancements and new features
-----------------------------

-  `copy-file <http://datalad.readthedocs.io/en/latest/generated/man/datalad-copy-file.html>`__
   now aborts if .git/ is in the target directory, adding to its
   existing .git/ safety checks.
   (`#5258 <https://github.com/datalad/datalad/issues/5258>`__)

0.13.6 (December 14, 2020) – .
==============================

.. _fixes-7:

Fixes
-----

-  An assortment of fixes for Windows compatibility.
   (`#5113 <https://github.com/datalad/datalad/issues/5113>`__)
   (`#5119 <https://github.com/datalad/datalad/issues/5119>`__)
   (`#5125 <https://github.com/datalad/datalad/issues/5125>`__)
   (`#5127 <https://github.com/datalad/datalad/issues/5127>`__)
   (`#5136 <https://github.com/datalad/datalad/issues/5136>`__)
   (`#5201 <https://github.com/datalad/datalad/issues/5201>`__)
   (`#5200 <https://github.com/datalad/datalad/issues/5200>`__)
   (`#5214 <https://github.com/datalad/datalad/issues/5214>`__)

-  Adding a subdataset on a system that defaults to using an adjusted
   branch (i.e. doesn’t support symlinks) didn’t properly set up the
   submodule URL if the source dataset was not in an adjusted state.
   (`#5127 <https://github.com/datalad/datalad/issues/5127>`__)

-  `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   failed to push to a remote that did not have an ``annex-uuid`` value
   in the local ``.git/config``.
   (`#5148 <https://github.com/datalad/datalad/issues/5148>`__)

-  The default renderer has been improved to avoid a spurious leading
   space, which led to the displayed path being incorrect in some cases.
   (`#5121 <https://github.com/datalad/datalad/issues/5121>`__)

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   showed an uninformative error message when asked to configure an
   unknown remote.
   (`#5146 <https://github.com/datalad/datalad/issues/5146>`__)

-  `drop <http://datalad.readthedocs.io/en/latest/generated/man/datalad-drop.html>`__
   confusingly relayed a suggestion from ``git annex drop`` to use
   ``--force``, an option that does not exist in ``datalad drop``.
   (`#5194 <https://github.com/datalad/datalad/issues/5194>`__)

-  `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   no longer offers user/password authentication because it is no longer
   supported by GitHub.
   (`#5218 <https://github.com/datalad/datalad/issues/5218>`__)

-  The internal command runner’s handling of the event loop has been
   tweaked to hopefully fix issues with runnning DataLad from IPython.
   (`#5106 <https://github.com/datalad/datalad/issues/5106>`__)

-  SSH cleanup wasn’t reliably triggered by the ORA special remote on
   failure, leading to a stall with a particular version of git-annex,
   8.20201103. (This is also resolved on git-annex’s end as of
   8.20201127.)
   (`#5151 <https://github.com/datalad/datalad/issues/5151>`__)

.. _enhancements-and-new-features-5:

Enhancements and new features
-----------------------------

-  The credential helper no longer asks the user to repeat tokens or AWS
   keys. (`#5219 <https://github.com/datalad/datalad/issues/5219>`__)

-  The new option ``datalad.locations.sockets`` controls where Datalad
   stores SSH sockets, allowing users to more easily work around file
   system and path length restrictions.
   (`#5238 <https://github.com/datalad/datalad/issues/5238>`__)

0.13.5 (October 30, 2020) – .
=============================

.. _fixes-8:

Fixes
-----

-  SSH connection handling has been reworked to fix cloning on Windows.
   A new configuration option, ``datalad.ssh.multiplex-connections``,
   defaults to false on Windows.
   (`#5042 <https://github.com/datalad/datalad/issues/5042>`__)

-  The ORA special remote and post-clone RIA configuration now provide
   authentication via DataLad’s credential mechanism and better handling
   of HTTP status codes.
   (`#5025 <https://github.com/datalad/datalad/issues/5025>`__)
   (`#5026 <https://github.com/datalad/datalad/issues/5026>`__)

-  By default, if a git executable is present in the same location as
   git-annex, DataLad modifies ``PATH`` when running git and git-annex
   so that the bundled git is used. This logic has been tightened to
   avoid unnecessarily adjusting the path, reducing the cases where the
   adjustment interferes with the local environment, such as special
   remotes in a virtual environment being masked by the system-wide
   variants.
   (`#5035 <https://github.com/datalad/datalad/issues/5035>`__)

-  git-annex is now consistently invoked as “git annex” rather than
   “git-annex” to work around failures on Windows.
   (`#5001 <https://github.com/datalad/datalad/issues/5001>`__)

-  `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   called ``git annex sync ...`` on plain git repositories.
   (`#5051 <https://github.com/datalad/datalad/issues/5051>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   in genernal doesn’t support registering multiple levels of untracked
   subdatasets, but it can now properly register nested subdatasets when
   all of the subdataset paths are passed explicitly (e.g.,
   ``datalad save -d. sub-a sub-a/sub-b``).
   (`#5049 <https://github.com/datalad/datalad/issues/5049>`__)

-  When called with ``--sidecar`` and ``--explicit``,
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   didn’t save the sidecar.
   (`#5017 <https://github.com/datalad/datalad/issues/5017>`__)

-  A couple of spots didn’t properly quote format fields when combining
   substrings into a format string.
   (`#4957 <https://github.com/datalad/datalad/issues/4957>`__)

-  The default credentials configured for ``indi-s3`` prevented
   anonymous access.
   (`#5045 <https://github.com/datalad/datalad/issues/5045>`__)

.. _enhancements-and-new-features-6:

Enhancements and new features
-----------------------------

-  Messages about suppressed similar results are now rate limited to
   improve performance when there are many similar results coming
   through quickly.
   (`#5060 <https://github.com/datalad/datalad/issues/5060>`__)

-  `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   can now be told to replace an existing sibling by passing
   ``--existing=replace``.
   (`#5008 <https://github.com/datalad/datalad/issues/5008>`__)

-  Progress bars now react to changes in the terminal’s width (requires
   tqdm 2.1 or later).
   (`#5057 <https://github.com/datalad/datalad/issues/5057>`__)

0.13.4 (October 6, 2020) – .
============================

.. _fixes-9:

Fixes
-----

-  Ephemeral clones mishandled bare repositories.
   (`#4899 <https://github.com/datalad/datalad/issues/4899>`__)

-  The post-clone logic for configuring RIA stores didn’t consider
   ``https://`` URLs.
   (`#4977 <https://github.com/datalad/datalad/issues/4977>`__)

-  DataLad custom remotes didn’t escape newlines in messages sent to
   git-annex.
   (`#4926 <https://github.com/datalad/datalad/issues/4926>`__)

-  The datalad-archives special remote incorrectly treated file names as
   percent-encoded.
   (`#4953 <https://github.com/datalad/datalad/issues/4953>`__)

-  The result handler didn’t properly escape “%” when constructing its
   message template.
   (`#4953 <https://github.com/datalad/datalad/issues/4953>`__)

-  In v0.13.0, the tailored rendering for specific subtypes of external
   command failures (e.g., “out of space” or “remote not available”) was
   unintentionally switched to the default rendering.
   (`#4966 <https://github.com/datalad/datalad/issues/4966>`__)

-  Various fixes and updates for the NDA authenticator.
   (`#4824 <https://github.com/datalad/datalad/issues/4824>`__)

-  The helper for getting a versioned S3 URL did not support anonymous
   access or buckets with “.” in their name.
   (`#4985 <https://github.com/datalad/datalad/issues/4985>`__)

-  Several issues with the handling of S3 credentials and token
   expiration have been addressed.
   (`#4927 <https://github.com/datalad/datalad/issues/4927>`__)
   (`#4931 <https://github.com/datalad/datalad/issues/4931>`__)
   (`#4952 <https://github.com/datalad/datalad/issues/4952>`__)

.. _enhancements-and-new-features-7:

Enhancements and new features
-----------------------------

-  A warning is now given if the detected Git is below v2.13.0 to let
   users that run into problems know that their Git version is likely
   the culprit.
   (`#4866 <https://github.com/datalad/datalad/issues/4866>`__)

-  A fix to
   `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   in v0.13.2 introduced a regression that surfaces when
   ``push.default`` is configured to “matching” and prevents the
   git-annex branch from being pushed. Note that, as part of the fix,
   the current branch is now always pushed even when it wouldn’t be
   based on the configured refspec or ``push.default`` value.
   (`#4896 <https://github.com/datalad/datalad/issues/4896>`__)

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__

   -  now allows spelling the empty string value of ``--since=`` as
      ``^`` for consistency with
      `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__.
      (`#4683 <https://github.com/datalad/datalad/issues/4683>`__)
   -  compares a revision given to ``--since=`` with ``HEAD`` rather
      than the working tree to speed up the operation.
      (`#4448 <https://github.com/datalad/datalad/issues/4448>`__)

-  `rerun <https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html>`__

   -  emits more INFO-level log messages.
      (`#4764 <https://github.com/datalad/datalad/issues/4764>`__)
   -  provides better handling of adjusted branches and aborts with a
      clear error for cases that are not supported.
      (`#5328 <https://github.com/datalad/datalad/issues/5328>`__)

-  The archives are handled with p7zip, if available, since DataLad
   v0.12.0. This implementation now supports .tgz and .tbz2 archives.
   (`#4877 <https://github.com/datalad/datalad/issues/4877>`__)

0.13.3 (August 28, 2020) – .
============================

.. _fixes-10:

Fixes
-----

-  Work around a Python bug that led to our asyncio-based command runner
   intermittently failing to capture the output of commands that exit
   very quickly.
   (`#4835 <https://github.com/datalad/datalad/issues/4835>`__)

-  `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   displayed an overestimate of the transfer size when multiple files
   pointed to the same key.
   (`#4821 <https://github.com/datalad/datalad/issues/4821>`__)

-  When
   `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   calls ``git annex addurl``, it catches and reports any failures
   rather than crashing. A change in v0.12.0 broke this handling in a
   particular case.
   (`#4817 <https://github.com/datalad/datalad/issues/4817>`__)

.. _enhancements-and-new-features-8:

Enhancements and new features
-----------------------------

-  The wrapper functions returned by decorators are now given more
   meaningful names to hopefully make tracebacks easier to digest.
   (`#4834 <https://github.com/datalad/datalad/issues/4834>`__)

0.13.2 (August 10, 2020) – .
============================

Deprecations
------------

-  The ``allow_quick`` parameter of ``AnnexRepo.file_has_content`` and
   ``AnnexRepo.is_under_annex`` is now ignored and will be removed in a
   later release. This parameter was only relevant for git-annex
   versions before 7.20190912.
   (`#4736 <https://github.com/datalad/datalad/issues/4736>`__)

.. _fixes-11:

Fixes
-----

-  Updates for compatibility with recent git and git-annex releases.
   (`#4746 <https://github.com/datalad/datalad/issues/4746>`__)
   (`#4760 <https://github.com/datalad/datalad/issues/4760>`__)
   (`#4684 <https://github.com/datalad/datalad/issues/4684>`__)

-  `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   didn’t sync the git-annex branch when ``--data=nothing`` was
   specified.
   (`#4786 <https://github.com/datalad/datalad/issues/4786>`__)

-  The ``datalad.clone.reckless`` configuration wasn’t stored in
   non-annex datasets, preventing the values from being inherited by
   annex subdatasets.
   (`#4749 <https://github.com/datalad/datalad/issues/4749>`__)

-  Running the post-update hook installed by ``create-sibling --ui``
   could overwrite web log files from previous runs in the unlikely
   event that the hook was executed multiple times in the same second.
   (`#4745 <https://github.com/datalad/datalad/issues/4745>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   inspected git’s standard error in a way that could cause an attribute
   error. (`#4775 <https://github.com/datalad/datalad/issues/4775>`__)

-  When cloning a repository whose ``HEAD`` points to a branch without
   commits,
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   tries to find a more useful branch to check out. It unwisely
   considered adjusted branches.
   (`#4792 <https://github.com/datalad/datalad/issues/4792>`__)

-  Since v0.12.0, ``SSHManager.close`` hasn’t closed connections when
   the ``ctrl_path`` argument was explicitly given.
   (`#4757 <https://github.com/datalad/datalad/issues/4757>`__)

-  When working in a dataset in which ``git annex init`` had not yet
   been called, the ``file_has_content`` and ``is_under_annex`` methods
   of ``AnnexRepo`` incorrectly took the “allow quick” code path on file
   systems that did not support it
   (`#4736 <https://github.com/datalad/datalad/issues/4736>`__)

Enhancements
------------

-  `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__
   now assigns version 4 (random) UUIDs instead of version 1 UUIDs that
   encode the time and hardware address.
   (`#4790 <https://github.com/datalad/datalad/issues/4790>`__)

-  The documentation for
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__
   now does a better job of describing the interaction between
   ``--dataset`` and ``PATH``.
   (`#4763 <https://github.com/datalad/datalad/issues/4763>`__)

-  The ``format_commit`` and ``get_hexsha`` methods of ``GitRepo`` have
   been sped up.
   (`#4807 <https://github.com/datalad/datalad/issues/4807>`__)
   (`#4806 <https://github.com/datalad/datalad/issues/4806>`__)

-  A better error message is now shown when the ``^`` or ``^.``
   shortcuts for ``--dataset`` do not resolve to a dataset.
   (`#4759 <https://github.com/datalad/datalad/issues/4759>`__)

-  A more helpful error message is now shown if a caller tries to
   download an ``ftp://`` link but does not have ``request_ftp``
   installed.
   (`#4788 <https://github.com/datalad/datalad/issues/4788>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   now tries harder to get up-to-date availability information after
   auto-enabling ``type=git`` special remotes.
   (`#2897 <https://github.com/datalad/datalad/issues/2897>`__)

0.13.1 (July 17, 2020) – .
==========================

.. _fixes-12:

Fixes
-----

-  Cloning a subdataset should inherit the parent’s
   ``datalad.clone.reckless`` value, but that did not happen when
   cloning via ``datalad get`` rather than ``datalad install`` or
   ``datalad clone``.
   (`#4657 <https://github.com/datalad/datalad/issues/4657>`__)

-  The default result renderer crashed when the result did not have a
   ``path`` key.
   (`#4666 <https://github.com/datalad/datalad/issues/4666>`__)
   (`#4673 <https://github.com/datalad/datalad/issues/4673>`__)

-  ``datalad push`` didn’t show information about ``git push`` errors
   when the output was not in the format that it expected.
   (`#4674 <https://github.com/datalad/datalad/issues/4674>`__)

-  ``datalad push`` silently accepted an empty string for ``--since``
   even though it is an invalid value.
   (`#4682 <https://github.com/datalad/datalad/issues/4682>`__)

-  Our JavaScript testing setup on Travis grew stale and has now been
   updated. (Thanks to Xiao Gui.)
   (`#4687 <https://github.com/datalad/datalad/issues/4687>`__)

-  The new class for running Git commands (added in v0.13.0) ignored any
   changes to the process environment that occurred after instantiation.
   (`#4703 <https://github.com/datalad/datalad/issues/4703>`__)

.. _enhancements-and-new-features-9:

Enhancements and new features
-----------------------------

-  ``datalad push`` now avoids unnecessary ``git push`` dry runs and
   pushes all refspecs with a single ``git push`` call rather than
   invoking ``git push`` for each one.
   (`#4692 <https://github.com/datalad/datalad/issues/4692>`__)
   (`#4675 <https://github.com/datalad/datalad/issues/4675>`__)

-  The readability of SSH error messages has been improved.
   (`#4729 <https://github.com/datalad/datalad/issues/4729>`__)

-  ``datalad.support.annexrepo`` avoids calling
   ``datalad.utils.get_linux_distribution`` at import time and caches
   the result once it is called because, as of Python 3.8, the function
   uses ``distro`` underneath, adding noticeable overhead.
   (`#4696 <https://github.com/datalad/datalad/issues/4696>`__)

   Third-party code should be updated to use ``get_linux_distribution``
   directly in the unlikely event that the code relied on the
   import-time call to ``get_linux_distribution`` setting the
   ``linux_distribution_name``, ``linux_distribution_release``, or
   ``on_debian_wheezy`` attributes in \`datalad.utils.

0.13.0 (June 23, 2020) – .
==========================

A handful of new commands, including ``copy-file``, ``push``, and
``create-sibling-ria``, along with various fixes and enhancements

.. _major-refactoring-and-deprecations-1:

Major refactoring and deprecations
----------------------------------

-  The ``no_annex`` parameter of
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__,
   which is exposed in the Python API but not the command line, is
   deprecated and will be removed in a later release. Use the new
   ``annex`` argument instead, flipping the value. Command-line callers
   that use ``--no-annex`` are unaffected.
   (`#4321 <https://github.com/datalad/datalad/issues/4321>`__)

-  ``datalad add``, which was deprecated in 0.12.0, has been removed.
   (`#4158 <https://github.com/datalad/datalad/issues/4158>`__)
   (`#4319 <https://github.com/datalad/datalad/issues/4319>`__)

-  The following ``GitRepo`` and ``AnnexRepo`` methods have been
   removed: ``get_changed_files``, ``get_missing_files``, and
   ``get_deleted_files``.
   (`#4169 <https://github.com/datalad/datalad/issues/4169>`__)
   (`#4158 <https://github.com/datalad/datalad/issues/4158>`__)

-  The ``get_branch_commits`` method of ``GitRepo`` and ``AnnexRepo``
   has been renamed to ``get_branch_commits_``.
   (`#3834 <https://github.com/datalad/datalad/issues/3834>`__)

-  The custom ``commit`` method of ``AnnexRepo`` has been removed, and
   ``AnnexRepo.commit`` now resolves to the parent method,
   ``GitRepo.commit``.
   (`#4168 <https://github.com/datalad/datalad/issues/4168>`__)

-  GitPython’s ``git.repo.base.Repo`` class is no longer available via
   the ``.repo`` attribute of ``GitRepo`` and ``AnnexRepo``.
   (`#4172 <https://github.com/datalad/datalad/issues/4172>`__)

-  ``AnnexRepo.get_corresponding_branch`` now returns ``None`` rather
   than the current branch name when a managed branch is not checked
   out. (`#4274 <https://github.com/datalad/datalad/issues/4274>`__)

-  The special UUID for git-annex web remotes is now available as
   ``datalad.consts.WEB_SPECIAL_REMOTE_UUID``. It remains accessible as
   ``AnnexRepo.WEB_UUID`` for compatibility, but new code should use
   ``consts.WEB_SPECIAL_REMOTE_UUID``
   (`#4460 <https://github.com/datalad/datalad/issues/4460>`__).

.. _fixes-13:

Fixes
-----

-  Widespread improvements in functionality and test coverage on Windows
   and crippled file systems in general.
   (`#4057 <https://github.com/datalad/datalad/issues/4057>`__)
   (`#4245 <https://github.com/datalad/datalad/issues/4245>`__)
   (`#4268 <https://github.com/datalad/datalad/issues/4268>`__)
   (`#4276 <https://github.com/datalad/datalad/issues/4276>`__)
   (`#4291 <https://github.com/datalad/datalad/issues/4291>`__)
   (`#4296 <https://github.com/datalad/datalad/issues/4296>`__)
   (`#4301 <https://github.com/datalad/datalad/issues/4301>`__)
   (`#4303 <https://github.com/datalad/datalad/issues/4303>`__)
   (`#4304 <https://github.com/datalad/datalad/issues/4304>`__)
   (`#4305 <https://github.com/datalad/datalad/issues/4305>`__)
   (`#4306 <https://github.com/datalad/datalad/issues/4306>`__)

-  ``AnnexRepo.get_size_from_key`` incorrectly handled file chunks.
   (`#4081 <https://github.com/datalad/datalad/issues/4081>`__)

-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   would too readily clobber existing paths when called with
   ``--existing=replace``. It now gets confirmation from the user before
   doing so if running interactively and unconditionally aborts when
   running non-interactively.
   (`#4147 <https://github.com/datalad/datalad/issues/4147>`__)

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   (`#4159 <https://github.com/datalad/datalad/issues/4159>`__)

   -  queried the incorrect branch configuration when updating non-annex
      repositories.
   -  didn’t account for the fact that the local repository can be
      configured as the upstream “remote” for a branch.

-  When the caller included ``--bare`` as a ``git init`` option,
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__
   crashed creating the bare repository, which is currently unsupported,
   rather than aborting with an informative error message.
   (`#4065 <https://github.com/datalad/datalad/issues/4065>`__)

-  The logic for automatically propagating the ‘origin’ remote when
   cloning a local source could unintentionally trigger a fetch of a
   non-local remote.
   (`#4196 <https://github.com/datalad/datalad/issues/4196>`__)

-  All remaining ``get_submodules()`` call sites that relied on the
   temporary compatibility layer added in v0.12.0 have been updated.
   (`#4348 <https://github.com/datalad/datalad/issues/4348>`__)

-  The custom result summary renderer for
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__,
   which was visible with ``--output-format=tailored``, displayed
   incorrect and confusing information in some cases. The custom
   renderer has been removed entirely.
   (`#4471 <https://github.com/datalad/datalad/issues/4471>`__)

-  The documentation for the Python interface of a command listed an
   incorrect default when the command overrode the value of command
   parameters such as ``result_renderer``.
   (`#4480 <https://github.com/datalad/datalad/issues/4480>`__)

.. _enhancements-and-new-features-10:

Enhancements and new features
-----------------------------

-  The default result renderer learned to elide a chain of results after
   seeing ten consecutive results that it considers similar, which
   improves the display of actions that have many results (e.g., saving
   hundreds of files).
   (`#4337 <https://github.com/datalad/datalad/issues/4337>`__)

-  The default result renderer, in addition to “tailored” result
   renderer, now triggers the custom summary renderer, if any.
   (`#4338 <https://github.com/datalad/datalad/issues/4338>`__)

-  The new command
   `create-sibling-ria <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-ria.html>`__
   provides support for creating a sibling in a `RIA
   store <http://handbook.datalad.org/en/latest/usecases/datastorage_for_institutions.html>`__.
   (`#4124 <https://github.com/datalad/datalad/issues/4124>`__)

-  DataLad ships with a new special remote, git-annex-remote-ora, for
   interacting with `RIA
   stores <http://handbook.datalad.org/en/latest/usecases/datastorage_for_institutions.html>`__
   and a new command
   `export-archive-ora <http://datalad.readthedocs.io/en/latest/generated/man/datalad-export-archive-ora.html>`__
   for exporting an archive from a local annex object store.
   (`#4260 <https://github.com/datalad/datalad/issues/4260>`__)
   (`#4203 <https://github.com/datalad/datalad/issues/4203>`__)

-  The new command
   `push <http://datalad.readthedocs.io/en/latest/generated/man/datalad-push.html>`__
   provides an alternative interface to
   `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   for pushing a dataset hierarchy to a sibling.
   (`#4206 <https://github.com/datalad/datalad/issues/4206>`__)
   (`#4581 <https://github.com/datalad/datalad/issues/4581>`__)
   (`#4617 <https://github.com/datalad/datalad/issues/4617>`__)
   (`#4620 <https://github.com/datalad/datalad/issues/4620>`__)

-  The new command
   `copy-file <http://datalad.readthedocs.io/en/latest/generated/man/datalad-copy-file.html>`__
   copies files and associated availability information from one dataset
   to another.
   (`#4430 <https://github.com/datalad/datalad/issues/4430>`__)

-  The command examples have been expanded and improved.
   (`#4091 <https://github.com/datalad/datalad/issues/4091>`__)
   (`#4314 <https://github.com/datalad/datalad/issues/4314>`__)
   (`#4464 <https://github.com/datalad/datalad/issues/4464>`__)

-  The tooling for linking to the `DataLad
   Handbook <http://handbook.datalad.org>`__ from DataLad’s
   documentation has been improved.
   (`#4046 <https://github.com/datalad/datalad/issues/4046>`__)

-  The ``--reckless`` parameter of
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   and
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   learned two new modes:

   -  “ephemeral”, where the .git/annex/ of the cloned repository is
      symlinked to the local source repository’s.
      (`#4099 <https://github.com/datalad/datalad/issues/4099>`__)
   -  “shared-{group|all|…}” that can be used to set up datasets for
      collaborative write access.
      (`#4324 <https://github.com/datalad/datalad/issues/4324>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__

   -  learned to handle dataset aliases in RIA stores when given a URL
      of the form ``ria+<protocol>://<storelocation>#~<aliasname>``.
      (`#4459 <https://github.com/datalad/datalad/issues/4459>`__)
   -  now checks ``datalad.get.subdataset-source-candidate-NAME`` to see
      if ``NAME`` starts with three digits, which is taken as a “cost”.
      Sources with lower costs will be tried first.
      (`#4619 <https://github.com/datalad/datalad/issues/4619>`__)

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   (`#4167 <https://github.com/datalad/datalad/issues/4167>`__)

   -  learned to disallow non-fast-forward updates when ``ff-only`` is
      given to the ``--merge`` option.
   -  gained a ``--follow`` option that controls how ``--merge``
      behaves, adding support for merging in the revision that is
      registered in the parent dataset rather than merging in the
      configured branch from the sibling.
   -  now provides a result record for merge events.

-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   now supports local paths as targets in addition to SSH URLs.
   (`#4187 <https://github.com/datalad/datalad/issues/4187>`__)

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   now

   -  shows a warning if the caller requests to delete a sibling that
      does not exist.
      (`#4257 <https://github.com/datalad/datalad/issues/4257>`__)
   -  phrases its warning about non-annex repositories in a less
      alarming way.
      (`#4323 <https://github.com/datalad/datalad/issues/4323>`__)

-  The rendering of command errors has been improved.
   (`#4157 <https://github.com/datalad/datalad/issues/4157>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   now

   -  displays a message to signal that the working tree is clean,
      making it more obvious that no results being rendered corresponds
      to a clean state.
      (`#4106 <https://github.com/datalad/datalad/issues/4106>`__)
   -  provides a stronger warning against using ``--to-git``.
      (`#4290 <https://github.com/datalad/datalad/issues/4290>`__)

-  `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   and
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   learned about scenarios where they could avoid unnecessary and
   expensive work.
   (`#4526 <https://github.com/datalad/datalad/issues/4526>`__)
   (`#4544 <https://github.com/datalad/datalad/issues/4544>`__)
   (`#4549 <https://github.com/datalad/datalad/issues/4549>`__)

-  Calling
   `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   without ``--recursive`` but with a path constraint within a
   subdataset (“/”) now traverses into the subdataset, as “/” would,
   restricting its report to “/”.
   (`#4235 <https://github.com/datalad/datalad/issues/4235>`__)

-  New option ``datalad.annex.retry`` controls how many times git-annex
   will retry on a failed transfer. It defaults to 3 and can be set to 0
   to restore the previous behavior.
   (`#4382 <https://github.com/datalad/datalad/issues/4382>`__)

-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__
   now warns when the specified dataset does not exist.
   (`#4331 <https://github.com/datalad/datalad/issues/4331>`__)

-  The ``repr`` and ``str`` output of the dataset and repo classes got a
   facelift.
   (`#4420 <https://github.com/datalad/datalad/issues/4420>`__)
   (`#4435 <https://github.com/datalad/datalad/issues/4435>`__)
   (`#4439 <https://github.com/datalad/datalad/issues/4439>`__)

-  The DataLad Singularity container now comes with p7zip-full.

-  DataLad emits a log message when the current working directory is
   resolved to a different location due to a symlink. This is now logged
   at the DEBUG rather than WARNING level, as it typically does not
   indicate a problem.
   (`#4426 <https://github.com/datalad/datalad/issues/4426>`__)

-  DataLad now lets the caller know that ``git annex init`` is scanning
   for unlocked files, as this operation can be slow in some
   repositories.
   (`#4316 <https://github.com/datalad/datalad/issues/4316>`__)

-  The ``log_progress`` helper learned how to set the starting point to
   a non-zero value and how to update the total of an existing progress
   bar, two features needed for planned improvements to how some
   commands display their progress.
   (`#4438 <https://github.com/datalad/datalad/issues/4438>`__)

-  The ``ExternalVersions`` object, which is used to check versions of
   Python modules and external tools (e.g., git-annex), gained an
   ``add`` method that enables DataLad extensions and other third-party
   code to include other programs of interest.
   (`#4441 <https://github.com/datalad/datalad/issues/4441>`__)

-  All of the remaining spots that use GitPython have been rewritten
   without it. Most notably, this includes rewrites of the ``clone``,
   ``fetch``, and ``push`` methods of ``GitRepo``.
   (`#4080 <https://github.com/datalad/datalad/issues/4080>`__)
   (`#4087 <https://github.com/datalad/datalad/issues/4087>`__)
   (`#4170 <https://github.com/datalad/datalad/issues/4170>`__)
   (`#4171 <https://github.com/datalad/datalad/issues/4171>`__)
   (`#4175 <https://github.com/datalad/datalad/issues/4175>`__)
   (`#4172 <https://github.com/datalad/datalad/issues/4172>`__)

-  When ``GitRepo.commit`` splits its operation across multiple calls to
   avoid exceeding the maximum command line length, it now amends to
   initial commit rather than creating multiple commits.
   (`#4156 <https://github.com/datalad/datalad/issues/4156>`__)

-  ``GitRepo`` gained a ``get_corresponding_branch`` method (which
   always returns None), allowing a caller to invoke the method without
   needing to check if the underlying repo class is ``GitRepo`` or
   ``AnnexRepo``.
   (`#4274 <https://github.com/datalad/datalad/issues/4274>`__)

-  A new helper function ``datalad.core.local.repo.repo_from_path``
   returns a repo class for a specified path.
   (`#4273 <https://github.com/datalad/datalad/issues/4273>`__)

-  New ``AnnexRepo`` method ``localsync`` performs a ``git annex sync``
   that disables external interaction and is particularly useful for
   propagating changes on an adjusted branch back to the main branch.
   (`#4243 <https://github.com/datalad/datalad/issues/4243>`__)

0.12.7 (May 22, 2020) – .
=========================

.. _fixes-14:

Fixes
-----

-  Requesting tailored output (``--output=tailored``) from a command
   with a custom result summary renderer produced repeated output.
   (`#4463 <https://github.com/datalad/datalad/issues/4463>`__)

-  A longstanding regression in argcomplete-based command-line
   completion for Bash has been fixed. You can enable completion by
   configuring a Bash startup file to run
   ``eval "$(register-python-argcomplete datalad)"`` or source DataLad’s
   ``tools/cmdline-completion``. The latter should work for Zsh as well.
   (`#4477 <https://github.com/datalad/datalad/issues/4477>`__)

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   didn’t prevent ``git-fetch`` from recursing into submodules, leading
   to a failure when the registered submodule was not present locally
   and the submodule did not have a remote named ‘origin’.
   (`#4560 <https://github.com/datalad/datalad/issues/4560>`__)

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__
   botched path handling when the file name format started with “./” and
   the call was made from a subdirectory of the dataset.
   (`#4504 <https://github.com/datalad/datalad/issues/4504>`__)

-  Double dash options in manpages were unintentionally escaped.
   (`#4332 <https://github.com/datalad/datalad/issues/4332>`__)

-  The check for HTTP authentication failures crashed in situations
   where content came in as bytes rather than unicode.
   (`#4543 <https://github.com/datalad/datalad/issues/4543>`__)

-  A check in ``AnnexRepo.whereis`` could lead to a type error.
   (`#4552 <https://github.com/datalad/datalad/issues/4552>`__)

-  When installing a dataset to obtain a subdataset,
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   confusingly displayed a message that described the containing dataset
   as “underneath” the subdataset.
   (`#4456 <https://github.com/datalad/datalad/issues/4456>`__)

-  A couple of Makefile rules didn’t properly quote paths.
   (`#4481 <https://github.com/datalad/datalad/issues/4481>`__)

-  With DueCredit support enabled (``DUECREDIT_ENABLE=1``), the query
   for metadata information could flood the output with warnings if
   datasets didn’t have aggregated metadata. The warnings are now
   silenced, with the overall failure of a
   `metadata <http://datalad.readthedocs.io/en/latest/generated/man/datalad-metadata.html>`__
   call logged at the debug level.
   (`#4568 <https://github.com/datalad/datalad/issues/4568>`__)

.. _enhancements-and-new-features-11:

Enhancements and new features
-----------------------------

-  The resource identifier helper learned to recognize URLs with
   embedded Git transport information, such as
   gcrypt::https://example.com.
   (`#4529 <https://github.com/datalad/datalad/issues/4529>`__)

-  When running non-interactively, a more informative error is now
   signaled when the UI backend, which cannot display a question, is
   asked to do so.
   (`#4553 <https://github.com/datalad/datalad/issues/4553>`__)

0.12.6 (April 23, 2020) – .
===========================

.. _major-refactoring-and-deprecations-2:

Major refactoring and deprecations
----------------------------------

-  The value of ``datalad.support.annexrep.N_AUTO_JOBS`` is no longer
   considered. The variable will be removed in a later release.
   (`#4409 <https://github.com/datalad/datalad/issues/4409>`__)

.. _fixes-15:

Fixes
-----

-  Staring with v0.12.0, ``datalad save`` recorded the current branch of
   a parent dataset as the ``branch`` value in the .gitmodules entry for
   a subdataset. This behavior is problematic for a few reasons and has
   been reverted.
   (`#4375 <https://github.com/datalad/datalad/issues/4375>`__)

-  The default for the ``--jobs`` option, “auto”, instructed DataLad to
   pass a value to git-annex’s ``--jobs`` equal to
   ``min(8, max(3, <number of CPUs>))``, which could lead to issues due
   to the large number of child processes spawned and file descriptors
   opened. To avoid this behavior, ``--jobs=auto`` now results in
   git-annex being called with ``--jobs=1`` by default. Configure the
   new option ``datalad.runtime.max-annex-jobs`` to control the maximum
   value that will be considered when ``--jobs='auto'``.
   (`#4409 <https://github.com/datalad/datalad/issues/4409>`__)

-  Various commands have been adjusted to better handle the case where a
   remote’s HEAD ref points to an unborn branch.
   (`#4370 <https://github.com/datalad/datalad/issues/4370>`__)

-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__

   -  learned to use the query as a regular expression that restricts
      the keys that are shown for ``--show-keys short``.
      (`#4354 <https://github.com/datalad/datalad/issues/4354>`__)
   -  gives a more helpful message when query is an invalid regular
      expression.
      (`#4398 <https://github.com/datalad/datalad/issues/4398>`__)

-  The code for parsing Git configuration did not follow Git’s behavior
   of accepting a key with no value as shorthand for key=true.
   (`#4421 <https://github.com/datalad/datalad/issues/4421>`__)

-  ``AnnexRepo.info`` needed a compatibility update for a change in how
   git-annex reports file names.
   (`#4431 <https://github.com/datalad/datalad/issues/4431>`__)

-  `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   did not gracefully handle a token that did not have the necessary
   permissions.
   (`#4400 <https://github.com/datalad/datalad/issues/4400>`__)

.. _enhancements-and-new-features-12:

Enhancements and new features
-----------------------------

-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   learned to use the query as a regular expression that restricts the
   keys that are shown for ``--show-keys short``.
   (`#4354 <https://github.com/datalad/datalad/issues/4354>`__)

-  ``datalad <subcommand>`` learned to point to the
   `datalad-container <https://github.com/datalad/datalad-container>`__
   extension when a subcommand from that extension is given but the
   extension is not installed.
   (`#4400 <https://github.com/datalad/datalad/issues/4400>`__)
   (`#4174 <https://github.com/datalad/datalad/issues/4174>`__)

0.12.5 (Apr 02, 2020) – a small step for datalad …
==================================================

Fix some bugs and make the world an even better place.

.. _fixes-16:

Fixes
-----

-  Our ``log_progress`` helper mishandled the initial display and step
   of the progress bar.
   (`#4326 <https://github.com/datalad/datalad/issues/4326>`__)

-  ``AnnexRepo.get_content_annexinfo`` is designed to accept
   ``init=None``, but passing that led to an error.
   (`#4330 <https://github.com/datalad/datalad/issues/4330>`__)

-  Update a regular expression to handle an output change in Git
   v2.26.0. (`#4328 <https://github.com/datalad/datalad/issues/4328>`__)

-  We now set ``LC_MESSAGES`` to ‘C’ while running git to avoid failures
   when parsing output that is marked for translation.
   (`#4342 <https://github.com/datalad/datalad/issues/4342>`__)

-  The helper for decoding JSON streams loaded the last line of input
   without decoding it if the line didn’t end with a new line, a
   regression introduced in the 0.12.0 release.
   (`#4361 <https://github.com/datalad/datalad/issues/4361>`__)

-  The clone command failed to git-annex-init a fresh clone whenever it
   considered to add the origin of the origin as a remote.
   (`#4367 <https://github.com/datalad/datalad/issues/4367>`__)

0.12.4 (Mar 19, 2020) – Windows?!
=================================

The main purpose of this release is to have one on PyPi that has no
associated wheel to enable a working installation on Windows
(`#4315 <https://github.com/datalad/datalad/issues/4315>`__).

.. _fixes-17:

Fixes
-----

-  The description of the ``log.outputs`` config switch did not keep up
   with code changes and incorrectly stated that the output would be
   logged at the DEBUG level; logging actually happens at a lower level.
   (`#4317 <https://github.com/datalad/datalad/issues/4317>`__)

0.12.3 (March 16, 2020) – .
===========================

Updates for compatibility with the latest git-annex, along with a few
miscellaneous fixes

.. _major-refactoring-and-deprecations-3:

Major refactoring and deprecations
----------------------------------

-  All spots that raised a ``NoDatasetArgumentFound`` exception now
   raise a ``NoDatasetFound`` exception to better reflect the situation:
   it is the *dataset* rather than the *argument* that is not found. For
   compatibility, the latter inherits from the former, but new code
   should prefer the latter.
   (`#4285 <https://github.com/datalad/datalad/issues/4285>`__)

.. _fixes-18:

Fixes
-----

-  Updates for compatibility with git-annex version 8.20200226.
   (`#4214 <https://github.com/datalad/datalad/issues/4214>`__)

-  ``datalad export-to-figshare`` failed to export if the generated
   title was fewer than three characters. It now queries the caller for
   the title and guards against titles that are too short.
   (`#4140 <https://github.com/datalad/datalad/issues/4140>`__)

-  Authentication was requested multiple times when git-annex launched
   parallel downloads from the ``datalad`` special remote.
   (`#4308 <https://github.com/datalad/datalad/issues/4308>`__)

-  At verbose logging levels, DataLad requests that git-annex display
   debugging information too. Work around a bug in git-annex that
   prevented that from happening.
   (`#4212 <https://github.com/datalad/datalad/issues/4212>`__)

-  The internal command runner looked in the wrong place for some
   configuration variables, including ``datalad.log.outputs``, resulting
   in the default value always being used.
   (`#4194 <https://github.com/datalad/datalad/issues/4194>`__)

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   failed when trying to publish to a git-lfs special remote for the
   first time.
   (`#4200 <https://github.com/datalad/datalad/issues/4200>`__)

-  ``AnnexRepo.set_remote_url`` is supposed to establish shared SSH
   connections but failed to do so.
   (`#4262 <https://github.com/datalad/datalad/issues/4262>`__)

.. _enhancements-and-new-features-13:

Enhancements and new features
-----------------------------

-  The message provided when a command cannot determine what dataset to
   operate on has been improved.
   (`#4285 <https://github.com/datalad/datalad/issues/4285>`__)

-  The “aws-s3” authentication type now allows specifying the host
   through “aws-s3_host”, which was needed to work around an
   authorization error due to a longstanding upstream bug.
   (`#4239 <https://github.com/datalad/datalad/issues/4239>`__)

-  The xmp metadata extractor now recognizes “.wav” files.

0.12.2 (Jan 28, 2020) – Smoothen the ride
=========================================

Mostly a bugfix release with various robustifications, but also makes
the first step towards versioned dataset installation requests.

.. _major-refactoring-and-deprecations-4:

Major refactoring and deprecations
----------------------------------

-  The minimum required version for GitPython is now 2.1.12.
   (`#4070 <https://github.com/datalad/datalad/issues/4070>`__)

.. _fixes-19:

Fixes
-----

-  The class for handling configuration values, ``ConfigManager``,
   inappropriately considered the current working directory’s dataset,
   if any, for both reading and writing when instantiated with
   ``dataset=None``. This misbehavior is fairly inaccessible through
   typical use of DataLad. It affects ``datalad.cfg``, the top-level
   configuration instance that should not consider repository-specific
   values. It also affects Python users that call ``Dataset`` with a
   path that does not yet exist and persists until that dataset is
   created. (`#4078 <https://github.com/datalad/datalad/issues/4078>`__)

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   saved the dataset when called with ``--merge``, which is unnecessary
   and risks committing unrelated changes.
   (`#3996 <https://github.com/datalad/datalad/issues/3996>`__)

-  Confusing and irrelevant information about Python defaults have been
   dropped from the command-line help.
   (`#4002 <https://github.com/datalad/datalad/issues/4002>`__)

-  The logic for automatically propagating the ‘origin’ remote when
   cloning a local source didn’t properly account for relative paths.
   (`#4045 <https://github.com/datalad/datalad/issues/4045>`__)

-  Various fixes to file name handling and quoting on Windows.
   (`#4049 <https://github.com/datalad/datalad/issues/4049>`__)
   (`#4050 <https://github.com/datalad/datalad/issues/4050>`__)

-  When cloning failed, error lines were not bubbled up to the user in
   some scenarios.
   (`#4060 <https://github.com/datalad/datalad/issues/4060>`__)

.. _enhancements-and-new-features-14:

Enhancements and new features
-----------------------------

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   (and thus
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__)

   -  now propagates the ``reckless`` mode from the superdataset when
      cloning a dataset into it.
      (`#4037 <https://github.com/datalad/datalad/issues/4037>`__)
   -  gained support for ``ria+<protocol>://`` URLs that point to
      `RIA <http://handbook.datalad.org/en/latest/usecases/datastorage_for_institutions.html>`__
      stores.
      (`#4022 <https://github.com/datalad/datalad/issues/4022>`__)
   -  learned to read “@version” from ``ria+`` URLs and install that
      version of a dataset
      (`#4036 <https://github.com/datalad/datalad/issues/4036>`__) and
      to apply URL rewrites configured through Git’s ``url.*.insteadOf``
      mechanism
      (`#4064 <https://github.com/datalad/datalad/issues/4064>`__).
   -  now copies ``datalad.get.subdataset-source-candidate-<name>``
      options configured within the superdataset into the subdataset.
      This is particularly useful for RIA data stores.
      (`#4073 <https://github.com/datalad/datalad/issues/4073>`__)

-  Archives are now (optionally) handled with 7-Zip instead of
   ``patool``. 7-Zip will be used by default, but ``patool`` will be
   used on non-Windows systems if the ``datalad.runtime.use-patool``
   option is set or the ``7z`` executable is not found.
   (`#4041 <https://github.com/datalad/datalad/issues/4041>`__)

0.12.1 (Jan 15, 2020) – Small bump after big bang
=================================================

Fix some fallout after major release.

.. _fixes-20:

Fixes
-----

-  Revert incorrect relative path adjustment to URLs in
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__.
   (`#3538 <https://github.com/datalad/datalad/issues/3538>`__)

-  Various small fixes to internal helpers and test to run on Windows
   (`#2566 <https://github.com/datalad/datalad/issues/2566>`__)
   (`#2534 <https://github.com/datalad/datalad/issues/2534>`__)

0.12.0 (Jan 11, 2020) – Krakatoa
================================

This release is the result of more than a year of development that
includes fixes for a large number of issues, yielding more robust
behavior across a wider range of use cases, and introduces major changes
in API and behavior. It is the first release for which extensive user
documentation is available in a dedicated `DataLad
Handbook <http://handbook.datalad.org>`__. Python 3 (3.5 and later) is
now the only supported Python flavor.

Major changes 0.12 vs 0.11
--------------------------

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   fully replaces
   `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__
   (which is obsolete now, and will be removed in a future release).

-  A new Git-annex aware
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   command enables detailed inspection of dataset hierarchies. The
   previously available
   `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   command has been adjusted to match
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   in argument semantics and behavior.

-  The ability to configure dataset procedures prior and after the
   execution of particular commands has been replaced by a flexible
   “hook” mechanism that is able to run arbitrary DataLad commands
   whenever command results are detected that match a specification.

-  Support of the Windows platform has been improved substantially.
   While performance and feature coverage on Windows still falls behind
   Unix-like systems, typical data consumer use cases, and standard
   dataset operations, such as
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__
   and
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__,
   are now working. Basic support for data provenance capture via
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   is also functional.

-  Support for Git-annex direct mode repositories has been removed,
   following the end of support in Git-annex itself.

-  The semantics of relative paths in command line arguments have
   changed. Previously, a call
   ``datalad save --dataset /tmp/myds some/relpath`` would have been
   interpreted as saving a file at ``/tmp/myds/some/relpath`` into
   dataset ``/tmp/myds``. This has changed to saving
   ``$PWD/some/relpath`` into dataset ``/tmp/myds``. More generally,
   relative paths are now always treated as relative to the current
   working directory, except for path arguments of
   `Dataset <http://docs.datalad.org/en/latest/generated/datalad.api.Dataset.html>`__
   class instance methods of the Python API. The resulting partial
   duplication of path specifications between path and dataset arguments
   is mitigated by the introduction of two special symbols that can be
   given as dataset argument: ``^`` and ``^.``, which identify the
   topmost superdataset and the closest dataset that contains the
   working directory, respectively.

-  The concept of a “core API” has been introduced. Commands situated in
   the module ``datalad.core`` (such as
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__,
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__,
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__,
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__,
   `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__)
   receive additional scrutiny regarding API and implementation, and are
   meant to provide longer-term stability. Application developers are
   encouraged to preferentially build on these commands.

Major refactoring and deprecations since 0.12.0rc6
--------------------------------------------------

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   has been incorporated into the growing core API. The public
   ``--alternative-source`` parameter has been removed, and a
   ``clone_dataset`` function with multi-source capabilities is provided
   instead. The ``--reckless`` parameter can now take literal mode
   labels instead of just beeing a binary flag, but backwards
   compatibility is maintained.

-  The ``get_file_content`` method of ``GitRepo`` was no longer used
   internally or in any known DataLad extensions and has been removed.
   (`#3812 <https://github.com/datalad/datalad/issues/3812>`__)

-  The function ``get_dataset_root`` has been replaced by
   ``rev_get_dataset_root``. ``rev_get_dataset_root`` remains as a
   compatibility alias and will be removed in a later release.
   (`#3815 <https://github.com/datalad/datalad/issues/3815>`__)

-  The ``add_sibling`` module, marked obsolete in v0.6.0, has been
   removed. (`#3871 <https://github.com/datalad/datalad/issues/3871>`__)

-  ``mock`` is no longer declared as an external dependency because we
   can rely on it being in the standard library now that our minimum
   required Python version is 3.5.
   (`#3860 <https://github.com/datalad/datalad/issues/3860>`__)

-  `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   now requires that directories be indicated with a trailing slash
   rather than interpreting a path as directory when it doesn’t exist.
   This avoids confusion that can result from typos and makes it
   possible to support directory targets that do not exist.
   (`#3854 <https://github.com/datalad/datalad/issues/3854>`__)

-  The ``dataset_only`` argument of the ``ConfigManager`` class is
   deprecated. Use ``source="dataset"`` instead.
   (`#3907 <https://github.com/datalad/datalad/issues/3907>`__)

-  The ``--proc-pre`` and ``--proc-post`` options have been removed, and
   configuration values for ``datalad.COMMAND.proc-pre`` and
   ``datalad.COMMAND.proc-post`` are no longer honored. The new result
   hook mechanism provides an alternative for ``proc-post`` procedures.
   (`#3963 <https://github.com/datalad/datalad/issues/3963>`__)

Fixes since 0.12.0rc6
---------------------

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   crashed when called with a detached HEAD. It now aborts with an
   informative message.
   (`#3804 <https://github.com/datalad/datalad/issues/3804>`__)

-  Since 0.12.0rc6 the call to
   `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   in
   `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   resulted in a spurious warning.
   (`#3877 <https://github.com/datalad/datalad/issues/3877>`__)

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   crashed if it encountered an annex repository that was marked as
   dead. (`#3892 <https://github.com/datalad/datalad/issues/3892>`__)

-  The update of
   `rerun <https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html>`__
   in v0.12.0rc3 for the rewritten
   `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   command didn’t account for a change in the output of ``diff``,
   leading to ``rerun --report`` unintentionally including unchanged
   files in its diff values.
   (`#3873 <https://github.com/datalad/datalad/issues/3873>`__)

-  In 0.12.0rc5
   `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   was updated to follow the new path handling logic, but its calls to
   AnnexRepo weren’t properly adjusted, resulting in incorrect path
   handling when the called from a dataset subdirectory.
   (`#3850 <https://github.com/datalad/datalad/issues/3850>`__)

-  `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   called ``git annex addurl`` in a way that failed to register a URL
   when its header didn’t report the content size.
   (`#3911 <https://github.com/datalad/datalad/issues/3911>`__)

-  With Git v2.24.0, saving new subdatasets failed due to a bug in that
   Git release.
   (`#3904 <https://github.com/datalad/datalad/issues/3904>`__)

-  With DataLad configured to stop on failure (e.g., specifying
   ``--on-failure=stop`` from the command line), a failing result record
   was not rendered.
   (`#3863 <https://github.com/datalad/datalad/issues/3863>`__)

-  Installing a subdataset yielded an “ok” status in cases where the
   repository was not yet in its final state, making it ineffective for
   a caller to operate on the repository in response to the result.
   (`#3906 <https://github.com/datalad/datalad/issues/3906>`__)

-  The internal helper for converting git-annex’s JSON output did not
   relay information from the “error-messages” field.
   (`#3931 <https://github.com/datalad/datalad/issues/3931>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   reported relative paths that were confusingly not relative to the
   current directory in some cases. It now always reports absolute
   paths. (`#3959 <https://github.com/datalad/datalad/issues/3959>`__)

-  `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   inappropriately reported files as deleted in some cases when ``to``
   was a value other than ``None``.
   (`#3999 <https://github.com/datalad/datalad/issues/3999>`__)

-  An assortment of fixes for Windows compatibility.
   (`#3971 <https://github.com/datalad/datalad/issues/3971>`__)
   (`#3974 <https://github.com/datalad/datalad/issues/3974>`__)
   (`#3975 <https://github.com/datalad/datalad/issues/3975>`__)
   (`#3976 <https://github.com/datalad/datalad/issues/3976>`__)
   (`#3979 <https://github.com/datalad/datalad/issues/3979>`__)

-  Subdatasets installed from a source given by relative path will now
   have this relative path used as ‘url’ in their .gitmodules record,
   instead of an absolute path generated by Git.
   (`#3538 <https://github.com/datalad/datalad/issues/3538>`__)

-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   will now correctly interpret ‘~/…’ paths as absolute path
   specifications.
   (`#3958 <https://github.com/datalad/datalad/issues/3958>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   mistakenly reported a directory as a procedure.
   (`#3793 <https://github.com/datalad/datalad/issues/3793>`__)

-  The cleanup for batched git-annex processes has been improved.
   (`#3794 <https://github.com/datalad/datalad/issues/3794>`__)
   (`#3851 <https://github.com/datalad/datalad/issues/3851>`__)

-  The function for adding a version ID to an AWS S3 URL doesn’t support
   URLs with an “s3://” scheme and raises a ``NotImplementedError``
   exception when it encounters one. The function learned to return a
   URL untouched if an “s3://” URL comes in with a version ID.
   (`#3842 <https://github.com/datalad/datalad/issues/3842>`__)

-  A few spots needed to be adjusted for compatibility with git-annex’s
   new ``--sameas``
   `feature <https://git-annex.branchable.com/tips/multiple_remotes_accessing_the_same_data_store/>`__,
   which allows special remotes to share a data store.
   (`#3856 <https://github.com/datalad/datalad/issues/3856>`__)

-  The ``swallow_logs`` utility failed to capture some log messages due
   to an incompatibility with Python 3.7.
   (`#3935 <https://github.com/datalad/datalad/issues/3935>`__)

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__

   -  crashed if ``--inherit`` was passed but the parent dataset did not
      have a remote with a matching name.
      (`#3954 <https://github.com/datalad/datalad/issues/3954>`__)
   -  configured the wrong pushurl and annexurl values in some cases.
      (`#3955 <https://github.com/datalad/datalad/issues/3955>`__)

Enhancements and new features since 0.12.0rc6
---------------------------------------------

-  By default, datasets cloned from local source paths will now get a
   configured remote for any recursively discoverable ‘origin’ sibling
   that is also available from a local path in order to maximize
   automatic file availability across local annexes.
   (`#3926 <https://github.com/datalad/datalad/issues/3926>`__)

-  The new `result hooks
   mechanism <http://handbook.datalad.org/en/latest/basics/101-145-hooks.html>`__
   allows callers to specify, via local Git configuration values,
   DataLad command calls that will be triggered in response to matching
   result records (i.e., what you see when you call a command with
   ``-f json_pp``).
   (`#3903 <https://github.com/datalad/datalad/issues/3903>`__)

-  The command interface classes learned to use a new ``_examples_``
   attribute to render documentation examples for both the Python and
   command-line API.
   (`#3821 <https://github.com/datalad/datalad/issues/3821>`__)

-  Candidate URLs for cloning a submodule can now be generated based on
   configured templates that have access to various properties of the
   submodule, including its dataset ID.
   (`#3828 <https://github.com/datalad/datalad/issues/3828>`__)

-  DataLad’s check that the user’s Git identity is configured has been
   sped up and now considers the appropriate environment variables as
   well. (`#3807 <https://github.com/datalad/datalad/issues/3807>`__)

-  The ``tag`` method of ``GitRepo`` can now tag revisions other than
   ``HEAD`` and accepts a list of arbitrary ``git tag`` options.
   (`#3787 <https://github.com/datalad/datalad/issues/3787>`__)

-  When ``get`` clones a subdataset and the subdataset’s HEAD differs
   from the commit that is registered in the parent, the active branch
   of the subdataset is moved to the registered commit if the registered
   commit is an ancestor of the subdataset’s HEAD commit. This handling
   has been moved to a more central location within ``GitRepo``, and now
   applies to any ``update_submodule(..., init=True)`` call.
   (`#3831 <https://github.com/datalad/datalad/issues/3831>`__)

-  The output of ``datalad -h`` has been reformatted to improve
   readability.
   (`#3862 <https://github.com/datalad/datalad/issues/3862>`__)

-  `unlock <http://datalad.readthedocs.io/en/latest/generated/man/datalad-unlock.html>`__
   has been sped up.
   (`#3880 <https://github.com/datalad/datalad/issues/3880>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   learned to provide and render more information about discovered
   procedures, including whether the procedure is overridden by another
   procedure with the same base name.
   (`#3960 <https://github.com/datalad/datalad/issues/3960>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   now (`#3817 <https://github.com/datalad/datalad/issues/3817>`__)

   -  records the active branch in the superdataset when registering a
      new subdataset.
   -  calls ``git annex sync`` when saving a dataset on an adjusted
      branch so that the changes are brought into the mainline branch.

-  `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   now aborts when its ``dataset`` argument points to a non-existent
   dataset. (`#3940 <https://github.com/datalad/datalad/issues/3940>`__)

-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__
   now

   -  reports the dataset ID if the current working directory is
      visiting a dataset.
      (`#3888 <https://github.com/datalad/datalad/issues/3888>`__)
   -  outputs entries deterministically.
      (`#3927 <https://github.com/datalad/datalad/issues/3927>`__)

-  The ``ConfigManager`` class

   -  learned to exclude ``.datalad/config`` as a source of
      configuration values, restricting the sources to standard Git
      configuration files, when called with ``source="local"``.
      (`#3907 <https://github.com/datalad/datalad/issues/3907>`__)
   -  accepts a value of “override” for its ``where`` argument to allow
      Python callers to more convenient override configuration.
      (`#3970 <https://github.com/datalad/datalad/issues/3970>`__)

-  Commands now accept a ``dataset`` value of “^.” as shorthand for “the
   dataset to which the current directory belongs”.
   (`#3242 <https://github.com/datalad/datalad/issues/3242>`__)

0.12.0rc6 (Oct 19, 2019) – some releases are better than the others
===================================================================

bet we will fix some bugs and make a world even a better place.

.. _major-refactoring-and-deprecations-5:

Major refactoring and deprecations
----------------------------------

-  DataLad no longer supports Python 2. The minimum supported version of
   Python is now 3.5.
   (`#3629 <https://github.com/datalad/datalad/issues/3629>`__)

-  Much of the user-focused content at http://docs.datalad.org has been
   removed in favor of more up to date and complete material available
   in the `DataLad Handbook <http://handbook.datalad.org>`__. Going
   forward, the plan is to restrict http://docs.datalad.org to technical
   documentation geared at developers.
   (`#3678 <https://github.com/datalad/datalad/issues/3678>`__)

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   used to allow the caller to specify which dataset(s) to update as a
   ``PATH`` argument or via the the ``--dataset`` option; now only the
   latter is supported. Path arguments only serve to restrict which
   subdataset are updated when operating recursively.
   (`#3700 <https://github.com/datalad/datalad/issues/3700>`__)

-  Result records from a
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   call no longer have a “state” key.
   (`#3746 <https://github.com/datalad/datalad/issues/3746>`__)

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   and
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   no longer support operating on independent hierarchies of datasets.
   (`#3700 <https://github.com/datalad/datalad/issues/3700>`__)
   (`#3746 <https://github.com/datalad/datalad/issues/3746>`__)

-  The
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   update in 0.12.0rc4 for the new path resolution logic broke the
   handling of inputs and outputs for calls from a subdirectory.
   (`#3747 <https://github.com/datalad/datalad/issues/3747>`__)

-  The ``is_submodule_modified`` method of ``GitRepo`` as well as two
   helper functions in gitrepo.py, ``kwargs_to_options`` and
   ``split_remote_branch``, were no longer used internally or in any
   known DataLad extensions and have been removed.
   (`#3702 <https://github.com/datalad/datalad/issues/3702>`__)
   (`#3704 <https://github.com/datalad/datalad/issues/3704>`__)

-  The ``only_remote`` option of ``GitRepo.is_with_annex`` was not used
   internally or in any known extensions and has been dropped.
   (`#3768 <https://github.com/datalad/datalad/issues/3768>`__)

-  The ``get_tags`` method of ``GitRepo`` used to sort tags by committer
   date. It now sorts them by the tagger date for annotated tags and the
   committer date for lightweight tags.
   (`#3715 <https://github.com/datalad/datalad/issues/3715>`__)

-  The ``rev_resolve_path`` substituted ``resolve_path`` helper.
   (`#3797 <https://github.com/datalad/datalad/issues/3797>`__)

.. _fixes-21:

Fixes
-----

-  Correctly handle relative paths in
   `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__.
   (`#3799 <https://github.com/datalad/datalad/issues/3799>`__)
   (`#3102 <https://github.com/datalad/datalad/issues/3102>`__)

-  Do not errorneously discover directory as a procedure.
   (`#3793 <https://github.com/datalad/datalad/issues/3793>`__)

-  Correctly extract version from manpage to trigger use of manpages for
   ``--help``.
   (`#3798 <https://github.com/datalad/datalad/issues/3798>`__)

-  The ``cfg_yoda`` procedure saved all modifications in the repository
   rather than saving only the files it modified.
   (`#3680 <https://github.com/datalad/datalad/issues/3680>`__)

-  Some spots in the documentation that were supposed appear as two
   hyphens were incorrectly rendered in the HTML output en-dashs.
   (`#3692 <https://github.com/datalad/datalad/issues/3692>`__)

-  `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__,
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__,
   and
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   treated paths as relative to the dataset even when the string form
   was given, violating the new path handling rules.
   (`#3749 <https://github.com/datalad/datalad/issues/3749>`__)
   (`#3777 <https://github.com/datalad/datalad/issues/3777>`__)
   (`#3780 <https://github.com/datalad/datalad/issues/3780>`__)

-  Providing the “^” shortcut to ``--dataset`` didn’t work properly when
   called from a subdirectory of a subdataset.
   (`#3772 <https://github.com/datalad/datalad/issues/3772>`__)

-  We failed to propagate some errors from git-annex when working with
   its JSON output.
   (`#3751 <https://github.com/datalad/datalad/issues/3751>`__)

-  With the Python API, callers are allowed to pass a string or list of
   strings as the ``cfg_proc`` argument to
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__,
   but the string form was mishandled.
   (`#3761 <https://github.com/datalad/datalad/issues/3761>`__)

-  Incorrect command quoting for SSH calls on Windows that rendered
   basic SSH-related functionality (e.g.,
   `sshrun <http://datalad.readthedocs.io/en/latest/generated/man/datalad-sshrun.html>`__)
   on Windows unusable.
   (`#3688 <https://github.com/datalad/datalad/issues/3688>`__)

-  Annex JSON result handling assumed platform-specific paths on Windows
   instead of the POSIX-style that is happening across all platforms.
   (`#3719 <https://github.com/datalad/datalad/issues/3719>`__)

-  ``path_is_under()`` was incapable of comparing Windows paths with
   different drive letters.
   (`#3728 <https://github.com/datalad/datalad/issues/3728>`__)

.. _enhancements-and-new-features-15:

Enhancements and new features
-----------------------------

-  Provide a collection of “public” ``call_git*`` helpers within GitRepo
   and replace use of “private” and less specific
   ``_git_custom_command`` calls.
   (`#3791 <https://github.com/datalad/datalad/issues/3791>`__)

-  `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   gained a ``--report-filetype``. Setting it to “raw” can give a
   performance boost for the price of no longer distinguishing symlinks
   that point to annexed content from other symlinks.
   (`#3701 <https://github.com/datalad/datalad/issues/3701>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   disables file type reporting by
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   to improve performance.
   (`#3712 <https://github.com/datalad/datalad/issues/3712>`__)

-  `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   (`#3743 <https://github.com/datalad/datalad/issues/3743>`__)

   -  now extends its result records with a ``contains`` field that
      lists which ``contains`` arguments matched a given subdataset.
   -  yields an ‘impossible’ result record when a ``contains`` argument
      wasn’t matched to any of the reported subdatasets.

-  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   now shows more readable output when cloning fails.
   (`#3775 <https://github.com/datalad/datalad/issues/3775>`__)

-  ``SSHConnection`` now displays a more informative error message when
   it cannot start the ``ControlMaster`` process.
   (`#3776 <https://github.com/datalad/datalad/issues/3776>`__)

-  If the new configuration option ``datalad.log.result-level`` is set
   to a single level, all result records will be logged at that level.
   If you’ve been bothered by DataLad’s double reporting of failures,
   consider setting this to “debug”.
   (`#3754 <https://github.com/datalad/datalad/issues/3754>`__)

-  Configuration values from ``datalad -c OPTION=VALUE ...`` are now
   validated to provide better errors.
   (`#3695 <https://github.com/datalad/datalad/issues/3695>`__)

-  `rerun <https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html>`__
   learned how to handle history with merges. As was already the case
   when cherry picking non-run commits, re-creating merges may results
   in conflicts, and ``rerun`` does not yet provide an interface to let
   the user handle these.
   (`#2754 <https://github.com/datalad/datalad/issues/2754>`__)

-  The ``fsck`` method of ``AnnexRepo`` has been enhanced to expose more
   features of the underlying ``git fsck`` command.
   (`#3693 <https://github.com/datalad/datalad/issues/3693>`__)

-  ``GitRepo`` now has a ``for_each_ref_`` method that wraps
   ``git for-each-ref``, which is used in various spots that used to
   rely on GitPython functionality.
   (`#3705 <https://github.com/datalad/datalad/issues/3705>`__)

-  Do not pretend to be able to work in optimized (``python -O``) mode,
   crash early with an informative message.
   (`#3803 <https://github.com/datalad/datalad/issues/3803>`__)

0.12.0rc5 (September 04, 2019) – .
==================================

Various fixes and enhancements that bring the 0.12.0 release closer.

.. _major-refactoring-and-deprecations-6:

Major refactoring and deprecations
----------------------------------

-  The two modules below have a new home. The old locations still exist
   as compatibility shims and will be removed in a future release.

   -  ``datalad.distribution.subdatasets`` has been moved to
      ``datalad.local.subdatasets``
      (`#3429 <https://github.com/datalad/datalad/issues/3429>`__)
   -  ``datalad.interface.run`` has been moved to
      ``datalad.core.local.run``
      (`#3444 <https://github.com/datalad/datalad/issues/3444>`__)

-  The ``lock`` method of ``AnnexRepo`` and the ``options`` parameter of
   ``AnnexRepo.unlock`` were unused internally and have been removed.
   (`#3459 <https://github.com/datalad/datalad/issues/3459>`__)

-  The ``get_submodules`` method of ``GitRepo`` has been rewritten
   without GitPython. When the new ``compat`` flag is true (the current
   default), the method returns a value that is compatible with the old
   return value. This backwards-compatible return value and the
   ``compat`` flag will be removed in a future release.
   (`#3508 <https://github.com/datalad/datalad/issues/3508>`__)

-  The logic for resolving relative paths given to a command has changed
   (`#3435 <https://github.com/datalad/datalad/issues/3435>`__). The new
   rule is that relative paths are taken as relative to the dataset only
   if a dataset *instance* is passed by the caller. In all other
   scenarios they’re considered relative to the current directory.

   The main user-visible difference from the command line is that using
   the ``--dataset`` argument does *not* result in relative paths being
   taken as relative to the specified dataset. (The undocumented
   distinction between “rel/path” and “./rel/path” no longer exists.)

   All commands under ``datalad.core`` and ``datalad.local``, as well as
   ``unlock`` and ``addurls``, follow the new logic. The goal is for all
   commands to eventually do so.

.. _fixes-22:

Fixes
-----

-  The function for loading JSON streams wasn’t clever enough to handle
   content that included a Unicode line separator like U2028.
   (`#3524 <https://github.com/datalad/datalad/issues/3524>`__)

-  When
   `unlock <http://datalad.readthedocs.io/en/latest/generated/man/datalad-unlock.html>`__
   was called without an explicit target (i.e., a directory or no paths
   at all), the call failed if any of the files did not have content
   present. (`#3459 <https://github.com/datalad/datalad/issues/3459>`__)

-  ``AnnexRepo.get_content_info`` failed in the rare case of a key
   without size information.
   (`#3534 <https://github.com/datalad/datalad/issues/3534>`__)

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   ignored ``--on-failure`` in its underlying call to
   `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__.
   (`#3470 <https://github.com/datalad/datalad/issues/3470>`__)

-  Calling
   `remove <http://datalad.readthedocs.io/en/latest/generated/man/datalad-remove.html>`__
   with a subdirectory displayed spurious warnings about the
   subdirectory files not existing.
   (`#3586 <https://github.com/datalad/datalad/issues/3586>`__)

-  Our processing of ``git-annex --json`` output mishandled info
   messages from special remotes.
   (`#3546 <https://github.com/datalad/datalad/issues/3546>`__)

-  `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__

   -  didn’t bypass the “existing subdataset” check when called with
      ``--force`` as of 0.12.0rc3
      (`#3552 <https://github.com/datalad/datalad/issues/3552>`__)
   -  failed to register the up-to-date revision of a subdataset when
      ``--cfg-proc`` was used with ``--dataset``
      (`#3591 <https://github.com/datalad/datalad/issues/3591>`__)

-  The base downloader had some error handling that wasn’t compatible
   with Python 3.
   (`#3622 <https://github.com/datalad/datalad/issues/3622>`__)

-  Fixed a number of Unicode py2-compatibility issues.
   (`#3602 <https://github.com/datalad/datalad/issues/3602>`__)

-  ``AnnexRepo.get_content_annexinfo`` did not properly chunk file
   arguments to avoid exceeding the command-line character limit.
   (`#3587 <https://github.com/datalad/datalad/issues/3587>`__)

.. _enhancements-and-new-features-16:

Enhancements and new features
-----------------------------

-  New command ``create-sibling-gitlab`` provides an interface for
   creating a publication target on a GitLab instance.
   (`#3447 <https://github.com/datalad/datalad/issues/3447>`__)

-  `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   (`#3429 <https://github.com/datalad/datalad/issues/3429>`__)

   -  now supports path-constrained queries in the same manner as
      commands like ``save`` and ``status``
   -  gained a ``--contains=PATH`` option that can be used to restrict
      the output to datasets that include a specific path.
   -  now narrows the listed subdatasets to those underneath the current
      directory when called with no arguments

-  `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__
   learned to accept a plain ``--annex`` (no value) as shorthand for
   ``--annex basic``.
   (`#3534 <https://github.com/datalad/datalad/issues/3534>`__)

-  The ``.dirty`` property of ``GitRepo`` and ``AnnexRepo`` has been
   sped up. (`#3460 <https://github.com/datalad/datalad/issues/3460>`__)

-  The ``get_content_info`` method of ``GitRepo``, used by ``status``
   and commands that depend on ``status``, now restricts its git calls
   to a subset of files, if possible, for a performance gain in
   repositories with many files.
   (`#3508 <https://github.com/datalad/datalad/issues/3508>`__)

-  Extensions that do not provide a command, such as those that provide
   only metadata extractors, are now supported.
   (`#3531 <https://github.com/datalad/datalad/issues/3531>`__)

-  When calling git-annex with ``--json``, we log standard error at the
   debug level rather than the warning level if a non-zero exit is
   expected behavior.
   (`#3518 <https://github.com/datalad/datalad/issues/3518>`__)

-  `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__
   no longer refuses to create a new dataset in the odd scenario of an
   empty .git/ directory upstairs.
   (`#3475 <https://github.com/datalad/datalad/issues/3475>`__)

-  As of v2.22.0 Git treats a sub-repository on an unborn branch as a
   repository rather than as a directory. Our documentation and tests
   have been updated appropriately.
   (`#3476 <https://github.com/datalad/datalad/issues/3476>`__)

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__
   learned to accept a ``--cfg-proc`` value and pass it to its
   ``create`` calls.
   (`#3562 <https://github.com/datalad/datalad/issues/3562>`__)

0.12.0rc4 (May 15, 2019) – the revolution is over
=================================================

With the replacement of the ``save`` command implementation with
``rev-save`` the revolution effort is now over, and the set of key
commands for local dataset operations (``create``, ``run``, ``save``,
``status``, ``diff``) is now complete. This new core API is available
from ``datalad.core.local`` (and also via ``datalad.api``, as any other
command).

.. _major-refactoring-and-deprecations-7:

Major refactoring and deprecations
----------------------------------

-  The ``add`` command is now deprecated. It will be removed in a future
   release.

.. _fixes-23:

Fixes
-----

-  Remove hard-coded dependencies on POSIX path conventions in SSH
   support code
   (`#3400 <https://github.com/datalad/datalad/issues/3400>`__)

-  Emit an ``add`` result when adding a new subdataset during
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   (`#3398 <https://github.com/datalad/datalad/issues/3398>`__)

-  SSH file transfer now actually opens a shared connection, if none
   exists yet
   (`#3403 <https://github.com/datalad/datalad/issues/3403>`__)

.. _enhancements-and-new-features-17:

Enhancements and new features
-----------------------------

-  ``SSHConnection`` now offers methods for file upload and dowload
   (``get()``, ``put()``. The previous ``copy()`` method only supported
   upload and was discontinued
   (`#3401 <https://github.com/datalad/datalad/issues/3401>`__)

0.12.0rc3 (May 07, 2019) – the revolution continues
===================================================

Continues API consolidation and replaces the ``create`` and ``diff``
command with more performant implementations.

.. _major-refactoring-and-deprecations-8:

Major refactoring and deprecations
----------------------------------

-  The previous ``diff`` command has been replaced by the diff variant
   from the
   `datalad-revolution <http://github.com/datalad/datalad-revolution>`__
   extension.
   (`#3366 <https://github.com/datalad/datalad/issues/3366>`__)

-  ``rev-create`` has been renamed to ``create``, and the previous
   ``create`` has been removed.
   (`#3383 <https://github.com/datalad/datalad/issues/3383>`__)

-  The procedure ``setup_yoda_dataset`` has been renamed to ``cfg_yoda``
   (`#3353 <https://github.com/datalad/datalad/issues/3353>`__).

-  The ``--nosave`` of ``addurls`` now affects only added content, not
   newly created subdatasets
   (`#3259 <https://github.com/datalad/datalad/issues/3259>`__).

-  ``Dataset.get_subdatasets`` (deprecated since v0.9.0) has been
   removed. (`#3336 <https://github.com/datalad/datalad/issues/3336>`__)

-  The ``.is_dirty`` method of ``GitRepo`` and ``AnnexRepo`` has been
   replaced by ``.status`` or, for a subset of cases, the ``.dirty``
   property.
   (`#3330 <https://github.com/datalad/datalad/issues/3330>`__)

-  ``AnnexRepo.get_status`` has been replaced by ``AnnexRepo.status``.
   (`#3330 <https://github.com/datalad/datalad/issues/3330>`__)

.. _fixes-24:

Fixes
-----

-  `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__

   -  reported on directories that contained only ignored files
      (`#3238 <https://github.com/datalad/datalad/issues/3238>`__)
   -  gave a confusing failure when called from a subdataset with an
      explicitly specified dataset argument and “.” as a path
      (`#3325 <https://github.com/datalad/datalad/issues/3325>`__)
   -  misleadingly claimed that the locally present content size was
      zero when ``--annex basic`` was specified
      (`#3378 <https://github.com/datalad/datalad/issues/3378>`__)

-  An informative error wasn’t given when a download provider was
   invalid. (`#3258 <https://github.com/datalad/datalad/issues/3258>`__)

-  Calling ``rev-save PATH`` saved unspecified untracked subdatasets.
   (`#3288 <https://github.com/datalad/datalad/issues/3288>`__)

-  The available choices for command-line options that take values are
   now displayed more consistently in the help output.
   (`#3326 <https://github.com/datalad/datalad/issues/3326>`__)

-  The new pathlib-based code had various encoding issues on Python 2.
   (`#3332 <https://github.com/datalad/datalad/issues/3332>`__)

.. _enhancements-and-new-features-18:

Enhancements and new features
-----------------------------

-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__
   now includes information about the Python version.
   (`#3255 <https://github.com/datalad/datalad/issues/3255>`__)

-  When operating in an annex repository, checking whether git-annex is
   available is now delayed until a call to git-annex is actually
   needed, allowing systems without git-annex to operate on annex
   repositories in a restricted fashion.
   (`#3274 <https://github.com/datalad/datalad/issues/3274>`__)

-  The ``load_stream`` on helper now supports auto-detection of
   compressed files.
   (`#3289 <https://github.com/datalad/datalad/issues/3289>`__)

-  ``create`` (formerly ``rev-create``)

   -  learned to be speedier by passing a path to ``status``
      (`#3294 <https://github.com/datalad/datalad/issues/3294>`__)
   -  gained a ``--cfg-proc`` (or ``-c``) convenience option for running
      configuration procedures (or more accurately any procedure that
      begins with “cfg\_”) in the newly created dataset
      (`#3353 <https://github.com/datalad/datalad/issues/3353>`__)

-  ``AnnexRepo.set_metadata`` now returns a list while
   ``AnnexRepo.set_metadata_`` returns a generator, a behavior which is
   consistent with the ``add`` and ``add_`` method pair.
   (`#3298 <https://github.com/datalad/datalad/issues/3298>`__)

-  ``AnnexRepo.get_metadata`` now supports batch querying of known annex
   files. Note, however, that callers should carefully validate the
   input paths because the batch call will silently hang if given
   non-annex files.
   (`#3364 <https://github.com/datalad/datalad/issues/3364>`__)

-  `status <http://datalad.readthedocs.io/en/latest/generated/man/datalad-status.html>`__

   -  now reports a “bytesize” field for files tracked by Git
      (`#3299 <https://github.com/datalad/datalad/issues/3299>`__)
   -  gained a new option ``eval_subdataset_state`` that controls how
      the subdataset state is evaluated. Depending on the information
      you need, you can select a less expensive mode to make ``status``
      faster.
      (`#3324 <https://github.com/datalad/datalad/issues/3324>`__)
   -  colors deleted files “red”
      (`#3334 <https://github.com/datalad/datalad/issues/3334>`__)

-  Querying repository content is faster due to batching of
   ``git cat-file`` calls.
   (`#3301 <https://github.com/datalad/datalad/issues/3301>`__)

-  The dataset ID of a subdataset is now recorded in the superdataset.
   (`#3304 <https://github.com/datalad/datalad/issues/3304>`__)

-  ``GitRepo.diffstatus``

   -  now avoids subdataset recursion when the comparison is not with
      the working tree, which substantially improves performance when
      diffing large dataset hierarchies
      (`#3314 <https://github.com/datalad/datalad/issues/3314>`__)
   -  got smarter and faster about labeling a subdataset as “modified”
      (`#3343 <https://github.com/datalad/datalad/issues/3343>`__)

-  ``GitRepo.get_content_info`` now supports disabling the file type
   evaluation, which gives a performance boost in cases where this
   information isn’t needed.
   (`#3362 <https://github.com/datalad/datalad/issues/3362>`__)

-  The XMP metadata extractor now filters based on file name to improve
   its performance.
   (`#3329 <https://github.com/datalad/datalad/issues/3329>`__)

0.12.0rc2 (Mar 18, 2019) – revolution!
======================================

.. _fixes-25:

Fixes
-----

-  ``GitRepo.dirty`` does not report on nested empty directories
   (`#3196 <https://github.com/datalad/datalad/issues/3196>`__).

-  ``GitRepo.save()`` reports results on deleted files.

.. _enhancements-and-new-features-19:

Enhancements and new features
-----------------------------

-  Absorb a new set of core commands from the datalad-revolution
   extension:

   -  ``rev-status``: like ``git status``, but simpler and working with
      dataset hierarchies
   -  ``rev-save``: a 2-in-1 replacement for save and add
   -  ``rev-create``: a ~30% faster create

-  JSON support tools can now read and write compressed files.

0.12.0rc1 (Mar 03, 2019) – to boldly go …
=========================================

.. _major-refactoring-and-deprecations-9:

Major refactoring and deprecations
----------------------------------

-  Discontinued support for git-annex direct-mode (also no longer
   supported upstream).

.. _enhancements-and-new-features-20:

Enhancements and new features
-----------------------------

-  Dataset and Repo object instances are now hashable, and can be
   created based on pathlib Path object instances

-  Imported various additional methods for the Repo classes to query
   information and save changes.

0.11.8 (Oct 11, 2019) – annex-we-are-catching-up
================================================

.. _fixes-26:

Fixes
-----

-  Our internal command runner failed to capture output in some cases.
   (`#3656 <https://github.com/datalad/datalad/issues/3656>`__)
-  Workaround in the tests around python in cPython >= 3.7.5 ‘;’ in the
   filename confusing mimetypes
   (`#3769 <https://github.com/datalad/datalad/issues/3769>`__)
   (`#3770 <https://github.com/datalad/datalad/issues/3770>`__)

.. _enhancements-and-new-features-21:

Enhancements and new features
-----------------------------

-  Prepared for upstream changes in git-annex, including support for the
   latest git-annex

   -  7.20190912 auto-upgrades v5 repositories to v7.
      (`#3648 <https://github.com/datalad/datalad/issues/3648>`__)
      (`#3682 <https://github.com/datalad/datalad/issues/3682>`__)
   -  7.20191009 fixed treatment of (larger/smaller)than in
      .gitattributes
      (`#3765 <https://github.com/datalad/datalad/issues/3765>`__)

-  The ``cfg_text2git`` procedure, as well the ``--text-no-annex``
   option of
   `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__,
   now configure .gitattributes so that empty files are stored in git
   rather than annex.
   (`#3667 <https://github.com/datalad/datalad/issues/3667>`__)

0.11.7 (Sep 06, 2019) – python2-we-still-love-you-but-…
=======================================================

Primarily bugfixes with some optimizations and refactorings.

.. _fixes-27:

Fixes
-----

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__

   -  now provides better handling when the URL file isn’t in the
      expected format.
      (`#3579 <https://github.com/datalad/datalad/issues/3579>`__)
   -  always considered a relative file for the URL file argument as
      relative to the current working directory, which goes against the
      convention used by other commands of taking relative paths as
      relative to the dataset argument.
      (`#3582 <https://github.com/datalad/datalad/issues/3582>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__

   -  hard coded “python” when formatting the command for non-executable
      procedures ending with “.py”. ``sys.executable`` is now used.
      (`#3624 <https://github.com/datalad/datalad/issues/3624>`__)
   -  failed if arguments needed more complicated quoting than simply
      surrounding the value with double quotes. This has been resolved
      for systems that support ``shlex.quote``, but note that on Windows
      values are left unquoted.
      (`#3626 <https://github.com/datalad/datalad/issues/3626>`__)

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   now displays an informative error message if a local path is given to
   ``--url`` but ``--name`` isn’t specified.
   (`#3555 <https://github.com/datalad/datalad/issues/3555>`__)

-  `sshrun <http://datalad.readthedocs.io/en/latest/generated/man/datalad-sshrun.html>`__,
   the command DataLad uses for ``GIT_SSH_COMMAND``, didn’t support all
   the parameters that Git expects it to.
   (`#3616 <https://github.com/datalad/datalad/issues/3616>`__)

-  Fixed a number of Unicode py2-compatibility issues.
   (`#3597 <https://github.com/datalad/datalad/issues/3597>`__)

-  `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   now will create leading directories of the output path if they do not
   exist (`#3646 <https://github.com/datalad/datalad/issues/3646>`__)

.. _enhancements-and-new-features-22:

Enhancements and new features
-----------------------------

-  The
   `annotate-paths <http://docs.datalad.org/en/latest/generated/man/datalad-annotate-paths.html>`__
   helper now caches subdatasets it has seen to avoid unnecessary calls.
   (`#3570 <https://github.com/datalad/datalad/issues/3570>`__)

-  A repeated configuration query has been dropped from the handling of
   ``--proc-pre`` and ``--proc-post``.
   (`#3576 <https://github.com/datalad/datalad/issues/3576>`__)

-  Calls to ``git annex find`` now use ``--in=.`` instead of the alias
   ``--in=here`` to take advantage of an optimization that git-annex (as
   of the current release, 7.20190730) applies only to the former.
   (`#3574 <https://github.com/datalad/datalad/issues/3574>`__)

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__
   now suggests close matches when the URL or file format contains an
   unknown field.
   (`#3594 <https://github.com/datalad/datalad/issues/3594>`__)

-  Shared logic used in the setup.py files of Datalad and its extensions
   has been moved to modules in the \_datalad_build_support/ directory.
   (`#3600 <https://github.com/datalad/datalad/issues/3600>`__)

-  Get ready for upcoming git-annex dropping support for direct mode
   (`#3631 <https://github.com/datalad/datalad/issues/3631>`__)

0.11.6 (Jul 30, 2019) – am I the last of 0.11.x?
================================================

Primarily bug fixes to achieve more robust performance

.. _fixes-28:

Fixes
-----

-  Our tests needed various adjustments to keep up with upstream changes
   in Travis and Git.
   (`#3479 <https://github.com/datalad/datalad/issues/3479>`__)
   (`#3492 <https://github.com/datalad/datalad/issues/3492>`__)
   (`#3493 <https://github.com/datalad/datalad/issues/3493>`__)

-  ``AnnexRepo.is_special_annex_remote`` was too selective in what it
   considered to be a special remote.
   (`#3499 <https://github.com/datalad/datalad/issues/3499>`__)

-  We now provide information about unexpected output when git-annex is
   called with ``--json``.
   (`#3516 <https://github.com/datalad/datalad/issues/3516>`__)

-  Exception logging in the ``__del__`` method of ``GitRepo`` and
   ``AnnexRepo`` no longer fails if the names it needs are no longer
   bound. (`#3527 <https://github.com/datalad/datalad/issues/3527>`__)

-  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__
   botched the construction of subdataset paths that were more than two
   levels deep and failed to create datasets in a reliable,
   breadth-first order.
   (`#3561 <https://github.com/datalad/datalad/issues/3561>`__)

-  Cloning a ``type=git`` special remote showed a spurious warning about
   the remote not being enabled.
   (`#3547 <https://github.com/datalad/datalad/issues/3547>`__)

.. _enhancements-and-new-features-23:

Enhancements and new features
-----------------------------

-  For calls to git and git-annex, we disable automatic garbage
   collection due to past issues with GitPython’s state becoming stale,
   but doing so results in a larger .git/objects/ directory that isn’t
   cleaned up until garbage collection is triggered outside of DataLad.
   Tests with the latest GitPython didn’t reveal any state issues, so
   we’ve re-enabled automatic garbage collection.
   (`#3458 <https://github.com/datalad/datalad/issues/3458>`__)

-  `rerun <https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html>`__
   learned an ``--explicit`` flag, which it relays to its calls to
   [run][[]]. This makes it possible to call ``rerun`` in a dirty
   working tree
   (`#3498 <https://github.com/datalad/datalad/issues/3498>`__).

-  The
   `metadata <http://datalad.readthedocs.io/en/latest/generated/man/datalad-metadata.html>`__
   command aborts earlier if a metadata extractor is unavailable.
   (`#3525 <https://github.com/datalad/datalad/issues/3525>`__)

0.11.5 (May 23, 2019) – stability is not overrated
==================================================

Should be faster and less buggy, with a few enhancements.

.. _fixes-29:

Fixes
-----

-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   (`#3318 <https://github.com/datalad/datalad/issues/3318>`__)

   -  Siblings are no longer configured with a post-update hook unless a
      web interface is requested with ``--ui``.
   -  ``git submodule update --init`` is no longer called from the
      post-update hook.
   -  If ``--inherit`` is given for a dataset without a superdataset, a
      warning is now given instead of raising an error.

-  The internal command runner failed on Python 2 when its ``env``
   argument had unicode values.
   (`#3332 <https://github.com/datalad/datalad/issues/3332>`__)
-  The safeguard that prevents creating a dataset in a subdirectory that
   already contains tracked files for another repository failed on Git
   versions before 2.14. For older Git versions, we now warn the caller
   that the safeguard is not active.
   (`#3347 <https://github.com/datalad/datalad/issues/3347>`__)
-  A regression introduced in v0.11.1 prevented
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   from committing changes under a subdirectory when the subdirectory
   was specified as a path argument.
   (`#3106 <https://github.com/datalad/datalad/issues/3106>`__)
-  A workaround introduced in v0.11.1 made it possible for
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   to do a partial commit with an annex file that has gone below the
   ``annex.largefiles`` threshold. The logic of this workaround was
   faulty, leading to files being displayed as typechanged in the index
   following the commit.
   (`#3365 <https://github.com/datalad/datalad/issues/3365>`__)
-  The resolve_path() helper confused paths that had a semicolon for SSH
   RIs. (`#3425 <https://github.com/datalad/datalad/issues/3425>`__)
-  The detection of SSH RIs has been improved.
   (`#3425 <https://github.com/datalad/datalad/issues/3425>`__)

.. _enhancements-and-new-features-24:

Enhancements and new features
-----------------------------

-  The internal command runner was too aggressive in its decision to
   sleep. (`#3322 <https://github.com/datalad/datalad/issues/3322>`__)
-  The “INFO” label in log messages now retains the default text color
   for the terminal rather than using white, which only worked well for
   terminals with dark backgrounds.
   (`#3334 <https://github.com/datalad/datalad/issues/3334>`__)
-  A short flag ``-R`` is now available for the ``--recursion-limit``
   flag, a flag shared by several subcommands.
   (`#3340 <https://github.com/datalad/datalad/issues/3340>`__)
-  The authentication logic for
   `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   has been revamped and now supports 2FA.
   (`#3180 <https://github.com/datalad/datalad/issues/3180>`__)
-  New configuration option ``datalad.ui.progressbar`` can be used to
   configure the default backend for progress reporting (“none”, for
   example, results in no progress bars being shown).
   (`#3396 <https://github.com/datalad/datalad/issues/3396>`__)
-  A new progress backend, available by setting datalad.ui.progressbar
   to “log”, replaces progress bars with a log message upon completion
   of an action.
   (`#3396 <https://github.com/datalad/datalad/issues/3396>`__)
-  DataLad learned to consult the `NO_COLOR <https://no-color.org/>`__
   environment variable and the new ``datalad.ui.color`` configuration
   option when deciding to color output. The default value, “auto”,
   retains the current behavior of coloring output if attached to a TTY
   (`#3407 <https://github.com/datalad/datalad/issues/3407>`__).
-  `clean <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clean.html>`__
   now removes annex transfer directories, which is useful for cleaning
   up failed downloads.
   (`#3374 <https://github.com/datalad/datalad/issues/3374>`__)
-  `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   no longer refuses to clone into a local path that looks like a URL,
   making its behavior consistent with ``git clone``.
   (`#3425 <https://github.com/datalad/datalad/issues/3425>`__)
-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__

   -  Learned to fall back to the ``dist`` package if ``platform.dist``,
      which has been removed in the yet-to-be-release Python 3.8, does
      not exist.
      (`#3439 <https://github.com/datalad/datalad/issues/3439>`__)
   -  Gained a ``--section`` option for limiting the output to specific
      sections and a ``--decor`` option, which currently knows how to
      format the output as GitHub’s ``<details>`` section.
      (`#3440 <https://github.com/datalad/datalad/issues/3440>`__)

0.11.4 (Mar 18, 2019) – get-ready
=================================

Largely a bug fix release with a few enhancements

Important
---------

-  0.11.x series will be the last one with support for direct mode of
   `git-annex <http://git-annex.branchable.com/>`__ which is used on
   crippled (no symlinks and no locking) filesystems. v7 repositories
   should be used instead.

.. _fixes-30:

Fixes
-----

-  Extraction of .gz files is broken without p7zip installed. We now
   abort with an informative error in this situation.
   (`#3176 <https://github.com/datalad/datalad/issues/3176>`__)

-  Committing failed in some cases because we didn’t ensure that the
   path passed to ``git read-tree --index-output=...`` resided on the
   same filesystem as the repository.
   (`#3181 <https://github.com/datalad/datalad/issues/3181>`__)

-  Some pointless warnings during metadata aggregation have been
   eliminated.
   (`#3186 <https://github.com/datalad/datalad/issues/3186>`__)

-  With Python 3 the LORIS token authenticator did not properly decode a
   response
   (`#3205 <https://github.com/datalad/datalad/issues/3205>`__).

-  With Python 3 downloaders unnecessarily decoded the response when
   getting the status, leading to an encoding error.
   (`#3210 <https://github.com/datalad/datalad/issues/3210>`__)

-  In some cases, our internal command Runner did not adjust the
   environment’s ``PWD`` to match the current working directory
   specified with the ``cwd`` parameter.
   (`#3215 <https://github.com/datalad/datalad/issues/3215>`__)

-  The specification of the pyliblzma dependency was broken.
   (`#3220 <https://github.com/datalad/datalad/issues/3220>`__)

-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   displayed an uninformative blank log message in some cases.
   (`#3222 <https://github.com/datalad/datalad/issues/3222>`__)

-  The logic for finding the location of the aggregate metadata DB
   anchored the search path incorrectly, leading to a spurious warning.
   (`#3241 <https://github.com/datalad/datalad/issues/3241>`__)

-  Some progress bars were still displayed when stdout and stderr were
   not attached to a tty.
   (`#3281 <https://github.com/datalad/datalad/issues/3281>`__)

-  Check for stdin/out/err to not be closed before checking for
   ``.isatty``.
   (`#3268 <https://github.com/datalad/datalad/issues/3268>`__)

.. _enhancements-and-new-features-25:

Enhancements and new features
-----------------------------

-  Creating a new repository now aborts if any of the files in the
   directory are tracked by a repository in a parent directory.
   (`#3211 <https://github.com/datalad/datalad/issues/3211>`__)

-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   learned to replace the ``{tmpdir}`` placeholder in commands with a
   temporary directory.
   (`#3223 <https://github.com/datalad/datalad/issues/3223>`__)

-  `duecredit <https://github.com/duecredit/duecredit>`__ support has
   been added for citing DataLad itself as well as datasets that an
   analysis uses.
   (`#3184 <https://github.com/datalad/datalad/issues/3184>`__)

-  The ``eval_results`` interface helper unintentionally modified one of
   its arguments.
   (`#3249 <https://github.com/datalad/datalad/issues/3249>`__)

-  A few DataLad constants have been added, changed, or renamed
   (`#3250 <https://github.com/datalad/datalad/issues/3250>`__):

   -  ``HANDLE_META_DIR`` is now ``DATALAD_DOTDIR``. The old name should
      be considered deprecated.
   -  ``METADATA_DIR`` now refers to ``DATALAD_DOTDIR/metadata`` rather
      than ``DATALAD_DOTDIR/meta`` (which is still available as
      ``OLDMETADATA_DIR``).
   -  The new ``DATASET_METADATA_FILE`` refers to
      ``METADATA_DIR/dataset.json``.
   -  The new ``DATASET_CONFIG_FILE`` refers to
      ``DATALAD_DOTDIR/config``.
   -  ``METADATA_FILENAME`` has been renamed to
      ``OLDMETADATA_FILENAME``.

0.11.3 (Feb 19, 2019) – read-me-gently
======================================

Just a few of important fixes and minor enhancements.

.. _fixes-31:

Fixes
-----

-  The logic for setting the maximum command line length now works
   around Python 3.4 returning an unreasonably high value for
   ``SC_ARG_MAX`` on Debian systems.
   (`#3165 <https://github.com/datalad/datalad/issues/3165>`__)

-  DataLad commands that are conceptually “read-only”, such as
   ``datalad ls -L``, can fail when the caller lacks write permissions
   because git-annex tries merging remote git-annex branches to update
   information about availability. DataLad now disables
   ``annex.merge-annex-branches`` in some common “read-only” scenarios
   to avoid these failures.
   (`#3164 <https://github.com/datalad/datalad/issues/3164>`__)

.. _enhancements-and-new-features-26:

Enhancements and new features
-----------------------------

-  Accessing an “unbound” dataset method now automatically imports the
   necessary module rather than requiring an explicit import from the
   Python caller. For example, calling ``Dataset.add`` no longer needs
   to be preceded by ``from datalad.distribution.add import Add`` or an
   import of ``datalad.api``.
   (`#3156 <https://github.com/datalad/datalad/issues/3156>`__)

-  Configuring the new variable ``datalad.ssh.identityfile`` instructs
   DataLad to pass a value to the ``-i`` option of ``ssh``.
   (`#3149 <https://github.com/datalad/datalad/issues/3149>`__)
   (`#3168 <https://github.com/datalad/datalad/issues/3168>`__)

0.11.2 (Feb 07, 2019) – live-long-and-prosper
=============================================

A variety of bugfixes and enhancements

.. _major-refactoring-and-deprecations-10:

Major refactoring and deprecations
----------------------------------

-  All extracted metadata is now placed under git-annex by default.
   Previously files smaller than 20 kb were stored in git.
   (`#3109 <https://github.com/datalad/datalad/issues/3109>`__)
-  The function ``datalad.cmd.get_runner`` has been removed.
   (`#3104 <https://github.com/datalad/datalad/issues/3104>`__)

.. _fixes-32:

Fixes
-----

-  Improved handling of long commands:

   -  The code that inspected ``SC_ARG_MAX`` didn’t check that the
      reported value was a sensible, positive number.
      (`#3025 <https://github.com/datalad/datalad/issues/3025>`__)
   -  More commands that invoke ``git`` and ``git-annex`` with file
      arguments learned to split up the command calls when it is likely
      that the command would fail due to exceeding the maximum supported
      length.
      (`#3138 <https://github.com/datalad/datalad/issues/3138>`__)

-  The ``setup_yoda_dataset`` procedure created a malformed
   .gitattributes line.
   (`#3057 <https://github.com/datalad/datalad/issues/3057>`__)
-  `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   unnecessarily tried to infer the dataset when ``--no-save`` was
   given. (`#3029 <https://github.com/datalad/datalad/issues/3029>`__)
-  `rerun <https://datalad.readthedocs.io/en/latest/generated/man/datalad-rerun.html>`__
   aborted too late and with a confusing message when a ref specified
   via ``--onto`` didn’t exist.
   (`#3019 <https://github.com/datalad/datalad/issues/3019>`__)
-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__:

   -  ``run`` didn’t preserve the current directory prefix (“./”) on
      inputs and outputs, which is problematic if the caller relies on
      this representation when formatting the command.
      (`#3037 <https://github.com/datalad/datalad/issues/3037>`__)
   -  Fixed a number of unicode py2-compatibility issues.
      (`#3035 <https://github.com/datalad/datalad/issues/3035>`__)
      (`#3046 <https://github.com/datalad/datalad/issues/3046>`__)
   -  To proceed with a failed command, the user was confusingly
      instructed to use ``save`` instead of ``add`` even though ``run``
      uses ``add`` underneath.
      (`#3080 <https://github.com/datalad/datalad/issues/3080>`__)

-  Fixed a case where the helper class for checking external modules
   incorrectly reported a module as unknown.
   (`#3051 <https://github.com/datalad/datalad/issues/3051>`__)
-  `add-archive-content <https://datalad.readthedocs.io/en/latest/generated/man/datalad-add-archive-content.html>`__
   mishandled the archive path when the leading path contained a
   symlink. (`#3058 <https://github.com/datalad/datalad/issues/3058>`__)
-  Following denied access, the credential code failed to consider a
   scenario, leading to a type error rather than an appropriate error
   message. (`#3091 <https://github.com/datalad/datalad/issues/3091>`__)
-  Some tests failed when executed from a ``git worktree`` checkout of
   the source repository.
   (`#3129 <https://github.com/datalad/datalad/issues/3129>`__)
-  During metadata extraction, batched annex processes weren’t properly
   terminated, leading to issues on Windows.
   (`#3137 <https://github.com/datalad/datalad/issues/3137>`__)
-  `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__
   incorrectly handled an “invalid repository” exception when trying to
   add a submodule.
   (`#3141 <https://github.com/datalad/datalad/issues/3141>`__)
-  Pass ``GIT_SSH_VARIANT=ssh`` to git processes to be able to specify
   alternative ports in SSH urls

.. _enhancements-and-new-features-27:

Enhancements and new features
-----------------------------

-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   learned to suggest closely matching keys if there are no hits.
   (`#3089 <https://github.com/datalad/datalad/issues/3089>`__)
-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__

   -  gained a ``--group`` option so that the caller can specify the
      file system group for the repository.
      (`#3098 <https://github.com/datalad/datalad/issues/3098>`__)
   -  now understands SSH URLs that have a port in them (i.e. the
      “ssh://[user@]host.xz[:port]/path/to/repo.git/” syntax mentioned
      in ``man git-fetch``).
      (`#3146 <https://github.com/datalad/datalad/issues/3146>`__)

-  Interface classes can now override the default renderer for
   summarizing results.
   (`#3061 <https://github.com/datalad/datalad/issues/3061>`__)
-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__:

   -  ``--input`` and ``--output`` can now be shortened to ``-i`` and
      ``-o``.
      (`#3066 <https://github.com/datalad/datalad/issues/3066>`__)
   -  Placeholders such as “{inputs}” are now expanded in the command
      that is shown in the commit message subject.
      (`#3065 <https://github.com/datalad/datalad/issues/3065>`__)
   -  ``interface.run.run_command`` gained an ``extra_inputs`` argument
      so that wrappers like
      `datalad-container <https://github.com/datalad/datalad-container>`__
      can specify additional inputs that aren’t considered when
      formatting the command string.
      (`#3038 <https://github.com/datalad/datalad/issues/3038>`__)
   -  “–” can now be used to separate options for ``run`` and those for
      the command in ambiguous cases.
      (`#3119 <https://github.com/datalad/datalad/issues/3119>`__)

-  The utilities ``create_tree`` and ``ok_file_has_content`` now support
   “.gz” files.
   (`#3049 <https://github.com/datalad/datalad/issues/3049>`__)
-  The Singularity container for 0.11.1 now uses
   `nd_freeze <https://github.com/neurodebian/neurodebian/blob/master/tools/nd_freeze>`__
   to make its builds reproducible.
-  A
   `publications <https://datalad.readthedocs.io/en/latest/publications.html>`__
   page has been added to the documentation.
   (`#3099 <https://github.com/datalad/datalad/issues/3099>`__)
-  ``GitRepo.set_gitattributes`` now accepts a ``mode`` argument that
   controls whether the .gitattributes file is appended to (default) or
   overwritten.
   (`#3115 <https://github.com/datalad/datalad/issues/3115>`__)
-  ``datalad --help`` now avoids using ``man`` so that the list of
   subcommands is shown.
   (`#3124 <https://github.com/datalad/datalad/issues/3124>`__)

0.11.1 (Nov 26, 2018) – v7-better-than-v6
=========================================

Rushed out bugfix release to stay fully compatible with recent
`git-annex <http://git-annex.branchable.com/>`__ which introduced v7 to
replace v6.

.. _fixes-33:

Fixes
-----

-  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__:
   be able to install recursively into a dataset
   (`#2982 <https://github.com/datalad/datalad/issues/2982>`__)
-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__:
   be able to commit/save changes whenever files potentially could have
   swapped their storage between git and annex
   (`#1651 <https://github.com/datalad/datalad/issues/1651>`__)
   (`#2752 <https://github.com/datalad/datalad/issues/2752>`__)
   (`#3009 <https://github.com/datalad/datalad/issues/3009>`__)
-  [aggregate-metadata][]:

   -  dataset’s itself is now not “aggregated” if specific paths are
      provided for aggregation
      (`#3002 <https://github.com/datalad/datalad/issues/3002>`__). That
      resolves the issue of ``-r`` invocation aggregating all
      subdatasets of the specified dataset as well
   -  also compare/verify the actual content checksum of aggregated
      metadata while considering subdataset metadata for re-aggregation
      (`#3007 <https://github.com/datalad/datalad/issues/3007>`__)

-  ``annex`` commands are now chunked assuming 50% “safety margin” on
   the maximal command line length. Should resolve crashes while
   operating ot too many files at ones
   (`#3001 <https://github.com/datalad/datalad/issues/3001>`__)
-  ``run`` sidecar config processing
   (`#2991 <https://github.com/datalad/datalad/issues/2991>`__)
-  no double trailing period in docs
   (`#2984 <https://github.com/datalad/datalad/issues/2984>`__)
-  correct identification of the repository with symlinks in the paths
   in the tests
   (`#2972 <https://github.com/datalad/datalad/issues/2972>`__)
-  re-evaluation of dataset properties in case of dataset changes
   (`#2946 <https://github.com/datalad/datalad/issues/2946>`__)
-  [text2git][] procedure to use ``ds.repo.set_gitattributes``
   (`#2974 <https://github.com/datalad/datalad/issues/2974>`__)
   (`#2954 <https://github.com/datalad/datalad/issues/2954>`__)
-  Switch to use plain ``os.getcwd()`` if inconsistency with env var
   ``$PWD`` is detected
   (`#2914 <https://github.com/datalad/datalad/issues/2914>`__)
-  Make sure that credential defined in env var takes precedence
   (`#2960 <https://github.com/datalad/datalad/issues/2960>`__)
   (`#2950 <https://github.com/datalad/datalad/issues/2950>`__)

.. _enhancements-and-new-features-28:

Enhancements and new features
-----------------------------

-  `shub://datalad/datalad:git-annex-dev <https://singularity-hub.org/containers/5663/view>`__
   provides a Debian buster Singularity image with build environment for
   `git-annex <http://git-annex.branchable.com/>`__.
   ``tools/bisect-git-annex`` provides a helper for running
   ``git bisect`` on git-annex using that Singularity container
   (`#2995 <https://github.com/datalad/datalad/issues/2995>`__)
-  Added ``.zenodo.json`` for better integration with Zenodo for
   citation
-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   now provides names and help messages with a custom renderer for
   (`#2993 <https://github.com/datalad/datalad/issues/2993>`__)
-  Documentation: point to
   `datalad-revolution <http://github.com/datalad/datalad-revolution>`__
   extension (prototype of the greater DataLad future)
-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__

   -  support injecting of a detached command
      (`#2937 <https://github.com/datalad/datalad/issues/2937>`__)

-  ``annex`` metadata extractor now extracts ``annex.key`` metadata
   record. Should allow now to identify uses of specific files etc
   (`#2952 <https://github.com/datalad/datalad/issues/2952>`__)
-  Test that we can install from http://datasets.datalad.org
-  Proper rendering of ``CommandError`` (e.g. in case of “out of space”
   error) (`#2958 <https://github.com/datalad/datalad/issues/2958>`__)

0.11.0 (Oct 23, 2018) – Soon-to-be-perfect
==========================================

`git-annex <http://git-annex.branchable.com/>`__ 6.20180913 (or later)
is now required - provides a number of fixes for v6 mode operations etc.

.. _major-refactoring-and-deprecations-11:

Major refactoring and deprecations
----------------------------------

-  ``datalad.consts.LOCAL_CENTRAL_PATH`` constant was deprecated in
   favor of ``datalad.locations.default-dataset``
   `configuration <http://docs.datalad.org/en/latest/config.html>`__
   variable (`#2835 <https://github.com/datalad/datalad/issues/2835>`__)

Minor refactoring
-----------------

-  ``"notneeded"`` messages are no longer reported by default results
   renderer
-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   no longer shows commit instructions upon command failure when
   ``explicit`` is true and no outputs are specified
   (`#2922 <https://github.com/datalad/datalad/issues/2922>`__)
-  ``get_git_dir`` moved into GitRepo
   (`#2886 <https://github.com/datalad/datalad/issues/2886>`__)
-  ``_gitpy_custom_call`` removed from GitRepo
   (`#2894 <https://github.com/datalad/datalad/issues/2894>`__)
-  ``GitRepo.get_merge_base`` argument is now called ``commitishes``
   instead of ``treeishes``
   (`#2903 <https://github.com/datalad/datalad/issues/2903>`__)

.. _fixes-34:

Fixes
-----

-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   should not leave the dataset in non-clean state
   (`#2858 <https://github.com/datalad/datalad/issues/2858>`__) and some
   other enhancements
   (`#2859 <https://github.com/datalad/datalad/issues/2859>`__)
-  Fixed chunking of the long command lines to account for decorators
   and other arguments
   (`#2864 <https://github.com/datalad/datalad/issues/2864>`__)
-  Progress bar should not crash the process on some missing progress
   information
   (`#2891 <https://github.com/datalad/datalad/issues/2891>`__)
-  Default value for ``jobs`` set to be ``"auto"`` (not ``None``) to
   take advantage of possible parallel get if in ``-g`` mode
   (`#2861 <https://github.com/datalad/datalad/issues/2861>`__)
-  `wtf <http://datalad.readthedocs.io/en/latest/generated/man/datalad-wtf.html>`__
   must not crash if ``git-annex`` is not installed etc
   (`#2865 <https://github.com/datalad/datalad/issues/2865>`__),
   (`#2865 <https://github.com/datalad/datalad/issues/2865>`__),
   (`#2918 <https://github.com/datalad/datalad/issues/2918>`__),
   (`#2917 <https://github.com/datalad/datalad/issues/2917>`__)
-  Fixed paths (with spaces etc) handling while reporting annex error
   output (`#2892 <https://github.com/datalad/datalad/issues/2892>`__),
   (`#2893 <https://github.com/datalad/datalad/issues/2893>`__)
-  ``__del__`` should not access ``.repo`` but ``._repo`` to avoid
   attempts for reinstantiation etc
   (`#2901 <https://github.com/datalad/datalad/issues/2901>`__)
-  Fix up submodule ``.git`` right in ``GitRepo.add_submodule`` to avoid
   added submodules being non git-annex friendly
   (`#2909 <https://github.com/datalad/datalad/issues/2909>`__),
   (`#2904 <https://github.com/datalad/datalad/issues/2904>`__)
-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   (`#2905 <https://github.com/datalad/datalad/issues/2905>`__)

   -  now will provide dataset into the procedure if called within
      dataset
   -  will not crash if procedure is an executable without ``.py`` or
      ``.sh`` suffixes

-  Use centralized ``.gitattributes`` handling while setting annex
   backend (`#2912 <https://github.com/datalad/datalad/issues/2912>`__)
-  ``GlobbedPaths.expand(..., full=True)`` incorrectly returned relative
   paths when called more than once
   (`#2921 <https://github.com/datalad/datalad/issues/2921>`__)

.. _enhancements-and-new-features-29:

Enhancements and new features
-----------------------------

-  Report progress on
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   when installing from “smart” git servers
   (`#2876 <https://github.com/datalad/datalad/issues/2876>`__)
-  Stale/unused ``sth_like_file_has_content`` was removed
   (`#2860 <https://github.com/datalad/datalad/issues/2860>`__)
-  Enhancements to
   `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   to operate on “improved” metadata layouts
   (`#2878 <https://github.com/datalad/datalad/issues/2878>`__)
-  Output of ``git annex init`` operation is now logged
   (`#2881 <https://github.com/datalad/datalad/issues/2881>`__)
-  New

   -  ``GitRepo.cherry_pick``
      (`#2900 <https://github.com/datalad/datalad/issues/2900>`__)
   -  ``GitRepo.format_commit``
      (`#2902 <https://github.com/datalad/datalad/issues/2902>`__)

-  `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   (`#2905 <https://github.com/datalad/datalad/issues/2905>`__)

   -  procedures can now recursively be discovered in subdatasets as
      well. The uppermost has highest priority
   -  Procedures in user and system locations now take precedence over
      those in datasets.

0.10.3.1 (Sep 13, 2018) – Nothing-is-perfect
============================================

Emergency bugfix to address forgotten boost of version in
``datalad/version.py``.

0.10.3 (Sep 13, 2018) – Almost-perfect
======================================

This is largely a bugfix release which addressed many (but not yet all)
issues of working with git-annex direct and version 6 modes, and
operation on Windows in general. Among enhancements you will see the
support of public S3 buckets (even with periods in their names), ability
to configure new providers interactively, and improved ``egrep`` search
backend.

Although we do not require with this release, it is recommended to make
sure that you are using a recent ``git-annex`` since it also had a
variety of fixes and enhancements in the past months.

.. _fixes-35:

Fixes
-----

-  Parsing of combined short options has been broken since DataLad
   v0.10.0. (`#2710 <https://github.com/datalad/datalad/issues/2710>`__)
-  The ``datalad save`` instructions shown by ``datalad run`` for a
   command with a non-zero exit were incorrectly formatted.
   (`#2692 <https://github.com/datalad/datalad/issues/2692>`__)
-  Decompression of zip files (e.g., through
   ``datalad add-archive-content``) failed on Python 3.
   (`#2702 <https://github.com/datalad/datalad/issues/2702>`__)
-  Windows:

   -  colored log output was not being processed by colorama.
      (`#2707 <https://github.com/datalad/datalad/issues/2707>`__)
   -  more codepaths now try multiple times when removing a file to deal
      with latency and locking issues on Windows.
      (`#2795 <https://github.com/datalad/datalad/issues/2795>`__)

-  Internal git fetch calls have been updated to work around a GitPython
   ``BadName`` issue.
   (`#2712 <https://github.com/datalad/datalad/issues/2712>`__),
   (`#2794 <https://github.com/datalad/datalad/issues/2794>`__)
-  The progess bar for annex file transferring was unable to handle an
   empty file.
   (`#2717 <https://github.com/datalad/datalad/issues/2717>`__)
-  ``datalad add-readme`` halted when no aggregated metadata was found
   rather than displaying a warning.
   (`#2731 <https://github.com/datalad/datalad/issues/2731>`__)
-  ``datalad rerun`` failed if ``--onto`` was specified and the history
   contained no run commits.
   (`#2761 <https://github.com/datalad/datalad/issues/2761>`__)
-  Processing of a command’s results failed on a result record with a
   missing value (e.g., absent field or subfield in metadata). Now the
   missing value is rendered as “N/A”.
   (`#2725 <https://github.com/datalad/datalad/issues/2725>`__).
-  A couple of documentation links in the “Delineation from related
   solutions” were misformatted.
   (`#2773 <https://github.com/datalad/datalad/issues/2773>`__)
-  With the latest git-annex, several known V6 failures are no longer an
   issue. (`#2777 <https://github.com/datalad/datalad/issues/2777>`__)
-  In direct mode, commit changes would often commit annexed content as
   regular Git files. A new approach fixes this and resolves a good
   number of known failures.
   (`#2770 <https://github.com/datalad/datalad/issues/2770>`__)
-  The reporting of command results failed if the current working
   directory was removed (e.g., after an unsuccessful ``install``).
   (`#2788 <https://github.com/datalad/datalad/issues/2788>`__)
-  When installing into an existing empty directory, ``datalad install``
   removed the directory after a failed clone.
   (`#2788 <https://github.com/datalad/datalad/issues/2788>`__)
-  ``datalad run`` incorrectly handled inputs and outputs for paths with
   spaces and other characters that require shell escaping.
   (`#2798 <https://github.com/datalad/datalad/issues/2798>`__)
-  Globbing inputs and outputs for ``datalad run`` didn’t work correctly
   if a subdataset wasn’t installed.
   (`#2796 <https://github.com/datalad/datalad/issues/2796>`__)
-  Minor (in)compatibility with git 2.19 - (no) trailing period in an
   error message now.
   (`#2815 <https://github.com/datalad/datalad/issues/2815>`__)

.. _enhancements-and-new-features-30:

Enhancements and new features
-----------------------------

-  Anonymous access is now supported for S3 and other downloaders.
   (`#2708 <https://github.com/datalad/datalad/issues/2708>`__)
-  A new interface is available to ease setting up new providers.
   (`#2708 <https://github.com/datalad/datalad/issues/2708>`__)
-  Metadata: changes to egrep mode search
   (`#2735 <https://github.com/datalad/datalad/issues/2735>`__)

   -  Queries in egrep mode are now case-sensitive when the query
      contains any uppercase letters and are case-insensitive otherwise.
      The new mode egrepcs can be used to perform a case-sensitive query
      with all lower-case letters.
   -  Search can now be limited to a specific key.
   -  Multiple queries (list of expressions) are evaluated using AND to
      determine whether something is a hit.
   -  A single multi-field query (e.g., ``pa*:findme``) is a hit, when
      any matching field matches the query.
   -  All matching key/value combinations across all (multi-field)
      queries are reported in the query_matched result field.
   -  egrep mode now shows all hits rather than limiting the results to
      the top 20 hits.

-  The documentation on how to format commands for ``datalad run`` has
   been improved.
   (`#2703 <https://github.com/datalad/datalad/issues/2703>`__)
-  The method for determining the current working directory on Windows
   has been improved.
   (`#2707 <https://github.com/datalad/datalad/issues/2707>`__)
-  ``datalad --version`` now simply shows the version without the
   license. (`#2733 <https://github.com/datalad/datalad/issues/2733>`__)
-  ``datalad export-archive`` learned to export under an existing
   directory via its ``--filename`` option.
   (`#2723 <https://github.com/datalad/datalad/issues/2723>`__)
-  ``datalad export-to-figshare`` now generates the zip archive in the
   root of the dataset unless ``--filename`` is specified.
   (`#2723 <https://github.com/datalad/datalad/issues/2723>`__)
-  After importing ``datalad.api``, ``help(datalad.api)`` (or
   ``datalad.api?`` in IPython) now shows a summary of the available
   DataLad commands.
   (`#2728 <https://github.com/datalad/datalad/issues/2728>`__)
-  Support for using ``datalad`` from IPython has been improved.
   (`#2722 <https://github.com/datalad/datalad/issues/2722>`__)
-  ``datalad wtf`` now returns structured data and reports the version
   of each extension.
   (`#2741 <https://github.com/datalad/datalad/issues/2741>`__)
-  The internal handling of gitattributes information has been improved.
   A user-visible consequence is that ``datalad create --force`` no
   longer duplicates existing attributes.
   (`#2744 <https://github.com/datalad/datalad/issues/2744>`__)
-  The “annex” metadata extractor can now be used even when no content
   is present.
   (`#2724 <https://github.com/datalad/datalad/issues/2724>`__)
-  The ``add_url_to_file`` method (called by commands like
   ``datalad download-url`` and ``datalad add-archive-content``) learned
   how to display a progress bar.
   (`#2738 <https://github.com/datalad/datalad/issues/2738>`__)

0.10.2 (Jul 09, 2018) – Thesecuriestever
========================================

Primarily a bugfix release to accommodate recent git-annex release
forbidding file:// and http://localhost/ URLs which might lead to
revealing private files if annex is publicly shared.

.. _fixes-36:

Fixes
-----

-  fixed testing to be compatible with recent git-annex (6.20180626)
-  `download-url <https://datalad.readthedocs.io/en/latest/generated/man/datalad-download-url.html>`__
   will now download to current directory instead of the top of the
   dataset

.. _enhancements-and-new-features-31:

Enhancements and new features
-----------------------------

-  do not quote ~ in URLs to be consistent with quote implementation in
   Python 3.7 which now follows RFC 3986
-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   support for user-configured placeholder values
-  documentation on native git-annex metadata support
-  handle 401 errors from LORIS tokens
-  ``yoda`` procedure will instantiate ``README.md``
-  ``--discover`` option added to
   `run-procedure <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run-procedure.html>`__
   to list available procedures

0.10.1 (Jun 17, 2018) – OHBM polish
===================================

The is a minor bugfix release.

.. _fixes-37:

Fixes
-----

-  Be able to use backports.lzma as a drop-in replacement for pyliblzma.
-  Give help when not specifying a procedure name in ``run-procedure``.
-  Abort early when a downloader received no filename.
-  Avoid ``rerun`` error when trying to unlock non-available files.

0.10.0 (Jun 09, 2018) – The Release
===================================

This release is a major leap forward in metadata support.

.. _major-refactoring-and-deprecations-12:

Major refactoring and deprecations
----------------------------------

-  Metadata

   -  Prior metadata provided by datasets under ``.datalad/meta`` is no
      longer used or supported. Metadata must be reaggregated using 0.10
      version
   -  Metadata extractor types are no longer auto-guessed and must be
      explicitly specified in ``datalad.metadata.nativetype`` config
      (could contain multiple values)
   -  Metadata aggregation of a dataset hierarchy no longer updates all
      datasets in the tree with new metadata. Instead, only the target
      dataset is updated. This behavior can be changed via the
      –update-mode switch. The new default prevents needless
      modification of (3rd-party) subdatasets.
   -  Neuroimaging metadata support has been moved into a dedicated
      extension: https://github.com/datalad/datalad-neuroimaging

-  Crawler

   -  moved into a dedicated extension:
      https://github.com/datalad/datalad-crawler

-  ``export_tarball`` plugin has been generalized to ``export_archive``
   and can now also generate ZIP archives.
-  By default a dataset X is now only considered to be a super-dataset
   of another dataset Y, if Y is also a registered subdataset of X.

.. _fixes-38:

Fixes
-----

A number of fixes did not make it into the 0.9.x series:

-  Dynamic configuration overrides via the ``-c`` option were not in
   effect.
-  ``save`` is now more robust with respect to invocation in
   subdirectories of a dataset.
-  ``unlock`` now reports correct paths when running in a dataset
   subdirectory.
-  ``get`` is more robust to path that contain symbolic links.
-  symlinks to subdatasets of a dataset are now correctly treated as a
   symlink, and not as a subdataset
-  ``add`` now correctly saves staged subdataset additions.
-  Running ``datalad save`` in a dataset no longer adds untracked
   content to the dataset. In order to add content a path has to be
   given, e.g. \ ``datalad save .``
-  ``wtf`` now works reliably with a DataLad that wasn’t installed from
   Git (but, e.g., via pip)
-  More robust URL handling in ``simple_with_archives`` crawler
   pipeline.

.. _enhancements-and-new-features-32:

Enhancements and new features
-----------------------------

-  Support for DataLad extension that can contribute API components from
   3rd-party sources, incl. commands, metadata extractors, and test case
   implementations. See
   https://github.com/datalad/datalad-extension-template for a demo
   extension.
-  Metadata (everything has changed!)

   -  Metadata extraction and aggregation is now supported for datasets
      and individual files.
   -  Metadata query via ``search`` can now discover individual files.
   -  Extracted metadata can now be stored in XZ compressed files, is
      optionally annexed (when exceeding a configurable size threshold),
      and obtained on demand (new configuration option
      ``datalad.metadata.create-aggregate-annex-limit``).
   -  Status and availability of aggregated metadata can now be reported
      via ``metadata --get-aggregates``
   -  New configuration option ``datalad.metadata.maxfieldsize`` to
      exclude too large metadata fields from aggregation.
   -  The type of metadata is no longer guessed during metadata
      extraction. A new configuration option
      ``datalad.metadata.nativetype`` was introduced to enable one or
      more particular metadata extractors for a dataset.
   -  New configuration option
      ``datalad.metadata.store-aggregate-content`` to enable the storage
      of aggregated metadata for dataset content (i.e. file-based
      metadata) in contrast to just metadata describing a dataset as a
      whole.

-  ``search`` was completely reimplemented. It offers three different
   modes now:

   -  ‘egrep’ (default): expression matching in a plain string version
      of metadata
   -  ‘textblob’: search a text version of all metadata using a fully
      featured query language (fast indexing, good for keyword search)
   -  ‘autofield’: search an auto-generated index that preserves
      individual fields of metadata that can be represented in a tabular
      structure (substantial indexing cost, enables the most detailed
      queries of all modes)

-  New extensions:

   -  `addurls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-addurls.html>`__,
      an extension for creating a dataset (and possibly subdatasets)
      from a list of URLs.
   -  export_to_figshare
   -  extract_metadata

-  add_readme makes use of available metadata
-  By default the wtf extension now hides sensitive information, which
   can be included in the output by passing ``--senstive=some`` or
   ``--senstive=all``.
-  Reduced startup latency by only importing commands necessary for a
   particular command line call.
-  `create <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create.html>`__:

   -  ``-d <parent> --nosave`` now registers subdatasets, when possible.
   -  ``--fake-dates`` configures dataset to use fake-dates

-  `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   now provides a way for the caller to save the result when a command
   has a non-zero exit status.
-  ``datalad rerun`` now has a ``--script`` option that can be used to
   extract previous commands into a file.
-  A DataLad Singularity container is now available on `Singularity
   Hub <https://singularity-hub.org/collections/667>`__.
-  More casts have been embedded in the `use case section of the
   documentation <http://docs.datalad.org/en/docs/usecases/index.html>`__.
-  ``datalad --report-status`` has a new value ‘all’ that can be used to
   temporarily re-enable reporting that was disable by configuration
   settings.

0.9.3 (Mar 16, 2018) – pi+0.02 release
======================================

Some important bug fixes which should improve usability

.. _fixes-39:

Fixes
-----

-  ``datalad-archives`` special remote now will lock on acquiring or
   extracting an archive - this allows for it to be used with -J flag
   for parallel operation
-  relax introduced in 0.9.2 demand on git being configured for datalad
   operation - now we will just issue a warning
-  ``datalad ls`` should now list “authored date” and work also for
   datasets in detached HEAD mode
-  ``datalad save`` will now save original file as well, if file was
   “git mv”ed, so you can now ``datalad run git mv old new`` and have
   changes recorded

.. _enhancements-and-new-features-33:

Enhancements and new features
-----------------------------

-  ``--jobs`` argument now could take ``auto`` value which would decide
   on # of jobs depending on the # of available CPUs. ``git-annex`` >
   6.20180314 is recommended to avoid regression with -J.
-  memoize calls to ``RI`` meta-constructor – should speed up operation
   a bit
-  ``DATALAD_SEED`` environment variable could be used to seed Python
   RNG and provide reproducible UUIDs etc (useful for testing and demos)

0.9.2 (Mar 04, 2018) – it is (again) better than ever
=====================================================

Largely a bugfix release with a few enhancements.

.. _fixes-40:

Fixes
-----

-  Execution of external commands (git) should not get stuck when lots
   of both stdout and stderr output, and should not loose remaining
   output in some cases
-  Config overrides provided in the command line (-c) should now be
   handled correctly
-  Consider more remotes (not just tracking one, which might be none)
   while installing subdatasets
-  Compatibility with git 2.16 with some changed behaviors/annotations
   for submodules
-  Fail ``remove`` if ``annex drop`` failed
-  Do not fail operating on files which start with dash (-)
-  URL unquote paths within S3, URLs and DataLad RIs (///)
-  In non-interactive mode fail if authentication/access fails
-  Web UI:

   -  refactored a little to fix incorrect listing of submodules in
      subdirectories
   -  now auto-focuses on search edit box upon entering the page

-  Assure that extracted from tarballs directories have executable bit
   set

.. _enhancements-and-new-features-34:

Enhancements and new features
-----------------------------

-  A log message and progress bar will now inform if a tarball to be
   downloaded while getting specific files (requires git-annex >
   6.20180206)
-  A dedicated ``datalad rerun`` command capable of rerunning entire
   sequences of previously ``run`` commands. **Reproducibility through
   VCS. Use ``run`` even if not interested in ``rerun``**
-  Alert the user if ``git`` is not yet configured but git operations
   are requested
-  Delay collection of previous ssh connections until it is actually
   needed. Also do not require ‘:’ while specifying ssh host
-  AutomagicIO: Added proxying of isfile, lzma.LZMAFile and io.open
-  Testing:

   -  added DATALAD_DATASETS_TOPURL=http://datasets-tests.datalad.org to
      run tests against another website to not obscure access stats
   -  tests run against temporary HOME to avoid side-effects
   -  better unit-testing of interactions with special remotes

-  CONTRIBUTING.md describes how to setup and use ``git-hub`` tool to
   “attach” commits to an issue making it into a PR
-  DATALAD_USE_DEFAULT_GIT env variable could be used to cause DataLad
   to use default (not the one possibly bundled with git-annex) git
-  Be more robust while handling not supported requests by annex in
   special remotes
-  Use of ``swallow_logs`` in the code was refactored away – less
   mysteries now, just increase logging level
-  ``wtf`` plugin will report more information about environment,
   externals and the system

0.9.1 (Oct 01, 2017) – “DATALAD!”(JBTM)
=======================================

Minor bugfix release

.. _fixes-41:

Fixes
-----

-  Should work correctly with subdatasets named as numbers of bool
   values (requires also GitPython >= 2.1.6)
-  Custom special remotes should work without crashing with git-annex >=
   6.20170924

0.9.0 (Sep 19, 2017) – isn’t it a lucky day even though not a Friday?
=====================================================================

.. _major-refactoring-and-deprecations-13:

Major refactoring and deprecations
----------------------------------

-  the ``files`` argument of
   `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   has been renamed to ``path`` to be uniform with any other command
-  all major commands now implement more uniform API semantics and
   result reporting. Functionality for modification detection of dataset
   content has been completely replaced with a more efficient
   implementation
-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   now features a ``--transfer-data`` switch that allows for a
   disambiguous specification of whether to publish data – independent
   of the selection which datasets to publish (which is done via their
   paths). Moreover,
   `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   now transfers data before repository content is pushed.

.. _fixes-42:

Fixes
-----

-  `drop <http://datalad.readthedocs.io/en/latest/generated/man/datalad-drop.html>`__
   no longer errors when some subdatasets are not installed
-  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   will no longer report nothing when a Dataset instance was given as a
   source argument, but rather perform as expected
-  `remove <http://datalad.readthedocs.io/en/latest/generated/man/datalad-remove.html>`__
   doesn’t remove when some files of a dataset could not be dropped
-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__

   -  no longer hides error during a repository push
   -  publish behaves “correctly” for ``--since=`` in considering only
      the differences the last “pushed” state
   -  data transfer handling while publishing with dependencies, to
      github

-  improved robustness with broken Git configuration
-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   should search for unicode strings correctly and not crash
-  robustify git-annex special remotes protocol handling to allow for
   spaces in the last argument
-  UI credentials interface should now allow to Ctrl-C the entry
-  should not fail while operating on submodules named with numerics
   only or by bool (true/false) names
-  crawl templates should not now override settings for ``largefiles``
   if specified in ``.gitattributes``

.. _enhancements-and-new-features-35:

Enhancements and new features
-----------------------------

-  **Exciting new feature**
   `run <http://datalad.readthedocs.io/en/latest/generated/man/datalad-run.html>`__
   command to protocol execution of an external command and rerun
   computation if desired. See
   `screencast <http://datalad.org/features.html#reproducible-science>`__
-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__
   now uses Git for detecting with sundatasets need to be inspected for
   potential changes, instead of performing a complete traversal of a
   dataset tree
-  `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__
   looks for changes relative to the last commited state of a dataset to
   discover files to add more efficiently
-  `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   can now report untracked files in addition to modified files
-  [uninstall][] will check itself whether a subdataset is properly
   registered in a superdataset, even when no superdataset is given in a
   call
-  `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   can now configure subdatasets for exclusion from recursive
   installation (``datalad-recursiveinstall`` submodule configuration
   property)
-  precrafted pipelines of [crawl][] now will not override
   ``annex.largefiles`` setting if any was set within ``.gitattribues``
   (e.g. by ``datalad create --text-no-annex``)
-  framework for screencasts: ``tools/cast*`` tools and sample cast
   scripts under ``doc/casts`` which are published at
   `datalad.org/features.html <http://datalad.org/features.html>`__
-  new `project YouTube
   channel <https://www.youtube.com/channel/UCB8-Zf7D0DSzAsREoIt0Bvw>`__
-  tests failing in direct and/or v6 modes marked explicitly

0.8.1 (Aug 13, 2017) – the best birthday gift
=============================================

Bugfixes

.. _fixes-43:

Fixes
-----

-  Do not attempt to
   `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   a not installed sub-dataset
-  In case of too many files to be specified for
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   or
   `copy_to <http://docs.datalad.org/en/latest/_modules/datalad/support/annexrepo.html?highlight=%22copy_to%22>`__,
   we will make multiple invocations of underlying git-annex command to
   not overfill command line
-  More robust handling of unicode output in terminals which might not
   support it

.. _enhancements-and-new-features-36:

Enhancements and new features
-----------------------------

-  Ship a copy of numpy.testing to facilitate [test][] without requiring
   numpy as dependency. Also allow to pass to command which test(s) to
   run
-  In
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   and
   `copy_to <http://docs.datalad.org/en/latest/_modules/datalad/support/annexrepo.html?highlight=%22copy_to%22>`__
   provide actual original requested paths, not the ones we deduced need
   to be transferred, solely for knowing the total

0.8.0 (Jul 31, 2017) – it is better than ever
=============================================

A variety of fixes and enhancements

.. _fixes-44:

Fixes
-----

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   would now push merged ``git-annex`` branch even if no other changes
   were done
-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   should be able to publish using relative path within SSH URI (git
   hook would use relative paths)
-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   should better tollerate publishing to pure git and ``git-annex``
   special remotes

.. _enhancements-and-new-features-37:

Enhancements and new features
-----------------------------

-  `plugin <http://datalad.readthedocs.io/en/latest/generated/man/datalad-plugin.html>`__
   mechanism came to replace
   `export <http://datalad.readthedocs.io/en/latest/generated/man/datalad-export.html>`__.
   See
   `export_tarball <http://docs.datalad.org/en/latest/generated/datalad.plugin.export_tarball.html>`__
   for the replacement of
   `export <http://datalad.readthedocs.io/en/latest/generated/man/datalad-export.html>`__.
   Now it should be easy to extend datalad’s interface with custom
   functionality to be invoked along with other commands.
-  Minimalistic coloring of the results rendering
-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__/``copy_to``
   got progress bar report now and support of ``--jobs``
-  minor fixes and enhancements to crawler (e.g. support of recursive
   removes)

0.7.0 (Jun 25, 2017) – when it works - it is quite awesome!
===========================================================

New features, refactorings, and bug fixes.

.. _major-refactoring-and-deprecations-14:

Major refactoring and deprecations
----------------------------------

-  `add-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add-sibling.html>`__
   has been fully replaced by the
   `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   command
-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__,
   and
   `unlock <http://datalad.readthedocs.io/en/latest/generated/man/datalad-unlock.html>`__
   have been re-written to support the same common API as most other
   commands

.. _enhancements-and-new-features-38:

Enhancements and new features
-----------------------------

-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   can now be used to query and configure a local repository by using
   the sibling name ``here``
-  `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   can now query and set annex preferred content configuration. This
   includes ``wanted`` (as previously supported in other commands), and
   now also ``required``
-  New
   `metadata <http://datalad.readthedocs.io/en/latest/generated/man/datalad-metadata.html>`__
   command to interface with datasets/files
   `meta-data <http://docs.datalad.org/en/latest/cmdline.html#meta-data-handling>`__
-  Documentation for all commands is now built in a uniform fashion
-  Significant parts of the documentation of been updated
-  Instantiate GitPython’s Repo instances lazily

.. _fixes-45:

Fixes
-----

-  API documentation is now rendered properly as HTML, and is easier to
   browse by having more compact pages
-  Closed files left open on various occasions (Popen PIPEs, etc)
-  Restored basic (consumer mode of operation) compatibility with
   Windows OS

0.6.0 (Jun 14, 2017) – German perfectionism
===========================================

This release includes a **huge** refactoring to make code base and
functionality more robust and flexible

-  outputs from API commands could now be highly customized. See
   ``--output-format``, ``--report-status``, ``--report-type``, and
   ``--report-type`` options for
   `datalad <http://docs.datalad.org/en/latest/generated/man/datalad.html>`__
   command.
-  effort was made to refactor code base so that underlying functions
   behave as generators where possible
-  input paths/arguments analysis was redone for majority of the
   commands to provide unified behavior

.. _major-refactoring-and-deprecations-15:

Major refactoring and deprecations
----------------------------------

-  ``add-sibling`` and ``rewrite-urls`` were refactored in favor of new
   `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   command which should be used for siblings manipulations
-  ‘datalad.api.alwaysrender’ config setting/support is removed in favor
   of new outputs processing

.. _fixes-46:

Fixes
-----

-  Do not flush manually git index in pre-commit to avoid “Death by the
   Lock” issue
-  Deployed by
   `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
   ``post-update`` hook script now should be more robust (tolerate
   directory names with spaces, etc.)
-  A variety of fixes, see `list of pull requests and issues
   closed <https://github.com/datalad/datalad/milestone/41?closed=1>`__
   for more information

.. _enhancements-and-new-features-39:

Enhancements and new features
-----------------------------

-  new
   `annotate-paths <http://docs.datalad.org/en/latest/generated/man/datalad-annotate-paths.html>`__
   plumbing command to inspect and annotate provided paths. Use
   ``--modified`` to summarize changes between different points in the
   history
-  new
   `clone <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clone.html>`__
   plumbing command to provide a subset (install a single dataset from a
   URL) functionality of
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
-  new
   `diff <http://datalad.readthedocs.io/en/latest/generated/man/datalad-diff.html>`__
   plumbing command
-  new
   `siblings <http://datalad.readthedocs.io/en/latest/generated/man/datalad-siblings.html>`__
   command to list or manipulate siblings
-  new
   `subdatasets <http://datalad.readthedocs.io/en/latest/generated/man/datalad-subdatasets.html>`__
   command to list subdatasets and their properties
-  `drop <http://datalad.readthedocs.io/en/latest/generated/man/datalad-drop.html>`__
   and
   `remove <http://datalad.readthedocs.io/en/latest/generated/man/datalad-remove.html>`__
   commands were refactored
-  ``benchmarks/`` collection of `Airspeed
   velocity <https://github.com/spacetelescope/asv/>`__ benchmarks
   initiated. See reports at http://datalad.github.io/datalad/
-  crawler would try to download a new url multiple times increasing
   delay between attempts. Helps to resolve problems with extended
   crawls of Amazon S3
-  `CRCNS <http://crcns.org>`__ crawler pipeline now also fetches and
   aggregates meta-data for the datasets from datacite
-  overall optimisations to benefit from the aforementioned refactoring
   and improve user-experience
-  a few stub and not (yet) implemented commands (e.g. ``move``) were
   removed from the interface
-  Web frontend got proper coloring for the breadcrumbs and some
   additional caching to speed up interactions. See
   http://datasets.datalad.org
-  Small improvements to the online documentation. See e.g. `summary of
   differences between
   git/git-annex/datalad <http://docs.datalad.org/en/latest/related.html#git-git-annex-datalad>`__

0.5.1 (Mar 25, 2017) – cannot stop the progress
===============================================

A bugfix release

.. _fixes-47:

Fixes
-----

-  `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__
   was forcing addition of files to annex regardless of settings in
   ``.gitattributes``. Now that decision is left to annex by default
-  ``tools/testing/run_doc_examples`` used to run doc examples as tests,
   fixed up to provide status per each example and not fail at once
-  ``doc/examples``

   -  `3rdparty_analysis_workflow.sh <http://docs.datalad.org/en/latest/generated/examples/3rdparty_analysis_workflow.html>`__
      was fixed up to reflect changes in the API of 0.5.0.

-  progress bars

   -  should no longer crash **datalad** and report correct sizes and
      speeds
   -  should provide progress reports while using Python 3.x

.. _enhancements-and-new-features-40:

Enhancements and new features
-----------------------------

-  ``doc/examples``

   -  `nipype_workshop_dataset.sh <http://docs.datalad.org/en/latest/generated/examples/nipype_workshop_dataset.html>`__
      new example to demonstrate how new super- and sub- datasets were
      established as a part of our datasets collection

0.5.0 (Mar 20, 2017) – it’s huge
================================

This release includes an avalanche of bug fixes, enhancements, and
additions which at large should stay consistent with previous behavior
but provide better functioning. Lots of code was refactored to provide
more consistent code-base, and some API breakage has happened. Further
work is ongoing to standardize output and results reporting
(`#1350 <https://github.com/datalad/datalad/issues/1350>`__)

Most notable changes
--------------------

-  requires `git-annex <http://git-annex.branchable.com/>`__ >=
   6.20161210 (or better even >= 6.20161210 for improved functionality)
-  commands should now operate on paths specified (if any), without
   causing side-effects on other dirty/staged files
-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__

   -  ``-a`` is deprecated in favor of ``-u`` or ``--all-updates`` so
      only changes known components get saved, and no new files
      automagically added
   -  ``-S`` does no longer store the originating dataset in its commit
      message

-  `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__

   -  can specify commit/save message with ``-m``

-  `add-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add-sibling.html>`__
   and
   `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__

   -  now take the name of the sibling (remote) as a ``-s`` (``--name``)
      option, not a positional argument
   -  ``--publish-depends`` to setup publishing data and code to
      multiple repositories (e.g. github + webserve) should now be
      functional see `this
      comment <https://github.com/datalad/datalad/issues/335#issuecomment-277240733>`__
   -  got ``--publish-by-default`` to specify what refs should be
      published by default
   -  got ``--annex-wanted``, ``--annex-groupwanted`` and
      ``--annex-group`` settings which would be used to instruct annex
      about preferred content.
      `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__
      then will publish data using those settings if ``wanted`` is set.
   -  got ``--inherit`` option to automagically figure out url/wanted
      and other git/annex settings for new remote sub-dataset to be
      constructed

-  `publish <http://datalad.readthedocs.io/en/latest/generated/man/datalad-publish.html>`__

   -  got ``--skip-failing`` refactored into ``--missing`` option which
      could use new feature of
      `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
      ``--inherit``

.. _fixes-48:

Fixes
-----

-  More consistent interaction through ssh - all ssh connections go
   through
   `sshrun <http://datalad.readthedocs.io/en/latest/generated/man/datalad-sshrun.html>`__
   shim for a “single point of authentication”, etc.
-  More robust
   `ls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-ls.html>`__
   operation outside of the datasets
-  A number of fixes for direct and v6 mode of annex

.. _enhancements-and-new-features-41:

Enhancements and new features
-----------------------------

-  New
   `drop <http://datalad.readthedocs.io/en/latest/generated/man/datalad-drop.html>`__
   and
   `remove <http://datalad.readthedocs.io/en/latest/generated/man/datalad-remove.html>`__
   commands
-  `clean <http://datalad.readthedocs.io/en/latest/generated/man/datalad-clean.html>`__

   -  got ``--what`` to specify explicitly what cleaning steps to
      perform and now could be invoked with ``-r``

-  ``datalad`` and ``git-annex-remote*`` scripts now do not use
   setuptools entry points mechanism and rely on simple import to
   shorten start up time
-  `Dataset <http://docs.datalad.org/en/latest/generated/datalad.api.Dataset.html>`__
   is also now using `Flyweight
   pattern <https://en.wikipedia.org/wiki/Flyweight_pattern>`__, so the
   same instance is reused for the same dataset
-  progressbars should not add more empty lines

Internal refactoring
--------------------

-  Majority of the commands now go through ``_prep`` for arguments
   validation and pre-processing to avoid recursive invocations

0.4.1 (Nov 10, 2016) – CA release
=================================

Requires now GitPython >= 2.1.0

.. _fixes-49:

Fixes
-----

-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__

   -  to not save staged files if explicit paths were provided

-  improved (but not yet complete) support for direct mode
-  `update <http://datalad.readthedocs.io/en/latest/generated/man/datalad-update.html>`__
   to not crash if some sub-datasets are not installed
-  do not log calls to ``git config`` to avoid leakage of possibly
   sensitive settings to the logs

.. _enhancements-and-new-features-42:

Enhancements and new features
-----------------------------

-  New `rfc822-compliant
   metadata <http://docs.datalad.org/en/latest/metadata.html#rfc822-compliant-meta-data>`__
   format
-  `save <http://datalad.readthedocs.io/en/latest/generated/man/datalad-save.html>`__

   -  -S to save the change also within all super-datasets

-  `add <http://datalad.readthedocs.io/en/latest/generated/man/datalad-add.html>`__
   now has progress-bar reporting
-  `create-sibling-github <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling-github.html>`__
   to create a :term:``sibling`` of a dataset on github
-  `OpenfMRI <http://openfmri.org>`__ crawler and datasets were enriched
   with URLs to separate files where also available from openfmri s3
   bucket (if upgrading your datalad datasets, you might need to run
   ``git annex enableremote datalad`` to make them available)
-  various enhancements to log messages
-  web interface

   -  populates “install” box first thus making UX better over slower
      connections

0.4 (Oct 22, 2016) – Paris is waiting
=====================================

Primarily it is a bugfix release but because of significant refactoring
of the
`install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
and
`get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
implementation, it gets a new minor release.

.. _fixes-50:

Fixes
-----

-  be able to
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
   or
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   while providing paths while being outside of a dataset
-  remote annex datasets get properly initialized
-  robust detection of outdated
   `git-annex <http://git-annex.branchable.com/>`__

.. _enhancements-and-new-features-43:

Enhancements and new features
-----------------------------

-  interface changes

   -  `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
      ``--recursion-limit=existing`` to not recurse into not-installed
      subdatasets
   -  `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
      ``-n`` to possibly install sub-datasets without getting any data
   -  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
      ``--jobs|-J`` to specify number of parallel jobs for annex
      `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
      call could use (ATM would not work when data comes from archives)

-  more (unit-)testing
-  documentation: see http://docs.datalad.org/en/latest/basics.html for
   basic principles and useful shortcuts in referring to datasets
-  various webface improvements: breadcrumb paths, instructions how to
   install dataset, show version from the tags, etc.

0.3.1 (Oct 1, 2016) – what a wonderful week
===========================================

Primarily bugfixes but also a number of enhancements and core
refactorings

.. _fixes-51:

Fixes
-----

-  do not build manpages and examples during installation to avoid
   problems with possibly previously outdated dependencies
-  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   can be called on already installed dataset (with ``-r`` or ``-g``)

.. _enhancements-and-new-features-44:

Enhancements and new features
-----------------------------

-  complete overhaul of datalad configuration settings handling (see
   `Configuration
   documentation <http://docs.datalad.org/config.html>`__), so majority
   of the environment. Now uses git format and stores persistent
   configuration settings under ``.datalad/config`` and local within
   ``.git/config`` variables we have used were renamed to match
   configuration names
-  `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__
   does not now by default upload web front-end
-  `export <http://datalad.readthedocs.io/en/latest/generated/man/datalad-export.html>`__
   command with a plug-in interface and ``tarball`` plugin to export
   datasets
-  in Python, ``.api`` functions with rendering of results in command
   line got a \_-suffixed sibling, which would render results as well in
   Python as well (e.g., using ``search_`` instead of ``search`` would
   also render results, not only output them back as Python objects)
-  `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__

   -  ``--jobs`` option (passed to ``annex get``) for parallel downloads
   -  total and per-download (with git-annex >= 6.20160923) progress
      bars (note that if content is to be obtained from an archive, no
      progress will be reported yet)

-  `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   ``--reckless`` mode option
-  `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__

   -  highlights locations and fieldmaps for better readability
   -  supports ``-d^`` or ``-d///`` to point to top-most or centrally
      installed meta-datasets
   -  “complete” paths to the datasets are reported now
   -  ``-s`` option to specify which fields (only) to search

-  various enhancements and small fixes to
   `meta-data <http://docs.datalad.org/en/latest/cmdline.html#meta-data-handling>`__
   handling,
   `ls <http://datalad.readthedocs.io/en/latest/generated/man/datalad-ls.html>`__,
   custom remotes, code-base formatting, downloaders, etc
-  completely switched to ``tqdm`` library (``progressbar`` is no longer
   used/supported)

0.3 (Sep 23, 2016) – winter is coming
=====================================

Lots of everything, including but not limited to

-  enhanced index viewer, as the one on http://datasets.datalad.org
-  initial new data providers support:
   `Kaggle <https://www.kaggle.com>`__,
   `BALSA <http://balsa.wustl.edu>`__,
   `NDA <http://data-archive.nimh.nih.gov>`__,
   `NITRC <https://www.nitrc.org>`__
-  initial `meta-data support and
   management <http://docs.datalad.org/en/latest/cmdline.html#meta-data-handling>`__
-  new and/or improved crawler pipelines for
   `BALSA <http://balsa.wustl.edu>`__, `CRCNS <http://crcns.org>`__,
   `OpenfMRI <http://openfmri.org>`__
-  refactored
   `install <http://datalad.readthedocs.io/en/latest/generated/man/datalad-install.html>`__
   command, now with separate
   `get <http://datalad.readthedocs.io/en/latest/generated/man/datalad-get.html>`__
-  some other commands renaming/refactoring (e.g.,
   `create-sibling <http://datalad.readthedocs.io/en/latest/generated/man/datalad-create-sibling.html>`__)
-  datalad
   `search <http://datalad.readthedocs.io/en/latest/generated/man/datalad-search.html>`__
   would give you an option to install datalad’s super-dataset under
   ~/datalad if ran outside of a dataset

0.2.3 (Jun 28, 2016) – busy OHBM
--------------------------------

New features and bugfix release

-  support of /// urls to point to http://datasets.datalad.org
-  variety of fixes and enhancements throughout

0.2.2 (Jun 20, 2016) – OHBM we are coming!
------------------------------------------

New feature and bugfix release

-  greately improved documentation
-  publish command API RFing allows for custom options to annex, and
   uses –to REMOTE for consistent with annex invocation
-  variety of fixes and enhancements throughout

0.2.1 (Jun 10, 2016)
--------------------

-  variety of fixes and enhancements throughout

0.2 (May 20, 2016)
==================

Major RFing to switch from relying on rdf to git native submodules etc

0.1 (Oct 14, 2015)
==================

Release primarily focusing on interface functionality including initial
publishing
{% if fullname == 'datalad.api' -%}
`{{ name }}`
=={%- for c in name %}={%- endfor %}
.. automodule:: datalad.api

.. currentmodule:: datalad.api

{% for item in members if not item.startswith('_') %}
`{{ item }}`
--{%- for c in item %}-{%- endfor %}

.. autofunction:: {{ item }}
{% endfor %}

{% else -%}
{{ fullname }}
{{ underline }}

.. automodule:: {{ fullname }}
   :members:
   :undoc-members:
   :show-inheritance:
{% endif %}
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_python_imports:

************************
Python import statements
************************

.. topic:: Specification scope and status

   This specification describes the current (albeit incomplete) implementation.

The following rules apply to any ``import`` statement in the code base:

- All imports *must* be absolute, unless they import individual pieces of an integrated code component that is only split across several source code files for technical or organizational reasons.

- Imports *must* be placed at the top of a source file, unless there is a
  specific reason not to do so (e.g., delayed import due to performance
  concerns, circular dependencies). If such a reason exists, it *must*
  be documented by a comment at the import statement.

- There *must* be no more than one import per line.

- Multiple individual imports from a single module *must* follow the pattern::

      from <module> import (
          symbol1,
          symbol2,
      )

  Individual imported symbols *should* be sorted alphabetically. The last symbol
  line *should* end with a comma.

- Imports from packages and modules *should* be grouped in categories like

  - Standard library packages

  - 3rd-party packages

  - Datalad core (absolute imports)

  - Datalad extensions
  
  - Datalad core ("local" relative imports)
  
  Sorting imports can be aided by https://github.com/PyCQA/isort (e.g. ``python -m isort -m3 --fgw 2 --tc <filename>``).



Examples
========

::

    from collections import OrderedDict
    import logging
    import os

    from datalad.utils import (
        bytes2human,
        ensure_list,
        ensure_unicode,
        get_dataset_root as gdr,
    )
    
 In the `datalad/submodule/tests/test_mod.py` test file demonstrating an "exception" to absolute imports
 rule where test files are accompanying corresponding files of the underlying module:: 
 
    import os
  
    from datalad.utils import ensure_list
    
    from ..mod import func1

    from datalad.tests.utils import assert_true
    
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_log_levels:

**********
Log levels
**********

.. topic:: Specification scope and status

   This specification provides a partial overview of the current
   implementation.


Log messages are emitted by a wide range of operations within DataLad. They are
categorized into distinct levels. While some levels have self-explanatory
descriptions (e.g. ``warning``, ``error``), others are less specific (e.g.
``info``, ``debug``).

Common principles
=================

Parenthical log message use the same level
  When log messages are used to indicate the start and end of an operation,
  both start and end message use the same log-level.

Use cases
=========

Command execution
-----------------

For the :class:`~datalad.cmd.WitlessRunner` and its protocols the following log levels are used:

- High-level execution -> ``debug``
- Process start/finish -> ``8``
- Threading and IO -> ``5``
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_dataset_argument:

********************
``dataset`` argument
********************

.. topic:: Specification scope and status

   This specification describes the current implementation.

All commands which operate on datasets have a ``dataset`` argument (``-d`` or
``--dataset`` for the :term:`CLI`) to identify a single dataset as the
context of an operation.
If ``--dataset`` argument is not provided, the context of an operation is command-specific.
For example, `clone` command will consider the :term:`dataset` which is being cloned to be the context.
But typically, a dataset which current working directory belongs to is the context of an operation.
In the latter case, if operation (e.g., `get`) does not find a dataset in current working directory, operation fails with an ``NoDatasetFound`` error.


Impact on relative path resolution
==================================

With one exception, the nature of a provided ``dataset`` argument does **not**
impact the interpretation of relative paths. Relative paths are always considered
to be relative to the process working directory.

The one exception to this rule is passing a ``Dataset`` object instance as
``dataset`` argument value in the Python API. In this, and only this, case, a
relative path is interpreted as relative to the root of the respective dataset.


Special values
==============

There are some pre-defined "shortcut" values for dataset arguments:

``^``
   Represents to the topmost superdataset that contains the dataset the current
   directory is part of.
``^.``
   Represents the root directory of the dataset the current directory is part of.
``///``
   Represents the "default" dataset located under `$HOME/datalad/`.


Use cases
=========

Save modification in superdataset hierarchy
-------------------------------------------

Sometimes it is convenient to work only in the context of a subdataset.
Executing a ``datalad save <subdataset content>`` will record changes to the
subdataset, but will leave existing superdatasets dirty, as the subdataset
state change will not be saved there. Using the ``dataset`` argument it is
possible to redefine the scope of the save operation. For example::

  datalad save -d^ <subdataset content>

will perform the exact same save operation in the subdataset, but additionally
save all subdataset state changes in all superdatasets until the root of a
dataset hierarchy. Except for the specification of the dataset scope there is
no need to adjust path arguments or change the working directory.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_threaded_runner:


****************
Threaded runner
****************

.. topic:: Specification scope and status

   This specification provides an overview over the current implementation of the subprocess runner that is used throughout datalad.

Threads
=======

Datalad often requires the execution of subprocesses. While subprocesses are executed, datalad, i.e. its main thread, should be able to read data from stdout and stderr of the subprocess as well as write data to stdin of the subprocess. This requires a way to efficiently multiplex reading from stdout and stderr of the subprocess as well as writing to stdin of the subprocess.

Since non-blocking IO and waiting on multiple sources (poll or select) differs vastly in terms of capabilities and API on different OSs, we decided to use blocking IO and threads to multiplex reading from different sources.

Generally we have a number of threads that might be created and executed, depending on the need for writing to stdin or reading from stdout or stderr. Each thread can read from either a single queue or a file descriptor. Reading is done blocking. Each thread can put data into multiple queues. This is used to transport data that was read as well as for signaling conditions like closed file descriptors.

Conceptually, there are the main thread and two different types of threads:

 - type 1: transport threads (1 thread per process I/O descriptor)
 - type 2: process waiting thread (1 thread)

Transport Threads
.................

Besides the main thread, there might be up to three additional threads to handle data transfer to ``stdin``, and from ``stdout`` and ``stderr``. Each of those threads copies data between queues and file descriptors in a tight loop. The stdin-thread reads from an input-queue, the stdout- and stderr-threads write to an output queue. Each thread signals its exit to a set of signal queues, which might be identical to the output queues.

The ``stdin``-thread reads data from a queue and writes it to the ``stdin``-file descriptor of the sub-process. If it reads ``None`` from the queue, it will exit. The thread will also exit, if an exit is requested by calling ``thread.request_exit()``, or if an error occurs during writing. In all cases it will enqueue a ``None`` to all its signal-queues.

The ``stdout``- and ``stderr``-threads read from the respective file descriptor and enqueue data into their output queue, unless the data has zero length (which indicates a closed descriptor). On a zero-length read they exit and enqueue ``None`` into their signal queues.

All queues are infinite. Nevertheless signaling is performed with a timeout of one 100 milliseconds in order to ensure that threads can exit.


Process Waiting Thread
......................

The process waiting thread waits for a given process to exit and enqueues an exit notification into it signal queues.



Main Thread
...........

There is a single queue, the ``output_queue``, on which the main thread waits, after all transport threads, and the process waiting thread are started. The ``output_queue`` is the signaling queue and the output queue of the stderr-thread and the stdout-thread. It is also the signaling queue of the stdin-thread, and it is the signaling queue for the process waiting threads.

The main thread waits on the ``output_queue`` for data or signals and handles them accordingly, i.e. calls data callbacks of the protocol if data arrives, and calls connection-related callbacks of the protocol if other signals arrive. If no messages arrive on the  ``output_queue``, the main thread blocks for 100ms. If it is unblocked, either by getting a message or due to elapsing of the 100ms, it will process timeouts. If the ``timeout``-parameter to the constructor was not ``None``, it will check the last time any of the monitored files (stdout and/or stderr) yielded data. If the time is larger than the specified timeout, it will call the ``tiemout`` method of the protocol instance. Due to this implementation, the resolution for timeouts is 100ms. The main thread handles the closing of ``stdin``-, ``stdout``-, and ``stderr``-file descriptors if all other threads have terminated and if ``output_queue`` is empty. These tasks are either performed in the method ``ThreadedRunner.run()`` or in a result generator that is returned by  ``ThreadedRunner.run()`` whenever ``send()`` is called on it.


Protocols
=========

Due to its history the runner implementation uses the interface defined in ``SubprocessProtocol`` (asyncio.protocols.SubprocessProtocol) (although the sub process protocol interface is defined in the asyncio libraries, the current thread-runner implementation does not make use of ``async``).

    - ``SubprocessProtocol.pipe_data_received(fd, data)``
    - ``SubprocessProtocol.pipe_connection_lost(fd, exc)``
    - ``SubprocessProtocol.process_exited()``

In addition the methods of ``BaseProtocol`` are called, i.e.:

    - ``BaseProtocol.connection_made(transport)``
    - ``BaseProtocol.connection_lost(exc)``


The datalad-provided protocol ``WitlessProtocol`` provides an additional callback:

    - ``WitlessProtocol.timeout(fd)``

The method ``timeout()`` will be called when the parameter ``timeout`` in ``WitlessRunner.run``, ``ThreadedRunner.run``, or ``run_command`` is set to a number specifying the desired timeout in seconds. If no data is received from ``stdin``, or ``stderr`` (if those are supposed to be captured), the method ``WitlessProtocol.timeout(fd)`` is called with ``fd`` set to the respective file number, e.g. 1, or 2. If ``WitlessProtocol.timeout(fd)`` returns ``True``, the file descriptor will be closed and the associated threads will exit.

The method ``WitlessProtocol.timeout(fd)`` is also called if stdout, stderr and stdin are closed and the process does not exit within the given interval. In this case ``fd`` is set to ``None``. If ``WitlessProtocol.timeout(fd)`` returns ``True`` the process is terminated.


Object and Generator Results
================================

If the protocol that is provided to ``run()`` does not inherit ``datalad.runner.protocol.GeneratorMixIn``, the final result that will be returned to the caller is determined by calling ``WitlessProtocol._prepare_result()``. Whatever object this method returns will be returned to the caller.

If the protocol that is provided to ``run()`` does inherit ``datalad.runner.protocol.GeneratorMixIn``, ``run()`` will return a ``Generator``. This generator will yield the elements that were sent to it in the protocol-implementation by calling ``GeneratorMixIn.send_result()`` in the order in which the method ``GeneratorMixIn.send_result()`` is called. For example, if ``GeneratorMixIn.send_result(43)`` is called, the generator will yield ``43``, and if ``GeneratorMixIn.send_result({"a": 123, "b": "some data"})`` is called, the generator will yield ``{"a": 123, "b": "some data"}``.

Internally the generator is implemented by keeping track of the process state and waiting in the ``output_queue`` once, when ``send`` is called on it.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_result_records:

**************
Result records
**************

.. topic:: Specification scope and status

   This specification describes the current implementation.

Result records are the standard return value format for all DataLad commands.
Each command invocation yields one or more result records. Result records are
routinely inspected throughout the code base, and are used to inform generic
error handling, as well as particular calling commands on how to proceed with
a specific operation.

The technical implementation of a result record is a Python dictionary.  This
dictionary must contain a number of mandatory fields/keys (see below). However,
an arbitrary number of additional fields may be added to a result record.

The ``get_status_dict()`` function simplifies the creation of result records.

.. note::
   Developers *must* compose result records with care! DataLad supports custom
   user-provided hook configurations that use result record fields to
   decide when to trigger a custom post-result operation. Such custom hooks
   rely on a persistent naming and composition of result record fields.
   Changes to result records, including field name changes, field value changes,
   but also timing/order of record emitting potentially break user set ups!


Mandatory fields
================

The following keys *must* be present in any result record. If any of these
keys is missing, DataLad's behavior is undefined.


``action``
----------

A string label identifying which type of operation a result is associated with.
Labels *must not* contain white space. They should be compact, and lower-cases,
and use ``_`` (underscore) to separate words in compound labels.

A result without an ``action`` label will not be processed and is discarded.


``path``
--------

A string with an *absolute* path describing the local entity a result is
associated with. Paths must be platform-specific (e.g., Windows paths on
Windows, and POSIX paths on other operating systems). When a result is about an
entity that has no meaningful relation to the local file system (e.g., a URL to
be downloaded), to ``path`` value should be determined with respect to the
potential impact of the result on any local entity (e.g., a URL downloaded
to a local file path, a local dataset modified based on remote information).


``status``
----------

This field indicates the nature of a result in terms of four categories, identified
by a string label.

- ``ok``: a standard, to-be-expected result
- ``notneeded``: an operation that was requested, but found to be unnecessary
  in order to achieve a desired goal
- ``impossible``: a requested operation cannot be performed, possibly because
  its preconditions are not met
- ``error``: an error occurred while performing an operation

Based on the ``status`` field, a result is categorized into *success* (``ok``,
``notneeded``) and *failure* (``impossible``, ``error``). Depending on the
``on_failure`` parameterization of a command call, any failure-result emitted
by a command can lead to an ``IncompleteResultsError`` being raised on command
exit, or a non-zero exit code on the command line. With ``on_failure='stop'``,
an operation is halted on the first failure and the command errors out
immediately, with ``on_failure='continue'`` an operation will continue despite
intermediate failures and the command only errors out at the very end, with
``on_failure='ignore'`` the command will not error even when failures occurred.
The latter mode can be used in cases where the initial status-characterization
needs to be corrected for the particular context of an operation (e.g., to
relabel expected and recoverable errors).


Common optional fields
======================

The following fields are not required, but can be used to enrich a result
record with additional information that improves its interpretability, or
triggers particular optional functionality in generic result processing.


``type``
--------

This field indicates the type of entity a result is associated with. This may
or may not be the type of the local entity identified by the ``path`` value.
The following values are common, and should be used in matching cases, but
arbitrary other values are supported too:

- ``dataset``: a DataLad dataset
- ``file``: a regular file
- ``directory``: a directory
- ``symlink``: a symbolic link
- ``key``: a git-annex key
- ``sibling``: a Dataset sibling or Git remote


``message``
-----------

A message providing additional human-readable information on the nature or
provenance of a result. Any non-``ok`` results *should* have a message providing
information on the rational of their status characterization.

A message can be a string or a tuple. In case of a tuple, the second item can
contain values for ``%``-expansion of the message string. Expansion is performed
only immediately prior to actually outputting the message, hence string formatting
runtime costs can be avoided this way, if a message is not actually shown.


``logger``
----------

If a result record has a ``message`` field, then a given `Logger` instance
(typically from ``logging.getLogger()``) will be used to automatically log
this message. The log channel/level is determined based on
``datalad.log.result-level`` configuration setting. By default, this is
the ``debug`` level. When set to ``match-status`` the log level is determined
based on the ``status`` field of a result record:

- ``debug`` for ``'ok'``, and ``'notneeded'`` results
- ``warning`` for ``'impossible'`` results
- ``error`` for ``'error'`` results

This feature should be used with care. Unconditional logging can lead to
confusing double-reporting when results rendered and also visibly logged.


``refds``
---------

This field can identify a path (using the same semantics and requirements as
the ``path`` field) to a reference dataset that represents the larger context
of an operation. For example, when recursively processing multiple files across
a number of subdatasets, a ``refds`` value may point to the common superdataset.
This value may influence, for example, how paths are rendered in user-output.


``parentds``
------------

This field can identify a path (using the same semantics and requirements as
the ``path`` field) to a dataset containing an entity.


``state``
---------

A string label categorizing the state of an entity. Common values are:

- ``clean``
- ``untracked``
- ``modified``
- ``deleted``
- ``absent``
- ``present``


``error-messages``
------------------

List of any error messages that were captured or produced while achieving a
result.


``exception``
-------------

An exception that occurred while achieving the reported result.


``exception_traceback``
-----------------------

A string with a traceback for the exception reported in ``exception``.


Additional fields observed "in the wild"
========================================

Given that arbitrary fields are supported in result records, it is impossible
to compose a comprehensive list of field names (keys). However, in order to
counteract needless proliferation, the following list describes fields that
have been observed in implementations. Developers are encouraged to preferably
use compatible names from this list, or extend the list for additional items.

In alphabetical order:

``bytesize``
  The size of an entity in bytes (integer).

``gitshasum``
  SHA1 of an entity (string)

``prev_gitshasum``
  SHA1 of a previous state of an entity (string)

``key``
  The git-annex key associated with a ``type``-``file`` entity.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_file_url_handling:

*****************
File URL handling
*****************

.. topic:: Specification scope and status

   This specification describes the current implementation.

Datalad datasets can record URLs for file content access as metadata. This is a
feature provided by git-annex and is available for any annexed file. DataLad
improves upon the git-annex functionality in two ways:

1. Support for a variety of (additional) protocols and authentication methods.

2. Support for special URLs pointing to individual files located in registered
   (annexed) archives, such as tarballs and ZIP files.

These additional features are available to all functionality that is processing
URLs, such as ``get``, ``addurls``, or ``download-url``.


Extensible protocol and authentication support
==============================================

DataLad ships with a dedicated implementation of an external `git-annex special
remote`_ named ``git-annex-remote-datalad``. This is a somewhat atypical special
remote, because it cannot receive files and store them, but only supports
read operations.

Specifically, it uses the ``CLAIMURL`` feature of the `external special remote
protocol`_ to take over processing of URLs with supported protocols in all
datasets that have this special remote configured and enabled.

This special remote is automatically configured and enabled in DataLad dataset
as a ``datalad`` remote, by commands that utilize its features, such as
``download-url``. Once enabled, DataLad (but also git-annex) is able to act on
additional protocols, such as ``s3://``, and the respective URLs can be given
directly to commands like ``git annex addurl``, or ``datalad download-url``.

Beyond additional protocol support, the ``datalad`` special remote also
interfaces with DataLad's :ref:`chap_design_credentials`. It can identify a
particular credential required for a given URL (based on something called a
"provider" configuration), ask for the credential or retrieve it from a
credential store, and supply it to the respective service in an appropriate
form. Importantly, this feature neither requires the necessary credential or
provider configuration to be encoded in a URL (where it would become part of
the git-annex metadata), nor to be committed to a dataset. Hence all
information that may depend on which entity is performing a URL request
and in what environment is completely separated from the location information
on a particular file content. This minimizes the required dataset maintenance
effort (when credentials change), and offers a clean separation of identity
and availability tracking vs. authentication management.


Indexing and access of archive content
======================================

Another `git-annex special remote`_, named
``git-annex-remote-datalad-archives``, is used to enable file content retrieval
from annexed archive files, such as tarballs and ZIP files. Its implementation
concept is closely related to the ``git-annex-remote-datalad``, described
above.  Its main difference is that it claims responsibility for a particular
type of "URL" (starting with ``dl+archive:``). These URLs encode the identity
of an archive file, in terms of its git-annex key name, and a relative path
inside this archive pointing to a particular file.

Like ``git-annex-remote-datalad``, only read operations are supported. When
a request to a ``dl+archive:`` "URL" is made, the special remote identifies
the archive file, if necessary obtains it at the precise version needed, and
extracts the respected file content from the archive at the correct location.

This special remote is automatically configured and enabled as
``datalad-archives`` by the ``add-archive-content`` command. This command
indexes annexed archives, extracts, and registers their content to a
dataset.  File content availability information is recorded in terms of the
``dl+archive:`` "URLs", which are put into the git-annex metadata on a file's
content.


.. _git-annex special remote: https://git-annex.branchable.com/special_remotes/
.. _external special remote protocol: https://git-annex.branchable.com/design/external_special_remote_protocol
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_batched_command:

*******************************
BatchedCommand and BatchedAnnex
*******************************

.. topic:: Specification scope and status

   This specification describes the new implementation of ``BatchedCommand`` and
   ``BatchedAnnex`` in ``datalad``.


Batched Command
===============

The class ``BatchedCommand`` (in ``datalad.cmd``), holds an instance of a running subprocess, allows to send requests to the subprocess over its stdin, and to receive responses from the subprocess over its stdout.

Requests can be provided to an instance of ``BatchedCommand`` by passing a single request or a list of requests to ``BatchCommand.__call__()``, i.e. by applying the function call-operator to an instance of ``BatchedCommand``. A request is either a string or a tuple of strings. In the latter case, the elements of the tuple will be joined by ``" "``. More than one request can be given by providing a list of requests, i.e. a list of strings or tuples. In this case, the return value will be a list with one response for every request.

``BatchedCommand`` will send each request that is sent to the subprocess as a single line, after terminating the line by ``"\n"``. After the request is sent, ``BatchedCommand`` calls an output-handler with stdout-ish (an object that provides a ``readline()``-function which operates on the stdout of the subprocess) of the subprocess as argument. The output-handler can be provided to the constructor. If no output-handler is provided, a default output-handler is used. The default output-handler reads a single output line on stdout, using ``io.IOBase.readline()``, and returns the ``rstrip()``-ed line.

The subprocess must at least emit one line of output per line of input in order to prevent the calling thread from blocking. In addition, the size of the output, i.e. the number of lines that the result consists of, must be discernible by the output-handler. That means, the subprocess must either return a fixed number of lines per input line, or it must indicate the end of a result in some other way, e.g. with an empty line.

Remark: In principle any output processing could be performed. But, if the output-handler blocks on stdout, the calling thread will be blocked. Due to the limited capabilities of the stdout-ish that is passed to the output-handler, the output-handler must rely on ``readline()`` to process the output of the subprocess. Together with the line-based request sending, ``BatchedCommand`` is geared towards supporting the batch processing modes of ``git`` and ``git-annex``. *This has to be taken into account when providing a custom output handler.*

When ``BatchedCommand.close()`` is called, stdin, stdout, and stderr of the subprocess are closed. This indicates the end of processing to the subprocess. Generally the subprocess is expected to exit shortly after that. ``BatchedCommand.close()`` will wait for the subprocess to end, if the configuration ``datalad.runtime.stalled-external`` is set to ``"wait"``. If the configuration ``datalad.runtime.stalled-external`` is set to ``"abandon"``, ``BatchedCommand.close()`` will return after "timeout" seconds if ``timeout`` was provided to ``BatchedCommand.__init__()``, otherwise it will return after 11 seconds. If a timeout occurred, the attribute ``wait_timed_out`` of the ``BatchedCommand`` instance will be set to ``True``. If ``exception_on_timeout=True`` is provided to ``BatchedCommand.__init__()``, a ``subprocess.TimeoutExpired`` exception will be raised on a timeout while waiting for the process. It is not safe to reused a ``BatchedCommand`` instance after such an exception was risen.

Stderr of the subprocess is gathered in a byte-string. Its content will be returned by ``BatchCommand.close()`` if the parameter ``return_stderr`` is ``True``.


Implementation details
......................

``BatchedCommand`` uses ``WitlessRunner`` with a protocol that has ``datalad.runner.protocol.GeneratorMixIn`` as a super-class. The protocol uses an output-handler to process data, if an output-handler was specified during construction of ``BatchedCommand``.

``BatchedCommand.close()`` queries the configuration key ``datalad.runtime.stalled-external`` to determine how to handle non-exiting processes (there is no killing, processes or process zombies might just linger around until the next reboot).

The current implementation of ``BatchedCommand`` can process a list of multiple requests at once, but it will collect all answers before returning a result. That means, if you send 1000 requests, ``BatchedCommand`` will return after having received 1000 responses.


BatchedAnnex
============
``BatchedAnnex`` is a subclass of ``BatchedCommand`` (which it actually doesn't have to be, it just adds git-annex specific parameters to the command and sets a specific output handler).

``BatchedAnnex`` provides a new output-handler if the constructor-argument ``json`` is ``True``. In this case, an output handler is used that reads a single line from stdout, strips the line and converts it into a json object, which is returned. If the stripped line is empty, an empty dictionary is returned.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_application_vs_libary_mode:

***************************************
Application-type vs. library-type usage
***************************************

.. topic:: Specification scope and status

   This specification describes the current implementation.

Historically, DataLad was implemented with the assumption of application-type
usage, i.e., a person using DataLad through any of its APIs. Consequently,
(error) messaging was primarily targeting humans, and usage advice focused on
interactive use. With the increasing utilization of DataLad as an
infrastructural component it was necessary to address use cases of library-type
or internal usage more explicitly.

DataLad continues to behave like a stand-alone application by default.

For internal use, Python and command-line APIs provide dedicated mode switches.

Library mode can be enabled by setting the boolean configuration setting
``datalad.runtime.librarymode`` **before the start of the DataLad process**.
From the command line, this can be done with the option
``-c datalad.runtime.librarymode=yes``, or any other means for setting
configuration. In an already running Python process, library mode can be
enabled by calling ``datalad.enable_libarymode()``. This should be done
immediately after importing the ``datalad`` package for maximum impact.

.. code-block:: python

   >>> import datalad
   >>> datalad.enable_libarymode()

In a Python session, library mode **cannot** be enabled reliably by just setting
the configuration flag **after** the ``datalad`` package was already imported.
The ``enable_librarymode()`` function must be used.

Moreover, with ``datalad.in_librarymode()`` a query utility is provided that
can be used throughout the code base for adjusting behavior according to the
usage scenario.

Switching back and forth between modes during the runtime of a process is not
supported.

A library mode setting is exported into the environment of the Python process.
By default, it will be inherited by all child-processes, such as dataset
procedure executions.


Library-mode implications
=========================

No Python API docs
  Generation of comprehensive doc-strings for all API commands is skipped. This
  speeds up ``import datalad.api`` by about 30%.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_exception_handling:

******************
Exception handling
******************

.. topic:: Specification scope and status

   This specification describes the current implementation target.


Catching exceptions
===================

Whenever we catch an exception in an ``except`` clause, the following rules
apply:

- unless we (re-)raise, the first line instantiates a
  :class:`~datalad.support.exceptions.CapturedException`::

      except Exception as e:
          ce = CapturedException(e)

  First, this ensures a low-level (8) log entry including the traceback of that
  exception. The depth of the included traceback can be limited by setting the
  ``datalad.exc.str.tb_limit`` config accordingly.

  Second, it deletes the frame stack references of the exception and keeps
  textual information only, in order to avoid circular references, where an
  object (whose method raised the exception) isn't going to be picked by the
  garbage collection. This can be particularly troublesome if that object holds
  a reference to a subprocess for example. However, it's not easy to see in what
  situation this would really be needed and we never need anything other than
  the textual information about what happened. Making the reference cleaning a
  general rule is easiest to write, maintain and review.

- if we raise, neither a log entry nor such a
  :class:`~datalad.support.exceptions.CapturedException` instance is to be
  created.
  Eventually, there will be a spot where that (re-)raised exception is caught.
  This then is the right place to log it. That log entry will have the
  traceback, there's no need to leave a trace by means of log messages!

- if we raise, but do not simply reraise that exact same exception, in order to
  change the exception class and/or its message, ``raise from`` must be used!::

      except SomeError as e:
          raise NewError("new message") from e

  This ensures that the original exception is properly registered as the cause
  for the exception via its ``__cause__`` attribute. Hence, the original
  exception's traceback will be part of the later on logged traceback of the new
  exception.


Messaging about an exception
============================

In addition to the auto-generated low-level log entry there might be a need to
create a higher-level log, a user message or a (result) dictionary that includes
information from that exception. While such messaging may use anything the
(captured) exception provides, please consider that "technical" details about an
exception are already auto-logged and generally not incredibly meaningful for
users.

For message creation :class:`~datalad.support.exceptions.CapturedException`
comes with a couple of ``format_*`` helper methods, its ``__str__`` provides a
short representation of the form ``ExceptionClass(message)`` and its
``__repr__`` the log form with a traceback tht is used for the auto-generated
log.

For result dictionaries :class:`~datalad.support.exceptions.CapturedException`
can be assigned to the field ``exception``. Currently, ``get_status_dict`` will
consider this field and create an additional field with a traceback string.
Hence, whether putting a captured exception into that field actually has an
effect depends on whether ``get_status_dict`` is subsequently used with that
dictionary. In the future such functionality may move into result renderers
instead, leaving the decision of what to do with the passed
:class:`~datalad.support.exceptions.CapturedException` to them. Therefore, even
if of no immediate effect, enhancing the result dicts accordingly makes sense
already, since it may be useful when using datalad via its python interface
already and provide instant benefits whenever the result rendering gets such an
upgrade.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_docstrings:

**********
Docstrings
**********

.. topic:: Specification scope and status

   This specification provides a partial overview of the current
   implementation.

Docstrings in DataLad source code are used and consumed in many ways. Besides
serving as documentation directly in the sources, they are also transformed
and rendered in various ways.

- Command line ``--help`` output
- Python's ``help()`` or IPython's ``?``
- Manpages
- Sphinx-rendered documentation for the Python API and the command line API

A common source docstring is transformed, amended and tuned specifically for
each consumption scenario.


Formatting overview and guidelines
==================================

Version information
-------------------

Additions, changes, or deprecation should be recorded in a docstring using the
standard Sphinx directives ``versionadded``, ``versionchanged``,
``deprecated``::

  .. deprecated:: 0.16
     The ``dryrun||--dryrun`` option will be removed in a future release, use
     the renamed ``dry_run||--dry-run`` option instead.


API-conditional docs
--------------------

The ``CMD`` and ``PY`` macros can be used to selectively include documentation
for specific APIs only::

  options to pass to :command:`git init`. [PY: Options can be given as a list
  of command line arguments or as a GitPython-style option dictionary PY][CMD:
  Any argument specified after the destination path of the repository will be
  passed to git-init as-is CMD].

For API-alternative command and argument specifications the following format
can be used::

  ``<python-api>||<cmdline-api``

where the double backticks are mandatory and ``<python-part>`` and
``<cmdline-part>`` represent the respective argument specification for each
API. In these specifications only valid argument/command names are allowed,
plus a comma character to list multiples, and the dot character to include an
ellipsis::

   ``github_organization||-g,--github-organization``

   ``create_sibling_...||create-sibling-...``


Reflow text
-----------

When automatic transformations negatively affect the presentation of a
docstring due to excessive removal of content, leaving "holes", the ``REFLOW``
macro can be used to enclose such segments, in order to reformat them
as the final processing step. Example::

  || REFLOW >>
  The API has been aligned with the some
  ``create_sibling_...||create-sibling-...`` commands of other GitHub-like
  services, such as GOGS, GIN, GitTea.<< REFLOW ||

The start macro must appear on a dedicated line.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_credentials:

*********************
Credential management
*********************

.. topic:: Specification scope and status

   This specification describes the current implementation.

Various components of DataLad need to be passed credentials to interact with services that require authentication. 
This includes downloading files, but also things like REST API usage.

Supported credential types include basic user/password combinations, access tokens, and a range of tailored solutions for particular services.
All credential type implementations are derived from a common :class:`Credential` base class.

Importantly, credentials must be identified by a name.
This name is a label that is often hard-coded in the program code of DataLad, any of its extensions, or specified in a dataset.

Given a credential ``name``, one or more credential ``component``\(s) (e.g., ``token``, ``username``, or ``password``) can be looked up by DataLad in at least two different locations.
These locations are tried in the following order, and the first successful lookup yields the final value.

1. A configuration item ``datalad.credential.<name>.<component>``.
   Such configuration items can be defined in any location supported by DataLad's configuration system.
   As with any other specification of configuration items, environment variables can be used to set or override credentials.
   Variable names take the form of ``DATALAD_CREDENTIAL_<NAME>_<COMPONENT>``, and standard replacement rules into configuration variable names apply.

2. DataLad uses the `keyring` package https://pypi.org/project/keyring to connect to any of its supported back-ends for setting or getting credentials.
   This provides support for credential storage on all major platforms, but also extensibility, providing 3rd-parties to implement and use specialized solutions.

When a credential is required for operation, but could not be obtained via any of the above approaches, DataLad can prompt for credentials in interactive terminal sessions.
Interactively entered credentials will be stored in the active credential store available via the ``keyring`` package.

When a credential value is known but invalid, the invalid value must be removed or replaced in the active credential store.
By setting the configuration flag ``datalad.credentials.force-ask``, DataLad can be instructed to force interactive credential re-entry to effectively override any store credential with a new value.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_pos_vs_kw_parameters:

********************************
Positional vs Keyword parameters
********************************

.. topic:: Specification scope and status

   This specification is a proposal, subject to review and further discussion.
   Technical preview was implemented in the `PR #6176 <https://github.com/datalad/datalad/pull/6176>`_.

Motivation
==========

Python allows for keyword arguments (arguments with default values) to be specified positionally.
That complicates addition or removal of new keyword arguments since such changes must account for their possible
positional use.
Moreover, in case of our Interface's, it contributes to inhomogeneity since when used in :term:`CLI`, all keyword
arguments
must be specified via non-positional ``--<option>``'s, whenever Python interface allows for them to be used
positionally.

Python 3 added possibility to use a ``*`` separator in the function definition to mandate that all keyword arguments
*after* it must be be used only via keyword (``<option>=<value>``) specification.
It is encouraged to use ``*`` to explicitly separate out positional from keyword arguments in majority of the cases,
and below we outline two major types of constructs.

Interfaces
==========

Subclasses of the :class:`~datalad.interface.base.Interface` provide specification and implementation for both
:term:`CLI` and Python API interfaces.
All new interfaces must separate all CLI ``--options`` from positional arguments using ``*`` in their ``__call__``
signature.

**Note:** that some positional arguments could still be optional (e.g., destination ``path`` for ``clone``),
and thus should be listed **before** ``*``, despite been defined as a keyword argument in the ``__call__`` signature.

A unit-test will be provided to guarantee such consistency between :term:`CLI` and Python interfaces.
Overall, exceptions to this rule could be only some old(er) interfaces.

Regular functions and methods
=============================

Use of ``*`` is encouraged for any function (or method) with keyword arguments.
Generally, ``*`` should come before the first keyword argument, but similarly to the Interfaces above, it is left to
the discretion of the developer to possibly allocate some (just few) arguments which could be used positionally if
specified... -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_drop:

***********************
Drop dataset components
***********************

.. topic:: Specification scope and status

   This specification is a proposal, subject to review and further discussion.
   It is now partially implemented in the `drop` command.

§1 The :command:`drop` command is the antagonist of :command:`get`. Whatever a
`drop` can do, should be undoable by a subsequent :command:`get` (given
unchanged remote availability).

§2 Like :command:`get`, :command:`drop` primarily operates on a mandatory path
specification (to discover relevant files and sudatasets to operate on).

§3 :command:`drop` has ``--what`` parameter that serves as an extensible
"mode-switch" to cover all relevant scenarios, like 'drop all file content in
the work-tree' (e.g. ``--what files``, default, `#5858
<https://github.com/datalad/datalad/issues/5858>`__), 'drop all keys from any
branch' (i.e. ``--what allkeys``, `#2328
<https://github.com/datalad/datalad/issues/2328>`__), but also '"drop" AKA
uninstall entire subdataset hierarchies' (e.g. ``--what all``), or drop
preferred content (``--what preferred-content``, `#3122
<https://github.com/datalad/datalad/issues/3122>`__).

§4 :command:`drop` prevents data loss by default (`#4750
<https://github.com/datalad/datalad/issues/4750>`__). Like :command:`get` it
features a ``--reckless`` "mode-switch" to disable some or all potentially slow
safety mechnism, i.e. 'key available in sufficient number of other remotes',
'main or all branches pushed to remote(s)' (`#1142
<https://github.com/datalad/datalad/issues/1142>`__), 'only check availability
of keys associated with the worktree, but not other branches'. "Reckless
operation" can be automatic, when following a reckless :command:`get` (`#4744
<https://github.com/datalad/datalad/issues/4744>`__).

§5 :command:`drop` properly manages annex lifetime information, e.g. by announcing
an annex as ``dead`` on removal of a repository (`#3887
<https://github.com/datalad/datalad/issues/3887>`__).

§6 Like :command:`get`, drop supports parallization `#1953
<https://github.com/datalad/datalad/issues/1953>`__ 

§7 `datalad drop` is not intended to be a comprehensive frontend to `git annex
drop` (e.g. limited support for e.g. `#1482
<https://github.com/datalad/datalad/issues/1482>`__ outside standard use cases
like `#2328 <https://github.com/datalad/datalad/issues/2328>`__).

.. note::
  It is understood that the current `uninstall` command is largely or
  completely made obsolete by this :command:`drop` concept.

§8 Given the development in `#5842
<https://github.com/datalad/datalad/issues/5842>`__  towards the complete
obsolescence of `remove` it becomes necessary to import one of its proposed
features:

§9 :command:`drop` should be able to recognize a botched attempt to delete a
dataset with a plain rm -rf, and act on it in a meaningful way, even if it is
just hinting at chmod + rm -rf.


Use cases
=========

The following use cases operate in the dataset hierarchy depicted below::

  super
  ├── dir
  │   ├── fileD1
  │   └── fileD2
  ├── fileS1
  ├── fileS2
  ├── subA
  │   ├── fileA
  │   ├── subsubC
  │   │   ├── fileC
  │   └── subsubD
  └── subB
      └── fileB

Unless explicitly stated, all command are assumed to be executed in the root of `super`.

- U1: ``datalad drop fileS1``

   Drops the file content of `file1` (as currently done by :command:`drop`)

- U2: ``datalad drop dir``

   Drop all file content in the directory (``fileD{1,2}``; as currently done by
   :command:`drop`

- U3: ``datalad drop subB``

   Drop all file content from the entire `subB` (`fileB`)

- U4: ``datalad drop subB --what all``

   Same as above (default ``--what files``), because it is not operating in the
   context of a superdataset (no automatic upward lookups). Possibly hint at
   next usage pattern).

- U5: ``datalad drop -d . subB --what all``

  Drop all from the superdataset under this path. I.e. drop all from the
  subdataset and drop the subdataset itself (AKA uninstall)

- U6: ``datalad drop subA --what all``

  Error: "``subA`` contains subdatasets, forgot --recursive?"

- U7: ``datalad drop -d . subA -r --what all``

  Drop all content from the subdataset (``fileA``) and its subdatasets
  (``fileC``), uninstall the subdataset (``subA``) and its subdatasets
  (``subsubC``, ``subsubD``)

- U8: ``datalad drop subA -r --what all``

  Same as above, but keep ``subA`` installed

- U9: ``datalad drop sub-A -r``

   Drop all content from the subdataset and its subdatasets (``fileA``,
   ``fileC``)

- U10: ``datalad drop . -r --what all``

  Drops all file content and subdatasets, but leaves the superdataset
  repository behind

- U11: ``datalad drop -d . subB``

  Does nothing and hints at alternative usage, see
  https://github.com/datalad/datalad/issues/5832#issuecomment-889656335

- U12: ``cd .. && datalad drop super/dir``

  Like :command:`get`, errors because the execution is not associated with a
  dataset. This avoids complexities, when the given `path`'s point to multiple
  (disjoint) datasets. It is understood that it could be done, but it is
  intentionally not done. `datalad -C super drop dir` or `datalad drop -d super
  super/dir` would work.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_standard_parameters:

*******************
Standard parameters
*******************

.. topic:: Specification scope and status

   This specification partially describes the current implementation, and partially is a proposal, subject to review and further discussion.

Several "standard parameters" are used in various DataLad commands.
Those standard parameters have an identical meaning across the commands they are used in.
Commands should ensure that they use those "standard parameters" where applicable and do not deviate from the common names nor the common meaning.

Currently used standard parameters are listed below, as well as suggestions on how to harmonize currently deviating standard parameters.
Deviations from the agreed upon list should be harmonized.
The parameters are listed in their command-line form, but similar names and descriptions apply to their Python form.

``-d``/``--dataset``
  A pointer to the dataset that a given command should operate on

``--dry-run``
  Display details about the command execution without actually running the command.

``-f``/``--force``
  Enforce the execution of a command, even when certain security checks would normally prevent this

``-J``/``--jobs``
  Number of parallel jobs to use.

``-m``/``--message``
  A commit message to attach to the saved change of a command execution.

``-r``/``--recursive``
  Perform an operation recursively across subdatasets

``-R``/``--recursion-limit``
  Limit recursion to a given amount of subdataset levels

``-s``/``--sibling-name`` [SUGGESTION]
  The identifier for a dataset sibling (remote)


Certain standard parameters will have their own design document.
Please refer to those documents for more in-depth information... -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_designpatterns:

**********************
Miscellaneous patterns
**********************

DataLad is the result of a distributed and collaborative development effort
over many years.  During this time the scope of the project has changed
multiple times. As a consequence, the API and employed technologies have been
adjusted repeatedly.  Depending on the age of a piece of code, a clear software
design is not always immediately visible. This section documents a few design
patterns that the project strives to adopt at present. Changes to existing code
and new contributions should follow these guidelines.


Generator methods in `Repo` classes
===================================

Substantial parts of DataLad are implemented to behave like Python generators
in order to be maximally responsive when processing long-running tasks. This
included methods of the core API classes
:class:`~datalad.support.gitrepo.GitRepo` and
:class:`~datalad.support.annexrepo.AnnexRepo`. By convention, such methods
carry a trailing `_` in their name. In some cases, sibling methods with the
same name, but without the trailing underscore are provided. These behave like
their generator-equivalent, but eventually return an iterable once processing
is fully completed.


Calls to Git commands
=====================

DataLad is built on Git, so calls to Git commands are a key element of the code
base. All such calls should be made through methods of the
:class:`~datalad.support.gitrepo.GitRepo` class.  This is necessary, as only
there it is made sure that Git operates under the desired conditions
(environment configuration, etc.).

For some functionality, for example querying and manipulating `gitattributes`,
dedicated methods are provided. However, in many cases simple one-off calls to
get specific information from Git, or trigger certain operations are needed.
For these purposes the :class:`~datalad.support.gitrepo.GitRepo` class provides
a set of convenience methods aiming to cover use cases requiring particular
return values:

- test success of a command:
  :meth:`~datalad.support.gitrepo.GitRepo.call_git_success`
- obtain `stdout` of a command:
  :meth:`~datalad.support.gitrepo.GitRepo.call_git`
- obtain a single output line:
  :meth:`~datalad.support.gitrepo.GitRepo.call_git_oneline`
- obtain items from output split by a separator:
  :meth:`~datalad.support.gitrepo.GitRepo.call_git_items_`

All these methods take care of raising appropriate exceptions when expected
conditions are not met. Whenever desired functionality can be achieved
using simple custom calls to Git via these methods, their use is preferred
over the implementation of additional, dedicated wrapper methods.

Command examples
================

Examples of Python and commandline invocations of DataLad's user-oriented
commands are defined in the class of the respective command as dictionaries
within `_examples_`:

.. code-block:: python

   _examples_ = [
    dict(text="""Create a dataset 'mydataset' in the current directory""",
         code_py="create(path='mydataset')",
         code_cmd="datalad create mydataset",
    dict(text="""Apply the text2git procedure upon creation of a dataset""",
         code_py="create(path='mydataset', cfg_proc='text2git')",
         code_cmd="datalad create -c text2git mydataset")
         ]

The formatting of code lines is preserved. Changes to existing examples and
new contributions should provide examples for Python and commandline API, as
well as a concise description.
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design:

******
Design
******

The chapter described command API principles and the design of particular
subsystems in DataLad.

.. toctree::
   :maxdepth: 2

   application_vs_library_mode
   file_url_handling
   result_records
   dataset_argument
   log_levels
   drop
   python_imports
   miscpatterns
   exception_handling
   credentials
   url_substitution
   threaded_runner
   batched_command
   standard_parameters
   pos_vs_kw_parameters
   docstrings
.. -*- mode: rst -*-
.. vi: set ft=rst sts=4 ts=4 sw=4 et tw=79:

.. _chap_design_url_substitution:

****************
URL substitution
****************

.. topic:: Specification scope and status

   This specification describes the current implementation. This implementation
   is covering URL substitution in ``clone`` only. A further extension to
   URL processing elsewhere is possible.

URL substitution is a transformation of a given URL using a set of
specifications. Such specification can be provided as configuration settings
(via all supported configuration sources). These configuration items must
follow the naming scheme ``datalad.clone.url-substitute.<label>``, where
``<label>`` is an arbitrary identifier.

A substitution specification is a string with a match and substitution
expression, each following Python's regular expression syntax.  Both
expressions are concatenated into a single string with an arbitrary delimiter
character. The delimiter is defined by prefixing the string with the delimiter.
Prefix and delimiter are stripped from the expressions before processing.
Example::

  ,^http://(.*)$,https://\\1

A particular configuration item can be defined multiple times (see examples
below) to form a substitution series. Substitutions in the same series will be
applied incrementally, in order of their definition. If the first substitution
expression does not match, the entire series will be ignored. However,
following a first positive match all further substitutions in a series are
processed, regardless whether intermediate expressions match or not.

Any number of substitution series can be configured. They will be considered in
no particular order. Consequently, it advisable to implement the first match
specification of any series as specific as possible, in order to prevent
undesired transformations.


Examples
========

Change the protocol component of a given URL in order to hand over further
processing to a dedicated Git remote helper. Specifically, the following
example converts Open Science Framework project URLs like
``https://osf.io/f5j3e/`` into ``osf://f5j3e``, a URL that can be handle by
``git-remote-osf``, the Git remote helper provided by the `datalad-osf
extension package <https://github.com/datalad/datalad-osf>`__::

  datalad.clone.url-substitute.osf = ,^https://osf.io/([^/]+)[/]*$,osf://\1

Here is a more complex examples with a series of substitutions. The first
expression ensures that only GitHub URLs are being processed. The associated
substitution disassembles the URL into its two only relevant components,
the organisation/user name, and the project name::

  datalad.clone.url-substitute.github = ,https?://github.com/([^/]+)/(.*)$,\1###\2

All other expressions in this series that are described below will only be considered
if the above expression matched.

The next two expressions in the series normalize URL components that maybe be
auto-generated by some DataLad functionality, e.g. subdataset location
candidate generation from directory names::

  # replace (back)slashes with a single dash
  datalad.clone.url-substitute.github = ,[/\\]+,-

  # replace with whitespace (URL-quoted or not) with a single underscore
  datalad.clone.url-substitute.github = ,\s+|(%2520)+|(%20)+,_

The final expression in the series is recombining the organization/user name
and project name components back into a complete URL::

  datalad.clone.url-substitute.github = ,([^#]+)###(.*),https://github.com/\1/\2
