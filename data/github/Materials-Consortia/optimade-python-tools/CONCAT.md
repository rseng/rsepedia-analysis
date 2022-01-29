# Changelog

## [v0.16.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.9) (2022-01-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.8...v0.16.9)

**Implemented enhancements:**

- Lower validator default read timeout and allow it to be customised [\#1051](https://github.com/Materials-Consortia/optimade-python-tools/pull/1051) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Dependabot not updating NumPy to 1.22 [\#1035](https://github.com/Materials-Consortia/optimade-python-tools/issues/1035)

**Security fixes:**

- Elastic search log4j vulnerability [\#1040](https://github.com/Materials-Consortia/optimade-python-tools/issues/1040)

**Closed issues:**

- Remove multiple "Update dependencies" entries in CHANGELOG generation [\#1038](https://github.com/Materials-Consortia/optimade-python-tools/issues/1038)
- Docs reference to `LarkParser` failing. [\#1037](https://github.com/Materials-Consortia/optimade-python-tools/issues/1037)

**Merged pull requests:**

- Attempt to fix syntax for actions workflow [\#1053](https://github.com/Materials-Consortia/optimade-python-tools/pull/1053) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1049](https://github.com/Materials-Consortia/optimade-python-tools/pull/1049) ([CasperWA](https://github.com/CasperWA))
- Update dependabot config and changelog generation [\#1048](https://github.com/Materials-Consortia/optimade-python-tools/pull/1048) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1044](https://github.com/Materials-Consortia/optimade-python-tools/pull/1044) ([CasperWA](https://github.com/CasperWA))
- Bump elasticsearch image version to avoid any log4j issues [\#1041](https://github.com/Materials-Consortia/optimade-python-tools/pull/1041) ([ml-evs](https://github.com/ml-evs))
- Make NumPy requirement py version-specific [\#1036](https://github.com/Materials-Consortia/optimade-python-tools/pull/1036) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1031](https://github.com/Materials-Consortia/optimade-python-tools/pull/1031) ([CasperWA](https://github.com/CasperWA))

## [v0.16.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.8) (2021-12-22)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.7...v0.16.8)

**Implemented enhancements:**

- Support for Python 3.10 [\#956](https://github.com/Materials-Consortia/optimade-python-tools/issues/956)

**Fixed bugs:**

- Overzealous validation of substring comparisons for chemical formula fields [\#1024](https://github.com/Materials-Consortia/optimade-python-tools/issues/1024)

**Closed issues:**

- Updating to pymatgen v2022+ [\#762](https://github.com/Materials-Consortia/optimade-python-tools/issues/762)

**Merged pull requests:**

- Update dependencies [\#1028](https://github.com/Materials-Consortia/optimade-python-tools/pull/1028) ([CasperWA](https://github.com/CasperWA))
- Add configurable field-specific validator overrides to set filter operators as optional [\#1025](https://github.com/Materials-Consortia/optimade-python-tools/pull/1025) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1023](https://github.com/Materials-Consortia/optimade-python-tools/pull/1023) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1017](https://github.com/Materials-Consortia/optimade-python-tools/pull/1017) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1008](https://github.com/Materials-Consortia/optimade-python-tools/pull/1008) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1004](https://github.com/Materials-Consortia/optimade-python-tools/pull/1004) ([CasperWA](https://github.com/CasperWA))
- Add Python 3.10 support [\#957](https://github.com/Materials-Consortia/optimade-python-tools/pull/957) ([ml-evs](https://github.com/ml-evs))

## [v0.16.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.7) (2021-11-21)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.6...v0.16.7)

**Implemented enhancements:**

- Automate dependency workflow further [\#958](https://github.com/Materials-Consortia/optimade-python-tools/issues/958)
- Stricter validation of chemical formulas in OpenAPI schema [\#708](https://github.com/Materials-Consortia/optimade-python-tools/issues/708)

**Fixed bugs:**

- `chemical_formula_anonymous` validator accepts incorrect proportion order if started with 1 [\#1002](https://github.com/Materials-Consortia/optimade-python-tools/issues/1002)
- Reinstate `typing-extensions` [\#999](https://github.com/Materials-Consortia/optimade-python-tools/issues/999)
- Updating permanent dependabot branch not working after updating dependencies [\#995](https://github.com/Materials-Consortia/optimade-python-tools/issues/995)
- Auto-merge dependabot PR-workflow not running [\#984](https://github.com/Materials-Consortia/optimade-python-tools/issues/984)

**Closed issues:**

- Update the auto-PR description for updating deps [\#988](https://github.com/Materials-Consortia/optimade-python-tools/issues/988)
- Versioned docs do not redirect all links correctly [\#977](https://github.com/Materials-Consortia/optimade-python-tools/issues/977)
- Missing support for timestamps/datetime in grammar [\#102](https://github.com/Materials-Consortia/optimade-python-tools/issues/102)

**Merged pull requests:**

- Fixed bug in check\_anonymous\_formula which caused `chemical_formula_anonymous = AB2` to pass validation. [\#1001](https://github.com/Materials-Consortia/optimade-python-tools/pull/1001) ([JPBergsma](https://github.com/JPBergsma))
- Use `diff` for checking PR body [\#1000](https://github.com/Materials-Consortia/optimade-python-tools/pull/1000) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#998](https://github.com/Materials-Consortia/optimade-python-tools/pull/998) ([CasperWA](https://github.com/CasperWA))
- Correct PR body comparison [\#996](https://github.com/Materials-Consortia/optimade-python-tools/pull/996) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#993](https://github.com/Materials-Consortia/optimade-python-tools/pull/993) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#991](https://github.com/Materials-Consortia/optimade-python-tools/pull/991) ([CasperWA](https://github.com/CasperWA))
- Update dependency auto-PR message [\#989](https://github.com/Materials-Consortia/optimade-python-tools/pull/989) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#987](https://github.com/Materials-Consortia/optimade-python-tools/pull/987) ([CasperWA](https://github.com/CasperWA))
- Stricter formula syntax [\#986](https://github.com/Materials-Consortia/optimade-python-tools/pull/986) ([merkys](https://github.com/merkys))
- Implement workflows for dependency updates [\#979](https://github.com/Materials-Consortia/optimade-python-tools/pull/979) ([CasperWA](https://github.com/CasperWA))
- Tidy up old grammars, add a development grammar for v1.2 and update filterparser tests [\#879](https://github.com/Materials-Consortia/optimade-python-tools/pull/879) ([ml-evs](https://github.com/ml-evs))

## [v0.16.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.6) (2021-10-19)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.5...v0.16.6)

**Merged pull requests:**

- Put docs release deployment in separate job [\#978](https://github.com/Materials-Consortia/optimade-python-tools/pull/978) ([CasperWA](https://github.com/CasperWA))

## [v0.16.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.5) (2021-10-18)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.4...v0.16.5)

**Closed issues:**

- 'elements\_ratios' model validator uses double-precision machine epsilon - could be relaxed [\#947](https://github.com/Materials-Consortia/optimade-python-tools/issues/947)
- Versioning in Docs [\#724](https://github.com/Materials-Consortia/optimade-python-tools/issues/724)

**Merged pull requests:**

- Update dependencies [\#974](https://github.com/Materials-Consortia/optimade-python-tools/pull/974) ([CasperWA](https://github.com/CasperWA))
- Fix option value for checkout in CD Docs workflow [\#972](https://github.com/Materials-Consortia/optimade-python-tools/pull/972) ([CasperWA](https://github.com/CasperWA))
- Correct default branch name to `master` [\#971](https://github.com/Materials-Consortia/optimade-python-tools/pull/971) ([CasperWA](https://github.com/CasperWA))
- Dependabot updates for v0.16.5 [\#964](https://github.com/Materials-Consortia/optimade-python-tools/pull/964) ([ml-evs](https://github.com/ml-evs))
- Automate versioned documentation [\#951](https://github.com/Materials-Consortia/optimade-python-tools/pull/951) ([CasperWA](https://github.com/CasperWA))
- Add JOSS citation [\#949](https://github.com/Materials-Consortia/optimade-python-tools/pull/949) ([ml-evs](https://github.com/ml-evs))
- Some validation QoL tweaks [\#948](https://github.com/Materials-Consortia/optimade-python-tools/pull/948) ([ml-evs](https://github.com/ml-evs))

## [v0.16.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.4) (2021-09-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.3...v0.16.4)

**Closed issues:**

- Code check fails because there is no valid version of jsmin [\#938](https://github.com/Materials-Consortia/optimade-python-tools/issues/938)
- Be properly compliant with the new pip resolver [\#625](https://github.com/Materials-Consortia/optimade-python-tools/issues/625)

**Merged pull requests:**

- Bump providers from `357c27b` to `fb05359` [\#945](https://github.com/Materials-Consortia/optimade-python-tools/pull/945) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump providers from `368f9f6` to `357c27b` [\#944](https://github.com/Materials-Consortia/optimade-python-tools/pull/944) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump providers from `91b51bd` to `368f9f6` [\#942](https://github.com/Materials-Consortia/optimade-python-tools/pull/942) ([dependabot[bot]](https://github.com/apps/dependabot))
- remove the dependency on mkdocs-minify because of issue \#938. [\#941](https://github.com/Materials-Consortia/optimade-python-tools/pull/941) ([JPBergsma](https://github.com/JPBergsma))
- Corrected command to call uvicorn server [\#937](https://github.com/Materials-Consortia/optimade-python-tools/pull/937) ([JPBergsma](https://github.com/JPBergsma))
- Use proper pip dependency resolver in publish workflow [\#935](https://github.com/Materials-Consortia/optimade-python-tools/pull/935) ([ml-evs](https://github.com/ml-evs))
- Dependency updates for v0.16.4 [\#901](https://github.com/Materials-Consortia/optimade-python-tools/pull/901) ([ml-evs](https://github.com/ml-evs))
- Add JOSS paper [\#804](https://github.com/Materials-Consortia/optimade-python-tools/pull/804) ([ml-evs](https://github.com/ml-evs))

## [v0.16.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.3) (2021-09-02)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.2...v0.16.3)

**Implemented enhancements:**

- Add validation that anonymous/reduced chemical formulae are in fact reduced [\#913](https://github.com/Materials-Consortia/optimade-python-tools/issues/913)

**Fixed bugs:**

- No error/warning when specifying a config file that does not exist [\#930](https://github.com/Materials-Consortia/optimade-python-tools/issues/930)
- Docker tests failing in CI: http://gh\_actions\_host no longer exists? [\#906](https://github.com/Materials-Consortia/optimade-python-tools/issues/906)
- Fix config file warnings when file is missing [\#931](https://github.com/Materials-Consortia/optimade-python-tools/pull/931) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Docs don't introduce the idea of "models" [\#910](https://github.com/Materials-Consortia/optimade-python-tools/issues/910)
- Docs don't mention anything about where to go for support [\#909](https://github.com/Materials-Consortia/optimade-python-tools/issues/909)
- `run.sh` does not appear to be available from the pip installation [\#904](https://github.com/Materials-Consortia/optimade-python-tools/issues/904)
- Missing guide for how to set up an implementation from existing database [\#176](https://github.com/Materials-Consortia/optimade-python-tools/issues/176)

**Merged pull requests:**

- Add tutorial-style guide on setting up an API [\#915](https://github.com/Materials-Consortia/optimade-python-tools/pull/915) ([ml-evs](https://github.com/ml-evs))
- Add validator to check whether anonymous and reduced formulae are reduced [\#914](https://github.com/Materials-Consortia/optimade-python-tools/pull/914) ([ml-evs](https://github.com/ml-evs))
- Clarify the "all models" documentation page [\#912](https://github.com/Materials-Consortia/optimade-python-tools/pull/912) ([ml-evs](https://github.com/ml-evs))
- Add more specific 'Getting Help' info to Contributing and README [\#911](https://github.com/Materials-Consortia/optimade-python-tools/pull/911) ([ml-evs](https://github.com/ml-evs))
- Bump Materials-Consortia/optimade-validator-action from 2.5.0 to 2.6.0 [\#907](https://github.com/Materials-Consortia/optimade-python-tools/pull/907) ([dependabot[bot]](https://github.com/apps/dependabot))
- Clarify installation methods by use-case [\#905](https://github.com/Materials-Consortia/optimade-python-tools/pull/905) ([ml-evs](https://github.com/ml-evs))
- Relax response top-level root validator [\#903](https://github.com/Materials-Consortia/optimade-python-tools/pull/903) ([CasperWA](https://github.com/CasperWA))
- Add integrated app docs, tweak other use case docs [\#883](https://github.com/Materials-Consortia/optimade-python-tools/pull/883) ([ml-evs](https://github.com/ml-evs))

## [v0.16.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.2) (2021-08-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.1...v0.16.2)

**Fixed bugs:**

- Provider fallbacks are still not working [\#896](https://github.com/Materials-Consortia/optimade-python-tools/issues/896)
- Fix provider fallbacks [\#897](https://github.com/Materials-Consortia/optimade-python-tools/pull/897) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Dependency updates for v0.16.2 [\#894](https://github.com/Materials-Consortia/optimade-python-tools/pull/894) ([ml-evs](https://github.com/ml-evs))
- Bump codecov/codecov-action from 2.0.1 to 2.0.2 [\#882](https://github.com/Materials-Consortia/optimade-python-tools/pull/882) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump codecov/codecov-action from 1.5.2 to 2.0.1 [\#878](https://github.com/Materials-Consortia/optimade-python-tools/pull/878) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.16.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.1) (2021-07-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.0...v0.16.1)

**Implemented enhancements:**

- Change MIME type to application/vnd.api+json where appropriate [\#875](https://github.com/Materials-Consortia/optimade-python-tools/issues/875)
- Minor corrections + use model aliases for `handle_response_fields()` [\#876](https://github.com/Materials-Consortia/optimade-python-tools/pull/876) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Wrong behaviour HAS ONLY query for MongoDB [\#810](https://github.com/Materials-Consortia/optimade-python-tools/issues/810)
- Correct the behaviour of HAS ONLY with MongoDB backend [\#861](https://github.com/Materials-Consortia/optimade-python-tools/pull/861) ([JPBergsma](https://github.com/JPBergsma))

**Merged pull requests:**

- Change default MIME type to "application/vnd.api+json" [\#877](https://github.com/Materials-Consortia/optimade-python-tools/pull/877) ([ml-evs](https://github.com/ml-evs))
- Update elements description to match specification [\#874](https://github.com/Materials-Consortia/optimade-python-tools/pull/874) ([ml-evs](https://github.com/ml-evs))

## [v0.16.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.0) (2021-07-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.5...v0.16.0)

**Closed issues:**

- Incoming model update \(new field: issue\_tracker\) [\#592](https://github.com/Materials-Consortia/optimade-python-tools/issues/592)

**Merged pull requests:**

- Add issue\_tracker field to provider model [\#593](https://github.com/Materials-Consortia/optimade-python-tools/pull/593) ([ml-evs](https://github.com/ml-evs))

## [v0.15.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.5) (2021-07-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.4...v0.15.5)

**Fixed bugs:**

- NOT filter operation of mongo query for complex expressions [\#79](https://github.com/Materials-Consortia/optimade-python-tools/issues/79)

**Closed issues:**

- Remove CI psycopg2-binary install when aiida-core\>1.6.3 [\#855](https://github.com/Materials-Consortia/optimade-python-tools/issues/855)
- Pytest fails at Setup environment for AiiDA [\#853](https://github.com/Materials-Consortia/optimade-python-tools/issues/853)
- Add timeout parameter to validator [\#681](https://github.com/Materials-Consortia/optimade-python-tools/issues/681)
- Add note in installation instructions about pulling submodule for providers [\#370](https://github.com/Materials-Consortia/optimade-python-tools/issues/370)

**Merged pull requests:**

- Update dependencies [\#872](https://github.com/Materials-Consortia/optimade-python-tools/pull/872) ([ml-evs](https://github.com/ml-evs))
- Add request --timeout parameter to validator [\#860](https://github.com/Materials-Consortia/optimade-python-tools/pull/860) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `fa25ed3` to `91b51bd` [\#858](https://github.com/Materials-Consortia/optimade-python-tools/pull/858) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update to AiiDA v1.6.4 and remove CI fix [\#857](https://github.com/Materials-Consortia/optimade-python-tools/pull/857) ([CasperWA](https://github.com/CasperWA))
- Temporary fix for CI tests with AiiDA [\#854](https://github.com/Materials-Consortia/optimade-python-tools/pull/854) ([CasperWA](https://github.com/CasperWA))
- Documentation tweaks [\#852](https://github.com/Materials-Consortia/optimade-python-tools/pull/852) ([JPBergsma](https://github.com/JPBergsma))
- Fix query negation in MongoDB [\#814](https://github.com/Materials-Consortia/optimade-python-tools/pull/814) ([JPBergsma](https://github.com/JPBergsma))

## [v0.15.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.4) (2021-06-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.3...v0.15.4)

**Implemented enhancements:**

- Missing documentation for new configuration methods [\#766](https://github.com/Materials-Consortia/optimade-python-tools/issues/766)

**Closed issues:**

- Add docs "use case" for the validator [\#841](https://github.com/Materials-Consortia/optimade-python-tools/issues/841)
- Use specific configuration file for Heroku deployment [\#738](https://github.com/Materials-Consortia/optimade-python-tools/issues/738)
- Potential submission to JOSS? [\#203](https://github.com/Materials-Consortia/optimade-python-tools/issues/203)
- Add more tests [\#104](https://github.com/Materials-Consortia/optimade-python-tools/issues/104)

**Merged pull requests:**

- Tweak configuration docs [\#851](https://github.com/Materials-Consortia/optimade-python-tools/pull/851) ([ml-evs](https://github.com/ml-evs))
- Add some more tutorial-style documentation [\#850](https://github.com/Materials-Consortia/optimade-python-tools/pull/850) ([ml-evs](https://github.com/ml-evs))
- Bump FastAPI version in setup.py [\#849](https://github.com/Materials-Consortia/optimade-python-tools/pull/849) ([CasperWA](https://github.com/CasperWA))
- Bump fastapi from 0.65.1 to 0.65.2 [\#848](https://github.com/Materials-Consortia/optimade-python-tools/pull/848) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.15.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.3) (2021-06-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.2...v0.15.3)

**Merged pull requests:**

- Update model descriptions following spec updates [\#847](https://github.com/Materials-Consortia/optimade-python-tools/pull/847) ([ml-evs](https://github.com/ml-evs))

## [v0.15.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.2) (2021-06-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.1...v0.15.2)

**Implemented enhancements:**

- Missing HTTP response codes in OpenAPI schema [\#763](https://github.com/Materials-Consortia/optimade-python-tools/issues/763)

**Merged pull requests:**

- Update response model information for routes [\#846](https://github.com/Materials-Consortia/optimade-python-tools/pull/846) ([CasperWA](https://github.com/CasperWA))
- Improve semver validation error messsage [\#845](https://github.com/Materials-Consortia/optimade-python-tools/pull/845) ([ml-evs](https://github.com/ml-evs))
- Bump codecov/codecov-action from 1.5.0 to 1.5.2 [\#843](https://github.com/Materials-Consortia/optimade-python-tools/pull/843) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.15.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.1) (2021-06-08)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.0...v0.15.1)

**Closed issues:**

- mongomock $size queries match all non-array fields for {$size: 1}, even nulls [\#807](https://github.com/Materials-Consortia/optimade-python-tools/issues/807)
- Allow custom headers to be specified for validation [\#790](https://github.com/Materials-Consortia/optimade-python-tools/issues/790)

**Merged pull requests:**

- Allow both Jinja2 v2 and v3 [\#838](https://github.com/Materials-Consortia/optimade-python-tools/pull/838) ([CasperWA](https://github.com/CasperWA))
- Update mongomock and remove test skip [\#836](https://github.com/Materials-Consortia/optimade-python-tools/pull/836) ([ml-evs](https://github.com/ml-evs))
- Add --headers argument to validator to allow passing e.g. API keys [\#806](https://github.com/Materials-Consortia/optimade-python-tools/pull/806) ([ml-evs](https://github.com/ml-evs))

## [v0.15.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.0) (2021-06-01)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.14.1...v0.15.0)

**Fixed bugs:**

- Provider fallbacks do not get used [\#829](https://github.com/Materials-Consortia/optimade-python-tools/issues/829)
- ParserError's should not return 500 HTTP status codes [\#812](https://github.com/Materials-Consortia/optimade-python-tools/issues/812)
- Fix provider fallback list [\#830](https://github.com/Materials-Consortia/optimade-python-tools/pull/830) ([ml-evs](https://github.com/ml-evs))
- Return 400 Bad Request \(not 500\) on filter parser errors, plus filterparser module facelift [\#813](https://github.com/Materials-Consortia/optimade-python-tools/pull/813) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- CI triggered by PRs does not test dep versions in setup.py [\#834](https://github.com/Materials-Consortia/optimade-python-tools/issues/834)
- Remove Django support for v0.15+ [\#832](https://github.com/Materials-Consortia/optimade-python-tools/issues/832)
- Move aliasing code to base transformer [\#743](https://github.com/Materials-Consortia/optimade-python-tools/issues/743)
- Missing optional fields are not returned as null when requested with response\_fields [\#516](https://github.com/Materials-Consortia/optimade-python-tools/issues/516)

**Merged pull requests:**

- Test all setup.py deps versions for every pull request, plus some deps updates [\#835](https://github.com/Materials-Consortia/optimade-python-tools/pull/835) ([ml-evs](https://github.com/ml-evs))
- Deprecate Python 3.6, remove Django and update dependencies/providers [\#828](https://github.com/Materials-Consortia/optimade-python-tools/pull/828) ([ml-evs](https://github.com/ml-evs))
- Update INSTALL docs [\#811](https://github.com/Materials-Consortia/optimade-python-tools/pull/811) ([ml-evs](https://github.com/ml-evs))
- Overhaul of filter transformers, mappers and response fields [\#797](https://github.com/Materials-Consortia/optimade-python-tools/pull/797) ([ml-evs](https://github.com/ml-evs))

## [v0.14.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.14.1) (2021-05-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.14.0...v0.14.1)

**Fixed bugs:**

- \[SECURITY\] Cycle secrets [\#777](https://github.com/Materials-Consortia/optimade-python-tools/issues/777)

**Closed issues:**

- Do not validate extension endpoints [\#793](https://github.com/Materials-Consortia/optimade-python-tools/issues/793)
- Verify that missing values are not returned in comparisons [\#792](https://github.com/Materials-Consortia/optimade-python-tools/issues/792)

**Merged pull requests:**

- Bump pydantic from 1.8.1 to 1.8.2 [\#805](https://github.com/Materials-Consortia/optimade-python-tools/pull/805) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update GH actions [\#803](https://github.com/Materials-Consortia/optimade-python-tools/pull/803) ([CasperWA](https://github.com/CasperWA))
- Handling null fields in the filtertransformer and validator [\#796](https://github.com/Materials-Consortia/optimade-python-tools/pull/796) ([ml-evs](https://github.com/ml-evs))
- Filter out extension endpoints before validation [\#794](https://github.com/Materials-Consortia/optimade-python-tools/pull/794) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `7a54843` to `fa25ed3` [\#791](https://github.com/Materials-Consortia/optimade-python-tools/pull/791) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump typing-extensions from 3.7.4.3 to 3.10.0.0 [\#789](https://github.com/Materials-Consortia/optimade-python-tools/pull/789) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#787](https://github.com/Materials-Consortia/optimade-python-tools/pull/787) ([CasperWA](https://github.com/CasperWA))
- Bump CharMixer/auto-changelog-action from v1.2 to v1.3 [\#778](https://github.com/Materials-Consortia/optimade-python-tools/pull/778) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump django from 3.1.7 to 3.1.8 [\#776](https://github.com/Materials-Consortia/optimade-python-tools/pull/776) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#773](https://github.com/Materials-Consortia/optimade-python-tools/pull/773) ([ml-evs](https://github.com/ml-evs))

## [v0.14.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.14.0) (2021-03-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.3...v0.14.0)

**Implemented enhancements:**

- Rename config variable use\_real\_mongo to something more general [\#742](https://github.com/Materials-Consortia/optimade-python-tools/issues/742)
- Custom configuration extensions & use standard pydantic way of loading config file [\#739](https://github.com/Materials-Consortia/optimade-python-tools/issues/739)
- Generalising collections and adding ElasticsearchCollection [\#660](https://github.com/Materials-Consortia/optimade-python-tools/pull/660) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Over-aggressive middleware to check versioned base URL [\#737](https://github.com/Materials-Consortia/optimade-python-tools/issues/737)
- Floating point comparisons should not be tested with the validator [\#735](https://github.com/Materials-Consortia/optimade-python-tools/issues/735)
- Mapper method `alias_of` extracts alias wrongly [\#667](https://github.com/Materials-Consortia/optimade-python-tools/issues/667)

**Closed issues:**

- Docs builds are not properly tested for each PR [\#747](https://github.com/Materials-Consortia/optimade-python-tools/issues/747)
- Remove SQLAlchemy version fix in CI with new AiiDA version [\#745](https://github.com/Materials-Consortia/optimade-python-tools/issues/745)

**Merged pull requests:**

- Update dependencies [\#760](https://github.com/Materials-Consortia/optimade-python-tools/pull/760) ([CasperWA](https://github.com/CasperWA))
- Fix CheckWronglyVersionedBaseUrls middleware \(for landing pages\) [\#752](https://github.com/Materials-Consortia/optimade-python-tools/pull/752) ([CasperWA](https://github.com/CasperWA))
- Deprecate Python 3.6 support, v0.14 last supported version [\#751](https://github.com/Materials-Consortia/optimade-python-tools/pull/751) ([CasperWA](https://github.com/CasperWA))
- Run full API docs invoke task for every PR [\#748](https://github.com/Materials-Consortia/optimade-python-tools/pull/748) ([ml-evs](https://github.com/ml-evs))
- Change aliasing method names in mapper and deprecate the old [\#746](https://github.com/Materials-Consortia/optimade-python-tools/pull/746) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `e2074e8` to `7a54843` [\#741](https://github.com/Materials-Consortia/optimade-python-tools/pull/741) ([dependabot[bot]](https://github.com/apps/dependabot))
- Config updates [\#740](https://github.com/Materials-Consortia/optimade-python-tools/pull/740) ([CasperWA](https://github.com/CasperWA))
- Disable all floating-point comparisons during validation [\#736](https://github.com/Materials-Consortia/optimade-python-tools/pull/736) ([ml-evs](https://github.com/ml-evs))
- Report user errors in filter as HTTP 400 Bad Request and not 501 Not Implemented [\#658](https://github.com/Materials-Consortia/optimade-python-tools/pull/658) ([markus1978](https://github.com/markus1978))

## [v0.13.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.3) (2021-03-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.2...v0.13.3)

**Fixed bugs:**

- Support `anyOf`, `allOf`, etc. standard OpenAPI fields [\#730](https://github.com/Materials-Consortia/optimade-python-tools/issues/730)
- Python 3.9 support invalid [\#728](https://github.com/Materials-Consortia/optimade-python-tools/issues/728)

**Merged pull requests:**

- Update dependencies [\#734](https://github.com/Materials-Consortia/optimade-python-tools/pull/734) ([CasperWA](https://github.com/CasperWA))
- Update pydantic to ~=1.8 [\#731](https://github.com/Materials-Consortia/optimade-python-tools/pull/731) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `da74513` to `e2074e8` [\#727](https://github.com/Materials-Consortia/optimade-python-tools/pull/727) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.13.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.2) (2021-03-01)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.1...v0.13.2)

**Implemented enhancements:**

- Improve validation of providers [\#723](https://github.com/Materials-Consortia/optimade-python-tools/issues/723)

**Merged pull requests:**

- Update dependencies [\#725](https://github.com/Materials-Consortia/optimade-python-tools/pull/725) ([CasperWA](https://github.com/CasperWA))

## [v0.13.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.1) (2021-02-23)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.0...v0.13.1)

**Fixed bugs:**

- Supported OPTIMADE \_\_api\_version\_\_ is incorrect in latest release [\#712](https://github.com/Materials-Consortia/optimade-python-tools/issues/712)

**Merged pull requests:**

- Bump OPTIMADE version [\#713](https://github.com/Materials-Consortia/optimade-python-tools/pull/713) ([ml-evs](https://github.com/ml-evs))

## [v0.13.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.0) (2021-02-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.9...v0.13.0)

**Closed issues:**

- Update species.mass model [\#630](https://github.com/Materials-Consortia/optimade-python-tools/issues/630)

**Merged pull requests:**

- Update species-\>mass field following specification change [\#631](https://github.com/Materials-Consortia/optimade-python-tools/pull/631) ([ml-evs](https://github.com/ml-evs))

## [v0.12.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.9) (2021-02-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.8...v0.12.9)

**Implemented enhancements:**

- Improve support for timestamp queries in MongoTransformer [\#590](https://github.com/Materials-Consortia/optimade-python-tools/pull/590) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Use Enums for pydantic model defaults instead of strings [\#683](https://github.com/Materials-Consortia/optimade-python-tools/issues/683)

**Closed issues:**

- When using `--as-type` in validator, one does not get a summary \(`--json` doesn't work\) [\#699](https://github.com/Materials-Consortia/optimade-python-tools/issues/699)
- Extension/import issue with mongo collection [\#682](https://github.com/Materials-Consortia/optimade-python-tools/issues/682)

**Merged pull requests:**

- Update dependencies [\#707](https://github.com/Materials-Consortia/optimade-python-tools/pull/707) ([CasperWA](https://github.com/CasperWA))
- Always print summary as last thing in validation [\#700](https://github.com/Materials-Consortia/optimade-python-tools/pull/700) ([CasperWA](https://github.com/CasperWA))
- Bump django from 3.1.5 to 3.1.6 [\#698](https://github.com/Materials-Consortia/optimade-python-tools/pull/698) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#697](https://github.com/Materials-Consortia/optimade-python-tools/pull/697) ([CasperWA](https://github.com/CasperWA))
- Fixes for new gateway implementation [\#684](https://github.com/Materials-Consortia/optimade-python-tools/pull/684) ([CasperWA](https://github.com/CasperWA))

## [v0.12.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.8) (2021-01-18)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.7...v0.12.8)

**Implemented enhancements:**

- Validate mandatory query field `structure_features` [\#678](https://github.com/Materials-Consortia/optimade-python-tools/issues/678)

**Fixed bugs:**

- Validator should not rely on `meta->data_available` [\#677](https://github.com/Materials-Consortia/optimade-python-tools/issues/677)
- Validator should not rely on SHOULD "meta" field "data\_returned" [\#675](https://github.com/Materials-Consortia/optimade-python-tools/issues/675)
- Validator: remove reliance on meta fields and check mandatory queries [\#676](https://github.com/Materials-Consortia/optimade-python-tools/pull/676) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Bump providers from `542ac0a` to `da74513` [\#679](https://github.com/Materials-Consortia/optimade-python-tools/pull/679) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.12.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.7) (2021-01-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.6...v0.12.7)

**Implemented enhancements:**

- Make content-type response checks on '/versions` endpoint optional [\#670](https://github.com/Materials-Consortia/optimade-python-tools/pull/670) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Publish workflow fails when no changes to api docs between versions [\#673](https://github.com/Materials-Consortia/optimade-python-tools/issues/673)
- /versions header `Content-Type` value should be granularized according to RFC requirements in validator [\#669](https://github.com/Materials-Consortia/optimade-python-tools/issues/669)
- Misleading error message from validator on failure from '/versions' [\#668](https://github.com/Materials-Consortia/optimade-python-tools/issues/668)
- Fix publishing workflow [\#674](https://github.com/Materials-Consortia/optimade-python-tools/pull/674) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Update codecov coverage config file [\#672](https://github.com/Materials-Consortia/optimade-python-tools/pull/672) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `fe5048b` to `542ac0a` [\#671](https://github.com/Materials-Consortia/optimade-python-tools/pull/671) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.12.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.6) (2021-01-08)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.5...v0.12.6)

**Implemented enhancements:**

- Create base transformer [\#286](https://github.com/Materials-Consortia/optimade-python-tools/issues/286)

**Fixed bugs:**

- Our models and validator are too strict [\#399](https://github.com/Materials-Consortia/optimade-python-tools/issues/399)
- Validator changes: always check unversioned '/versions' and handle rich HTML pages [\#665](https://github.com/Materials-Consortia/optimade-python-tools/pull/665) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Add more prominent link to rendered docs [\#628](https://github.com/Materials-Consortia/optimade-python-tools/issues/628)
- Review the required properties of StructureResourceAttributes in openapi.json [\#198](https://github.com/Materials-Consortia/optimade-python-tools/issues/198)

**Merged pull requests:**

- Added GitHub CODEOWNERS [\#664](https://github.com/Materials-Consortia/optimade-python-tools/pull/664) ([ml-evs](https://github.com/ml-evs))
- Robustness improvements to validator [\#659](https://github.com/Materials-Consortia/optimade-python-tools/pull/659) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#655](https://github.com/Materials-Consortia/optimade-python-tools/pull/655) ([CasperWA](https://github.com/CasperWA))
- Bugfixes for elasticsearch filtertransformer comparision operators. [\#648](https://github.com/Materials-Consortia/optimade-python-tools/pull/648) ([markus1978](https://github.com/markus1978))
- Update dependencies [\#647](https://github.com/Materials-Consortia/optimade-python-tools/pull/647) ([ml-evs](https://github.com/ml-evs))
- Added "root\_path" config parameter for FastAPI apps [\#634](https://github.com/Materials-Consortia/optimade-python-tools/pull/634) ([markus1978](https://github.com/markus1978))
- Bump providers from `2673be6` to `fe5048b` [\#633](https://github.com/Materials-Consortia/optimade-python-tools/pull/633) ([dependabot[bot]](https://github.com/apps/dependabot))
- Updated README and moved some files to top-level [\#629](https://github.com/Materials-Consortia/optimade-python-tools/pull/629) ([ml-evs](https://github.com/ml-evs))
- insert reading of default optimade\_config.json in example run script run.sh [\#627](https://github.com/Materials-Consortia/optimade-python-tools/pull/627) ([rartino](https://github.com/rartino))
- Create template filtertransformer BaseTransformer [\#287](https://github.com/Materials-Consortia/optimade-python-tools/pull/287) ([ml-evs](https://github.com/ml-evs))

## [v0.12.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.5) (2020-12-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.4...v0.12.5)

**Closed issues:**

- PyPI publishing build is broken by latest pip [\#624](https://github.com/Materials-Consortia/optimade-python-tools/issues/624)
- Empty endpoints raise errors on validation [\#622](https://github.com/Materials-Consortia/optimade-python-tools/issues/622)
- Frequency of updating online docs [\#452](https://github.com/Materials-Consortia/optimade-python-tools/issues/452)

**Merged pull requests:**

- Fix PyPI publishing in CI [\#623](https://github.com/Materials-Consortia/optimade-python-tools/pull/623) ([ml-evs](https://github.com/ml-evs))
- Change validation error to warning on empty endpoints [\#621](https://github.com/Materials-Consortia/optimade-python-tools/pull/621) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#620](https://github.com/Materials-Consortia/optimade-python-tools/pull/620) ([CasperWA](https://github.com/CasperWA))
- Upstream fixes from specification [\#611](https://github.com/Materials-Consortia/optimade-python-tools/pull/611) ([ml-evs](https://github.com/ml-evs))
- Minor fixes for the validator [\#610](https://github.com/Materials-Consortia/optimade-python-tools/pull/610) ([ml-evs](https://github.com/ml-evs))
- Dependency updates [\#607](https://github.com/Materials-Consortia/optimade-python-tools/pull/607) ([ml-evs](https://github.com/ml-evs))
- include LICENSE in pip Package [\#594](https://github.com/Materials-Consortia/optimade-python-tools/pull/594) ([jan-janssen](https://github.com/jan-janssen))
- Relax models to allow for all SHOULD fields to be None [\#560](https://github.com/Materials-Consortia/optimade-python-tools/pull/560) ([ml-evs](https://github.com/ml-evs))
- Python 3.9 support [\#558](https://github.com/Materials-Consortia/optimade-python-tools/pull/558) ([ml-evs](https://github.com/ml-evs))
- ReadTheDocs configuration file \(v2\) [\#485](https://github.com/Materials-Consortia/optimade-python-tools/pull/485) ([CasperWA](https://github.com/CasperWA))

## [v0.12.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.4) (2020-11-16)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.3...v0.12.4)

**Merged pull requests:**

- Minor fixes for versions endpoint validation [\#591](https://github.com/Materials-Consortia/optimade-python-tools/pull/591) ([ml-evs](https://github.com/ml-evs))
- Add --minimal/--page\_limit validator options and remove old code [\#571](https://github.com/Materials-Consortia/optimade-python-tools/pull/571) ([ml-evs](https://github.com/ml-evs))

## [v0.12.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.3) (2020-11-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.2...v0.12.3)

**Fixed bugs:**

- GITHUB\_TOKEN not useful for changelog action [\#587](https://github.com/Materials-Consortia/optimade-python-tools/issues/587)
- Hill notation wrong \(still\) [\#585](https://github.com/Materials-Consortia/optimade-python-tools/issues/585)
- Hill notation validation turning around C and H [\#581](https://github.com/Materials-Consortia/optimade-python-tools/issues/581)

**Closed issues:**

- Make structure "deformity" tests more robust [\#583](https://github.com/Materials-Consortia/optimade-python-tools/issues/583)
- Incomplete output of optimade-validator [\#568](https://github.com/Materials-Consortia/optimade-python-tools/issues/568)

**Merged pull requests:**

- Use special release PAT for CHANGELOG generation action [\#588](https://github.com/Materials-Consortia/optimade-python-tools/pull/588) ([CasperWA](https://github.com/CasperWA))
- Check for carbon in elements for Hill [\#586](https://github.com/Materials-Consortia/optimade-python-tools/pull/586) ([CasperWA](https://github.com/CasperWA))
- Added better expected error messages to deformity tests [\#584](https://github.com/Materials-Consortia/optimade-python-tools/pull/584) ([ml-evs](https://github.com/ml-evs))
- Fix Hill ordering validation [\#582](https://github.com/Materials-Consortia/optimade-python-tools/pull/582) ([CasperWA](https://github.com/CasperWA))
- Bump mkdocs-material from 6.1.0 to 6.1.2 [\#580](https://github.com/Materials-Consortia/optimade-python-tools/pull/580) ([dependabot[bot]](https://github.com/apps/dependabot))
- Moved CONFIG import so it does not get triggered when just importing mapper [\#569](https://github.com/Materials-Consortia/optimade-python-tools/pull/569) ([ml-evs](https://github.com/ml-evs))

## [v0.12.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.2) (2020-10-31)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.1...v0.12.2)

**Implemented enhancements:**

- Add convenience method for adding all required middleware [\#536](https://github.com/Materials-Consortia/optimade-python-tools/issues/536)
- Add model validators and regexp for chemical formulae fields [\#547](https://github.com/Materials-Consortia/optimade-python-tools/pull/547) ([ml-evs](https://github.com/ml-evs))
- Validator improvements [\#515](https://github.com/Materials-Consortia/optimade-python-tools/pull/515) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- 'Chosen entry had no value for ...' when property is not requested [\#514](https://github.com/Materials-Consortia/optimade-python-tools/issues/514)
- Fix Species validators and error messages [\#561](https://github.com/Materials-Consortia/optimade-python-tools/pull/561) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Chemical symbols D and T [\#570](https://github.com/Materials-Consortia/optimade-python-tools/issues/570)
- Push back dependabot to monthly updates [\#567](https://github.com/Materials-Consortia/optimade-python-tools/issues/567)
- Spurious validation errors in Structure-\>Species [\#559](https://github.com/Materials-Consortia/optimade-python-tools/issues/559)
- Chemical formulae are not properly validated on model creation [\#546](https://github.com/Materials-Consortia/optimade-python-tools/issues/546)

**Merged pull requests:**

- Update dependencies [\#578](https://github.com/Materials-Consortia/optimade-python-tools/pull/578) ([CasperWA](https://github.com/CasperWA))
- Bump CasperWA/push-protected from v1 to v2.1.0 [\#573](https://github.com/Materials-Consortia/optimade-python-tools/pull/573) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update deps [\#566](https://github.com/Materials-Consortia/optimade-python-tools/pull/566) ([ml-evs](https://github.com/ml-evs))
- Improve handling of MongoDB ObjectID [\#557](https://github.com/Materials-Consortia/optimade-python-tools/pull/557) ([ml-evs](https://github.com/ml-evs))
- Update deps [\#556](https://github.com/Materials-Consortia/optimade-python-tools/pull/556) ([ml-evs](https://github.com/ml-evs))
- Updated dependencies [\#551](https://github.com/Materials-Consortia/optimade-python-tools/pull/551) ([ml-evs](https://github.com/ml-evs))
- Update dependencies - remove black as direct dependency [\#545](https://github.com/Materials-Consortia/optimade-python-tools/pull/545) ([CasperWA](https://github.com/CasperWA))
- Added convenience variables for middleware and exception handlers [\#537](https://github.com/Materials-Consortia/optimade-python-tools/pull/537) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#531](https://github.com/Materials-Consortia/optimade-python-tools/pull/531) ([ml-evs](https://github.com/ml-evs))

## [v0.12.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.1) (2020-09-24)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.0...v0.12.1)

**Implemented enhancements:**

- Move entry schemas to separate submodule [\#511](https://github.com/Materials-Consortia/optimade-python-tools/pull/511) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator should allow implementations to return "501 Not Implemented" for unsupported filters [\#518](https://github.com/Materials-Consortia/optimade-python-tools/issues/518)
- Landing page wrong URL  [\#371](https://github.com/Materials-Consortia/optimade-python-tools/issues/371)

**Merged pull requests:**

- This should ensure requirements\*.txt are tested [\#527](https://github.com/Materials-Consortia/optimade-python-tools/pull/527) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#526](https://github.com/Materials-Consortia/optimade-python-tools/pull/526) ([CasperWA](https://github.com/CasperWA))
- Fix landing page URL [\#519](https://github.com/Materials-Consortia/optimade-python-tools/pull/519) ([shyamd](https://github.com/shyamd))
- Update dependencies [\#510](https://github.com/Materials-Consortia/optimade-python-tools/pull/510) ([ml-evs](https://github.com/ml-evs))
- Fixing typo `validatated` -\> `validated` [\#506](https://github.com/Materials-Consortia/optimade-python-tools/pull/506) ([merkys](https://github.com/merkys))
- Make validator respond to KeyboardInterrupts [\#505](https://github.com/Materials-Consortia/optimade-python-tools/pull/505) ([ml-evs](https://github.com/ml-evs))
- Add support levels to validator config [\#503](https://github.com/Materials-Consortia/optimade-python-tools/pull/503) ([ml-evs](https://github.com/ml-evs))
- Enable JSON response from the validator [\#502](https://github.com/Materials-Consortia/optimade-python-tools/pull/502) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#501](https://github.com/Materials-Consortia/optimade-python-tools/pull/501) ([CasperWA](https://github.com/CasperWA))

## [v0.12.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.0) (2020-09-11)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.11.0...v0.12.0)

**Fixed bugs:**

- Missing field descriptions in schema for Species-\>name and Person-\>name [\#492](https://github.com/Materials-Consortia/optimade-python-tools/issues/492)
- "type" field not marked as required for derived entry resource models [\#479](https://github.com/Materials-Consortia/optimade-python-tools/issues/479)
- OpenAPI validations fails due to incorrect type of "dimension\_types" [\#478](https://github.com/Materials-Consortia/optimade-python-tools/issues/478)
- Have fallbacks for retrieving providers list [\#450](https://github.com/Materials-Consortia/optimade-python-tools/issues/450)
- Commit only when necessary [\#495](https://github.com/Materials-Consortia/optimade-python-tools/pull/495) ([CasperWA](https://github.com/CasperWA))
- Fix field optonality inconsistency in schema [\#482](https://github.com/Materials-Consortia/optimade-python-tools/pull/482) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator message for wrong version [\#493](https://github.com/Materials-Consortia/optimade-python-tools/issues/493)
- Validator should validate versions endpoint [\#491](https://github.com/Materials-Consortia/optimade-python-tools/issues/491)
- List of providers not included in `/links` endpoint for index meta-database [\#454](https://github.com/Materials-Consortia/optimade-python-tools/issues/454)
- Validate bad version URLs responding with 553 Version Not Supported [\#427](https://github.com/Materials-Consortia/optimade-python-tools/issues/427)
- Nonexistent property 'list' in validator tests [\#423](https://github.com/Materials-Consortia/optimade-python-tools/issues/423)
- Test `data_returned` [\#402](https://github.com/Materials-Consortia/optimade-python-tools/issues/402)
- AiiDA tests only run on Python 3.8 in CI [\#401](https://github.com/Materials-Consortia/optimade-python-tools/issues/401)
- Links under top-level 'links' may be objects [\#394](https://github.com/Materials-Consortia/optimade-python-tools/issues/394)
- Suggestion: use absolute imports in app code to allow re-use [\#298](https://github.com/Materials-Consortia/optimade-python-tools/issues/298)
- Update mongomock requirement when next released [\#207](https://github.com/Materials-Consortia/optimade-python-tools/issues/207)
- error when browsing OpenAPI docs [\#192](https://github.com/Materials-Consortia/optimade-python-tools/issues/192)

**Merged pull requests:**

- Don't report untracked and ignored files [\#496](https://github.com/Materials-Consortia/optimade-python-tools/pull/496) ([CasperWA](https://github.com/CasperWA))
- Improved error message for bad version returning 553 [\#494](https://github.com/Materials-Consortia/optimade-python-tools/pull/494) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#490](https://github.com/Materials-Consortia/optimade-python-tools/pull/490) ([CasperWA](https://github.com/CasperWA))
- Allow Link objects for pagination [\#484](https://github.com/Materials-Consortia/optimade-python-tools/pull/484) ([ml-evs](https://github.com/ml-evs))
- Absolute imports [\#483](https://github.com/Materials-Consortia/optimade-python-tools/pull/483) ([CasperWA](https://github.com/CasperWA))
- Validate OpenAPI specification in CI [\#481](https://github.com/Materials-Consortia/optimade-python-tools/pull/481) ([ml-evs](https://github.com/ml-evs))
- Update types to align with OpenAPI [\#480](https://github.com/Materials-Consortia/optimade-python-tools/pull/480) ([CasperWA](https://github.com/CasperWA))
- Update dependencies and pre-commit [\#477](https://github.com/Materials-Consortia/optimade-python-tools/pull/477) ([CasperWA](https://github.com/CasperWA))
- Unpin CI Python version for AiiDA tests [\#472](https://github.com/Materials-Consortia/optimade-python-tools/pull/472) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#471](https://github.com/Materials-Consortia/optimade-python-tools/pull/471) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#466](https://github.com/Materials-Consortia/optimade-python-tools/pull/466) ([CasperWA](https://github.com/CasperWA))
- Provider list fallback and list of providers in both servers' `/links`-endpoints [\#455](https://github.com/Materials-Consortia/optimade-python-tools/pull/455) ([CasperWA](https://github.com/CasperWA))
- SHOULD/MUST/OPTIONAL fields in models [\#453](https://github.com/Materials-Consortia/optimade-python-tools/pull/453) ([ml-evs](https://github.com/ml-evs))
- Validator overhaul [\#417](https://github.com/Materials-Consortia/optimade-python-tools/pull/417) ([ml-evs](https://github.com/ml-evs))

## [v0.11.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.11.0) (2020-08-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.10.0...v0.11.0)

**Implemented enhancements:**

- Use logging more thoroughly throughout the code base [\#242](https://github.com/Materials-Consortia/optimade-python-tools/issues/242)
- Implement `warnings` [\#105](https://github.com/Materials-Consortia/optimade-python-tools/issues/105)

**Fixed bugs:**

- Heroku is failing - raising OSError when making LOGS\_DIR [\#448](https://github.com/Materials-Consortia/optimade-python-tools/issues/448)
- `/versions` endpoint content-type parameter "header=present" is provided in the wrong place [\#418](https://github.com/Materials-Consortia/optimade-python-tools/issues/418)
- Publish workflow cannot push to protected branch [\#341](https://github.com/Materials-Consortia/optimade-python-tools/issues/341)
- Fix circular dep and extra permission error in logs [\#436](https://github.com/Materials-Consortia/optimade-python-tools/pull/436) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- log\_dir option in config is unused [\#435](https://github.com/Materials-Consortia/optimade-python-tools/issues/435)
- Allow all types of JSON API relationships [\#429](https://github.com/Materials-Consortia/optimade-python-tools/issues/429)
- OPTIMADE version badge was not bumped on 1.0 release [\#415](https://github.com/Materials-Consortia/optimade-python-tools/issues/415)
- Add `api_hint` query parameter [\#392](https://github.com/Materials-Consortia/optimade-python-tools/issues/392)
- Return 553 for wrongly versioned base URLs [\#391](https://github.com/Materials-Consortia/optimade-python-tools/issues/391)
- Private/dunder methods incorrectly documented in mkdocs [\#365](https://github.com/Materials-Consortia/optimade-python-tools/issues/365)
- Configuration documentation [\#310](https://github.com/Materials-Consortia/optimade-python-tools/issues/310)
- Improve handling of sorting in MongoDB backend [\#276](https://github.com/Materials-Consortia/optimade-python-tools/issues/276)

**Merged pull requests:**

- Catch OSError instead of PermissionError when making log dir [\#449](https://github.com/Materials-Consortia/optimade-python-tools/pull/449) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#447](https://github.com/Materials-Consortia/optimade-python-tools/pull/447) ([CasperWA](https://github.com/CasperWA))
- Bump mkdocstrings from 0.12.1 to 0.12.2 and mkdocs-material from 5.5.0 to 5.5.2 [\#440](https://github.com/Materials-Consortia/optimade-python-tools/pull/440) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump uvicorn from 0.11.5 to 0.11.7 [\#433](https://github.com/Materials-Consortia/optimade-python-tools/pull/433) ([dependabot[bot]](https://github.com/apps/dependabot))
- Introduce logging [\#432](https://github.com/Materials-Consortia/optimade-python-tools/pull/432) ([CasperWA](https://github.com/CasperWA))
- New middleware to catch any `OptimadeWarning`s [\#431](https://github.com/Materials-Consortia/optimade-python-tools/pull/431) ([CasperWA](https://github.com/CasperWA))
- Auto-generate API reference in docs and an overhaul [\#430](https://github.com/Materials-Consortia/optimade-python-tools/pull/430) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `52027b1` to `9712dd8` [\#428](https://github.com/Materials-Consortia/optimade-python-tools/pull/428) ([dependabot[bot]](https://github.com/apps/dependabot))
- Cleanup config files [\#426](https://github.com/Materials-Consortia/optimade-python-tools/pull/426) ([CasperWA](https://github.com/CasperWA))
- Update more unittest tests to pytest [\#425](https://github.com/Materials-Consortia/optimade-python-tools/pull/425) ([CasperWA](https://github.com/CasperWA))
- Sorting on unknown properties: returning Bad Request when appropriate [\#424](https://github.com/Materials-Consortia/optimade-python-tools/pull/424) ([ml-evs](https://github.com/ml-evs))
- Minor CI updates [\#422](https://github.com/Materials-Consortia/optimade-python-tools/pull/422) ([CasperWA](https://github.com/CasperWA))
- Add `api_hint` query parameter [\#421](https://github.com/Materials-Consortia/optimade-python-tools/pull/421) ([CasperWA](https://github.com/CasperWA))
- Implement 553 Version Not Supported [\#420](https://github.com/Materials-Consortia/optimade-python-tools/pull/420) ([CasperWA](https://github.com/CasperWA))
- Fix incorrect placement of header=present in versions endpoint [\#419](https://github.com/Materials-Consortia/optimade-python-tools/pull/419) ([ml-evs](https://github.com/ml-evs))
- Bump optimade-version.json to 1.0.0 [\#416](https://github.com/Materials-Consortia/optimade-python-tools/pull/416) ([ml-evs](https://github.com/ml-evs))
- Use optimade-validator-action v2 [\#413](https://github.com/Materials-Consortia/optimade-python-tools/pull/413) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `a96d424` to `52027b1` [\#389](https://github.com/Materials-Consortia/optimade-python-tools/pull/389) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.10.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.10.0) (2020-07-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.8...v0.10.0)

**Implemented enhancements:**

- Move tests to pytest system from unittest [\#270](https://github.com/Materials-Consortia/optimade-python-tools/issues/270)

**Fixed bugs:**

- Fix /vMAJOR/info in index server [\#414](https://github.com/Materials-Consortia/optimade-python-tools/pull/414) ([CasperWA](https://github.com/CasperWA))

**Closed issues:**

- Validation of 'structures' type crashes [\#397](https://github.com/Materials-Consortia/optimade-python-tools/issues/397)
- Validator verbosity levels need more detailed description [\#396](https://github.com/Materials-Consortia/optimade-python-tools/issues/396)
- Validator treats top-level 'included' array as mandatory [\#393](https://github.com/Materials-Consortia/optimade-python-tools/issues/393)
- \(Un\)versioned URLs [\#379](https://github.com/Materials-Consortia/optimade-python-tools/issues/379)

**Merged pull requests:**

- Update dependencies [\#412](https://github.com/Materials-Consortia/optimade-python-tools/pull/412) ([CasperWA](https://github.com/CasperWA))
- Bump pydantic from 1.5.1 to 1.6.1 [\#405](https://github.com/Materials-Consortia/optimade-python-tools/pull/405) ([dependabot[bot]](https://github.com/apps/dependabot))
- Temporarily run AiiDA tests on Python 3.8 only [\#400](https://github.com/Materials-Consortia/optimade-python-tools/pull/400) ([ml-evs](https://github.com/ml-evs))
- Make the example for --as\_type more similar to a real use case [\#398](https://github.com/Materials-Consortia/optimade-python-tools/pull/398) ([merkys](https://github.com/merkys))
- Fix some validator-specific crashes [\#395](https://github.com/Materials-Consortia/optimade-python-tools/pull/395) ([ml-evs](https://github.com/ml-evs))
- Use pytest instead of unittest [\#390](https://github.com/Materials-Consortia/optimade-python-tools/pull/390) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#388](https://github.com/Materials-Consortia/optimade-python-tools/pull/388) ([CasperWA](https://github.com/CasperWA))

## [v0.9.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.8) (2020-07-03)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.7...v0.9.8)

**Implemented enhancements:**

- Set implementation version in config by default [\#385](https://github.com/Materials-Consortia/optimade-python-tools/pull/385) ([CasperWA](https://github.com/CasperWA))

**Merged pull requests:**

- Update models, endpoints and responses to 1.0.0 [\#380](https://github.com/Materials-Consortia/optimade-python-tools/pull/380) ([ml-evs](https://github.com/ml-evs))

## [v0.9.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.7) (2020-06-28)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.6...v0.9.7)

## [v0.9.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.6) (2020-06-28)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.5...v0.9.6)

**Fixed bugs:**

- Fix publish workflow - final\(TM\) fix [\#378](https://github.com/Materials-Consortia/optimade-python-tools/pull/378) ([CasperWA](https://github.com/CasperWA))

## [v0.9.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.5) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.4...v0.9.5)

**Implemented enhancements:**

- Use new action for publishing [\#377](https://github.com/Materials-Consortia/optimade-python-tools/pull/377) ([CasperWA](https://github.com/CasperWA))

## [v0.9.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.4) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.3...v0.9.4)

## [v0.9.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.3) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.2...v0.9.3)

**Merged pull requests:**

- Fix version issues in the publish workflow [\#376](https://github.com/Materials-Consortia/optimade-python-tools/pull/376) ([shyamd](https://github.com/shyamd))
- Bump providers from `732593a` to `a96d424` [\#368](https://github.com/Materials-Consortia/optimade-python-tools/pull/368) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.9.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.2) (2020-06-25)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.1...v0.9.2)

**Fixed bugs:**

- Heroku cannot handle submodules when deploying via GitHub [\#373](https://github.com/Materials-Consortia/optimade-python-tools/issues/373)

**Closed issues:**

- Updates to models \(new OPTIONAL `type` field under `properties`\) [\#345](https://github.com/Materials-Consortia/optimade-python-tools/issues/345)
- Add aggregatation fields to links model [\#344](https://github.com/Materials-Consortia/optimade-python-tools/issues/344)
- Updates to models \(nperiodic\_dimensions\) [\#343](https://github.com/Materials-Consortia/optimade-python-tools/issues/343)
- Updates to models \(changing unknown atoms\) [\#342](https://github.com/Materials-Consortia/optimade-python-tools/issues/342)
- Improvements/fixes for openapi.json [\#332](https://github.com/Materials-Consortia/optimade-python-tools/issues/332)
- Update to v1.0.0-rc.1 [\#329](https://github.com/Materials-Consortia/optimade-python-tools/issues/329)
- Decouple updates in providers repo [\#311](https://github.com/Materials-Consortia/optimade-python-tools/issues/311)
- RST not rendering with mkdocs [\#307](https://github.com/Materials-Consortia/optimade-python-tools/issues/307)

**Merged pull requests:**

- Retrieve providers list if no submodule is found [\#374](https://github.com/Materials-Consortia/optimade-python-tools/pull/374) ([CasperWA](https://github.com/CasperWA))
- Update default implementation information [\#372](https://github.com/Materials-Consortia/optimade-python-tools/pull/372) ([shyamd](https://github.com/shyamd))
- Bump spec version to 1.0.0-rc.2 [\#367](https://github.com/Materials-Consortia/optimade-python-tools/pull/367) ([ml-evs](https://github.com/ml-evs))
- Dependabot updates: numpy, mkdocs-material, mkdocstrings, requests [\#364](https://github.com/Materials-Consortia/optimade-python-tools/pull/364) ([ml-evs](https://github.com/ml-evs))
- Merge all Dependabot updates [\#353](https://github.com/Materials-Consortia/optimade-python-tools/pull/353) ([shyamd](https://github.com/shyamd))
- Update model descriptions and openapi.json for 1.0.0-rc2 [\#351](https://github.com/Materials-Consortia/optimade-python-tools/pull/351) ([ml-evs](https://github.com/ml-evs))
- Update models according to changes during CECAM 2020 meeting [\#350](https://github.com/Materials-Consortia/optimade-python-tools/pull/350) ([ml-evs](https://github.com/ml-evs))
- Decouple changes in providers repo [\#312](https://github.com/Materials-Consortia/optimade-python-tools/pull/312) ([shyamd](https://github.com/shyamd))

## [v0.9.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.1) (2020-06-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.0...v0.9.1)

## [v0.9.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.0) (2020-06-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.8.1...v0.9.0)

**Implemented enhancements:**

- Breaking up the python tools into seperable packages [\#255](https://github.com/Materials-Consortia/optimade-python-tools/issues/255)
- Run both servers as standard [\#238](https://github.com/Materials-Consortia/optimade-python-tools/issues/238)

**Fixed bugs:**

- Non-running CI job [\#331](https://github.com/Materials-Consortia/optimade-python-tools/issues/331)
- Special species "X" not tested for non-disordered structures [\#304](https://github.com/Materials-Consortia/optimade-python-tools/issues/304)
- Standardize timezone of datetime responses [\#288](https://github.com/Materials-Consortia/optimade-python-tools/issues/288)
- Queries on aliased/provider fields are broken for nested properties [\#282](https://github.com/Materials-Consortia/optimade-python-tools/issues/282)
- General exceptions not being put into response [\#281](https://github.com/Materials-Consortia/optimade-python-tools/issues/281)
- Issue with CIF export [\#271](https://github.com/Materials-Consortia/optimade-python-tools/issues/271)
- Type-cast inputs for general Error [\#280](https://github.com/Materials-Consortia/optimade-python-tools/pull/280) ([CasperWA](https://github.com/CasperWA))

**Security fixes:**

- \[Security\] Bump django from 3.0.4 to 3.0.7 in /.github/workflows [\#291](https://github.com/Materials-Consortia/optimade-python-tools/pull/291) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))

**Closed issues:**

- Update links resources [\#299](https://github.com/Materials-Consortia/optimade-python-tools/issues/299)
- Need to set up mkdocs [\#289](https://github.com/Materials-Consortia/optimade-python-tools/issues/289)
- Need to add custom schema entries for unit/sortable \(and eventually type\) [\#278](https://github.com/Materials-Consortia/optimade-python-tools/issues/278)
- /info/\<entry-endpoint\> missing `sortable` key under each property [\#273](https://github.com/Materials-Consortia/optimade-python-tools/issues/273)
- Make CI linting more useful [\#269](https://github.com/Materials-Consortia/optimade-python-tools/issues/269)
- \[PR SPECIFIC\] Reminder: Validator test pinned to specific commit [\#268](https://github.com/Materials-Consortia/optimade-python-tools/issues/268)
- Validator does not check that pagination links work [\#265](https://github.com/Materials-Consortia/optimade-python-tools/issues/265)
- available\_api\_versions is not correctly validated [\#261](https://github.com/Materials-Consortia/optimade-python-tools/issues/261)
- Implementation model should allow for any URL type in `source_url` [\#260](https://github.com/Materials-Consortia/optimade-python-tools/issues/260)
- Extra structure endpoints in the api specification @ odbx [\#259](https://github.com/Materials-Consortia/optimade-python-tools/issues/259)
- Wrong response structure at info endpoint @ cod [\#258](https://github.com/Materials-Consortia/optimade-python-tools/issues/258)
- Missing base url for api's docs @ materialscloud [\#257](https://github.com/Materials-Consortia/optimade-python-tools/issues/257)
- Handling of KNOWN in mongo backend [\#254](https://github.com/Materials-Consortia/optimade-python-tools/issues/254)
- `None` values in `lattice_vectors` [\#170](https://github.com/Materials-Consortia/optimade-python-tools/issues/170)
- Make sure that the PyPI distribution works [\#143](https://github.com/Materials-Consortia/optimade-python-tools/issues/143)
- Move run.sh to a python file to be environment-agnostic [\#81](https://github.com/Materials-Consortia/optimade-python-tools/issues/81)

**Merged pull requests:**

- Another fix for release pipeline [\#355](https://github.com/Materials-Consortia/optimade-python-tools/pull/355) ([shyamd](https://github.com/shyamd))
- Fix publish workflow [\#354](https://github.com/Materials-Consortia/optimade-python-tools/pull/354) ([CasperWA](https://github.com/CasperWA))
- Fix publish workflow [\#352](https://github.com/Materials-Consortia/optimade-python-tools/pull/352) ([CasperWA](https://github.com/CasperWA))
- Update publish workflow [\#340](https://github.com/Materials-Consortia/optimade-python-tools/pull/340) ([shyamd](https://github.com/shyamd))
- Remove test publish action [\#338](https://github.com/Materials-Consortia/optimade-python-tools/pull/338) ([shyamd](https://github.com/shyamd))
- Fix 'publish\_TestPyPI' CI job [\#337](https://github.com/Materials-Consortia/optimade-python-tools/pull/337) ([CasperWA](https://github.com/CasperWA))
- Specify versions for all setup.py deps [\#336](https://github.com/Materials-Consortia/optimade-python-tools/pull/336) ([CasperWA](https://github.com/CasperWA))
- Represent the datetime objects as UTC in RFC3339 format [\#333](https://github.com/Materials-Consortia/optimade-python-tools/pull/333) ([fekad](https://github.com/fekad))
- dependamat: Bump \<package\_name\> v x.y.z to vx.y.\(z+1\) [\#330](https://github.com/Materials-Consortia/optimade-python-tools/pull/330) ([ml-evs](https://github.com/ml-evs))
- Bump fastapi from 0.53.1 to 0.56.0 [\#324](https://github.com/Materials-Consortia/optimade-python-tools/pull/324) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump pydantic from 1.4 to 1.5.1 [\#320](https://github.com/Materials-Consortia/optimade-python-tools/pull/320) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update links resources [\#306](https://github.com/Materials-Consortia/optimade-python-tools/pull/306) ([CasperWA](https://github.com/CasperWA))
- Add special species for adapters testing [\#305](https://github.com/Materials-Consortia/optimade-python-tools/pull/305) ([CasperWA](https://github.com/CasperWA))
- Clean Up Build Environment [\#301](https://github.com/Materials-Consortia/optimade-python-tools/pull/301) ([shyamd](https://github.com/shyamd))
- Enable CI failures for linting [\#300](https://github.com/Materials-Consortia/optimade-python-tools/pull/300) ([ml-evs](https://github.com/ml-evs))
- Adding jarvis-tools structures [\#297](https://github.com/Materials-Consortia/optimade-python-tools/pull/297) ([knc6](https://github.com/knc6))
- Update Docs [\#295](https://github.com/Materials-Consortia/optimade-python-tools/pull/295) ([shyamd](https://github.com/shyamd))
- Setup MKDocs for Documentation [\#294](https://github.com/Materials-Consortia/optimade-python-tools/pull/294) ([shyamd](https://github.com/shyamd))
- Fix filters on nested provider/aliased fields [\#285](https://github.com/Materials-Consortia/optimade-python-tools/pull/285) ([ml-evs](https://github.com/ml-evs))
- Use heroku-shields instead of heroku-badge [\#284](https://github.com/Materials-Consortia/optimade-python-tools/pull/284) ([CasperWA](https://github.com/CasperWA))
- Add OPTIMADE logo to badge by extending JSON [\#283](https://github.com/Materials-Consortia/optimade-python-tools/pull/283) ([CasperWA](https://github.com/CasperWA))
- Add null check to mongo filtertransformer for KNOWN/UNKNOWN filters [\#279](https://github.com/Materials-Consortia/optimade-python-tools/pull/279) ([ml-evs](https://github.com/ml-evs))
- Add `sortable=True` to all properties [\#274](https://github.com/Materials-Consortia/optimade-python-tools/pull/274) ([CasperWA](https://github.com/CasperWA))
- Make \_atom\_site\_label unique in CIF generation [\#272](https://github.com/Materials-Consortia/optimade-python-tools/pull/272) ([CasperWA](https://github.com/CasperWA))
- Not so quick fix to allow "/" at end of validator URL, plus fixes and tests for --as\_type [\#267](https://github.com/Materials-Consortia/optimade-python-tools/pull/267) ([ml-evs](https://github.com/ml-evs))
- Check pagination links-\>next with validator [\#266](https://github.com/Materials-Consortia/optimade-python-tools/pull/266) ([ml-evs](https://github.com/ml-evs))
- Relax HTTP URL constraints on meta-\>implementation-\>source\_url field. [\#262](https://github.com/Materials-Consortia/optimade-python-tools/pull/262) ([ml-evs](https://github.com/ml-evs))
- Validate lattice\_vectors for all null or all float [\#171](https://github.com/Materials-Consortia/optimade-python-tools/pull/171) ([CasperWA](https://github.com/CasperWA))

## [v0.8.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.8.1) (2020-04-25)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.8.0...v0.8.1)

**Fixed bugs:**

- Pip install missing some files [\#252](https://github.com/Materials-Consortia/optimade-python-tools/issues/252)

**Merged pull requests:**

- v0.8.1 hotfix [\#256](https://github.com/Materials-Consortia/optimade-python-tools/pull/256) ([ml-evs](https://github.com/ml-evs))
- Fix 252 missing landing page [\#253](https://github.com/Materials-Consortia/optimade-python-tools/pull/253) ([shyamd](https://github.com/shyamd))

## [v0.8.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.8.0) (2020-04-22)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.7.1...v0.8.0)

**Implemented enhancements:**

- Switch to pydantic's BaseSettings for the config file? [\#152](https://github.com/Materials-Consortia/optimade-python-tools/issues/152)
- Use services for testing/updating dependencies? [\#96](https://github.com/Materials-Consortia/optimade-python-tools/issues/96)
- Remove query constraints for /links-endpoint [\#244](https://github.com/Materials-Consortia/optimade-python-tools/pull/244) ([CasperWA](https://github.com/CasperWA))
- Add adapters - Base design + 'structures' \(+ 'references'... sort of\) [\#241](https://github.com/Materials-Consortia/optimade-python-tools/pull/241) ([CasperWA](https://github.com/CasperWA))
- Add dependabot and last commit date badges [\#237](https://github.com/Materials-Consortia/optimade-python-tools/pull/237) ([CasperWA](https://github.com/CasperWA))
- Add mongo length operator functionality with length aliases [\#222](https://github.com/Materials-Consortia/optimade-python-tools/pull/222) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Use Path.home\(\) instead of ~ in default config path values [\#245](https://github.com/Materials-Consortia/optimade-python-tools/issues/245)

**Closed issues:**

- Have Dependabot take care of various requirements.txt files as well [\#249](https://github.com/Materials-Consortia/optimade-python-tools/issues/249)
- Remove commented out GH Action job `deps_clean-install` [\#247](https://github.com/Materials-Consortia/optimade-python-tools/issues/247)
- Local testing fails without default config [\#239](https://github.com/Materials-Consortia/optimade-python-tools/issues/239)
- Release only when pushing to master [\#229](https://github.com/Materials-Consortia/optimade-python-tools/issues/229)
- Do we need `server.cfg`? [\#134](https://github.com/Materials-Consortia/optimade-python-tools/issues/134)
- Implement LENGTH in query [\#86](https://github.com/Materials-Consortia/optimade-python-tools/issues/86)

**Merged pull requests:**

- Up to v0.8.0 [\#251](https://github.com/Materials-Consortia/optimade-python-tools/pull/251) ([CasperWA](https://github.com/CasperWA))
- Remove old commented GH Action job [\#250](https://github.com/Materials-Consortia/optimade-python-tools/pull/250) ([CasperWA](https://github.com/CasperWA))
- Use Path.home\(\) instead of `~` [\#246](https://github.com/Materials-Consortia/optimade-python-tools/pull/246) ([CasperWA](https://github.com/CasperWA))
- Fix path in default config [\#243](https://github.com/Materials-Consortia/optimade-python-tools/pull/243) ([ml-evs](https://github.com/ml-evs))
- Fixes Local Tests [\#240](https://github.com/Materials-Consortia/optimade-python-tools/pull/240) ([shyamd](https://github.com/shyamd))
- Revert "Fix github actions for non-release tags" [\#236](https://github.com/Materials-Consortia/optimade-python-tools/pull/236) ([shyamd](https://github.com/shyamd))
- Enable filtering on relationships with mongo  [\#234](https://github.com/Materials-Consortia/optimade-python-tools/pull/234) ([ml-evs](https://github.com/ml-evs))
- Update filter examples and validate optional cases [\#227](https://github.com/Materials-Consortia/optimade-python-tools/pull/227) ([ml-evs](https://github.com/ml-evs))
- Switch from config init to BaseSettings [\#226](https://github.com/Materials-Consortia/optimade-python-tools/pull/226) ([shyamd](https://github.com/shyamd))

## [v0.7.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.7.1) (2020-03-16)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.7.0...v0.7.1)

**Closed issues:**

- Fix all capitalisation of OPTIMADE [\#232](https://github.com/Materials-Consortia/optimade-python-tools/issues/232)
- Remove validator action from README [\#230](https://github.com/Materials-Consortia/optimade-python-tools/issues/230)

**Merged pull requests:**

- Fix github actions for non-release tags [\#235](https://github.com/Materials-Consortia/optimade-python-tools/pull/235) ([shyamd](https://github.com/shyamd))
- Update OPTIMADE capitalisation [\#233](https://github.com/Materials-Consortia/optimade-python-tools/pull/233) ([ml-evs](https://github.com/ml-evs))
- Update mentions of action in readme [\#231](https://github.com/Materials-Consortia/optimade-python-tools/pull/231) ([ml-evs](https://github.com/ml-evs))

## [v0.7.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.7.0) (2020-03-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.6.0...v0.7.0)

**Implemented enhancements:**

- Validate all non-optional :filter: examples from the spec [\#213](https://github.com/Materials-Consortia/optimade-python-tools/pull/213) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Some mandatory filter examples from spec do not work [\#217](https://github.com/Materials-Consortia/optimade-python-tools/issues/217)
- Add txt-files in optimade.validator.data to MANIFEST [\#225](https://github.com/Materials-Consortia/optimade-python-tools/pull/225) ([CasperWA](https://github.com/CasperWA))
- Handle arbitrary nested NOT/AND/OR in queries [\#221](https://github.com/Materials-Consortia/optimade-python-tools/pull/221) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator only validates what we have working, not what is required by the spec [\#182](https://github.com/Materials-Consortia/optimade-python-tools/issues/182)

**Merged pull requests:**

- v0.7.0 release [\#228](https://github.com/Materials-Consortia/optimade-python-tools/pull/228) ([ml-evs](https://github.com/ml-evs))
- Remove GH Action to validate OPTiMaDe instances [\#224](https://github.com/Materials-Consortia/optimade-python-tools/pull/224) ([CasperWA](https://github.com/CasperWA))
- Codecov-action supports token-less uploads [\#220](https://github.com/Materials-Consortia/optimade-python-tools/pull/220) ([CasperWA](https://github.com/CasperWA))
- Update django requirement from \>=2.2.9,~=2.2 to \>=2.2,\<4.0 [\#219](https://github.com/Materials-Consortia/optimade-python-tools/pull/219) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))
- Update elasticsearch-dsl requirement from ~=6.4 to \>=6.4,\<8.0 [\#218](https://github.com/Materials-Consortia/optimade-python-tools/pull/218) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))

## [v0.6.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.6.0) (2020-03-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.5.0...v0.6.0)

**Implemented enhancements:**

- Possibly add CORS middleware [\#159](https://github.com/Materials-Consortia/optimade-python-tools/issues/159)
- Add debug flag to server [\#130](https://github.com/Materials-Consortia/optimade-python-tools/issues/130)
- Make validator GitHub Action [\#191](https://github.com/Materials-Consortia/optimade-python-tools/pull/191) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- meta/query/representation value not cutting off version properly [\#199](https://github.com/Materials-Consortia/optimade-python-tools/issues/199)
- URL for providers.json from Materials-Consortia has changed [\#186](https://github.com/Materials-Consortia/optimade-python-tools/issues/186)
- Relationships don't work when "/" present in id [\#181](https://github.com/Materials-Consortia/optimade-python-tools/issues/181)
- Redirect middleware not hitting single-entry endpoints [\#174](https://github.com/Materials-Consortia/optimade-python-tools/issues/174)

**Closed issues:**

- /info/ reports wrong url under available\_api\_versions [\#215](https://github.com/Materials-Consortia/optimade-python-tools/issues/215)
- Query parameters not handled correctly [\#208](https://github.com/Materials-Consortia/optimade-python-tools/issues/208)
- Test for AvailableApiVersion is correct for the wrong reasons [\#204](https://github.com/Materials-Consortia/optimade-python-tools/issues/204)
- Drop '/optimade' from paths in openapi.json [\#197](https://github.com/Materials-Consortia/optimade-python-tools/issues/197)
- heroku is failing [\#185](https://github.com/Materials-Consortia/optimade-python-tools/issues/185)
- List properties and HAS \_ operators missing [\#98](https://github.com/Materials-Consortia/optimade-python-tools/issues/98)
- Checklist for OPTiMaDe v0.10.1 [\#29](https://github.com/Materials-Consortia/optimade-python-tools/issues/29)

**Merged pull requests:**

- Removed /optimade/ prefix in info response [\#216](https://github.com/Materials-Consortia/optimade-python-tools/pull/216) ([ml-evs](https://github.com/ml-evs))
- Self load data [\#212](https://github.com/Materials-Consortia/optimade-python-tools/pull/212) ([shyamd](https://github.com/shyamd))
- Update tests for available\_api\_versions [\#211](https://github.com/Materials-Consortia/optimade-python-tools/pull/211) ([CasperWA](https://github.com/CasperWA))
- Up to v0.6.0 [\#210](https://github.com/Materials-Consortia/optimade-python-tools/pull/210) ([CasperWA](https://github.com/CasperWA))
- Update handling of include parameter \(and other query parameters\) [\#209](https://github.com/Materials-Consortia/optimade-python-tools/pull/209) ([CasperWA](https://github.com/CasperWA))
- Skip HAS ONLY test if mongomock version \<= 3.19.0 [\#206](https://github.com/Materials-Consortia/optimade-python-tools/pull/206) ([ml-evs](https://github.com/ml-evs))
- Test mandatory queries in validator [\#205](https://github.com/Materials-Consortia/optimade-python-tools/pull/205) ([ml-evs](https://github.com/ml-evs))
- Fix include query parameter [\#202](https://github.com/Materials-Consortia/optimade-python-tools/pull/202) ([CasperWA](https://github.com/CasperWA))
- Fix meta.query.representation and remove /optimade in base URLs [\#201](https://github.com/Materials-Consortia/optimade-python-tools/pull/201) ([CasperWA](https://github.com/CasperWA))
- Use mongo for CI [\#196](https://github.com/Materials-Consortia/optimade-python-tools/pull/196) ([ml-evs](https://github.com/ml-evs))
- \(Cosmetic\) updates to models [\#195](https://github.com/Materials-Consortia/optimade-python-tools/pull/195) ([CasperWA](https://github.com/CasperWA))
- Add CORSMiddleware [\#194](https://github.com/Materials-Consortia/optimade-python-tools/pull/194) ([CasperWA](https://github.com/CasperWA))
- Add "debug mode" [\#190](https://github.com/Materials-Consortia/optimade-python-tools/pull/190) ([CasperWA](https://github.com/CasperWA))
- Use https://provider.optimade.org/providers.json [\#187](https://github.com/Materials-Consortia/optimade-python-tools/pull/187) ([CasperWA](https://github.com/CasperWA))
- Fix errors parsing IDs that contain slashes [\#183](https://github.com/Materials-Consortia/optimade-python-tools/pull/183) ([ml-evs](https://github.com/ml-evs))
- Added default mongo implementations for HAS ALL/ANY/ONLY [\#173](https://github.com/Materials-Consortia/optimade-python-tools/pull/173) ([ml-evs](https://github.com/ml-evs))

## [v0.5.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.5.0) (2020-02-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.4.0...v0.5.0)

**Implemented enhancements:**

- Implement a landing page for requests to the base URL [\#169](https://github.com/Materials-Consortia/optimade-python-tools/issues/169)

**Fixed bugs:**

- 'minor' and 'patch' versioned base URL prefixes are wrong [\#177](https://github.com/Materials-Consortia/optimade-python-tools/issues/177)

**Closed issues:**

- Handle `include` standard JSON API query parameter [\#94](https://github.com/Materials-Consortia/optimade-python-tools/issues/94)

**Merged pull requests:**

- Bump to v0.5.0 [\#179](https://github.com/Materials-Consortia/optimade-python-tools/pull/179) ([CasperWA](https://github.com/CasperWA))
- Correctly create optional versioned base URLs [\#178](https://github.com/Materials-Consortia/optimade-python-tools/pull/178) ([CasperWA](https://github.com/CasperWA))
- Make mapper aliases configurable [\#175](https://github.com/Materials-Consortia/optimade-python-tools/pull/175) ([ml-evs](https://github.com/ml-evs))
- Add landing page at base URL [\#172](https://github.com/Materials-Consortia/optimade-python-tools/pull/172) ([ml-evs](https://github.com/ml-evs))
- Implement `include` query parameter [\#163](https://github.com/Materials-Consortia/optimade-python-tools/pull/163) ([CasperWA](https://github.com/CasperWA))
- Add docker for index meta-database [\#140](https://github.com/Materials-Consortia/optimade-python-tools/pull/140) ([CasperWA](https://github.com/CasperWA))

## [v0.4.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.4.0) (2020-02-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.4...v0.4.0)

**Implemented enhancements:**

- switch to pipenv? [\#37](https://github.com/Materials-Consortia/optimade-python-tools/issues/37)
- Reorder tests [\#162](https://github.com/Materials-Consortia/optimade-python-tools/pull/162) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Server app intermingles [\#161](https://github.com/Materials-Consortia/optimade-python-tools/issues/161)
- `response_fields` not working [\#154](https://github.com/Materials-Consortia/optimade-python-tools/issues/154)

**Closed issues:**

- Change `page_page` to `page_number` [\#165](https://github.com/Materials-Consortia/optimade-python-tools/issues/165)
- Add schema-relevant parameters to query parameters [\#164](https://github.com/Materials-Consortia/optimade-python-tools/issues/164)
- Alias optimade/structures/ to optimade/structure [\#128](https://github.com/Materials-Consortia/optimade-python-tools/issues/128)
- Minor changes to specification v0.10.1-develop [\#115](https://github.com/Materials-Consortia/optimade-python-tools/issues/115)
- Update models with new levels of REQUIRED response properties [\#114](https://github.com/Materials-Consortia/optimade-python-tools/issues/114)
- Constraining list/array types in the schema [\#55](https://github.com/Materials-Consortia/optimade-python-tools/issues/55)

**Merged pull requests:**

- Bump to v0.4.0 [\#168](https://github.com/Materials-Consortia/optimade-python-tools/pull/168) ([CasperWA](https://github.com/CasperWA))
- Describe query parameters in OpenAPI schema [\#166](https://github.com/Materials-Consortia/optimade-python-tools/pull/166) ([CasperWA](https://github.com/CasperWA))
- Redirect slashed URLs [\#160](https://github.com/Materials-Consortia/optimade-python-tools/pull/160) ([CasperWA](https://github.com/CasperWA))
- New REQUIRED level properties [\#153](https://github.com/Materials-Consortia/optimade-python-tools/pull/153) ([CasperWA](https://github.com/CasperWA))

## [v0.3.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.4) (2020-02-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.3...v0.3.4)

**Implemented enhancements:**

- Include `develop` or not? Default branch? - Create INSTALL.md [\#136](https://github.com/Materials-Consortia/optimade-python-tools/issues/136)

**Fixed bugs:**

- Excepting non-existent exception [\#129](https://github.com/Materials-Consortia/optimade-python-tools/issues/129)

**Closed issues:**

- disable serving API under /v0.10 and /v0.10.0 by default? [\#122](https://github.com/Materials-Consortia/optimade-python-tools/issues/122)
- PyPI release checklist [\#67](https://github.com/Materials-Consortia/optimade-python-tools/issues/67)

**Merged pull requests:**

- Bump to v0.3.4 [\#158](https://github.com/Materials-Consortia/optimade-python-tools/pull/158) ([CasperWA](https://github.com/CasperWA))
- Fix heroku badge [\#157](https://github.com/Materials-Consortia/optimade-python-tools/pull/157) ([ml-evs](https://github.com/ml-evs))
- Move installation instructions [\#156](https://github.com/Materials-Consortia/optimade-python-tools/pull/156) ([ml-evs](https://github.com/ml-evs))
- Update base URLs [\#155](https://github.com/Materials-Consortia/optimade-python-tools/pull/155) ([CasperWA](https://github.com/CasperWA))
- Extend OpenAPI/spec description [\#151](https://github.com/Materials-Consortia/optimade-python-tools/pull/151) ([CasperWA](https://github.com/CasperWA))
- Non Local Mongo [\#150](https://github.com/Materials-Consortia/optimade-python-tools/pull/150) ([shyamd](https://github.com/shyamd))

## [v0.3.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.3) (2020-01-24)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.2...v0.3.3)

**Fixed bugs:**

- Lark files not being distributed [\#141](https://github.com/Materials-Consortia/optimade-python-tools/issues/141)

**Closed issues:**

- Tests fail with lark-parser\>=0.8 [\#146](https://github.com/Materials-Consortia/optimade-python-tools/issues/146)

**Merged pull requests:**

- Updated lark-parser to 0.8.1 [\#149](https://github.com/Materials-Consortia/optimade-python-tools/pull/149) ([ml-evs](https://github.com/ml-evs))
- Split eager and standard tests to avoid unnecessary badge of shame [\#148](https://github.com/Materials-Consortia/optimade-python-tools/pull/148) ([ml-evs](https://github.com/ml-evs))
- Bump to v0.3.3 [\#147](https://github.com/Materials-Consortia/optimade-python-tools/pull/147) ([CasperWA](https://github.com/CasperWA))
- Fix root\_validator issues with optional fields and made meta optional [\#145](https://github.com/Materials-Consortia/optimade-python-tools/pull/145) ([ml-evs](https://github.com/ml-evs))
- Handle `JSONDecodeError`s in validator [\#144](https://github.com/Materials-Consortia/optimade-python-tools/pull/144) ([ml-evs](https://github.com/ml-evs))

## [v0.3.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.2) (2020-01-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.1...v0.3.2)

**Implemented enhancements:**

- Add base URL to configuration file [\#135](https://github.com/Materials-Consortia/optimade-python-tools/pull/135) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Fix `load_from_json` [\#137](https://github.com/Materials-Consortia/optimade-python-tools/pull/137) ([CasperWA](https://github.com/CasperWA))

**Merged pull requests:**

- Make sure relevant package data is included in distributions [\#142](https://github.com/Materials-Consortia/optimade-python-tools/pull/142) ([CasperWA](https://github.com/CasperWA))
- Add database page limit [\#139](https://github.com/Materials-Consortia/optimade-python-tools/pull/139) ([CasperWA](https://github.com/CasperWA))

## [v0.3.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.1) (2020-01-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.0...v0.3.1)

**Merged pull requests:**

- Update requirements [\#138](https://github.com/Materials-Consortia/optimade-python-tools/pull/138) ([CasperWA](https://github.com/CasperWA))

## [v0.3.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.0) (2020-01-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.2...v0.3.0)

**Implemented enhancements:**

- Implement optional `implementation` in top-level meta response [\#117](https://github.com/Materials-Consortia/optimade-python-tools/issues/117)
- Create "special" index meta-database server [\#100](https://github.com/Materials-Consortia/optimade-python-tools/issues/100)
- Implement relationships in server [\#71](https://github.com/Materials-Consortia/optimade-python-tools/issues/71)
- Add missing /references endpoint to server [\#69](https://github.com/Materials-Consortia/optimade-python-tools/issues/69)
- Automatically publish version tags to PyPI via GH Actions [\#107](https://github.com/Materials-Consortia/optimade-python-tools/pull/107) ([CasperWA](https://github.com/CasperWA))
- Using routers [\#99](https://github.com/Materials-Consortia/optimade-python-tools/pull/99) ([CasperWA](https://github.com/CasperWA))
- Add relationships functionality [\#91](https://github.com/Materials-Consortia/optimade-python-tools/pull/91) ([ml-evs](https://github.com/ml-evs))
- Added external API validator based on our pydantic models [\#74](https://github.com/Materials-Consortia/optimade-python-tools/pull/74) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- The invoke task `update-openapijson` is incomplete [\#123](https://github.com/Materials-Consortia/optimade-python-tools/issues/123)
- Django vulnerability [\#108](https://github.com/Materials-Consortia/optimade-python-tools/issues/108)

**Closed issues:**

- info endpoint duplicated? [\#120](https://github.com/Materials-Consortia/optimade-python-tools/issues/120)
- Commented-out validator [\#111](https://github.com/Materials-Consortia/optimade-python-tools/issues/111)
- FastAPI v0.44.0 supports pydantic \> 1.0.0 [\#101](https://github.com/Materials-Consortia/optimade-python-tools/issues/101)
- Server is missing /links endpoint [\#89](https://github.com/Materials-Consortia/optimade-python-tools/issues/89)
- Make sure all validators are tested [\#87](https://github.com/Materials-Consortia/optimade-python-tools/issues/87)
- The `sortable` field must be added to models [\#84](https://github.com/Materials-Consortia/optimade-python-tools/issues/84)
- Package structure [\#72](https://github.com/Materials-Consortia/optimade-python-tools/issues/72)
- Possibly make /info/{endpoint} dynamic [\#70](https://github.com/Materials-Consortia/optimade-python-tools/issues/70)
- setuptools package with server as "extra" [\#62](https://github.com/Materials-Consortia/optimade-python-tools/issues/62)
- use examples from specs as resources [\#57](https://github.com/Materials-Consortia/optimade-python-tools/issues/57)
- httptools dependency has build issues on GCC/Linux [\#54](https://github.com/Materials-Consortia/optimade-python-tools/issues/54)
- Lark grammar file for v0.9.8 [\#50](https://github.com/Materials-Consortia/optimade-python-tools/issues/50)
- type is missing in response [\#43](https://github.com/Materials-Consortia/optimade-python-tools/issues/43)
- Enforce use of autoformatter  [\#33](https://github.com/Materials-Consortia/optimade-python-tools/issues/33)
- switch license to MIT [\#28](https://github.com/Materials-Consortia/optimade-python-tools/issues/28)
- write a lark JSONTransformer / JSONdecoder [\#26](https://github.com/Materials-Consortia/optimade-python-tools/issues/26)
- server.jsonapi has no additionalProperties=false  [\#23](https://github.com/Materials-Consortia/optimade-python-tools/issues/23)
- server.jsonapi has no patternProperties  [\#22](https://github.com/Materials-Consortia/optimade-python-tools/issues/22)
- Developer-friendly pre-commit openapi.json visual diff [\#21](https://github.com/Materials-Consortia/optimade-python-tools/issues/21)
- add JSON schema API [\#12](https://github.com/Materials-Consortia/optimade-python-tools/issues/12)
- generate static documentation on github from openapi.json [\#9](https://github.com/Materials-Consortia/optimade-python-tools/issues/9)
- test how to generate a client from the openapi.json [\#8](https://github.com/Materials-Consortia/optimade-python-tools/issues/8)
- come up with suggested toolchain for validating existing optimade API against openapi.json [\#7](https://github.com/Materials-Consortia/optimade-python-tools/issues/7)
- add travis test that checks openapi.json is valid OpenAPI spec [\#6](https://github.com/Materials-Consortia/optimade-python-tools/issues/6)
- add 2 examples of how to include documentation in python classes [\#5](https://github.com/Materials-Consortia/optimade-python-tools/issues/5)
- add one-line command to update openapi.json [\#4](https://github.com/Materials-Consortia/optimade-python-tools/issues/4)

**Merged pull requests:**

- Fixed CI readme badge [\#133](https://github.com/Materials-Consortia/optimade-python-tools/pull/133) ([ml-evs](https://github.com/ml-evs))
- Add meta.description to BaseRelationshipResource [\#131](https://github.com/Materials-Consortia/optimade-python-tools/pull/131) ([CasperWA](https://github.com/CasperWA))
- Added homepage attribute to LinksResource [\#127](https://github.com/Materials-Consortia/optimade-python-tools/pull/127) ([ml-evs](https://github.com/ml-evs))
- Updated structure models and validators [\#126](https://github.com/Materials-Consortia/optimade-python-tools/pull/126) ([ml-evs](https://github.com/ml-evs))
- Minor change to fallback server.cfg [\#125](https://github.com/Materials-Consortia/optimade-python-tools/pull/125) ([ml-evs](https://github.com/ml-evs))
- Update local OpenAPI schemes prior to copying [\#124](https://github.com/Materials-Consortia/optimade-python-tools/pull/124) ([CasperWA](https://github.com/CasperWA))
- Update OpenAPI tags [\#121](https://github.com/Materials-Consortia/optimade-python-tools/pull/121) ([CasperWA](https://github.com/CasperWA))
- A few fixes related to usage as a library [\#119](https://github.com/Materials-Consortia/optimade-python-tools/pull/119) ([ml-evs](https://github.com/ml-evs))
- Add implementation to top-level meta response [\#118](https://github.com/Materials-Consortia/optimade-python-tools/pull/118) ([CasperWA](https://github.com/CasperWA))
- Add heroku deployment scripts [\#116](https://github.com/Materials-Consortia/optimade-python-tools/pull/116) ([ltalirz](https://github.com/ltalirz))
- Reorganize package [\#113](https://github.com/Materials-Consortia/optimade-python-tools/pull/113) ([CasperWA](https://github.com/CasperWA))
- Introduce grammar v0.10.1 [\#112](https://github.com/Materials-Consortia/optimade-python-tools/pull/112) ([CasperWA](https://github.com/CasperWA))
- Update to pydantic v1 [\#110](https://github.com/Materials-Consortia/optimade-python-tools/pull/110) ([CasperWA](https://github.com/CasperWA))
- Minimum requirement of django v2.2.8 [\#109](https://github.com/Materials-Consortia/optimade-python-tools/pull/109) ([CasperWA](https://github.com/CasperWA))
- Index meta-database [\#103](https://github.com/Materials-Consortia/optimade-python-tools/pull/103) ([CasperWA](https://github.com/CasperWA))
- restrict pydantic version [\#97](https://github.com/Materials-Consortia/optimade-python-tools/pull/97) ([ltalirz](https://github.com/ltalirz))
- Add /links [\#95](https://github.com/Materials-Consortia/optimade-python-tools/pull/95) ([CasperWA](https://github.com/CasperWA))
- Fix data\_returned and data\_available [\#93](https://github.com/Materials-Consortia/optimade-python-tools/pull/93) ([CasperWA](https://github.com/CasperWA))
- Use GitHub Actions for CI [\#92](https://github.com/Materials-Consortia/optimade-python-tools/pull/92) ([ml-evs](https://github.com/ml-evs))
- Remove inappropriate lint messages [\#90](https://github.com/Materials-Consortia/optimade-python-tools/pull/90) ([CasperWA](https://github.com/CasperWA))
- Fix dependencies [\#88](https://github.com/Materials-Consortia/optimade-python-tools/pull/88) ([CasperWA](https://github.com/CasperWA))
- Add sortable field to EntryInfoProperty model [\#85](https://github.com/Materials-Consortia/optimade-python-tools/pull/85) ([CasperWA](https://github.com/CasperWA))
- Validate illegal fields are not present under attributes and relationships [\#83](https://github.com/Materials-Consortia/optimade-python-tools/pull/83) ([CasperWA](https://github.com/CasperWA))
- Add references endpoint [\#78](https://github.com/Materials-Consortia/optimade-python-tools/pull/78) ([CasperWA](https://github.com/CasperWA))
- fix travis build [\#77](https://github.com/Materials-Consortia/optimade-python-tools/pull/77) ([ltalirz](https://github.com/ltalirz))
- Fix manual verification of elements\_ratios [\#76](https://github.com/Materials-Consortia/optimade-python-tools/pull/76) ([CasperWA](https://github.com/CasperWA))
- add automatic PyPI deployment [\#75](https://github.com/Materials-Consortia/optimade-python-tools/pull/75) ([ltalirz](https://github.com/ltalirz))
- Remove reference to `"all"` endpoint and rename collections submodule [\#73](https://github.com/Materials-Consortia/optimade-python-tools/pull/73) ([ml-evs](https://github.com/ml-evs))
- Updates to README and docs for v0.10.0 [\#68](https://github.com/Materials-Consortia/optimade-python-tools/pull/68) ([ml-evs](https://github.com/ml-evs))
- Adding grammar for v0.10.0 [\#66](https://github.com/Materials-Consortia/optimade-python-tools/pull/66) ([fekad](https://github.com/fekad))
- Schema updates and fixes relative to the v0.10.0 spec [\#65](https://github.com/Materials-Consortia/optimade-python-tools/pull/65) ([ml-evs](https://github.com/ml-evs))
- Break requirements down on per backend basis [\#64](https://github.com/Materials-Consortia/optimade-python-tools/pull/64) ([ml-evs](https://github.com/ml-evs))
- 0.10.0 grammer, elasticsearch transformer, setuptools extra [\#63](https://github.com/Materials-Consortia/optimade-python-tools/pull/63) ([markus1978](https://github.com/markus1978))
- Added a Lark to Django Query converter [\#61](https://github.com/Materials-Consortia/optimade-python-tools/pull/61) ([tachyontraveler](https://github.com/tachyontraveler))
- Some minor fixes [\#60](https://github.com/Materials-Consortia/optimade-python-tools/pull/60) ([ml-evs](https://github.com/ml-evs))
- Added codecov to CI [\#59](https://github.com/Materials-Consortia/optimade-python-tools/pull/59) ([ml-evs](https://github.com/ml-evs))
- Enforce black via `pre-commit` tool [\#53](https://github.com/Materials-Consortia/optimade-python-tools/pull/53) ([dwinston](https://github.com/dwinston))
- Update setup.py and version [\#51](https://github.com/Materials-Consortia/optimade-python-tools/pull/51) ([dwinston](https://github.com/dwinston))
- /structure/info endpoint [\#49](https://github.com/Materials-Consortia/optimade-python-tools/pull/49) ([fawzi](https://github.com/fawzi))
- add constrained list type [\#48](https://github.com/Materials-Consortia/optimade-python-tools/pull/48) ([dwinston](https://github.com/dwinston))
- Refactored into submodules and added test data [\#47](https://github.com/Materials-Consortia/optimade-python-tools/pull/47) ([ml-evs](https://github.com/ml-evs))
- Update structure endpoint to pre-alpha 0.10 spec [\#45](https://github.com/Materials-Consortia/optimade-python-tools/pull/45) ([ltalirz](https://github.com/ltalirz))
- Adding Resource Links [\#44](https://github.com/Materials-Consortia/optimade-python-tools/pull/44) ([tpurcell90](https://github.com/tpurcell90))
- Reblacken [\#42](https://github.com/Materials-Consortia/optimade-python-tools/pull/42) ([ml-evs](https://github.com/ml-evs))
- Documented json [\#41](https://github.com/Materials-Consortia/optimade-python-tools/pull/41) ([tpurcell90](https://github.com/tpurcell90))
- fix example output [\#40](https://github.com/Materials-Consortia/optimade-python-tools/pull/40) ([dwinston](https://github.com/dwinston))
- use jsonapi better at top level, add error response [\#36](https://github.com/Materials-Consortia/optimade-python-tools/pull/36) ([fawzi](https://github.com/fawzi))
- add JSONTransformer [\#35](https://github.com/Materials-Consortia/optimade-python-tools/pull/35) ([dwinston](https://github.com/dwinston))
- switch to MIT license [\#34](https://github.com/Materials-Consortia/optimade-python-tools/pull/34) ([ltalirz](https://github.com/ltalirz))
- Updated entry definitions and renamed Response classes [\#32](https://github.com/Materials-Consortia/optimade-python-tools/pull/32) ([ml-evs](https://github.com/ml-evs))
- update readme [\#31](https://github.com/Materials-Consortia/optimade-python-tools/pull/31) ([ltalirz](https://github.com/ltalirz))
- Seperated Links from JSON API into its own file [\#30](https://github.com/Materials-Consortia/optimade-python-tools/pull/30) ([tpurcell90](https://github.com/tpurcell90))
- simplify schema update [\#27](https://github.com/Materials-Consortia/optimade-python-tools/pull/27) ([ltalirz](https://github.com/ltalirz))
- add openapi\_diff to travis [\#25](https://github.com/Materials-Consortia/optimade-python-tools/pull/25) ([ltalirz](https://github.com/ltalirz))
- Json api add [\#24](https://github.com/Materials-Consortia/optimade-python-tools/pull/24) ([tpurcell90](https://github.com/tpurcell90))
- Added JSON diff test [\#20](https://github.com/Materials-Consortia/optimade-python-tools/pull/20) ([ml-evs](https://github.com/ml-evs))
- info endpoint [\#19](https://github.com/Materials-Consortia/optimade-python-tools/pull/19) ([fawzi](https://github.com/fawzi))
- adding run.sh script to start webserver [\#18](https://github.com/Materials-Consortia/optimade-python-tools/pull/18) ([fawzi](https://github.com/fawzi))
- error response [\#17](https://github.com/Materials-Consortia/optimade-python-tools/pull/17) ([fawzi](https://github.com/fawzi))
- Links can be strings [\#16](https://github.com/Materials-Consortia/optimade-python-tools/pull/16) ([fawzi](https://github.com/fawzi))
- response should be either many \(list\) or one \(object\), not an union [\#15](https://github.com/Materials-Consortia/optimade-python-tools/pull/15) ([fawzi](https://github.com/fawzi))
- reorg models [\#14](https://github.com/Materials-Consortia/optimade-python-tools/pull/14) ([dwinston](https://github.com/dwinston))
- Update the OptimadeMetaResponse to development schema [\#13](https://github.com/Materials-Consortia/optimade-python-tools/pull/13) ([ml-evs](https://github.com/ml-evs))
- add openapi spec validator [\#10](https://github.com/Materials-Consortia/optimade-python-tools/pull/10) ([ltalirz](https://github.com/ltalirz))
- fix test data download [\#3](https://github.com/Materials-Consortia/optimade-python-tools/pull/3) ([ltalirz](https://github.com/ltalirz))
- \[WIP\] Mongoconverter [\#1](https://github.com/Materials-Consortia/optimade-python-tools/pull/1) ([wuxiaohua1011](https://github.com/wuxiaohua1011))

## [v0.1.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.2) (2018-06-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.1...v0.1.2)

## [v0.1.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.1) (2018-06-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.0...v0.1.1)

## [v0.1.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.0) (2018-06-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/6680fb28a60ec4ff43a303b7b4dbf41e159a25b6...v0.1.0)



\* *This Changelog was automatically generated by [github_changelog_generator](https://github.com/github-changelog-generator/github-changelog-generator)*
<img width="10%" align="left" src="images/optimade_logo_180x180.svg">

# OPTIMADE Python tools

| Latest release | Build status | Activity |
|:--------------:|:------------:|:--------:|
| [![PyPI Version](https://img.shields.io/pypi/v/optimade?logo=pypi&logoColor=white)](https://pypi.org/project/optimade/)<br>[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/optimade?logo=python&logoColor=white)](https://pypi.org/project/optimade/)<br>[![OPTIMADE](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/optimade-version.json)](https://github.com/Materials-Consortia/OPTIMADE/) | [![Build Status](https://img.shields.io/github/workflow/status/Materials-Consortia/optimade-python-tools/CI%20tests?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/actions?query=branch%3Amaster+)<br>[![codecov](https://img.shields.io/codecov/c/github/Materials-Consortia/optimade-python-tools?logo=codecov&logoColor=white&token=UJAtmqkZZO)](https://codecov.io/gh/Materials-Consortia/optimade-python-tools)<br>[![Heroku App Status](https://heroku-shields.herokuapp.com/optimade??logo=heroku)](https://optimade.herokuapp.com) | [![Commit Activity](https://img.shields.io/github/commit-activity/m/Materials-Consortia/optimade-python-tools?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/pulse)<br>[![Last Commit](https://img.shields.io/github/last-commit/Materials-Consortia/optimade-python-tools/master?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/commits/master)<br>[![Contributors](https://badgen.net/github/contributors/Materials-Consortia/optimade-python-tools?icon=github)](https://github.com/Materials-Consortia/optimade-python-tools/graphs/contributors) |

The aim of OPTIMADE is to develop a common API, compliant with the [JSON:API 1.0](http://jsonapi.org/format/1.0/) specification.
This is to enable interoperability among databases that contain calculated properties of existing and hypothetical materials.

This repository contains a library of tools for implementing and consuming [OPTIMADE](https://www.optimade.org) APIs using Python.
Server implementations can make use of the supported MongoDB (v4) and Elasticsearch (v6) database backends, or plug in a custom backend implementation.
The package also contains a server validator tool, which may be called from the shell (`optimade-validator`) or used as a GitHub Action from [optimade-validator-action](https://github.com/Materials-Consortia/optimade-validator-action).

The release history and changelog can be found in [the changelog](CHANGELOG.md).

## Documentation

This document, guides, and the full module API documentation can be found online at [https://optimade.org/optimade-python-tools](https://optimade.org/optimade-python-tools).
In particular, documentation of the OPTIMADE API response data models (implemented here with [pydantic](https://github.com/samuelcolvin/pydantic)) can be found online under [OPTIMADE Data Models](https://optimade.org/optimade-python-tools/all_models).

## Installation

Detailed installation instructions for different use cases (e.g., using the library or running a server) can be found in [the installation documentation](INSTALL.md).

The latest stable version of this package can be obtained from [PyPI](https://pypi.org/project/optimade) `pip install optimade`.
The latest development version of this package can be installed from the master branch of this repository `git clone https://github.com/Materials-Consortia/optimade-python-tools`.

## Supported OPTIMADE versions

Each release of the `optimade` package from this repository only targets one version of the OPTIMADE specification, summarised in the table below.

| OPTIMADE API version | `optimade` version |
|:--------------------:|:------------------:|
| [v1.0.0](https://github.com/Materials-Consortia/OPTIMADE/blob/v1.0.0/optimade.rst) | [v0.12.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.9) |
| [v1.1.0](https://github.com/Materials-Consortia/OPTIMADE/blob/v1.1.0/optimade.rst) | [v0.16.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.0) |

## Contributing and Getting Help

All development of this package (bug reports, suggestions, feedback and pull requests) occurs in the [optimade-python-tools GitHub repository](https://github.com/Materials-Consortia/optimade-python-tools).
Contribution guidelines and tips for getting help can be found in the [contributing notes](CONTRIBUTING.md).

## How to cite

If you use this package to access or serve OPTIMADE data, we kindly request that you consider citing the following:

- Andersen *et al.*, OPTIMADE, an API for exchanging materials data, *Sci. Data* **8**, 217 (2021) [10.1038/s41597-021-00974-z](https://doi.org/10.1038/s41597-021-00974-z)
- Evans *et al.*, optimade-python-tools: a Python library for serving and consuming materials data via OPTIMADE APIs. *Journal of Open Source Software*, **6**(65), 3458 (2021) [10.21105/joss.03458](https://doi.org/10.21105/joss.03458)

## Links

- [OPTIMADE Specification](https://github.com/Materials-Consortia/OPTIMADE/blob/develop/optimade.rst), the human-readable specification that this library is based on.
- [optimade-validator-action](https://github.com/Materials-Consortia/optimade-validator-action), a GitHub action that can be used to validate implementations from a URL (using the validator from this repo).
- [OpenAPI](https://github.com/OAI/OpenAPI-Specification), the machine-readable format used to specify the OPTIMADE API in [`openapi.json`](openapi/openapi.json) and [`index_openapi.json`](openapi/index_openapi.json).
- [Interactive documentation](https://petstore.swagger.io/?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/openapi/openapi.json) generated from [`openapi.json`](openapi/openapi.json) (see also [interactive JSON editor](https://editor.swagger.io/?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/openapi/openapi.json)).
- [pydantic](https://pydantic-docs.helpmanual.io/), the library used for generating the OpenAPI schema from [Python models](https://www.optimade.org/optimade-python-tools/all_models/).
- [FastAPI](https://fastapi.tiangolo.com/), the framework used for generating the reference implementation expressed by the [`openapi.json`](openapi/openapi.json) specification.
- [lark](https://github.com/lark-parser/lark), the library used to parse the filter language in OPTIMADE queries.
# Contributing and getting help

If you run into any problems using this package, or if you have a question, suggestion or feedback, then please raise an issue on [GitHub](https://github.com/Materials-Consortia/optimade-python-tools/issues/new).

The [Materials Consortia](https://github.com/Materials-Consortia) is very open to contributions across all of its packages.
This may be anything from simple feedback and raising [new issues](https://github.com/Materials-Consortia/optimade-python-tools/issues/new) to creating [new PRs](https://github.com/Materials-Consortia/optimade-python-tools/compare).

If you are interested in contributing but don't know where to begin, some issues have been marked with the [good first issue](https://github.com/Materials-Consortia/optimade-python-tools/labels/good%20first%20issue) label, typically where an isolated enhancement has a concrete suggestion.
Simply add a comment under an issue if you are interested in tackling it!

Recommendations for setting up a development environment for this package can be found in the [Installation instructions](https://www.optimade.org/optimade-python-tools/INSTALL/#full-development-installation).

More broadly, if you would like to ask questions or contact the consortium about creating an OPTIMADE implementation for a new database, then please read the relevant "get involved" section on the [OPTIMADE website](https://www.optimade.org/#get-involved).
# Installation

This package can be installed from PyPI, or by cloning the repository, depending on your use-case.

1. To use the `optimade` Python package as a library, (e.g., using the models for validation, parsing filters with the grammar, or using the command-line tool `optimade-validator` tool), it is recommended that you install the latest release of the package from PyPI with `pip install optimade`.
2. If you want to run, use or modify the reference server implementation, then it is recommended that you clone this repository and install it from your local files (with `pip install .`, or `pip install -e .` for an editable installation).

## The index meta-database

This package may be used to setup and run an [OPTIMADE index meta-database](https://github.com/Materials-Consortia/OPTIMADE/blob/develop/optimade.rst#index-meta-database).
Clone this repository and install the package locally with `pip install -e .[server]`.

There is a built-in index meta-database set up to populate a `mongomock` in-memory database with resources from a static `json` file containing the `child` resources you, as a database provider, want to serve under this index meta-database.
The location of that `json` file is controllable using the `index_links_path` property of the configuration or setting via the environment variable `optimade_index_links_path`.

Running the index meta-database is then as simple as writing `./run.sh index` in a terminal from the root of this package.
You can find it at the base URL: <http://localhost:5001/v1>.

Here is an example of how it may look to start your server:

```sh
:~$ export OPTIMADE_CONFIG_FILE=/home/optimade_server/config.json
:~$ ./path/to/optimade/run.sh index
```

## Full development installation

The dependencies of this package can be found in `setup.py` with their latest supported versions.
By default, a minimal set of requirements are installed to work with the filter language and the `pydantic` models.
After cloning the repository, the install mode `server` (i.e. `pip install .[server]`) is sufficient to run a `uvicorn` server using the `mongomock` backend (or MongoDB with `pymongo`, if present).
The suite of development and testing tools are installed with via the install modes `dev` and `testing`.
There are additionally two backend-specific install modes, `elastic` and `mongo`, as well as the `all` mode, which installs all dependencies.
All contributed Python code, must use the [black](https://github.com/ambv/black) code formatter, and must pass the [flake8](http://flake8.pycqa.org/en/latest/) linter that is run automatically on all PRs.

```sh
# Clone this repository to your computer
git clone git@github.com:Materials-Consortia/optimade-python-tools.git
cd optimade-python-tools

# Ensure a Python>=3.7 (virtual) environment (example below using Anaconda/Miniconda)
conda create -n optimade python=3.7
conda activate optimade

# Install package and dependencies in editable mode (including "dev" requirements).
pip install -e ".[dev]"

# Optional: Retrieve the list of OPTIMADE providers. (Without this submodule, some of the tests will fail because "providers.json" cannot be found.)
git submodule update --init

# Run the tests with pytest
py.test

# Install pre-commit environment (e.g., auto-formats code on `git commit`)
pre-commit install

# Optional: Install MongoDB (and set `database_backend = mongodb`)
# Below method installs in conda environment and
# - starts server in background
# - ensures and uses ~/dbdata directory to store data
conda install -c anaconda mongodb
mkdir -p ~/dbdata && mongod --dbpath ~/dbdata --syslog --fork

# Start a development server (auto-reload on file changes at http://localhost:5000
# You can also execute ./run.sh
uvicorn optimade.server.main:app --reload --port 5000

# View auto-generated docs
open http://localhost:5000/docs
# View Open API Schema
open http://localhost:5000/openapi.json
```

When developing, you can run both the server and an index meta-database server at the same time (from two separate terminals).
Running the following:

```shell
./run.sh index
# or
uvicorn optimade.server.main_index:app --reload --port 5001
```

will run the index meta-database server at <http://localhost:5001/v1>.

## Testing specific backends

In order to run the test suite for a specific backend, the
`OPTIMADE_DATABASE_BACKEND` [environment variable (or config
option)](https://www.optimade.org/optimade-python-tools/configuration/) can be
set to one of `'mongodb'`, `'mongomock'` or `'elastic'` (see
[`ServerConfig.database_backend`][optimade.server.config.ServerConfig.database_backend]).
Tests for the two "real" database backends, MongoDB and Elasticsearch, require a writable, temporary database to be accessible.

The easiest way to deploy these databases and run the tests is with Docker, as shown below.
[Docker installation instructions](https://docs.docker.com/engine/install/) will depend on your system; on Linux, the `docker` commands below may need to be prepended with `sudo`, depending on your distribution.
These commands should be run from a local optimade-python-tools directory.

The following command starts a local Elasticsearch v6 instance, runs the test suite, then stops and deletes the containers (required as the tests insert some data):
```shell
docker run -d --name elasticsearch_test -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" elasticsearch:6.8.22 \
&& sleep 10 \
&& OPTIMADE_DATABASE_BACKEND="elastic" py.test; \
docker container stop elasticsearch_test; docker container rm elasticsearch_test
```

The following command starts a local MongoDB instance, runs the test suite, then stops and deletes the containers:
```shell
docker run -d --name mongo_test -p 27017:27017 -d mongo:4.4.6 \
&& OPTIMADE_DATABASE_BACKEND="mongodb" py.test; \
docker container stop mongo_test; docker container rm mongo_test
```
---
title: '`optimade-python-tools`: a Python library for serving and consuming materials data via OPTIMADE APIs'
tags:
  - Python
  - REST API
  - JSON:API
  - OPTIMADE API
  - crystallography
  - density-functional theory
  - ab initio
  - materials discovery
  - databases
authors:
  - name: Matthew L. Evans^[corresponding author]
    orcid:  0000-0002-1182-9098
    affiliation: "1, 2"
  - name: Casper W. Andersen^[co-first author]
    orcid: 0000-0002-2547-155X
    affiliation: 3
  - name: Shyam Dwaraknath
    orcid: 0000-0003-0289-2607
    affiliation: 4
  - name: Markus Scheidgen
    orcid: 0000-0002-8038-2277
    affiliation: "5, 6"
  - name: Ádám Fekete
    orcid: 0000-0002-6263-897X
    affiliation: "1, 8, 9"
  - name: Donald Winston
    orcid: 0000-0002-8424-0604
    affiliation: "4, 7"
affiliations:
 - name: Institut de la Matière Condensée et des Nanosciences, Université catholique de Louvain, Chemin des Étoiles 8, Louvain-la-Neuve 1348, Belgium
   index: 1
 - name: Theory of Condensed Matter Group, Cavendish Laboratory, University of Cambridge, J. J. Thomson Avenue, Cambridge, CB3 0HE, United Kingdom
   index: 2
 - name: Theory and Simulation of Materials (THEOS), Faculté des Sciences et Techniques de l'Ingénieur, École Polytechnique Fédérale de Lausanne, CH-1015 Lausanne, Switzerland
   index: 3
 - name: Lawrence Berkeley National Laboratory, Berkeley, CA, USA
   index: 4
 - name: Fritz-Haber-Institut der Max-Planck-Gesellschaft, Faradayweg 4-6, 14195, Berlin, Germany
   index: 5
 - name: Humboldt-Universität zu Berlin, Institut für Physik and IRIS Adlershof, 12489 Berlin, Germany
   index: 6
 - name: Polyneme LLC, New York, NY, USA
   index: 7
 - name: Department of Physics, King's College London, Strand, London WC2R 2LS, United Kingdom
   index: 8
 - name: Department of Physics and Namur Institute of Structured Materials, University of Namur, Rue de Bruxelles 51, 5000 Namur, Belgium
   index: 9

date: June 2021
bibliography: paper.bib
---

# Summary

In recent decades, improvements in algorithms, hardware, and theory have enabled crystalline materials to be studied computationally at the atomistic level with great accuracy and speed.
To enable dissemination, reproducibility, and reuse, many digital crystal structure databases have been created and curated, ready for comparison with existing infrastructure that stores structural characterizations (e.g., diffraction) of real crystals.
Each database will typically have a bespoke, stateless, web-based Application Programming Interface (API); users can submit a query via specially-crafted URLs.
Such esoteric and specialized APIs incur maintenance and usability costs upon both the data providers and consumers, who may not be software specialists.

The [OPTIMADE API](https://optimade.org) specification [@andersen2021optimade; @OPTIMADE_spec], released in July 2020, aimed to reduce these costs by designing a common API for use across a consortium of collaborating materials databases and beyond.
Whilst based on the robust JSON:API standard [@JSONAPI], the OPTIMADE API specification presents several domain-specific features and requirements that can be tricky to implement for non-specialist teams.
The repository presented here, `optimade-python-tools`, provides a modular reference server implementation and a set of associated tools to accelerate the development process for data providers, toolmakers and end-users.

# Statement of need

In order to accommodate existing materials database APIs, the OPTIMADE API specification allows for flexibility in the specific data served, but enforces a simple yet domain-specific filter language on well-defined resources.
However, this flexibility could be daunting to database providers, likely acting to increase the barrier to hosting an OPTIMADE API.
`optimade-python-tools` aims to catalyse the creation of APIs from existing and new data sources by providing a configurable and modular reference server implementation for hosting materials data in an OPTIMADE-compliant way.
The repository hosts the `optimade` Python package, which leverages the modern Python libraries pydantic [@pydantic] and FastAPI [@FastAPI] to specify the data models and API routes defined in the OPTIMADE API specification, additionally providing a schema following the OpenAPI format [@OpenAPI].
As this package was developed concomitantly with the OPTIMADE specification itself, the present authors are not aware of any other generic packages with similar functionality.
Two storage back-ends are supported out of the box, with full filter support for databases that employ the popular [MongoDB](https://www.mongodb.com) or [Elasticsearch](https://elastic.co) frameworks.

# Functionality

The modular functionality of `optimade` can be broken down by the different stages of a user query to the reference server.
Consider the following query URL to an OPTIMADE API, which should filter for any crystal structures in the database with a composition that consists of any three elements in a 1:1:1 ratio:

```
https://example.org/v1/structures?filter=chemical_formula_anonymous="ABC"
```

1. After routing the query to the appropriate `/structures/` endpoint adhering to version `v1` of the specification, the filter string `chemical_formula_anonymous="ABC"` is tokenized and parsed into an abstract tree by a `FilterParser` object using the Lark parsing library [@Lark] against the formal grammar defined by the specification.
2. The abstract tree is then transformed by a `FilterTransformer` object into a database query specific to the configured back-end for the server.
This transformation can include aliasing and custom transformations such that the underlying database format can be accommodated.
3. The results from the database query are then de-serialized by `EntryResourceMapper` objects into the OPTIMADE-defined data models and finally re-serialized into JSON before being served to the user over HTTP.

Beyond this query functionality, the package also provides:

- A fuzzy implementation validator that performs HTTP queries against remote or local OPTIMADE APIs, with test queries and expected responses generated dynamically based on the data served at the introspective `/info/` endpoint of the API implementation.
- Entry "adapters" that can convert between OPTIMADE-compliant entries and the data models of popular Python libraries used widely in the materials science community: `pymatgen` [@pymatgen], ASE [@ASE], AiiDA [@AiiDA], and JARVIS [@JARVIS].

# Use cases

The package is currently used in production by three major data providers for materials science data:

- The Materials Project [@MaterialsProject] uses `optimade-python-tools` alongside their existing API [@MAPI] and MongoDB database, providing access to highly-curated density-functional theory calculations across all known inorganic materials.
`optimade-python-tools` handles filter parsing, database query generation and response validation by running the reference server implementation with minimal configuration.
- NOMAD [@nomad] uses `optimade-python-tools` as a library to extend its existing web app with OPTIMADE API routes.
 It uses the Elasticsearch implementation to filter millions of structures from published first-principles calculations provided by users and other projects.
NOMAD also uses the filtering module in its own API to expose the OPTIMADE filter language in the user-centric web interface search bar.
NOMAD uses a released version of `optimade-python-tools` and all necessary customization can be realized via configuration and sub-classing.
- Materials Cloud [@MaterialsCloud] uses `optimade-python-tools` as a library to provide an OPTIMADE API entry to archived computational materials studies, created with the AiiDA [@AiiDA] Python framework and published through their archive.
In this case, each individual study and archive entry has its own database and separate API entry.
The Python classes within the `optimade` package have been extended to make use of AiiDA and its underlying [PostgreSQL](https://postgresql.org) storage engine.
- The `optimade.adapters` module from the `optimade-python-tools` library is used in a graphical web client hosted on Materials Cloud [@MaterialsCloudClient].
It allows users to query OPTIMADE API implementations using user-friendly widgets as well as raw filter strings.
The client uses the registry of known OPTIMADE providers to allow easy switching between databases.
The crystal structures returned can be inspected visually and either downloaded in formats provided by conversion functions in the `optimade.adapters` module, or used seamlessly within other Materials Cloud web tools, where the structure is automatically validated and transferred in the background, partly using the `optimade.adapters` module.

# Acknowledgements

All authors acknowledge contributions and feedback from other members of the OPTIMADE consortium, with special thanks to Michael Wu, Leopold Talirz, Thomas Purcell, Abhijith Gopakumar, Andrius Merkys and Fawzi Mohamed for their direct contributions to the `optimade` package.
M.L.E. would like to acknowledge the EPSRC Centre for Doctoral Training in Computational Methods for Materials Science for funding under grant number EP/L015552/1 and support from the European Union's Horizon 2020 research and innovation program under the European Union's Grant agreement No. 951786 (NOMAD CoE).
C.W.A. acknowledges financial support by the MARKETPLACE project, which is funded by Horizon 2020 under H2020-NMBP-25-2017 call with Grant agreement number: 760173 as well as the National Centres of Competence in Research (NCCR) Materials' revolution: Computational Design and Discovery of Novel Materials (MARVEL) created by the Swiss National Science Foundation (SNSF).
S.D. acknowledges financial support by the U.S. Department of Energy, Office of Science, Office of Basic Energy Sciences, Materials Sciences and Engineering Division under Contract No. DE-AC02-05-CH11231 (Materials Project program KC23MP).
M.S. acknowledges support from the European Union's Horizon 2020 research and innovation program under the European Union's Grant agreement No. 676580 (NoMaD) and No. 951786 (NOMAD CoE) as well as financial support from the Max Planck research network on big-data-driven materials science (BiGmax).
A.F. acknowledges support from the Communauté française de Belgique through the SURFASCOPE project (ARC 19/24-102).
# OPTIMADE Data Models

This page provides documentation for the `optimade.models` submodule, where all the OPTIMADE (and JSON:API)-defined data models are located.

For example, the three OPTIMADE entry types, `structures`, `references` and `links`, are defined primarily through the corresponding attribute models:

- [`StructureResourceAttributes`](#optimade.models.structures.StructureResourceAttributes)
- [`ReferenceResourceAttributes`](#optimade.models.references.ReferenceResourceAttributes)
- [`LinksResourceAttributes`](#optimade.models.links.LinksResourceAttributes)

As well as validating data types when creating instances of these models, this package defines several OPTIMADE-specific validators that ensure consistency between fields (e.g., the value of `nsites` matches the number of positions provided in `cartesian_site_positions`).

::: optimade.models
<img width="10%" align="left" src="images/optimade_logo_180x180.svg">

# OPTIMADE Python tools

| Latest release | Build status | Activity |
|:--------------:|:------------:|:--------:|
| [![PyPI Version](https://img.shields.io/pypi/v/optimade?logo=pypi&logoColor=white)](https://pypi.org/project/optimade/)<br>[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/optimade?logo=python&logoColor=white)](https://pypi.org/project/optimade/)<br>[![OPTIMADE](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/optimade-version.json)](https://github.com/Materials-Consortia/OPTIMADE/) | [![Build Status](https://img.shields.io/github/workflow/status/Materials-Consortia/optimade-python-tools/CI%20tests?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/actions?query=branch%3Amaster+)<br>[![codecov](https://img.shields.io/codecov/c/github/Materials-Consortia/optimade-python-tools?logo=codecov&logoColor=white&token=UJAtmqkZZO)](https://codecov.io/gh/Materials-Consortia/optimade-python-tools)<br>[![Heroku App Status](https://heroku-shields.herokuapp.com/optimade??logo=heroku)](https://optimade.herokuapp.com) | [![Commit Activity](https://img.shields.io/github/commit-activity/m/Materials-Consortia/optimade-python-tools?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/pulse)<br>[![Last Commit](https://img.shields.io/github/last-commit/Materials-Consortia/optimade-python-tools/master?logo=github)](https://github.com/Materials-Consortia/optimade-python-tools/commits/master)<br>[![Contributors](https://badgen.net/github/contributors/Materials-Consortia/optimade-python-tools?icon=github)](https://github.com/Materials-Consortia/optimade-python-tools/graphs/contributors) |

The aim of OPTIMADE is to develop a common API, compliant with the [JSON:API 1.0](http://jsonapi.org/format/1.0/) specification.
This is to enable interoperability among databases that contain calculated properties of existing and hypothetical materials.

This repository contains a library of tools for implementing and consuming [OPTIMADE](https://www.optimade.org) APIs using Python.
Server implementations can make use of the supported MongoDB (v4) and Elasticsearch (v6) database backends, or plug in a custom backend implementation.
The package also contains a server validator tool, which may be called from the shell (`optimade-validator`) or used as a GitHub Action from [optimade-validator-action](https://github.com/Materials-Consortia/optimade-validator-action).

The release history and changelog can be found in [the changelog](CHANGELOG.md).

## Documentation

This document, guides, and the full module API documentation can be found online at [https://optimade.org/optimade-python-tools](https://optimade.org/optimade-python-tools).
In particular, documentation of the OPTIMADE API response data models (implemented here with [pydantic](https://github.com/samuelcolvin/pydantic)) can be found online under [OPTIMADE Data Models](https://optimade.org/optimade-python-tools/all_models).

## Installation

Detailed installation instructions for different use cases (e.g., using the library or running a server) can be found in [the installation documentation](INSTALL.md).

The latest stable version of this package can be obtained from [PyPI](https://pypi.org/project/optimade) `pip install optimade`.
The latest development version of this package can be installed from the master branch of this repository `git clone https://github.com/Materials-Consortia/optimade-python-tools`.

## Supported OPTIMADE versions

Each release of the `optimade` package from this repository only targets one version of the OPTIMADE specification, summarised in the table below.

| OPTIMADE API version | `optimade` version |
|:--------------------:|:------------------:|
| [v1.0.0](https://github.com/Materials-Consortia/OPTIMADE/blob/v1.0.0/optimade.rst) | [v0.12.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.9) |
| [v1.1.0](https://github.com/Materials-Consortia/OPTIMADE/blob/v1.1.0/optimade.rst) | [v0.16.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.0) |

## Contributing and Getting Help

All development of this package (bug reports, suggestions, feedback and pull requests) occurs in the [optimade-python-tools GitHub repository](https://github.com/Materials-Consortia/optimade-python-tools).
Contribution guidelines and tips for getting help can be found in the [contributing notes](CONTRIBUTING.md).

## How to cite

If you use this package to access or serve OPTIMADE data, we kindly request that you consider citing the following:

- Andersen *et al.*, OPTIMADE, an API for exchanging materials data, *Sci. Data* **8**, 217 (2021) [10.1038/s41597-021-00974-z](https://doi.org/10.1038/s41597-021-00974-z)
- Evans *et al.*, optimade-python-tools: a Python library for serving and consuming materials data via OPTIMADE APIs. *Journal of Open Source Software*, **6**(65), 3458 (2021) [10.21105/joss.03458](https://doi.org/10.21105/joss.03458)

## Links

- [OPTIMADE Specification](https://github.com/Materials-Consortia/OPTIMADE/blob/develop/optimade.rst), the human-readable specification that this library is based on.
- [optimade-validator-action](https://github.com/Materials-Consortia/optimade-validator-action), a GitHub action that can be used to validate implementations from a URL (using the validator from this repo).
- [OpenAPI](https://github.com/OAI/OpenAPI-Specification), the machine-readable format used to specify the OPTIMADE API in [`openapi.json`](openapi/openapi.json) and [`index_openapi.json`](openapi/index_openapi.json).
- [Interactive documentation](https://petstore.swagger.io/?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/openapi/openapi.json) generated from [`openapi.json`](openapi/openapi.json) (see also [interactive JSON editor](https://editor.swagger.io/?url=https://raw.githubusercontent.com/Materials-Consortia/optimade-python-tools/master/openapi/openapi.json)).
- [pydantic](https://pydantic-docs.helpmanual.io/), the library used for generating the OpenAPI schema from [Python models](https://www.optimade.org/optimade-python-tools/all_models/).
- [FastAPI](https://fastapi.tiangolo.com/), the framework used for generating the reference implementation expressed by the [`openapi.json`](openapi/openapi.json) specification.
- [lark](https://github.com/lark-parser/lark), the library used to parse the filter language in OPTIMADE queries.
# Changelog

## [v0.16.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.9) (2022-01-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.8...v0.16.9)

**Implemented enhancements:**

- Lower validator default read timeout and allow it to be customised [\#1051](https://github.com/Materials-Consortia/optimade-python-tools/pull/1051) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Dependabot not updating NumPy to 1.22 [\#1035](https://github.com/Materials-Consortia/optimade-python-tools/issues/1035)

**Security fixes:**

- Elastic search log4j vulnerability [\#1040](https://github.com/Materials-Consortia/optimade-python-tools/issues/1040)

**Closed issues:**

- Remove multiple "Update dependencies" entries in CHANGELOG generation [\#1038](https://github.com/Materials-Consortia/optimade-python-tools/issues/1038)
- Docs reference to `LarkParser` failing. [\#1037](https://github.com/Materials-Consortia/optimade-python-tools/issues/1037)

**Merged pull requests:**

- Attempt to fix syntax for actions workflow [\#1053](https://github.com/Materials-Consortia/optimade-python-tools/pull/1053) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1049](https://github.com/Materials-Consortia/optimade-python-tools/pull/1049) ([CasperWA](https://github.com/CasperWA))
- Update dependabot config and changelog generation [\#1048](https://github.com/Materials-Consortia/optimade-python-tools/pull/1048) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1044](https://github.com/Materials-Consortia/optimade-python-tools/pull/1044) ([CasperWA](https://github.com/CasperWA))
- Bump elasticsearch image version to avoid any log4j issues [\#1041](https://github.com/Materials-Consortia/optimade-python-tools/pull/1041) ([ml-evs](https://github.com/ml-evs))
- Make NumPy requirement py version-specific [\#1036](https://github.com/Materials-Consortia/optimade-python-tools/pull/1036) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1031](https://github.com/Materials-Consortia/optimade-python-tools/pull/1031) ([CasperWA](https://github.com/CasperWA))

## [v0.16.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.8) (2021-12-22)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.7...v0.16.8)

**Implemented enhancements:**

- Support for Python 3.10 [\#956](https://github.com/Materials-Consortia/optimade-python-tools/issues/956)

**Fixed bugs:**

- Overzealous validation of substring comparisons for chemical formula fields [\#1024](https://github.com/Materials-Consortia/optimade-python-tools/issues/1024)

**Closed issues:**

- Updating to pymatgen v2022+ [\#762](https://github.com/Materials-Consortia/optimade-python-tools/issues/762)

**Merged pull requests:**

- Update dependencies [\#1028](https://github.com/Materials-Consortia/optimade-python-tools/pull/1028) ([CasperWA](https://github.com/CasperWA))
- Add configurable field-specific validator overrides to set filter operators as optional [\#1025](https://github.com/Materials-Consortia/optimade-python-tools/pull/1025) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#1023](https://github.com/Materials-Consortia/optimade-python-tools/pull/1023) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1017](https://github.com/Materials-Consortia/optimade-python-tools/pull/1017) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1008](https://github.com/Materials-Consortia/optimade-python-tools/pull/1008) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#1004](https://github.com/Materials-Consortia/optimade-python-tools/pull/1004) ([CasperWA](https://github.com/CasperWA))
- Add Python 3.10 support [\#957](https://github.com/Materials-Consortia/optimade-python-tools/pull/957) ([ml-evs](https://github.com/ml-evs))

## [v0.16.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.7) (2021-11-21)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.6...v0.16.7)

**Implemented enhancements:**

- Automate dependency workflow further [\#958](https://github.com/Materials-Consortia/optimade-python-tools/issues/958)
- Stricter validation of chemical formulas in OpenAPI schema [\#708](https://github.com/Materials-Consortia/optimade-python-tools/issues/708)

**Fixed bugs:**

- `chemical_formula_anonymous` validator accepts incorrect proportion order if started with 1 [\#1002](https://github.com/Materials-Consortia/optimade-python-tools/issues/1002)
- Reinstate `typing-extensions` [\#999](https://github.com/Materials-Consortia/optimade-python-tools/issues/999)
- Updating permanent dependabot branch not working after updating dependencies [\#995](https://github.com/Materials-Consortia/optimade-python-tools/issues/995)
- Auto-merge dependabot PR-workflow not running [\#984](https://github.com/Materials-Consortia/optimade-python-tools/issues/984)

**Closed issues:**

- Update the auto-PR description for updating deps [\#988](https://github.com/Materials-Consortia/optimade-python-tools/issues/988)
- Versioned docs do not redirect all links correctly [\#977](https://github.com/Materials-Consortia/optimade-python-tools/issues/977)
- Missing support for timestamps/datetime in grammar [\#102](https://github.com/Materials-Consortia/optimade-python-tools/issues/102)

**Merged pull requests:**

- Fixed bug in check\_anonymous\_formula which caused `chemical_formula_anonymous = AB2` to pass validation. [\#1001](https://github.com/Materials-Consortia/optimade-python-tools/pull/1001) ([JPBergsma](https://github.com/JPBergsma))
- Use `diff` for checking PR body [\#1000](https://github.com/Materials-Consortia/optimade-python-tools/pull/1000) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#998](https://github.com/Materials-Consortia/optimade-python-tools/pull/998) ([CasperWA](https://github.com/CasperWA))
- Correct PR body comparison [\#996](https://github.com/Materials-Consortia/optimade-python-tools/pull/996) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#993](https://github.com/Materials-Consortia/optimade-python-tools/pull/993) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#991](https://github.com/Materials-Consortia/optimade-python-tools/pull/991) ([CasperWA](https://github.com/CasperWA))
- Update dependency auto-PR message [\#989](https://github.com/Materials-Consortia/optimade-python-tools/pull/989) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#987](https://github.com/Materials-Consortia/optimade-python-tools/pull/987) ([CasperWA](https://github.com/CasperWA))
- Stricter formula syntax [\#986](https://github.com/Materials-Consortia/optimade-python-tools/pull/986) ([merkys](https://github.com/merkys))
- Implement workflows for dependency updates [\#979](https://github.com/Materials-Consortia/optimade-python-tools/pull/979) ([CasperWA](https://github.com/CasperWA))
- Tidy up old grammars, add a development grammar for v1.2 and update filterparser tests [\#879](https://github.com/Materials-Consortia/optimade-python-tools/pull/879) ([ml-evs](https://github.com/ml-evs))

## [v0.16.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.6) (2021-10-19)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.5...v0.16.6)

**Merged pull requests:**

- Put docs release deployment in separate job [\#978](https://github.com/Materials-Consortia/optimade-python-tools/pull/978) ([CasperWA](https://github.com/CasperWA))

## [v0.16.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.5) (2021-10-18)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.4...v0.16.5)

**Closed issues:**

- 'elements\_ratios' model validator uses double-precision machine epsilon - could be relaxed [\#947](https://github.com/Materials-Consortia/optimade-python-tools/issues/947)
- Versioning in Docs [\#724](https://github.com/Materials-Consortia/optimade-python-tools/issues/724)

**Merged pull requests:**

- Update dependencies [\#974](https://github.com/Materials-Consortia/optimade-python-tools/pull/974) ([CasperWA](https://github.com/CasperWA))
- Fix option value for checkout in CD Docs workflow [\#972](https://github.com/Materials-Consortia/optimade-python-tools/pull/972) ([CasperWA](https://github.com/CasperWA))
- Correct default branch name to `master` [\#971](https://github.com/Materials-Consortia/optimade-python-tools/pull/971) ([CasperWA](https://github.com/CasperWA))
- Dependabot updates for v0.16.5 [\#964](https://github.com/Materials-Consortia/optimade-python-tools/pull/964) ([ml-evs](https://github.com/ml-evs))
- Automate versioned documentation [\#951](https://github.com/Materials-Consortia/optimade-python-tools/pull/951) ([CasperWA](https://github.com/CasperWA))
- Add JOSS citation [\#949](https://github.com/Materials-Consortia/optimade-python-tools/pull/949) ([ml-evs](https://github.com/ml-evs))
- Some validation QoL tweaks [\#948](https://github.com/Materials-Consortia/optimade-python-tools/pull/948) ([ml-evs](https://github.com/ml-evs))

## [v0.16.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.4) (2021-09-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.3...v0.16.4)

**Closed issues:**

- Code check fails because there is no valid version of jsmin [\#938](https://github.com/Materials-Consortia/optimade-python-tools/issues/938)
- Be properly compliant with the new pip resolver [\#625](https://github.com/Materials-Consortia/optimade-python-tools/issues/625)

**Merged pull requests:**

- Bump providers from `357c27b` to `fb05359` [\#945](https://github.com/Materials-Consortia/optimade-python-tools/pull/945) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump providers from `368f9f6` to `357c27b` [\#944](https://github.com/Materials-Consortia/optimade-python-tools/pull/944) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump providers from `91b51bd` to `368f9f6` [\#942](https://github.com/Materials-Consortia/optimade-python-tools/pull/942) ([dependabot[bot]](https://github.com/apps/dependabot))
- remove the dependency on mkdocs-minify because of issue \#938. [\#941](https://github.com/Materials-Consortia/optimade-python-tools/pull/941) ([JPBergsma](https://github.com/JPBergsma))
- Corrected command to call uvicorn server [\#937](https://github.com/Materials-Consortia/optimade-python-tools/pull/937) ([JPBergsma](https://github.com/JPBergsma))
- Use proper pip dependency resolver in publish workflow [\#935](https://github.com/Materials-Consortia/optimade-python-tools/pull/935) ([ml-evs](https://github.com/ml-evs))
- Dependency updates for v0.16.4 [\#901](https://github.com/Materials-Consortia/optimade-python-tools/pull/901) ([ml-evs](https://github.com/ml-evs))
- Add JOSS paper [\#804](https://github.com/Materials-Consortia/optimade-python-tools/pull/804) ([ml-evs](https://github.com/ml-evs))

## [v0.16.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.3) (2021-09-02)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.2...v0.16.3)

**Implemented enhancements:**

- Add validation that anonymous/reduced chemical formulae are in fact reduced [\#913](https://github.com/Materials-Consortia/optimade-python-tools/issues/913)

**Fixed bugs:**

- No error/warning when specifying a config file that does not exist [\#930](https://github.com/Materials-Consortia/optimade-python-tools/issues/930)
- Docker tests failing in CI: http://gh\_actions\_host no longer exists? [\#906](https://github.com/Materials-Consortia/optimade-python-tools/issues/906)
- Fix config file warnings when file is missing [\#931](https://github.com/Materials-Consortia/optimade-python-tools/pull/931) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Docs don't introduce the idea of "models" [\#910](https://github.com/Materials-Consortia/optimade-python-tools/issues/910)
- Docs don't mention anything about where to go for support [\#909](https://github.com/Materials-Consortia/optimade-python-tools/issues/909)
- `run.sh` does not appear to be available from the pip installation [\#904](https://github.com/Materials-Consortia/optimade-python-tools/issues/904)
- Missing guide for how to set up an implementation from existing database [\#176](https://github.com/Materials-Consortia/optimade-python-tools/issues/176)

**Merged pull requests:**

- Add tutorial-style guide on setting up an API [\#915](https://github.com/Materials-Consortia/optimade-python-tools/pull/915) ([ml-evs](https://github.com/ml-evs))
- Add validator to check whether anonymous and reduced formulae are reduced [\#914](https://github.com/Materials-Consortia/optimade-python-tools/pull/914) ([ml-evs](https://github.com/ml-evs))
- Clarify the "all models" documentation page [\#912](https://github.com/Materials-Consortia/optimade-python-tools/pull/912) ([ml-evs](https://github.com/ml-evs))
- Add more specific 'Getting Help' info to Contributing and README [\#911](https://github.com/Materials-Consortia/optimade-python-tools/pull/911) ([ml-evs](https://github.com/ml-evs))
- Bump Materials-Consortia/optimade-validator-action from 2.5.0 to 2.6.0 [\#907](https://github.com/Materials-Consortia/optimade-python-tools/pull/907) ([dependabot[bot]](https://github.com/apps/dependabot))
- Clarify installation methods by use-case [\#905](https://github.com/Materials-Consortia/optimade-python-tools/pull/905) ([ml-evs](https://github.com/ml-evs))
- Relax response top-level root validator [\#903](https://github.com/Materials-Consortia/optimade-python-tools/pull/903) ([CasperWA](https://github.com/CasperWA))
- Add integrated app docs, tweak other use case docs [\#883](https://github.com/Materials-Consortia/optimade-python-tools/pull/883) ([ml-evs](https://github.com/ml-evs))

## [v0.16.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.2) (2021-08-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.1...v0.16.2)

**Fixed bugs:**

- Provider fallbacks are still not working [\#896](https://github.com/Materials-Consortia/optimade-python-tools/issues/896)
- Fix provider fallbacks [\#897](https://github.com/Materials-Consortia/optimade-python-tools/pull/897) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Dependency updates for v0.16.2 [\#894](https://github.com/Materials-Consortia/optimade-python-tools/pull/894) ([ml-evs](https://github.com/ml-evs))
- Bump codecov/codecov-action from 2.0.1 to 2.0.2 [\#882](https://github.com/Materials-Consortia/optimade-python-tools/pull/882) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump codecov/codecov-action from 1.5.2 to 2.0.1 [\#878](https://github.com/Materials-Consortia/optimade-python-tools/pull/878) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.16.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.1) (2021-07-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.16.0...v0.16.1)

**Implemented enhancements:**

- Change MIME type to application/vnd.api+json where appropriate [\#875](https://github.com/Materials-Consortia/optimade-python-tools/issues/875)
- Minor corrections + use model aliases for `handle_response_fields()` [\#876](https://github.com/Materials-Consortia/optimade-python-tools/pull/876) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Wrong behaviour HAS ONLY query for MongoDB [\#810](https://github.com/Materials-Consortia/optimade-python-tools/issues/810)
- Correct the behaviour of HAS ONLY with MongoDB backend [\#861](https://github.com/Materials-Consortia/optimade-python-tools/pull/861) ([JPBergsma](https://github.com/JPBergsma))

**Merged pull requests:**

- Change default MIME type to "application/vnd.api+json" [\#877](https://github.com/Materials-Consortia/optimade-python-tools/pull/877) ([ml-evs](https://github.com/ml-evs))
- Update elements description to match specification [\#874](https://github.com/Materials-Consortia/optimade-python-tools/pull/874) ([ml-evs](https://github.com/ml-evs))

## [v0.16.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.16.0) (2021-07-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.5...v0.16.0)

**Closed issues:**

- Incoming model update \(new field: issue\_tracker\) [\#592](https://github.com/Materials-Consortia/optimade-python-tools/issues/592)

**Merged pull requests:**

- Add issue\_tracker field to provider model [\#593](https://github.com/Materials-Consortia/optimade-python-tools/pull/593) ([ml-evs](https://github.com/ml-evs))

## [v0.15.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.5) (2021-07-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.4...v0.15.5)

**Fixed bugs:**

- NOT filter operation of mongo query for complex expressions [\#79](https://github.com/Materials-Consortia/optimade-python-tools/issues/79)

**Closed issues:**

- Remove CI psycopg2-binary install when aiida-core\>1.6.3 [\#855](https://github.com/Materials-Consortia/optimade-python-tools/issues/855)
- Pytest fails at Setup environment for AiiDA [\#853](https://github.com/Materials-Consortia/optimade-python-tools/issues/853)
- Add timeout parameter to validator [\#681](https://github.com/Materials-Consortia/optimade-python-tools/issues/681)
- Add note in installation instructions about pulling submodule for providers [\#370](https://github.com/Materials-Consortia/optimade-python-tools/issues/370)

**Merged pull requests:**

- Update dependencies [\#872](https://github.com/Materials-Consortia/optimade-python-tools/pull/872) ([ml-evs](https://github.com/ml-evs))
- Add request --timeout parameter to validator [\#860](https://github.com/Materials-Consortia/optimade-python-tools/pull/860) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `fa25ed3` to `91b51bd` [\#858](https://github.com/Materials-Consortia/optimade-python-tools/pull/858) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update to AiiDA v1.6.4 and remove CI fix [\#857](https://github.com/Materials-Consortia/optimade-python-tools/pull/857) ([CasperWA](https://github.com/CasperWA))
- Temporary fix for CI tests with AiiDA [\#854](https://github.com/Materials-Consortia/optimade-python-tools/pull/854) ([CasperWA](https://github.com/CasperWA))
- Documentation tweaks [\#852](https://github.com/Materials-Consortia/optimade-python-tools/pull/852) ([JPBergsma](https://github.com/JPBergsma))
- Fix query negation in MongoDB [\#814](https://github.com/Materials-Consortia/optimade-python-tools/pull/814) ([JPBergsma](https://github.com/JPBergsma))

## [v0.15.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.4) (2021-06-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.3...v0.15.4)

**Implemented enhancements:**

- Missing documentation for new configuration methods [\#766](https://github.com/Materials-Consortia/optimade-python-tools/issues/766)

**Closed issues:**

- Add docs "use case" for the validator [\#841](https://github.com/Materials-Consortia/optimade-python-tools/issues/841)
- Use specific configuration file for Heroku deployment [\#738](https://github.com/Materials-Consortia/optimade-python-tools/issues/738)
- Potential submission to JOSS? [\#203](https://github.com/Materials-Consortia/optimade-python-tools/issues/203)
- Add more tests [\#104](https://github.com/Materials-Consortia/optimade-python-tools/issues/104)

**Merged pull requests:**

- Tweak configuration docs [\#851](https://github.com/Materials-Consortia/optimade-python-tools/pull/851) ([ml-evs](https://github.com/ml-evs))
- Add some more tutorial-style documentation [\#850](https://github.com/Materials-Consortia/optimade-python-tools/pull/850) ([ml-evs](https://github.com/ml-evs))
- Bump FastAPI version in setup.py [\#849](https://github.com/Materials-Consortia/optimade-python-tools/pull/849) ([CasperWA](https://github.com/CasperWA))
- Bump fastapi from 0.65.1 to 0.65.2 [\#848](https://github.com/Materials-Consortia/optimade-python-tools/pull/848) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.15.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.3) (2021-06-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.2...v0.15.3)

**Merged pull requests:**

- Update model descriptions following spec updates [\#847](https://github.com/Materials-Consortia/optimade-python-tools/pull/847) ([ml-evs](https://github.com/ml-evs))

## [v0.15.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.2) (2021-06-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.1...v0.15.2)

**Implemented enhancements:**

- Missing HTTP response codes in OpenAPI schema [\#763](https://github.com/Materials-Consortia/optimade-python-tools/issues/763)

**Merged pull requests:**

- Update response model information for routes [\#846](https://github.com/Materials-Consortia/optimade-python-tools/pull/846) ([CasperWA](https://github.com/CasperWA))
- Improve semver validation error messsage [\#845](https://github.com/Materials-Consortia/optimade-python-tools/pull/845) ([ml-evs](https://github.com/ml-evs))
- Bump codecov/codecov-action from 1.5.0 to 1.5.2 [\#843](https://github.com/Materials-Consortia/optimade-python-tools/pull/843) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.15.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.1) (2021-06-08)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.15.0...v0.15.1)

**Closed issues:**

- mongomock $size queries match all non-array fields for {$size: 1}, even nulls [\#807](https://github.com/Materials-Consortia/optimade-python-tools/issues/807)
- Allow custom headers to be specified for validation [\#790](https://github.com/Materials-Consortia/optimade-python-tools/issues/790)

**Merged pull requests:**

- Allow both Jinja2 v2 and v3 [\#838](https://github.com/Materials-Consortia/optimade-python-tools/pull/838) ([CasperWA](https://github.com/CasperWA))
- Update mongomock and remove test skip [\#836](https://github.com/Materials-Consortia/optimade-python-tools/pull/836) ([ml-evs](https://github.com/ml-evs))
- Add --headers argument to validator to allow passing e.g. API keys [\#806](https://github.com/Materials-Consortia/optimade-python-tools/pull/806) ([ml-evs](https://github.com/ml-evs))

## [v0.15.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.15.0) (2021-06-01)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.14.1...v0.15.0)

**Fixed bugs:**

- Provider fallbacks do not get used [\#829](https://github.com/Materials-Consortia/optimade-python-tools/issues/829)
- ParserError's should not return 500 HTTP status codes [\#812](https://github.com/Materials-Consortia/optimade-python-tools/issues/812)
- Fix provider fallback list [\#830](https://github.com/Materials-Consortia/optimade-python-tools/pull/830) ([ml-evs](https://github.com/ml-evs))
- Return 400 Bad Request \(not 500\) on filter parser errors, plus filterparser module facelift [\#813](https://github.com/Materials-Consortia/optimade-python-tools/pull/813) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- CI triggered by PRs does not test dep versions in setup.py [\#834](https://github.com/Materials-Consortia/optimade-python-tools/issues/834)
- Remove Django support for v0.15+ [\#832](https://github.com/Materials-Consortia/optimade-python-tools/issues/832)
- Move aliasing code to base transformer [\#743](https://github.com/Materials-Consortia/optimade-python-tools/issues/743)
- Missing optional fields are not returned as null when requested with response\_fields [\#516](https://github.com/Materials-Consortia/optimade-python-tools/issues/516)

**Merged pull requests:**

- Test all setup.py deps versions for every pull request, plus some deps updates [\#835](https://github.com/Materials-Consortia/optimade-python-tools/pull/835) ([ml-evs](https://github.com/ml-evs))
- Deprecate Python 3.6, remove Django and update dependencies/providers [\#828](https://github.com/Materials-Consortia/optimade-python-tools/pull/828) ([ml-evs](https://github.com/ml-evs))
- Update INSTALL docs [\#811](https://github.com/Materials-Consortia/optimade-python-tools/pull/811) ([ml-evs](https://github.com/ml-evs))
- Overhaul of filter transformers, mappers and response fields [\#797](https://github.com/Materials-Consortia/optimade-python-tools/pull/797) ([ml-evs](https://github.com/ml-evs))

## [v0.14.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.14.1) (2021-05-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.14.0...v0.14.1)

**Fixed bugs:**

- \[SECURITY\] Cycle secrets [\#777](https://github.com/Materials-Consortia/optimade-python-tools/issues/777)

**Closed issues:**

- Do not validate extension endpoints [\#793](https://github.com/Materials-Consortia/optimade-python-tools/issues/793)
- Verify that missing values are not returned in comparisons [\#792](https://github.com/Materials-Consortia/optimade-python-tools/issues/792)

**Merged pull requests:**

- Bump pydantic from 1.8.1 to 1.8.2 [\#805](https://github.com/Materials-Consortia/optimade-python-tools/pull/805) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update GH actions [\#803](https://github.com/Materials-Consortia/optimade-python-tools/pull/803) ([CasperWA](https://github.com/CasperWA))
- Handling null fields in the filtertransformer and validator [\#796](https://github.com/Materials-Consortia/optimade-python-tools/pull/796) ([ml-evs](https://github.com/ml-evs))
- Filter out extension endpoints before validation [\#794](https://github.com/Materials-Consortia/optimade-python-tools/pull/794) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `7a54843` to `fa25ed3` [\#791](https://github.com/Materials-Consortia/optimade-python-tools/pull/791) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump typing-extensions from 3.7.4.3 to 3.10.0.0 [\#789](https://github.com/Materials-Consortia/optimade-python-tools/pull/789) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#787](https://github.com/Materials-Consortia/optimade-python-tools/pull/787) ([CasperWA](https://github.com/CasperWA))
- Bump CharMixer/auto-changelog-action from v1.2 to v1.3 [\#778](https://github.com/Materials-Consortia/optimade-python-tools/pull/778) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump django from 3.1.7 to 3.1.8 [\#776](https://github.com/Materials-Consortia/optimade-python-tools/pull/776) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#773](https://github.com/Materials-Consortia/optimade-python-tools/pull/773) ([ml-evs](https://github.com/ml-evs))

## [v0.14.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.14.0) (2021-03-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.3...v0.14.0)

**Implemented enhancements:**

- Rename config variable use\_real\_mongo to something more general [\#742](https://github.com/Materials-Consortia/optimade-python-tools/issues/742)
- Custom configuration extensions & use standard pydantic way of loading config file [\#739](https://github.com/Materials-Consortia/optimade-python-tools/issues/739)
- Generalising collections and adding ElasticsearchCollection [\#660](https://github.com/Materials-Consortia/optimade-python-tools/pull/660) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Over-aggressive middleware to check versioned base URL [\#737](https://github.com/Materials-Consortia/optimade-python-tools/issues/737)
- Floating point comparisons should not be tested with the validator [\#735](https://github.com/Materials-Consortia/optimade-python-tools/issues/735)
- Mapper method `alias_of` extracts alias wrongly [\#667](https://github.com/Materials-Consortia/optimade-python-tools/issues/667)

**Closed issues:**

- Docs builds are not properly tested for each PR [\#747](https://github.com/Materials-Consortia/optimade-python-tools/issues/747)
- Remove SQLAlchemy version fix in CI with new AiiDA version [\#745](https://github.com/Materials-Consortia/optimade-python-tools/issues/745)

**Merged pull requests:**

- Update dependencies [\#760](https://github.com/Materials-Consortia/optimade-python-tools/pull/760) ([CasperWA](https://github.com/CasperWA))
- Fix CheckWronglyVersionedBaseUrls middleware \(for landing pages\) [\#752](https://github.com/Materials-Consortia/optimade-python-tools/pull/752) ([CasperWA](https://github.com/CasperWA))
- Deprecate Python 3.6 support, v0.14 last supported version [\#751](https://github.com/Materials-Consortia/optimade-python-tools/pull/751) ([CasperWA](https://github.com/CasperWA))
- Run full API docs invoke task for every PR [\#748](https://github.com/Materials-Consortia/optimade-python-tools/pull/748) ([ml-evs](https://github.com/ml-evs))
- Change aliasing method names in mapper and deprecate the old [\#746](https://github.com/Materials-Consortia/optimade-python-tools/pull/746) ([ml-evs](https://github.com/ml-evs))
- Bump providers from `e2074e8` to `7a54843` [\#741](https://github.com/Materials-Consortia/optimade-python-tools/pull/741) ([dependabot[bot]](https://github.com/apps/dependabot))
- Config updates [\#740](https://github.com/Materials-Consortia/optimade-python-tools/pull/740) ([CasperWA](https://github.com/CasperWA))
- Disable all floating-point comparisons during validation [\#736](https://github.com/Materials-Consortia/optimade-python-tools/pull/736) ([ml-evs](https://github.com/ml-evs))
- Report user errors in filter as HTTP 400 Bad Request and not 501 Not Implemented [\#658](https://github.com/Materials-Consortia/optimade-python-tools/pull/658) ([markus1978](https://github.com/markus1978))

## [v0.13.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.3) (2021-03-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.2...v0.13.3)

**Fixed bugs:**

- Support `anyOf`, `allOf`, etc. standard OpenAPI fields [\#730](https://github.com/Materials-Consortia/optimade-python-tools/issues/730)
- Python 3.9 support invalid [\#728](https://github.com/Materials-Consortia/optimade-python-tools/issues/728)

**Merged pull requests:**

- Update dependencies [\#734](https://github.com/Materials-Consortia/optimade-python-tools/pull/734) ([CasperWA](https://github.com/CasperWA))
- Update pydantic to ~=1.8 [\#731](https://github.com/Materials-Consortia/optimade-python-tools/pull/731) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `da74513` to `e2074e8` [\#727](https://github.com/Materials-Consortia/optimade-python-tools/pull/727) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.13.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.2) (2021-03-01)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.1...v0.13.2)

**Implemented enhancements:**

- Improve validation of providers [\#723](https://github.com/Materials-Consortia/optimade-python-tools/issues/723)

**Merged pull requests:**

- Update dependencies [\#725](https://github.com/Materials-Consortia/optimade-python-tools/pull/725) ([CasperWA](https://github.com/CasperWA))

## [v0.13.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.1) (2021-02-23)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.13.0...v0.13.1)

**Fixed bugs:**

- Supported OPTIMADE \_\_api\_version\_\_ is incorrect in latest release [\#712](https://github.com/Materials-Consortia/optimade-python-tools/issues/712)

**Merged pull requests:**

- Bump OPTIMADE version [\#713](https://github.com/Materials-Consortia/optimade-python-tools/pull/713) ([ml-evs](https://github.com/ml-evs))

## [v0.13.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.13.0) (2021-02-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.9...v0.13.0)

**Closed issues:**

- Update species.mass model [\#630](https://github.com/Materials-Consortia/optimade-python-tools/issues/630)

**Merged pull requests:**

- Update species-\>mass field following specification change [\#631](https://github.com/Materials-Consortia/optimade-python-tools/pull/631) ([ml-evs](https://github.com/ml-evs))

## [v0.12.9](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.9) (2021-02-10)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.8...v0.12.9)

**Implemented enhancements:**

- Improve support for timestamp queries in MongoTransformer [\#590](https://github.com/Materials-Consortia/optimade-python-tools/pull/590) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Use Enums for pydantic model defaults instead of strings [\#683](https://github.com/Materials-Consortia/optimade-python-tools/issues/683)

**Closed issues:**

- When using `--as-type` in validator, one does not get a summary \(`--json` doesn't work\) [\#699](https://github.com/Materials-Consortia/optimade-python-tools/issues/699)
- Extension/import issue with mongo collection [\#682](https://github.com/Materials-Consortia/optimade-python-tools/issues/682)

**Merged pull requests:**

- Update dependencies [\#707](https://github.com/Materials-Consortia/optimade-python-tools/pull/707) ([CasperWA](https://github.com/CasperWA))
- Always print summary as last thing in validation [\#700](https://github.com/Materials-Consortia/optimade-python-tools/pull/700) ([CasperWA](https://github.com/CasperWA))
- Bump django from 3.1.5 to 3.1.6 [\#698](https://github.com/Materials-Consortia/optimade-python-tools/pull/698) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update dependencies [\#697](https://github.com/Materials-Consortia/optimade-python-tools/pull/697) ([CasperWA](https://github.com/CasperWA))
- Fixes for new gateway implementation [\#684](https://github.com/Materials-Consortia/optimade-python-tools/pull/684) ([CasperWA](https://github.com/CasperWA))

## [v0.12.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.8) (2021-01-18)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.7...v0.12.8)

**Implemented enhancements:**

- Validate mandatory query field `structure_features` [\#678](https://github.com/Materials-Consortia/optimade-python-tools/issues/678)

**Fixed bugs:**

- Validator should not rely on `meta->data_available` [\#677](https://github.com/Materials-Consortia/optimade-python-tools/issues/677)
- Validator should not rely on SHOULD "meta" field "data\_returned" [\#675](https://github.com/Materials-Consortia/optimade-python-tools/issues/675)
- Validator: remove reliance on meta fields and check mandatory queries [\#676](https://github.com/Materials-Consortia/optimade-python-tools/pull/676) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Bump providers from `542ac0a` to `da74513` [\#679](https://github.com/Materials-Consortia/optimade-python-tools/pull/679) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.12.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.7) (2021-01-15)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.6...v0.12.7)

**Implemented enhancements:**

- Make content-type response checks on '/versions` endpoint optional [\#670](https://github.com/Materials-Consortia/optimade-python-tools/pull/670) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Publish workflow fails when no changes to api docs between versions [\#673](https://github.com/Materials-Consortia/optimade-python-tools/issues/673)
- /versions header `Content-Type` value should be granularized according to RFC requirements in validator [\#669](https://github.com/Materials-Consortia/optimade-python-tools/issues/669)
- Misleading error message from validator on failure from '/versions' [\#668](https://github.com/Materials-Consortia/optimade-python-tools/issues/668)
- Fix publishing workflow [\#674](https://github.com/Materials-Consortia/optimade-python-tools/pull/674) ([ml-evs](https://github.com/ml-evs))

**Merged pull requests:**

- Update codecov coverage config file [\#672](https://github.com/Materials-Consortia/optimade-python-tools/pull/672) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `fe5048b` to `542ac0a` [\#671](https://github.com/Materials-Consortia/optimade-python-tools/pull/671) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.12.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.6) (2021-01-08)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.5...v0.12.6)

**Implemented enhancements:**

- Create base transformer [\#286](https://github.com/Materials-Consortia/optimade-python-tools/issues/286)

**Fixed bugs:**

- Our models and validator are too strict [\#399](https://github.com/Materials-Consortia/optimade-python-tools/issues/399)
- Validator changes: always check unversioned '/versions' and handle rich HTML pages [\#665](https://github.com/Materials-Consortia/optimade-python-tools/pull/665) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Add more prominent link to rendered docs [\#628](https://github.com/Materials-Consortia/optimade-python-tools/issues/628)
- Review the required properties of StructureResourceAttributes in openapi.json [\#198](https://github.com/Materials-Consortia/optimade-python-tools/issues/198)

**Merged pull requests:**

- Added GitHub CODEOWNERS [\#664](https://github.com/Materials-Consortia/optimade-python-tools/pull/664) ([ml-evs](https://github.com/ml-evs))
- Robustness improvements to validator [\#659](https://github.com/Materials-Consortia/optimade-python-tools/pull/659) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#655](https://github.com/Materials-Consortia/optimade-python-tools/pull/655) ([CasperWA](https://github.com/CasperWA))
- Bugfixes for elasticsearch filtertransformer comparision operators. [\#648](https://github.com/Materials-Consortia/optimade-python-tools/pull/648) ([markus1978](https://github.com/markus1978))
- Update dependencies [\#647](https://github.com/Materials-Consortia/optimade-python-tools/pull/647) ([ml-evs](https://github.com/ml-evs))
- Added "root\_path" config parameter for FastAPI apps [\#634](https://github.com/Materials-Consortia/optimade-python-tools/pull/634) ([markus1978](https://github.com/markus1978))
- Bump providers from `2673be6` to `fe5048b` [\#633](https://github.com/Materials-Consortia/optimade-python-tools/pull/633) ([dependabot[bot]](https://github.com/apps/dependabot))
- Updated README and moved some files to top-level [\#629](https://github.com/Materials-Consortia/optimade-python-tools/pull/629) ([ml-evs](https://github.com/ml-evs))
- insert reading of default optimade\_config.json in example run script run.sh [\#627](https://github.com/Materials-Consortia/optimade-python-tools/pull/627) ([rartino](https://github.com/rartino))
- Create template filtertransformer BaseTransformer [\#287](https://github.com/Materials-Consortia/optimade-python-tools/pull/287) ([ml-evs](https://github.com/ml-evs))

## [v0.12.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.5) (2020-12-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.4...v0.12.5)

**Closed issues:**

- PyPI publishing build is broken by latest pip [\#624](https://github.com/Materials-Consortia/optimade-python-tools/issues/624)
- Empty endpoints raise errors on validation [\#622](https://github.com/Materials-Consortia/optimade-python-tools/issues/622)
- Frequency of updating online docs [\#452](https://github.com/Materials-Consortia/optimade-python-tools/issues/452)

**Merged pull requests:**

- Fix PyPI publishing in CI [\#623](https://github.com/Materials-Consortia/optimade-python-tools/pull/623) ([ml-evs](https://github.com/ml-evs))
- Change validation error to warning on empty endpoints [\#621](https://github.com/Materials-Consortia/optimade-python-tools/pull/621) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#620](https://github.com/Materials-Consortia/optimade-python-tools/pull/620) ([CasperWA](https://github.com/CasperWA))
- Upstream fixes from specification [\#611](https://github.com/Materials-Consortia/optimade-python-tools/pull/611) ([ml-evs](https://github.com/ml-evs))
- Minor fixes for the validator [\#610](https://github.com/Materials-Consortia/optimade-python-tools/pull/610) ([ml-evs](https://github.com/ml-evs))
- Dependency updates [\#607](https://github.com/Materials-Consortia/optimade-python-tools/pull/607) ([ml-evs](https://github.com/ml-evs))
- include LICENSE in pip Package [\#594](https://github.com/Materials-Consortia/optimade-python-tools/pull/594) ([jan-janssen](https://github.com/jan-janssen))
- Relax models to allow for all SHOULD fields to be None [\#560](https://github.com/Materials-Consortia/optimade-python-tools/pull/560) ([ml-evs](https://github.com/ml-evs))
- Python 3.9 support [\#558](https://github.com/Materials-Consortia/optimade-python-tools/pull/558) ([ml-evs](https://github.com/ml-evs))
- ReadTheDocs configuration file \(v2\) [\#485](https://github.com/Materials-Consortia/optimade-python-tools/pull/485) ([CasperWA](https://github.com/CasperWA))

## [v0.12.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.4) (2020-11-16)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.3...v0.12.4)

**Merged pull requests:**

- Minor fixes for versions endpoint validation [\#591](https://github.com/Materials-Consortia/optimade-python-tools/pull/591) ([ml-evs](https://github.com/ml-evs))
- Add --minimal/--page\_limit validator options and remove old code [\#571](https://github.com/Materials-Consortia/optimade-python-tools/pull/571) ([ml-evs](https://github.com/ml-evs))

## [v0.12.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.3) (2020-11-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.2...v0.12.3)

**Fixed bugs:**

- GITHUB\_TOKEN not useful for changelog action [\#587](https://github.com/Materials-Consortia/optimade-python-tools/issues/587)
- Hill notation wrong \(still\) [\#585](https://github.com/Materials-Consortia/optimade-python-tools/issues/585)
- Hill notation validation turning around C and H [\#581](https://github.com/Materials-Consortia/optimade-python-tools/issues/581)

**Closed issues:**

- Make structure "deformity" tests more robust [\#583](https://github.com/Materials-Consortia/optimade-python-tools/issues/583)
- Incomplete output of optimade-validator [\#568](https://github.com/Materials-Consortia/optimade-python-tools/issues/568)

**Merged pull requests:**

- Use special release PAT for CHANGELOG generation action [\#588](https://github.com/Materials-Consortia/optimade-python-tools/pull/588) ([CasperWA](https://github.com/CasperWA))
- Check for carbon in elements for Hill [\#586](https://github.com/Materials-Consortia/optimade-python-tools/pull/586) ([CasperWA](https://github.com/CasperWA))
- Added better expected error messages to deformity tests [\#584](https://github.com/Materials-Consortia/optimade-python-tools/pull/584) ([ml-evs](https://github.com/ml-evs))
- Fix Hill ordering validation [\#582](https://github.com/Materials-Consortia/optimade-python-tools/pull/582) ([CasperWA](https://github.com/CasperWA))
- Bump mkdocs-material from 6.1.0 to 6.1.2 [\#580](https://github.com/Materials-Consortia/optimade-python-tools/pull/580) ([dependabot[bot]](https://github.com/apps/dependabot))
- Moved CONFIG import so it does not get triggered when just importing mapper [\#569](https://github.com/Materials-Consortia/optimade-python-tools/pull/569) ([ml-evs](https://github.com/ml-evs))

## [v0.12.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.2) (2020-10-31)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.1...v0.12.2)

**Implemented enhancements:**

- Add convenience method for adding all required middleware [\#536](https://github.com/Materials-Consortia/optimade-python-tools/issues/536)
- Add model validators and regexp for chemical formulae fields [\#547](https://github.com/Materials-Consortia/optimade-python-tools/pull/547) ([ml-evs](https://github.com/ml-evs))
- Validator improvements [\#515](https://github.com/Materials-Consortia/optimade-python-tools/pull/515) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- 'Chosen entry had no value for ...' when property is not requested [\#514](https://github.com/Materials-Consortia/optimade-python-tools/issues/514)
- Fix Species validators and error messages [\#561](https://github.com/Materials-Consortia/optimade-python-tools/pull/561) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Chemical symbols D and T [\#570](https://github.com/Materials-Consortia/optimade-python-tools/issues/570)
- Push back dependabot to monthly updates [\#567](https://github.com/Materials-Consortia/optimade-python-tools/issues/567)
- Spurious validation errors in Structure-\>Species [\#559](https://github.com/Materials-Consortia/optimade-python-tools/issues/559)
- Chemical formulae are not properly validated on model creation [\#546](https://github.com/Materials-Consortia/optimade-python-tools/issues/546)

**Merged pull requests:**

- Update dependencies [\#578](https://github.com/Materials-Consortia/optimade-python-tools/pull/578) ([CasperWA](https://github.com/CasperWA))
- Bump CasperWA/push-protected from v1 to v2.1.0 [\#573](https://github.com/Materials-Consortia/optimade-python-tools/pull/573) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update deps [\#566](https://github.com/Materials-Consortia/optimade-python-tools/pull/566) ([ml-evs](https://github.com/ml-evs))
- Improve handling of MongoDB ObjectID [\#557](https://github.com/Materials-Consortia/optimade-python-tools/pull/557) ([ml-evs](https://github.com/ml-evs))
- Update deps [\#556](https://github.com/Materials-Consortia/optimade-python-tools/pull/556) ([ml-evs](https://github.com/ml-evs))
- Updated dependencies [\#551](https://github.com/Materials-Consortia/optimade-python-tools/pull/551) ([ml-evs](https://github.com/ml-evs))
- Update dependencies - remove black as direct dependency [\#545](https://github.com/Materials-Consortia/optimade-python-tools/pull/545) ([CasperWA](https://github.com/CasperWA))
- Added convenience variables for middleware and exception handlers [\#537](https://github.com/Materials-Consortia/optimade-python-tools/pull/537) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#531](https://github.com/Materials-Consortia/optimade-python-tools/pull/531) ([ml-evs](https://github.com/ml-evs))

## [v0.12.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.1) (2020-09-24)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.12.0...v0.12.1)

**Implemented enhancements:**

- Move entry schemas to separate submodule [\#511](https://github.com/Materials-Consortia/optimade-python-tools/pull/511) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator should allow implementations to return "501 Not Implemented" for unsupported filters [\#518](https://github.com/Materials-Consortia/optimade-python-tools/issues/518)
- Landing page wrong URL  [\#371](https://github.com/Materials-Consortia/optimade-python-tools/issues/371)

**Merged pull requests:**

- This should ensure requirements\*.txt are tested [\#527](https://github.com/Materials-Consortia/optimade-python-tools/pull/527) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#526](https://github.com/Materials-Consortia/optimade-python-tools/pull/526) ([CasperWA](https://github.com/CasperWA))
- Fix landing page URL [\#519](https://github.com/Materials-Consortia/optimade-python-tools/pull/519) ([shyamd](https://github.com/shyamd))
- Update dependencies [\#510](https://github.com/Materials-Consortia/optimade-python-tools/pull/510) ([ml-evs](https://github.com/ml-evs))
- Fixing typo `validatated` -\> `validated` [\#506](https://github.com/Materials-Consortia/optimade-python-tools/pull/506) ([merkys](https://github.com/merkys))
- Make validator respond to KeyboardInterrupts [\#505](https://github.com/Materials-Consortia/optimade-python-tools/pull/505) ([ml-evs](https://github.com/ml-evs))
- Add support levels to validator config [\#503](https://github.com/Materials-Consortia/optimade-python-tools/pull/503) ([ml-evs](https://github.com/ml-evs))
- Enable JSON response from the validator [\#502](https://github.com/Materials-Consortia/optimade-python-tools/pull/502) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#501](https://github.com/Materials-Consortia/optimade-python-tools/pull/501) ([CasperWA](https://github.com/CasperWA))

## [v0.12.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.12.0) (2020-09-11)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.11.0...v0.12.0)

**Fixed bugs:**

- Missing field descriptions in schema for Species-\>name and Person-\>name [\#492](https://github.com/Materials-Consortia/optimade-python-tools/issues/492)
- "type" field not marked as required for derived entry resource models [\#479](https://github.com/Materials-Consortia/optimade-python-tools/issues/479)
- OpenAPI validations fails due to incorrect type of "dimension\_types" [\#478](https://github.com/Materials-Consortia/optimade-python-tools/issues/478)
- Have fallbacks for retrieving providers list [\#450](https://github.com/Materials-Consortia/optimade-python-tools/issues/450)
- Commit only when necessary [\#495](https://github.com/Materials-Consortia/optimade-python-tools/pull/495) ([CasperWA](https://github.com/CasperWA))
- Fix field optonality inconsistency in schema [\#482](https://github.com/Materials-Consortia/optimade-python-tools/pull/482) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator message for wrong version [\#493](https://github.com/Materials-Consortia/optimade-python-tools/issues/493)
- Validator should validate versions endpoint [\#491](https://github.com/Materials-Consortia/optimade-python-tools/issues/491)
- List of providers not included in `/links` endpoint for index meta-database [\#454](https://github.com/Materials-Consortia/optimade-python-tools/issues/454)
- Validate bad version URLs responding with 553 Version Not Supported [\#427](https://github.com/Materials-Consortia/optimade-python-tools/issues/427)
- Nonexistent property 'list' in validator tests [\#423](https://github.com/Materials-Consortia/optimade-python-tools/issues/423)
- Test `data_returned` [\#402](https://github.com/Materials-Consortia/optimade-python-tools/issues/402)
- AiiDA tests only run on Python 3.8 in CI [\#401](https://github.com/Materials-Consortia/optimade-python-tools/issues/401)
- Links under top-level 'links' may be objects [\#394](https://github.com/Materials-Consortia/optimade-python-tools/issues/394)
- Suggestion: use absolute imports in app code to allow re-use [\#298](https://github.com/Materials-Consortia/optimade-python-tools/issues/298)
- Update mongomock requirement when next released [\#207](https://github.com/Materials-Consortia/optimade-python-tools/issues/207)
- error when browsing OpenAPI docs [\#192](https://github.com/Materials-Consortia/optimade-python-tools/issues/192)

**Merged pull requests:**

- Don't report untracked and ignored files [\#496](https://github.com/Materials-Consortia/optimade-python-tools/pull/496) ([CasperWA](https://github.com/CasperWA))
- Improved error message for bad version returning 553 [\#494](https://github.com/Materials-Consortia/optimade-python-tools/pull/494) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#490](https://github.com/Materials-Consortia/optimade-python-tools/pull/490) ([CasperWA](https://github.com/CasperWA))
- Allow Link objects for pagination [\#484](https://github.com/Materials-Consortia/optimade-python-tools/pull/484) ([ml-evs](https://github.com/ml-evs))
- Absolute imports [\#483](https://github.com/Materials-Consortia/optimade-python-tools/pull/483) ([CasperWA](https://github.com/CasperWA))
- Validate OpenAPI specification in CI [\#481](https://github.com/Materials-Consortia/optimade-python-tools/pull/481) ([ml-evs](https://github.com/ml-evs))
- Update types to align with OpenAPI [\#480](https://github.com/Materials-Consortia/optimade-python-tools/pull/480) ([CasperWA](https://github.com/CasperWA))
- Update dependencies and pre-commit [\#477](https://github.com/Materials-Consortia/optimade-python-tools/pull/477) ([CasperWA](https://github.com/CasperWA))
- Unpin CI Python version for AiiDA tests [\#472](https://github.com/Materials-Consortia/optimade-python-tools/pull/472) ([ml-evs](https://github.com/ml-evs))
- Update dependencies [\#471](https://github.com/Materials-Consortia/optimade-python-tools/pull/471) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#466](https://github.com/Materials-Consortia/optimade-python-tools/pull/466) ([CasperWA](https://github.com/CasperWA))
- Provider list fallback and list of providers in both servers' `/links`-endpoints [\#455](https://github.com/Materials-Consortia/optimade-python-tools/pull/455) ([CasperWA](https://github.com/CasperWA))
- SHOULD/MUST/OPTIONAL fields in models [\#453](https://github.com/Materials-Consortia/optimade-python-tools/pull/453) ([ml-evs](https://github.com/ml-evs))
- Validator overhaul [\#417](https://github.com/Materials-Consortia/optimade-python-tools/pull/417) ([ml-evs](https://github.com/ml-evs))

## [v0.11.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.11.0) (2020-08-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.10.0...v0.11.0)

**Implemented enhancements:**

- Use logging more thoroughly throughout the code base [\#242](https://github.com/Materials-Consortia/optimade-python-tools/issues/242)
- Implement `warnings` [\#105](https://github.com/Materials-Consortia/optimade-python-tools/issues/105)

**Fixed bugs:**

- Heroku is failing - raising OSError when making LOGS\_DIR [\#448](https://github.com/Materials-Consortia/optimade-python-tools/issues/448)
- `/versions` endpoint content-type parameter "header=present" is provided in the wrong place [\#418](https://github.com/Materials-Consortia/optimade-python-tools/issues/418)
- Publish workflow cannot push to protected branch [\#341](https://github.com/Materials-Consortia/optimade-python-tools/issues/341)
- Fix circular dep and extra permission error in logs [\#436](https://github.com/Materials-Consortia/optimade-python-tools/pull/436) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- log\_dir option in config is unused [\#435](https://github.com/Materials-Consortia/optimade-python-tools/issues/435)
- Allow all types of JSON API relationships [\#429](https://github.com/Materials-Consortia/optimade-python-tools/issues/429)
- OPTIMADE version badge was not bumped on 1.0 release [\#415](https://github.com/Materials-Consortia/optimade-python-tools/issues/415)
- Add `api_hint` query parameter [\#392](https://github.com/Materials-Consortia/optimade-python-tools/issues/392)
- Return 553 for wrongly versioned base URLs [\#391](https://github.com/Materials-Consortia/optimade-python-tools/issues/391)
- Private/dunder methods incorrectly documented in mkdocs [\#365](https://github.com/Materials-Consortia/optimade-python-tools/issues/365)
- Configuration documentation [\#310](https://github.com/Materials-Consortia/optimade-python-tools/issues/310)
- Improve handling of sorting in MongoDB backend [\#276](https://github.com/Materials-Consortia/optimade-python-tools/issues/276)

**Merged pull requests:**

- Catch OSError instead of PermissionError when making log dir [\#449](https://github.com/Materials-Consortia/optimade-python-tools/pull/449) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#447](https://github.com/Materials-Consortia/optimade-python-tools/pull/447) ([CasperWA](https://github.com/CasperWA))
- Bump mkdocstrings from 0.12.1 to 0.12.2 and mkdocs-material from 5.5.0 to 5.5.2 [\#440](https://github.com/Materials-Consortia/optimade-python-tools/pull/440) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump uvicorn from 0.11.5 to 0.11.7 [\#433](https://github.com/Materials-Consortia/optimade-python-tools/pull/433) ([dependabot[bot]](https://github.com/apps/dependabot))
- Introduce logging [\#432](https://github.com/Materials-Consortia/optimade-python-tools/pull/432) ([CasperWA](https://github.com/CasperWA))
- New middleware to catch any `OptimadeWarning`s [\#431](https://github.com/Materials-Consortia/optimade-python-tools/pull/431) ([CasperWA](https://github.com/CasperWA))
- Auto-generate API reference in docs and an overhaul [\#430](https://github.com/Materials-Consortia/optimade-python-tools/pull/430) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `52027b1` to `9712dd8` [\#428](https://github.com/Materials-Consortia/optimade-python-tools/pull/428) ([dependabot[bot]](https://github.com/apps/dependabot))
- Cleanup config files [\#426](https://github.com/Materials-Consortia/optimade-python-tools/pull/426) ([CasperWA](https://github.com/CasperWA))
- Update more unittest tests to pytest [\#425](https://github.com/Materials-Consortia/optimade-python-tools/pull/425) ([CasperWA](https://github.com/CasperWA))
- Sorting on unknown properties: returning Bad Request when appropriate [\#424](https://github.com/Materials-Consortia/optimade-python-tools/pull/424) ([ml-evs](https://github.com/ml-evs))
- Minor CI updates [\#422](https://github.com/Materials-Consortia/optimade-python-tools/pull/422) ([CasperWA](https://github.com/CasperWA))
- Add `api_hint` query parameter [\#421](https://github.com/Materials-Consortia/optimade-python-tools/pull/421) ([CasperWA](https://github.com/CasperWA))
- Implement 553 Version Not Supported [\#420](https://github.com/Materials-Consortia/optimade-python-tools/pull/420) ([CasperWA](https://github.com/CasperWA))
- Fix incorrect placement of header=present in versions endpoint [\#419](https://github.com/Materials-Consortia/optimade-python-tools/pull/419) ([ml-evs](https://github.com/ml-evs))
- Bump optimade-version.json to 1.0.0 [\#416](https://github.com/Materials-Consortia/optimade-python-tools/pull/416) ([ml-evs](https://github.com/ml-evs))
- Use optimade-validator-action v2 [\#413](https://github.com/Materials-Consortia/optimade-python-tools/pull/413) ([CasperWA](https://github.com/CasperWA))
- Bump providers from `a96d424` to `52027b1` [\#389](https://github.com/Materials-Consortia/optimade-python-tools/pull/389) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.10.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.10.0) (2020-07-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.8...v0.10.0)

**Implemented enhancements:**

- Move tests to pytest system from unittest [\#270](https://github.com/Materials-Consortia/optimade-python-tools/issues/270)

**Fixed bugs:**

- Fix /vMAJOR/info in index server [\#414](https://github.com/Materials-Consortia/optimade-python-tools/pull/414) ([CasperWA](https://github.com/CasperWA))

**Closed issues:**

- Validation of 'structures' type crashes [\#397](https://github.com/Materials-Consortia/optimade-python-tools/issues/397)
- Validator verbosity levels need more detailed description [\#396](https://github.com/Materials-Consortia/optimade-python-tools/issues/396)
- Validator treats top-level 'included' array as mandatory [\#393](https://github.com/Materials-Consortia/optimade-python-tools/issues/393)
- \(Un\)versioned URLs [\#379](https://github.com/Materials-Consortia/optimade-python-tools/issues/379)

**Merged pull requests:**

- Update dependencies [\#412](https://github.com/Materials-Consortia/optimade-python-tools/pull/412) ([CasperWA](https://github.com/CasperWA))
- Bump pydantic from 1.5.1 to 1.6.1 [\#405](https://github.com/Materials-Consortia/optimade-python-tools/pull/405) ([dependabot[bot]](https://github.com/apps/dependabot))
- Temporarily run AiiDA tests on Python 3.8 only [\#400](https://github.com/Materials-Consortia/optimade-python-tools/pull/400) ([ml-evs](https://github.com/ml-evs))
- Make the example for --as\_type more similar to a real use case [\#398](https://github.com/Materials-Consortia/optimade-python-tools/pull/398) ([merkys](https://github.com/merkys))
- Fix some validator-specific crashes [\#395](https://github.com/Materials-Consortia/optimade-python-tools/pull/395) ([ml-evs](https://github.com/ml-evs))
- Use pytest instead of unittest [\#390](https://github.com/Materials-Consortia/optimade-python-tools/pull/390) ([CasperWA](https://github.com/CasperWA))
- Update dependencies [\#388](https://github.com/Materials-Consortia/optimade-python-tools/pull/388) ([CasperWA](https://github.com/CasperWA))

## [v0.9.8](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.8) (2020-07-03)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.7...v0.9.8)

**Implemented enhancements:**

- Set implementation version in config by default [\#385](https://github.com/Materials-Consortia/optimade-python-tools/pull/385) ([CasperWA](https://github.com/CasperWA))

**Merged pull requests:**

- Update models, endpoints and responses to 1.0.0 [\#380](https://github.com/Materials-Consortia/optimade-python-tools/pull/380) ([ml-evs](https://github.com/ml-evs))

## [v0.9.7](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.7) (2020-06-28)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.6...v0.9.7)

## [v0.9.6](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.6) (2020-06-28)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.5...v0.9.6)

**Fixed bugs:**

- Fix publish workflow - final\(TM\) fix [\#378](https://github.com/Materials-Consortia/optimade-python-tools/pull/378) ([CasperWA](https://github.com/CasperWA))

## [v0.9.5](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.5) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.4...v0.9.5)

**Implemented enhancements:**

- Use new action for publishing [\#377](https://github.com/Materials-Consortia/optimade-python-tools/pull/377) ([CasperWA](https://github.com/CasperWA))

## [v0.9.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.4) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.3...v0.9.4)

## [v0.9.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.3) (2020-06-26)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.2...v0.9.3)

**Merged pull requests:**

- Fix version issues in the publish workflow [\#376](https://github.com/Materials-Consortia/optimade-python-tools/pull/376) ([shyamd](https://github.com/shyamd))
- Bump providers from `732593a` to `a96d424` [\#368](https://github.com/Materials-Consortia/optimade-python-tools/pull/368) ([dependabot[bot]](https://github.com/apps/dependabot))

## [v0.9.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.2) (2020-06-25)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.1...v0.9.2)

**Fixed bugs:**

- Heroku cannot handle submodules when deploying via GitHub [\#373](https://github.com/Materials-Consortia/optimade-python-tools/issues/373)

**Closed issues:**

- Updates to models \(new OPTIONAL `type` field under `properties`\) [\#345](https://github.com/Materials-Consortia/optimade-python-tools/issues/345)
- Add aggregatation fields to links model [\#344](https://github.com/Materials-Consortia/optimade-python-tools/issues/344)
- Updates to models \(nperiodic\_dimensions\) [\#343](https://github.com/Materials-Consortia/optimade-python-tools/issues/343)
- Updates to models \(changing unknown atoms\) [\#342](https://github.com/Materials-Consortia/optimade-python-tools/issues/342)
- Improvements/fixes for openapi.json [\#332](https://github.com/Materials-Consortia/optimade-python-tools/issues/332)
- Update to v1.0.0-rc.1 [\#329](https://github.com/Materials-Consortia/optimade-python-tools/issues/329)
- Decouple updates in providers repo [\#311](https://github.com/Materials-Consortia/optimade-python-tools/issues/311)
- RST not rendering with mkdocs [\#307](https://github.com/Materials-Consortia/optimade-python-tools/issues/307)

**Merged pull requests:**

- Retrieve providers list if no submodule is found [\#374](https://github.com/Materials-Consortia/optimade-python-tools/pull/374) ([CasperWA](https://github.com/CasperWA))
- Update default implementation information [\#372](https://github.com/Materials-Consortia/optimade-python-tools/pull/372) ([shyamd](https://github.com/shyamd))
- Bump spec version to 1.0.0-rc.2 [\#367](https://github.com/Materials-Consortia/optimade-python-tools/pull/367) ([ml-evs](https://github.com/ml-evs))
- Dependabot updates: numpy, mkdocs-material, mkdocstrings, requests [\#364](https://github.com/Materials-Consortia/optimade-python-tools/pull/364) ([ml-evs](https://github.com/ml-evs))
- Merge all Dependabot updates [\#353](https://github.com/Materials-Consortia/optimade-python-tools/pull/353) ([shyamd](https://github.com/shyamd))
- Update model descriptions and openapi.json for 1.0.0-rc2 [\#351](https://github.com/Materials-Consortia/optimade-python-tools/pull/351) ([ml-evs](https://github.com/ml-evs))
- Update models according to changes during CECAM 2020 meeting [\#350](https://github.com/Materials-Consortia/optimade-python-tools/pull/350) ([ml-evs](https://github.com/ml-evs))
- Decouple changes in providers repo [\#312](https://github.com/Materials-Consortia/optimade-python-tools/pull/312) ([shyamd](https://github.com/shyamd))

## [v0.9.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.1) (2020-06-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.9.0...v0.9.1)

## [v0.9.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.9.0) (2020-06-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.8.1...v0.9.0)

**Implemented enhancements:**

- Breaking up the python tools into seperable packages [\#255](https://github.com/Materials-Consortia/optimade-python-tools/issues/255)
- Run both servers as standard [\#238](https://github.com/Materials-Consortia/optimade-python-tools/issues/238)

**Fixed bugs:**

- Non-running CI job [\#331](https://github.com/Materials-Consortia/optimade-python-tools/issues/331)
- Special species "X" not tested for non-disordered structures [\#304](https://github.com/Materials-Consortia/optimade-python-tools/issues/304)
- Standardize timezone of datetime responses [\#288](https://github.com/Materials-Consortia/optimade-python-tools/issues/288)
- Queries on aliased/provider fields are broken for nested properties [\#282](https://github.com/Materials-Consortia/optimade-python-tools/issues/282)
- General exceptions not being put into response [\#281](https://github.com/Materials-Consortia/optimade-python-tools/issues/281)
- Issue with CIF export [\#271](https://github.com/Materials-Consortia/optimade-python-tools/issues/271)
- Type-cast inputs for general Error [\#280](https://github.com/Materials-Consortia/optimade-python-tools/pull/280) ([CasperWA](https://github.com/CasperWA))

**Security fixes:**

- \[Security\] Bump django from 3.0.4 to 3.0.7 in /.github/workflows [\#291](https://github.com/Materials-Consortia/optimade-python-tools/pull/291) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))

**Closed issues:**

- Update links resources [\#299](https://github.com/Materials-Consortia/optimade-python-tools/issues/299)
- Need to set up mkdocs [\#289](https://github.com/Materials-Consortia/optimade-python-tools/issues/289)
- Need to add custom schema entries for unit/sortable \(and eventually type\) [\#278](https://github.com/Materials-Consortia/optimade-python-tools/issues/278)
- /info/\<entry-endpoint\> missing `sortable` key under each property [\#273](https://github.com/Materials-Consortia/optimade-python-tools/issues/273)
- Make CI linting more useful [\#269](https://github.com/Materials-Consortia/optimade-python-tools/issues/269)
- \[PR SPECIFIC\] Reminder: Validator test pinned to specific commit [\#268](https://github.com/Materials-Consortia/optimade-python-tools/issues/268)
- Validator does not check that pagination links work [\#265](https://github.com/Materials-Consortia/optimade-python-tools/issues/265)
- available\_api\_versions is not correctly validated [\#261](https://github.com/Materials-Consortia/optimade-python-tools/issues/261)
- Implementation model should allow for any URL type in `source_url` [\#260](https://github.com/Materials-Consortia/optimade-python-tools/issues/260)
- Extra structure endpoints in the api specification @ odbx [\#259](https://github.com/Materials-Consortia/optimade-python-tools/issues/259)
- Wrong response structure at info endpoint @ cod [\#258](https://github.com/Materials-Consortia/optimade-python-tools/issues/258)
- Missing base url for api's docs @ materialscloud [\#257](https://github.com/Materials-Consortia/optimade-python-tools/issues/257)
- Handling of KNOWN in mongo backend [\#254](https://github.com/Materials-Consortia/optimade-python-tools/issues/254)
- `None` values in `lattice_vectors` [\#170](https://github.com/Materials-Consortia/optimade-python-tools/issues/170)
- Make sure that the PyPI distribution works [\#143](https://github.com/Materials-Consortia/optimade-python-tools/issues/143)
- Move run.sh to a python file to be environment-agnostic [\#81](https://github.com/Materials-Consortia/optimade-python-tools/issues/81)

**Merged pull requests:**

- Another fix for release pipeline [\#355](https://github.com/Materials-Consortia/optimade-python-tools/pull/355) ([shyamd](https://github.com/shyamd))
- Fix publish workflow [\#354](https://github.com/Materials-Consortia/optimade-python-tools/pull/354) ([CasperWA](https://github.com/CasperWA))
- Fix publish workflow [\#352](https://github.com/Materials-Consortia/optimade-python-tools/pull/352) ([CasperWA](https://github.com/CasperWA))
- Update publish workflow [\#340](https://github.com/Materials-Consortia/optimade-python-tools/pull/340) ([shyamd](https://github.com/shyamd))
- Remove test publish action [\#338](https://github.com/Materials-Consortia/optimade-python-tools/pull/338) ([shyamd](https://github.com/shyamd))
- Fix 'publish\_TestPyPI' CI job [\#337](https://github.com/Materials-Consortia/optimade-python-tools/pull/337) ([CasperWA](https://github.com/CasperWA))
- Specify versions for all setup.py deps [\#336](https://github.com/Materials-Consortia/optimade-python-tools/pull/336) ([CasperWA](https://github.com/CasperWA))
- Represent the datetime objects as UTC in RFC3339 format [\#333](https://github.com/Materials-Consortia/optimade-python-tools/pull/333) ([fekad](https://github.com/fekad))
- dependamat: Bump \<package\_name\> v x.y.z to vx.y.\(z+1\) [\#330](https://github.com/Materials-Consortia/optimade-python-tools/pull/330) ([ml-evs](https://github.com/ml-evs))
- Bump fastapi from 0.53.1 to 0.56.0 [\#324](https://github.com/Materials-Consortia/optimade-python-tools/pull/324) ([dependabot[bot]](https://github.com/apps/dependabot))
- Bump pydantic from 1.4 to 1.5.1 [\#320](https://github.com/Materials-Consortia/optimade-python-tools/pull/320) ([dependabot[bot]](https://github.com/apps/dependabot))
- Update links resources [\#306](https://github.com/Materials-Consortia/optimade-python-tools/pull/306) ([CasperWA](https://github.com/CasperWA))
- Add special species for adapters testing [\#305](https://github.com/Materials-Consortia/optimade-python-tools/pull/305) ([CasperWA](https://github.com/CasperWA))
- Clean Up Build Environment [\#301](https://github.com/Materials-Consortia/optimade-python-tools/pull/301) ([shyamd](https://github.com/shyamd))
- Enable CI failures for linting [\#300](https://github.com/Materials-Consortia/optimade-python-tools/pull/300) ([ml-evs](https://github.com/ml-evs))
- Adding jarvis-tools structures [\#297](https://github.com/Materials-Consortia/optimade-python-tools/pull/297) ([knc6](https://github.com/knc6))
- Update Docs [\#295](https://github.com/Materials-Consortia/optimade-python-tools/pull/295) ([shyamd](https://github.com/shyamd))
- Setup MKDocs for Documentation [\#294](https://github.com/Materials-Consortia/optimade-python-tools/pull/294) ([shyamd](https://github.com/shyamd))
- Fix filters on nested provider/aliased fields [\#285](https://github.com/Materials-Consortia/optimade-python-tools/pull/285) ([ml-evs](https://github.com/ml-evs))
- Use heroku-shields instead of heroku-badge [\#284](https://github.com/Materials-Consortia/optimade-python-tools/pull/284) ([CasperWA](https://github.com/CasperWA))
- Add OPTIMADE logo to badge by extending JSON [\#283](https://github.com/Materials-Consortia/optimade-python-tools/pull/283) ([CasperWA](https://github.com/CasperWA))
- Add null check to mongo filtertransformer for KNOWN/UNKNOWN filters [\#279](https://github.com/Materials-Consortia/optimade-python-tools/pull/279) ([ml-evs](https://github.com/ml-evs))
- Add `sortable=True` to all properties [\#274](https://github.com/Materials-Consortia/optimade-python-tools/pull/274) ([CasperWA](https://github.com/CasperWA))
- Make \_atom\_site\_label unique in CIF generation [\#272](https://github.com/Materials-Consortia/optimade-python-tools/pull/272) ([CasperWA](https://github.com/CasperWA))
- Not so quick fix to allow "/" at end of validator URL, plus fixes and tests for --as\_type [\#267](https://github.com/Materials-Consortia/optimade-python-tools/pull/267) ([ml-evs](https://github.com/ml-evs))
- Check pagination links-\>next with validator [\#266](https://github.com/Materials-Consortia/optimade-python-tools/pull/266) ([ml-evs](https://github.com/ml-evs))
- Relax HTTP URL constraints on meta-\>implementation-\>source\_url field. [\#262](https://github.com/Materials-Consortia/optimade-python-tools/pull/262) ([ml-evs](https://github.com/ml-evs))
- Validate lattice\_vectors for all null or all float [\#171](https://github.com/Materials-Consortia/optimade-python-tools/pull/171) ([CasperWA](https://github.com/CasperWA))

## [v0.8.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.8.1) (2020-04-25)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.8.0...v0.8.1)

**Fixed bugs:**

- Pip install missing some files [\#252](https://github.com/Materials-Consortia/optimade-python-tools/issues/252)

**Merged pull requests:**

- v0.8.1 hotfix [\#256](https://github.com/Materials-Consortia/optimade-python-tools/pull/256) ([ml-evs](https://github.com/ml-evs))
- Fix 252 missing landing page [\#253](https://github.com/Materials-Consortia/optimade-python-tools/pull/253) ([shyamd](https://github.com/shyamd))

## [v0.8.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.8.0) (2020-04-22)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.7.1...v0.8.0)

**Implemented enhancements:**

- Switch to pydantic's BaseSettings for the config file? [\#152](https://github.com/Materials-Consortia/optimade-python-tools/issues/152)
- Use services for testing/updating dependencies? [\#96](https://github.com/Materials-Consortia/optimade-python-tools/issues/96)
- Remove query constraints for /links-endpoint [\#244](https://github.com/Materials-Consortia/optimade-python-tools/pull/244) ([CasperWA](https://github.com/CasperWA))
- Add adapters - Base design + 'structures' \(+ 'references'... sort of\) [\#241](https://github.com/Materials-Consortia/optimade-python-tools/pull/241) ([CasperWA](https://github.com/CasperWA))
- Add dependabot and last commit date badges [\#237](https://github.com/Materials-Consortia/optimade-python-tools/pull/237) ([CasperWA](https://github.com/CasperWA))
- Add mongo length operator functionality with length aliases [\#222](https://github.com/Materials-Consortia/optimade-python-tools/pull/222) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Use Path.home\(\) instead of ~ in default config path values [\#245](https://github.com/Materials-Consortia/optimade-python-tools/issues/245)

**Closed issues:**

- Have Dependabot take care of various requirements.txt files as well [\#249](https://github.com/Materials-Consortia/optimade-python-tools/issues/249)
- Remove commented out GH Action job `deps_clean-install` [\#247](https://github.com/Materials-Consortia/optimade-python-tools/issues/247)
- Local testing fails without default config [\#239](https://github.com/Materials-Consortia/optimade-python-tools/issues/239)
- Release only when pushing to master [\#229](https://github.com/Materials-Consortia/optimade-python-tools/issues/229)
- Do we need `server.cfg`? [\#134](https://github.com/Materials-Consortia/optimade-python-tools/issues/134)
- Implement LENGTH in query [\#86](https://github.com/Materials-Consortia/optimade-python-tools/issues/86)

**Merged pull requests:**

- Up to v0.8.0 [\#251](https://github.com/Materials-Consortia/optimade-python-tools/pull/251) ([CasperWA](https://github.com/CasperWA))
- Remove old commented GH Action job [\#250](https://github.com/Materials-Consortia/optimade-python-tools/pull/250) ([CasperWA](https://github.com/CasperWA))
- Use Path.home\(\) instead of `~` [\#246](https://github.com/Materials-Consortia/optimade-python-tools/pull/246) ([CasperWA](https://github.com/CasperWA))
- Fix path in default config [\#243](https://github.com/Materials-Consortia/optimade-python-tools/pull/243) ([ml-evs](https://github.com/ml-evs))
- Fixes Local Tests [\#240](https://github.com/Materials-Consortia/optimade-python-tools/pull/240) ([shyamd](https://github.com/shyamd))
- Revert "Fix github actions for non-release tags" [\#236](https://github.com/Materials-Consortia/optimade-python-tools/pull/236) ([shyamd](https://github.com/shyamd))
- Enable filtering on relationships with mongo  [\#234](https://github.com/Materials-Consortia/optimade-python-tools/pull/234) ([ml-evs](https://github.com/ml-evs))
- Update filter examples and validate optional cases [\#227](https://github.com/Materials-Consortia/optimade-python-tools/pull/227) ([ml-evs](https://github.com/ml-evs))
- Switch from config init to BaseSettings [\#226](https://github.com/Materials-Consortia/optimade-python-tools/pull/226) ([shyamd](https://github.com/shyamd))

## [v0.7.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.7.1) (2020-03-16)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.7.0...v0.7.1)

**Closed issues:**

- Fix all capitalisation of OPTIMADE [\#232](https://github.com/Materials-Consortia/optimade-python-tools/issues/232)
- Remove validator action from README [\#230](https://github.com/Materials-Consortia/optimade-python-tools/issues/230)

**Merged pull requests:**

- Fix github actions for non-release tags [\#235](https://github.com/Materials-Consortia/optimade-python-tools/pull/235) ([shyamd](https://github.com/shyamd))
- Update OPTIMADE capitalisation [\#233](https://github.com/Materials-Consortia/optimade-python-tools/pull/233) ([ml-evs](https://github.com/ml-evs))
- Update mentions of action in readme [\#231](https://github.com/Materials-Consortia/optimade-python-tools/pull/231) ([ml-evs](https://github.com/ml-evs))

## [v0.7.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.7.0) (2020-03-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.6.0...v0.7.0)

**Implemented enhancements:**

- Validate all non-optional :filter: examples from the spec [\#213](https://github.com/Materials-Consortia/optimade-python-tools/pull/213) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- Some mandatory filter examples from spec do not work [\#217](https://github.com/Materials-Consortia/optimade-python-tools/issues/217)
- Add txt-files in optimade.validator.data to MANIFEST [\#225](https://github.com/Materials-Consortia/optimade-python-tools/pull/225) ([CasperWA](https://github.com/CasperWA))
- Handle arbitrary nested NOT/AND/OR in queries [\#221](https://github.com/Materials-Consortia/optimade-python-tools/pull/221) ([ml-evs](https://github.com/ml-evs))

**Closed issues:**

- Validator only validates what we have working, not what is required by the spec [\#182](https://github.com/Materials-Consortia/optimade-python-tools/issues/182)

**Merged pull requests:**

- v0.7.0 release [\#228](https://github.com/Materials-Consortia/optimade-python-tools/pull/228) ([ml-evs](https://github.com/ml-evs))
- Remove GH Action to validate OPTiMaDe instances [\#224](https://github.com/Materials-Consortia/optimade-python-tools/pull/224) ([CasperWA](https://github.com/CasperWA))
- Codecov-action supports token-less uploads [\#220](https://github.com/Materials-Consortia/optimade-python-tools/pull/220) ([CasperWA](https://github.com/CasperWA))
- Update django requirement from \>=2.2.9,~=2.2 to \>=2.2,\<4.0 [\#219](https://github.com/Materials-Consortia/optimade-python-tools/pull/219) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))
- Update elasticsearch-dsl requirement from ~=6.4 to \>=6.4,\<8.0 [\#218](https://github.com/Materials-Consortia/optimade-python-tools/pull/218) ([dependabot-preview[bot]](https://github.com/apps/dependabot-preview))

## [v0.6.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.6.0) (2020-03-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.5.0...v0.6.0)

**Implemented enhancements:**

- Possibly add CORS middleware [\#159](https://github.com/Materials-Consortia/optimade-python-tools/issues/159)
- Add debug flag to server [\#130](https://github.com/Materials-Consortia/optimade-python-tools/issues/130)
- Make validator GitHub Action [\#191](https://github.com/Materials-Consortia/optimade-python-tools/pull/191) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- meta/query/representation value not cutting off version properly [\#199](https://github.com/Materials-Consortia/optimade-python-tools/issues/199)
- URL for providers.json from Materials-Consortia has changed [\#186](https://github.com/Materials-Consortia/optimade-python-tools/issues/186)
- Relationships don't work when "/" present in id [\#181](https://github.com/Materials-Consortia/optimade-python-tools/issues/181)
- Redirect middleware not hitting single-entry endpoints [\#174](https://github.com/Materials-Consortia/optimade-python-tools/issues/174)

**Closed issues:**

- /info/ reports wrong url under available\_api\_versions [\#215](https://github.com/Materials-Consortia/optimade-python-tools/issues/215)
- Query parameters not handled correctly [\#208](https://github.com/Materials-Consortia/optimade-python-tools/issues/208)
- Test for AvailableApiVersion is correct for the wrong reasons [\#204](https://github.com/Materials-Consortia/optimade-python-tools/issues/204)
- Drop '/optimade' from paths in openapi.json [\#197](https://github.com/Materials-Consortia/optimade-python-tools/issues/197)
- heroku is failing [\#185](https://github.com/Materials-Consortia/optimade-python-tools/issues/185)
- List properties and HAS \_ operators missing [\#98](https://github.com/Materials-Consortia/optimade-python-tools/issues/98)
- Checklist for OPTiMaDe v0.10.1 [\#29](https://github.com/Materials-Consortia/optimade-python-tools/issues/29)

**Merged pull requests:**

- Removed /optimade/ prefix in info response [\#216](https://github.com/Materials-Consortia/optimade-python-tools/pull/216) ([ml-evs](https://github.com/ml-evs))
- Self load data [\#212](https://github.com/Materials-Consortia/optimade-python-tools/pull/212) ([shyamd](https://github.com/shyamd))
- Update tests for available\_api\_versions [\#211](https://github.com/Materials-Consortia/optimade-python-tools/pull/211) ([CasperWA](https://github.com/CasperWA))
- Up to v0.6.0 [\#210](https://github.com/Materials-Consortia/optimade-python-tools/pull/210) ([CasperWA](https://github.com/CasperWA))
- Update handling of include parameter \(and other query parameters\) [\#209](https://github.com/Materials-Consortia/optimade-python-tools/pull/209) ([CasperWA](https://github.com/CasperWA))
- Skip HAS ONLY test if mongomock version \<= 3.19.0 [\#206](https://github.com/Materials-Consortia/optimade-python-tools/pull/206) ([ml-evs](https://github.com/ml-evs))
- Test mandatory queries in validator [\#205](https://github.com/Materials-Consortia/optimade-python-tools/pull/205) ([ml-evs](https://github.com/ml-evs))
- Fix include query parameter [\#202](https://github.com/Materials-Consortia/optimade-python-tools/pull/202) ([CasperWA](https://github.com/CasperWA))
- Fix meta.query.representation and remove /optimade in base URLs [\#201](https://github.com/Materials-Consortia/optimade-python-tools/pull/201) ([CasperWA](https://github.com/CasperWA))
- Use mongo for CI [\#196](https://github.com/Materials-Consortia/optimade-python-tools/pull/196) ([ml-evs](https://github.com/ml-evs))
- \(Cosmetic\) updates to models [\#195](https://github.com/Materials-Consortia/optimade-python-tools/pull/195) ([CasperWA](https://github.com/CasperWA))
- Add CORSMiddleware [\#194](https://github.com/Materials-Consortia/optimade-python-tools/pull/194) ([CasperWA](https://github.com/CasperWA))
- Add "debug mode" [\#190](https://github.com/Materials-Consortia/optimade-python-tools/pull/190) ([CasperWA](https://github.com/CasperWA))
- Use https://provider.optimade.org/providers.json [\#187](https://github.com/Materials-Consortia/optimade-python-tools/pull/187) ([CasperWA](https://github.com/CasperWA))
- Fix errors parsing IDs that contain slashes [\#183](https://github.com/Materials-Consortia/optimade-python-tools/pull/183) ([ml-evs](https://github.com/ml-evs))
- Added default mongo implementations for HAS ALL/ANY/ONLY [\#173](https://github.com/Materials-Consortia/optimade-python-tools/pull/173) ([ml-evs](https://github.com/ml-evs))

## [v0.5.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.5.0) (2020-02-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.4.0...v0.5.0)

**Implemented enhancements:**

- Implement a landing page for requests to the base URL [\#169](https://github.com/Materials-Consortia/optimade-python-tools/issues/169)

**Fixed bugs:**

- 'minor' and 'patch' versioned base URL prefixes are wrong [\#177](https://github.com/Materials-Consortia/optimade-python-tools/issues/177)

**Closed issues:**

- Handle `include` standard JSON API query parameter [\#94](https://github.com/Materials-Consortia/optimade-python-tools/issues/94)

**Merged pull requests:**

- Bump to v0.5.0 [\#179](https://github.com/Materials-Consortia/optimade-python-tools/pull/179) ([CasperWA](https://github.com/CasperWA))
- Correctly create optional versioned base URLs [\#178](https://github.com/Materials-Consortia/optimade-python-tools/pull/178) ([CasperWA](https://github.com/CasperWA))
- Make mapper aliases configurable [\#175](https://github.com/Materials-Consortia/optimade-python-tools/pull/175) ([ml-evs](https://github.com/ml-evs))
- Add landing page at base URL [\#172](https://github.com/Materials-Consortia/optimade-python-tools/pull/172) ([ml-evs](https://github.com/ml-evs))
- Implement `include` query parameter [\#163](https://github.com/Materials-Consortia/optimade-python-tools/pull/163) ([CasperWA](https://github.com/CasperWA))
- Add docker for index meta-database [\#140](https://github.com/Materials-Consortia/optimade-python-tools/pull/140) ([CasperWA](https://github.com/CasperWA))

## [v0.4.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.4.0) (2020-02-06)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.4...v0.4.0)

**Implemented enhancements:**

- switch to pipenv? [\#37](https://github.com/Materials-Consortia/optimade-python-tools/issues/37)
- Reorder tests [\#162](https://github.com/Materials-Consortia/optimade-python-tools/pull/162) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Server app intermingles [\#161](https://github.com/Materials-Consortia/optimade-python-tools/issues/161)
- `response_fields` not working [\#154](https://github.com/Materials-Consortia/optimade-python-tools/issues/154)

**Closed issues:**

- Change `page_page` to `page_number` [\#165](https://github.com/Materials-Consortia/optimade-python-tools/issues/165)
- Add schema-relevant parameters to query parameters [\#164](https://github.com/Materials-Consortia/optimade-python-tools/issues/164)
- Alias optimade/structures/ to optimade/structure [\#128](https://github.com/Materials-Consortia/optimade-python-tools/issues/128)
- Minor changes to specification v0.10.1-develop [\#115](https://github.com/Materials-Consortia/optimade-python-tools/issues/115)
- Update models with new levels of REQUIRED response properties [\#114](https://github.com/Materials-Consortia/optimade-python-tools/issues/114)
- Constraining list/array types in the schema [\#55](https://github.com/Materials-Consortia/optimade-python-tools/issues/55)

**Merged pull requests:**

- Bump to v0.4.0 [\#168](https://github.com/Materials-Consortia/optimade-python-tools/pull/168) ([CasperWA](https://github.com/CasperWA))
- Describe query parameters in OpenAPI schema [\#166](https://github.com/Materials-Consortia/optimade-python-tools/pull/166) ([CasperWA](https://github.com/CasperWA))
- Redirect slashed URLs [\#160](https://github.com/Materials-Consortia/optimade-python-tools/pull/160) ([CasperWA](https://github.com/CasperWA))
- New REQUIRED level properties [\#153](https://github.com/Materials-Consortia/optimade-python-tools/pull/153) ([CasperWA](https://github.com/CasperWA))

## [v0.3.4](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.4) (2020-02-04)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.3...v0.3.4)

**Implemented enhancements:**

- Include `develop` or not? Default branch? - Create INSTALL.md [\#136](https://github.com/Materials-Consortia/optimade-python-tools/issues/136)

**Fixed bugs:**

- Excepting non-existent exception [\#129](https://github.com/Materials-Consortia/optimade-python-tools/issues/129)

**Closed issues:**

- disable serving API under /v0.10 and /v0.10.0 by default? [\#122](https://github.com/Materials-Consortia/optimade-python-tools/issues/122)
- PyPI release checklist [\#67](https://github.com/Materials-Consortia/optimade-python-tools/issues/67)

**Merged pull requests:**

- Bump to v0.3.4 [\#158](https://github.com/Materials-Consortia/optimade-python-tools/pull/158) ([CasperWA](https://github.com/CasperWA))
- Fix heroku badge [\#157](https://github.com/Materials-Consortia/optimade-python-tools/pull/157) ([ml-evs](https://github.com/ml-evs))
- Move installation instructions [\#156](https://github.com/Materials-Consortia/optimade-python-tools/pull/156) ([ml-evs](https://github.com/ml-evs))
- Update base URLs [\#155](https://github.com/Materials-Consortia/optimade-python-tools/pull/155) ([CasperWA](https://github.com/CasperWA))
- Extend OpenAPI/spec description [\#151](https://github.com/Materials-Consortia/optimade-python-tools/pull/151) ([CasperWA](https://github.com/CasperWA))
- Non Local Mongo [\#150](https://github.com/Materials-Consortia/optimade-python-tools/pull/150) ([shyamd](https://github.com/shyamd))

## [v0.3.3](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.3) (2020-01-24)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.2...v0.3.3)

**Fixed bugs:**

- Lark files not being distributed [\#141](https://github.com/Materials-Consortia/optimade-python-tools/issues/141)

**Closed issues:**

- Tests fail with lark-parser\>=0.8 [\#146](https://github.com/Materials-Consortia/optimade-python-tools/issues/146)

**Merged pull requests:**

- Updated lark-parser to 0.8.1 [\#149](https://github.com/Materials-Consortia/optimade-python-tools/pull/149) ([ml-evs](https://github.com/ml-evs))
- Split eager and standard tests to avoid unnecessary badge of shame [\#148](https://github.com/Materials-Consortia/optimade-python-tools/pull/148) ([ml-evs](https://github.com/ml-evs))
- Bump to v0.3.3 [\#147](https://github.com/Materials-Consortia/optimade-python-tools/pull/147) ([CasperWA](https://github.com/CasperWA))
- Fix root\_validator issues with optional fields and made meta optional [\#145](https://github.com/Materials-Consortia/optimade-python-tools/pull/145) ([ml-evs](https://github.com/ml-evs))
- Handle `JSONDecodeError`s in validator [\#144](https://github.com/Materials-Consortia/optimade-python-tools/pull/144) ([ml-evs](https://github.com/ml-evs))

## [v0.3.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.2) (2020-01-20)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.1...v0.3.2)

**Implemented enhancements:**

- Add base URL to configuration file [\#135](https://github.com/Materials-Consortia/optimade-python-tools/pull/135) ([CasperWA](https://github.com/CasperWA))

**Fixed bugs:**

- Fix `load_from_json` [\#137](https://github.com/Materials-Consortia/optimade-python-tools/pull/137) ([CasperWA](https://github.com/CasperWA))

**Merged pull requests:**

- Make sure relevant package data is included in distributions [\#142](https://github.com/Materials-Consortia/optimade-python-tools/pull/142) ([CasperWA](https://github.com/CasperWA))
- Add database page limit [\#139](https://github.com/Materials-Consortia/optimade-python-tools/pull/139) ([CasperWA](https://github.com/CasperWA))

## [v0.3.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.1) (2020-01-17)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.3.0...v0.3.1)

**Merged pull requests:**

- Update requirements [\#138](https://github.com/Materials-Consortia/optimade-python-tools/pull/138) ([CasperWA](https://github.com/CasperWA))

## [v0.3.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.3.0) (2020-01-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.2...v0.3.0)

**Implemented enhancements:**

- Implement optional `implementation` in top-level meta response [\#117](https://github.com/Materials-Consortia/optimade-python-tools/issues/117)
- Create "special" index meta-database server [\#100](https://github.com/Materials-Consortia/optimade-python-tools/issues/100)
- Implement relationships in server [\#71](https://github.com/Materials-Consortia/optimade-python-tools/issues/71)
- Add missing /references endpoint to server [\#69](https://github.com/Materials-Consortia/optimade-python-tools/issues/69)
- Automatically publish version tags to PyPI via GH Actions [\#107](https://github.com/Materials-Consortia/optimade-python-tools/pull/107) ([CasperWA](https://github.com/CasperWA))
- Using routers [\#99](https://github.com/Materials-Consortia/optimade-python-tools/pull/99) ([CasperWA](https://github.com/CasperWA))
- Add relationships functionality [\#91](https://github.com/Materials-Consortia/optimade-python-tools/pull/91) ([ml-evs](https://github.com/ml-evs))
- Added external API validator based on our pydantic models [\#74](https://github.com/Materials-Consortia/optimade-python-tools/pull/74) ([ml-evs](https://github.com/ml-evs))

**Fixed bugs:**

- The invoke task `update-openapijson` is incomplete [\#123](https://github.com/Materials-Consortia/optimade-python-tools/issues/123)
- Django vulnerability [\#108](https://github.com/Materials-Consortia/optimade-python-tools/issues/108)

**Closed issues:**

- info endpoint duplicated? [\#120](https://github.com/Materials-Consortia/optimade-python-tools/issues/120)
- Commented-out validator [\#111](https://github.com/Materials-Consortia/optimade-python-tools/issues/111)
- FastAPI v0.44.0 supports pydantic \> 1.0.0 [\#101](https://github.com/Materials-Consortia/optimade-python-tools/issues/101)
- Server is missing /links endpoint [\#89](https://github.com/Materials-Consortia/optimade-python-tools/issues/89)
- Make sure all validators are tested [\#87](https://github.com/Materials-Consortia/optimade-python-tools/issues/87)
- The `sortable` field must be added to models [\#84](https://github.com/Materials-Consortia/optimade-python-tools/issues/84)
- Package structure [\#72](https://github.com/Materials-Consortia/optimade-python-tools/issues/72)
- Possibly make /info/{endpoint} dynamic [\#70](https://github.com/Materials-Consortia/optimade-python-tools/issues/70)
- setuptools package with server as "extra" [\#62](https://github.com/Materials-Consortia/optimade-python-tools/issues/62)
- use examples from specs as resources [\#57](https://github.com/Materials-Consortia/optimade-python-tools/issues/57)
- httptools dependency has build issues on GCC/Linux [\#54](https://github.com/Materials-Consortia/optimade-python-tools/issues/54)
- Lark grammar file for v0.9.8 [\#50](https://github.com/Materials-Consortia/optimade-python-tools/issues/50)
- type is missing in response [\#43](https://github.com/Materials-Consortia/optimade-python-tools/issues/43)
- Enforce use of autoformatter  [\#33](https://github.com/Materials-Consortia/optimade-python-tools/issues/33)
- switch license to MIT [\#28](https://github.com/Materials-Consortia/optimade-python-tools/issues/28)
- write a lark JSONTransformer / JSONdecoder [\#26](https://github.com/Materials-Consortia/optimade-python-tools/issues/26)
- server.jsonapi has no additionalProperties=false  [\#23](https://github.com/Materials-Consortia/optimade-python-tools/issues/23)
- server.jsonapi has no patternProperties  [\#22](https://github.com/Materials-Consortia/optimade-python-tools/issues/22)
- Developer-friendly pre-commit openapi.json visual diff [\#21](https://github.com/Materials-Consortia/optimade-python-tools/issues/21)
- add JSON schema API [\#12](https://github.com/Materials-Consortia/optimade-python-tools/issues/12)
- generate static documentation on github from openapi.json [\#9](https://github.com/Materials-Consortia/optimade-python-tools/issues/9)
- test how to generate a client from the openapi.json [\#8](https://github.com/Materials-Consortia/optimade-python-tools/issues/8)
- come up with suggested toolchain for validating existing optimade API against openapi.json [\#7](https://github.com/Materials-Consortia/optimade-python-tools/issues/7)
- add travis test that checks openapi.json is valid OpenAPI spec [\#6](https://github.com/Materials-Consortia/optimade-python-tools/issues/6)
- add 2 examples of how to include documentation in python classes [\#5](https://github.com/Materials-Consortia/optimade-python-tools/issues/5)
- add one-line command to update openapi.json [\#4](https://github.com/Materials-Consortia/optimade-python-tools/issues/4)

**Merged pull requests:**

- Fixed CI readme badge [\#133](https://github.com/Materials-Consortia/optimade-python-tools/pull/133) ([ml-evs](https://github.com/ml-evs))
- Add meta.description to BaseRelationshipResource [\#131](https://github.com/Materials-Consortia/optimade-python-tools/pull/131) ([CasperWA](https://github.com/CasperWA))
- Added homepage attribute to LinksResource [\#127](https://github.com/Materials-Consortia/optimade-python-tools/pull/127) ([ml-evs](https://github.com/ml-evs))
- Updated structure models and validators [\#126](https://github.com/Materials-Consortia/optimade-python-tools/pull/126) ([ml-evs](https://github.com/ml-evs))
- Minor change to fallback server.cfg [\#125](https://github.com/Materials-Consortia/optimade-python-tools/pull/125) ([ml-evs](https://github.com/ml-evs))
- Update local OpenAPI schemes prior to copying [\#124](https://github.com/Materials-Consortia/optimade-python-tools/pull/124) ([CasperWA](https://github.com/CasperWA))
- Update OpenAPI tags [\#121](https://github.com/Materials-Consortia/optimade-python-tools/pull/121) ([CasperWA](https://github.com/CasperWA))
- A few fixes related to usage as a library [\#119](https://github.com/Materials-Consortia/optimade-python-tools/pull/119) ([ml-evs](https://github.com/ml-evs))
- Add implementation to top-level meta response [\#118](https://github.com/Materials-Consortia/optimade-python-tools/pull/118) ([CasperWA](https://github.com/CasperWA))
- Add heroku deployment scripts [\#116](https://github.com/Materials-Consortia/optimade-python-tools/pull/116) ([ltalirz](https://github.com/ltalirz))
- Reorganize package [\#113](https://github.com/Materials-Consortia/optimade-python-tools/pull/113) ([CasperWA](https://github.com/CasperWA))
- Introduce grammar v0.10.1 [\#112](https://github.com/Materials-Consortia/optimade-python-tools/pull/112) ([CasperWA](https://github.com/CasperWA))
- Update to pydantic v1 [\#110](https://github.com/Materials-Consortia/optimade-python-tools/pull/110) ([CasperWA](https://github.com/CasperWA))
- Minimum requirement of django v2.2.8 [\#109](https://github.com/Materials-Consortia/optimade-python-tools/pull/109) ([CasperWA](https://github.com/CasperWA))
- Index meta-database [\#103](https://github.com/Materials-Consortia/optimade-python-tools/pull/103) ([CasperWA](https://github.com/CasperWA))
- restrict pydantic version [\#97](https://github.com/Materials-Consortia/optimade-python-tools/pull/97) ([ltalirz](https://github.com/ltalirz))
- Add /links [\#95](https://github.com/Materials-Consortia/optimade-python-tools/pull/95) ([CasperWA](https://github.com/CasperWA))
- Fix data\_returned and data\_available [\#93](https://github.com/Materials-Consortia/optimade-python-tools/pull/93) ([CasperWA](https://github.com/CasperWA))
- Use GitHub Actions for CI [\#92](https://github.com/Materials-Consortia/optimade-python-tools/pull/92) ([ml-evs](https://github.com/ml-evs))
- Remove inappropriate lint messages [\#90](https://github.com/Materials-Consortia/optimade-python-tools/pull/90) ([CasperWA](https://github.com/CasperWA))
- Fix dependencies [\#88](https://github.com/Materials-Consortia/optimade-python-tools/pull/88) ([CasperWA](https://github.com/CasperWA))
- Add sortable field to EntryInfoProperty model [\#85](https://github.com/Materials-Consortia/optimade-python-tools/pull/85) ([CasperWA](https://github.com/CasperWA))
- Validate illegal fields are not present under attributes and relationships [\#83](https://github.com/Materials-Consortia/optimade-python-tools/pull/83) ([CasperWA](https://github.com/CasperWA))
- Add references endpoint [\#78](https://github.com/Materials-Consortia/optimade-python-tools/pull/78) ([CasperWA](https://github.com/CasperWA))
- fix travis build [\#77](https://github.com/Materials-Consortia/optimade-python-tools/pull/77) ([ltalirz](https://github.com/ltalirz))
- Fix manual verification of elements\_ratios [\#76](https://github.com/Materials-Consortia/optimade-python-tools/pull/76) ([CasperWA](https://github.com/CasperWA))
- add automatic PyPI deployment [\#75](https://github.com/Materials-Consortia/optimade-python-tools/pull/75) ([ltalirz](https://github.com/ltalirz))
- Remove reference to `"all"` endpoint and rename collections submodule [\#73](https://github.com/Materials-Consortia/optimade-python-tools/pull/73) ([ml-evs](https://github.com/ml-evs))
- Updates to README and docs for v0.10.0 [\#68](https://github.com/Materials-Consortia/optimade-python-tools/pull/68) ([ml-evs](https://github.com/ml-evs))
- Adding grammar for v0.10.0 [\#66](https://github.com/Materials-Consortia/optimade-python-tools/pull/66) ([fekad](https://github.com/fekad))
- Schema updates and fixes relative to the v0.10.0 spec [\#65](https://github.com/Materials-Consortia/optimade-python-tools/pull/65) ([ml-evs](https://github.com/ml-evs))
- Break requirements down on per backend basis [\#64](https://github.com/Materials-Consortia/optimade-python-tools/pull/64) ([ml-evs](https://github.com/ml-evs))
- 0.10.0 grammer, elasticsearch transformer, setuptools extra [\#63](https://github.com/Materials-Consortia/optimade-python-tools/pull/63) ([markus1978](https://github.com/markus1978))
- Added a Lark to Django Query converter [\#61](https://github.com/Materials-Consortia/optimade-python-tools/pull/61) ([tachyontraveler](https://github.com/tachyontraveler))
- Some minor fixes [\#60](https://github.com/Materials-Consortia/optimade-python-tools/pull/60) ([ml-evs](https://github.com/ml-evs))
- Added codecov to CI [\#59](https://github.com/Materials-Consortia/optimade-python-tools/pull/59) ([ml-evs](https://github.com/ml-evs))
- Enforce black via `pre-commit` tool [\#53](https://github.com/Materials-Consortia/optimade-python-tools/pull/53) ([dwinston](https://github.com/dwinston))
- Update setup.py and version [\#51](https://github.com/Materials-Consortia/optimade-python-tools/pull/51) ([dwinston](https://github.com/dwinston))
- /structure/info endpoint [\#49](https://github.com/Materials-Consortia/optimade-python-tools/pull/49) ([fawzi](https://github.com/fawzi))
- add constrained list type [\#48](https://github.com/Materials-Consortia/optimade-python-tools/pull/48) ([dwinston](https://github.com/dwinston))
- Refactored into submodules and added test data [\#47](https://github.com/Materials-Consortia/optimade-python-tools/pull/47) ([ml-evs](https://github.com/ml-evs))
- Update structure endpoint to pre-alpha 0.10 spec [\#45](https://github.com/Materials-Consortia/optimade-python-tools/pull/45) ([ltalirz](https://github.com/ltalirz))
- Adding Resource Links [\#44](https://github.com/Materials-Consortia/optimade-python-tools/pull/44) ([tpurcell90](https://github.com/tpurcell90))
- Reblacken [\#42](https://github.com/Materials-Consortia/optimade-python-tools/pull/42) ([ml-evs](https://github.com/ml-evs))
- Documented json [\#41](https://github.com/Materials-Consortia/optimade-python-tools/pull/41) ([tpurcell90](https://github.com/tpurcell90))
- fix example output [\#40](https://github.com/Materials-Consortia/optimade-python-tools/pull/40) ([dwinston](https://github.com/dwinston))
- use jsonapi better at top level, add error response [\#36](https://github.com/Materials-Consortia/optimade-python-tools/pull/36) ([fawzi](https://github.com/fawzi))
- add JSONTransformer [\#35](https://github.com/Materials-Consortia/optimade-python-tools/pull/35) ([dwinston](https://github.com/dwinston))
- switch to MIT license [\#34](https://github.com/Materials-Consortia/optimade-python-tools/pull/34) ([ltalirz](https://github.com/ltalirz))
- Updated entry definitions and renamed Response classes [\#32](https://github.com/Materials-Consortia/optimade-python-tools/pull/32) ([ml-evs](https://github.com/ml-evs))
- update readme [\#31](https://github.com/Materials-Consortia/optimade-python-tools/pull/31) ([ltalirz](https://github.com/ltalirz))
- Seperated Links from JSON API into its own file [\#30](https://github.com/Materials-Consortia/optimade-python-tools/pull/30) ([tpurcell90](https://github.com/tpurcell90))
- simplify schema update [\#27](https://github.com/Materials-Consortia/optimade-python-tools/pull/27) ([ltalirz](https://github.com/ltalirz))
- add openapi\_diff to travis [\#25](https://github.com/Materials-Consortia/optimade-python-tools/pull/25) ([ltalirz](https://github.com/ltalirz))
- Json api add [\#24](https://github.com/Materials-Consortia/optimade-python-tools/pull/24) ([tpurcell90](https://github.com/tpurcell90))
- Added JSON diff test [\#20](https://github.com/Materials-Consortia/optimade-python-tools/pull/20) ([ml-evs](https://github.com/ml-evs))
- info endpoint [\#19](https://github.com/Materials-Consortia/optimade-python-tools/pull/19) ([fawzi](https://github.com/fawzi))
- adding run.sh script to start webserver [\#18](https://github.com/Materials-Consortia/optimade-python-tools/pull/18) ([fawzi](https://github.com/fawzi))
- error response [\#17](https://github.com/Materials-Consortia/optimade-python-tools/pull/17) ([fawzi](https://github.com/fawzi))
- Links can be strings [\#16](https://github.com/Materials-Consortia/optimade-python-tools/pull/16) ([fawzi](https://github.com/fawzi))
- response should be either many \(list\) or one \(object\), not an union [\#15](https://github.com/Materials-Consortia/optimade-python-tools/pull/15) ([fawzi](https://github.com/fawzi))
- reorg models [\#14](https://github.com/Materials-Consortia/optimade-python-tools/pull/14) ([dwinston](https://github.com/dwinston))
- Update the OptimadeMetaResponse to development schema [\#13](https://github.com/Materials-Consortia/optimade-python-tools/pull/13) ([ml-evs](https://github.com/ml-evs))
- add openapi spec validator [\#10](https://github.com/Materials-Consortia/optimade-python-tools/pull/10) ([ltalirz](https://github.com/ltalirz))
- fix test data download [\#3](https://github.com/Materials-Consortia/optimade-python-tools/pull/3) ([ltalirz](https://github.com/ltalirz))
- \[WIP\] Mongoconverter [\#1](https://github.com/Materials-Consortia/optimade-python-tools/pull/1) ([wuxiaohua1011](https://github.com/wuxiaohua1011))

## [v0.1.2](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.2) (2018-06-14)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.1...v0.1.2)

## [v0.1.1](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.1) (2018-06-13)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/v0.1.0...v0.1.1)

## [v0.1.0](https://github.com/Materials-Consortia/optimade-python-tools/tree/v0.1.0) (2018-06-05)

[Full Changelog](https://github.com/Materials-Consortia/optimade-python-tools/compare/6680fb28a60ec4ff43a303b7b4dbf41e159a25b6...v0.1.0)



\* *This Changelog was automatically generated by [github_changelog_generator](https://github.com/github-changelog-generator/github-changelog-generator)*
# Contributing and getting help

If you run into any problems using this package, or if you have a question, suggestion or feedback, then please raise an issue on [GitHub](https://github.com/Materials-Consortia/optimade-python-tools/issues/new).

The [Materials Consortia](https://github.com/Materials-Consortia) is very open to contributions across all of its packages.
This may be anything from simple feedback and raising [new issues](https://github.com/Materials-Consortia/optimade-python-tools/issues/new) to creating [new PRs](https://github.com/Materials-Consortia/optimade-python-tools/compare).

If you are interested in contributing but don't know where to begin, some issues have been marked with the [good first issue](https://github.com/Materials-Consortia/optimade-python-tools/labels/good%20first%20issue) label, typically where an isolated enhancement has a concrete suggestion.
Simply add a comment under an issue if you are interested in tackling it!

Recommendations for setting up a development environment for this package can be found in the [Installation instructions](https://www.optimade.org/optimade-python-tools/INSTALL/#full-development-installation).

More broadly, if you would like to ask questions or contact the consortium about creating an OPTIMADE implementation for a new database, then please read the relevant "get involved" section on the [OPTIMADE website](https://www.optimade.org/#get-involved).
# Configuration

Since the server implementation is built with [FastAPI](https://fastapi.tiangolo.com/), which uses [pydantic](https://pydantic-docs.helpmanual.io/), the configuration is based on pydantic's [Setting management](https://pydantic-docs.helpmanual.io/usage/settings/).
This way of handling configuration options supports various different approaches to configure the server.
We recommend either or a combination of the following:

1. Create a JSON or YAML configuration file with an implementation's complete configuration in the default location [DEFAULT_CONFIG_FILE_PATH][optimade.server.config.DEFAULT_CONFIG_FILE_PATH] or specify its location with the `OPTIMADE_CONFIG_FILE` environment variable.
2. Set environment variables prefixed with `OPTIMADE_` or `optimade_`.
3. Create a custom [`ServerConfig`][optimade.server.config.ServerConfig] object with the desired settings directly.
4. Load settings from a secret file (see [pydantic documentation](https://pydantic-docs.helpmanual.io/usage/settings/#secret-support) for more information).

## The JSON configuration file

The main way of configuring the OPTIMADE server is by creating a configuration JSON file.

An example of one that works with the example implementation can be found in [`optimade_config.json`](static/optimade_config.json):

=== "Configuration file for the default OPTIMADE server"
    ```json
    --8<-- "optimade_config.json"
    ```

## Environment variables

In order for the implementation to know where your configuration JSON file is located, you can set an environment variable `OPTIMADE_CONFIG_FILE` with either the value of the _absolute_ path to the configuration file or the _relative_ path to the file from the current working directory of where the server is run.

This variable is actually an extension of the configuration option `config_file`.
By default, the server will try to load a JSON file called `.optimade.json` located in your home folder (or equivalent).

Here the generally recognized environment variable prefix becomes evident, namely `OPTIMADE_` or `optimade_`.
Hence, you can set (or overwrite) any configuration option from the server's defaults or a value read from the configuration JSON by setting an environment variable named `OPTIMADE_<configuration_option>`.

## Custom configuration options

One can extend the current list of configuration options by sub-classing [`ServerConfig`][optimade.server.config.ServerConfig] and adding configuration options as attributes with values of `Field` (`pydantic.field`).
Any attribute type will be validated through `pydantic` as is the case for all of the regular configuration options.

This is useful for, e.g., custom database backends, if one wants to utilize the general server configuration setup implemented in `optimade` to declare specific database information.
It can also be useful if one wishes to extend and build upon the general `optimade` server with new endpoints and routes.

Remember to instantiate an instance of the sub-class, which can be imported and used in your application.

## List of configuration options

See [`config.py`][optimade.server.config.ServerConfig] for a complete list of configuration options.

The following configuration file represents the default values for all configuration options:

=== "Default values for all configuration options"
    ```json
    --8<-- "docs/static/default_config.json"
    ```
# Installation

This package can be installed from PyPI, or by cloning the repository, depending on your use-case.

1. To use the `optimade` Python package as a library, (e.g., using the models for validation, parsing filters with the grammar, or using the command-line tool `optimade-validator` tool), it is recommended that you install the latest release of the package from PyPI with `pip install optimade`.
2. If you want to run, use or modify the reference server implementation, then it is recommended that you clone this repository and install it from your local files (with `pip install .`, or `pip install -e .` for an editable installation).

## The index meta-database

This package may be used to setup and run an [OPTIMADE index meta-database](https://github.com/Materials-Consortia/OPTIMADE/blob/develop/optimade.rst#index-meta-database).
Clone this repository and install the package locally with `pip install -e .[server]`.

There is a built-in index meta-database set up to populate a `mongomock` in-memory database with resources from a static `json` file containing the `child` resources you, as a database provider, want to serve under this index meta-database.
The location of that `json` file is controllable using the `index_links_path` property of the configuration or setting via the environment variable `optimade_index_links_path`.

Running the index meta-database is then as simple as writing `./run.sh index` in a terminal from the root of this package.
You can find it at the base URL: <http://localhost:5001/v1>.

Here is an example of how it may look to start your server:

```sh
:~$ export OPTIMADE_CONFIG_FILE=/home/optimade_server/config.json
:~$ ./path/to/optimade/run.sh index
```

## Full development installation

The dependencies of this package can be found in `setup.py` with their latest supported versions.
By default, a minimal set of requirements are installed to work with the filter language and the `pydantic` models.
After cloning the repository, the install mode `server` (i.e. `pip install .[server]`) is sufficient to run a `uvicorn` server using the `mongomock` backend (or MongoDB with `pymongo`, if present).
The suite of development and testing tools are installed with via the install modes `dev` and `testing`.
There are additionally two backend-specific install modes, `elastic` and `mongo`, as well as the `all` mode, which installs all dependencies.
All contributed Python code, must use the [black](https://github.com/ambv/black) code formatter, and must pass the [flake8](http://flake8.pycqa.org/en/latest/) linter that is run automatically on all PRs.

```sh
# Clone this repository to your computer
git clone git@github.com:Materials-Consortia/optimade-python-tools.git
cd optimade-python-tools

# Ensure a Python>=3.7 (virtual) environment (example below using Anaconda/Miniconda)
conda create -n optimade python=3.7
conda activate optimade

# Install package and dependencies in editable mode (including "dev" requirements).
pip install -e ".[dev]"

# Optional: Retrieve the list of OPTIMADE providers. (Without this submodule, some of the tests will fail because "providers.json" cannot be found.)
git submodule update --init

# Run the tests with pytest
py.test

# Install pre-commit environment (e.g., auto-formats code on `git commit`)
pre-commit install

# Optional: Install MongoDB (and set `database_backend = mongodb`)
# Below method installs in conda environment and
# - starts server in background
# - ensures and uses ~/dbdata directory to store data
conda install -c anaconda mongodb
mkdir -p ~/dbdata && mongod --dbpath ~/dbdata --syslog --fork

# Start a development server (auto-reload on file changes at http://localhost:5000
# You can also execute ./run.sh
uvicorn optimade.server.main:app --reload --port 5000

# View auto-generated docs
open http://localhost:5000/docs
# View Open API Schema
open http://localhost:5000/openapi.json
```

When developing, you can run both the server and an index meta-database server at the same time (from two separate terminals).
Running the following:

```shell
./run.sh index
# or
uvicorn optimade.server.main_index:app --reload --port 5001
```

will run the index meta-database server at <http://localhost:5001/v1>.

## Testing specific backends

In order to run the test suite for a specific backend, the
`OPTIMADE_DATABASE_BACKEND` [environment variable (or config
option)](https://www.optimade.org/optimade-python-tools/configuration/) can be
set to one of `'mongodb'`, `'mongomock'` or `'elastic'` (see
[`ServerConfig.database_backend`][optimade.server.config.ServerConfig.database_backend]).
Tests for the two "real" database backends, MongoDB and Elasticsearch, require a writable, temporary database to be accessible.

The easiest way to deploy these databases and run the tests is with Docker, as shown below.
[Docker installation instructions](https://docs.docker.com/engine/install/) will depend on your system; on Linux, the `docker` commands below may need to be prepended with `sudo`, depending on your distribution.
These commands should be run from a local optimade-python-tools directory.

The following command starts a local Elasticsearch v6 instance, runs the test suite, then stops and deletes the containers (required as the tests insert some data):
```shell
docker run -d --name elasticsearch_test -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" elasticsearch:6.8.22 \
&& sleep 10 \
&& OPTIMADE_DATABASE_BACKEND="elastic" py.test; \
docker container stop elasticsearch_test; docker container rm elasticsearch_test
```

The following command starts a local MongoDB instance, runs the test suite, then stops and deletes the containers:
```shell
docker run -d --name mongo_test -p 27017:27017 -d mongo:4.4.6 \
&& OPTIMADE_DATABASE_BACKEND="mongodb" py.test; \
docker container stop mongo_test; docker container rm mongo_test
```
# warnings

::: optimade.adapters.warnings
# exceptions

::: optimade.adapters.exceptions
# logger

::: optimade.adapters.logger
# base

::: optimade.adapters.base
# adapter

::: optimade.adapters.references.adapter
# adapter

::: optimade.adapters.structures.adapter
# utils

::: optimade.adapters.structures.utils
# jarvis

::: optimade.adapters.structures.jarvis
# aiida

::: optimade.adapters.structures.aiida
# pymatgen

::: optimade.adapters.structures.pymatgen
# proteindatabank

::: optimade.adapters.structures.proteindatabank
# cif

::: optimade.adapters.structures.cif
# ase

::: optimade.adapters.structures.ase
# utils

::: optimade.validator.utils
# validator

::: optimade.validator.validator
# config

::: optimade.validator.config
# lark_parser

::: optimade.filterparser.lark_parser
# warnings

::: optimade.server.warnings
# exception_handlers

::: optimade.server.exception_handlers
# query_params

::: optimade.server.query_params
# middleware

::: optimade.server.middleware
# exceptions

::: optimade.server.exceptions
# config

::: optimade.server.config
# main_index

::: optimade.server.main_index
# logger

::: optimade.server.logger
# main

::: optimade.server.main
# schemas

::: optimade.server.schemas
# structures

::: optimade.server.mappers.structures
# entries

::: optimade.server.mappers.entries
# references

::: optimade.server.mappers.references
# links

::: optimade.server.mappers.links
# structures

::: optimade.server.routers.structures
# references

::: optimade.server.routers.references
# utils

::: optimade.server.routers.utils
# links

::: optimade.server.routers.links
# landing

::: optimade.server.routers.landing
# versions

::: optimade.server.routers.versions
# index_info

::: optimade.server.routers.index_info
# info

::: optimade.server.routers.info
# elasticsearch

::: optimade.server.entry_collections.elasticsearch
# entry_collections

::: optimade.server.entry_collections.entry_collections
# mongo

::: optimade.server.entry_collections.mongo
# baseinfo

::: optimade.models.baseinfo
    rendering:
      show_if_no_docstring: true
# structures

::: optimade.models.structures
    rendering:
      show_if_no_docstring: true
# entries

::: optimade.models.entries
    rendering:
      show_if_no_docstring: true
# references

::: optimade.models.references
    rendering:
      show_if_no_docstring: true
# utils

::: optimade.models.utils
    rendering:
      show_if_no_docstring: true
# links

::: optimade.models.links
    rendering:
      show_if_no_docstring: true
# index_metadb

::: optimade.models.index_metadb
    rendering:
      show_if_no_docstring: true
# jsonapi

::: optimade.models.jsonapi
    rendering:
      show_if_no_docstring: true
# responses

::: optimade.models.responses
    rendering:
      show_if_no_docstring: true
# optimade_json

::: optimade.models.optimade_json
    rendering:
      show_if_no_docstring: true
# elasticsearch

::: optimade.filtertransformers.elasticsearch
# mongo

::: optimade.filtertransformers.mongo
# base_transformer

::: optimade.filtertransformers.base_transformer
# Validation of OPTIMADE APIs

`optimade-python-tools` contains tools for validating external OPTIMADE implementations that may be helpful for all OPTIMADE providers.
The validator is dynamic and fuzzy, in that the tested filters are generated based on *random* entries served by the API, and the description of the API provided at the `/info` endpoint.

The validator is implemented in the [`optimade.validator`][optimade.validator.validator] submodule, but the two main entry points are:

1. The `optimade-validator` script, which is installed alongside the package.
2. The [`optimade-validator-action`](https://github.com/Materials-Consortia/optimade-validator-action) which allows the validator to be used as a GitHub Action.

To run the script, simply provide an OPTIMADE URL to the script at the command-line.
You can use the following to validate the Heroku deployment of our reference server:

```shell
$ optimade-validator https://optimade.herokuapp.com/
```

Several additional options can be found under the `--help` flag, with the most important being `-v/-vvvv` to set the verbosity, `--index` to validate OPTIMADE index meta-databases and `--json` to receive the validation results as JSON document for programmatic use.

```shell
$ optimade-validator --help
usage: optimade-validator [-h] [-v] [-j] [-t AS_TYPE] [--index]
                          [--skip-optional] [--fail-fast] [-m]
                          [--page_limit PAGE_LIMIT]
                          [--headers HEADERS]
                          [base_url]

Tests OPTIMADE implementations for compliance with the optimade-python-tools models.

- To test an entire implementation (at say example.com/optimade/v1) for all required/available endpoints:

    $ optimade-validator http://example.com/optimade/v1

- To test a particular response of an implementation against a particular model:

    $ optimade-validator http://example.com/optimade/v1/structures/id=1234 --as-type structure

- To test a particular response of an implementation against a particular model:

    $ optimade-validator http://example.com/optimade/v1/structures --as-type structures


positional arguments:
  base_url              The base URL of the OPTIMADE
                        implementation to point at, e.g.
                        'http://example.com/optimade/v1' or
                        'http://localhost:5000/v1'

optional arguments:
  -h, --help            show this help message and exit
  -v, --verbosity       Increase the verbosity of the output.
                        (-v: warning, -vv: info, -vvv: debug)
  -j, --json            Only a JSON summary of the validator
                        results will be printed to stdout.
  -t AS_TYPE, --as-type AS_TYPE
                        Validate the request URL with the
                        provided type, rather than scanning the
                        entire implementation e.g. optimade-
                        validator `http://example.com/optimade/v1
                        /structures/0 --as-type structure`
  --index               Flag for whether the specified OPTIMADE
                        implementation is an Index meta-database
                        or not.
  --skip-optional       Flag for whether the skip the tests of
                        optional features.
  --fail-fast           Whether to exit on first test failure.
  -m, --minimal         Run only a minimal test set.
  --page_limit PAGE_LIMIT
                        Alter the requested page limit for some
                        tests.
  --headers HEADERS     Additional HTTP headers to use for each
                        request, specified as a JSON object.
```
# Example use cases

## Serving a single database

The [Materials Project](https://materialsproject.org) uses `optimade-python-tools` alongside their existing API and MongoDB database, providing [OPTIMADE-compliant access](https://optimade.materialsproject.org) to highly-curated density-functional theory calculations across all known inorganic materials.

`optimade-python-tools` handles filter parsing, database query generation and response validation by running the reference server implementation with minimal configuration.

[*odbx*](https://odbx.science), a small database of results from crystal structure prediction calculations, follows a similar approach.
This implementation is open source, available on GitHub at [ml-evs/odbx.science](https://github.com/ml-evs/odbx.science).

## Serving multiple databases

[Materials Cloud](https://materialscloud.org) uses `optimade-python-tools` as a library to provide an OPTIMADE API entry to archived computational materials studies, created with the [AiiDA](https://aiida.net) Python framework and published through their archive.
In this case, each individual study and archive entry has its own database and separate API entry.
The Python classes within the `optimade` package have been extended to make use of AiiDA and its underlying [PostgreSQL](https://postgresql.org) storage engine.

Details of this implementation can be found on GitHub at [aiidateam/aiida-optimade](https://github.com/aiidateam/aiida-optimade).

## Extending an existing API

[NOMAD](https://nomad-lab.eu/) uses `optimade-python-tools` as a library to add OPTIMADE API endpoints to an existing web app.
Their implementation uses the Elasticsearch database backend to filter on millions of structures from aggregated first-principles calculations provided by their users and partners.
NOMAD also uses the package to implement a GUI search bar that accepts the OPTIMADE filter language.
NOMAD uses the release versions of the `optimade-python-tools` package, performing all customisation via configuration and sub-classing.
The NOMAD OPTIMADE API implementation is available in the [NOMAD FAIR GitLab repository](https://gitlab.mpcdf.mpg.de/nomad-lab/nomad-FAIR).

This use case is demonstrated in the example [Integrate OPTIMADE with an existing web application](integrated.md).
# Setting up an OPTIMADE API

These notes describe how to set up and customize an OPTIMADE API based on the reference server in this package for some existing crystal structure data.

To follow this guide, you will need to have a working development installation, as described in the [installation instructions](../INSTALL.md#full-development-installation).
Complete examples of APIs that use this package are described in the [Use Cases](./use_cases.md) section.

## Setting up the database

The `optimade` reference server requires a data source per OPTIMADE entry type (`structures`, `references`, `links`).
In the simplest case, these can be configured as named MongoDB collections with a defined MongoDB URI and database name (see below), but they can also be set up as custom subclasses of [`EntryCollection`][optimade.server.entry_collections.entry_collections.EntryCollection] that could simply read from a static file.
In the reference server, these data sources, or collections, are created in the submodule for the corresponding routers/endpoints.

Here, we shall use the built-in MongoDB collections for each entry type, by simply specifying the appropriate options in the [configuration](../configuration.md), namely [`"database_backend": "mongodb"`][optimade.server.config.ServerConfig.database_backend], [`"mongo_uri": "mongodb://localhost:27017"`][optimade.server.config.ServerConfig.mongo_uri], [`"mongo_database": "optimade"`][optimade.server.config.ServerConfig.mongo_database] and the collection names for each entry type ([`"structures_collection": "structures"`][optimade.server.config.ServerConfig.structures_collection] etc.).
These notes will now assume that you have a MongoDB instance running and you have created a database that matches your [`"mongo_database"`][optimade.server.config.ServerConfig.mongo_database] config option.

If you disable inserting test data (with the [`"insert_test_data": false`][optimade.server.config.ServerConfig.insert_test_data] configuration option), you can test your API/database connection by running the web server with `uvicorn optimade.server.main:app --port 5000` and visiting the (hopefully empty) structures endpoint at `localhost:5000/v1/structures` (or your chosen base URL).

!!! note
    As of version v0.16, the other supported database backend is Elasticsearch.
    If you are interested in using another backend, or would like it to be supported in the `optimade` package, please raise an issue on [GitHub](https://github.com/Materials-Consortia/optimade-python-tools/issues/new) and visit the notes on implementing new [filter transformers](./filtering.md#developing-new-filter-transformers).

## Mapping non-OPTIMADE data

There are two ways to work with data that does not exactly match the OPTIMADE specification, both of which require configuring a subclass of [`BaseResourceMapper`][optimade.server.mappers.entries.BaseResourceMapper] that converts your stored data format into an OPTIMADE-compliant entry.
The two options are:

- Use the mapper to dynamically convert the data stored in the database, and the filters on that data, to an OPTIMADE format when responding to API requests.
- Apply the mapper to your entries before ingestion and use it to create a secondary database that stores the converted entries (e.g., normalized data), or equivalently, adding all the required OPTIMADE fields inside the existing entries (e.g., denormalized data)

The main consideration when choosing these options is not necessarily how closely your data matches the OPTIMADE format, but instead how readily the OPTIMADE filtering of that document can be mapped into the corresponding database query.
This could require writing or extending the [`BaseFilterTransformer`][optimade.filtertransformers.base_transformer.BaseTransformer] class, which takes an OPTIMADE filter string and converts it into a backend-specific query.

For example, if your database stores chemical formulae with extraneous "1"'s, e.g., SiO<sub>2</sub> is represented as `"Si1O2"`, then the incoming OPTIMADE filter (which asserts that elements must be alphabetical, and "1"'s must be omitted) for `chemical_formula_reduced="O2Si"` will also need to be transformed so that the corresponding database query matches the stored string, which in this case can be done easily.
Instead, if you are storing chemical formulae as an unreduced count per simulation cell, e.g., `"Si4O8"`, then it is impossible to remap the filter `chemical_formula_reduced="O2Si"` such that it matches all structures with the correct formula unit (e.g., `"SiO2"`, `"Si2O4"`, ...).
This would then instead require option 2 above, namely either the addition of auxiliary fields that store the correct (or mappable) OPTIMADE format in the database, or the creation of a secondary database that returns the pre-converted structures.

In the simplest case, the mapper classes can be used to define aliases between fields in the database and the OPTIMADE field name; these can be configured via the [`aliases`][optimade.server.config.ServerConfig.aliases] option as a dictionary mapping stored in a dictionary under the appropriate endpoint name, e.g. `"aliases": {"structures": {"chemical_formula_reduced": "my_chem_form"}}`, or defined as part of a custom mapper class.

In either option, you should now be able to insert your data into the corresponding MongoDB (or otherwise) collection.

## Serving custom fields/properties

According to the OPTIMADE specification, any field not standardized in the specification must be prefixed with an appropriate "provider prefix" (e.g., "`_aflow`" for [AFLOW](https://aflow.org) and "`_cod`" for [COD](https://crystallography.net)).
This prefix is intended to be unique across all [OPTIMADE providers](https://github.com/Materials-Consortia/providers) to enable filters to work across different implementations.
The prefix can be set in the [configuration](../configuration.md) as part of the [`provider`][optimade.server.config.ServerConfig.provider] option.

Once the prefix has been set, custom fields can be listed by endpoint in the [`provider_fields`][optimade.server.config.ServerConfig.provider_fields] configuration option.
Filters that use the prefixed form of these fields will then be passed through to the underlying database without the prefix, and then the prefix will be reinstated in the response.

!!! warning
    This config-only approach does not provide any way of **describing** the underlying field (via `description`), its type, or any potential physical units, and the field will not be added to the corresponding entry info endpoint (e.g., `/info/structures`).
    For this, you will need to follow the more complicated method below, under [More advanced usage](#more_advanced_usage).

### More advanced usage

It is recommended that you provide a description, type and unit for each custom field that can be returned at the corresponding `/info/<entry_type>` endpoint.
To do this, the underlying `EntryResourceAttributes` model will need to be sub-classed, the pydantic fields added to that class, and the server adjusted to make use of those models in responses.
In this case, it may be easier to write a custom endpoint for your entry type, that copies the existing reference endpoint.

Your custom model will need to be registered in three places:

1. The data collection.
1. The resource mapper class used by the collection.
1. The [`ENTRY_INFO_SCHEMAS`][optimade.server.schemas.ENTRY_INFO_SCHEMAS] dictionary.

Finally, the model must be instructed to use the prefixed (aliased) fields when generating its schemas.

Pulling all of this together:

```python
from optimade.server.schemas import ENTRY_INFO_SCHEMAS
from optimade.models import (
    StructureResource, StructureResourceAttributes, OptimadeField
)


class MyStructureResourceAttributes(StructureResourceAttributes):
     my_custom_field: str = OptimadeField(
        "default value",
        description="This is a custom field",
    )

    class Config:
        """Add a pydantic `Config` that defines the alias generator,
        based on our configured `provider_fields`.

        """
        @classmethod
        def alias_generator(cls, name: str) -> str:
            if name in CONFIG.provider_fields.get("structures", []):
                return f"_{CONFIG.provider.prefix}_{name}"
            return name


class MyStructureResource(StructureResource):
    attributes: MyStructureResourceAttributes


ENTRY_INFO_SCHEMAS["structures"] = MyStructureResource.schema
```

Currently, the reference server is not flexible enough to use custom response classes via configuration only (there is an open issue tracking this [#929](https://github.com/Materials-Consortia/optimade-python-tools/issues/929)), so instead the code will need to be forked and modified for your implementation.

## Validating your implementation

With the database collections, mappers, aliases and provider configured, you can try running the web server (with e.g., `uvicorn optimade.server.main:app`, if your app is in the same file as the reference server) and validating it as an OPTIMADE API, following the [validation guide](./validation.md).

## Registering as a provider

If you host your API at a persistent URL, you should consider registering as an OPTIMADE provider, which will add you to the federated list used by users and clients to discover data.
Instructions for how to do this can be found at in the [Materials-Consortia/providers](https://github.com/Materials-Consortia/providers) repository.
# Filter parsing and transforming

One of the aims of this package is to integrate with existing databases and APIs, and as such your particular backend may not have a supported filter transformer.
This guide will briefly outline how to parse OPTIMADE filter strings into database or API-specific queries.

## Parsing OPTIMADE filter strings

The [`LarkParser`][optimade.filterparser.lark_parser.LarkParser] class will take an OPTIMADE filter string, supplied by the user, and parse it into a `lark.Tree` instance.

Example use:

```python
from optimade.filterparser import LarkParser

p = LarkParser(version=(1, 0, 0))
tree = p.parse("nelements<3")
print(tree)
```

```shell
Tree('filter', [Tree('expression', [Tree('expression_clause', [Tree('expression_phrase', [Tree('comparison', [Tree('property_first_comparison', [Tree('property', [Token('IDENTIFIER', 'nelements')]), Tree('value_op_rhs', [Token('OPERATOR', '<'), Tree('value', [Tree('number', [Token('SIGNED_INT', '3')])])])])])])])])])
```

```python
print(tree.pretty())
```

```shell
filter
  expression
    expression_clause
      expression_phrase
        comparison
          property_first_comparison
            property	nelements
            value_op_rhs
              <
              value
                number	3
```

```python
tree = p.parse('_mp_bandgap > 5.0 AND _cod_molecular_weight < 350')
print(tree.pretty())
```

```shell
filter
  expression
    expression_clause
      expression_phrase
        comparison
          property_first_comparison
            property	_mp_bandgap
            value_op_rhs
              >
              value
                number	5.0
      expression_phrase
        comparison
          property_first_comparison
            property	_cod_molecular_weight
            value_op_rhs
              <
              value
                number	350
```

## Flow for parsing a user-supplied filter and converting to a backend query

After the [`LarkParser`][optimade.filterparser.lark_parser.LarkParser] has turned the filter string into a `lark.Tree`, it is fed to a `lark.Transformer` instance, which transforms the 'lark.Tree' into a backend-specific representation of the query.
For example, [`MongoTransformer`][optimade.filtertransformers.mongo.MongoTransformer] will turn the tree into something useful for a MongoDB backend:

```python
# Example: Converting to MongoDB Query Syntax
from optimade.filtertransformers.mongo import MongoTransformer

transformer = MongoTransformer()

tree = p.parse('_mp_bandgap > 5.0 AND _cod_molecular_weight < 350')
query = transformer.transform(tree)
print(query)
```

```json
{
    "$and": [
        {"_mp_bandgap": {"$gt": 5.0}},
        {"_cod_molecular_weight": {"$lt": 350.0}}
    ]
}
```


## Developing new filter transformers

In order to support a new backend, you will need to create a new filter transformer that inherits from the [`BaseTransformer`][optimade.filtertransformers.base_transformer.BaseTransformer].
This transformer will need to override the methods that match the particular grammatical constructs in the Lark grammar in order to construct a query.
Two examples can be found within `optimade-python-tools`, one for MongoDB ([`MongoTransformer`][optimade.filtertransformers.mongo.MongoTransformer]) and one for Elasticsearch ([`ElasticTransformer`][optimade.filtertransformers.elasticsearch.ElasticTransformer]).

In some cases, you may also need to extend the base [`EntryCollection`][optimade.server.entry_collections.entry_collections.EntryCollection], the class that receives the transformed filter as an argument to its private `._run_db_query()` method.
This class handles the connections to the underlying database, formatting of the response in an OPTIMADE format, and other API features such as sorting and pagination.
Again, the examples for MongoDB ([`MongoCollection`][optimade.server.entry_collections.mongo.MongoCollection]) and Elasticsearch ([`ElasticCollection`][optimade.server.entry_collections.elasticsearch.ElasticCollection]) should be helpful.

If you would like to contribute your new filter transformer back to the package, please raise an issue to signal your intent (in case someone else is already working on this).
Adding a transformer requires the following:

1. A new submodule (`.py` file) in the `optimade/filtertransformers` folder containing an implementation of the transformer object that extends `optimade.filtertransformers.base_transformer.BaseTransformer`.
2. Any additional Python requirements must be optional and provided as a separate "`extra_requires`" entry in `setup.py` and in the `requirements.txt` file.
3. Tests in `optimade/filtertransformers/tests` that are skipped if the required packages fail to import.
# Integrate OPTIMADE with an existing web application

The `optimade` package can be used to create a standalone web application that serves the OPTIMADE API based on a pre-configured MongoDB backend.
In this document, we are going to use `optimade` differently and use it to add an OPTIMADE API implementation alongside an existing API that employs an Elasticsearch storage layer.

Let's assume we already have a *FastAPI* application that runs an unrelated web service, and that we use an Elasticsearch backend that contains all structure data, but not necessarily in a form that OPTIMADE expects.

## Providing the `optimade` configuration

`optimade` can read its configuration from a JSON file.
It uses the `OPTIMADE_CONFIG_FILE` environment variable (or a default path) to find the config file.
If you run `optimade` code inside another application, you might want to provide this config file as part of the source code and not via environment variables.
Let's say you have a file `optimade_config.json` as part of the Python module that you use to create your OPTIMADE API.

!!! tip
    You can find more detailed information about configuring the `optimade` server in the [Configuration](../configuration.md) section.

Before importing any `optimade` modules, you can set the `OPTIMADE_CONFIG_FILE` environment variable to refer to your config file:

```python
import os
from pathlib import Path

os.environ['OPTIMADE_CONFIG_FILE'] = str(Path(__file__).parent / "optimade_config.json")
```

## Customize the [`EntryCollection`][optimade.server.entry_collections.entry_collections.EntryCollection] implementation

Let's assume that your Elasticsearch backend stores structure data in a different enough manner that you need to provide your own custom implementation.
The following code customizes the [`EntryCollection`][optimade.server.entry_collections.entry_collections.EntryCollection] class for structures, whilst keeping the default MongoDB-based implementation (using [`MongoCollection`][optimade.server.entry_collections.mongo.MongoCollection]) for all other entry types.

```python
from optimade.server.routers import structures

structures.structures_coll = MyElasticsearchStructureCollection()
```

You can imagine that `MyElasticsearchStructureCollection` either sub-classes the default `optimade` Elasticsearch implementation ([`ElasticsearchCollection`][optimade.server.entry_collections.elasticsearch.ElasticCollection]) or sub-classes [`EntryCollection`][optimade.server.entry_collections.entry_collections.EntryCollection], depending on how deeply you need to customize the default `optimade` behavior.

## Mounting the OPTIMADE Python tools *FastAPI* app into an existing *FastAPI* app

Let's assume you have an existing *FastAPI* app `my_app`.
It already implements a few routers under certain path prefixes, and now you want to add an OPTIMADE implementation under the path prefix `/optimade`.

First, you have to set the `root_path` in the `optimade` configuration, so that the app expects all requests to be prefixed with `/optimade`.

Second, you simply mount the `optimade` app into your existing app `my_app`:

```python
from optimade.server.config import CONFIG

CONFIG.root_path = "/optimade"

from optimade.server import main as optimade

optimade.add_major_version_base_url(optimade.app)
my_app.mount("/optimade", optimade.app)
```

!!! tip
    In the example above, we imported `CONFIG` before `main` so that our config was loaded before app creation.
    To avoid the need for this, the `root_path` can be set in your JSON config file, passed as an environment variable, or declared in a custom Python module (see [Configuration](../configuration.md)).

See also the *FastAPI* documentation on [sub-applications](https://fastapi.tiangolo.com/advanced/sub-applications/).

Now, if you run `my_app`, it will still serve all its routers as before and in addition it will also serve all OPTIMADE routes under `/optimade/` and the versioned URLs `/optimade/v1/`.
