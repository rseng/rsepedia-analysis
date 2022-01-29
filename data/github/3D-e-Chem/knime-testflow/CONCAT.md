# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
The file is formatted as described on http://keepachangelog.com/.

## [Unreleased]

### Changes

- Requires KNIME 4.0 [#7](https://github.com/3D-e-Chem/knime-testflow/issues/7)

## [1.0.2] - 2017-03-07

### Changes

- Requires KNIME 3.3

### Fixes

- Disabling check log message will leak log into next test (#5)
- NoSuchMethodError setMountpointRoot (#4)

## [1.0.1] - 2017-01-31

### Fixes

- Workaround for OSX SWT issues (#3)

## [1.0.0] - 2016-06-17

Initial release.

### Added

* Testing 2 workflows with several check permutations (#1)
Test Knime workflows from a Junit test.

[![Build Status](https://travis-ci.org/3D-e-Chem/knime-testflow.svg?branch=master)](https://travis-ci.org/3D-e-Chem/knime-testflow)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/ba09652161144d9abbe4827fd16bbaec)](https://www.codacy.com/app/3D-e-Chem/knime-testflow?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=3D-e-Chem/knime-testflow&amp;utm_campaign=Badge_Grade)
[![DOI](https://zenodo.org/badge/doi/10.5281/zenodo.55805.svg)](http://dx.doi.org/10.5281/zenodo.55805)

The Knime Testing Framework can run a test workflow either:
* Inside Knime, if you right-click on a workflow in your local workspace, you can select "Run as workflow test".
* From the command line, using `knime -application org.knime.testing.NGTestflowRunner -root <workflow dir>`.

This repo gives you another option run a test workflow inside of a Junit @Test method declaration.

This project uses [Eclipse Tycho](https://www.eclipse.org/tycho/) to perform build steps.

# Usage

Using the plugin requires several steps.

## 1. Add repository

This plugin is available in the `https://3d-e-chem.github.io/updates` update site.

To make use of in a Tycho based project add to the `<repositories>` tag of the `pom.xml` file the following:
```
<repository>
    <id>3d-e-chem</id>
    <layout>p2</layout>
    <url>https://3d-e-chem.github.io/updates</url>
</repository>
```

## 2. Add dependency to tests

In the `Require-Bundle` attribute of the `META-INF/MANIFEST.MF` of the tests module add
```
nl.esciencecenter.e3dchem.knime.testing.plugin;bundle-version="[1.0.0,2.0.0)",
org.knime.testing;bundle-version="[4.0.0,5.0.0)",
```

## 3. Add test workflow

Create a test workflow as described in the "Testing Framework" manual that you get when you install the "KNIME Testing Framework" (look in plugins/org.knime.testing_x.y.z/doc/Regression Tests.pdf).

Place the workflow as a directory inside the `src/knime/` directory of the tests module.

## 4. Add test

Create a new test class and inside the class put the following:
```
@Rule
public ErrorCollector collector = new ErrorCollector();
private TestFlowRunner runner;

@Before
public void setUp() {
    TestrunConfiguration runConfiguration = new TestrunConfiguration();
    runner = new TestFlowRunner(collector, runConfiguration);
}

@Test
public void test_simple() throws IOException, InvalidSettingsException, CanceledExecutionException,
        UnsupportedWorkflowVersionException, LockFailedException, InterruptedException {
    File workflowDir = new File("src/knime/my-workflow-test");
    runner.runTestWorkflow(workflowDir);
}
```

This will test the workflow put in `src/knime/my-workflow-test` in the previous step.

This will run minimal checks, to check more configure `runConfiguration` object.  
For example add some more checks by adding 
```
runConfiguration.setTestDialogs(true);
runConfiguration.setReportDeprecatedNodes(true);
runConfiguration.setCheckMemoryLeaks(true);
```

## 5. Run tests

```
mvn verify
```

The test results can be found in the `T E S T S` section of the standard output.

## 6. Add GUI testing to Travis-CI.

As you might have noticed during the previouse step, running test will quickly show some dialogs and windows.
To show graphical user elements an X-server is required, sadly Travis-CI does not run an X-server. 
A temporary X-server can be run with Xvfb, which is luckily available on all Travis-CI environments.

Prepend `xvfb-run` before the `mvn verify` command in the `.travis.yml` file.

For example
```
script: xvfb-run mvn verify -B
```

# Build

```
mvn verify
```

An Eclipse update site will be made in `p2/target/repository` repository.
The update site can be used to perform a local installation.

# Development

Steps to get development environment setup based on https://github.com/knime/knime-sdk-setup#sdk-setup:

1. Install Java 8
2. Install Eclipse for [RCP and RAP developers](https://www.eclipse.org/downloads/packages/release/2018-12/r/eclipse-ide-rcp-and-rap-developers)
3. Configure Java 8 inside Eclipse Window > Preferences > Java > Installed JREs
4. Import this repo as an Existing Maven project
5. Activate target platform by going to Window > Preferences > Plug-in Development > Target Platform and check the `KNIME Analytics Platform (4.0) - nl.esciencecenter.e3dchem.knime.testing.targetplatform/KNIME-AP-4.0.target` target definition.

During import the Tycho Eclipse providers must be installed.

# New release

1. Update versions in pom files with `mvn org.eclipse.tycho:tycho-versions-plugin:set-version -DnewVersion=<version>-SNAPSHOT` command.
2. Commit and push changes
3. Create package with `mvn package`, will create update site in `p2/target/repository`
4. Append new release to an update site
  1. Make clone of an update site repo
  2. Append release to the update site with `mvn install -Dtarget.update.site=<path to update site>`
5. Commit and push changes in this repo and update site repo.

