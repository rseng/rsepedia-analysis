# Change Log

## [v0.3.0](https://github.com/recipy/recipy/tree/v0.3.0) (2016-09-13)
[Full Changelog](https://github.com/recipy/recipy/compare/v0.2.3...v0.3.0)

This is a major release with a number of new features and lots of tidying up. It has been released in time for PyConUK, where we will be sprinting on recipy and hopefully developing even more new features and fixing more bugs.

The major new features are:

 * A hash is now computed for all input and output files, and this is used for searching for files. Practically, this means that you can create `graph.png`, send it to a colleague, and even if they send back a file called `RobinsGraph.png` you can still look up how it was created, as the hash will be the same even if the filename isn't.
 * Hashes are computed by default, rather than hashing having to be switched on in the config file (Note: this is a change from the pre-release versions, as hashing is now working effectively).
 * `recipy search` has gained a new option: `--filepath` (or `-p`) which forces searching based on the file and path rather than the hash.
 * There is now a beta-quality `recipy.open` function that can be used instead of Python's standard `open` function to log files that are opened. (The reason why the standard Python `open` command can't be wrapped by recipy is a topic for a blog post, I think). There are some known issues with this function, but it should work for many use cases.
 * Versions of the libraries that are patched by recipy are now stored with the run
 * There is now a `recipy annotate` command that can be used to add notes to a 'recipy run'.
 * There is a new tab in the GUI that lists all of modules and functions that are patched by this version of recipy.
 * The terminal output from the `recipy` command is now far more easy to read, as it uses different font styles for headers and the actual content
 * The GUI has been updated to handle and display all of the new information we're storing.
 * Lots of bugs have been fixed (far more than listed in the issues below, as not all were recorded officially via issues) and the `recipy search` command is far more robust now.

 Huge thanks must also go to [mikej888](https://github.com/mikej888) who has been working on recipy through the [Software Sustainability Institute](http://www.software.ac.uk) and has reported many bugs and give lots of ideas for improvements. He is currently working on a test suite that should make the whole package far more robust.

**Implemented enhancements:**

- Store hash for input/output files [\#25](https://github.com/recipy/recipy/issues/25)
- Hook into standard `open` etc in pure python [\#44](https://github.com/recipy/recipy/issues/44)
- Store data diffs [\#107](https://github.com/recipy/recipy/issues/107)
- Format terminal output nicely [\#102](https://github.com/recipy/recipy/issues/102)
- How should GUI display hashes? [\#100](https://github.com/recipy/recipy/issues/100)
- Store versions of libraries used [\#87](https://github.com/recipy/recipy/issues/87)
- Use colours in terminal output of CLI [\#86](https://github.com/recipy/recipy/issues/86)
- Add 'annotate' functionality [\#69](https://github.com/recipy/recipy/issues/69)

**Fixed bugs:**

- Output to JSON with no results found fails with IndexError [\#142](https://github.com/recipy/recipy/issues/142)
- recipy command fails on Windows due to blessings use [\#141](https://github.com/recipy/recipy/issues/141)
- Issues with searching by ID [\#124](https://github.com/recipy/recipy/issues/124)
- latest\_run throws an error if there are no runs in the database [\#118](https://github.com/recipy/recipy/issues/118)
- tinydb error [\#92](https://github.com/recipy/recipy/issues/92)
- Recipy not able to find a file  [\#89](https://github.com/recipy/recipy/issues/89)
- Add blessings to setup.py [\#140](https://github.com/recipy/recipy/issues/140)
- Diffs can't be applied as patches because no filename info is included [\#133](https://github.com/recipy/recipy/issues/133)
- matplotlib isn't a dependency [\#115](https://github.com/recipy/recipy/issues/115)
- TypeError: 'Query' object is not callable when searching for runs in GUI [\#112](https://github.com/recipy/recipy/issues/112)

**Closed issues:**

- Split the GUI into a separate PyPI package? [\#117](https://github.com/recipy/recipy/issues/117)
- Use with R via doit? [\#114](https://github.com/recipy/recipy/issues/114)
- Source code not PEP-8 [\#110](https://github.com/recipy/recipy/issues/110)
- Python 3 compatibility: deprecated modules [\#109](https://github.com/recipy/recipy/issues/109)

**Merged pull requests:**

- adds check for $EDITOR for the annotate command [\#97](https://github.com/recipy/recipy/pull/97) ([cash](https://github.com/cash))
- updated readme to be consistent about location of recipy console script [\#96](https://github.com/recipy/recipy/pull/96) ([cash](https://github.com/cash))
- add test\_requirements.txt for installing unit test dependencies [\#94](https://github.com/recipy/recipy/pull/94) ([cash](https://github.com/cash))
- adds tinydb-serialization to requirements.txt for manually installing deps [\#93](https://github.com/recipy/recipy/pull/93) ([cash](https://github.com/cash))


## [v0.2.3](https://github.com/recipy/recipy/tree/HEAD)
[Full Changelog](https://github.com/recipy/recipy/compare/v0.2.3...HEAD)

This is a very small release to fix recipy to work with the latest version of TinyDB (v3.0.0).

**Fixed bugs:**

- TinyDB import error [\#92](https://github.com/recipy/recipy/issues/92)

## [v0.2.1](https://github.com/recipy/recipy/tree/v0.2.1)
[Full Changelog](https://github.com/recipy/recipy/compare/v0.2.0...HEAD)

Minor bug-fix release.

**Fixed bugs:**

- Patching of Pillow does not work due to [\#52](https://github.com/recipy/recipy/issues/52), which caused matplotlib imports to fail. Fixed by removing Pillow support for the moment

## [v0.2.0](https://github.com/recipy/recipy/tree/v0.2.0) (2015-09-21)
[Full Changelog](https://github.com/recipy/recipy/compare/v0.1.0...v0.2.0)

This is the first new release of recipy since its debut at EuroSciPy 2015. Sorry for the delay in getting this out, life has been rather chaotic for all of the members of the recipy team. We want to say a huge thank you to all of the people who have submitted issues and pull requests: we couldn't have done this without you!

Full details are below, but the major new features include logging command-line arguments, exporting runs to JSON, adding more configuration options, and adding more commands to the GUI. A number of bugs have also been fixed.

**Implemented enhancements:**

- Add export option to the GUI [\#68](https://github.com/recipy/recipy/issues/68)
- Add 'latest' option to CLI [\#59](https://github.com/recipy/recipy/issues/59)
- Make DB path configurable [\#57](https://github.com/recipy/recipy/issues/57)
- Configuration option not to track input or output files [\#56](https://github.com/recipy/recipy/issues/56)
- Add search by id [\#49](https://github.com/recipy/recipy/issues/49)
- Configuration file should be ~/.recipy/recipyrc [\#46](https://github.com/recipy/recipy/issues/46)
- Add 'quiet' option to config file to stop display of 'recipy run inserted xyz' message [\#43](https://github.com/recipy/recipy/issues/43)
- Add export of individual runs, or the whole database [\#50](https://github.com/recipy/recipy/issues/50)
- Log command-line arguments [\#47](https://github.com/recipy/recipy/issues/47)

**Fixed bugs:**

- Loading a numpy file gives error message, only after recipy import. [\#83](https://github.com/recipy/recipy/issues/83)
- Config file reading doesn't work as specified in the docs [\#64](https://github.com/recipy/recipy/issues/64)
- recipyCommon is not a package [\#62](https://github.com/recipy/recipy/issues/62)

**Closed issues:**

- GUI port number consistency [\#80](https://github.com/recipy/recipy/issues/80)
- README.md example is inconsistent [\#78](https://github.com/recipy/recipy/issues/78)
- Make UTC explicit in interfaces [\#75](https://github.com/recipy/recipy/issues/75)
- Tidy up inputs/outputs lists when there are none of them [\#74](https://github.com/recipy/recipy/issues/74)
- Add release notes [\#67](https://github.com/recipy/recipy/issues/67)
- Make proper documentation about config file options [\#60](https://github.com/recipy/recipy/issues/60)
- GUI Internal Server Error when viewing run with no git metadata [\#42](https://github.com/recipy/recipy/issues/42)
- keeping track of in- and outfile versions [\#41](https://github.com/recipy/recipy/issues/41)
- Logging of parameters [\#40](https://github.com/recipy/recipy/issues/40)
- Convert to use TinyDB [\#39](https://github.com/recipy/recipy/issues/39)
- Add 'command' functionality to recipy-cmd [\#37](https://github.com/recipy/recipy/issues/37)
- Create text index on all fields in database [\#35](https://github.com/recipy/recipy/issues/35)
- Upload to PyPI [\#23](https://github.com/recipy/recipy/issues/23)
- Add recipy-cmd function to create DB and set text index [\#10](https://github.com/recipy/recipy/issues/10)

**Merged pull requests:**

- GUI port number [\#81](https://github.com/recipy/recipy/pull/81) ([sjdenny](https://github.com/sjdenny))
- Add python highlighting in the README.md [\#79](https://github.com/recipy/recipy/pull/79) ([musically-ut](https://github.com/musically-ut))
- fixed typo in modulename of PatchPillow [\#73](https://github.com/recipy/recipy/pull/73) ([MichielCottaar](https://github.com/MichielCottaar))
- Command args [\#71](https://github.com/recipy/recipy/pull/71) ([oew1v07](https://github.com/oew1v07))
- Allow running recipy as a module \(python -m recipy\) [\#66](https://github.com/recipy/recipy/pull/66) ([kynan](https://github.com/kynan))
- Nibabel support [\#51](https://github.com/recipy/recipy/pull/51) ([MichielCottaar](https://github.com/MichielCottaar))

## [v0.1.0](https://github.com/recipy/recipy/tree/v0.1.0) (2015-08-16)

This is the first public release of recipy. The changes listed below are compared to early pre-release versions.

**Fixed bugs:**

- Not all pandas output methods are wrapped [\#11](https://github.com/recipy/recipy/issues/11)

**Closed issues:**

- Make compatible with Python 3 \(and still with Python 2\) [\#34](https://github.com/recipy/recipy/issues/34)
- Add configuration file support [\#33](https://github.com/recipy/recipy/issues/33)
- Add working setup.py [\#31](https://github.com/recipy/recipy/issues/31)
- Deal with git status of multiple files [\#29](https://github.com/recipy/recipy/issues/29)
- Use PyMongo v3? [\#28](https://github.com/recipy/recipy/issues/28)
- Add fuzzy searching option to recipy-cmd [\#27](https://github.com/recipy/recipy/issues/27)
- Make recipy-cmd work to search for full paths when given file in current directory [\#26](https://github.com/recipy/recipy/issues/26)
- Store diffs between script as executed and latest git commit [\#24](https://github.com/recipy/recipy/issues/24)
- Add hooks for scikit-learn [\#16](https://github.com/recipy/recipy/issues/16)
- Add hooks for scikit-image [\#15](https://github.com/recipy/recipy/issues/15)
- Add hooks for scikit-image [\#14](https://github.com/recipy/recipy/issues/14)
- Add hooks for PIL and Pillow [\#13](https://github.com/recipy/recipy/issues/13)
- Get more details on the environment and store in the DB [\#9](https://github.com/recipy/recipy/issues/9)
- Add a way to share entries from the DB [\#8](https://github.com/recipy/recipy/issues/8)
- Add interfaces from more languages and other tools [\#7](https://github.com/recipy/recipy/issues/7)
- Keep track of git commit version of script [\#6](https://github.com/recipy/recipy/issues/6)
- Add automatic installation for various OS's [\#5](https://github.com/recipy/recipy/issues/5)
- Add proper setup.py file [\#4](https://github.com/recipy/recipy/issues/4)
- Add hooks for more commonly-used Python modules [\#3](https://github.com/recipy/recipy/issues/3)
- Store full paths to files [\#2](https://github.com/recipy/recipy/issues/2)
- Make recipy-cmd search with fuzzy matches [\#1](https://github.com/recipy/recipy/issues/1)
# recipy

## What is it and who cares?
Imagine the situation: You’ve written some wonderful Python code which produces a beautiful graph as an output. You save that graph, naturally enough, as `graph.png`. You run the code a couple of times, each time making minor modifications. You come back to it the next week/month/year. Do you know how you created that graph? What input data? What version of your code? If you’re anything like me then the answer will often, frustratingly, be “no”. Of course, you then waste lots of time trying to work out how you created it, or even give up and never use it in that journal paper that will win you a Nobel Prize…

This talk will introduce ReciPy (from *recipe* and *python*), a Python module that will save you from this situation! (Although it can’t guarantee that your resulting paper will win a Nobel Prize!) With the addition of a single line of code to the top of your Python files, ReciPy will log each run of your code to a database, keeping track of the input files, output files and the version of your code, and then let you query this database to find out how you actually did create `graph.png`.

## Installation:
The easiest way to install is by simply running

    pip install recipy

Alternatively, you can clone this repository and run:

	python setup.py install

If you want to install the dependencies manually (they should be installed automatically if you're following the instructions above) then run:

	pip install -r requirements.txt

You can upgrade from a previous release by running:

	pip install -U recipy

To find out what has changed since the last release, see the [changelog](https://github.com/recipy/recipy/blob/master/CHANGELOG.md)

**Note:** Previous (unreleased) versions of recipy required MongoDB to be installed and set up manually. This is no longer required, as a pure Python database (TinyDB) is used instead. Also, the GUI is now integrated fully into recipy and does not require installing separately.

## Usage
Simply add the following line to the top of your Python script:

``` python
import recipy
```

Note that this **must** be the **very top** line of your script, before you import anything else.

Then just run your script as usual, and all of the data will be logged into the TinyDB database (don't worry, the database is automatically created if needed). You can then use the `recipy` script to quickly query the database to find out what run of your code produced what output file. So, for example, if you run some code like this:

``` python
import recipy
import numpy

arr = numpy.arange(10)
arr = arr + 500

numpy.save('test.npy', arr)
```

(Note the addition of `import recipy` at the beginning of script - but there are no other changes from a standard script)

Alternatively, run an unmodified script with `python -m recipy SCRIPT [ARGS ...]` to enable recipy logging. This invokes recipy's module entry point, which takes care of import recipy for you, before running your script.

it will produce an output called `test.npy`. To find out the details of the run which created this file you can search using

    recipy search test.npy

and it will display information like the following:

    Created by robin on 2015-05-25 19:00:15.631000
	Ran /Users/robin/code/recipy/example_script.py using /usr/local/opt/python/bin/python2.7
	Git: commit 91a245e5ea82f33ae58380629b6586883cca3ac4, in repo /Users/robin/code/recipy, with origin git@github.com:recipy/recipy.git
	Environment: Darwin-14.3.0-x86_64-i386-64bit, python 2.7.9 (default, Feb 10 2015, 03:28:08)
	Inputs:

	Outputs:
	  /Users/robin/code/recipy/test.npy

An alternative way to view this is to use the GUI. Just run `recipy gui` and a browser window will open with an interface that you can use to search all of your recipy 'runs':

![Screenshot of GUI](http://rtwilson.com/images/RecipyGUI.png)

If you want to log inputs and outputs of files read or written with built-in open, you need to do a little more work. Either use `recipy.open` (only requires `import recipy` at the top of your script), or add `from recipy import open` and just use `open`.
This workaround is required, because many libraries use built-in open internally, and you only want to record the files you explicitly opened yourself.

If you use Python 2, you can pass an `encoding` parameter to `recipy.open`. In this case `codecs` is used to open the file with proper encoding.

Once you've got some runs in your database, you can 'annotate' these runs with any notes that you want to keep about them. This can be particularly useful for recording which runs worked well, or particular problems you ran into. This can be done from the 'details' page in the GUI, or by running

	recipy annotate

which will open an editor to allow you to write notes that will be attached to the run. These will then be viewable via the command-line and the GUI when searching for runs.

There are other features in the command-line interface too: `recipy --help` to see the other options. You can view diffs, see all runs that created a file with a given name, search based on ids, show the latest entry and more:

	recipy - a frictionless provenance tool for Python

	Usage:
	  recipy search [options] <outputfile>
	  recipy latest [options]
	  recipy gui [options]
	  recipy annotate [<idvalue>]
	  recipy (-h | --help)
	  recipy --version

	Options:
	  -h --help     Show this screen
	  --version     Show version
	  -a --all      Show all results (otherwise just latest result given)
	  -f --fuzzy    Use fuzzy searching on filename
	  -r --regex    Use regex searching on filename
	  -i --id       Search based on (a fragment of) the run ID
	  -v --verbose  Be verbose
	  -d --diff     Show diff
	  -j --json     Show output as JSON
	  --no-browser  Do not open browser window
	  --debug       Turn on debugging mode

## Configuration
Recipy stores all of its configuration and the database itself in `~/.recipy`. Recipy's  main configuration file is inside this folder, called `recipyrc`. The configuration file format is very simple, and is based on Windows INI files - and having a configuration file is completely optional: the defaults will work fine with no configuration file.

An example configuration is:

	[ignored metadata]
	diff

	[general]
	debug

This simply instructs recipy not to save `git diff` information when it records metadata about a run, and also to print debug messages (which can be really useful if you're trying to work out why certain functions aren't patched). At the moment, the only possible options are:

 * `[general]`
	 * `debug` - print debug messages
 	 * `editor = vi` - Configure the default text editor that will be used when recipy needs you to type in a message. Use notepad if on Windows, for example
	 * `quiet` - don't print any messages
	 * `port` - specify port to use for the GUI
 *  `[data]`
	 * `file_diff_outputs` - store diff between the old output and new output file, if the output file exists before the script is executed
 *  `[database]`
 	 * `path = /path/to/file.json` - set the path to the database file
 * `[ignored metadata]`
	 * `diff` - don't store the output of `git diff` in the metadata for a recipy run
	 * `git` - don't store anything relating to git (origin, commit, repo etc) in the metadata for a recipy run
     * `input_hashes` - don't compute and store SHA-1 hashes of input files
     * `output_hashes` - don't compute and store SHA-1 hashes of output files
 * `[ignored inputs]`
 	 * List any module here (eg. `numpy`) to instruct recipy *not* to record inputs from this module, or `all` to ignore inputs from all modules
 * `[ignored outputs]`
 	 * List any module here (eg. `numpy`) to instruct recipy *not* to record outputs from this module, or `all` to ignore outputs from all modules	 

By default all metadata is stored (ie. no metadata is ignored) and debug messages are not shown. A `.recipyrc` file in the current directory takes precedence over the `~/.recipy/recipyrc` file, allowing per-project configurations to be easily handled.

**Note:** No default configuration file is provided with recipy, so if you wish to configure anything you will need to create a properly-formatted file yourself.

## How it works
When you import recipy it adds a number of classes to `sys.meta_path`. These are then used by Python as part of the importing procedure for modules. The classes that we add are classes derived from `PatchImporter`, often using the easier interface provided by `PatchSimple`, which allow us to wrap functions that do input/output in a function that calls recipy first to log the information.

Generally, most of the complexity is hidden away in `PatchImporter` and `PatchSimple` (plus `utils.py`), so the actual code to wrap a module, such as `numpy` is fairly simple:

``` python
# Inherit from PatchSimple
class PatchNumpy(PatchSimple):
    # Specify the full name of the module
    modulename = 'numpy'

    # List functions that are involved in input/output
    # these can be anything that can go after "modulename."
    # so they could be something like "pyplot.savefig" for example
    input_functions = ['genfromtxt', 'loadtxt', 'load', 'fromfile']
    output_functions = ['save', 'savez', 'savez_compressed', 'savetxt']

    # Define the functions that will be used to wrap the input/output
    # functions.
    # In this case we are calling the log_input function to log it to the DB
    # and we are giving it the 0th argument from the function (because all of
    # the functions above take the filename as the 0th argument), and telling
    # it that it came from numpy.
    input_wrapper = create_wrapper(log_input, 0, 'numpy')
    output_wrapper = create_wrapper(log_output, 0, 'numpy')
```

A class like this must be implemented for each module whose input/output needs logging. At the moment the following input and output functions are patched:

Patched modules
===============

This table lists the modules recipy has patches for, and the input and output functions that are patched.

<table>
<colgroup>
<col width="33%" />
<col width="33%" />
<col width="33%" />
</colgroup>
<thead>
<tr class="header">
<th align="left">Module</th>
<th align="left">Input functions</th>
<th align="left">Output functions</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left"><code>pandas</code></td>
<td align="left"><code>read_csv</code>, <code>read_table</code>, <code>read_excel</code>, <code>read_hdf</code>, <code>read_pickle</code>, <code>read_stata</code>, <code>read_msgpack</code></td>
<td align="left"><code>DataFrame.to_csv</code>, <code>DataFrame.to_excel</code>, <code>DataFrame.to_hdf</code>, <code>DataFrame.to_msgpack</code>, <code>DataFrame.to_stata</code>, <code>DataFrame.to_pickle</code>, <code>Panel.to_excel</code>, <code>Panel.to_hdf</code>, <code>Panel.to_msgpack</code>, <code>Panel.to_pickle</code>, <code>Series.to_csv</code>, <code>Series.to_hdf</code>, <code>Series.to_msgpack</code>, <code>Series.to_pickle</code></td>
</tr>
<tr class="even">
<td align="left"><code>matplotlib.pyplot</code></td>
<td align="left"></td>
<td align="left"><code>savefig</code></td>
</tr>
<tr class="odd">
<td align="left"><code>numpy</code></td>
<td align="left"><code>genfromtxt</code>, <code>loadtxt</code>, <code>fromfile</code></td>
<td align="left"><code>save</code>, <code>savez</code>, <code>savez_compressed</code>, <code>savetxt</code></td>
</tr>
<tr class="even">
<td align="left"><code>lxml.etree</code></td>
<td align="left"><code>parse</code>, <code>iterparse</code></td>
<td align="left"></td>
</tr>
<tr class="odd">
<td align="left"><code>bs4</code></td>
<td align="left"><code>BeautifulSoup</code></td>
<td align="left"></td>
</tr>
<tr class="even">
<td align="left"><code>gdal</code></td>
<td align="left"><code>Open</code></td>
<td align="left"><code>Driver.Create</code>, <code>Driver.CreateCopy</code></td>
</tr>
<tr class="odd">
<td align="left"><code>sklearn</code></td>
<td align="left"><code>datasets.load_svmlight_file</code></td>
<td align="left"><code>datasets.dump_svmlight_file</code></td>
</tr>
<tr class="even">
<td align="left"><code>nibabel</code></td>
<td align="left"><code>nifti1.Nifti1Image.from_filename</code>, <code>nifti2.Nifti2Image.from_filename</code>, <code>freesurfer.mghformat.MGHImage.from_filename</code>, <code>spm99analyze.Spm99AnalyzeImage.from_filename</code>, <code>minc1.Minc1Image.from_filename</code>, <code>minc2.Minc2Image.from_filename</code>, <code>analyze.AnalyzeImage.from_filename</code>, <code>parrec.PARRECImage.from_filename</code>, <code>spm2analyze.Spm2AnalyzeImage.from_filename</code></td>
<td align="left"><code>nifti1.Nifti1Image.to_filename</code>, <code>nifti2.Nifti2Image.to_filename</code>, <code>freesurfer.mghformat.MGHImage.to_filename</code>, <code>spm99analyze.Spm99AnalyzeImage.to_filename</code>, <code>minc1.Minc1Image.to_filename</code>, <code>minc2.Minc2Image.to_filename</code>, <code>analyze.AnalyzeImage.to_filename</code>, <code>parrec.PARRECImage.to_filename</code>, <code>spm2analyze.Spm2AnalyzeImage.to_filename</code></td>
</tr>
</tbody>
</table>

However, the code example above shows how easy it is to write a class to wrap a new module - so please feel free to submit a Pull Request to make recipy work with your favourite scientific modules!

## Test framework

recipy's test framework is in `integration_test`. The test framework has been designed to run under both Python 2.7+ and Python 3+. For more information see [recipy test framework](./docs/TestFramework.md).

The test framework is run on the following platforms:

* Travis CI: [![Integration test status image](https://travis-ci.org/recipy/recipy.svg)](https://travis-ci.org/recipy/recipy)
* AppVeyor: [![Build status](https://ci.appveyor.com/api/projects/status/irvathkx02yigjfn?svg=true
)](https://ci.appveyor.com/project/recipy/recipy)

##########################
User Manual
##########################

With the addition of a single line of code to the top of a Python script, recipy
logs each run of your code to a database, keeping track of the input files, output
files and the version of your code. It then lets you query this database to
help you to recall the exact steps you took to create a certain output file.

Logging Provenance Information
==============================

To log provenance information, simply add the following line to the top of your
code:

.. code-block:: python

   import recipy

Note that this **must** be the **very top** line of your script, before you
import anything else.

Then just run your script as usual, and all of the data will be logged into the
TinyDB database (don't worry, the database is automatically created if needed).
You can then use the ``recipy`` command to quickly query the database to find out
what run of your code produced what output file. So, for example, if you run some
code like this:

.. code-block:: python

   import recipy
   import numpy

   arr = numpy.arange(10)
   arr = arr + 500

   numpy.save('test.npy', arr)

(Note the addition of ``import recipy`` at the beginning of script - but there
are no other changes from a standard script.)

Alternatively, run an unmodified script with ``python -m recipy SCRIPT [ARGS ...]``
to enable recipy logging. This invokes recipy's module entry point, which takes
care of import recipy for you, before running your script.

Retrieving Information about Runs
=================================

it will produce an output called ``test.npy``. To find out the details of the
run which created this file you can search using

.. code-block:: sh

   recipy search test.npy


and it will display information like the following:

.. code-block:: sh

   Created by robin on 2015-05-25 19:00:15.631000
   Ran /Users/robin/code/recipy/example_script.py using /usr/local/opt/python/bin/python2.7
   Git: commit 91a245e5ea82f33ae58380629b6586883cca3ac4, in repo /Users/robin/code/recipy, with origin git@github.com:recipy/recipy.git
   Environment: Darwin-14.3.0-x86_64-i386-64bit, python 2.7.9 (default, Feb 10 2015, 03:28:08)
   Inputs:

   Outputs:
     /Users/robin/code/recipy/test.npy

An alternative way to view this is to use the GUI. Just run ``recipy gui`` and
a browser window will open with an interface that you can use to search all of
your recipy 'runs':

.. image:: http://rtwilson.com/images/RecipyGUI.png
   :target: http://rtwilson.com/images/RecipyGUI.png
   :alt: Screenshot of GUI

Logging Files Using Built-In Open
=================================

If you want to log inputs and outputs of files read or written with built-in
open, you need to do a little more work. Either use ``recipy.open``
(only requires ``import recipy`` at the top of your script), or add
``from recipy import open`` and just use ``open``.

This workaround is required, because many libraries use built-in open internally,
and you only want to record the files you explicitly opened yourself.

If you use Python 2, you can pass an ``encoding`` parameter to ``recipy.open``.
In this case :mod:`codecs` is used to open the file with proper encoding.

Annotating Runs
===============

Once you've got some runs in your database, you can 'annotate' these runs with
any notes that you want to keep about them. This can be particularly useful for
recording which runs worked well, or particular problems you ran into. This can
be done from the 'details' page in the GUI, or by running

.. code-block:: sh

   recipy annotate [run-id]


which will open an editor to allow you to write notes that will be attached to
the run. These will then be viewable via the command-line and the GUI when
searching for runs.

Saving Custom Values
====================

In your script, you can also add custom key-value pairs to the run:

.. code-block:: python

   recipy.log_values(key='value')
   recipy.log_values({'key': 'value'})


Please note that, at the moment, `these values are not displayed in the CLI or
in the GUI <https://github.com/recipy/recipy/issues/202>`_.

Command Line Interface
======================

There are other features in the command-line interface too: ``recipy --help``
to see the other options. You can view diffs, see all runs that created a file
with a given name, search based on ids, show the latest entry and more:

.. code-block:: sh

   recipy - a frictionless provenance tool for Python

   Usage:
     recipy search [options] <outputfile>
     recipy latest [options]
     recipy gui [options]
     recipy annotate [<idvalue>]
     recipy pm [--format <rst|plain>]
     recipy (-h | --help)
     recipy --version

   Options:
     -h --help        Show this screen
     --version        Show version
     -p --filepath    Search based on filepath rather than hash
     -f --fuzzy       Use fuzzy searching on filename
     -r --regex       Use regex searching on filename
     -i --id          Search based on (a fragment of) the run ID
     -a --all         Show all results (otherwise just latest result given)
     -v --verbose     Be verbose
     -d --diff        Show diff
     -j --json        Show output as JSON
     --no-browser     Do not open browser window
     --debug          Turn on debugging mode

Configuration
=============

By default, recipy stores all of its configuration and the database itself in
``~/.recipy``. Recipy's  main configuration file is inside this folder, called
``recipyrc``. The configuration file format is very simple, and is based on
Windows INI files - and having a configuration file is completely optional:
the defaults will work fine with no configuration file.

An example configuration is:

.. code-block:: sh

   [ignored metadata]
   diff

   [general]
   debug


This simply instructs recipy not to save ``git diff`` information when it
records metadata about a run, and also to print debug messages (which can be
really useful if you're trying to work out why certain functions aren't
patched). At the moment, the only possible options are:

* ``[general]``

  * ``debug`` - print debug messages
  * ``editor = vi`` - Configure the default text editor that will be used when
    recipy needs you to type in a message. Use notepad if on Windows, for example
  * ``quiet`` - don't print any messages
  * ``port`` - specify port to use for the GUI

* ``[data]``

  * ``file_diff_outputs`` - store diff between the old output and new output
    file, if the output file exists before the script is executed

* ``[database]``

  * ``path = /path/to/file.json`` - set the path to the database file

* ``[ignored metadata]``

  * ``diff`` - don't store the output of ``git diff`` in the metadata for a
    recipy run
  * ``git`` - don't store anything relating to git (origin, commit, repo etc)
    in the metadata for a recipy run
  * ``input_hashes`` - don't compute and store SHA-1 hashes of input files
  * ``output_hashes`` - don't compute and store SHA-1 hashes of output files

* ``[ignored inputs]``

  * List any module here (eg. ``numpy``\ ) to instruct recipy *not* to record
    inputs from this module, or ``all`` to ignore inputs from all modules

* ``[ignored outputs]``

  * List any module here (eg. ``numpy``\ ) to instruct recipy *not* to record
    outputs from this module, or ``all`` to ignore outputs from all modules

By default all metadata is stored (ie. no metadata is ignored) and debug messages
are not shown. A ``.recipyrc`` file in the current directory takes precedence over
the ``~/.recipy/recipyrc`` file, allowing per-project configurations to be easily
handled.

**Note:** No default configuration file is provided with recipy, so if you wish
to configure anything you will need to create a properly-formatted file yourself.
****************
Test Framework
****************

recipy's test framework is in ``integration_test``. The test framework
has been designed to run under both Python 2.7+ and Python 3+.

Running Tests with py.test
==========================

The tests use the `py.test <http://pytest.org>`_ test framework, and are
run using its ``py.test`` command. Useful ``py.test`` flags and
command-line options include:


* ``-v``\ : increase verbosity. This shows the name each test function
  that is run.
* ``-s``\ : show any output to standard output.
* ``-rs``\ : show extra test summary information for tests that were
  skipped.
* ``--junit-xml=report.xml``\ : create a JUnit-style test report in the
  file ``report.xml``.

Running General Tests
=====================

To run tests of recipy's command-line functions, run:

.. code-block:: console

   py.test -v integration_test/test_recipy.py

To run tests of recipy's ``.recipyrc`` configuration file, run:

.. code-block:: console

   py.test -v integration_test/test_recipyrc.py

To run tests of recipy invocations using ``python -m recipy script.py``\ ,
run:

.. code-block:: console

   py.test -v integration_test/test_m_flag.py

Running Package-Specific Tests
==============================

To run tests that check recipy logs information about packages it has
been configured to log, run:

.. code-block:: console

   py.test -v -rs integration_test/test_packages.py

**Note:** this assumes that all the packages have been installed.

To run a single test, provide the name of the test, for example:

.. code-block:: console

   $ py.test -v -rs integration_test/test_packages.py::test_scripts\[run_numpy_py_loadtxt\]

**Note:** ``[`` and ``]`` need to be prefixed by ``\``.

Package-specific tests use a test configuration file located in
``integration_test/config/test_packages.yml``.

You can specify a different test configuration file using a
``RECIPY_TEST_CONFIG`` environment variable. For example:

.. code-block:: console

   RECIPY_TEST_CONFIG=test_my_package.yml py.test -v -rs \
       integration_test/test_packages.py

**Note:** the above command should be typed on one line, omitting ``\``.

For Windows, run:

.. code-block:: console

   set RECIPY_TEST_CONFIG=test_my_package.yml
   py.test -v -rs integration_test\test_packages.py

Test Configuration File
-----------------------

The test configuration file is written in `YAML <http://yaml.org/>`_
(YAML Ain't Markup Language). YAML syntax is:


* ``---`` indicates the start of a document.
* ``:`` denotes a dictionary. ``:`` must be followed by a space.
* ``-`` denotes a list.

The test configuration file has format:

.. code-block:: console

   ---
   script: SCRIPT
   standalone: True | False
   libraries: [LIBRARY, LIBRARY, ... ]
   test_cases:
   - libraries: [LIBRARY, LIBRARY, ... ]
     arguments: [..., ..., ...]
     inputs: [INPUT, INPUT, ...]
     outputs: [OUTPUT, OUTPUT, ...]
   - libraries: [LIBRARY, LIBRARY, ... ]
     arguments: [..., ..., ...]
     inputs: [INPUT, INPUT, ...]
     outputs: [OUTPUT, OUTPUT, ...]
     skip: "Known issue with recipy"
     skip_py_version: [3.4, ...]
   - etc
   ---
   script: SCRIPT
   etc

Each script to be tested is defined by:


* ``SCRIPT``\ : Python script, with a relative or absolute path. For
  recipy sample scripts (see below), the script is assumed to be in a
  sub-directory ``integration_test/packages``.
* ``standalone``\ : is the script a standalone script? If ``False``\ , or if
  omitted, then the script is assumed to be a recipy sample script
  (see below).
* ``libraries``\ : A list of zero or more Python libraries used by the
  script, which are expected to be logged by recipy when the script
  is run regardless of arguments (i.e. any libraries common to all
  test cases). If none, then this can be omitted.

Each script also has one or more test cases, each of which defines:


* ``libraries``\ : A list of zero or more Python libraries used by the
  script, which are expected to be logged by recipy when the script
  is run with the given arguments for this test case. If none, then
  this can be omitted.
* ``arguments``\ : A list of arguments to be passed to the script. If
  none, then this can be omitted.
* ``inputs``\ : A list of zero or more input files which the script will
  read, and which are expected to be logged by recipy when running
  the script with the arguments. If none, then this can be omitted.
* ``outputs``\ : A list of zero or more output files which the script
  will write, and which are expected to be logged by recipy when
  running the script with the arguments. If none, then this can be
  omitted.
* ``skip``\ : An optional value. If present this test case is marked as
  skipped. The value is the reason for skipping the test case.
* ``skip_py_version``\ : An optional value. If present this test case is marked
   as skipped if the current Python version is in the list of values. Should
   be used when a patched library does not support a Python version that is
   supported by recipy.

For example:

.. code-block:: console

   ---
   script: run_numpy.py
   libraries: [numpy]
   test_cases:
   - arguments: [loadtxt]
     inputs: [input.csv]
   - arguments: [savetxt]
     outputs: [output.csv]
   - arguments: [load_and_save_txt]
     inputs: [input.csv]
     outputs: [output.csv]
   ---
   script: "/home/users/user/run_my_script.py"
   standalone: True
   test_cases:
   - arguments: [ ]
     libraries: [ numpy ]
     outputs: [ data.csv ]
   ---
   script: run_nibabel.py
   libraries: [ nibabel ]
   test_cases:
   - arguments: [ analyze_from_filename ]
     inputs: [ analyze_image ]
   - arguments: [ analyze_to_filename ]
     outputs: [ out_analyze_image ]
   - arguments: [ minc1_from_filename ]
     inputs: [ minc1_image ]
   - arguments: [ minc1_to_filename ]
     outputs: [ out_minc1_image ]
     skip: "nibabel.minc1.Minc1Image.to_filename raises NotImplementedError"

There may be a number of entries for a single script, if desired. For
example:

.. code-block:: console

   ---
   script: run_numpy.py
   libraries: [numpy]
   test_cases:
   - arguments: [loadtxt]
     inputs: [input.csv]
   - arguments: [savetxt]
     outputs: [output.csv]
   ---
   script: run_numpy.py
   libraries: [numpy]
   test_cases:
   - arguments: [load_and_save_txt]
     inputs: [input.csv]
     outputs: [output.csv]

It is up to you to ensure the ``library``\ , ``input`` and ``output`` file
names record the libraries, input and output files used by the
associated script, and which you expect to be logged by recipy.

Comments can be added to configuration files, prefixed by ``#``\ , for
example:

.. code-block:: console

   # This is a comment

Issues
======

The sample scripts in ``integration_tests/packages`` may fail to run
with older versions of third-party packages. Known package versions
that can cause failures are listed in `Package versioning
problems <./PackageVersionFailures.md>`_.

Certain third-party packages gave rise to issues, when attempting to
configure the test framework for these. The packages and issues, and
how the test framework has been configured to currently skip these are
described in `recipy and third-party package issues <./Issues.md>`_.

How the Test Framework uses Test Configuration Files
====================================================

A test configuration file is used to auto-generate test functions for
each test case using py.test's support for
`parameterization <http://doc.pytest.org/en/latest/parametrize.html>`_.

In the first example above, 8 test functions are created, 3 for
``run_numpy.py`` and 1 for ``run_my_scripy.py`` and 4 for ``run_nibabel.py``
(of which 1 is marked to be skipped. In the second example, 3 test
functions are created, 2 for the first group of ``run_numpy.py`` test
cases and 1 for the second group.

Test function names are auto-generated according to the following
template:

.. code-block:: console

   test_scripts[SCRIPT_ARGUMENTS]

where ``SCRIPT`` is the ``script`` value and arguments the ``argument``
values, concatenated using underscores (\ ``_``\ ) and with all forward
slashes, backslashes, colons, semi-colons and spaces also replaced by
``_``. For example, ``test_scripts[run_nibabel_py_analyze_from_filename]``.

The test framework runs the script with its arguments as follows. For
recipy sample scripts:

.. code-block:: console

   python -m integration_test.package.SCRIPT ARGUMENTS

For scripts marked ``standalone: True``\ :

.. code-block:: console

   python SCRIPT ARGUMENTS

Once the script has run, the test framework carries out the following
checks on the recipy database:


* There is only one new run in the database i.e. number of logs has
  increased by 1.
* ``script`` refers to the same file as the ``script``.
* ``command_args`` matches the test case's ``arguments``.
* ``libraries`` matches all the test case's ``libraries`` and all the
  ``libraries`` common to all test case's for a script, and the
  recorded versions match the versions used when ``script`` was run.
* ``inputs`` match the test case's ``inputs`` (in terms of local file
  names).
* ``outputs`` match test case's ``outputs`` (in terms of local file
  names).
* ``date`` is a valid date.
* ``exit_date`` is a valid date and is <= ``date``.
* ``command`` holds the current Python interpreter.
* ``environment`` holds the operating system and version of the current
  Python interpreter.
* ``author`` holds the current user.
* ``description`` is empty.

Recipy Sample Scripts
---------------------

``integration_test/packages`` has a collection of package-specific
scripts. Each script corresponds to one package logged by recipy. Each
script has a function that invokes each of the input/output functions
of a specific package logged by recipy. For example, ``run_numpy.py``
has functions that invoke:


* ``numpy.loadtxt``
* ``numpy.savetxt``
* ``numpy.fromfile``
* ``numpy.save``
* ``numpy.savez``
* ``numpy.savez_compressed``
* ``numpy.genfromtxt``

Each function is expected to invoke input and/or output functions
using one or more functions which recipy can log.

Input and output files are the responsibility of each script
itself. It can either create its own input and output files, or
read these in from somewhere (but it is not expected that the caller
(i.e. the test framework) create these.

Each test class has access to its own directory, via a
``self.current_dir`` field. It can use this to access any files it
needs within the current directory or, by convention, within a
sub-directory of ``data`` (for example ``run_numpy.py`` assumes a
``data/numpy`` sub-directory).

These scripts consist of classes that inherit from
``integration_test.base.Base`` which provides sub-classes with a simple
command-line interface which takes a function name as argument and, if
that function is provided by the script's class (and takes no
arguments beyond ``self``\ ), invokes that function. For example:

.. code-block:: console

   python SCRIPT.py FUNCTION

A sample script can be run as follows:

.. code-block:: console

   python -m integration_test.packages.SCRIPT FUNCTION

For example:

.. code-block:: console

   python -m integration_test.packages.run_numpy loadtxt

``test_packages.py`` assumes that if it is given a relative path to a
script, then that script is in ``integration_test/packages`` and will
create this form of invocation.

**Running scripts as modules**

Note that the script needs to be specified as a module that is run as
a script (the ``-m`` flag). Running it directly as a script e.g.

.. code-block:: console

   $ python integration_test/packages/run_numpy.py loadtxt

will fail:

.. code-block:: console

   Traceback (most recent call last):
     File "integration_test/packages/run_numpy.py", line 17, in <module>
       from integration_test.packages.base import Base
   ImportError: No module named 'integration_test.packages'

For the technical detail of why this is so, please see `Execution of
Python code with -m option or
not <http://stackoverflow.com/questions/22241420/execution-of-python-code-with-m-option-or-not>`_.

Providing a Test Configuration for Any Script
=============================================

A recipy test configuration can be written for any script that uses
recipy. For example, to test a script that uses ``numpy.loadtxt`` you
could write a configuration file which specifies:


* Full path to your script.
* Command-line arguments to be passed to your script.
* Libraries you expect to be logged by recipy.
* Local names of input files you expect to be logged by recipy.
* Local names of output files you expect to be logged by recipy.

For example, ``my_tests.yml``\ :

.. code-block:: console

   ---
   script: "/home/ubuntu/sample/run_numpy.py"
   standalone: True
   test_cases:
   - arguments: [ "/home/ubuntu/data/data.csv",
                  "/home/ubuntu/data/out.csv" ]
     libraries: [ numpy ]
     inputs: [ data.csv ]
     outputs: [ out.csv ]

You can run this as follows:

.. code-block:: console

   RECIPY_TEST_CONFIG=my_tests.yml py.test -v -rs \
       integration_test/test_packages.py

The output might look like

.. code-block:: console

   ============================= test session starts ==============================
   platform linux2 -- Python 2.7.12, pytest-2.9.2, py-1.4.31, pluggy-0.3.1 -- /home/ubuntu/anaconda2/bin/python
   cachedir: .cache
   rootdir: /home/ubuntu/recipy, inifile:
   collected 1 items

   integration_test/test_packages.py::test_scripts[run_numpy_py__home_ubuntu_data_data_csv__home_ubuntu_data_out_csv] PASSED

   =========================== 1 passed in 4.39 seconds ===========================

If using Anaconda and Git Bash on Windows, the file might look like:

.. code-block:: console

   ---
   script: "c:/Users/mjj/Local\ Documents/sample/run_numpy.py"
   standalone: True
   test_cases:
   - arguments: [ "c:/Users/mjj/Local\ Documents/data/data.csv",
                  "c:/Users/mjj/Local\ Documents/data/out.csv" ]
     libraries: [ numpy ]
     inputs: [ data.csv ]
     outputs: [ out.csv ]

Note the escaped spaces in the path.

Test Framework Limitations
==========================

The test framework does not support filtering tests depending upon
which versions of packages are being tested e.g. specific versions of
numpy or matplotlib. The test framework is designed to run tests
within a single execution environment: a Python interpreter and a set
of installed libraries.

If wishing to test different versions of packages then this could be
done by:


* Writing a Python script that invokes input/output functions of that
  package.
* Writing a test configuration file that just runs that script.
* Setting up a test environment (e.g. as part of a Travis CI or
  AppVeyor configuration file) that installs the specific package and
  runs ``py.test integration_test/test_packages.py`` using the
  test configuration file.

``test_recipy.py`` does not validate whether multiple test results are
returned by ``recipy search -i``.

Recipy Dependencies
===================

The test framework has no dependencies on any other part of the recipy
repository: it uses recipy as if it were a stand-alone tool and
queries the recipy database directly.
.. _database_schema:

Database Schema
===============

.. list-table::
   :header-rows: 1

   * - Field
     - Description
     - Example
   * - ``unique_id``
     -
     - ff129fee-c0d9-47cc-b0a9-997a530230a8
   * - ``author``
     - Username of the one running the script
     -
   * - ``description``
     - Currently not used?
     -
   * - ``inputs``
     - List of input files (can be empty)
     -
   * - ``outputs``
     - List of output files (can be empty)
     -
   * - ``script``
     - Path to the script that was run
     -
   * - ``command``
     - Path to the python executable that was used to run the script
     - ``/usr/bin/python2.7``
   * - ``environment``
     - Information about the environment in which the script was run (including Python version)
     - Linux-4.15.0-29-generic-x86_64-with-debian-stretch-sid python 2.7.15 -Anaconda, Inc.- (default, May 1 2018, 23:32:55)
   * - ``date``
     - Date and time the script was run (in UTC)
     -
   * - ``command_args``
     - Command line arguments of the script (can be empty)
     -
   * - ``warnings``
     - Warnings issued during execution of the script (if any)
     -
   * - ``exception``
     - Exception raised during execution of the script (if any)
     -
   * - ``libraries``
     - Names and versions of patched libraries used in this script
     -
   * - ``notes``
     - Notes added by the user by running ``recipy annotate`` or in the gui
     -
   * - ``custom_values``
     - Values logged explicitly in the script by calling ``recipy.log_values({'name': 'value'})``
     -
   * - ``gitrepo``\ /\ ``svnrepo``
     -
     -
   * - ``gitcommit``\ /\ ``svncommit``
     - Current git/svn commit
     -
   * - ``gitorigin``
     -
     -
   * - ``diff``
     -
     -
##########################
For Developers
##########################

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   how_does_it_work
   creating_patches
   databaseSchema
   TestFramework
#######################
Known Problems
#######################

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Issues
   PackageVersionFailures
###################
Installation
###################

The easiest way to install is by simply running:

.. code-block:: sh

    pip install recipy

Alternatively, you can clone this repository and run:

.. code-block:: sh

	python setup.py install

If you want to install the dependencies manually (they should be installed
automatically if you're following the instructions above) then run:

.. code-block:: sh

	pip install -r requirements.txt

You can upgrade from a previous release by running:

.. code-block:: sh

	pip install -U recipy

To find out what has changed since the last release, see the
`changelog <https://github.com/recipy/recipy/blob/master/CHANGELOG.md>`_

**Note:** Previous (unreleased) versions of recipy required MongoDB to be
installed and set up manually. This is no longer required, as a pure Python
database (TinyDB) is used instead. Also, the GUI is now integrated fully into
recipy and does not require installing separately.
How does it work?
=================

For each run of your script, recipy logs information about which script was run,
the command line arguments, Python version, git or svn repo, warnings, errors,
and inputs and outputs (see :ref:`database_schema` for a complete overview).
Gathering most of this information is straightforward using the Python Standard
Library and specialized packages (such as,
`GitPython <https://gitpython.readthedocs.io/en/stable/>`_
or `svn <https://github.com/dsoprea/PySvn>`_).
Automatically logging inputs and outputs, recipy's most important feature, is
more complicated.

When a Python module that reads or writes files is imported, recipy wraps
methods for reading and writing files to log file paths to the database.
To make this happen, recipy contains a patch for every library that reads or
writes files.  When you import recipy, the patches are added to
:data:`sys.meta_path` so they can be used to wrap a module's functions that read
or write files when it is imported. This is why `import recipy` should be called
before importing other modules.

Currently, recipy contains patches for many popular (scientific) libraries,
including :mod:`numpy` and :mod:`pandas`. For an overview see :ref:`Patched Modules`.
Is your favourite module missing? Have a look at :ref:`Creating Patches`.
We are looking forward to your pull request!
.. recipy documentation master file, created by
   sphinx-quickstart on Sun Oct 14 20:23:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the recipy documentation!
====================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   user_manual
   patched_modules
   developer
   known_problems



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
******************************
Package Versioning Problems
******************************

The sample scripts in ``integration_tests/packages`` may fail to run
with older versions of third-party packages. Known package versions
that can cause failures are listed below.

These can arise due to differences in versions of packages (especially
their APIs) installed by ``conda``\ , ``pip``\ , ``easy_install`` or within
Python-related packages installed via ``apt-get`` or ``yum`` under Linux.

NiBabel ``AttributeError``
===========================

.. code-block:: console

   $ python -m integration_test.packages.run_nibabel nifti2_from_filename
       data = nib.Nifti2Image.from_filename(file_name)
   AttributeError: 'module' object has no attribute 'Nifti2Image'

   $ python -m integration_test.packages.run_nibabel nifti2_to_filename
       img = nib.Nifti2Image(self.get_data(), self,get_affine())
   AttributeError: 'module' object has no attribute 'Nifti2Image'

   $ python -m integration_test.packages.run_nibabel minc1_from_filename
       data = nib.minc1.Minc1Image.from_filename(file_name)
   AttributeError: 'module' object has no attribute 'minc1'

   $ python -m integration_test.packages.run_nibabel minc1_to_filename
       img = nib.minc1.Minc1Image(self.get_data(), np.eye(4))
   AttributeError: 'module' object has no attribute 'minc1'

   $ python -m integration_test.packages.run_nibabel minc2_from_filename
       data = nib.minc2.Minc2Image.from_filename(file_name)
   AttributeError: 'module' object has no attribute 'minc2'

   $ python -m integration_test.packages.run_nibabel minc2_to_filename
       img = nib.minc2.Minc2Image(self.get_data(), np.eye(4))
   AttributeError: 'module' object has no attribute 'minc2'

   $ python -m integration_test.packages.run_nibabel parrec_from_filename
       data = nib.parrec.PARRECImage.from_filename(file_name)
   AttributeError: 'module' object has no attribute 'parrec'

   $ python -m integration_test.packages.run_nibabel parrec_to_filename
       img = nib.parrec.PARRECImage.from_filename(file_name)
   AttributeError: 'module' object has no attribute 'parrec'
   `

Fails on nibabel 1.2.2 (bundled in Ubuntu package ``python-nibabel``\ )
due to change in package API.

Succeeds on 2.0.2.

PIL ``AttributeError``
========================

.. code-block:: console

   $ python -m integration_test.packages.run_pil image_open
     File "/usr/lib/python2.7/dist-packages/PIL/Image.py", line 528, in __getattr__
       raise AttributeError(name)
   AttributeError: __exit__

   $ python -m integration_test.packages.run_pil image_save
   ...as above...

Fails on PIL/pillow 2.3.0 (bundled in Ubuntu package ``python-pillow``\ )
due to change in package API.

Succeeds on 3.2.0+.

pandas ``TypeError``
======================

.. code-block:: console

   $ python -m integration_test.packages.run_pandas read_excel
   TypeError: read_excel() takes exactly 2 arguments (1 given)

   $ python3 -m integration_test.packages.run_pandas read_excel
   TypeError: read_excel() missing 1 required positional argument: 'sheetname'

   $ python -m integration_test.packages.run_pandas read_hdf
   TypeError: read_hdf() takes exactly 2 arguments (1 given)

   $ python3 -m integration_test.packages.run_pandas read_hdf
   TypeError: read_hdf() missing 1 required positional argument: 'key'

Fails on pandas 0.13.1 (bundled in Ubuntu package ``python-pandas``\ ) due
to change in package API.

Succeeds on 0.18.1.

pandas ``ImportError``
========================

.. code-block:: console

   $ python -m integration_test.packages.run_pandas read_pickle
   ImportError: No module named indexes.base

   $ python3 -m integration_test.packages.run_pandas read_pickle
   ImportError: No module named 'pandas.indexes'

   During handling of the above exception, another exception occurred:
   ...

Fails on pandas 0.13.1 (bundled in Ubuntu package ``python-pandas``\ ) due
to change in package API.

Succeeds on 0.18.1.

pandas ``ValueError``
========================

.. code-block:: console

   $ python -m integration_test.packages.run_pandas read_msgpack
   ValueError: Unpack failed: error = -1

   $ python3 -m integration_test.packages.run_pandas read_msgpack
   ...as above...

Fails on pandas 0.13.1 (bundled in Ubuntu package ``python-pandas``\ ) due
to change in file format. The scripts work if run using data files
created by pandas 0.13.1.

Succeeds on 0.18.1.

skimage ``NameError``
=======================

.. code-block:: console

   $ python -m integration_test.packages.run_skimage io_load_sift
   NameError: name 'file' is not defined

   $ python -m integration_test.packages.run_skimage io_load_surf
   NameError: name 'file' is not defined

Fails on Python 3 as ``file()`` built-in removed in Python -M
Integration_Test.Packages.3 (see
`Builtins <https://docs.python.org/release/3.0/whatsnew/3.0.html#builtins>`_\ ).

skimage ``ImportError``
=========================

``run_skimage.py`` examples fail with:

.. code-block:: console

   Traceback (most recent call last):
     File "run_skimage.py", line 13, in <module>
       from skimage import external
   ImportError: cannot import name 'external'

Commenting out:

.. code-block:: console

   from skimage import external

allows non-external.tifffile examples to run.

Fails on skimage 0.9.3 (bundled in Ubuntu package ``python-skimage``\ )
due to change in package API.

Succeeds on 0.12.3.
################################
Patched Modules
################################

This table lists the modules recipy has patches for, and the input and output
functions that are patched.

Are you missing a patch for your favourite Python package? Learn how to
`create a new patch <#Creating Patches>`_. We are looking forward to your pull request!

=====================  =================================================  ===============================================
modulename             input_functions                                    output_functions
=====================  =================================================  ===============================================
``pandas``             ``read_csv``,                                      ``DataFrame.to_csv``,
                       ``read_table``,                                    ``DataFrame.to_excel``,
                       ``read_excel``,                                    ``DataFrame.to_hdf``,
                       ``read_hdf``,                                      ``DataFrame.to_msgpack``,
                       ``read_pickle``,                                   ``DataFrame.to_stata``,
                       ``read_stata``,                                    ``DataFrame.to_pickle``,
                       ``read_msgpack``                                   ``Panel.to_excel``,
                                                                          ``Panel.to_hdf``,
                                                                          ``Panel.to_msgpack``,
                                                                          ``Panel.to_pickle``,
                                                                          ``Series.to_csv``,
                                                                          ``Series.to_hdf``,
                                                                          ``Series.to_msgpack``,
                                                                          ``Series.to_pickle``
``matplotlib.pyplot``                                                     ``savefig``
``numpy``              ``genfromtxt``,                                    ``save``,
                       ``loadtxt``,                                       ``savez``,
                       ``fromfile``                                       ``savez_compressed``,
                                                                          ``savetxt``
``lxml.etree``         ``parse``,
                       ``iterparse``
``bs4``                ``BeautifulSoup``
``gdal``               ``Open``                                           ``Driver.Create``,
                                                                          ``Driver.CreateCopy``
``sklearn``            ``datasets.load_svmlight_file``                    ``datasets.dump_svmlight_file``
``nibabel``            ``nifti1.Nifti1Image.from_filename``,              ``nifti1.Nifti1Image.to_filename``,
                       ``nifti2.Nifti2Image.from_filename``,              ``nifti2.Nifti2Image.to_filename``,
                       ``freesurfer.mghformat.MGHImage.from_filename``,   ``freesurfer.mghformat.MGHImage.to_filename``,
                       ``spm99analyze.Spm99AnalyzeImage.from_filename``,  ``spm99analyze.Spm99AnalyzeImage.to_filename``,
                       ``minc1.Minc1Image.from_filename``,                ``minc1.Minc1Image.to_filename``,
                       ``minc2.Minc2Image.from_filename``,                ``minc2.Minc2Image.to_filename``,
                       ``analyze.AnalyzeImage.from_filename``,            ``analyze.AnalyzeImage.to_filename``,
                       ``parrec.PARRECImage.from_filename``,              ``parrec.PARRECImage.to_filename``,
                       ``spm2analyze.Spm2AnalyzeImage.from_filename``     ``spm2analyze.Spm2AnalyzeImage.to_filename``
``tifffile``           ``imread``                                         ``imsave``
``imageio``            ``core.functions.get_reader``,                     ``core.functions.get_writer``
                       ``core.functions.read``
``netCDF4``            ``Dataset``                                        ``Dataset``
``xarray``             ``open_dataset``,                                  ``Dataset.to_netcdf``,
                       ``open_mfdataset``,                                ``DataArray.to_netcdf``
                       ``open_rasterio``,
                       ``open_dataarray``
``iris``               ``iris.load``,                                     ``iris.save``
                       ``iris.load_cube``,
                       ``iris.load_cubes``,
                       ``iris.load_raw``
=====================  =================================================  ===============================================

To generate this table do:

.. code-block:: sh

   recipy pm --format rst
****************************************
Recipy and Third-Party Package Issues
****************************************

Issues encountered during test framework development and how the test
framework has been configured to get round these issues.

No information logged by recipy (recipy code is commented-out)
================================================================

All these tests are marked as skipped.

``PIL`` operations
-------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_pil image_open
   recipy latest

.. code-block:: console

   Run ID: cf695a6c-e9f0-4f4b-896f-9fa64127e2e8
   Created by ubuntu on 2016-10-26 11:59:04 UTC
   Ran /home/ubuntu/recipy/integration_test/packages/run_pil.py using /home/ubuntu/anaconda2/bin/python
   Using command-line arguments: image_open
   Git: commit 5dc58a3cf3432441c83fd9899768eb9b63583208, in repo /home/ubuntu/recipy, with origin https://mikej888@github.com/mikej888/recipy
   Environment: Linux-3.19.0-25-generic-x86_64-with-debian-jessie-sid, python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40)
   Libraries: recipy v0.3.0
   Inputs: none
   Outputs: none

``skimage`` operations
-----------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_skimage io_imread
   recipy latest

.. code-block:: console

   Run ID: 10c803f7-3741-4008-8548-2f8d7ba4462c
   Created by ubuntu on 2016-10-26 11:57:07 UTC
   Ran /home/ubuntu/recipy/integration_test/packages/run_skimage.py using /home/ubuntu/anaconda2/bin/python
   Using command-line arguments: io_imread
   Git: commit 5dc58a3cf3432441c83fd9899768eb9b63583208, in repo /home/ubuntu/recipy, with origin https://mikej888@github.com/mikej888/recipy
   Environment: Linux-3.19.0-25-generic-x86_64-with-debian-jessie-sid, python 2.7.12 |Anaconda custom (64-bit)| (default, Jul  2 2016, 17:42:40)
   Libraries: recipy v0.3.0
   Inputs: none
   Outputs: none

Inaccurate Information Logged by Recipy
=========================================

None.

Bugs Arising within Recipy Logging
====================================

All these tests are marked as skipped.

``recipy.open``
-----------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_python open

Under Python 3 this fails with:

.. code-block:: console

   recipy run inserted, with ID 9abe4113-5158-4dcb-8fe8-afd1b1f6505e

   Traceback (most recent call last):
     File "c:\Users\mjj\AppData\Local\Continuum\Anaconda3\lib\runpy.py", line 170,in _run_module_as_main
       "__main__", mod_spec)
     File "c:\Users\mjj\AppData\Local\Continuum\Anaconda3\lib\runpy.py", line 85, in _run_code
       exec(code, run_globals)
     File "c:\Users\mjj\Local Documents\recipy\recipy\integration_test\packages\run_python.py", line 45, in <module>
       PythonSample().invoke(sys.argv)
     File "c:\Users\mjj\Local Documents\recipy\recipy\integration_test\packages\base.py", line 57, in invoke
       function()
     File "c:\Users\mjj\Local Documents\recipy\recipy\integration_test\packages\run_python.py", line 39, in open
       with recipy.open('out.txt', 'w') as f:
     File "c:\Users\mjj\Local Documents\recipy\recipy\recipy\utils.py", line 20, in open
       mode = kwargs['mode']
   KeyError: 'mode'

Under Python 2 this fails with:

.. code-block:: console

   recipy run inserted, with ID 5d80b88b-0d56-428d-b9e0-d95eca423044

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_python.py", line 42, in <module>
       python_sample.invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_python.py", line 35, in open
       with recipy.open('out.txt', 'w') as f:
     File "recipy/utils.py", line 35, in open
       log_output(args[0], 'recipy.open')
     File "recipy/log.py", line 153, in log_output
       db.update(append("libraries", get_version(source), no_duplicates=True), eids=[RUN_ID])
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 377, in update
       cond, eids
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 230, in process_elements
       data = self._read()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 277, in _read
       return self._storage.read()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 31, in read
       raw_data = (self._storage.read() or {})[self._table_name]
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb_serialization/__init__.py", line 139, in read
       data = self.storage.read()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/storages.py", line 93, in read
       self._handle.seek(0, 2)
   ValueError: I/O operation on closed file

``bs4.beautifulsoup.prettify``
--------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_bs4 beautifulsoup

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_bs4.py", line 53, in <module>
       Bs4Sample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_bs4.py", line 49, in beautifulsoup
       print((soup.prettify()))
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/bs4/element.py", line 1160, in prettify
       return self.decode(True, formatter=formatter)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/bs4/__init__.py", line 439, in decode
       return prefix + super(BeautifulSoup, self).decode(
   TypeError: super() argument 1 must be type, not FunctionWrapper

``pandas.Panel.to_excel``
---------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_pandas panel_to_excel

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_pandas.py", line 355, in <module>
       PandasSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_pandas.py", line 195, in panel_to_excel
       panel.to_excel(file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/pandas/core/panel.py", line 460, in to_excel
       df.to_excel(writer, name, **kwargs)
     File "recipyCommon/utils.py", line 90, in f
       function(args[arg_loc], source)
     File "recipy/log.py", line 139, in log_output
       filename = os.path.abspath(filename)
     File "/home/ubuntu/anaconda2/lib/python2.7/posixpath.py", line 360, in abspath
       if not isabs(path):
     File "/home/ubuntu/anaconda2/lib/python2.7/posixpath.py", line 54, in isabs
       return s.startswith('/')
   AttributeError: '_XlwtWriter' object has no attribute 'startswith'

``nibabel.minc2.Minc2Image.from_filename``
------------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_nibabel minc2_from_filename

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 302, in <module>
       NibabelSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 143, in minc2_from_filename
       data = nib.minc2.Minc2Image.from_filename(file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 699, in from_filename
       return klass.from_file_map(file_map)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/minc1.py", line 299, in from_file_map
       minc_file = Minc1File(netcdf_file(fobj))
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/externals/netcdf.py", line 230, in __init__
       self._read()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/externals/netcdf.py", line 513, in _read
       self.filename)
   TypeError: Error: None is not a valid NetCDF 3 file

``nibabel.Nifti2Image.from_filename``
--------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_nibabel nifti2_from_filename

.. code-block:: console

   sizeof_hdr should be 348; set sizeof_hdr to 348
   data code 0 not supported; not attempting fix
   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 302, in <module>
       NibabelSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 182, in nifti2_from_filename
       data = nib.Nifti2Image.from_filename(file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/keywordonly.py", line 16, in wrapper
       return func(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/analyze.py", line 986, in from_filename
       return klass.from_file_map(file_map, mmap=mmap)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/keywordonly.py", line 16, in wrapper
       return func(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/analyze.py", line 947, in from_file_map
       header = klass.header_class.from_fileobj(hdrf)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/nifti1.py", line 594, in from_fileobj
       hdr = klass(raw_str, endianness, check)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/nifti1.py", line 577, in __init__
       check)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/analyze.py", line 252, in __init__
       super(AnalyzeHeader, self).__init__(binaryblock, endianness, check)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/wrapstruct.py", line 176, in __init__
       self.check_fix()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/wrapstruct.py", line 361, in check_fix
       report.log_raise(logger, error_level)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/batteryrunners.py", line 275, in log_raise
       raise self.error(self.problem_msg)
   nibabel.spatialimages.HeaderDataError: data code 0 not supported

``sklearn.load_svmlight_file`` and ``sklearn.dump_svmlight_file``
-------------------------------------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_sklearn load_svmlight_file

Under Python 3 this fails with:

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda3/lib/python3.5/runpy.py", line 184, in _run_module_as_main
       "__main__", mod_spec)
     File "/home/ubuntu/anaconda3/lib/python3.5/runpy.py", line 85, in _run_code
       exec(code, run_globals)
     File "/home/ubuntu/recipy/integration_test/packages/run_sklearn.py", line 16, in <module>
       from sklearn import datasets
     File "<frozen importlib._bootstrap>", line 969, in _find_and_load
     File "<frozen importlib._bootstrap>", line 958, in _find_and_load_unlocked
     File "<frozen importlib._bootstrap>", line 664, in _load_unlocked
     File "<frozen importlib._bootstrap>", line 634, in _load_backward_compatible
     File "/home/ubuntu/recipy/recipy/PatchImporter.py", line 52, in load_module
       mod = self.patch(mod)
     File "/home/ubuntu/recipy/recipy/PatchSimple.py", line 25, in patch
       patch_function(mod, f, self.input_wrapper)
     File "/home/ubuntu/recipy/recipyCommon/utils.py", line 82, in patch_function
       setattr(mod, old_f_name, recursive_getattr(mod, function))
     File "/home/ubuntu/recipy/recipyCommon/utils.py", line 54, in recursive_getattr
       prev_part = getattr(prev_part, part)
   AttributeError: module 'sklearn' has no attribute 'datasets'

Under Python 2 this fails with:

.. code-block:: console

   Traceback (most recent call last):
     File "recipy/log.py", line 165, in log_exception
       db.update({"exception": exception}, eids=[RUN_ID])
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 382, in update
       cond, eids
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 235, in process_elements
       func(data, eid)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 381, in <lambda>
       lambda data, eid: data[eid].update(fields),
   KeyError: 316

   Original exception was:
   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_sklearn.py", line 16, in <module>
       from sklearn import datasets
     File "recipy/PatchImporter.py", line 52, in load_module
       mod = self.patch(mod)
     File "recipy/PatchSimple.py", line 25, in patch
       patch_function(mod, f, self.input_wrapper)
     File "recipyCommon/utils.py", line 82, in patch_function
       setattr(mod, old_f_name, recursive_getattr(mod, function))
     File "recipyCommon/utils.py", line 54, in recursive_getattr
       prev_part = getattr(prev_part, part)
   AttributeError: 'module' object has no attribute 'datasets'
   Error in atexit._run_exitfuncs:
   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/atexit.py", line 24, in _run_exitfuncs
       func(*targs, **kargs)
     File "recipy/log.py", line 244, in hash_outputs
       for filename in run.get('outputs')]
   AttributeError: 'NoneType' object has no attribute 'get'
   Error in atexit._run_exitfuncs:
   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/atexit.py", line 24, in _run_exitfuncs
       func(*targs, **kargs)
     File "recipy/log.py", line 231, in log_exit
       db.update({'exit_date': exit_date}, eids=[RUN_ID])
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 382, in update
       cond, eids
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 235, in process_elements
       func(data, eid)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 381, in <lambda>
       lambda data, eid: data[eid].update(fields),
   KeyError: 316
   Error in sys.exitfunc:
   Error in sys.excepthook:
   Traceback (most recent call last):
     File "recipy/log.py", line 165, in log_exception
       db.update({"exception": exception}, eids=[RUN_ID])
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 382, in update
       cond, eids
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 235, in process_elements
       func(data, eid)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 381, in <lambda>
       lambda data, eid: data[eid].update(fields),
   KeyError: 316

   Original exception was:
   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/atexit.py", line 24, in _run_exitfuncs
       func(*targs, **kargs)
     File "recipy/log.py", line 231, in log_exit
       db.update({'exit_date': exit_date}, eids=[RUN_ID])
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 382, in update
       cond, eids
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 235, in process_elements
       func(data, eid)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/tinydb/database.py", line 381, in <lambda>
       lambda data, eid: data[eid].update(fields),
   KeyError: 316

Operations Not Implemented by Packages
========================================

All these tests are marked as skipped.

``nibabel.minc1.Minc1Image.to_filename``
-----------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_nibabel minc1_to_filename

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 302, in <module>
       NibabelSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 134, in minc1_to_filename
       img.to_filename(file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 781, in to_filename
       self.to_file_map()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 790, in to_file_map
       raise NotImplementedError
   NotImplementedError

``nibabel.minc2.Minc2Image.to_filename``
-----------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_nibabel minc2_to_filename

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 302, in <module>
       NibabelSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 154, in minc2_to_filename
       img.to_filename(file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 781, in to_filename
       self.to_file_map()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 790, in to_file_map
       raise NotImplementedError
   NotImplementedError

``nibabel.parrec.PARRECImage.to_filename``
-------------------------------------------

To replicate:

.. code-block:: console

   python -m integration_test.packages.run_nibabel parrec_to_filename

.. code-block:: console

   Traceback (most recent call last):
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 174, in _run_module_as_main
       "__main__", fname, loader, pkg_name)
     File "/home/ubuntu/anaconda2/lib/python2.7/runpy.py", line 72, in _run_code
       exec code in run_globals
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 302, in <module>
       NibabelSample().invoke(sys.argv)
     File "integration_test/packages/base.py", line 57, in invoke
       function()
     File "/home/ubuntu/recipy/integration_test/packages/run_nibabel.py", line 217, in parrec_to_filename
       img.to_filename(par_file_name)
     File "recipyCommon/utils.py", line 91, in f
       return wrapped(*args, **kwargs)
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 781, in to_filename
       self.to_file_map()
     File "/home/ubuntu/anaconda2/lib/python2.7/site-packages/nibabel/spatialimages.py", line 790, in to_file_map
       raise NotImplementedError
   NotImplementedError

Using py.test and recipy
------------------------

An issue that does not affect the test framework, but may affect
future test development is that recipy and py.test do not
integrate. For example, given test_sample.py:

.. code-block:: python

   class TestSample:

       def test_sample(self):
           pass

Running:

.. code-block:: console

   py.test test_sample.py

gives:

.. code-block:: console

   ============================= test session starts =============================
   platform win32 -- Python 3.5.1, pytest-3.0.2, py-1.4.31, pluggy-0.3.1
   rootdir: c:\Users\mjj\Local Documents, inifile:
   collected 1 items

   test_sample.py .

   ========================== 1 passed in 0.02 seconds ===========================

Adding:

.. code-block:: python

   import recipy

Running py.test gives:

.. code-block:: console

   ============================= test session starts =============================
   platform win32 -- Python 3.5.1, pytest-3.0.2, py-1.4.31, pluggy-0.3.1
   rootdir: c:\Users\mjj\Local Documents, inifile:
   collected 0 items / 1 errors

   =================================== ERRORS ====================================
   _______________________ ERROR collecting test_sample.py _______________________
   ..\appdata\local\continuum\anaconda3\lib\site-packages\_pytest\python.py:209: in fget
       return self._obj
   E   AttributeError: 'Module' object has no attribute '_obj'

   During handling of the above exception, another exception occurred:
   test_sample.py:1: in <module>
       import recipy
   ..\appdata\local\continuum\anaconda3\lib\site-packages\recipy-0.3.0-py3.5.egg\recipy\__init__.py:12: in <module>
       log_init()
   ..\appdata\local\continuum\anaconda3\lib\site-packages\recipy-0.3.0-py3.5.egg\recipy\log.py:74: in log_init
       add_git_info(run, scriptpath)
   ..\appdata\local\continuum\anaconda3\lib\site-packages\recipy-0.3.0-py3.5.egg\recipyCommon\version_control.py:30: in add_git_info
       repo = Repo(scriptpath, search_parent_directories=True)
   ..\appdata\local\continuum\anaconda3\lib\site-packages\gitpython-2.0.8-py3.5.egg\git\repo\base.py:139: in __init__
       raise NoSuchPathError(epath)
   E   git.exc.NoSuchPathError: c:\Users\mjj\AppData\Local\Continuum\Anaconda3\Scripts\py.test
   !!!!!!!!!!!!!!!!!!! Interrupted: 1 errors during collection !!!!!!!!!!!!!!!!!!!
   =========================== 1 error in 4.55 seconds ===========================
   Error in atexit._run_exitfuncs:
   Traceback (most recent call last):
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\recipy-0.3.0-py3.5.egg\recipy\log.py", line 242, in hash_outputs
       run = db.get(eid=RUN_ID)
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\tinydb-3.2.1-py3.5.egg\tinydb\database.py", line 432, in get
   TypeError: unhashable type: 'dict'
   Error in atexit._run_exitfuncs:
   Traceback (most recent call last):
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\recipy-0.3.0-py3.5.egg\recipy\log.py", line 231, in log_exit
       db.update({'exit_date': exit_date}, eids=[RUN_ID])
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\tinydb-3.2.1-py3.5.egg\tinydb\database.py", line 382, in update
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\tinydb-3.2.1-py3.5.egg\tinydb\database.py", line 235, in process_elements
     File "c:\users\mjj\appdata\local\continuum\anaconda3\lib\site-packages\tinydb-3.2.1-py3.5.egg\tinydb\database.py", line 381, in <lambda>
   TypeError: unhashable type: 'dict'
Creating Patches
================

Patches are derived from :class:`~recipy.PatchImporter.PatchImporter`, often
using the easier interface provided by :class:`~recipy.PatchSimple.PatchSimple`.
To create a patch, you need to specify the module name, the input and output
functions, and the argument referencing the file path.

Simple patches
***************

Because most of the complexity is hidden away, the actual code to wrap a module
is fairly simple. Using :class:`~recipy.PatchSimple.PatchSimple`, the patch
for `numpy <http://www.numpy.org>`_ looks like:

.. code-block:: python
  :linenos:

  from .PatchSimple import PatchSimple
  from .log import log_input, log_output, add_module_to_db
  from recipyCommon.utils import create_wrapper

  # Inherit from PatchSimple
  class PatchNumpy(PatchSimple):
    # Specify the full name of the module
    modulename = 'numpy'

    # List functions that are involved in input/output
    # these can be anything that can go after "modulename."
    # so they could be something like "pyplot.savefig" for example
    input_functions = ['genfromtxt', 'loadtxt', 'load', 'fromfile']
    output_functions = ['save', 'savez', 'savez_compressed', 'savetxt']

    # Define the functions that will be used to wrap the input/output
    # functions.
    # In this case we are calling the log_input function to log it to
    # the DB and we are giving it the 0th argument from the function
    # (because all of the functions above take the filename as the
    # 0th argument), and telling it that it came from numpy.
    input_wrapper = create_wrapper(log_input, 0, 'numpy')
    output_wrapper = create_wrapper(log_output, 0, 'numpy')

    # Add the module to the database, so we can list it.
    add_module_to_db(modulename, input_functions, output_functions)

First, we import the required functionality (lines 1-3). In addition to
:class:`~recipy.PatchSimple.PatchSimple`, we need functions that do the actual
logging (:meth:`~recipy.log.log_input` and :meth:`~recipy.log.log_output`),
write details about the module to the database
(:meth:`~recipy.log.add_module_to_db`), and a function to create wrappers
(:meth:`~recipyCommon.utils.create_wrapper`). Next, we define a class
for the patch that inherits from :class:`~recipy.PatchSimple.PatchSimple` (line
6). In this class, we need to define the module name (line 8), and the input and
output functions (line 13 and 14).
Because the patches rely on class attributes, it is important to use the
variable names specified in the example code.
Then, we need to define the input and output
wrapper (line 22 and 23). The wrappers are created using a predefined function,
that needs as inputs a logging function (i.e., :meth:`~recipy.log.log_input`
or :meth:`~recipy.log.log_output`), the index of the argument that specifies
the file path in the input or output function, and the name of the module.
Finally, the module is added to the database (line 26).

Enabling a patch
****************

Patch objects for specific modules can be found in
:mod:`~recipy.PatchBaseScientific` and :mod:`~recipy.PatchScientific`. To enable
a new patch, one more step is required; the constructor of the new object should
be called inside :meth:`~recipyCommon.utils.multiple_insert`. This function is
called at the bottom of :mod:`~recipy.PatchBaseScientific` and
:mod:`~recipy.PatchScientific`.

Patching more complex modules
*****************************

Most of the patched modules are based on :class:`~recipy.PatchSimple.PatchSimple`.
However, some modules require more complexity. Some modules, e.g.,
`netcdf4-python <http://unidata.github.io/netcdf4-python/>`_, use a file open like
method for reading and writing files; the method called is the same, and whether
the file is read or written depends on arguments:

.. code-block:: python
  :linenos:

  from netCDF4 import Dataset

  ds = Dataset('test.nc', 'r')  # read test.nc
  ds = Dataset('test.nc', 'w')  # write test.nc

For this situation a different patch type and wrapper creation function are
available: :class:`~recipy.PatchFileOpenLike` in combination with
:meth:`~recipy.log.create_argument_wrapper`.

The complete code example (:class:`~recipy.PatchNetCDF4`):

.. code-block:: python
  :linenos:

  from .PatchFileOpenLike import PatchFileOpenLike
  from .log import log_input, log_output, add_module_to_db
  from recipyCommon.utils import create_argument_wrapper

  class PatchNetCDF4(PatchFileOpenLike):
    modulename = 'netCDF4'

    functions = ['Dataset']

    wrapper = create_argument_wrapper(log_input, log_output, 0, 'mode', 'ra',
                                      'aw', 'r', 'netCDF4')

    add_module_to_db(modulename, functions, functions)

Another instance where it isn't possible to use
:class:`~recipy.PatchSimple.PatchSimple` is when not all input or output
functions have the file path in the same position. For example,
`xarray <http://xarray.pydata.org/>`_
has three functions for writing files, i.e., :meth:`~xarray.Dataset.to_netcdf`,
:meth:`~xarray.DataArray.to_netcdf`, and :func:`~xarray.save_mfdataset`.
For :meth:`~xarray.Dataset.to_netcdf` and :meth:`~xarray.DataArray.to_netcdf`,
the file paths are argument 0, as in the previous examples. However, the
argument for the file paths of :func:`~xarray.save_mfdataset` (it is a method
for writing multiple files at once) is 1.
With :class:`~recipy.PatchSimple.PatchSimple` there is no way to represent this.

A patch class that allows specifying separate wrappers for different functions
is :class:`~recipy.PatchMultipleWrappers.PatchMultipleWrappers`. Using this
class involves defining a :class:`~recipy.PatchMultipleWrappers.WrapperList` to
which inputs and outputs can be added. As can be seen in the following code
example (:class:`~recipy.PatchXarray`), you can specify a wrapper for a list
of functions (line 13 and 14) or for one function (line 15).

.. code-block:: python
  :linenos:

  from .PatchMultipleWrappers import PatchMultipleWrappers, WrapperList
  from .log import log_input, log_output, add_module_to_db

  class PatchXarray(PatchMultipleWrappers):
    modulename = 'xarray'

    wrappers = WrapperList()

    input_functions = ['open_dataset', 'open_mfdataset', 'open_rasterio',
                       'open_dataarray']
    output_functions = ['Dataset.to_netcdf', 'DataArray.to_netcdf']

    wrappers.add_inputs(input_functions, log_input, 0, modulename)
    wrappers.add_outputs(output_functions, log_output, 0, modulename)
    wrappers.add_outputs('save_mfdataset', log_output, 1, modulename)

    add_module_to_db(modulename, input_functions, output_functions)

Writing tests
*************

If you make a new patch, please include tests (your pull request won't be
accepted without them)! Recipy has a testing framework that checks
whether inputs and outputs are actually logged when a function is called
(see :ref:`Test Framework` for more details).

If you create a patch for a module, follow these steps to create tests:

1. Prepare small data files for testing, and create a directory
   ``integration_test/packages/data/<module name>`` containing these files.

   * Small means kilobytes!

2. Create a test script for the patch. The name of the script
   should be ``run_<module name>.py``.

   * It is probably easiest to copy one of the
     existing scripts in ``integration_test/packages/``, so you can reuse the
     set up (it is pretty self-explanatory).

3. Write a test method for each input/output method.

   * Be sure to add docstrings!

4. Add the test configuration.

   * Open ``integration_test/config/test_packages.yml``
   * Add a new section by typing:

   .. code-block:: sh

     ---
     script: run_<module name>.py
     libraries: [ modulename ]
     test_cases:

   * Add a test case for each method in ``run_<module name>.py``, specifying
     ``arguments``, a list containing the name of the test method,
     ``inputs``, a list of the input files this method needs (if any), and
     ``outputs``, a list of the output files this method  creates (if any):

   .. code-block:: sh

     - arguments: [ <test method name> ]
       inputs: [<input names>]
       outputs: [<output names>]

   * Specified input files should exist in
     ``integration_test/packages/data/<module name>``
   * Test cases can skipped by adding the line: ``skip: "<reason for skipping>"``
   * If you need to skip a test, please add a description of the problem plus
     how to reproduce it in the :ref:`Recipy and Third-Party Package Issues`
     section
   * If a test case should be skipped in certain Python versions, add
     ``skip_py_version: [ <python versions> ]``

5. Run the tests by typing: ``py.test -v integration_test/``
