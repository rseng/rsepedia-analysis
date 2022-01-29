# CHANGELOG

This is a manually generated log to track changes to the repository for each release.
Each section should include general headers such as **Implemented enhancements**
and **Merged pull requests**. Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes
 - migration guidance
 - changed behaviour

The versions coincide with releases on pip.

## [0.2.x](https://github.com/vsoch/gridtest/tree/master) (0.0.x)
 - code cleanup and JoSS review (0.0.15)
 - adding machine learning examples for grids (0.0.14)
 - refactoring to use Grid class, json-tricks export (0.0.13)
 - adding grid exports, and variables section (0.0.12)
 - adding "tests" level to config files (0.0.11)
 - adding grids section in test file, and function definition for args
 - first alpha release with basic grid test generation (0.0.1)
 - skeleton release (0.0.0)
# GridTest

[![PyPI version](https://badge.fury.io/py/gridtest.svg)](https://badge.fury.io/py/gridtest)
[![DOI](https://zenodo.org/badge/256346804.svg)](https://zenodo.org/badge/latestdoi/256346804)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02284/status.svg)](https://doi.org/10.21105/joss.02284)

Simple grid parameterization and testing setup for Python functions and modules.
See [Documentation](https://vsoch.github.io/gridtest/) to get started.

![docs/assets/img/logo/gridtest.gif](https://raw.githubusercontent.com/vsoch/gridtest/master/docs/assets/img/logo/gridtest.gif)

## Overview 

GridTest is a library that specializes in generating [parameter grids](https://vsoch.github.io/gridtest/#parameterization). The grids are most obviously used for testing, but can extend to other use cases.
In the context of testing, GridTest makes it easy to discover functions,
classes, and arguments for your python scripts or modules, and then generate
a template for you to easily populate. Outside of testing, you can define
grids that are version controlled, programatically defined with functions,
and easy to interact with from the command line or Python interpreter.
You might be interested in GridTest if you need:

   - low overhead tests for Python scripts and small packages
   - to generate input data for reproducible computations

To learn more, it's recommended to reference the [documentation](https://vsoch.github.io/gridtest/),
take a look at the [getting started](https://vsoch.github.io/gridtest/getting-started/index.html) pages,
or browse one of the many [tutorials](https://vsoch.github.io/gridtest/tutorials/index.html) available.

 * Free software: MPL 2.0 License

## Support

If you have any questions or requests for examples or tutorials, please don't hesitate
to [open an issue](https://github.com/vsoch/gridtest/issues).

## Contributing

Please see the [documentation contributing guide](https://vsoch.github.io/gridtest/contributing/index.html)
for details on how to contribute to documentation or code, or the GitHub [CONTRIBUTING.md](.github/CONTRIBUTING.md) 
for a list of checks when opening a pull request.

## Known Issues 

The following are known to not work, and development will depend on how useful
the average user will assess each of these points. The developer @vsoch has not
added them yet because she doesn't think them overall useful.

 - support for system libraries (e.g., sys) or anything without a filename in site-packages
---
title: 'GridTest: testing and metrics collection for Python'
tags:
  - testing
  - metrics
  - continuous integration
authors:
 - name: Vanessa Sochat
   orcid: 0000-0002-4387-3819
   affiliation: 1
affiliations:
 - name: Stanford University Research Computing
   index: 1
date: 12 May 2020
bibliography: paper.bib
---

## Summary

Creating reproducible testing and parameterization is a consistent challenge
in research because it is hard to do [@Koteska2015ScientificST]. 
Researchers often don't have the bandwidth to think about tests, and
consequently, creating tools for this use case is often overlooked.
GridTest [@gridtest] is a library that specializes in generating parameter grids. The grids
are most obviously used for testing, but can extend to other use cases.
In the context of testing, GridTest makes it easy to discover functions,
classes, and arguments for Python scripts or modules, and then generate
a template to easily populate. Outside of testing, grids can be created
that are version controlled, programatically defined with functions,
and easy to interact with from the command line or Python interpreter.
Grids can be used with tests that can further be parameterized
and configured to collect metrics for each case run. Both grid and test 
specifications are stored in a simple YAML configuration that the library helps to generate,
and features include interactive debugging, interactive report generation,
and provided metrics (Python decorators) that can assist with research.

## Background

While several scientific libraries [@sklearn; @gridregular] make it possible to
generate parameter grids within code, they require a substantial list of
numerical library dependencies that might be overkill for the user's needs,
and further, they don't allow for representation of grids outside of the code.
Additionally, these libraries do not natively allow for interactive debugging,
collection of metrics around grids, or report generation. Another subset
of grid generation libraries are specifically intended for hyperparameter tuning
[@keras; @h2o; @Willkoehrsen2018-hm], and are thus packaged
alongside machine learning libraries. While these libraries are rich and hugely
useful for their intended purposes, none of them present a domain-agnostic,
simple grid definition that can be defined outside of the programming language (e.g., R, Python)
code that uses it. The landscape is missing a library that places grids alongside
code, and can define them without needing it. GridTest offers this ability,
and further, introduces a new paradigm that grids might be shared between
code bases as their own entity, and are not required to be embedded within it.

### Grids as First Class Citizens

Parameters always come as a second thought when writing tests, and this is
why they are commonly applied as decorators. The author of this software
realized that she might want to define just sets of parameters that expand
into matrices that can be useful across many use cases. This makes
the grids "first class citizens." For example, instead of a top to bottom
script that loops over some set of datasets, parameters, and algorithms, 
the user could define grids to generate each in a modular fashion. This
is explained in detail for the [clustering grids](https://github.com/vsoch/gridtest/tree/master/examples/clustering-grids) example derived from scikit-learn. As another example, the user might
just want to parameterize some set of inputs to randomly generate a cohort.
This example is detailed [in another tutorial here](https://vsoch.github.io/gridtest/tutorials/samplegrid/). The overall idea is simple. The current practice is generally to write parameterizations alongside code, whether that means nested for loops or decorators for testing.
GridTest allows for this same kind of functionality, but storing the parameterization
alongside the code and not embedded with it. This makes it easy to change grids
or tests without touching the code.

## Concepts

### Testing

A **gridtest**: is one that is run over a grid of parameter settings. Each test
can include a set of argument specifications, and optionally mapping these arguments
to functions so they can be programatically defined. 
A grid can be inline to the test (if not used elsewhere) or defined globally and shared.
For an example of command line usage, the reader is directed to the ["How does it work"](https://vsoch.github.io/gridtest/getting-started/index.html) section in the Getting Started guide.


### Parameterization

A **grid** is a global definition of a parameter matrix. A user can define arguments,
and optionally functions to run to be mapped to arguments. Grids are generated
on demand, meaning when the user iterates over a Grid object, so that no large lists are stored in memory.
Grids can be put to many uses. The user might share a repository that only defines grids that people
can use across many different kinds of machine learning models, possibly to collect metrics
and compare different analysis strategies being used. An introduction to grids
is available [here](https://vsoch.github.io/gridtest/getting-started/grids/).

### Metrics

A **metric** is a Python decorator that is paired with a test that will measure some
attribute of a test. For example:
 - the user might run a function across a grid of arguments, and then measure the time that each combination takes (the metric), and generate a report for inspection.
 - the user might be doing text processing and having functions to parse text. Each function might be run over a grid of sentences and counts, and for each result, the number of unique and total words is counted (metrics).

Metrics are fully described in the [metrics](https://vsoch.github.io/gridtest/getting-started/metrics/) section
of the documentation.


## Use Cases

GridTest has use cases well beyond testing, because parameterization is used
widely across data science, and version control for reproducibility of those
parameters is essential for reproducible, sustainable science. The following
set of examples are good use cases.

### 1. Parameter Grids

GridTest extends the traditional definition [@sklearn-tutorial] of a grid to include:

 - [generating random samples](https://vsoch.github.io/gridtest/tutorials/samplegrid/)
 - [loading grids via a GridRunner](https://vsoch.github.io/gridtest/getting-started/grids/index.html) class separate from the application's Python code.
 - generating grids as they are needed (meaning as an iterator)
 - previewing grids on the command line before using them
 - generating content of grids via external functions, and optionally unwrapping list values

Grids are generated on demand for more efficient memory allocation, and can be extended to any use case that requires some combination of arguments, and optionally functions to run to be mapped to arguments. See the section on the concept of a grid for more detail. 


### 2. Capturing Metrics

How long does a function take to run when provided parameter X as one value, versus
another? By way of allowing the user to specify one or more metrics alongside tests,
they can easily capture metrics (Python decorators to functions to test)
to output in an interactive report. For example, if the user writes a test that runs
a machine learning algorithm across a grid of datasets and algorithms, they can easily
add a metric to record the time that each takes, and save this result to a file.
GridTest provides a standard set of [decorators](https://vsoch.github.io/gridtest/tutorials/decorators/index.html) for
ease of use, and the user is also free to write their own functions to collect
metrics.

### 3. Generating Reports

It can be handy to save results to a data file (e.g., results.json) or generate
an interactive report for GitHub pages. GridTest allows for this by way of the
`--save` or `--save-web` flags. An example web report is shown in Figure 1. Any grid can also be exported
to JSON for archiving in a repository or extension to other custom visualizations.

![Figure 1. An example GridTest web report](report.png)

An interactive, live report is available to view [here](https://vsoch.github.io/gridtest/templates/report/),
and more information about reports and export formats is provided [here](https://vsoch.github.io/gridtest/getting-started/results/index.html). More domain-specific reports can be developed as requested.

### 4. Debugging

What programmer hasn't been in the scenario of running a group of tests,
and then having some fail? What can be done? The user might start an interactive
shell, import what is needed, and try to reproduce, or they can turn up verbosity
and add a bunch of print statements to figure out what is going on. GridTest makes
this much easier with its `--interactive` mode, which allows the user to 
shell into an interactive session right before the function is run to allow for debugging. 
A detailed walkthrough of debugging is provided [here](https://vsoch.github.io/gridtest/getting-started/debugging/).

### 5. Running Reproducible Tests

When the user writes tests for a file, local, or system module, they store them in
a YAML file that is stored alongside the code, and can be tested with CI.
The YAML file can have grids of parameters defined so the user can easily test many
different combinations.

### 6. Knowing the tests to write

Whether tests are written during development or at the end of it, the user typically needs to look through files
to know the function names and arguments that need to be tested. GridTest solves
this problem by way of discovery - the user can give it a module name, a file name, or
an entire directory with Python files, and it will generate a template
to fill in that already includes arguments and functions.  
Once tests are written, the user can run GridTest with the `--check` feature
to find newly added functions in the code. For more details about creating, checking, and updating
tests, see the [testing](https://vsoch.github.io/gridtest/getting-started/testing/index.html)
documentation.

In summary, GridTest:

 1. Lets the user define grids to be generated programatically, version controlled, and used for multiple purposes
 2. Allows measuring metrics alongside tests
 3. Stores tests in a YAML file that can be stored in version control
 4. Generates data exports and interactive reports for results
 5. Provides an easy way to interactively debug
 6. Helps to discover the tests that need to be written, and creates a template to fill in
 7. Makes it easy to define and interact with expanded parameter grids

## Conclusion

GridTest aims to be a general tool for data scientists and research software engineers
alike. If there is a need to create a collection of grids, regardless of being used
for testing or another use case, GridTest can solve this problem. By way of saving
to YAML, GridTest can represent grids separately from code.
By way of providing JSON export with interactive web
reports, GridTest can be used in a continuous integration setup to generate a report
for some tests or metrics of interest. The author had amazing fun creating this library,
and is excited for its potential generalizability to support many different kinds
of research tasks. For more examples, tutorials, and details, see the official documentation at https://vsoch.github.io/gridtest [@gridtest-docs].

# References
* GridTest version:
* Python version:
* Operating System:

### Description

Describe what you were trying to get done.
Tell us what happened, what went wrong, and what you expected to happen.

### What I Did

```
Paste the command(s) you ran and the output.
If there was a crash, please include the traceback here.
```
# Contributor's Agreement

This code is licensed under the MPL 2.0 [LICENSE](LICENSE). This contributing
guide provides detailed instructions and a checklist for pull requests.
Please also note contributing documentation available [here](https://vsoch.github.io/gridtest/contributing/index.html).

# Contributing

When contributing to GridTest, it is important to properly communicate the
gist of the contribution. If it is a simple code or editorial fix, simply
explaining this within the GitHub Pull Request (PR) will suffice. But if this
is a larger fix or enhancement, it should be first discussed with the project
team.

Please note we have a code of conduct, described below. Please follow it in
all your interactions with the project members and users.

## Pull Request Process

1. All PRs should be against the master branch.
2. Follow the existing code style precedent. We use black for linting, and you
   can install it `pip install black` and then run `black gridtest` to format the code.
   This will be tested.
3. Test your PR locally, and provide the steps necessary to test for the
   reviewers.
4. The project's default copyright and header have been included in any new
   source files.
5. All (major) changes to GridTest must be documented in the docs folder.
   If your PR changes a core functionality, please 
   include clear description of the changes in your PR so that the docs 
   can be updated, or better, submit another PR to update the docs directly.
6. Update the CHANGELOG.md. If necessary, update the README.md.
7. The pull request will be reviewed by others, and the final merge must be
   done by the project lead, @vsoch (or approved by her).


# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

### Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project leader @vsoch. All
complaints will be reviewed and investigated and will result in a response
that is deemed necessary and appropriate to the circumstances. The project
team is obligated to maintain confidentiality with regard to the reporter of
an incident. Further details of specific enforcement policies may be posted
separately.

Project maintainers, contributors and users who do not follow or enforce the
Code of Conduct in good faith may face temporary or permanent repercussions 
with their involvement in the project as determined by the project's leader(s).

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Examples

## Grids

GridTest can be used just to generate grids that are useful for however you please! For example, we can define
and then load a grids.yml file in our Python code to parameterize tests or other functions.

 - [grids](grids): is an example of defining different grids for your use in a grids.yml file.

## Testing

 - [basic](basic): shows writing a grid test file for a basic script
 - [read-write](read-write): shows how to use template substitutes `{% tmp_path %}` and `{% tmp_dir %}` for creating temporary files and directories for tests.
 - [is-true-false](is-true-false): example tests for using the `istrue` and `isfalse` conditions.
 - [package](package): generate for a python package you've already installed
 - [class](class): run grid tests for a car class

## Metric Collection

 - [metrics](metrics): use decorators to measure metrics across a grid of tests (eg., timeit)
 - [custom-decorator](custom-decorator): writing a custom decorator to count words over a grid of parameters
 - [grid-function](grid-function) gridtest can derive input parameters from one or more functions.


## Grids

 - [grids](grids): an introduction to generating and using grids.
 - [sample-grid](sample-grid): using gridtest to generate grids of samples.
 - [clustering-grids](clustering-grids): using GridTest to easily run the [clustering example for scikit-learn](https://scikit-learn.org/stable/auto_examples/cluster/plot_linkage_comparison.html#sphx-glr-auto-examples-cluster-plot-linkage-comparison-py)

Also see the [tutorials](https://vsoch.github.io/gridtest/tutorials/) section of the web documentation for more details.
# istrue and isfalse

This is-true-false example will show how to generate and run a grid test
to use the istrue and isfalse conditional checks.

## Install

After you've installed gridtest:

```bash
git clone git@github.com:vsoch/gridtest
cd gridtest
pip install -e .
```

or

```bash
pip install gridtest
```

## Generate

then you can cd into this folder, and test generating a gridtest file for the
[truefalse.py](truefalse.py) included here:

```bash
$ gridtest generate truefalse.py gridtest.yml
Extracting add from truefalse
Extracting add_with_type from truefalse
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the [gridtest.yml](gridtest.yml) will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case.

## Customize

You can then open the file in a text editor, and add arguments to each.
If your library uses typing, the typing will be checked at testing time,
and it's not specified here. You'll generally want to fill in args for
each testing condition (or leave null for None). For example, we might want to 
change:

```yaml
  script.add:
    args:
    - one: null
    - two: null
```

to instead be:

```yaml
  script.add:
    args:
    - one: 1
    - two: 2
```

To test adding 1+2. 

### Return Types

Since we are primarily interested with the `istrue` and `isfalse` return
types, let's just look at examples for each of those.

**istrue**

istrue is used when we want to check if something is True.
You usually would want to refer to an input or output variable:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    istrue: isinstance({% returns %}, int)
```

**isfalse**

or you might want the opposite, isfalse:


```yaml
  script.add:
  - args:
      one: 1
      two: 2
    isfalse: not isinstance({% returns %}, float)
```

This means that we can edit our script from this:

```yaml
truefalse:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/is-true-false/truefalse.py
  tests:
    truefalse.add:
    - args:
        one: null
        two: null
    truefalse.add_with_type:
    - args:
        one: null
        two: null
```

to be something more reasonable to test:

```yaml
truefalse:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/is-true-false/truefalse.py
  tests:
    truefalse.add:
    - args:
        one: 1.0
        two: 2
      istrue: "isinstance({{ result }}, float)"
      isfalse: "isinstance({{ result }}, int)"
    truefalse.add_with_type:
    - args:
        one: 1
        two: 2
      returns: 3
    - args:
        one: 1.0
        two: 2
      raises: TypeError
```

For typing, given that a function uses typing, that will be tested. For example,
the last function "add_with_type" would raise a TypeError if we give it a float.
This is why we have added a test case for it. Finally, the template strings `{{ result }}`
and `{{ returns }}` are spceial cases. Returns references what you specify your
test to return, and result specifies the result that is actually returned.
We use "result" above because we didn't do any checks for return values for
the same function.

## Test

Finally, you'll have your test file, and an environment where you want to
test. You can run tests like this:

```bash
$ gridtest test gridtest.yml
[4/4] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
truefalse.add.0                success                        istrue isinstance(3.0, float) isfalse isinstance(3.0, int)
truefalse.add.1                success                        equals 1+2                    
truefalse.add_with_type.0      success                        returns 3                     
truefalse.add_with_type.1      success                        raises TypeError              

4/4 tests passed
```

Or since gridtest.yml is the default, just do:

```bash
$ gridtest test
```
# Grid Function

This is a derivative of the [custom-decorator](../custom-decorator)
example, except we derive inputs for the grid from a custom function.

### Create your Functions

Originally we started with functions to take some text input and parse it 
(`multiply_sentences`).

```python
# These are functions in my script

def multiply_sentence(sentence, count):
    return sentence * count
```

But since we want to test having a function generate inputs for our test function,
let's change this up a bit. Instead of generating a sentence of some length,
let's return the ascii characters for a pokemon. The input that we need for the
function below is the pokemon id (pid).

```python
# These are functions in my script

from pokemon.master import get_pokemon, catch_em_all

def generate_pokemon(pid):
    """Generate a pokemon based on a particular identifier. This is excessive
       because we could just use get_pokemon() to return a random pokemon, but
       we are taking in the pid (pokemon id) as an example for providing a function
       to generate input arguments for a gridtest.
    """
    catch = get_pokemon(pid=pid)
    catch_id = list(catch.keys())[0]
    return catch[catch_id]['ascii']
```

How will we generate that?

### Create your Yaml

Remember that for the example in the [custom-decorator](../custom-decorator), 
we just provide a list of sentences for our args:

```yaml
  script.multiply_sentence:
  - metrics:
    - '@script.countwords'
    args:
      count:
        list: [1, 5, 10]
      sentence:
        list:
          - "He ran for the hills."
          - "Skiddery-a rinky dinky dinky, skittery rinky doo."
          - "You are my sunshine, my only sunshine."
```

This seems reasonable for a dummy example, or small set of inputs, but likely 
isn't reasonable if we want to more programatically generated inputs.
For example, for our current pokemon function we would need to hard
code a specific set of ids (as a list):

```
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  tests:
    script.generate_pokemon:
    - metrics:
      - '@script.countwords'
      args:
        pid: 
          list: [1, 2, 3]
```

or some range:

```
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  tests:
    script.generate_pokemon:
      args:
        pid: 
          min: 1
          max: 3
```

What we really want to do is write a function that randomly generates pokemon ids for us.
Let's create another function that we want to run many times to generate
pokemon ids for this function:

```python
def get_pokemon_id():
    """Return a random pokemon id"""
    return random.choice(list(catch_em_all()), 1)
```

A more realistic use case would be having some numerical optimization function
that requires a distribution (an array of numbers) to be generated for some input.
We'd want to generate that array with a function over some grid of parameters.

### Add Arguments

We basically need to point the parameter "pid" to be generated from the function. What
does that look like? First we need to add the function under "grids" - a grid
is a general parameterization that can be used to populate sections of your testing
file. We will give it a name, generate_pids for "generate pokemon ids."

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  grids:
    generate_pids:
      functions: 
        pid: script.get_pokemon_id
...
```

The function doesn't have any input arguments, so we can just specify it's name
under functions, and notice that we are mapping it to populate the variable "pid."
But not we run into another issue - how do we tell the grid to be run 10 times
(each to randomly generate a pokemon id?) To do this, we add a "count" to our
grid:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  grids:
    generate_pids:
      count: 10
      functions: 
        pid: 
          func: script.get_pokemon_id
```

Then we need to specify our generate_pokemon function to use it, and we do this
by adding the reference to the grid "generate_pids" under `script.generate_pokemon`:

```yaml
  tests:
    script.generate_pokemon:
      - metrics:
          - "@script.uniquechars"
        grid: generate_pids
```

The above says that "for script.generate_pokemon test, use the grid named
"generate_pids" to parameterize and populate input parameters. We can then
pop up to the grids section to see that we will generate pid by running
a function 10 times.

> Why can't I specify the arguments alongside the test?

You can! The above recipe is actually defining a global grid to use, which
might be optimal if you want to share grids between tests. But you can
just as easily define an inline grid. That would look like this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  tests:
    script.generate_pokemon:
      - metrics:
          - "@script.uniquechars"
        count: 10
        functions: 
          pid: 
            func: script.get_pokemon_id
```

And is a more succinct (but possibly redundant) method to run the same thing.
See [gridtest-inline.yml](gridtest-inline.yml) for this example, and to run:

```bash
$ gridtest test gridtest-inline.yml
```

### Add a Metric

Finally, let's create a metric to simply take the ascii (the result of `generate_pokemon`
and count the number of unique characters:

```python
def uniquechars(func):
    """this is a simple example of a custom decorator - We run the generate
       pokemon function and print the number of unique characters in the ascii.
    """
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        chars = len(set(result))
        print(f"@script.uniquechars {chars} unique chars")
        return result

    return wrapper
```

The final template (with the metric) looks like this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  grids:
    generate_pids:
      func: script.get_pokemon_id
      count: 10

  tests:
    script.generate_pokemon:
      - metrics:
          - "@script.uniquechars"
        grid:
          pid: generate_pids
```

### Running Tests

Let's run the tests! We should see a count of unique characters for each pokemon generated:

```bash
$ gridtest test
[10/10] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.generate_pokemon.0      success                                                      
script.generate_pokemon.1      success                                                      
script.generate_pokemon.2      success                                                      
script.generate_pokemon.3      success                                                      
script.generate_pokemon.4      success                                                      
script.generate_pokemon.5      success                                                      
script.generate_pokemon.6      success                                                      
script.generate_pokemon.7      success                                                      
script.generate_pokemon.8      success                                                      
script.generate_pokemon.9      success                                                      

________________________________________________________________________________________________________________________
script.generate_pokemon.0      @script.uniquechars            11 unique chars               
script.generate_pokemon.1      @script.uniquechars            11 unique chars               
script.generate_pokemon.2      @script.uniquechars            11 unique chars               
script.generate_pokemon.3      @script.uniquechars            11 unique chars               
script.generate_pokemon.4      @script.uniquechars            11 unique chars               
script.generate_pokemon.5      @script.uniquechars            11 unique chars               
script.generate_pokemon.6      @script.uniquechars            11 unique chars               
script.generate_pokemon.7      @script.uniquechars            11 unique chars               
script.generate_pokemon.8      @script.uniquechars            11 unique chars               
script.generate_pokemon.9      @script.uniquechars            11 unique chars               

10/10 tests passed
```

Spoiler alert - they all use the same characters! Notice that we've also successfully
run the unique id generator function 10 times to result in 10 tests.
See the [gridtest.yml](gridtest.yml) for the full test recipe.

### Save to File

If we want to save the complete results, we can add `--save` with a filename:

```bash
$ gridtest test --save results.json
```

You can see the example [results.json](results.json) in this folder, and notice
how the result of each run is the ascii for a complete pokemon (that would
render nicely if printed to the screen).
# Metrics

The most simple case of creating a GridTest is to do so by running a set of 
metrics (which are actually Python decorators) for a testing
file, and then setting ranges of variables to run within tests. For example,
let's start with the [is-true-false](../is-true-false/) example that has functions
for addition, and customize the file to have metric.

## Adding Metrics

The first thing to do is to add a metrics section. A metric can be any
custom python decorator that you can easily import. As a good starting example,
we can use a custom timeit provided in `gridtest.decorators`.
We would add this to the functions we want to measure:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/optimize/script.py

  tests:
    script.add:
      # metrics are each a decorator, we first look to gridtest.decorators then external import
    - metrics:
        - "@timeit"
```

The `@` indicates that we are using a decorator, which is the current convention
in case other kinds of metrics might be used. Once we do this, when we run gridtest
with the file, it will by default detect the metric, and then run the functions
through it.

When you run the test, you'll see the results for metric decorators at the bottom.



```bash
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        istrue isinstance(self.result, float) isfalse isinstance(self.result, int)               
script.add.1                   success                        equals 1+2                                                   
script.add_with_type.0         success                        returns 3                                                    
script.add_with_type.1         success                        raises TypeError                                             

________________________________________________________________________________________________________________________
script.add.0                   @timeit                        0.00 ms                       

4/4 tests passed
```

## Adding Arguments

Great - so we've measured time for one function. What if we want to measure the time for
a function, but across a parameter grid? We might want to adopt our recipe to allow for this:

```yaml
    args:
      one:
        min: 0
        max: 5
        list: [10, 15]
    args:
      two: 2
```

but this won't produce a very interesting result, because the result is just adding
two numbers. Let's try something that will work with out timeit function, namely
write a function that will sleep for some number of input seconds. We would then want to see
that the timeout output increases to match the input seconds.

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/optimize/script.py
  tests:
    script.gotosleep:
    - metrics:
      - '@timeit'
      args:
        seconds:
          list: [10, 15]
          max: 5
          min: 0
```

Before we run the test, let's talk about the formatting of the yaml.
Notice that we've moved the "one" argument up from the "args" section into the grid
section. This tells gridtest that we want to run a grid of tests for some 
parameterization of our argument "one." In the example above, we want to include
a range from 0 to 5 (default increment is by 1) and also add the one off values of 10 and 15 
(provided in list). This gives us the following allowable keys for a grid
parameter:

**min** and **max** are to be used when specifying a range. When unset, **by** 
would be 1. If you want to decrease, set a negative value for by.

**list** is for when you want to include a list of values, even in addition to a
range already specified as in the example above.

## Run the GridTest

Now let's run the test!

```bash
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.gotosleep.0             success                                                      
script.gotosleep.1             success                                                      
script.gotosleep.2             success                                                      
script.gotosleep.3             success                                                      
script.gotosleep.4             success                                                      
script.gotosleep.5             success                                                      
script.gotosleep.6             success                                                      
script.add.0                   success                        istrue isinstance(self.result, float) isfalse isinstance(self.result, int)
script.add.1                   success                        equals 1+2                    
script.add_with_type.0         success                        returns 3                     
script.add_with_type.1         success                        raises TypeError              

________________________________________________________________________________________________________________________
script.gotosleep.0             @timeit                        0.01 ms                       
script.gotosleep.1             @timeit                        1000.77 ms                    
script.gotosleep.2             @timeit                        2001.84 ms                    
script.gotosleep.3             @timeit                        3001.72 ms                    
script.gotosleep.4             @timeit                        4003.81 ms                    
script.gotosleep.5             @timeit                        10009.59 ms                   
script.gotosleep.6             @timeit                        15013.49 ms                   

11/11 tests passed
```

Awesome! We've run a grid of tests over the different values for seconds, and have 
reported the total time taken via the timeit decorator. See the [gridtest.yml](gridtest.yml)
for the full test recipe.
# Clustering Grids

This small tutorial will show how to use Gridtest to generate grids to run the
[scikit-learn](https://scikit-learn.org/stable/auto_examples/cluster/plot_linkage_comparison.html#sphx-glr-auto-examples-cluster-plot-linkage-comparison-py) cluster plot linkage comparison.
There are three parts:

## Dataset and Algorithms Grids

[1-dataset-algorithms-grids](1-dataset-algorithms-grids/) walks through
converting the original script to have functions that populate dataset and
algorithms grids that are parameterized to run over a single function
to plot one dataset. This example serves to show the basic steps and thinking
to convert a top to bottom script to a simple grid design.

## Modular Grids

[2-modular-grids](2-modular-grids) takes the original design, and improves
upon it by making the datasets and algorithms completely modular, meaning
that we can run each dataset/algorithm combination as a separate task, and
add metrics to compare performance.

# Modular Clustering Grids

We are starting with the [1-dataset-algorithms-grids](../1-dataset-algorithms-grids)
example and modifying it to be completely modular - meaning that the algorithms
and datasets are completely separated, and the the logic looks like this:

```
for every dataset and parameter combination
    for every algorithm and parameter combination
       generate a plot along with metrics
```

The idea would be that each grid in our [grids.yml](grids.yml) is modular, so
I might use either analysis or dataset for another analysis, or write
different tests that manipulate the grids. Let's start with updating our
functions.

## Updating the Grids

### Algorithm Generation

Let's start with the function to generate algorithms. If you remember from
the previous example, the function took in a dataset as a dependency. But
now we want to eliminate that, so we will instead input just the number of clusters.

```bash
def generate_algorithms(n_clusters=2):
    """Here we have updated the original generate_algorithms function to separate
       the dataset from it, so it can be run in parallel, and the number of clusters
       variable varied.
    """
    ward = cluster.AgglomerativeClustering(
        n_clusters=n_clusters, linkage='ward')
    complete = cluster.AgglomerativeClustering(
        n_clusters=n_clusters, linkage='complete')
    average = cluster.AgglomerativeClustering(
        n_clusters=n_clusters, linkage='average')

    return [
        ('Average Linkage', average),
        ('Complete Linkage', complete),
        ('Ward Linkage', ward),
    ]
```

And then we can vary the parameter for the number of clusters. The following says to
run the function for value of n_clusters, and then "unwrap" the result. This should
mean that the variable "algorithms" is a list of 3 algorithms x 2 args, for a total of
6:

```bash
$ gridtest gridview grids.yml generate_algorithms --arg algorithms
[('Average Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='average', memory=None,
            n_clusters=2, pooling_func=<function mean at 0x7f3705f50cb0>)), ('Complete Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='complete', memory=None,
            n_clusters=2, pooling_func=<function mean at 0x7f3705f50cb0>)), ('Ward Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='ward', memory=None, n_clusters=2,
            pooling_func=<function mean at 0x7f3705f50cb0>)), ('Average Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='average', memory=None,
            n_clusters=5, pooling_func=<function mean at 0x7f3705f50cb0>)), ('Complete Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='complete', memory=None,
            n_clusters=5, pooling_func=<function mean at 0x7f3705f50cb0>)), ('Ward Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='ward', memory=None, n_clusters=5,
            pooling_func=<function mean at 0x7f3705f50cb0>))]
```
```bash
$ gridtest gridview grids.yml generate_algorithms --arg algorithms --count
Variable algorithms has length 6.
```

Great! Now let's review our datasets.

### Dataset Generation

We actually don't need to update the datasets grid, because it already generates
a list of datasets that we can use. As a reminder, here is the grid definition,
which also uses unwrap for the variable:

```yaml
    generate_datasets:
      functions:
        datasets:
           func: analysis.generate_datasets
           unwrap: true
```

And we can confirm that 6 are generated:

```bash
$ gridtest gridview grids.yml generate_datasets --arg datasets --count
Variable datasets has length 6.
```

However we will remove the n_clusters (and params provided with the dataset)
since they aren't used anymore. The new function just returns the dataset
with it's name. The list of tuples will be unwrapped so each dataset is treated
as a single unit in the parameterization.

```python
def generate_datasets(n_samples=1500, factor=0.5, noise=0.05):
    """This generate_datasets function is simply taking the original 
       (top to bottom) style code, and converting into a function to
       return datasets. This function can be provided to a grid, and then
       each dataset will be run across some number of variables to produce
       a grid
    """
    noisy_circles = datasets.make_circles(n_samples=n_samples, factor=factor,
                                          noise=noise)
    noisy_moons = datasets.make_moons(n_samples=n_samples, noise=.05)
    blobs = datasets.make_blobs(n_samples=n_samples, random_state=8)
    no_structure = np.random.rand(n_samples, 2), None

    # Anisotropicly distributed data
    random_state = 170
    X, y = datasets.make_blobs(n_samples=n_samples, random_state=random_state)
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    aniso = (X_aniso, y)

    # blobs with varied variances
    varied = datasets.make_blobs(n_samples=n_samples,
                                 cluster_std=[1.0, 2.5, 0.5],
                                 random_state=random_state)

    # name, dataset
    return [
        ("circles", noisy_circles),
        ("moons", noisy_moons),
        ("varied", varied),
        ("aniso", aniso),
        ("blobs", blobs),
        ("no-structure", no_structure)]
```

### Plot Generation

We are almost done! We now want to run a test that will iterate over the list
of algorithms, datasets, and any additional arguments that we provide. First
we can make a grid to parameterize the two:

```yaml
    generate_inputs:
      ref:
        algorithm: generate_algorithms.algorithms
        dataset: generate_datasets.datasets
```


Since we know the list of algorithms produces 6 (3 algorithms over 2 cluster sizes)
and there are 6 datasets, we should have a total of 36 in the parameterized result:

```bash
$ gridtest gridview grids.yml generate_inputs --count
36 argument sets produced.
```

And we do! Now we can update the analysis.generate_plot function (note I've
renamed it from generate_plots to generate_plot since we are just running
one dataset for one algorithm now!) It just points to the `generate_inputs` grid - and we updated this from
previously when it pointed to the generate_algorithms grid, since we had not
separated out the algorithms.

```yaml
    # To be used in this test.
    analysis.generate_plot:
      - grid: generate_inputs
```

### Testing the Test

Let's run the test! Actually, first let's run an interactive test so we can see
the arguments being handed to the function, and debug or develop if necessary:

```bash
$ gridtest test grids.yml --interactive
[analysis.generate_plot:1/36] ||----------------------------------|   2.8% 

Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function generate_plot at 0x7f5e3dbc6b00>
 module: analysis
   args: {'algorithm': ('Average Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='average', memory=None,
            n_clusters=2, pooling_func=<function mean at 0x7f5e4c047cb0>)), 'dataset': ('circles', (array([[-0.67799938, -0.69875698],
       [ 0.93143746,  0.19139133],
       [ 0.54829131, -0.00601715],
       ...,
       [-0.34518816, -0.35804797],
       [ 0.01719727, -0.94513802],
       [ 0.91377877, -0.59884164]]), array([0, 0, 1, ..., 1, 0, 0])), {'n_clusters': 2})}
returns: None

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
```

The args look as we expect, a dictionary (key word arguments or kwargs) that
will be exploded into the function, which we also see 

```bash
Variables
   func: <function generate_plot at 0x7f5e3dbc6b00>
```

We can see `How to test` in the terminal as well, and quit with Control+D when we finish.

### Running the Test

```bash
$ gridtest test grids.yml
$ gridtest test grids.yml
[36/36] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
analysis.generate_plot.0       success                                                      
analysis.generate_plot.1       success                                                      
analysis.generate_plot.2       success                                                      
analysis.generate_plot.3       success                                                      
analysis.generate_plot.4       success                                                      
analysis.generate_plot.5       success                                                      
analysis.generate_plot.6       success                                                      
analysis.generate_plot.7       success                                                      
analysis.generate_plot.8       success                                                      
analysis.generate_plot.9       success                                                      
analysis.generate_plot.10      success                                                      
analysis.generate_plot.11      success                                                      
analysis.generate_plot.12      success                                                      
analysis.generate_plot.13      success                                                      
analysis.generate_plot.14      success                                                      
analysis.generate_plot.15      success                                                      
analysis.generate_plot.16      success                                                      
analysis.generate_plot.17      success                                                      
analysis.generate_plot.18      success                                                      
analysis.generate_plot.19      success                                                      
analysis.generate_plot.20      success                                                      
analysis.generate_plot.21      success                                                      
analysis.generate_plot.22      success                                                      
analysis.generate_plot.23      success                                                      
analysis.generate_plot.24      success                                                      
analysis.generate_plot.25      success                                                      
analysis.generate_plot.26      success                                                      
analysis.generate_plot.27      success                                                      
analysis.generate_plot.28      success                                                      
analysis.generate_plot.29      success                                                      
analysis.generate_plot.30      success                                                      
analysis.generate_plot.31      success                                                      
analysis.generate_plot.32      success
analysis.generate_plot.33      success                                                      
analysis.generate_plot.34      success                                                      
analysis.generate_plot.35      success                                                      
```

The images generated are shown in the present working directory here. Here is one example:

![no-structure-complete-linkage.png](no-structure-complete-linkage.png)

### Adding Metrics

This is where it gets a little more interesting! We can add [metrics](https://vsoch.github.io/gridtest/getting-started/metrics/) to compare across the runs. As an easy example, let's time each run - the @timeit function
is provded by gridtest. 


If you write your own metric or want to use one (a decorator) that already exists, you can
just reference it as the importable path.

```yaml
    # To be used in this test.
    analysis.generate_plot:
      - grid: generate_inputs
        metrics: ["@timeit"]
```

And we now get a time for each algorithm and dataset generation:

```bash
$ gridtest test grids.yml[36/36] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
analysis.generate_plot.0       success                                                      
analysis.generate_plot.1       success                                                      
analysis.generate_plot.2       success                                                      
analysis.generate_plot.3       success                                                      
analysis.generate_plot.4       success                                                      
analysis.generate_plot.5       success                                                      
analysis.generate_plot.6       success                                                      
analysis.generate_plot.7       success                                                      
analysis.generate_plot.8       success                                                      
analysis.generate_plot.9       success                                                      
analysis.generate_plot.10      success                                                      
analysis.generate_plot.11      success                                                      
analysis.generate_plot.12      success                                                      
analysis.generate_plot.13      success                                                      
analysis.generate_plot.14      success                                                      
analysis.generate_plot.15      success                                                      
analysis.generate_plot.16      success                                                      
analysis.generate_plot.17      success                                                      
analysis.generate_plot.18      success                                                      
analysis.generate_plot.19      success                                                      
analysis.generate_plot.20      success                                                      
analysis.generate_plot.21      success                                                      
analysis.generate_plot.22      success                                                      
analysis.generate_plot.23      success                                                      
analysis.generate_plot.24      success                                                      
analysis.generate_plot.25      success                                                      
analysis.generate_plot.26      success                                                      
analysis.generate_plot.27      success                                                      
analysis.generate_plot.28      success                                                      
analysis.generate_plot.29      success                                                      
analysis.generate_plot.30      success                                                      
analysis.generate_plot.31      success                                                      
analysis.generate_plot.32      success                                                      
analysis.generate_plot.33      success                                                      
analysis.generate_plot.34      success                                                      
analysis.generate_plot.35      success                                                      

________________________________________________________________________________________________________________________
analysis.generate_plot.0       @timeit                        4836.47 ms                    
analysis.generate_plot.1       @timeit                        4832.80 ms                    
analysis.generate_plot.2       @timeit                        3736.94 ms                    
analysis.generate_plot.3       @timeit                        4345.53 ms                    
analysis.generate_plot.4       @timeit                        4861.36 ms                    
analysis.generate_plot.5       @timeit                        3906.09 ms                    
analysis.generate_plot.6       @timeit                        4704.65 ms                    
analysis.generate_plot.7       @timeit                        5144.60 ms                    
analysis.generate_plot.8       @timeit                        4780.60 ms                    
analysis.generate_plot.9       @timeit                        4818.56 ms                    
analysis.generate_plot.10      @timeit                        4808.86 ms                    
analysis.generate_plot.11      @timeit                        5685.36 ms                    
analysis.generate_plot.12      @timeit                        4822.39 ms                    
analysis.generate_plot.13      @timeit                        5762.99 ms                    
analysis.generate_plot.14      @timeit                        5397.26 ms                    
analysis.generate_plot.15      @timeit                        5839.27 ms                    
analysis.generate_plot.16      @timeit                        5394.22 ms                    
analysis.generate_plot.17      @timeit                        2114.58 ms                    
analysis.generate_plot.18      @timeit                        1972.74 ms                    
analysis.generate_plot.19      @timeit                        1901.97 ms                    
analysis.generate_plot.20      @timeit                        1706.98 ms                    
analysis.generate_plot.21      @timeit                        2070.13 ms                    
analysis.generate_plot.22      @timeit                        2010.58 ms                    
analysis.generate_plot.23      @timeit                        2021.10 ms                    
analysis.generate_plot.24      @timeit                        1966.76 ms                    
analysis.generate_plot.25      @timeit                        1980.24 ms                    
analysis.generate_plot.26      @timeit                        2005.79 ms                    
analysis.generate_plot.27      @timeit                        1383.97 ms                    
analysis.generate_plot.28      @timeit                        1851.91 ms                    
analysis.generate_plot.29      @timeit                        1732.63 ms                    
analysis.generate_plot.30      @timeit                        1798.06 ms                    
analysis.generate_plot.31      @timeit                        1417.78 ms                    
analysis.generate_plot.32      @timeit                        1510.49 ms                    
analysis.generate_plot.33      @timeit                        1465.70 ms                    
analysis.generate_plot.34      @timeit                        1325.80 ms                    
analysis.generate_plot.35      @timeit                        1467.37 ms                    

36/36 tests passed
```

We then might just want to save metrics:

```bash
$ gridtest test grids.yml --save-metrics metrics.json
```

For the full results, most of the time json wouldn't be feasible, so you
can save to pickle:

```bash
$ gridtest test grids.yml --save results.pkl
```

### Summary

Now that we've made our grids modular, it would be very easy to use the datasets
OR the algorithms for another purpose. We could even read them into Python
for interactive use, or provide the grids.yml and script for others
to more easily develop other functions that use them.
# Clustering Grids

This small tutorial will show how to use Gridtest to generate grids to run the
[scikit-learn](https://scikit-learn.org/stable/auto_examples/cluster/plot_linkage_comparison.html#sphx-glr-auto-examples-cluster-plot-linkage-comparison-py) cluster plot linkage comparison.

## Developing the Grids

When you look at the original code, it's a bit "top to bottom" (yes I'm avoiding
calling it spaghetti code!) The first thing we want to do is think about how
we would turn this into neater functions, and further, have those functions
easy to plug into a grid. 

### Dataset Generation

The first section is easy - we see a bunch of datasets being produced, and then we want to run those
datasets across a set of models. This is perfect! We can write a `generate_datasets`
function that will take in some input parameters and be used to generate a 
grid. Here is what that looks like - this was largely adding intendation
to the already existing code.

```python
def generate_datasets(n_samples=1500, factor=0.5, noise=0.05):
    """This generate_datasets function is simply taking the original 
       (top to bottom) style code, and converting into a function to
       return datasets. This function can be provided to a grid, and then
       each dataset will be run across some number of variables to produce
       a grid
    """
    noisy_circles = datasets.make_circles(n_samples=n_samples, factor=factor,
                                          noise=noise)
    noisy_moons = datasets.make_moons(n_samples=n_samples, noise=.05)
    blobs = datasets.make_blobs(n_samples=n_samples, random_state=8)
    no_structure = np.random.rand(n_samples, 2), None

    # Anisotropicly distributed data
    random_state = 170
    X, y = datasets.make_blobs(n_samples=n_samples, random_state=random_state)
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    aniso = (X_aniso, y)

    # blobs with varied variances
    varied = datasets.make_blobs(n_samples=n_samples,
                                 cluster_std=[1.0, 2.5, 0.5],
                                 random_state=random_state)

    return [
        (noisy_circles, {'n_clusters': 2}),
        (noisy_moons, {'n_clusters': 2}),
        (varied, {'n_neighbors': 2}),
        (aniso, {'n_neighbors': 2}),
        (blobs, {}),
        (no_structure, {})]
```

At this point we have enough for a very simple grid. Let's write that now:

```yaml
analysis:
  grids:
    generate_datasets:
      functions:
        datasets: analysis.generate_datasets
```

In the above, we've created a grid called `generate_datasets`. We have 
define one variables that is generated by the function `analysis.generate_datasets`.
From the above we can see this function produces a list of datasets, and this list
will be piped into the dataset variable. We can check this:

```bash
gridtest gridview grids.yml generate_datasets --count
1 argument sets produced.
```
```bash
$ gridtest gridview grids.yml generate_datasets                                                                                                                      
{'datasets': [((array([[-0.67799938, -0.69875698],
       [ 0.93143746,  0.19139133],
       [ 0.54829131, -0.00601715],
       ...,
       [-0.34518816, -0.35804797],
       [ 0.01719727, -0.94513802],
       [ 0.91377877, -0.59884164]]), array([0, 0, 1, ..., 1, 0, 0])), {'n_clusters': 2}), ((array([[ 0.49627131, -0.34275349],
       [-0.16629956,  0.92234209],
       [ 0.71895601,  0.66529038],
       ...,
       [ 1.90950927,  0.02989686],
       [ 0.54623069, -0.36003133],
       [ 0.04090016,  0.37069297]]), array([1, 0, 0, ..., 1, 1, 1])), {'n_clusters': 2}), ((array([[ -6.11119721,   1.47153062],
       [ -7.49665361,   0.9134251 ],
       [-10.84489837,  -7.55352273],
       ...,
       [  1.64990343,  -0.20117787],
       [  0.79230661,   0.60868888],
       [  1.91226342,   0.25327399]]), array([1, 1, 0, ..., 2, 2, 2])), {'n_neighbors': 2}), ((array([[-3.37561542,  3.63236314],
       [-3.61882807,  3.78627892],
       [-3.48552993,  0.46412084],
       ...,
       [ 1.17962827, -1.54262502],
       [-0.49738132,  0.78227797],
       [ 1.13089877, -1.13033403]]), array([1, 1, 0, ..., 2, 2, 2])), {'n_neighbors': 2}), ((array([[ 5.86749807,  8.17715188],
       [ 5.61369982,  9.93295527],
       [ 7.22508428, 10.44886194],
       ...,
       [ 7.73674097, 10.82855388],
       [-4.61701094, -9.64855983],
       [-3.48640175, -9.25766922]]), array([0, 0, 0, ..., 0, 2, 2])), {}), ((array([[0.59945663, 0.24694133],
       [0.5173267 , 0.57255303],
       [0.55229185, 0.40567924],
       ...,
       [0.8384347 , 0.52906874],
       [0.84228843, 0.11517496],
       [0.91963613, 0.22592146]]), None), {})]}
```

But actually, we would want to unwrap this list so that each of the datasets
can be parsed by some other function. Let's do that now by adding "unwrap" set
to true to our grid:


```yaml
analysis:
  grids:
    generate_datasets:
      functions:
        datasets:
           func: analysis.generate_datasets
           unwrap: true
```

Do we generate 6 separate datasets?

```bash
$ gridtest gridview grids.yml generate_datasets
{'datasets': ((array([[-0.67799938, -0.69875698],
       [ 0.93143746,  0.19139133],
       [ 0.54829131, -0.00601715],
       ...,
       [-0.34518816, -0.35804797],
       [ 0.01719727, -0.94513802],
       [ 0.91377877, -0.59884164]]), array([0, 0, 1, ..., 1, 0, 0])), {'n_clusters': 2})}
{'datasets': ((array([[ 0.49627131, -0.34275349],
       [-0.16629956,  0.92234209],
       [ 0.71895601,  0.66529038],
       ...,
       [ 1.90950927,  0.02989686],
       [ 0.54623069, -0.36003133],
       [ 0.04090016,  0.37069297]]), array([1, 0, 0, ..., 1, 1, 1])), {'n_clusters': 2})}
{'datasets': ((array([[ -6.11119721,   1.47153062],
       [ -7.49665361,   0.9134251 ],
       [-10.84489837,  -7.55352273],
       ...,
       [  1.64990343,  -0.20117787],
       [  0.79230661,   0.60868888],
       [  1.91226342,   0.25327399]]), array([1, 1, 0, ..., 2, 2, 2])), {'n_neighbors': 2})}
{'datasets': ((array([[-3.37561542,  3.63236314],
       [-3.61882807,  3.78627892],
       [-3.48552993,  0.46412084],
       ...,
       [ 1.17962827, -1.54262502],
       [-0.49738132,  0.78227797],
       [ 1.13089877, -1.13033403]]), array([1, 1, 0, ..., 2, 2, 2])), {'n_neighbors': 2})}
{'datasets': ((array([[ 5.86749807,  8.17715188],
       [ 5.61369982,  9.93295527],
       [ 7.22508428, 10.44886194],
       ...,
       [ 7.73674097, 10.82855388],
       [-4.61701094, -9.64855983],
       [-3.48640175, -9.25766922]]), array([0, 0, 0, ..., 0, 2, 2])), {})}
{'datasets': ((array([[0.59945663, 0.24694133],
       [0.5173267 , 0.57255303],
       [0.55229185, 0.40567924],
       ...,
       [0.8384347 , 0.52906874],
       [0.84228843, 0.11517496],
       [0.91963613, 0.22592146]]), None), {})}
```
```bash
gridtest gridview grids.yml generate_datasets --count
6 argument sets produced.
```

If we wanted to add a variable for each generation of a dataset,
we would add it under args. As an example, let's add a variable `n_samples`
to vary the number of samples. Let's say we want to do each of 1500 and
3000. We can add an arg that will be added to the parameterization:

```python
analysis:
  grids:

    generate_datasets:
      args:
        n_samples: [1500, 3000]
      functions:
        dataset:
          func: analysis.generate_datasets
          unwrap: true
        algorithms: analysis.generate_algorithms
```

If this works as expected, we would want 6 datasets (the result of unwrap) x2 (one
for each of the number of samples (12). We can ask explicitly to see the variable "dataset"
with `--arg` (truncated for readability):

```bash
gridtest gridview grids.yml generate_datasets --arg datasets
[((array([[-0.67799938, -0.69875698],
       [ 0.93143746,  0.19139133],
       [ 0.54829131, -0.00601715],
       ...,
       [-0.34518816, -0.35804797],
       [ 0.01719727, -0.94513802],
       [ 0.91377877, -0.59884164]]), array([0, 0, 1, ..., 1, 0, 0])), {'n_clusters': 2}), ((array([[ 0.49627131, -0.34275349],
       [-0.16629956,  0.92234209],
       [ 0.71895601,  0.66529038],
       ...,
       ...
       ...,
       [0.07990429, 0.71819444],
       [0.38310211, 0.44144583],
       [0.75244099, 0.22064084]]), None), {})]

```
```bash
$ gridtest gridview grids.yml generate_datasets --arg datasets --count
Variable datasets has length 12.
```

We won't use this for this simple example, but will add additional args 
later in this example. Let's start with our original 6 datasets, and
generate algorithms for each.

### Algorithm Generation

The algorithms are generated from the datasets, and we can again write a simple
function to do this:

We naturally want to do the same thing for the set of algorithms that are looped
over. In this case, each algorithm takes a dataset as input.

```python
def get_algorithms(dataset):
    """this was previously in the nested for loop, and returns a tuple of
       algorithms. Since params is derived from a dataset (the index 1) we
       take a dataset as input and then grab the params from it.
    """
    params = dataset[1]
    ward = cluster.AgglomerativeClustering(
        n_clusters=params['n_clusters'], linkage='ward')
    complete = cluster.AgglomerativeClustering(
        n_clusters=params['n_clusters'], linkage='complete')
    average = cluster.AgglomerativeClustering(
        n_clusters=params['n_clusters'], linkage='average')

    return (
        ('Average Linkage', average),
        ('Complete Linkage', complete),
        ('Ward Linkage', ward),
    )
```

The design is really up to you - in this case I chose to pair each dataset
with generation of 3 algorithms since the original code was modeled in that
form. However you could also have each algorithm generation be a separate function.
Let's add a grid that will run this algorithm for each dataset. What set of
datasets? Let's just make a reference to the 6 datasets that we just made.

```yaml
    generate_algorithms:
      ref:
        dataset: generate_datasets.datasets
      functions:
        algorithms: analysis.generate_algorithms
```

And we would expect the 6 datasets to be run through the function once,
generating 6 items, each with algorithms added:

```bash
$ gridtest gridview grids.yml generate_algorithms 
{'dataset': ('blobs', (array([[ 8.26573117,  0.56159687],
       [ 8.70688264, 10.68949521],
       [ 9.32813486,  0.17474761],
       ...,
       [ 8.33646324,  0.0705084 ],
       [ 7.16968586,  8.66745049],
       [ 8.01727721,  0.21072996]]), array([1, 0, 1, ..., 1, 0, 1])), {}), 'algorithms': [('Average Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='average', memory=None,
            n_clusters=3, pooling_func=<function mean at 0x7f049be74cb0>)), ('Complete Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='complete', memory=None,
            n_clusters=3, pooling_func=<function mean at 0x7f049be74cb0>)), ('Ward Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='ward', memory=None, n_clusters=3,
            pooling_func=<function mean at 0x7f049be74cb0>))]}
{'dataset': ('no-structure', (array([[0.25114494, 0.74550472],
       [0.40554956, 0.34608139],
       [0.72326419, 0.8238374 ],
       ...,
       [0.07990429, 0.71819444],
       [0.38310211, 0.44144583],
       [0.75244099, 0.22064084]]), None), {}), 'algorithms': [('Average Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='average', memory=None,
            n_clusters=3, pooling_func=<function mean at 0x7f049be74cb0>)), ('Complete Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='complete', memory=None,
            n_clusters=3, pooling_func=<function mean at 0x7f049be74cb0>)), ('Ward Linkage', AgglomerativeClustering(affinity='euclidean', compute_full_tree='auto',
            connectivity=None, linkage='ward', memory=None, n_clusters=3,
            pooling_func=<function mean at 0x7f049be74cb0>))]}
       ...
```

The count is 6 as we expect.

```bash
$ gridtest gridview grids.yml generate_algorithms --count
6 argument sets produced.
```

### Test Generation

We need our final function that takes in two arguments, a dataset and a list
of algorithms, and generates an output image. The function is basically
the rest of th eoriginal script wrapped in a function, and the test in the
file looks like this:

```yaml
  tests:
    analysis.generate_plots:
      - grid: generate_algorithms
```

In the above, we tell the `analysis.generate_plots` function that takes
as input a dataset and set of algorithms to use the grid `generate_algorithms`
which has parameter grids with each of `dataset` and `algorithms` defined.
We can then run tests, and see an output png generated for each dataset,
across all algorithms.

```
gridtest test grids.yml                                                                                                 
[6/6] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
analysis.generate_plots.0      success                                                      
analysis.generate_plots.1      success                                                      
analysis.generate_plots.2      success                                                      
analysis.generate_plots.3      success                                                      
analysis.generate_plots.4      success                                                      
analysis.generate_plots.5      success                                                      
```

Here is an example!

![linkage-no-structure.png](linkage-no-structure.png)

And our final grids.yml looks like this:

```yaml
analysis:
  grids:
   
    generate_datasets:
      functions:
        datasets:
           func: analysis.generate_datasets
           unwrap: true

    generate_algorithms:
      ref:
        dataset: generate_datasets.datasets
      functions:
        algorithms: analysis.generate_algorithms
        
  tests:
    analysis.generate_plots:
      - grid: generate_algorithms
```

But the goal of gridtest would be to make these runs modular, meaning that we
could record metrics for each dataset and each algorithm. The example
here served to show a simple thought process to go from a single
top to bottom script to a simple grid. The [next example](../2-modular-grids) will take this a step
further to make it optimized in terms of modularity.
# Custom Decorator

If you want to run a gridtest that uses a custom metric, you can easily
do this by defining your own decorator. For example, let's say we have a function
to do some kind of text processing. It takes some number of inputs, and
returns raw text. We would then want to count the unique words in the raw text.
Let's go!

### Create your Functions

Let's first write our functions. We will write a simple function to take
some text input and parse it (`multiply_sentences`) and a decorator
to run any function that returns a string of text, and count the unique
words (`countwords`)

```python
# These are functions in my script

def multiply_sentence(sentence, count):
    return sentence * count

def countwords(func):
    """this is a simple example of a custom decorator - the idea would be that
       the function we are decorating returns some texty value, and we split
       this value by a blank space and then count the number of tokens (words).
    """
    def counter(*args, **kwargs):
        result = func(*args, **kwargs)
        words = len(set(result.split(' ')))
        print(f"@script.countwords {words} words")
        return result

    return counter
```
An important note about the decorator - it needs to be importable, meaning either
the module is already on your Python path, or it's included somewhere in the library
that you are testing. Also note that in order for gridtest to parse the result, you
need to print something to stdout that is prefixed **exactly** with the name of
the decorator defined under metrics. E.g., if we changed `script.countwords` to just
`countwords` the result wouldn't be properly parsed, because gridtest is looking
for the the first.


### Generate your Template

Let's generate a simple template that we can fill in to include a grid. We can
first preview it:

```bash
$ gridtest generate script.py

script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.countwords:
    - args:
        func: null
    script.multiply_sentence:
    - args:
        count: null
        sentence: null
```

We don't want a test for the decorator, so we will write this to file, and remove it.

```bash
$ gridtest generate script.py gridtest.yml

script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - args:
        count: null
        sentence: null
```

Next, let's better refine our arguments.

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - metrics:
      - '@script.countwords'
      args:
        count: [1, 5, 10]
        sentence:
          - "He ran for the hills."
          - "Skiddery-a rinky dinky dinky, skittery rinky doo."
          - "You are my sunshine, my only sunshine."
```

Just for your FYI - if you had wanted to have some set of arguments shared
between tests, you could have defined them as a named grid:

```yaml
grids:
  script_inputs:
    args:
      count: [1, 5, 10]
      sentence:
        - "He ran for the hills."
        - "Skiddery-a rinky dinky dinky, skittery rinky doo."
        - "You are my sunshine, my only sunshine."
```

and then instead pointed to it for your test:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - metrics:
      - '@script.countwords'
      grid: script_inputs
```

By default, a grid is generated at the test creation time. However, if you have 
a grid shared by many functions that you want to calculate once and cache,
just set the cache variable to true:


```yaml
grids:
  script_inputs:
    cache: true
    args:
      count: [1, 5, 10]
      sentence:
        - "He ran for the hills."
        - "Skiddery-a rinky dinky dinky, skittery rinky doo."
        - "You are my sunshine, my only sunshine."
```

Regardless of how we specify our grid (globally or inline) the grid says that
for each sentence under the list of sentences, we will run the function
`multiply_sentence` with counts of 1,5, and 10. This would come down to 3x3 or 9 total tests.

### Running Tests

Let's run the tests! We should see a count of words for each function.

```bash
$ gridtest test
[9/9] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.multiply_sentence.0     success                                                      
script.multiply_sentence.1     success                                                      
script.multiply_sentence.2     success                                                      
script.multiply_sentence.3     success                                                      
script.multiply_sentence.4     success                                                      
script.multiply_sentence.5     success                                                      
script.multiply_sentence.6     success                                                      
script.multiply_sentence.7     success                                                      
script.multiply_sentence.8     success                                                      

________________________________________________________________________________________________________________________
script.multiply_sentence.0     @script.countwords             5 words                       
script.multiply_sentence.1     @script.countwords             7 words                       
script.multiply_sentence.2     @script.countwords             7 words                       
script.multiply_sentence.3     @script.countwords             21 words                      
script.multiply_sentence.4     @script.countwords             31 words                      
script.multiply_sentence.5     @script.countwords             31 words                      
script.multiply_sentence.6     @script.countwords             41 words                      
script.multiply_sentence.7     @script.countwords             61 words                      
script.multiply_sentence.8     @script.countwords             61 words                      

9/9 tests passed
```

Awesome! We've run a grid of tests over a grid of inputs, and measured some metric with
our custom decorators. See the [gridtest.yml](gridtest.yml)
for the full test recipe.

### Save to File

If we want to save the complete results, we can add `--save` with a filename:

```bash
$ gridtest test --save results.json
```

You can see the example [results.json](results.json) in this folder.
# Class Testing

The example here [car.py](car.py) provides a script file with a simple class definition.
We can first generate a testing config for preview - if we use gridtest generate without
an output file, it will print to the screen:


```bash
$ gridtest generate car.py
Extracting Car from car
Extracting honk from car
Extracting lights from car
```
```yaml
car:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/class/car.py
  tests:
    car.Car:
    - args:
        color: red
        lights: false
        wheels: 4
    car.Car.honk:
    - args:
        self: null
    car.Car.lights:
    - args:
        self: null
```

We can then write to file:

```bash
$ gridtest generate car.py gridtest.yml
```

Notice that we would want to be able to define an instance of a class to be
used for the subsequent testing functions of the class. To do this, take the 
section you want to use, for example this one:


```yaml
  car.Car:
  - args:
      color: red
      lights: false
      wheels: 4
```

and add an "instance" key to it.


```yaml
  car.Car:
  - instance: thisone
    args:
      color: red
      lights: false
      wheels: 4
```


We now want to reference "thisone" as the instance to
use. Just update the "self" variable in each of the class test cases.

```yaml
car:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/class/car.py
  tests:
    car.Car:
    - args:
        color: notacolor
        lights: false
        wheels: 4
      raises: ColorException
    - instance: thisone
      args:
        color: red
        lights: false
        wheels: 4
      isinstance: car.Car 
    car.Car.honk:
    - args:
        self: "{{ instance.thisone }}"
    car.Car.switch_lights:
    - args:
        self: "{{ instance.thisone }}"
```

And then the instance of the Car named as instance "this one" (the second block)
will be used for those tests. 

### Testing

Now, we can run tests! Since we've named the testing file `gridtest.yml` we can
just run:

```bash
$ gridtest test
[4/4] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
car.Car.0                      success                        raises ColorException         
car.Car.1                      success                        isinstance Car                
car.Car.honk.0                 success                                                      
car.Car.switch_lights.0        success                                                      

4/4 tests passed
```

This is a very basic usage for a class, and we 
expect more complex cases to be written up when they are determined.
# Results Interface

This is an example that starts with the [custom-decorator](../custom-decorator)
example and shows how to generate an interactive interface with results. We've
added an additional function and metric decorator to produce a more interesting
result:

```bash
$ gridtest test
[12/12] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.multiply_sentence.0     success                                                      
script.multiply_sentence.1     success                                                      
script.multiply_sentence.2     success                                                      
script.multiply_sentence.3     success                                                      
script.multiply_sentence.4     success                                                      
script.multiply_sentence.5     success                                                      
script.multiply_sentence.6     success                                                      
script.multiply_sentence.7     success                                                      
script.multiply_sentence.8     success                                                      
script.unique_sentence.0       success                                                      
script.unique_sentence.1       success                                                      
script.unique_sentence.2       success                                                      

________________________________________________________________________________________________________________________
script.multiply_sentence.0     @script.countwords             5 words                       
script.multiply_sentence.0     @script.countletters           17 letters                    
script.multiply_sentence.1     @script.countwords             7 words                       
script.multiply_sentence.1     @script.countletters           43 letters                    
script.multiply_sentence.2     @script.countwords             7 words                       
script.multiply_sentence.2     @script.countletters           32 letters                    
script.multiply_sentence.3     @script.countwords             21 words                      
script.multiply_sentence.3     @script.countletters           85 letters                    
script.multiply_sentence.4     @script.countwords             31 words                      
script.multiply_sentence.4     @script.countletters           215 letters                   
script.multiply_sentence.5     @script.countwords             31 words                      
script.multiply_sentence.5     @script.countletters           160 letters                   
script.multiply_sentence.6     @script.countwords             41 words                      
script.multiply_sentence.6     @script.countletters           170 letters                   
script.multiply_sentence.7     @script.countwords             61 words                      
script.multiply_sentence.7     @script.countletters           430 letters                   
script.multiply_sentence.8     @script.countwords             61 words                      
script.multiply_sentence.8     @script.countletters           320 letters                   
script.unique_sentence.0       @script.countwords             5 words                       
script.unique_sentence.0       @script.countletters           17 letters                    
script.unique_sentence.1       @script.countwords             6 words                       
script.unique_sentence.1       @script.countletters           38 letters                    
script.unique_sentence.2       @script.countwords             6 words                       
script.unique_sentence.2       @script.countletters           30 letters                    

12/12 tests passed
```

## Run Tests

Generating the interface is simple! Remember that if we wanted to save the raw
results, we could add `--save` with a filename:

```bash
$ gridtest test --save results.json
```

You can do the same for the interface, but instead provide the --save-web
flag. If you don't provide an argument to this flag, it will generate
a folder for you in `/tmp` with files that can be added to a web server.
If you do provide a folder path, the same folder will be renamed
to be there (so make sure it doesn't exist). The following command
will save pages to [web](web) in the present working directory.

```bash
$ gridtest test --save-web web/
```

You should then be able to put these static files on GitHub pages, or just cd
into the folder and run a webserver:

```bash
python -m http.server 9999
```

And see the content at [http://localhost:9999](http://localhost:9999). This
is a fairly simple results template that lets you select functions in the left
columns, and then see specific tests in the right table

![img/gridtest.png](img/gridtest.png)

and mouse over the test name to see the output, error, and metrics recorded.

![img/detail.png](img/detail.png)
# Package Testing

You don't necessarily need to write tests just for local files or modules!
Gridtest also works to generate tests for modules installed with your Python.
This provides a quick example for that. First, generate a testing file:

```bash
gridtest generate requests gridtest.yml --skip-classes
```

Note that we just want to test the simple requests functions, so we are 
skipping classes. If I wanted to include private functions:

```bash
$ gridtest generate --include-private --skip-classes requests gridtest.yml
```

In the case of requests, I don't really want to test more of the complex
functionlity, but just the functions that I might need to use. Thus,  
I can then open the file and just delete those chunks. I'd then
update the file to customize it with my tests, meaning that I change this:

```yaml
requests:
  filename: /home/vanessa/anaconda3/lib/python3.7/site-packages/requests/__init__.py
  tests:
    _internal_utils.to_native_string:
    - args:
        encoding: null
        string: ascii
    _internal_utils.unicode_is_ascii:
    - args:
        u_string: null
    adapters._basic_auth_str:
    - args:
        password: null
        username: null
    adapters.extract_cookies_to_jar:
    - args:
        jar: null
        request: null
        response: null
    adapters.extract_zipped_paths:
    - args:
        path: null
    adapters.get_auth_from_url:
    - args:
        url: null
    adapters.get_encoding_from_headers:
    - args:
        headers: null
    adapters.prepend_scheme_if_needed:
    - args:
        new_scheme: null
        url: null
    adapters.select_proxy:
    - args:
        proxies: null
        url: null
    adapters.urldefragauth:
    - args:
        url: null
    auth.parse_dict_header:
    - args:
        value: null
    cookies._copy_cookie_jar:
    - args:
        jar: null
    cookies.cookiejar_from_dict:
    - args:
        cookie_dict: null
        cookiejar: true
        overwrite: null
    cookies.create_cookie:
    - args:
        name: null
        value: null
    cookies.get_cookie_header:
    - args:
        jar: null
        request: null
    cookies.merge_cookies:
    - args:
        cookiejar: null
        cookies: null
    cookies.morsel_to_cookie:
    - args:
        morsel: null
    cookies.remove_cookie_by_name:
    - args:
        cookiejar: null
        domain: null
        name: null
        path: null
    hooks.default_hooks:
    - args: {}
    hooks.dispatch_hook:
    - args:
        hook_data: null
        hooks: null
        key: null
    models.check_header_validity:
    - args:
        header: null
    models.guess_filename:
    - args:
        obj: null
    models.guess_json_utf:
    - args:
        data: null
    models.iter_slices:
    - args:
        slice_length: null
        string: null
    models.parse_header_links:
    - args:
        value: null
    models.requote_uri:
    - args:
        uri: null
    models.stream_decode_response_unicode:
    - args:
        iterator: null
        r: null
    models.super_len:
    - args:
        o: null
    models.to_key_val_list:
    - args:
        value: null
    requests._check_cryptography:
    - args:
        cryptography_version: null
    requests.check_compatibility:
    - args:
        chardet_version: null
        urllib3_version: null
    requests.delete:
    - args:
        url: null
    requests.get:
    - args:
        params: null
        url: null
    requests.head:
    - args:
        url: null
    requests.options:
    - args:
        url: null
    requests.patch:
    - args:
        data: null
        url: null
    requests.post:
    - args:
        data: null
        json: null
        url: null
    requests.put:
    - args:
        data: null
        url: null
    requests.request:
    - args:
        method: null
        url: null
    requests.session:
    - args: {}
    sessions.default_headers:
    - args: {}
    sessions.get_environ_proxies:
    - args:
        no_proxy: null
        url: null
    sessions.get_netrc_auth:
    - args:
        raise_errors: null
        url: false
    sessions.merge_hooks:
    - args:
        dict_class: null
        request_hooks: &id001 !!python/name:collections.OrderedDict ''
        session_hooks: null
    sessions.merge_setting:
    - args:
        dict_class: null
        request_setting: *id001
        session_setting: null
    sessions.rewind_body:
    - args:
        prepared_request: null
    sessions.should_bypass_proxies:
    - args:
        no_proxy: null
        url: null
    status_codes._init:
    - args: {}
    utils._parse_content_type_header:
    - args:
        header: null
    utils.add_dict_to_cookiejar:
    - args:
        cj: null
        cookie_dict: null
    utils.address_in_network:
    - args:
        ip: null
        net: null
    utils.default_user_agent:
    - args:
        name: python-requests
    utils.dict_from_cookiejar:
    - args:
        cj: null
    utils.dict_to_sequence:
    - args:
        d: null
    utils.dotted_netmask:
    - args:
        mask: null
    utils.from_key_val_list:
    - args:
        value: null
    utils.get_encodings_from_content:
    - args:
        content: null
    utils.get_unicode_from_response:
    - args:
        r: null
    utils.is_ipv4_address:
    - args:
        string_ip: null
    utils.is_valid_cidr:
    - args:
        string_network: null
    utils.parse_list_header:
    - args:
        value: null
    utils.unquote_header_value:
    - args:
        is_filename: null
        value: false
    utils.unquote_unreserved:
    - args:
        uri: null
```

to this:

```yaml
requests:
  filename: /home/vanessa/anaconda3/lib/python3.7/site-packages/requests/__init__.py
  tests:
    requests.api.get:
    - args:
        params: null
        url: https://google.com
      isinstance: Response
    requests.api.head:
    - args:
        url: https://google.com
      istrue: "self.result.status_code == 301"
    requests.api.options:
    - args:
        url: https://google.com
```

This recipe also shows a good example of how to check for an instance type.
The string "Response" should be given if I check the `type(result).__name__`
for the result. Also notice the `istrue` statement isn't targeting a `{{ result }}`
template that can be converted to a string and evaluated, but rather the GridTest
instance directly (self.result) and more specifically, the status_code attribute.

## Running Tests

And then run your tests:

```bash
$ gridtest test
[3/3] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
requests.api.get.0             success                        isinstance Response           
requests.api.head.0            success                        istrue self.result.status_code == 301
requests.api.options.0         success                                                      

3/3 tests passed
```
# Read and Write Example

This simple example will show you how to use the functions `{% tmp_path %}`
and `{% tmp_dir %}` to create and then reference a temporary file path or
directory, respectively.

## Install

After you've installed gridtest:

```bash
git clone git@github.com:vsoch/gridtest
cd gridtest
pip install -e .
```
or

```bash
pip install gridtest
```

## Generate

We generated the initial template for the script [temp.py](temp.py)
included here like this:

```bash
$ gridtest generate temp.py gridtest.yml
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the [gridtest.yml](gridtest.yml) will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case.

## Customize

You can then open the file in a text editor, and add arguments to each.
If your library uses typing, the typing will be checked at testing time,
and it's not specified here. You'll generally want to fill in args for
each testing condition (or leave null for None). For more detail about
templating, see [the templating documentation](https://vsoch.github.io/gridtest/getting-started/templating). 
We ultimately updated the template to include the following:

```yaml
temp:
  filename: temp.py
  tests:
    temp.create_directory:
    - args:
        dirname: "{% tmp_dir %}"
      returns: "{{ args.dirname }}"
    temp.write_file:
    - args:
        filename: "{% tmp_path %}"
      exists: "{{ args.filename }}"
```

The above recipe says that we want to test the function `create_directory`
in [temp.py](temp.py), and we want gridtest to generate a temporary directory
(the variable `{% tmp_dir %}` and then test that whatever name is generated
(`{{ args.dirname }}`) is returned by the function. For the function write_file,
the input "filename" will have a randomly generated temporary file created
for it with `{% tmp_path %}` that will want to ensure exists. Both will be 
cleaned up at the completion of the test.

## Test

We can run the tests as follows:

```bash
$ gridtest test gridtest.yml
[2/2] |===================================| 100.0%
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
temp.create_directory.0        success                        returns /tmp/gridtest-dir.ok6x5kiy
temp.write_file.0              success                        exists /tmp/gridtest-file-flp_cmi4

2/2 tests passed
```

And the directory mentioned, `/tmp/gridtest-dir.j0aio0l6` will be cleaned up
upon completion. If we don't want to clean it up, we can add `--no-cleanup`:

```bash
$ gridtest test gridtest.yml --no-cleanup
[2/2] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
temp.create_directory.0        success                        returns /tmp/gridtest-dir.n1devo3f
temp.write_file.0              success                        exists /tmp/gridtest-file-iqff2y__

2/2 tests passed
``` 

And then the directory generated would still exist after the run:

```bash
$ ls -l /tmp/gridtest-dir.e1c4gbr8/
total 0
```

Also, since gridtest.yml is the default, you don't need to specify it to
find the file in the present working directory:

```bash
gridtest gridtest.yml
```

The same is true for the testing file.
# Sample Grid

Gridtest can easily be used to generate random samplies for some number of inputs,
where each input is returned via a function as a list of options to select from.

## Write your Functions

Let's start by writing a set of functions. Each of these will return a list
of attributes that we might want to parameterize. For our first example,
we will generate cohorts with all possible combinations. Let's create functions
to return colors, ages, shapes, and animals.

```python
# Example functions to generate lists to parameterize

import colorsys
import random

def generate_rgb_color():
    """return a randomly selected color
    """
    N = int(random.choice(range(0,100)))
    HSV_tuple = (N*1.0/N, 0.5, 0.5)
    return colorsys.hsv_to_rgb(*HSV_tuple)

def generate_shape():
    return random.choice(["square", "triangle", "circle", "ellipsis", "rectangle", "octagon"])
        

def generate_age():
    return random.choice(range(0,100))

def generate_animal():
    return random.choice(["dog", "cat", "bird", "cow", "chicken"])
```

It's fairly simple - each one spits out a random value. You could also
imagine loading data from some other source.

## Generate the Grids

We next want to generate grids for each! This is fairly simple to do too -
each grid is named based on what it selects, and uses the appropriate function
to be run:

```yaml
script:
  grids:

    # Each grid below generates one randomly selected value
    select_color:
      functions:
        color: script.generate_rgb_color

    select_shape:
      functions:
        shape: script.generate_shape

    select_animal:
      functions:
        animal: script.generate_animal

    select_age:
      functions:
        age: script.generate_age
```
We can also verify that each one generates what we would expect with `gridview`:

```bash
$ gridtest gridview cohort-grids.yml select_color
{'color': (0.5, 0.25, 0.25)}

$ gridtest gridview cohort-grids.yml select_animal
{'animal': 'bird'}

$ gridtest gridview cohort-grids.yml select_age
{'age': 47}

$ gridtest gridview cohort-grids.yml select_shape
{'shape': 'triangle'}
```

But that's not very interesting! Let's instead create a full set of data
for a single sample of a cohort:

```yaml
    generate_cohort:
      functions:
        age: script.generate_age
        animal: script.generate_animal
        shape: script.generate_shape
        color: script.generate_rgb_color
```

and generate it:

```yaml
$ gridtest gridview cohort-grids.yml generate_cohort
{'age': 54, 'animal': 'dog', 'shape': 'ellipsis', 'color': (0.5, 0.25, 0.25)}
```

We could run it again to get a different result, but why not just add a
count for the number of samples that we need?

```yaml
    generate_cohort:
      count: 10
      functions:
        age: script.generate_age
        animal: script.generate_animal
        shape: script.generate_shape
        color: script.generate_rgb_color
```

```bash
{'age': 71, 'animal': 'dog', 'shape': 'triangle', 'color': (0.5, 0.25, 0.25)}
{'age': 18, 'animal': 'bird', 'shape': 'circle', 'color': (0.5, 0.25, 0.25)}
{'age': 17, 'animal': 'bird', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 24, 'animal': 'chicken', 'shape': 'rectangle', 'color': (0.5, 0.25, 0.25)}
{'age': 14, 'animal': 'bird', 'shape': 'octagon', 'color': (0.5, 0.25, 0.25)}
{'age': 82, 'animal': 'cow', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 4, 'animal': 'cat', 'shape': 'rectangle', 'color': (0.5, 0.25, 0.25)}
{'age': 55, 'animal': 'cat', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 48, 'animal': 'bird', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 84, 'animal': 'chicken', 'shape': 'triangle', 'color': (0.5, 0.25, 0.25)}
```

We can see that we have ten results:

```bash
$ gridtest gridview cohort-grids.yml generate_cohort --count
10 argument sets produced.
```

Great! Let's save this to file. We don't need gridtest anymore, we can just keep
our parameters for use.

```bash
$ gridtest gridview cohort-grids.yml generate_cohort --export samples.json
```

The results are saved to json.

```json
[
    {
        "age": 88,
        "animal": "bird",
        "shape": "ellipsis",
        "color": [
            0.5,
            0.25,
            0.25
        ]
    },
    {
        "age": 94,
        "animal": "cat",
        "shape": "rectangle",
        "color": [
            0.5,
            0.25,
            0.25
        ]
    },
...
]
```

Take a look at [cohort-grids.yml](cohort-grids.yml) for the grids.

## Combination Grids

The above works, but it calls the same functions many times. We would do much
better to, for some number of samples that we need, run each function once with
that number (returning a list of lists) and then unwrap it into a grid. For this
first example, we want all possible combinations of the parameters. Let's again
write our functions, and this time, we return a list:

```python
def generate_rgb_color(N=10):
    return generate_rgb_color() * N

def generate_shapes(N):
    return generate_shape() * N
        
def generate_ages(N=10):
    return generate_age() * N

def generate_animals(N=10):
    return generate_animal() * N
```

We can now use "unwrap" to make sure that each argument is represented on its own.
This way, we can parameterize them each over a larger grid. We can check to see this
in the output:

```bash
$ gridtest gridview grids.yml select_shapes
{'shapes': ['octagon']}
{'shapes': ['circle']}
{'shapes': ['rectangle']}
{'shapes': ['square']}
{'shapes': ['square']}
{'shapes': ['rectangle']}
{'shapes': ['ellipsis']}
{'shapes': ['octagon']}
{'shapes': ['ellipsis']}
{'shapes': ['rectangle']}

$ gridtest gridview grids.yml select_colors
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}

$ gridtest gridview grids.yml select_ages
{'ages': [59]}
{'ages': [44]}
{'ages': [52]}
{'ages': [78]}
{'ages': [78]}
{'ages': [12]}
{'ages': [82]}
{'ages': [40]}
{'ages': [64]}
{'ages': [66]}

$ gridtest gridview grids.yml select_animals
{'animals': ['chicken']}
{'animals': ['cow']}
{'animals': ['cat']}
{'animals': ['chicken']}
{'animals': ['cow']}
{'animals': ['cow']}
{'animals': ['cat']}
{'animals': ['dog']}
{'animals': ['chicken']}
{'animals': ['cow']}
```

Finally, let's make our grid that will generate many combinations of each. Just
a heads up ... it's every possible combination, so we will have 10 x 10 x 10 x 10...
10,000!

```bash
$ gridtest gridview grids.yml generate_cohort --count
10000 argument sets produced.
```
```bash
$ gridtest gridview grids.yml generate_cohort --count
...
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
```

And you could save this to file again:

```bash
$ gridtest gridview grids.yml generate_cohort --save combinations.json
```

Take a look at [grids.yml](grids.yml) for the grids. This grid could then be
input into some testing function for further use.

## Overview

Hopefully you can see the two use cases here - the first is to generate a
specific set of samples, given some input functions that generate values.
The second is to produce all possible combinations, and we take
advantage of unwrapping lists. For this dummy example, we are actually
calling the same function many times for the generation of lists of lists.
However, you could imagine having a single computationally intensive function
or maybe an API call that generates a list that you want to use the items
across another grid.
# Grids

GridTest isn't just for testing! In fact, you can write a file of grid
specifications that can be loaded and used for parameterization, or your own
custom functions. Let's take a look at an example [grids.yml](grids.yml)
for how to do that, and [grids-with-function.yml](grids-with-function.yml)
for plugging those parameters into functions.

## Loading via a GridRunner

Let's say that we have a grids.yml file with grids defined. We can load with
gridtest easily via the GridRunner:

```python
from gridtest.main.test import GridRunner
runner = GridRunner('grids.yml')
# [gridtest|grids.yml]
```

We can easily generate the grids as follows:

```python
runner.get_grids()
{'generate_empty': [grid|generate_empty],
 'generate_matrix': [grid|generate_matrix],
 'generate_lists_matrix': [grid|generate_lists_matrix],
 'generate_by_min_max': [grid|generate_by_min_max],
 'generate_by_min_max_twovars': [grid|generate_by_min_max_twovars]}
```

## Loading as a Grid

A Grid is a first class citizen, so you can also generate without
the GridRunner, either via a json object or yaml you've loaded:

```python
from gridtest.main.grids import Grid

grid = Grid(name="mygrid", params={"args": {"one": [1, 11, 111], "two": 2}})
[grid|mygrid]
```

### Cached Loading

By default, grids are generated on demand, meaning that the argument sets (optionally
derived by functions you have provided under params) are generated at this point in time.
For example:

```
for argset in grid:
    print(argset)

{'one': 1, 'two': 2}
{'one': 11, 'two': 2}
{'one': 111, 'two': 2}
```

This also means that, by default, your grid.argsets will be empty when you
instantiate it.

```python
grid = Grid(name="mygrid", params={"args": {"one": [1, 11, 111], "two": 2}})
grid.argsets
[]
```

This is done intentionally because it can sometimes take a lot of memory to store
a long list of argument sets. However, if you want to generate the argument sets
on init, just set "cache" to True in your params:

```python
grid = Grid(name="mygrid", params={"cache": True, "args": {"one": [1, 11, 111], "two": 2}})
```

You should then be able to see the argsets ready for your use!

```python
grid.argsets
[{'one': 1, 'two': 2}, {'one': 11, 'two': 2}, {'one': 111, 'two': 2}]
```

### Grid Functions

What if you want to generate a long list of items for an argument, and it would
be crazy to type them out? This is what the functions specification for a grid
is for. You would provide your functions, each associated with an argument, under
the "functions" section. As an example, let's say we want to run random.choice,
and select from a list defined under the argument "choices":

```python
import random
grid = Grid(name="mygrid", params={"functions": {"pid": random.choice}, "args": {"seq": [[1,2,3,4,5,6,7]]}})
```

The above would look like this in yaml:

```yaml
grids:
  mygrid:
    args:
      seq: [[1,2,3,4,5,6,7]]
    functions:
      pid: random.choice
```

And we would generate our argument sets:

```python
list(grid)
[{'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 2}]
```

This also points out another important point - if you were to provide a list for an
argument, it would be parameterized. If you want the entire list treated as one variable,
then define it within another list as shown above. The above grid says "generate the argument
pid by using random.choice to select from seq." How does the grid know to match seq to random.choice?
It's a known keyword argument! Another insight here is that although we define a sequence argument, we don't actually
need it, we really only care about the result (pid). But let's say we want to generate
10 x a pid, how do we do that? We add a count to the params:

```python
grid = Grid(name="mygrid", params={"count": 10, "functions": {"pid": random.choice}, "args": {"seq": [[1,2,3,4,5,6,7]]}})
list(grid)

[{'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 5},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 3},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 7},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 5},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 7},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 4},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 3}]
```

These values are being generated by an iterator and not saved to any single list in memory,
so although it seems very rundant to see the same list more than once, it shouldn't
serve significant issues with respect to memory.


## Viewing on the Command Line

You can also preview the grids generated from the command line. First you might
want to list all grids defined in a file:

```bash
$ gridtest gridview grids.yml
generate_empty
generate_matrix
generate_lists_matrix
generate_by_min_max
generate_by_min_max_twovars
```

and then inspect a specific named grid:

```bash
$ gridtest gridview grids.yml generate_empty
{}
{}
{}
{}
{}
{}
{}
{}
{}
{}
```

If you want to export your entire grid to json, possibly to read in as-is
later and use, add --export:

```bash
$  gridtest gridview grids.yml generate_by_min_max --export exported.json 
```
```bash
cat exported.json
[
    {
        "x": 0.0
    },
    {
        "x": 2.0
    },
    {
        "x": 4.0
    },
    {
        "x": 6.0
    },
    {
        "x": 8.0
    }
]
```

That's not a really interesting grid, but you get the gist.
Each of the grids listed above will be explained below in more detail.

## Writing Grids

Let's start with the most basic of grids, which are those that don't import any
special functions. 

### Header

The header should have a named section for the grids that you want to define
(e.g. `mygrids`) and then a `grids` index, under which we will write each
of our named grids. Here is the start:

```
mygrids:
  grids:
  ...
```

Then each grid is added as a named section below that:

```
mygrids:
  grids:

    # A grid that will generate 10 empty set of arguments
    generate_empty:
      count: 10
```

The full content of this file will be written to the [grids.yml](grids.yml)
example here. Each example grid is discussed below. For each example,
you can preview the grid in the terminal with `gridtest gridview` or
obtain the grid by instantiating the `GridRunner` as shown above.

### Empty

It could be that you need to generate empty lists of arguments for a function,
in which case the `count` variable will be useful to you. Here is
how to specify a grid that will generate 10 empty set of arguments

```yaml
generate_empty:
  count: 10
```

The result comes out to be:

```bash
$ gridtest gridview grids.yml generate_empty
{}
{}
{}
{}
{}
{}
{}
{}
{}
{}
```

### Parameterize Variables

Let's say that we have two variables, x and y, and we want to generate a grid
of all possible combinations for a listing of each. That would look like this:

```yaml
    generate_matrix:
      args:
        x: [1, 2, 3] 
        y: [1, 2, 3]
```

And the resulting grid will have 3x3 or 9 total combinations of x and y:

```bash
$ gridtest gridview grids.yml generate_matrix
{'x': 1, 'y': 1}
{'x': 1, 'y': 2}
{'x': 1, 'y': 3}
{'x': 2, 'y': 1}
{'x': 2, 'y': 2}
{'x': 2, 'y': 3}
{'x': 3, 'y': 1}
{'x': 3, 'y': 2}
{'x': 3, 'y': 3}
```

### Parameterize Lists

If you want to do similar but instead have a list of values be paramaterized, just 
specify a list of lists instead.

```yaml
    generate_lists_matrix:
      args:
        x: [[1, 2, 3], [4, 5, 6]] 
        y: [[1, 2, 3], [4, 5, 6]] 
```

The result will have 2x2 or 4 entries:

```bash
$ gridtest gridview grids.yml generate_lists_matrix
{'x': [1, 2, 3], 'y': [1, 2, 3]}
{'x': [1, 2, 3], 'y': [4, 5, 6]}
{'x': [4, 5, 6], 'y': [1, 2, 3]}
{'x': [4, 5, 6], 'y': [4, 5, 6]}
```

Here is an easier way to check the count:

```bash
$ gridtest gridview grids.yml generate_lists_matrix --count
4 argument sets produced.
```

### Range of Values

For most use cases, you'll want to generate a list of values over a range. You
can do that with **min** and **max** and (optionally) **by** that defaults to 1.

```yaml
    generate_by_min_max:
      args:
        x:
          min: 0
          max: 10
          by: 2
```

```bash
$ gridtest gridview grids.yml generate_by_min_max
{'x': 0.0}
{'x': 2.0}
{'x': 4.0}
{'x': 6.0}
{'x': 8.0
```

Here is an example with two variables:

```yaml
    generate_by_min_max_twovars:
      args:
        x:
          min: 0
          max: 10
          by: 2
        y:
          min: 10
          max: 20
          by: 2
```

```bash
$ gridtest gridview grids.yml generate_by_min_max_twovars 
{'x': 0.0, 'y': 10.0}
{'x': 0.0, 'y': 12.0}
{'x': 0.0, 'y': 14.0}
{'x': 0.0, 'y': 16.0}
{'x': 0.0, 'y': 18.0}
{'x': 2.0, 'y': 10.0}
{'x': 2.0, 'y': 12.0}
{'x': 2.0, 'y': 14.0}
{'x': 2.0, 'y': 16.0}
{'x': 2.0, 'y': 18.0}
{'x': 4.0, 'y': 10.0}
{'x': 4.0, 'y': 12.0}
{'x': 4.0, 'y': 14.0}
{'x': 4.0, 'y': 16.0}
{'x': 4.0, 'y': 18.0}
{'x': 6.0, 'y': 10.0}
{'x': 6.0, 'y': 12.0}
{'x': 6.0, 'y': 14.0}
{'x': 6.0, 'y': 16.0}
{'x': 6.0, 'y': 18.0}
{'x': 8.0, 'y': 10.0}
{'x': 8.0, 'y': 12.0}
{'x': 8.0, 'y': 14.0}
{'x': 8.0, 'y': 16.0}
{'x': 8.0, 'y': 18.0}
```

Logically there are the previous number of tests, but squared.

```bash
$ gridtest gridview grids.yml generate_by_min_max_twovars --count
25 argument sets produced.
```

## Grids with Functions

## What does it mean to use a function?

You can map functions to parametrs in a grid, meaning that the values for
the parameters will be generated by the function. For example,
let's use the function `random.choice` to dynamically generic a grid of parameters.
The grid below will call `random.choice` ten times (count is set to 10) 
across the sequence of values `[1, 2, 3]`, which is an input key word argument
to random choice.

```yaml
...
  grids:
    random_choice:
      count: 10
      functions: 
        pid: random.choice
      args:
         seq: [[1, 2, 3]]
```

Notice that the sequence argument input is a list of lists, and this is because we want the entire list to be treated as an argument. If we run this from it's respective file, we get a result with 10
argument sets:

```bash
$ gridtest gridview grids-with-function.yml random_choice
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 1}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 1}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 3}
```


## Grids with Custom Functions

If you are using a script for any of your grids that isn't a system installed
module, then (akin to a standard gridtest) it needs to be included under a section 
header that is named by the relevant module, and that includes the filename to import.
For example, the file [script.py](script.py) in the present working directory has
a function, `get_pokemon_id` that I want to use. Here is how I'd write the recipe:

script:
  filename: script.py 
  grids:

    # A grid that will generate 10 random values using a custom function
    generate_pids:
      functions: 
        pid: script.get_pokemon_id
      count: 10
```

Now let's run it!

```bash
$ gridtest gridview grids-with-function.yml generate_pids
{'pid': '125'}
{'pid': '780'}
{'pid': '508'}
{'pid': '566'}
{'pid': '803'}
{'pid': '513'}
{'pid': '854'}
{'pid': '405'}
{'pid': '639'}
{'pid': '353'}
```

Akin to the previous example, since we've provided a function to pass our grid
arguments into, the results are returned. This is how gridtest can
use a grid specified under a test to generate a list of values for an argument.

## Grids with Unwrapped Functions

Let's say that we have some function that produces a list of lists, and we want to
use that as a list to be parameterized for another grid. First, here is our
function:

```python
```

And here is how we might define it in the grid:
```yaml
    # Generate a list of lists, intended to be unwrapped and used for another grid
    unwrapped_grid:
      functions:
        numbers: 
          func: script.generate_numbers
          unwrap: true
```

The "unwrap" serves to take the list, and unwrap it so that each result can
be used separately. If we generate the grid, the parameterization is done with
any additional arguments. Since we don't have any, we get a list of 10 inputs, 
each of length 10

```bash
$ gridtest gridview grids-with-function.yml unwrapped_grid 
{'numbers': [82, 82, 82, 82, 82, 82, 82, 82, 82, 82]}
{'numbers': [52, 52, 52, 52, 52, 52, 52, 52, 52, 52]}
{'numbers': [41, 41, 41, 41, 41, 41, 41, 41, 41, 41]}
{'numbers': [25, 25, 25, 25, 25, 25, 25, 25, 25, 25]}
{'numbers': [43, 43, 43, 43, 43, 43, 43, 43, 43, 43]}
{'numbers': [64, 64, 64, 64, 64, 64, 64, 64, 64, 64]}
{'numbers': [43, 43, 43, 43, 43, 43, 43, 43, 43, 43]}
{'numbers': [39, 39, 39, 39, 39, 39, 39, 39, 39, 39]}
{'numbers': [27, 27, 27, 27, 27, 27, 27, 27, 27, 27]}
{'numbers': [64, 64, 64, 64, 64, 64, 64, 64, 64, 64]}
```

We can also ask to just look at the numbers parameter, as it was defined
before we produced the parameter matrix above (it's a list of lists)

```bash
$ gridtest gridview grids-with-function.yml unwrapped_grid --arg numbers --count
Variable numbers has length 10.

$ gridtest gridview grids-with-function.yml unwrapped_grid --arg numbers
[[32, 32, 32, 32, 32, 32, 32, 32, 32, 32], [70, 70, 70, 70, 70, 70, 70, 70, 70, 70], [75, 75, 75, 75, 75, 75, 75, 75, 75, 75], [86, 86, 86, 86, 86, 86, 86, 86, 86, 86], [86, 86, 86, 86, 86, 86, 86, 86, 86, 86], [50, 50, 50, 50, 50, 50, 50, 50, 50, 50], [92, 92, 92, 92, 92, 92, 92, 92, 92, 92], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [80, 80, 80, 80, 80, 80, 80, 80, 80, 80], [35, 35, 35, 35, 35, 35, 35, 35, 35, 35]]
```

We can also add another parameter (e.g., defining two for length doubles the results)

```yaml
    # Generate a list of lists, intended to be unwrapped and used for another grid
    unwrapped_grid:
      args:
        length: [10, 20]
      functions:
        numbers: 
          func: script.generate_numbers
          unwrap: true
```

Regardless, we could then reference this variable with `ref` for another grid.  For
example, here we want to 

```yaml
    sum_unwrapped_numbers:
      ref:
        numbers: unwrapped_grid.numbers
      functions: 
        total: script.dosum
```

And then we generate a grid with the result.

```bash
$ gridtest gridview grids-with-function.yml sum_unwrapped_numbers
{'numbers': [90, 90, 90, 90, 90, 90, 90, 90, 90, 90], 'total': 900}
{'numbers': [72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 'total': 720}
{'numbers': [88, 88, 88, 88, 88, 88, 88, 88, 88, 88], 'total': 880}
{'numbers': [21, 21, 21, 21, 21, 21, 21, 21, 21, 21], 'total': 210}
{'numbers': [37, 37, 37, 37, 37, 37, 37, 37, 37, 37], 'total': 370}
{'numbers': [96, 96, 96, 96, 96, 96, 96, 96, 96, 96], 'total': 960}
{'numbers': [40, 40, 40, 40, 40, 40, 40, 40, 40, 40], 'total': 400}
{'numbers': [37, 37, 37, 37, 37, 37, 37, 37, 37, 37], 'total': 370}
{'numbers': [26, 26, 26, 26, 26, 26, 26, 26, 26, 26], 'total': 260}
{'numbers': [2, 2, 2, 2, 2, 2, 2, 2, 2, 2], 'total': 20}
```
# Basic Example

This basic example will show how to generate and run a grid test.

## Install

After you've installed gridtest:

```bash
git clone git@github.com:vsoch/gridtest
cd gridtest
pip install -e .
```

or

```bash
pip install gridtest
```

## Generate

then you can cd into this folder, and test generating a gridtest file for the
[script.py](script.py) included here:

```bash
$ gridtest generate script.py gridtest.yml
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the [gridtest.yml](gridtest.yml) will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case.

## Customize

You can then open the file in a text editor, and add arguments to each.
If your library uses typing, the typing will be checked at testing time,
and it's not specified here. You'll generally want to fill in args for
each testing condition (or leave null for None). For example, we might want to 
change:

```yaml
  script.add:
    args:
    - one: null
    - two: null
```

to instead be:

```yaml
  script.add:
    args:
    - one: 1
    - two: 2
```

To test adding 1+2. 

### Input Types

Gridtest can support multiple different kinds of input types that correspond
with common use cases. For example, many tests want to use a temporary
file as input. We might do something like:


```yaml
  script.write:
    args:
    - name: "dinosaur"
    - outfile: {% tmp_file %}
    returns: {{ args.outfile }}
```

In the above, the double `{{}}` refers to a variable, in this case which is
an argument name. The `{% %}` refers to a function that is known natively
to grid test. This syntax is based on [jinja2](https://jinja.palletsprojects.com/en/2.11.x/).
You could also check that the file exists (and it might not be returned).

```yaml
  script.write:
    args:
    - name: vanessa
    - outfile: {% tmp_file %}
    exists: {{ args.outfile }}
```

### Return Types

For basic testing, there are typically a few obvious cases we want to test for:

 - returns: meaning that a function returns a specific value
 - raises: meaning that the function raises an error
 - exists: meaning that some output file is deemed to exist.

If you see another simple testing case that you want added, please 
[open an issue](https://github.com/vsoch/gridtest/issues). Highly complex testing needs
probably should use a more substantial testing library. Let's look through how each of
these examples might be used for our add function.

**returns**

If we want to ensure that the correct value is returned, we would do:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    returns: 3
```

**raises**

If we wanted to ensure that an exception was raised, we would do:

```yaml
  script.add:
  - args:
      one: 1
      two: null
    raises: TypeError
```

**istrue**

istrue is used when we want to check if something is True.
You usually would want to refer to an input or output variable:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    istrue: isinstance({% returns %}, int)
```

**isfalse**

or you might want the opposite, isfalse:


```yaml
  script.add:
  - args:
      one: 1
      two: 2
    isfalse: not isinstance({% returns %}, float)
```

**equals**

or you might want to evaluate a statement. In the example below, we want to 
make sure an attribute of the value returned is equal to 200.

```yaml
  requests.get:
  - args:
      url: https://google.com
      data: null
      json: null
    eval: {% returns %}.status_code == 200
```

This is different from providing an explicit value.

**success**

You might just want to run something, and be sure that the success status is False.
For example, if you give the wrong type as input to a function, it will by default
exit with an error. If we update the test config from this:

```yaml
  script.hello_with_type:
  - args:
      name: 1
```

to this (to indicate that we expect failure):

```yaml
  script.hello_with_type:
  - args:
      name: 1
  success: false
```

the tests will pass!

**exists**

And finally, if our function saved a file, we'd want to check for that like this.
The following example checks that whatever filename is returned does exist:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    exists: {% returns %}
```

A previous example showed how you could reference a specific input variable,
"outfile" existing. Since we also used the function `{% tmp_file %}` this output
file will be cleaned up after the fact.

```yaml
  script.write:
    args:
    - name: vanessa
    - outfile: {% tmp_file %}
    exists: {{ args.outfile }}
```


This means that we can edit our script from this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/basic/script.py
  tests:
    script.add:
    - args:
        one: null
        two: null
    script.add_with_type:
    - args:
        one: null
        two: null
    script.hello:
    - args:
        name: null
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: null
```

to be something more reasonable to test:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/basic/script.py
  tests:
    script.add:
    - args:
        one: 1
        two: 2
      returns: 3
    - args:
        one: 1
        two: null
      raises: TypeError
    script.add_with_type:
    - args:
        one: 1
        two: 2
      returns: 3
    script.hello:
    - args:
        name: Vanessa
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: 1
```

For typing, given that a function uses typing, that will be tested. For example,
the last function "hello_with_type" will not be successful.

## Test

Finally, you'll have your test file, and an environment where you want to
test. You can run tests like this:

```bash
$ gridtest test gridtest.yml
```
Or since gridtest.yml is the default, just do:

```bash
$ gridtest test
```

And here is an (under development) snapshot of what a result currently looks like
(and this particular shot was run in serial).

![img/failed.png](img/failed.png)

And here is an example of when all tests pass:

![img/success.png](img/success.png)


### Verbose

You can print a little more output about successes or failures with `--verbose`

```bash
$ gridtest test --verbose
[6/6] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        Exception: TypeError raised as desired. raises TypeError
script.add_with_type.0         success                        returns 3                     
script.hello.0                 success                        hello Vanessa!                
script.hello_with_default.0    success                        hello Dinosaur!               
script.hello_with_type.0       success                        success key set to false, expected failure.

```

Or you can filter to a regular expression (pattern) to only run a subset of
tests:

```bash
$ gridtest test --pattern script.add 
[3/3] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        raises TypeError              
script.add_with_type.0         success                        returns 3                     

3/3 tests passed
```

These are of course very simple test cases - we don't have any
classes, custom types, or matrices of tests. These will be developed
and discussed in other examples.c
---
title: GridTest
permalink: /
---

> What is GridTest?

GridTest is a library that specializes in generating parameter grids. The grids
are most obviously used for testing, but can extend to other use cases.
In the context of testing, GridTest makes it easy to discover functions,
classes, and arguments for your python scripts or modules, and then generate
a template for you to easily populate. Outside of testing, you can define
grids that are version controlled, programatically defined with functions,
and easy to interact with from the command line or Python interpreter.
Gridtest treats grids as first class citizens, and not as "the other thing
I need when I write tests."

> How might I use gridtest?

GridTest has use cases well beyond testing, because parameterization is used
widely across data science, and version control for reproducibility of those
parameterization is essential for reproducible, sustainable science.
This means that you might use GridTest for...

<a id="testing">
### Testing

A **gridtest**: is one that is run over a grid of parameter settings. Each test
can include a set of argument specifications, and optionally mapping these arguments
to functions so they can be programatically defined. 
A grid can be inline to the test (if not used elsewhere) or defined globally and shared.

<a id="grids">
### Parameterization

A **grid** is a global definition of a parameter matrix. While you are probably familiar
with a [traditional](https://scikit-learn.org/stable/modules/generated/sklearn.model_selection.ParameterGrid.html) definition of a grid, GridTest extends the idea and functionality of grids to include:

 - [generating random samples](https://vsoch.github.io/gridtest/tutorials/samplegrid/)
 - [loading grids via a GridRunner](https://vsoch.github.io/gridtest/getting-started/grids/index.html#loading-via-a-gridrunner) class separate from any Python code.
 - generating grids as you go (meaning as an iterator)
 - previewing grids on the command line before you use them
 - generating content of grids via external functions, and optionally unwrapping list values

Grids are generated on demand, meaning when you iterate over a grid object so that they are more
optimal to use because we don't save any single, large list to memory.
You can generally extend a grid to any use case that requires some combination of define arguments, 
and optionally functions to run to be mapped to arguments.
Grids do not have to be used in testing! You might share a repository that only defines grids that people
can use across many different kinds of machine learning models, likely to run metrics.

<a id="metrics">
### Metrics

A **metric** is a Python decorator that is paired with a test that will measure some
attribute of a test. For example:
   - you might run a function across a grid of arguments, and then measure the time that each combination takes (the metric), and generate a report for inspection.
   - you might be doing text processing and having functions to parse text. Each function might be run over a grid of sentences and counts, and for each result, we want to count the number of unique words, and total words (metrics). This is the [interface example](examples/interface).

Take a look at the [examples](examples) folder or the [documentation](https://vsoch.github.io/gridtest) for getting started. An example report is available to view [here](https://vsoch.github.io/gridtest/templates/report/),
and as we get more real world use cases, the report templates and data export options will be expanded
to use and visualize them beautifully. Please [open an issue](https://github.com/vsoch/gridtest/issues) 
if you have a use case that @vsoch can help with!

<a id="kind-of-projects">

> What kind of projects is GridTest for?

If you manage a large open source project with a diverse community, and you
are looking for a testing framework, you probably don't need GridTest. 
<a href="https://docs.pytest.org/en/latest/" target="_blank">pytest</a> is
really good at this. GridTest is designed for smaller or individual-run
Python projects that don't have the bandwidth to invest time into writing tests, as it
makes it easy to generate and fill in a template to run tests,
and run them during continuous integration.

<a id="just-for-testing">
> Is gridtest just for testing?

As mentioned above, one of the powerful features of GridTest is the 
ability to define [grids](https://vsoch.github.io/gridtest/getting-started/grids/index.html) 
that can be used within gridtests, but also outside of them. For example, you can easily define sets of named
parameterizations over sets of variables, and even map arguments to functions. 
The result can be an expanded list of arguments, or of results.

> What makes GridTest easy to use?

GridTest nicely manages small annoyances for writing tests, or generating grids of parameters.
The author has listed these in the order that she finds interesting:

<a id="grids-as-first-class-citizens">
### 1. Grids as First Class Citizens

Parameters always come as a second thought when writing tests, and this is
why they are commonly applied as decorators. The author of this software
realized that she might want to define just sets of parameters that expand
into matrices that can be useful across many use cases. This makes
the grids first class citizens.

<a id="capturing-metrics">
### 2. Capturing Metrics

How long does your function take when you provide parameter X as one value, versus
another? By way of allowing you to specify one or more metrics alongside tests,
you can easily capture metrics (Python decorators to your functions to test)
to output in an interactive report.

<a id="generating-reports">
### 3. Generating Reports

If you need to save results to a data file (e.g., results.json) or generate
an interactive report for GitHub pages, this is easy to do do with running
Gridtest with the `--save` or `--save-web` flags. An example report is 
available to view [here](https://vsoch.github.io/gridtest/templates/report/),
and as we get more real world use cases, the report templates and data export 
options can be expanded to use and visualize them beautifully. Please [open an issue](https://github.com/vsoch/gridtest/issues) if you have a use case that @vsoch can help with!

<a id="debugging">
### 4. Debugging

What programmer hasn't been in the scenario of running a group of tests,
and then having some fail? What do you do in that case? You can start an interactive
shell, import what you need, and try to reproduce, or you can turn up verbosity
and add a bunch of print statements to figure out what is going on. GridTest makes
this much easier with it's `--interactive` mode, which will let you simply
shell into an interpreter right before the function is run, and let you debug 
away.

<a id="reproducible-tests">
### 5. Running Reproducible Tests

When you write tests for a file, local, or system module, you store them in
a yaml file that is stored alongside the code, and can be tested with CI.
The yaml file can have grids of parameters defined so you can easily test many
different combinations.

<a id="knowing-tests">
### 6. Knowing the tests to write

Whether you write as you go or at the end, you have to look back at your files
to know the functions names and arguments that need to be tested. GridTest solves
this problem by way of discovery - give it a module name, a file name, or
an entire directory with Python files, and it will generate a template for you
to easily fill in that already includes arguments and functions. 

<a id="knowing-new-tests">
### 7. Knowing new tests to write

Okay great, so you've already written your tests. What if you add a function,
and haven't written tests for it yet? GridTest can tell you this too with it's
`--check` feature. It will let you run it against your previously generated file
and tell you exactly the functions that need to be added. Then remove `--check`
and it will add them.


In summary, GridTest:

 1. Let's you define grids to be generated programatically, version controlled, and used for multiple purposes
 2. Allows measuring of metrics alongside tests
 3. Stores tests in a yaml file that can be stored in version control
 4. Generates data exports and interactive reports for results
 5. Provides an easier way to interactively debug
 6. Helps you to discover the tests that you need to write, and creates a template to fill in
 7. Makes it easy to define and interact with expanded parameter grids


> Where do I go from here?

A good place to start is the [getting started]({{ site.baseurl }}/getting-started/) page,
which has links for getting started with writing tests, running tests, and many examples.
---
title: Not Found
permalink: /404.html
sitemap: false
---

This page doesn't exist!
---
title:
type: major
---

This release introduces

**Features:**

*

**Fixes:**

*
---
title:
category:
order: 1
---
---
title: Contributing
category: Contributing
permalink: /contributing/index.html
order: 1
---

It's so great that you want to contribute! There are several ways to contribute.

 - [Contribute to Documentation]({{ site.baseurl }}/contributing/docs/)
 - [Contribute Code]({{ site.baseurl }}/contributing/code/)
 - Ask a question, report a bug, or suggest a feature by opening up an [issue](https://github.com/{{ site.repo }}/issues)

Whether you ask a question, help to update the code with a bug fix, or
simply want to say hello, your contribution is appreciated.

---
title: Contributing to Documentation
category: Contributing
order: 3
---

It's so great that you want to contribute! The documentation here includes information
about using and developing {{ site.title }}, and they are hosted on Github, meaning that you
can easily contribute via a [pull request](https://help.github.com/articles/about-pull-requests/).

## Getting Started

### Installing Dependencies

Initially (on OS X), you will need to setup [Brew](http://brew.sh/) which is a 
package manager for OS X and [Git](https://git-scm.com/). To install Brew and Git, 
run the following commands:

```bash
/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
brew install git
```

If you are on Debian/Ubuntu, then you can easily install git with `apt-get`

```bash
apt-get update && apt-get install -y git
```

### Fork the repo

To contribute to the web based documentation, you should obtain a GitHub account and *fork* the <a href="https://github.com/{{ site.repo }}" target="_blank">{{ site.title }} Documentation</a> repository by clicking the *fork* button on the top right of the page. Once forked, you will want to clone the fork of the repo to your computer. Let's say my GitHub username is *meatball*:

```bash
git clone https://github.com/meatball/{{ site.reponame }}
cd {{ site.reponame }}/
```

### Install a local Jekyll server
This step is required if you want to render your work locally before committing the changes. This is highly recommended to ensure that your changes will render properly and will be accepted.

```bash
brew install ruby
gem install jekyll
gem install bundler
bundle install
```

The documentation is located in the "docs" subfolder, so after cloning the repository,
you can change the directory to there and then run the jekyll serve.
Specifically, you can see the site locally by running the server with jekyll:

```bash
cd docs
bundle exec jekyll serve
```

This will make the site viewable at <a href="http://localhost:4000/{{ site.title }}/" target="_blank">http://localhost:4000/{{ site.title }}/</a>.
---
title: Github Contribution
category: Contributing
order: 3
---

To contribute to code, you should first *fork* the <a href="https://www.github.com/{{ site.repo }}" target="_blank">{{ site.title }}</a> repository by clicking the *fork* button on the top right of the page. Once forked, you will want to clone the fork of the repository to your computer:

```bash
git clone https://github.com/<username>/{{ site.reponame }}
cd {{ site.reponame }}/
```

The main python module, gridtest, is in the top level folder. You should checkout a branch,
push the branch to your remote, and when you are ready, open a pull request against
the master branch of {{ site.repo }}.

```bash
git checkout -b add/my-feature
git commit -a -m 'adding my new feature!'
git push origin add/my-feature
```
---
title: Introduction
category: Getting Started
permalink: /getting-started/index.html
order: 1
---

You should first [install]({{ site.baseurl }}/install/) gridtest.
This will place the executable `gridtest` in your bin folder, which is the client
for generating and running tests. 

## Getting Started

### Introduction

 - [How does it work?](#how-does-it-work): How does gridtest work?
 - [Concepts](#concepts): What are common gridtest concepts?

### Writing Tests

 - [Templates](templates/): for test yaml files, including function and argument substitution
 - [Grids](grids/): can be used to programatically generate inputs for tests, or outside of testing when you want to parameterize some values, optionally mapping args to one or more functions.
 - [Metrics](metrics/): decorators to measure metrics across a grid of tests.

### Running Tests

 - [Environment](environment/): variables to change defaults for gridtest behavior
 - [Debugging](#debugging): is easy with the interactive interpreter.
 - [Testing](testing/): via continuous integration, or checking if tests need updating.
 - [Python](python/): interacting with a GridRunner from within Python
 - [Results](results/) are available via a json export or interactive web report.

<a id="#how-does-it-work">
### How does it work?

GridTest takes advantage of the <a href="https://docs.python.org/3/library/inspect.html">inspect
module</a> in Python so that you can provide any Python file, local, or system installed module
and extract not only a list of functions and classes, but also arguments. By way of providing
a file, module name, or folder path, GridTest can first generate for you a "test template," 
which comes down to a yaml file that you can easily fill in and customize for your tests.
For example, I might generate a test yaml file for my script like this:

```bash
$ gridtest generate script.py gridtest.yml
```

And you would then open the file and customize the default template to your liking.
The testing file can then be run, in the same way, using gridtest:

```bash
$ gridtest test gridtest.yml
```

Or since gridtest.yml is the default, just leave it out to find the file in
the present working directory:

```bash
$ gridtest test
```

The tests are run in parallel with multiprocessing to be efficient, but if you want
you can specify running in serial:

```bash
$ gridtest test --serial
```

Along with testing the conditions that you specify in the testing file, GridTest will also run type checking
if you have defined types in your code. You can also check if you need to add
more tests:

```bash
# exit code is 0 if all tests written, 1 otherwise
$ gridtest check gridtest.yml

No new tests to add!✨ 🥑️ ✨
```

And if there are new tests (meaning new functions or classes in the file that aren't
tested) you can update your testing file with those missing:

```bash
$ gridtest update gridtest.yml
```

If you don't care about testing, you can use GridTest to generate a yaml specification
of parameterizations. Let's say we have a set of grids defined in grids.yml.
We can list the named grids that are defined in the file:

```bash
$ gridtest gridview grids.yml --list
generate_empty
generate_matrix
generate_lists_matrix
generate_by_min_max
generate_by_min_max_twovars
```

And then either print them all to the screen:

```bash
$ gridtest gridview grids.yml
```

or just print a specific grid we found with `--list`.

```bash
$ gridtest gridview grids.yml generate_by_min_max
{'x': 0.0}
{'x': 2.0}
{'x': 4.0}
{'x': 6.0}
{'x': 8.0}
```

For larger grids, it's nice to be able to get a count. The "twovars"
variant of the function above adds another variable (also with 5 values)
so we get 5x5:

```bash
$ $ gridtest gridview grids.yml generate_by_min_max_twovars --count
25 lists produced.
```

The rest of this getting started guide will review overall functionality. 
You should look at [tutorials]({{ site.baseurl }}/tutorials/) 
if you have never done this before, and would like to
walk through basic examples.

<a id="concepts">
### What are gridtest concepts

<a id="gridtest">
#### GridTest

A gridtest is a yaml file that is executed using the software described here,
GridTest. It's called a GridTest because of the ability to define grids of tests.
For example, let's say that we want to measure the time it takes to run a function
over a set of variables. We do this by way of adding an optimization, which
comes down to a decorator function that will measure the value and return it to
gridtest. We then might specify a range of numbers for one or more function
variables. This will ultimately generate a grid of tests for the function,
each with a set of results from the function itself and decorators. 
GridTest also will make it easy to run over a grid of different environments,
although this is not developed yet.

<a id="gridtest.yml">
#### gridtest.yml

GridTest's main convention is that  it will look for a yaml file, `gridtest.yml` to run tests by default. This means
that if you write your testing file in the root of a repository with Python
software like:

```
setup.py
requirements.txt
gridtest.yml
mymodule/
  subfolder/
  __init__.py
```

You can easily run `gridtest test` in the root of that folder to discover the
file. This is similar to a Dockerfile, or Singularity recipe file.

<a id="grids">
#### grids

Gridtest does parameterization by way of grids, which is a section of the gridtest.yml
(alongside tests) that has definitions for one or more named grids. For example,
if we wanted to use a function (`script.get_pokemon_id`) to generate a list of
values for the argument `pid` for a test called `script.generate_pokemon`, we might
write a recipe like this. 

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/grid-function/script.py
  grids:
    generate_pids:
      func: script.get_pokemon_id
      count: 10

  tests:
    script.generate_pokemon:
        grid:
          pid: generate_pids
```

The reference to "generate_pids" for the argument "pid" is referenced under
grids, and we know to use the `scripts.get_pokemon_id`, run 10 times, to
generate it. We could also define arguments for the parameterization. For
a simple example, see the [grid-function](https://github.com/vsoch/gridtest/tree/master/examples/grid-function) 
example, or read more in the getting started [templates](https://vsoch.github.io/gridtest/getting-started/templates/index.html) guide.

<a id="#debugging">
### Debugging

For most testing frameworks, when you hit an error there is a flash of red across the
screen, and at best you can stop execution on the first error (e.g., with pytest you would
add the -x flag) or add print statements to the test to see what might be going on. 
You might start an interactive shell to import modules and functions needed to debug,
and then need to exit and run again to re-run the test. Instead of these approaches,
gridtest gives you an ability to run a test with `--interactive`, meaning that you
can open an interactive shell at the onset of each test.

```python
$ gridtest test examples/basic/script-tests.yml --interactive
[script.add:1/6] |=====|-----------------------------|  16.7% 

Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function add at 0x7fe4c0a44200>
 module: <module 'script' from '/home/vanessa/Desktop/Code/gridtest/examples/basic/script.py'>
   args: {'one': 1, 'two': 2}
returns: 3

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:                                                                     
```

In the example above, we show that the function, module, and arguments for
the test are loaded, and you are shown how to run the tests. If you have added any
decorators (optimizations to measure) they will be applied already to the variable 
func. By way of having the interactive terminal, you can of course interact with functions and variables
directly, and debug what might be the issue for a test. In the case of a file
with multiple tests (the typical case) You can also specify the name of the test you want
to interact with:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add
```

For the above, this would interact with all tests that start with script.add. If you
want to limit to the script.add module, you might want to do:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add.
```

Or a specific indexed text for the module:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add.0
```

For a more detailed debugging example, see the [debugging]({{ site.baseurl }}/getting-started/debugging/)
documentation.

## Licenses

This code is licensed under the Mozilla, version 2.0 or later [LICENSE](LICENSE).

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Testing
category: Getting Started
permalink: /getting-started/testing/index.html
order: 6
---

You likely want to integrate GridTest into your favorite continuous integration
(CI) service and there are many good ways to do that.

## GridTest Linting

If you add a new function to your testing library, you probably
want to ensure that you've written a test for it. Akin to black for code
formatting, gridtest provides a gridtest check command
to check if a file has all tests written. In the case that all functions
are represented, the return code is 0 (and the CI will pass):

```bash
$ gridtest check tests/modules/temp-tests.yml 

No new tests to add!✨ 🥑️ ✨

# Return code is 0
$ echo $?
0
```

But now let's say that we update the module associated with the testing file
and add a new function.

```python
def new_function():
    print("No test for me!")
```

And then we run checks again:


```bash
$ gridtest check tests/modules/temp-tests.yml 

New sections to add:
tests.modules.temp.new_function
```

Since we have a new test (and missed it) the return code is exit with an error (1):

```bash
echo $?
1
```

If we wanted to skip this test, we could do that by defining a skip pattern:

```bash
$ gridtest check tests/modules/temp-tests.yml --skip-patterns new_function

No new tests to add!✨ 🥑️ ✨
```

And the return code is again 0 (meaning it would pass in CI).

```bash
$ echo $?
0
```

## GridTest Updating

Let's take this same file, and say that we've checked it, and found that a new test
needs to be added. How do we add the function to the template? You can run
the update command:

```bash
$ gridtest update tests/modules/temp-tests.yml 
Adding function temp.new_function
Writing [gridtest|temp-tests.yml] to tests/modules/temp-tests.yml
```

Update is similar to generate, but we already have the filename (for the file
or module) added to the testing file and don't need to specify it again.
The test file is updated with the new function (which is rather boring, it doesn't
have any arguments)

```yaml
temp:
  filename: /home/vanessa/Desktop/Code/gridtest/tests/modules/temp.py
  tests:
    temp.create_directory:
    - args:
        dirname: '{% raw %}{% tmp_dir %}{% endraw %}'
    temp.new_function:
    - args: {}
    temp.write_file:
    - args:
        filename: '{% raw %}{% tmp_path %}{% endraw %}'
```

## Continuous Integration Recipes

Gridtest has templates available (CI services added on request) 
for particular continuous integration providers. Generally, you can use the 
above `gridest check <testfile>.yml` to check for needing to write tests, 
and then the standard grid runner (with any options you need) to run the tests, 
after installing gridtest. For example, in a run statement you might do:

```bash
pip install gridtest
gridtest run testfile.yml
```

### GitHub Workflows

An example GitHub workflow is provided at [github-workflow](https://github.com/vsoch/gridtest/blob/master/.github/workflows/github-workflow-example.yml), and run with this repository. It will run the grid of
tests for the [interface] example:

![img/github-workflow.png](../img/github-workflow.png)

 and then upload the results report as an artifact.

![img/github-workflow-artifact.png](../img/github-workflow-artifact.png)

You could imagine also committing the static files to a docs folder, and then
opening a pull request (or pushing directly) to update GitHub pages. If you
want some help setting this up, please don't be afraid to [reach out](https://github.com/{{ site.repo }}/issues).
---
title: Debugging
category: Getting Started
permalink: /getting-started/debugging/index.html
order: 2
---

### Why Debugging?

For most testing frameworks, when you hit an error there is a flash of red across the
screen, and at best you can stop execution on the first error (e.g., with pytest you would
add the -x flag) or add print statements to the test to see what might be going on. 
You might start an interactive shell to import modules and functions needed to debug,
and then need to exit and run again to re-run the test. Instead of these approaches,
gridtest gives you an ability to run a test with `--interactive`, meaning that you
can open an interactive shell at the onset of each test.

### An example Python script

Let's say that we have this script, called "script.py"

```python
# These are functions in my script
# Typing is here, so Python 

def add(one, two):
    """add will add two numbers, one and two. There is no typing here"""
    return one + two

def add_with_type(one: int, two: int) -> int:
    """add_with_type will add two numbers, one and two, with typing for ints."""
    return one + two

def hello(name):
    """print hello to a name, with no typing"""
    print(f"hello {name}!")

def hello_with_default(name="Dinosaur"):
    """print hello to a name with a default"""
    print(f"hello {name}!")

def hello_with_type(name: str) -> None:
    """print hello to a name, with typing"""
    print(f"hello {name}!")
```

And from it we've produced this testing file "script-tests.yml" in the same directory:

```yaml
script:
  filename: script.py
  tests:
    script.add:
    - args:
        one: 1
        two: 2
      returns: 3
    - args:
        one: 1
        two: null
      raises: TypeError
    script.add_with_type:
    - args:
        one: 1
        two: 2
      returns: 3
    script.hello:
    - args:
        name: Vanessa
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: 1
      success: false
```

### Running Tests

We would run it with gridtest as follows:

```bash
$ gridtest test script-tests.yml 
[6/6] |===================================| 100.0% 
success: script.add.0 returns 3 
success: script.add.1 raises TypeError 
success: script.add_with_type.0 returns 3 
success: script.hello.0 
success: script.hello_with_default.0 
success: script.hello_with_type.0 
6/6 tests passed
```

Now let's say there is an error in a script. Let's randomly raise an exception:

```python
def hello(name):
    """print hello to a name, with no typing"""
    raise Exception('ruhroh')
```

If we run tests again, we see a failure with an unexpected exception:

```bash
$ gridtest test script-tests.yml 
[6/6] |===================================| 100.0% 
success: script.add.0 returns 3 
success: script.add.1 raises TypeError 
success: script.add_with_type.0 returns 3 
failure: script.hello.0 ruhroh Unexpected Exception: Exception.
success: script.hello_with_default.0 
success: script.hello_with_type.0 
5/6 tests passed
```

### Adding Interactivity

If we add the `--interactive` flag, it's going to allow us to cycle through
*every single test* and press Control+d to jump to the next test:

```python
$ gridtest test script-tests.yml --interactive
[script.add:1/6] |=====|-----------------------------|  16.7% 

Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function add at 0x7fe4c0a44200>
 module: <module 'script' from '/home/vanessa/Desktop/Code/gridtest/examples/basic/script.py'>
   args: {'one': 1, 'two': 2}
returns: 3

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:                                                                     
```

But that's not really what we want - we know the failing test is `script.hello.0`
so let's run the tests, but only this particular test for interactive:

```python
$ gridtest test script-tests.yml --interactive --name script.hello

[script.add:1/6] |=====|-----------------------------|  16.7% False
[script.add:2/6] |===========|-----------------------|  33.3% False
[script.add_with_type:3/6] |=================|-----------------|  50.0% False
[script.hello:4/6] |=======================|-----------|  66.7% True


Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function hello at 0x7fc118e8c680>
 module: <module 'script' from '/home/vanessa/Desktop/Code/gridtest/examples/basic/script.py'>
   args: {'name': 'Vanessa'}
returns: None

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:
```

In the example above, we show that the function, module, and arguments for
the test are loaded, and you are shown how to run the tests. For example,
here is how we would inspect arguments and then test typing:

```python
In [1]: args                                                                                                                                 
Out[1]: {'name': 'Vanessa'}

In [2]: passed, error = test_types(func, args, returns)                                                                                      

In [3]: passed                                                                                                                               
Out[3]: True

In [4]: error                                                                                                                                
Out[4]: []
```

And then run the actual test to trigger the full error:

```python
In [6]: result = func(**args)                                                                                                                
hello Vanessa!
---------------------------------------------------------------------------
Exception                                 Traceback (most recent call last)
~/Desktop/Code/gridtest/gridtest/main/helpers.py in <module>
----> 1 result = func(**args)

~/Desktop/Code/gridtest/examples/basic/script.py in hello(name)
     18     """print hello to a name, with no typing"""
     19     print(f"hello {name}!")
---> 20     raise Exception('ruhroh')
     21 
     22 def hello_with_default(name="Dinosaur"):

Exception: ruhroh
```

At this point, you can interact with your arguments or the function to debug further.

### Specifying Names

In the case of a file with multiple tests (the typical case) You can also specify the name of the test you want
to interact with:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add
```

For the above, this would interact with all tests that start with script.add. If you
want to limit to the script.add module, you might want to do:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add.
```

Or a specific indexed text for the module:

```python
$ gridtest test examples/basic/script-tests.yml --interactive --name script.add.0
```

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Environment
category: Getting Started
permalink: /getting-started/environment/index.html
order: 5
---

Gridtest supports export of a few environment variables that can drive
testing or other functionality.

## Gridtest Workers

By default, the number of workers will be the number of processes (nproc) multiplied
by 2 plus 1. This isn't some hard defined standard, but what was found to work
relatively well. You can change this by exporting the variable `GRIDTEST_WORKERS`:

```bash
export GRIDTEST_WORKERS=2
```

Note that workers are only relevant for tests run with multiprocessing (the default)
so if you add the `--serial` flag, or if you are using `--interactive` mode (which
also requires running in serial) this variable will not be relevant.

## GridTest Shell

If you use the `gridtest shell` mode to interactively create a gridtest running
in Python, it will default to looking for an ipython installation, and then check
for standard python, and then bpython. If you want to change the default,
just export the string for the python that you want:

```bash
# one of ipython, python, or bpython
export GRIDTEST_SHELL=ipython
```

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Templates
category: Getting Started
permalink: /getting-started/templates/index.html
order: 3
---

The power of your recipe comes down to your template. Knowing how input
and return types work is thus essential for writing awesome templates! The basic
recipe format has a module name as the first key, and then a list of [tests](#testing-types),
a filename, and [grids](#grids).

```yaml
script:
  filename: script.py
  tests:
     ...
  grids:
     ...
```

Each of these sections will be discussed below.

<a id="testing-types">
## Testing Types

<a id="basic-test">
### Basic Test

The default (most basic kind of test) that gridtest can run is to take some function
name, a set of arguments, and then just test that it runs:

```yaml
  script.write:
    args:
    - name: "dinosaur"
```

In the above snippet, if `script.write` runs with the
input "dinosaur" as name, this test will be successful. The next type is what
gives gridtest it's name, the "grid" specification. This is when we don't want
to define a single argument, but some set of parameters for gridtest to iterate
over. As an example, let's say that we have a function that takes an input, `seconds`, and sleeps
for that many seconds. Our default might start like this:

```yaml
  script.gotosleep:
    args:
      - seconds: 1
```

The above would run one test, and sleep for 1 second. But that's not really so
useful. We would want to define a range of values between 0 and 5, 
and then a few explicit higher values, 10 and 11. How would that look?

<a id="grid-test">
### Grid Test

We can accomplish better parameterization by using a grid in our test. That
looks like this:

```yaml
  script.gotosleep:
    grid:
      seconds:
        max: 5
        min: 0
        list: [10, 15]
```

In the example above, the previous argument (arg) has been moved to a section
called "grid" to indicate to gridtest that this is a grid of parameters to run
over, and not a one off value. Under grid we have an equivalent entry for
seconds, but this time, we define a min (0), a max (5), and a list to include
10 and 15. Gridtest would parse this to run tests over our sleep function
for all of these argument sets:

```python
[{'seconds': 0}, {'seconds': 1}, {'seconds': 2}, {'seconds': 3}, {'seconds': 4}, {'seconds': 10}, {'seconds': 15}]
```

And if we had other grid parameters defined, we'd build a matrix over them too.
Single values can remain in the args section, but **you are not allowed to have
the same parameter defined under both args and grid**. To be explicit about the
grid section:

**min** and **max** are to be used when specifying a range. When unset, **by** 
would be 1. If you want to decrease, set a negative value for by. You can assume
the values are going into range like `range(min, max, by)`.

**list** is for when you want to include a list of values, even in addition to a
range already specified as in the example above.

An interactive result report is planned to better illustrate the output here.

## Input Types

Gridtest can support multiple different kinds of input types that correspond
with common use cases. For example, many tests want to use a temporary
file as input. We might do something like:

```yaml
  script.write:
    args:
    - name: "dinosaur"
    - outfile: {% raw %}{% tmp_path %}{% endraw %}
    returns: {% raw %}{{ args.outfile }}{% endraw %}
```

In the above, the double `{% raw %}{{}}{% endraw %}` refers to a variable, in this case which is
an argument name. The `{% raw %}{% %}{% endraw %}` refers to a function that is known natively
to grid test. This syntax is based on [jinja2](https://jinja.palletsprojects.com/en/2.11.x/).
You could also check that the file exists (and it might not be returned).

```yaml
  script.write:
    args:
    - name: vanessa
    - outfile: {% raw %}{% tmp_path %}{% endraw %}
    exists: {% raw %}{{ args.outfile }}{% endraw %}
```

### Input Variables

An input variable is distinguished by being in the format `{% raw %}{{ <name> }}{% endraw %}`
where `<name>` would be some named variable, referenced below.

| Name        | Description | Syntax | Example |
|-------------|-------------|--------|---------|
| args.<name> | Refer to a named argument under args. |`{% raw %}{{ args.<name> }}{% endraw %}` | `{% raw %}{{ args.name }}{% endraw %}` |
| returns | Refer to the returns value you defined. |`{% raw %}{{ returns }}{% endraw %}` | `{% raw %}{{ returns }}{% endraw %}` |

If you want to refer to a results object just refer to the GridTest instance (self) directly. 
For example, to test that a requests.response_code is equal to a certain value, we can do:

```yaml
  requests.api.head:
  - args:
      url: https://google.com
    istrue: "self.result.status_code == 301"
```



### Function Variables

In addition to input variables, gridtest provides a few functions to make testing easier. Functions
are distinguished based on being in the format `{% raw %}{% <func> %}{% endraw %}`,
where "func" refers to the name of the function. Gridtest currently supports the following functions:

| Name        | Description | Syntax |
|-------------|-------------|--------|
| tmp_dir | Create (and cleanup) a temporary directory for testing. |`{% raw %}{% tmp_dir %}{% endraw %}` |
| tmp_path | Create (and cleanup) a temporary filename for testing. |`{% raw %}{% tmp_path %}{% endraw %}`|

For specifics about tmp_dir and tmp_path, see the [temp tutorial]({{ site.baseurl }}/tutorials/temp/).


### Return Types

For basic testing, there are typically a few obvious cases we want to test for:

 - returns: meaning that a function returns a specific value
 - raises: meaning that the function raises an error
 - exists: meaning that some output file is deemed to exist.
 - success: boolean to indicate if success (no error) is desired.
 - istrue: determine if a custom evaluated statement is true
 - isfalse: determine if a custom evaluated statement is false
 - equals: determine if a custom evaluated statement is equal to the result

If you see another simple testing case that you want added, please 
[open an issue](https://github.com/vsoch/gridtest/issues). Highly complex testing needs
probably should use a more substantial testing library. Let's look through how each of
these examples might be used for our add function.

**returns**

If we want to ensure that the correct value is returned, we would do:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    returns: 3
```

**raises**

If we wanted to ensure that an exception was raised, we would do:

```yaml
  script.add:
  - args:
      one: 1
      two: null
    raises: TypeError
```

**istrue**

istrue is used when we want to check if something is True.
You usually would want to refer to an input or output variable:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    istrue: isinstance({% raw %}{% returns %}{% endraw %}, int)
```

**isfalse**

or you might want the opposite, isfalse:


```yaml
  script.add:
  - args:
      one: 1
      two: 2
    isfalse: not isinstance({% raw %}{% returns %}{% endraw %}, float)
```

**equals**

or you might want to evaluate a statement. In the example below, we want to 
make sure an attribute of the value returned is equal to 200.

```yaml
  requests.get:
  - args:
      url: https://google.com
      data: null
      json: null
    eval: {% raw %}{% returns %}{% endraw %}.status_code == 200
```

This is different from providing an explicit value.

**success**

You might just want to run something, and be sure that the success status is False.
For example, if you give the wrong type as input to a function, it will by default
exit with an error:

```bash
$ gridtest test examples/basic/script-tests.yml 
[script.hello_with_type:6/6] |===================================| 100.0% 3% 
success: script.add.0 returns 3 
success: script.add.1 raises TypeError 
success: script.add_with_type.0 returns 3 
success: script.hello.0 
success: script.hello_with_default.0 
failure: script.hello_with_type.0 TypeError name (1) is <class 'int'>, should be <class 'str'>
```

However, if we update the test config from this:

```yaml
  script.hello_with_type:
  - args:
      name: 1
```

to this (to indicate that we expect failure):

```yaml
  script.hello_with_type:
  - args:
      name: 1
  success: false
```

the tests will pass!

**exists**

And finally, if our function saved a file, we'd want to check for that like this.
The following example checks that whatever filename is returned does exist:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    exists: {% raw %}{% returns %}{% endraw %}
```

A previous example showed how you could reference a specific input variable,
"outfile" existing. Since we also used the function `{% raw %}{% tmp_path %}{% endraw %}` this output
file will be cleaned up after the fact.

```yaml
  script.write:
    args:
    - name: vanessa
    - outfile: {% raw %}{% tmp_path %}{% endraw %}
    exists: {% raw %}{{ args.outfile }}{% endraw %}
```

**isinstance**

To check that a result is of a particular instance type, you can use `isinstance`
and provide the name of the class that would be returned as a string with `type(self.result).__name__`.
For example, to test if a custom Car instance is of type car.Car, we would do:

```yaml
  car.Car:
  - instance: thisone
    args:
      color: red
      lights: false
      wheels: 4
    isinstance: car.Car 
```

### Grids

If you've already read about [grids](../grids), you know that grids can be defined
in the context of testing to run tests across a grid of input arguments. A good
example of adding a grid is provided in the 
[metrics example](https://vsoch.github.io/gridtest/getting-started/metrics/#adding-a-grid).



You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Grids
category: Getting Started
permalink: /getting-started/grids/index.html
order: 6
---

GridTest isn't just for testing! In fact, you can write a file of grid
specifications that can be loaded and used for parameterization, or your own
custom functions. The description here will follow example the files in the
[grids](https://github.com/vsoch/gridtest/tree/master/examples/grids) examples folder.
Here is a peek at the top of this file, where we define two grids:

```yaml
mygrids:
  grids:

    # A grid that will generate 10 empty set of arguments
    generate_empty:
      count: 10

    # A grid that will generate each cross of x and y (total of 9)
    generate_matrix:
      args:
        x: [1, 2, 3] 
        y: [1, 2, 3] 
```

## Loading via a GridRunner

Let's say that we have a grids.yml file with grids defined. We can load with
gridtest easily via the GridRunner:

```python
from gridtest.main.test import GridRunner
runner = GridRunner('grids.yml')
# [gridtest|grids.yml]
```

We can easily generate the grids as follows:

```python
runner.get_grids()
{'generate_empty': [grid|generate_empty],
 'generate_matrix': [grid|generate_matrix],
 'generate_lists_matrix': [grid|generate_lists_matrix],
 'generate_by_min_max': [grid|generate_by_min_max],
 'generate_by_min_max_twovars': [grid|generate_by_min_max_twovars]}
```

## Loading as a Grid

A Grid is a first class citizen, so you can also generate without
the GridRunner, either via a json object or yaml you've loaded:

```python
from gridtest.main.grids import Grid

grid = Grid(name="mygrid", params={"args": {"one": [1, 11, 111], "two": 2}})
[grid|mygrid]
```

### Cached Loading

By default, grids are generated on demand, meaning that the argument sets (optionally
derived by functions you have provided under params) are generated at this point in time.
For example:

```
for argset in grid:
    print(argset)

{'one': 1, 'two': 2}
{'one': 11, 'two': 2}
{'one': 111, 'two': 2}
```

This also means that, by default, your grid.argsets will be empty when you
instantiate it.

```python
grid = Grid(name="mygrid", params={"args": {"one": [1, 11, 111], "two": 2}})
grid.argsets
[]
```

This is done intentionally because it can sometimes take a lot of memory to store
a long list of argument sets. However, if you want to generate the argument sets
on init, just set "cache" to True in your params:

```python
grid = Grid(name="mygrid", params={"cache": True, "args": {"one": [1, 11, 111], "two": 2}})
```

You should then be able to see the argsets ready for your use!

```python
grid.argsets
[{'one': 1, 'two': 2}, {'one': 11, 'two': 2}, {'one': 111, 'two': 2}]
```

### Grid Functions

What if you want to generate a long list of items for an argument, and it would
be crazy to type them out? This is what the functions specification for a grid
is for. You would provide your functions, each associated with an argument, under
the "functions" section. As an example, let's say we want to run random.choice,
and select from a list defined under the argument "choices":

```python
import random
grid = Grid(name="mygrid", params={"functions": {"pid": random.choice}, "args": {"seq": [[1,2,3,4,5,6,7]]}})
```

The above would look like this in yaml:

```yaml
grids:
  mygrid:
    args:
      seq: [[1,2,3,4,5,6,7]]
    functions:
      pid: random.choice
```

And we would generate our argument sets:

```python
list(grid)
[{'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 2}]
```

This also points out another important point - if you were to provide a list for an
argument, it would be parameterized. If you want the entire list treated as one variable,
then define it within another list as shown above. The above grid says "generate the argument
pid by using random.choice to select from seq." How does the grid know to match seq to random.choice?
It's a known keyword argument! Another insight here is that although we define a sequence argument, we don't actually
need it, we really only care about the result (pid). But let's say we want to generate
10 x a pid, how do we do that? We add a count to the params:

```python
grid = Grid(name="mygrid", params={"count": 10, "functions": {"pid": random.choice}, "args": {"seq": [[1,2,3,4,5,6,7]]}})
list(grid)

[{'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 5},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 3},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 7},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 5},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 7},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 4},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 1},
 {'seq': [1, 2, 3, 4, 5, 6, 7], 'pid': 3}]
```

These values are being generated by an iterator and not saved to any single list in memory,
so although it seems very rundant to see the same list more than once, it shouldn't
serve significant issues with respect to memory.


## Viewing on the Command Line

You can also preview the grids generated from the command line. First you might
want to list all grids defined in a file:

```bash
$ gridtest gridview grids.yml
generate_empty
generate_matrix
generate_lists_matrix
generate_by_min_max
generate_by_min_max_twovars
```

and then inspect a specific named grid:

```bash
$ gridtest gridview grids.yml generate_empty
{}
{}
{}
{}
{}
{}
{}
{}
{}
{}
```

If you want to export your entire grid to json, possibly to read in as-is
later and use, add --export:

```bash
$  gridtest gridview grids.yml generate_by_min_max --export exported.json 
```
```bash
cat exported.json
[
    {
        "x": 0.0
    },
    {
        "x": 2.0
    },
    {
        "x": 4.0
    },
    {
        "x": 6.0
    },
    {
        "x": 8.0
    }
]
```

That's not a really interesting grid, but you get the gist.
Each of the grids listed above will be explained below in more detail.

## Writing Grids

Let's start with the most basic of grids, which are those that don't import any
special functions. 

### Header

The header should have a named section for the grids that you want to define
(e.g. `mygrids`) and then a `grids` index, under which we will write each
of our named grids. Here is the start:

```
mygrids:
  grids:
  ...
```

Then each grid is added as a named section below that:

```
mygrids:
  grids:

    # A grid that will generate 10 empty set of arguments
    generate_empty:
      count: 10
```

The full content of this file will be written to the [grids.yml](grids.yml)
example here. Each example grid is discussed below. For each example,
you can preview the grid in the terminal with `gridtest gridview` or
obtain the grid by instantiating the `GridRunner` as shown above.

### Empty

It could be that you need to generate empty lists of arguments for a function,
in which case the `count` variable will be useful to you. Here is
how to specify a grid that will generate 10 empty set of arguments

```yaml
generate_empty:
  count: 10
```

The result comes out to be:

```bash
$ gridtest gridview grids.yml generate_empty
{}
{}
{}
{}
{}
{}
{}
{}
{}
{}
```

### Parameterize Variables

Let's say that we have two variables, x and y, and we want to generate a grid
of all possible combinations for a listing of each. That would look like this:

```yaml
    generate_matrix:
      args:
        x: [1, 2, 3] 
        y: [1, 2, 3]
```

And the resulting grid will have 3x3 or 9 total combinations of x and y:

```bash
$ gridtest gridview grids.yml generate_matrix
{'x': 1, 'y': 1}
{'x': 1, 'y': 2}
{'x': 1, 'y': 3}
{'x': 2, 'y': 1}
{'x': 2, 'y': 2}
{'x': 2, 'y': 3}
{'x': 3, 'y': 1}
{'x': 3, 'y': 2}
{'x': 3, 'y': 3}
```

### Parameterize Lists

If you want to do similar but instead have a list of values be paramaterized, just 
specify a list of lists instead.

```yaml
    generate_lists_matrix:
      args:
        x: [[1, 2, 3], [4, 5, 6]] 
        y: [[1, 2, 3], [4, 5, 6]] 
```

The result will have 2x2 or 4 entries:

```bash
$ gridtest gridview grids.yml generate_lists_matrix
{'x': [1, 2, 3], 'y': [1, 2, 3]}
{'x': [1, 2, 3], 'y': [4, 5, 6]}
{'x': [4, 5, 6], 'y': [1, 2, 3]}
{'x': [4, 5, 6], 'y': [4, 5, 6]}
```

Here is an easier way to check the count:

```bash
$ gridtest gridview grids.yml generate_lists_matrix --count
4 argument sets produced.
```

### Range of Values

For most use cases, you'll want to generate a list of values over a range. You
can do that with **min** and **max** and (optionally) **by** that defaults to 1.

```yaml
    generate_by_min_max:
      args:
        x:
          min: 0
          max: 10
          by: 2
```

```bash
$ gridtest gridview grids.yml generate_by_min_max
{'x': 0.0}
{'x': 2.0}
{'x': 4.0}
{'x': 6.0}
{'x': 8.0
```

Here is an example with two variables:

```yaml
    generate_by_min_max_twovars:
      args:
        x:
          min: 0
          max: 10
          by: 2
        y:
          min: 10
          max: 20
          by: 2
```

```bash
$ gridtest gridview grids.yml generate_by_min_max_twovars 
{'x': 0.0, 'y': 10.0}
{'x': 0.0, 'y': 12.0}
{'x': 0.0, 'y': 14.0}
{'x': 0.0, 'y': 16.0}
{'x': 0.0, 'y': 18.0}
{'x': 2.0, 'y': 10.0}
{'x': 2.0, 'y': 12.0}
{'x': 2.0, 'y': 14.0}
{'x': 2.0, 'y': 16.0}
{'x': 2.0, 'y': 18.0}
{'x': 4.0, 'y': 10.0}
{'x': 4.0, 'y': 12.0}
{'x': 4.0, 'y': 14.0}
{'x': 4.0, 'y': 16.0}
{'x': 4.0, 'y': 18.0}
{'x': 6.0, 'y': 10.0}
{'x': 6.0, 'y': 12.0}
{'x': 6.0, 'y': 14.0}
{'x': 6.0, 'y': 16.0}
{'x': 6.0, 'y': 18.0}
{'x': 8.0, 'y': 10.0}
{'x': 8.0, 'y': 12.0}
{'x': 8.0, 'y': 14.0}
{'x': 8.0, 'y': 16.0}
{'x': 8.0, 'y': 18.0}
```

Logically there are the previous number of tests, but squared.

```bash
$ gridtest gridview grids.yml generate_by_min_max_twovars --count
25 argument sets produced.
```

## Grids with Functions

## What does it mean to use a function?

You can map functions to parametrs in a grid, meaning that the values for
the parameters will be generated by the function. For example,
let's use the function `random.choice` to dynamically generic a grid of parameters.
The grid below will call `random.choice` ten times (count is set to 10) 
across the sequence of values `[1, 2, 3]`, which is an input key word argument
to random choice.

```yaml
...
  grids:
    random_choice:
      count: 10
      functions: 
        pid: random.choice
      args:
         seq: [[1, 2, 3]]
```

Notice that the sequence argument input is a list of lists, and this is because we want the entire list to be treated as an argument. If we run this from it's respective file, we get a result with 10
argument sets:

```bash
$ gridtest gridview grids-with-function.yml random_choice
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 1}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 1}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 3}
{'seq': [1, 2, 3], 'pid': 2}
{'seq': [1, 2, 3], 'pid': 3}
```


## Grids with Custom Functions

If you are using a script for any of your grids that isn't a system installed
module, then (akin to a standard gridtest) it needs to be included under a section 
header that is named by the relevant module, and that includes the filename to import.
For example, the file [script.py](script.py) in the present working directory has
a function, `get_pokemon_id` that I want to use. Here is how I'd write the recipe:

script:
  filename: script.py 
  grids:

    # A grid that will generate 10 random values using a custom function
    generate_pids:
      functions: 
        pid: script.get_pokemon_id
      count: 10
```

Now let's run it!

```bash
$ gridtest gridview grids-with-function.yml generate_pids
{'pid': '125'}
{'pid': '780'}
{'pid': '508'}
{'pid': '566'}
{'pid': '803'}
{'pid': '513'}
{'pid': '854'}
{'pid': '405'}
{'pid': '639'}
{'pid': '353'}
```

Akin to the previous example, since we've provided a function to pass our grid
arguments into, the results are returned. This is how gridtest can
use a grid specified under a test to generate a list of values for an argument.

## Grids with Unwrapped Functions

Let's say that we have some function that produces a list of lists, and we want to
use that as a list to be parameterized for another grid. First, here is our
function:

```python
```

And here is how we might define it in the grid:
```yaml
    # Generate a list of lists, intended to be unwrapped and used for another grid
    unwrapped_grid:
      functions:
        numbers: 
          func: script.generate_numbers
          unwrap: true
```

The "unwrap" serves to take the list, and unwrap it so that each result can
be used separately. If we generate the grid, the parameterization is done with
any additional arguments. Since we don't have any, we get a list of 10 inputs, 
each of length 10

```bash
$ gridtest gridview grids-with-function.yml unwrapped_grid 
{'numbers': [82, 82, 82, 82, 82, 82, 82, 82, 82, 82]}
{'numbers': [52, 52, 52, 52, 52, 52, 52, 52, 52, 52]}
{'numbers': [41, 41, 41, 41, 41, 41, 41, 41, 41, 41]}
{'numbers': [25, 25, 25, 25, 25, 25, 25, 25, 25, 25]}
{'numbers': [43, 43, 43, 43, 43, 43, 43, 43, 43, 43]}
{'numbers': [64, 64, 64, 64, 64, 64, 64, 64, 64, 64]}
{'numbers': [43, 43, 43, 43, 43, 43, 43, 43, 43, 43]}
{'numbers': [39, 39, 39, 39, 39, 39, 39, 39, 39, 39]}
{'numbers': [27, 27, 27, 27, 27, 27, 27, 27, 27, 27]}
{'numbers': [64, 64, 64, 64, 64, 64, 64, 64, 64, 64]}
```

We can also ask to just look at the numbers parameter, as it was defined
before we produced the parameter matrix above (it's a list of lists)

```bash
$ gridtest gridview grids-with-function.yml unwrapped_grid --arg numbers --count
Variable numbers has length 10.

$ gridtest gridview grids-with-function.yml unwrapped_grid --arg numbers
[[32, 32, 32, 32, 32, 32, 32, 32, 32, 32], [70, 70, 70, 70, 70, 70, 70, 70, 70, 70], [75, 75, 75, 75, 75, 75, 75, 75, 75, 75], [86, 86, 86, 86, 86, 86, 86, 86, 86, 86], [86, 86, 86, 86, 86, 86, 86, 86, 86, 86], [50, 50, 50, 50, 50, 50, 50, 50, 50, 50], [92, 92, 92, 92, 92, 92, 92, 92, 92, 92], [1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [80, 80, 80, 80, 80, 80, 80, 80, 80, 80], [35, 35, 35, 35, 35, 35, 35, 35, 35, 35]]
```

We can also add another parameter (e.g., defining two for length doubles the results)

```yaml
    # Generate a list of lists, intended to be unwrapped and used for another grid
    unwrapped_grid:
      args:
        length: [10, 20]
      functions:
        numbers: 
          func: script.generate_numbers
          unwrap: true
```

Regardless, we could then reference this variable with `ref` for another grid.  For
example, here we want to 

```yaml
    sum_unwrapped_numbers:
      ref:
        numbers: unwrapped_grid.numbers
      functions: 
        total: script.dosum
```

And then we generate a grid with the result.

```bash
$ gridtest gridview grids-with-function.yml sum_unwrapped_numbers
{'numbers': [90, 90, 90, 90, 90, 90, 90, 90, 90, 90], 'total': 900}
{'numbers': [72, 72, 72, 72, 72, 72, 72, 72, 72, 72], 'total': 720}
{'numbers': [88, 88, 88, 88, 88, 88, 88, 88, 88, 88], 'total': 880}
{'numbers': [21, 21, 21, 21, 21, 21, 21, 21, 21, 21], 'total': 210}
{'numbers': [37, 37, 37, 37, 37, 37, 37, 37, 37, 37], 'total': 370}
{'numbers': [96, 96, 96, 96, 96, 96, 96, 96, 96, 96], 'total': 960}
{'numbers': [40, 40, 40, 40, 40, 40, 40, 40, 40, 40], 'total': 400}
{'numbers': [37, 37, 37, 37, 37, 37, 37, 37, 37, 37], 'total': 370}
{'numbers': [26, 26, 26, 26, 26, 26, 26, 26, 26, 26], 'total': 260}
{'numbers': [2, 2, 2, 2, 2, 2, 2, 2, 2, 2], 'total': 20}
```

## What to do with Grids?

You might just be using grids inline to go with your [tests](../testing/). However,
grids are useful to many things:

 - create input spaces to test different machine learning models
 - keep parameterizations under version control
 - create a parameter space to export to json to use for a non-Python use case
 - and of course testing to create a parameter space

Next you might want to read about how gridtest grids can be 
[used in tests](https://vsoch.github.io/gridtest/getting-started/metrics/#adding-arguments).
---
title: Metrics
category: Getting Started
permalink: /getting-started/metrics/index.html
order: 6
---

### Metrics

Gridtest includes a suite of decorators that can be added to gridtest.yml
test yaml files in order to measure a benchmark or metric.  These include:

| Name   | Description                            | Usage    |
|--------|----------------------------------------|----------|
| timeit | print time (ms) for function execution | @timeit  |
| result | print function result as a metric | @result  |
| length | calculate length of a result, None if not relevant | @length  |

This namespace of decorators will be looked for in the `gridtest.decorators`
module and you don't need to specify this path. If you define a custom decorator, 
you can simply define the module and function to import (e.g., `@script.mydecorator`).

<a id="an-example-decorator">
### An Example Decorator

As an example, let's take a look at using the (likely familiar) timeit decorator,
which we can find in `gridtest.decorators`:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/optimize/script.py
  tests:
    script.add:
      # metrics are each a decorator, we first look to gridtest.decorators then external import
    - metrics:
        - "@timeit"
      args:
        one: 1.0
        two: 2
      istrue: "isinstance(self.result, float)"
      isfalse: "isinstance(self.result, int)"
```

The recipe above is very simple - we are testing a single function, script.add, and
we've defined one metric, a decorator called "timeit" to measure the total time.
Note that metrics belong on the level of the test, since it's likely we don't want
to use any single decorator across all tests.

<a id="attributes-of-a-metric-decorator">
### Attributes of a Metric Decorator

Decorators can only work together (meaning multiple applied to the same function,
and collected for the same single run) given that:

 1. they don't interfere with the function input and return of output
 2. they don't add significant processing / computational needs
 3. they must print their name (identifier) to stdout on a single line, followed by their result

Given these three criteria are met, we can apply multiple decorators to one
function run, and easily collect output based on parsing the stdout. We can then
even generate a report for the run that shows the different metrics. 

<a id="running-with-a-decorator">
### Running with a Decorator

So let's first try running with a simple metric to record time. You'll notice that the metrics
is printed in a second table to the screen!

```bash
$ gridtest test
[4/4] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        istrue isinstance(self.result, float) isfalse isinstance(self.result, int)
script.add.1                   success                        equals 1+2                    
script.add_with_type.0         success                        returns 3                     
script.add_with_type.1         success                        raises TypeError              

________________________________________________________________________________________________________________________
script.add.0                   @timeit                        0.00 ms                       

4/4 tests passed
```

<a id="adding-arguments">
### Adding Arguments

Great - so we've measured time for one function. What if we want to measure the time for
a function, but across a parameter grid? We might want to adopt our recipe to allow for this:

```yaml
    args:
      one:
        min: 0
        max: 5
        list: [10, 15]
    args:
      two: 2
```

but this won't produce a very interesting result, because the result is just adding
two numbers. Let's try something that will work with out timeit function, namely
write a function that will sleep for some number of input seconds. 
Our goal is to see that the timeout output increases to match the input seconds.
That might look like this:

```python
from time import sleep
def gotosleep(seconds):
    """sleep for whatever specified number of seconds are provided"""
    sleep(seconds)
```

We would then use `gridtest check` to detect the new function:

```bash
$ gridtest check test.yml 

New sections to add:
script.script.gotosleep
```

And add the template to work with.

```bash
$ gridtest update gridtest.yml
Adding function script.gotosleep
Writing [gridtest|test.yml] to test.yml
```

And then fill in the template to add the metrics `timeit` and also define a grid
of parameters to run it over. That might look like this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/optimize/script.py
  script.gotosleep:
  - metrics:
    - '@timeit'
    args:
      seconds:
        list: [10, 15]
        max: 5
        min: 0
```

Before we run the test, let's talk about the formatting of the yaml.
Notice that we've moved the "one" argument up from the "args" section into the grid
section. This tells gridtest that we want to run a grid of tests for some 
parameterization of our argument "one." In the example above, we want to include
a range from 0 to 5 (default increment is by 1) and also add the one off values of 10 and 15 
(provided in list). This gives us the following allowable keys for a grid
parameter:

**min** and **max** are to be used when specifying a range. When unset, **by** 
would be 1. If you want to decrease, set a negative value for by.

**list** is for when you want to include a list of values, even in addition to a
range already specified as in the example above.

> Did you know that an arguments section is just an inline grid? You can also define this section globally and share between tests, see [grids](../grids) to learn more.

<a id="run-the-gridtest">
## Run the GridTest

Now let's run the test!

```bash
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.gotosleep.0             success                                                      
script.gotosleep.1             success                                                      
script.gotosleep.2             success                                                      
script.gotosleep.3             success                                                      
script.gotosleep.4             success                                                      
script.gotosleep.5             success                                                      
script.gotosleep.6             success                                                      
script.add.0                   success                        istrue isinstance(self.result, float) isfalse isinstance(self.result, int)
script.add.1                   success                        equals 1+2                    
script.add_with_type.0         success                        returns 3                     
script.add_with_type.1         success                        raises TypeError              

________________________________________________________________________________________________________________________
script.gotosleep.0             @timeit                        0.01 ms                       
script.gotosleep.1             @timeit                        1000.77 ms                    
script.gotosleep.2             @timeit                        2001.84 ms                    
script.gotosleep.3             @timeit                        3001.72 ms                    
script.gotosleep.4             @timeit                        4003.81 ms                    
script.gotosleep.5             @timeit                        10009.59 ms                   
script.gotosleep.6             @timeit                        15013.49 ms                   

11/11 tests passed
```

Awesome! We've run a grid of tests over the different values for seconds, and have 
reported the total time taken via the timeit decorator. See the [gridtest.yml](gridtest.yml)
for the full test recipe.

A more interactive results view will be developed, along with more real world examples for 
using a decorator, and custom decorator.

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Python
category: Getting Started
permalink: /getting-started/python/index.html
order: 4
---

GridTest can be interacted with from within a python shell. We do this
by way of the GridRunner.

## GridTest Shell

You can start a gridtest shell with a run by doing the following. In the command
below, we don't have a particular test file in mind to load, so we just issue
the `gridtest shell` command:

```bash
$ gridtest shell

Gridtest Interactive Shell
runner = GridRunner('tests.yml')
Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:      
```

Note that it's detected not having a testing file, and showed us how to load
a generic one:

```python
runner = GridRunner('tests.yml')
```

In fact, GridRunner is already loaded in the workspace, so we could just use it
on (an existing) test file.

```python
runner = GridRunner("temp-tests.yml")

runner                                                                  
# [gridtest|temp-tests.yml]
```

If you want to have this preloaded (and the runner available) just provide the
test file back when you originally do `gridtest shell`:

```bash
$ gridtest shell temp-tests.yml

Gridtest Interactive Shell
testfile: temp-tests.yml
  runner: [gridtest|temp-tests.yml]
Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.
```

You'll notice that now the interpreter shows us that we have the runner loaded,
and the testfile defined! Either way, once you are here, the easiest thing to do is just run the tests.

```python
runner.run()
runner.run(parallel=False)
```

And you can look at the function docstrings to see how to customize this command.
At this point we might want to get the tests, which are each of type
GridTest. This is a dictionary with keys corresponding to the test name,
and values the instantiation of a GridTest.

```python
tests = runner.get_tests()                                                      
{'temp.create_directory.0': [test|temp.create_directory],
 'temp.write_file.0': [test|temp.write_file]}
```

Let's say we wanted to interact with the first, we would inspect it as follows:

```python
test = tests['temp.create_directory.0']
test                                                                   
[test|temp.create_directory]
```

You can inspect parameters, run the test, or directly interact with many of the individual
functions for the class:

```python
test.params                                                            
{'args': {'dirname': '/tmp/gridtest-dir.ho18j_29'}}

test.run()
test.result
test.name                                                              
# 'temp.create_directory'
```

Please [open an issue](https://github.com/{{ site.repo }}/issues) if you would like specific help for using the classes
interactively.

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Results
category: Getting Started
permalink: /getting-started/results/index.html
order: 7
---

 - [Export](#export): as a json file for your own usage.
 - [Web Report](#report): an interactive web report to share more easily.


<a id="export">
## Export via Json

GridTest currently provide a simple way to export a json file of results. As an
example, for the [custom-decorator](https://github.com/vsoch/gridtest/tree/master/examples/custom-decorator) example, we can run tests and generate a results.json
file as follows:

```bash
$ gridtest test --save results.json
[9/9] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.multiply_sentence.0     success                                                      
script.multiply_sentence.1     success                                                      
script.multiply_sentence.2     success                                                      
script.multiply_sentence.3     success                                                      
script.multiply_sentence.4     success                                                      
script.multiply_sentence.5     success                                                      
script.multiply_sentence.6     success                                                      
script.multiply_sentence.7     success                                                      
script.multiply_sentence.8     success                                                      

________________________________________________________________________________________________________________________
script.multiply_sentence.0     @script.countwords             5 words                       
script.multiply_sentence.1     @script.countwords             7 words                       
script.multiply_sentence.2     @script.countwords             7 words                       
script.multiply_sentence.3     @script.countwords             21 words                      
script.multiply_sentence.4     @script.countwords             31 words                      
script.multiply_sentence.5     @script.countwords             31 words                      
script.multiply_sentence.6     @script.countwords             41 words                      
script.multiply_sentence.7     @script.countwords             61 words                      
script.multiply_sentence.8     @script.countwords             61 words                      

9/9 tests passed
```

The results file is a list of results, each a dictionary of attributes for one of
the tests (meaning that the tests above would produce a list of nine entities). Here
is an example of one:

```json
[
    {
        "name": "script.multiply_sentence.8",
        "function": "script.multiply_sentence",
        "filename": "/home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py",
        "out": [],
        "err": [],
        "result": "You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.You are my sunshine, my only sunshine.",
        "params": {
            "metrics": [
                "@script.countwords"
            ],
            "grid": {
                "count": {
                    "list": [
                        1,
                        5,
                        10
                    ]
                },
                "sentence": {
                    "list": [
                        "He ran for the hills.",
                        "Skiddery-a rinky dinky dinky, skittery rinky doo.",
                        "You are my sunshine, my only sunshine."
                    ]
                }
            },
            "args": {
                "count": 10,
                "sentence": "You are my sunshine, my only sunshine."
            }
        },
        "raises": null,
        "success": true,
        "metrics": {
            "@script.countwords": [
                "61 words"
            ]
        },
        "module": "script"
    }
]
```

This might be useful, for example, to see that our function isn't putting a space between
string combinations, so we might be counting words incorrectly. Oh no! Thank goodness it's just a
dummy example.

<a id="report">
## Web Reports

GridTest allows for generation of static html reports to go along with tests.

### Run Tests

As an example, we can start with the [interface](https://github.com/vsoch/gridtest/tree/master/examples/interface)
example, which is an extended example of the above. Running tests looks like this:

```bash
$ gridtest test
[12/12] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.multiply_sentence.0     success                                                      
script.multiply_sentence.1     success                                                      
script.multiply_sentence.2     success                                                      
script.multiply_sentence.3     success                                                      
script.multiply_sentence.4     success                                                      
script.multiply_sentence.5     success                                                      
script.multiply_sentence.6     success                                                      
script.multiply_sentence.7     success                                                      
script.multiply_sentence.8     success                                                      
script.unique_sentence.0       success                                                      
script.unique_sentence.1       success                                                      
script.unique_sentence.2       success                                                      

________________________________________________________________________________________________________________________
script.multiply_sentence.0     @script.countwords             5 words                       
script.multiply_sentence.0     @script.countletters           17 letters                    
script.multiply_sentence.1     @script.countwords             7 words                       
script.multiply_sentence.1     @script.countletters           43 letters                    
script.multiply_sentence.2     @script.countwords             7 words                       
script.multiply_sentence.2     @script.countletters           32 letters                    
script.multiply_sentence.3     @script.countwords             21 words                      
script.multiply_sentence.3     @script.countletters           85 letters                    
script.multiply_sentence.4     @script.countwords             31 words                      
script.multiply_sentence.4     @script.countletters           215 letters                   
script.multiply_sentence.5     @script.countwords             31 words                      
script.multiply_sentence.5     @script.countletters           160 letters                   
script.multiply_sentence.6     @script.countwords             41 words                      
script.multiply_sentence.6     @script.countletters           170 letters                   
script.multiply_sentence.7     @script.countwords             61 words                      
script.multiply_sentence.7     @script.countletters           430 letters                   
script.multiply_sentence.8     @script.countwords             61 words                      
script.multiply_sentence.8     @script.countletters           320 letters                   
script.unique_sentence.0       @script.countwords             5 words                       
script.unique_sentence.0       @script.countletters           17 letters                    
script.unique_sentence.1       @script.countwords             6 words                       
script.unique_sentence.1       @script.countletters           38 letters                    
script.unique_sentence.2       @script.countwords             6 words                       
script.unique_sentence.2       @script.countletters           30 letters                    

12/12 tests passed
```

As you can see, when you have multiple tests with more than one metric, the output
to the terminal can get crowded. We can help that with report generation.

### Generate Web Report

Generating the interface is simple! Remember that if we wanted to save the raw
results, we could add `--save` with a filename:

```bash
$ gridtest test --save results.json
```

You can do the same for the interface, but instead provide the --save-web
flag. If you don't provide an argument to this flag, it will generate
a folder for you in `/tmp` with files that can be added to a web server.
If you do provide a folder path, the same folder will be renamed
to be there (so make sure it doesn't exist). The following command
will save pages to a subfolder called "web" in the present working directory,
which should not exist.

```bash
$ gridtest test --save-web web/
```

You'll notice in the folder that we've generated some basic webby files (.html, .js. and .css)
and also the same results.json that will populate the interface.

```bash
$ tree web
web/
├── gridtest.css
├── gridtest.js
├── index.html
└── results.json
```

You should be able to put these static files on GitHub pages, or just cd
into the folder and run a webserver:

```bash
python -m http.server 9999
```

And see the content at `http://localhost:9999`. This
is a fairly simple results template that lets you select functions in the left
columns, and then see specific tests in the right table

![img/gridtest.png](../img/gridtest.png)

and mouse over the test name to see the output, error, and metrics recorded.

![img/detail.png](../img/detail.png)

The report generated above is available for viewing [here]({{ site.baseurl }}/templates/report/)

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Installation
category: Installation
permalink: /install/index.html
order: 1
---

<a id="install">
## Install

GridTest can be installed natively (python 3 recommended) with pip:

```bash
pip install gridtest
```

or via conda-forge:

```bash
conda install --channel conda-forge gridtest
```

or you can clone and install from source:

```bash
$ git clone https://github.com/vsoch/gridtest
$ cd gridtest
$ python setup.py install
```

or

```bash
$ pip install -e .
```

When you have installed GridTest, there will be an executable "gridtest"
placed in your bin folder:

```bash
which gridtest
/home/vanessa/anaconda3/bin/gridtest
```

and you should be able to run the executable and see the usage:

```bash
$ gridtest

GridTest Python v0.0.0
usage: gridtest [-h] [--version] {version,test,generate} ...

Python Grid Testing

optional arguments:
  -h, --help            show this help message and exit
  --version             suppress additional output.

actions:
  actions for gridtest

  {version,test,generate}
                        gridtest actions
    version             show software version
    test                run a grid test.
    generate            generate a grid test yaml file.
```

<a id="running-tests">
## Running Tests

Once you've installed gridtest, you can run the test suite with pytest, and
install a dependency for testing, the pokemon library:

```bash
pip install pokemon
pytest -sv tests/*py
```

The test suite is also run during continuous integration for GitHub actions,
and will run on pull requests if you don't want to run these commands locally.


If you have any questions or issues, please [open an issue]({{ site.repo }}/issues).
---
title: Temporary Tutorial
category: Tutorials
permalink: /tutorials/temp/index.html
order: 3
---

## Temporary Paths

In this tutorial we will show you how to use the functions `{% raw %}{% tmp_path %}{% endraw %}` and `{% raw %}{% tmp_dir %}{% endraw %}` to  create and then reference a temporary file path or directory, respectively.
If you want to see how to generate tests you might see the [basic](basic) tutorial instead.
If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

### Generate Testing Template

Let's say that we have this script, called "temp.py"

```python
# Basic functions for testing custom functions in format {% raw %}{% tmp_path %}{% endraw %}

import os


def write_file(filename):
    """write a file with some random nonsense"""
    with open(filename, "w") as filey:
        filey.write("I heard there was an octupus living in that Christmas tree.")


def create_directory(dirname):
    """create a directory named according to input variable dirname"""
    if not os.path.exists(dirname):
        os.mkdir(dirname)
    return dirname
```

Notice that the first function writes a temporary file that comes from
an input, and the second creates a directory and returns it. For
this gridtest, we will want to test these functions. We thus
generated the initial template for the script temp.py shown above
like this:

```bash
$ gridtest generate temp.py gridtest.yml
Extracting write_file from temp
Extracting create_directory from temp
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the "gridtest.yml" file will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case.

### Customize

You can then open the file in a text editor, and add arguments to each.
If your library uses typing, the typing will be checked at testing time,
and it's not specified here. You'll generally want to fill in args for
each testing condition (or leave null for None). For more detail about
templating, see [the templating documentation](https://vsoch.github.io/gridtest/getting-started/templating). 
We ultimately updated the template to include the following:

```yaml
temp:
  filename: temp.py
  tests:
    temp.create_directory:
    - args:
        dirname: "{% raw %}{% tmp_dir %}{% endraw %}"
      returns: "{% raw %}{{ args.dirname }}{% endraw %}"
    temp.write_file:
    - args:
        filename: "{% raw %}{% tmp_path %}{% endraw %}"
      exists: "{% raw %}{{ args.filename }}{% endraw %}"
```

The above recipe says that we want to test the function `create_directory`
temp.py, and we want gridtest to generate a temporary directory
(the variable `{% raw %}{% tmp_dir %}{% endraw %}` and then test that whatever name is generated
(`{% raw %}{{ args.dirname }}{% endraw %}`) is returned by the function. For the function write_file,
the input "filename" will have a randomly generated temporary file created
for it with `{% raw %}{% tmp_path %}{% endraw %}`, and we will test that it exists by referencing it
with `{% raw %}{{ args.filename }}{% endraw %}`. Both will also be cleaned up at the completion of the test.

### Temp Variables

For each of `tmp_path` and `tmp_dir`, you can optionally define a boolean to cleanup,
or a prefix. That might look like this:

```yaml
"{% raw %}{% tmp_dir cleanup=False %}{% endraw %}"
"{% raw %}{% tmp_path prefix=mytest %}{% endraw %}"
```

See the [function helpers](/gridtest/api/source/gridtest.html#module-gridtest.func) module
for more details.

## Test

Once we have our testing file and the original script, we can run the tests as follows:

```bash
$ gridtest test gridtest.yml 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
temp.create_directory.0        success                        returns /tmp/gridtest-dir._cfykliv
temp.write_file.0              success                        exists /tmp/gridtest-file-ns9yrb5g

2/2 tests passed
```

Or since gridtest.yml is the default, just leave it out to find the file in
the present working directory:

```bash
$ gridtest test
```

And the directory mentioned, `/tmp/gridtest-dir.j0aio0l6` will be cleaned up
upon completion. If we don't want to clean it up, we can add `--no-cleanup`:

```bash
$ gridtest test gridtest.yml --no-cleanup
[2/2] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
temp.create_directory.0        success                        returns /tmp/gridtest-dir.ntgs4pp3
temp.write_file.0              success                        exists /tmp/gridtest-file-x9hw6l12

2/2 tests passed
``` 

And then the directory generated would still exist after the run:

```bash
$ ls -l /tmp/gridtest-dir.ntgs4pp3
total 0
```

The same would be true for the testing file.

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Tutorials
category: Tutorials
permalink: /tutorials/index.html
order: 1
---

## GridTest Examples

The following tutorials are meant to show you basic through more advanced use cases for
using GridTest.

### Generate and Test

 - [Basic Tutorial](basic/): Generating tests for a single Python file with functions
 - [Temp Tutorial](temp/): Create temporary files and folders for testing, optionally without cleanup
 - [Boolean Logic Tutorial](boolean/): Use `istrue` and `isfalse` to check custom boolean logic
 - [Class Tutorial](class/): Writing tests for Python classes
 - [Package Tutorial](package/): Writing tests for a system package, requests
 - [Decorators Tutorial](decorators/) Run a grid of tests for a gridtest or custom decorator.

## Grids

 - [Sample Grid](samplegrid/): creating grids by selecting randomly via separate functions.

## Real World Examples

 - [pySchrodinger](https://github.com/researchapps/pySchrodinger): a walk through of creating a grid to run over a Schrodinger class to generate metrics and results for inspection.
 - [clustering-grids](https://github.com/vsoch/gridtest/tree/master/examples/clustering-grids): an example of taking a scikit-learn tutorial of algorithms and datasets, and splitting into grids to collect metrics for each.

If you want to request or suggest a tutorial, please [open an issue](https://github.com/{{ site.repo }}/issues).
---
title: Boolean Logic Tutorial
category: Tutorials
permalink: /tutorials/boolean/index.html
order: 4
---

## Boolean Logic

If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

### Write Functions

Let's say we start with these functions, and save them to a file called truefalse.py

```python
# These are functions in my script
# Typing is here, so Python 

def add(one, two):
    """add will add two numbers, one and two. There is no typing here"""
    return one + two

def add_with_type(one: int, two: int) -> int:
    """add_with_type will add two numbers, one and two, with typing for ints."""
    return one + two
```

We would first generate a gridtest template like the following:

```bash
$ gridtest generate truefalse.py gridtest.yml
Extracting add from truefalse
Extracting add_with_type from truefalse
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the file "gridtest.yml" will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case.

### Return Types

Since we are primarily interested with the `istrue` and `isfalse` return
types, let's just look at examples for each of those.

**istrue**

istrue is used when we want to check if something is True.
You usually would want to refer to an input or output variable:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    istrue: isinstance({% raw %}{{{ result }}{% endraw %}, int)
```

**isfalse**

or you might want the opposite, isfalse:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    isfalse: not isinstance({% raw %}{{{ result }}{% endraw %}, float)
```

**equals**

Equals is similar to returns, but implies some custom code or logic that
needs to be evaluated first:

```yaml
  script.add:
  - args:
      one: 1
      two: 2
    equals: 1+2
```

### Customize

You can then open the file in a text editor, and add arguments to each.
Since in this tutorial we want to test boolean logic, I'll show you what
my final testing yaml file looks like. I went from this starting template:

```yaml
truefalse:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/is-true-false/truefalse.py
  tests:
    truefalse.add:
    - args:
       one: null
        two: null
    truefalse.add_with_type:
    - args:
        one: null
        two: null
```

to be something more reasonable to test:

```yaml
truefalse:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/is-true-false/truefalse.py
  tests:
    truefalse.add:
    - args:
        one: 1.0
        two: 2
      istrue: "isinstance({% raw %}{{ result }}{% endraw %}, float)"
      isfalse: "isinstance({% raw %}{{ result }}{% endraw %}, int)"
    - args:
        one: 1.0
        two: 2
      equals: 1+2
    truefalse.add_with_type:
    - args:
        one: 1
        two: 2
      returns: 3
    - args:
        one: 1.0
        two: 2
      raises: TypeError
```

Notice that we are using istrue and isfalse conditional checks.
For typing, given that a function uses typing, that will be tested. For example,
the last function "add_with_type" would raise a TypeError if we give it a float.
This is why we have added a test case for it. Finally, the template strings `{% raw %}{{ result }}{% endraw %}`
and `{% raw %}{{ returns }}{% endraw %}` are spceial cases. Returns references what you specify your
test to return, and result specifies the result that is actually returned.
We use "result" above because we didn't do any checks for return values for
the same function.

## Test

Finally, you'll have your test file, and an environment where you want to
test. You can run tests like this:

```bash
$ gridtest test gridtest.yml 
[4/4] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
truefalse.add.0                success                        istrue isinstance(3.0, float) isfalse isinstance(3.0, int)
truefalse.add.1                success                        equals 1+2                    
truefalse.add_with_type.0      success                        returns 3                     
truefalse.add_with_type.1      success                        raises TypeError              

4/4 tests passed
```

Or since gridtest.yml is the default, just leave it out to find the file in
the present working directory:

```bash
$ gridtest test
```

You might next want to browse other [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Basic Tutorial
category: Tutorials
permalink: /tutorials/basic/index.html
order: 2
---

## A Basic Script

This tutorial will walk through starting with a minimal Python script (with functions
and typing) and:

 1. Generate a test template
 2. Customize the template
 3. Run tests.

If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

### Generate Testing Template

Let's say that we have this script, called "script.py"

```python
# These are functions in my script
# Typing is here, so Python 

def add(one, two):
    """add will add two numbers, one and two. There is no typing here"""
    return one + two

def add_with_type(one: int, two: int) -> int:
    """add_with_type will add two numbers, one and two, with typing for ints."""
    return one + two

def hello(name):
    """print hello to a name, with no typing"""
    print(f"hello {name}!")

def hello_with_default(name="Dinosaur"):
    """print hello to a name with a default"""
    print(f"hello {name}!")

def hello_with_type(name: str) -> None:
    """print hello to a name, with typing"""
    print(f"hello {name}!")
```

The first step is to generate your testing template, which needs to be done once,
and will generate a yaml file that you can fill in as a template.

```bash
$ gridtest generate script.py gridtest.yml2
Extracting add from script
Extracting add_with_type from script
Extracting hello from script
Extracting hello_with_default from script
Extracting hello_with_type from script
```

The first argument is the input for the generate command, and this can be
a filename, a folder name (that might contain multiple scripts) or a python
module string (.e.g, requests.get). The second argument is the gridtest
output file that will be produced with your tests. After you finish,
the output template file (gridtest.yml) will have a list of tests that
you can add values for. You can delete sections that aren't relevant, or copy
paste new entries to each list for another testing case. The base template
looks like this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/basic/script.py
  tests:
    script.add:
    - args:
        one: null
        two: null
    script.add_with_type:
    - args:
        one: null
        two: null
    script.hello:
    - args:
        name: null
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: null
```

As you can see, each function is discovered, and arguments (and any defaults)
are provided.

### Customize

You can then open the file in a text editor, and add arguments to each.
If your library uses typing, the typing will be checked at testing time,
and it's not specified here. You'll generally want to fill in args for
each testing condition (or leave null for None). For example, we might want to 
change:

```yaml
  script.add:
    args:
    - one: null
    - two: null
```

to instead be:

```yaml
  script.add:
    args:
    - one: 1
    - two: 2
```

To test adding 1+2. You should look at the complete options to [templates]({{ site.baseurl }}/getting-started/templates/)
for details on...


This means that we can edit our script from this:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/basic/script.py
  tests:
    script.add:
    - args:
        one: null
        two: null
    script.add_with_type:
    - args:
        one: null
        two: null
    script.hello:
    - args:
        name: null
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: null
```

to be something more reasonable to test:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/basic/script.py
  tests:
    script.add:
    - args:
        one: 1
        two: 2
      returns: 3
    - args:
        one: 1
        two: null
      raises: TypeError
    script.add_with_type:
    - args:
        one: 1
        two: 2
      returns: 3
    script.hello:
    - args:
        name: Vanessa
    script.hello_with_default:
    - args:
        name: Dinosaur
    script.hello_with_type:
    - args:
        name: 1
```

For typing, given that a function uses typing, that will be tested. For example,
the last function "hello_with_type" will not be successful.

## Test

Finally, you'll have your test file, and an environment where you want to
test. You can run tests like this:

```bash
$ gridtest test gridtest.yml
[6/6] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        raises TypeError              
script.add_with_type.0         success                        returns 3                     
script.hello.0                 success                                                      
script.hello_with_default.0    success                                                      
script.hello_with_type.0       success                                                      

6/6 tests passed
```

And since gridtest looks for the gridtest.yml, you can just do:

```bash
$ gridtest test
```

### Verbose

You can print a little more output about successes or failures with `--verbose`

```bash
$ gridtest test --verbose 
[6/6] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        Exception: TypeError raised as desired. raises TypeError
script.add_with_type.0         success                        returns 3                     
script.hello.0                 success                        hello Vanessa!                
script.hello_with_default.0    success                        hello Dinosaur!               
script.hello_with_type.0       success                        success key set to false, expected failure.

6/6 tests passed
```

Or you can filter to a regular expression (pattern) to only run a subset of
tests:

```bash
$ gridtest test --pattern script.add
[3/3] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        raises TypeError              
script.add_with_type.0         success                        returns 3                     

3/3 tests passed
```

### Debugging

Now let's say there is an error in a script. Let's randomly raise an exception:

```python
def hello(name):
    """print hello to a name, with no typing"""
    raise Exception('ruhroh')
```

If we run tests again, we see a failure with an unexpected exception:

```bash
$ gridtest test gridtest.yml 
[6/6] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.add.0                   success                        returns 3                     
script.add.1                   success                        raises TypeError              
script.add_with_type.0         success                        returns 3                     
script.hello.0                 failure                        ruhroh Unexpected Exception: Exception.                                                     
script.hello_with_default.0    success                                                      
script.hello_with_type.0       success                                                      

5/6 tests passed
```

### Adding Interactivity

If we add the `--interactive` flag, it's going to allow us to cycle through
*every single test* and press Control+d to jump to the next test:

```python
$ gridtest test gridtest.yml --interactive
[script.add:1/6] |=====|-----------------------------|  16.7% 

Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function add at 0x7fe4c0a44200>
 module: <module 'script' from '/home/vanessa/Desktop/Code/gridtest/examples/basic/script.py'>
   args: {'one': 1, 'two': 2}
returns: 3

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:                                                                     
```

But that's not really what we want - we know the failing test is `script.hello.0`
so let's run the tests, but only this particular test for interactive:

```python
$ gridtest test gridtest.yml --interactive --name script.hello

[script.add:1/6] |=====|-----------------------------|  16.7% False
[script.add:2/6] |===========|-----------------------|  33.3% False
[script.add_with_type:3/6] |=================|-----------------|  50.0% False
[script.hello:4/6] |=======================|-----------|  66.7% True


Gridtest interactive mode! Press Control+D to cycle to next test.

Variables
   func: <function hello at 0x7fc118e8c680>
 module: <module 'script' from '/home/vanessa/Desktop/Code/gridtest/examples/basic/script.py'>
   args: {'name': 'Vanessa'}
returns: None

How to test
passed, error = test_types(func, args, returns)
result = func(**args)

Python 3.7.4 (default, Aug 13 2019, 20:35:49) 
Type 'copyright', 'credits' or 'license' for more information
IPython 7.8.0 -- An enhanced Interactive Python. Type '?' for help.

In [1]:
```

In the example above, we show that the function, module, and arguments for
the test are loaded, and you are shown how to run the tests. For example,
here is how we would inspect arguments and then test typing:

```python
In [1]: args                                                                                                                                 
Out[1]: {'name': 'Vanessa'}

In [2]: passed, error = test_types(func, args, returns)                                                                                      

In [3]: passed                                                                                                                               
Out[3]: True

In [4]: error                                                                                                                                
Out[4]: []
```

And then run the actual test to trigger the full error:

```python
In [6]: result = func(**args)                                                                                                                
hello Vanessa!
---------------------------------------------------------------------------
Exception                                 Traceback (most recent call last)
~/Desktop/Code/gridtest/gridtest/main/helpers.py in <module>
----> 1 result = func(**args)

~/Desktop/Code/gridtest/examples/basic/script.py in hello(name)
     18     """print hello to a name, with no typing"""
     19     print(f"hello {name}!")
---> 20     raise Exception('ruhroh')
     21 
     22 def hello_with_default(name="Dinosaur"):

Exception: ruhroh
```

At this point, you can interact with your arguments or the function to debug further.

### Specifying Names

In the case of a file with multiple tests (the typical case) You can also specify the name of the test you want
to interact with:

```python
$ gridtest test --interactive --name script.add
```

For the above, this would interact with all tests that start with script.add. If you
want to limit to the script.add module, you might want to do:

```python
$ gridtest test --interactive --name script.add.
```

Or a specific indexed text for the module:

```python
$ gridtest test --interactive --name script.add.0
```

You might next want to browse [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Class Tutorial
category: Tutorials
permalink: /tutorials/class/index.html
order: 5
---

## Classes

If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

### Write Functions

Let's say we start with these functions, and save them to a file called car.py

```python
# Example car class to run tests for

valid_colors = ["red", "black", "white", "blue"]

# An example of raising a custom exception
class ColorException(Exception):
    pass

class WheelsException(Exception):
    pass


class Car:

    def __init__(self, wheels:int=4, color:str="red", lights:bool=False):
        """a new car must have an even number of wheels, and be a valid color
        """
        if color not in valid_colors:
            raise ColorException
        if wheels % 2 != 0:
            raise WheelsException
        self.wheels = wheels
        self.color = color
        self.lights = lights

    def __str__(self):
        return (f"[car][wheels:{self.wheels}][color:{self.color}]")
    def __repr__(self):
        return self.__str__()

    def honk(self):
        print("Honk, Honk!")

    @property
    def axels(self) -> int:
        return int(wheels/2)

    def switch_lights(self) -> bool:
        self.lights = not self.lights
```

We can first generate a testing config for preview - if we use gridtest generate without
an output file, it will print to the screen:


```bash
$ gridtest generate car.py
Extracting Car from car
Extracting honk from car
Extracting lights from car

car:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/class/car.py
  tests:
    car.Car:
    - args:
        color: red
        lights: false
        wheels: 4
    car.Car.honk:
    - args:
        self: null
    car.Car.lights:
    - args:
        self: null
```

We can then write to file, and we'll use the default name that gridtest can easily discover.

```bash
$ gridtest generate car.py gridtest.yml
```

Notice that we would want to be able to define an instance of a class to be
used for the subsequent testing functions of the class. To do this, take the 
section you want to use, for example this one:


```yaml
  car.Car:
  - args:
      color: red
      lights: false
      wheels: 4
```

and add an "instance" key to it.


```yaml
  car.Car:
  - instance: thisone
    args:
      color: red
      lights: false
      wheels: 4
```


We now want to reference "thisone" as the instance to
use. Just update the "self" variable in each of the class test cases.

```yaml
car:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/class/car.py
  tests:
    car.Car:
    - args:
        color: notacolor
        lights: false
        wheels: 4
      raises: ColorException
    - instance: thisone
      args:
        color: red
        lights: false
        wheels: 4
      isinstance: car.Car 
    car.Car.honk:
    - args:
        self: "{% raw %}{{ instance.thisone }}{% endraw %}"
    car.Car.switch_lights:
    - args:
        self: "{% raw %}{{ instance.thisone }}{% endraw %}"
```

And then the instance of the Car named as instance "this one" (the second block)
will be used for those tests. This is a very basic usage for a class, and we 
expect more complex cases to be written up when they are determined.

### Testing

Now, we can run tests! Since we've named the testing file `gridtest.yml` we can
just run:

```bash
$ gridtest test
[4/4] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
car.Car.0                      success                        raises ColorException         
car.Car.1                      success                        isinstance Car                
car.Car.honk.0                 success                                                      
car.Car.switch_lights.0        success                                                      

4/4 tests passed
```

You might next want to browse other [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Decorators
category: Tutorials
permalink: /tutorials/decorators/index.html
order: 6
---

## Decorators

If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

### GridTest Decorators

GridTest provides a suite of it's own decorators to measure [metrics]({{ site.baseurl }}/getting-started/metrics/).
For a simple example of using the gridtest timeit decorator, see [this tutorial](https://github.com/vsoch/gridtest/tree/master/examples/metrics). For a custom decorator example, you can see the [tutorial here](https://github.com/vsoch/gridtest/tree/master/examples/custom-decorator/) or keep reading.

### Custom Decorator

If you want to run a gridtest that uses a custom metric, you can easily
do this by defining your own decorator. For example, let's say we have a function
to do some kind of text processing. It takes some number of inputs, and
returns raw text. We would then want to count the unique words in the raw text.
Let's go!

#### Create your Functions

Let's first write our functions. We will write a simple function to take
some text input and parse it (`multiply_sentences`) and a decorator
to run any function that returns a string of text, and count the unique
words (`countwords`)

```python
# These are functions in my script

def multiply_sentence(sentence, count):
    return sentence * count

def countwords(func):
    """this is a simple example of a custom decorator - the idea would be that
       the function we are decorating returns some texty value, and we split
       this value by a blank space and then count the number of tokens (words).
    """
    def counter(*args, **kwargs):
        result = func(*args, **kwargs)
        words = len(set(result.split(' ')))
        print(f"@script.countwords {words} words")
        return result

    return counter
```
An important note about the decorator - it needs to be importable, meaning either
the module is already on your Python path, or it's included somewhere in the library
that you are testing. Also note that in order for gridtest to parse the result, you
need to print something to stdout that is prefixed **exactly** with the name of
the decorator defined under metrics. E.g., if we changed `script.countwords` to just
`countwords` the result wouldn't be properly parsed, because gridtest is looking
for the the first.

### Generate your Template

Let's generate a simple template that we can fill in to include a grid. We can
first preview it:

```bash
$ gridtest generate script.py

script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.countwords:
    - args:
        func: null
    script.multiply_sentence:
    - args:
        count: null
        sentence: null
```

We don't want a test for the decorator, so we will write this to file, and remove it.

```bash
$ gridtest generate script.py gridtest.yml

script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - args:
        count: null
        sentence: null
```

Next, let's better refine our arguments.

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - metrics:
      - '@script.countwords'
      args:
        count: [1, 5, 10]
        sentence:
          - "He ran for the hills."
          - "Skiddery-a rinky dinky dinky, skittery rinky doo."
          - "You are my sunshine, my only sunshine."
```

Just for your FYI - if you had wanted to have some set of arguments shared
between tests, you could have defined them as a named grid:

```yaml
grids:
  script_inputs:
    args:
      count: [1, 5, 10]
      sentence:
        - "He ran for the hills."
        - "Skiddery-a rinky dinky dinky, skittery rinky doo."
        - "You are my sunshine, my only sunshine."
```

and then instead pointed to it for your test:

```yaml
script:
  filename: /home/vanessa/Desktop/Code/gridtest/examples/custom-decorator/script.py
  tests:
    script.multiply_sentence:
    - metrics:
      - '@script.countwords'
      grid: script_inputs
```

By default, a grid is generated at the test creation time. However, if you have 
a grid shared by many functions that you want to calculate once and cache,
just set the cache variable to true:


```yaml
grids:
  script_inputs:
    cache: true
    args:
      count: [1, 5, 10]
      sentence:
        - "He ran for the hills."
        - "Skiddery-a rinky dinky dinky, skittery rinky doo."
        - "You are my sunshine, my only sunshine."
```

Regardless of how we specify our grid (globally or inline) the grid says that
for each sentence under the list of sentences, we will run the function
`multiply_sentence` with counts of 1,5, and 10. This would come down to 3x3 or 9 total tests.


#### Running Tests

Let's run the tests! We should see a count of words for each function.

```bash
$ gridtest test
[9/9] |===================================| 100.0% 
Name                           Status                         Summary                       
________________________________________________________________________________________________________________________
script.multiply_sentence.0     success                                                      
script.multiply_sentence.1     success                                                      
script.multiply_sentence.2     success                                                      
script.multiply_sentence.3     success                                                      
script.multiply_sentence.4     success                                                      
script.multiply_sentence.5     success                                                      
script.multiply_sentence.6     success                                                      
script.multiply_sentence.7     success                                                      
script.multiply_sentence.8     success                                                      

________________________________________________________________________________________________________________________
script.multiply_sentence.0     @script.countwords             5 words                       
script.multiply_sentence.1     @script.countwords             7 words                       
script.multiply_sentence.2     @script.countwords             7 words                       
script.multiply_sentence.3     @script.countwords             21 words                      
script.multiply_sentence.4     @script.countwords             31 words                      
script.multiply_sentence.5     @script.countwords             31 words                      
script.multiply_sentence.6     @script.countwords             41 words                      
script.multiply_sentence.7     @script.countwords             61 words                      
script.multiply_sentence.8     @script.countwords             61 words                      

9/9 tests passed
```

If you don't see the result of the decorator under the summary, likely you didn't write
the stdout print within the function to match the name that you supplied under metrics.
For this example, this should be `@script.countwords` in both cases.
And that's it! We've run a grid of tests over a grid of inputs, and measured some metric with
our custom decorators. 

You might next want to browse other [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Sample Grid Tutorial
category: Tutorials
permalink: /tutorials/samplegrid/index.html
order: 7
---

Gridtest can easily be used to generate random samplies for some number of inputs,
where each input is returned via a function as a list of options to select from.

## Write your Functions

Let's start by writing a set of functions. Each of these will return a list
of attributes that we might want to parameterize. For our first example,
we will generate cohorts with all possible combinations. Let's create functions
to return colors, ages, shapes, and animals.

```python
# Example functions to generate lists to parameterize

import colorsys
import random

def generate_rgb_color():
    """return a randomly selected color
    """
    N = int(random.choice(range(0,100)))
    HSV_tuple = (N*1.0/N, 0.5, 0.5)
    return colorsys.hsv_to_rgb(*HSV_tuple)

def generate_shape():
    return random.choice(["square", "triangle", "circle", "ellipsis", "rectangle", "octagon"])
        

def generate_age():
    return random.choice(range(0,100))

def generate_animal():
    return random.choice(["dog", "cat", "bird", "cow", "chicken"])
```

It's fairly simple - each one spits out a random value. You could also
imagine loading data from some other source.

## Generate the Grids

We next want to generate grids for each! This is fairly simple to do too -
each grid is named based on what it selects, and uses the appropriate function
to be run:

```yaml
script:
  grids:

    # Each grid below generates one randomly selected value
    select_color:
      functions:
        color: script.generate_rgb_color

    select_shape:
      functions:
        shape: script.generate_shape

    select_animal:
      functions:
        animal: script.generate_animal

    select_age:
      functions:
        age: script.generate_age
```

We can also verify that each one generates what we would expect with `gridview`:

```bash
$ gridtest gridview cohort-grids.yml select_color
{'color': (0.5, 0.25, 0.25)}

$ gridtest gridview cohort-grids.yml select_animal
{'animal': 'bird'}

$ gridtest gridview cohort-grids.yml select_age
{'age': 47}

$ gridtest gridview cohort-grids.yml select_shape
{'shape': 'triangle'}
```

But that's not very interesting! Let's instead create a full set of data
for a single sample of a cohort:

```yaml
    generate_cohort:
      functions:
        age: script.generate_age
        animal: script.generate_animal
        shape: script.generate_shape
        color: script.generate_rgb_color
```

and generate it:

```yaml
$ gridtest gridview cohort-grids.yml generate_cohort
{'age': 54, 'animal': 'dog', 'shape': 'ellipsis', 'color': (0.5, 0.25, 0.25)}
```

We could run it again to get a different result, but why not just add a
count for the number of samples that we need?

```yaml
    generate_cohort:
      count: 10
      functions:
        age: script.generate_age
        animal: script.generate_animal
        shape: script.generate_shape
        color: script.generate_rgb_color
```

```bash
{'age': 71, 'animal': 'dog', 'shape': 'triangle', 'color': (0.5, 0.25, 0.25)}
{'age': 18, 'animal': 'bird', 'shape': 'circle', 'color': (0.5, 0.25, 0.25)}
{'age': 17, 'animal': 'bird', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 24, 'animal': 'chicken', 'shape': 'rectangle', 'color': (0.5, 0.25, 0.25)}
{'age': 14, 'animal': 'bird', 'shape': 'octagon', 'color': (0.5, 0.25, 0.25)}
{'age': 82, 'animal': 'cow', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 4, 'animal': 'cat', 'shape': 'rectangle', 'color': (0.5, 0.25, 0.25)}
{'age': 55, 'animal': 'cat', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 48, 'animal': 'bird', 'shape': 'square', 'color': (0.5, 0.25, 0.25)}
{'age': 84, 'animal': 'chicken', 'shape': 'triangle', 'color': (0.5, 0.25, 0.25)}
```

We can see that we have ten results:

```bash
$ gridtest gridview cohort-grids.yml generate_cohort --count
10 argument sets produced.
```

Great! Let's save this to file. We don't need gridtest anymore, we can just keep
our parameters for use.

```bash
$ gridtest gridview cohort-grids.yml generate_cohort --export samples.json
```

The results are saved to json.

```json
[
    {
        "age": 88,
        "animal": "bird",
        "shape": "ellipsis",
        "color": [
            0.5,
            0.25,
            0.25
        ]
    },
    {
        "age": 94,
        "animal": "cat",
        "shape": "rectangle",
        "color": [
            0.5,
            0.25,
            0.25
        ]
    },
...
]
```

## Combination Grids

The above works, but it calls the same functions many times. We would do much
better to, for some number of samples that we need, run each function once with
that number (returning a list of lists) and then unwrap it into a grid. For this
first example, we want all possible combinations of the parameters. Let's again
write our functions, and this time, we return a list:

```python
def generate_rgb_color(N=10):
    return generate_rgb_color() * N

def generate_shapes(N):
    return generate_shape() * N
        
def generate_ages(N=10):
    return generate_age() * N

def generate_animals(N=10):
    return generate_animal() * N
```

We can now use "unwrap" to make sure that each argument is represented on its own.
This way, we can parameterize them each over a larger grid. We can check to see this
in the output:

```bash
$ gridtest gridview grids.yml select_shapes
{'shapes': ['octagon']}
{'shapes': ['circle']}
{'shapes': ['rectangle']}
{'shapes': ['square']}
{'shapes': ['square']}
{'shapes': ['rectangle']}
{'shapes': ['ellipsis']}
{'shapes': ['octagon']}
{'shapes': ['ellipsis']}
{'shapes': ['rectangle']}

$ gridtest gridview grids.yml select_colors
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}
{'colors': (0.5, 0.25, 0.25)}

$ gridtest gridview grids.yml select_ages
{'ages': [59]}
{'ages': [44]}
{'ages': [52]}
{'ages': [78]}
{'ages': [78]}
{'ages': [12]}
{'ages': [82]}
{'ages': [40]}
{'ages': [64]}
{'ages': [66]}

$ gridtest gridview grids.yml select_animals
{'animals': ['chicken']}
{'animals': ['cow']}
{'animals': ['cat']}
{'animals': ['chicken']}
{'animals': ['cow']}
{'animals': ['cow']}
{'animals': ['cat']}
{'animals': ['dog']}
{'animals': ['chicken']}
{'animals': ['cow']}
```

Finally, let's make our grid that will generate many combinations of each. Just
a heads up ... it's every possible combination, so we will have 10 x 10 x 10 x 10...
10,000!

```bash
$ gridtest gridview grids.yml generate_cohort --count
10000 argument sets produced.
```
```bash
$ gridtest gridview grids.yml generate_cohort --count
...
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['octagon'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['circle'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
{'age': [56], 'animal': ['cow'], 'shape': ['ellipsis'], 'color': (0.5, 0.25, 0.25)}
```

And you could save this to file again:

```bash
$ gridtest gridview grids.yml generate_cohort --save combinations.json
```

## Overview

Hopefully you can see the two use cases here - the first is to generate a
specific set of samples, given some input functions that generate values.
The second is to produce all possible combinations, and we take
advantage of unwrapping lists. For this dummy example, we are actually
calling the same function many times for the generation of lists of lists.
However, you could imagine having a single computationally intensive function
or maybe an API call that generates a list that you want to use the items
across another grid. You can see all files for these grids and scripts
in the repository [here](https://github.com/vsoch/gridtest/tree/master/examples/sample-grid).

You might next want to browse other [tutorials]({{ site.baseurl }}/tutorials/) available.
---
title: Package Tutorial
category: Tutorials
permalink: /tutorials/package/index.html
order: 5
---

This tutorial covers writing gridtests for a specific package of interest,
such as the `requests` module that might already be installed in your python
site-packages. If you haven't [installed]({{ site.baseurl }}/install/) gridtest, you should do this first.

# Package Testing

You don't necessarily need to write tests just for local files or modules!
Gridtest also works to generate tests for modules installed with your Python.
This provides a quick example for that. First, generate a testing file:

```bash
gridtest generate requests gridtest.yml --skip-classes
```

Note that we just want to test the simple requests functions, so we are 
skipping classes. If I wanted to include private functions:

```bash
$ gridtest generate --include-private --skip-private requests gridtest.yml
```

In the case of requests, I don't really want to test more of the complex
functionlity, but just the functions that I might need to use. Thus,  
I can then open the file and just delete those chunks. I'd then
update the file to customize it with my tests, meaning that I change this:

```yaml
requests:
  filename: /home/vanessa/anaconda3/lib/python3.7/site-packages/requests/__init__.py
  tests:
    requests.adapters.extract_cookies_to_jar:
    - args:
        jar: null
        request: null
        response: null
    requests.adapters.extract_zipped_paths:
    - args:
        path: null
    requests.adapters.get_auth_from_url:
    - args:
        url: null
    requests.adapters.get_encoding_from_headers:
    - args:
        headers: null
    requests.adapters.prepend_scheme_if_needed:
    - args:
        new_scheme: null
        url: null
    requests.adapters.select_proxy:
    - args:
        proxies: null
        url: null
    requests.adapters.urldefragauth:
    - args:
        url: null
    requests.api.delete:
    - args:
        url: null
    requests.api.get:
    - args:
        params: null
        url: null
    requests.api.head:
    - args:
        url: null
    requests.api.options:
    - args:
        url: null
    requests.api.patch:
    - args:
        data: null
        url: null
    requests.api.post:
    - args:
        data: null
        json: null
        url: null
    requests.api.put:
    - args:
        data: null
        url: null
    requests.api.request:
    - args:
        method: null
        url: null
    requests.auth.extract_cookies_to_jar:
    - args:
        jar: null
        request: null
        response: null
    requests.auth.parse_dict_header:
    - args:
        value: null
    requests.auth.to_native_string:
    - args:
        encoding: null
        string: ascii
    requests.cookies.cookiejar_from_dict:
    - args:
        cookie_dict: null
        cookiejar: true
        overwrite: null
    requests.cookies.create_cookie:
    - args:
        name: null
        value: null
    requests.cookies.extract_cookies_to_jar:
    - args:
        jar: null
        request: null
        response: null
    requests.cookies.get_cookie_header:
    - args:
        jar: null
        request: null
    requests.cookies.merge_cookies:
    - args:
        cookiejar: null
        cookies: null
    requests.cookies.morsel_to_cookie:
    - args:
        morsel: null
    requests.cookies.remove_cookie_by_name:
    - args:
        cookiejar: null
        domain: null
        name: null
        path: null
    requests.cookies.to_native_string:
    - args:
        encoding: null
        string: ascii
    requests.hooks.default_hooks:
    - args: {}
    requests.hooks.dispatch_hook:
    - args:
        hook_data: null
        hooks: null
        key: null
    requests.models.check_header_validity:
    - args:
        header: null
    requests.models.cookiejar_from_dict:
    - args:
        cookie_dict: null
        cookiejar: true
        overwrite: null
    requests.models.default_hooks:
    - args: {}
    requests.models.get_auth_from_url:
    - args:
        url: null
    requests.models.get_cookie_header:
    - args:
        jar: null
        request: null
    requests.models.guess_filename:
    - args:
        obj: null
    requests.models.guess_json_utf:
    - args:
        data: null
    requests.models.iter_slices:
    - args:
        slice_length: null
        string: null
    requests.models.parse_header_links:
    - args:
        value: null
    requests.models.requote_uri:
    - args:
        uri: null
    requests.models.stream_decode_response_unicode:
    - args:
        iterator: null
        r: null
    requests.models.super_len:
    - args:
        o: null
    requests.models.to_key_val_list:
    - args:
        value: null
    requests.models.to_native_string:
    - args:
        encoding: null
        string: ascii
    requests.models.unicode_is_ascii:
    - args:
        u_string: null
    requests.requests.check_compatibility:
    - args:
        chardet_version: null
        urllib3_version: null
    requests.requests.delete:
    - args:
        url: null
    requests.requests.get:
    - args:
        params: null
        url: null
    requests.requests.head:
    - args:
        url: null
    requests.requests.options:
    - args:
        url: null
    requests.requests.patch:
    - args:
        data: null
        url: null
    requests.requests.post:
    - args:
        data: null
        json: null
        url: null
    requests.requests.put:
    - args:
        data: null
        url: null
    requests.requests.request:
    - args:
        method: null
        url: null
    requests.requests.session:
    - args: {}
    requests.sessions.cookiejar_from_dict:
    - args:
        cookie_dict: null
        cookiejar: true
        overwrite: null
    requests.sessions.default_headers:
    - args: {}
    requests.sessions.default_hooks:
    - args: {}
    requests.sessions.dispatch_hook:
    - args:
        hook_data: null
        hooks: null
        key: null
    requests.sessions.extract_cookies_to_jar:
    - args:
        jar: null
        request: null
        response: null
    requests.sessions.get_auth_from_url:
    - args:
        url: null
    requests.sessions.get_environ_proxies:
    - args:
        no_proxy: null
        url: null
    requests.sessions.get_netrc_auth:
    - args:
        raise_errors: null
        url: false
    requests.sessions.merge_cookies:
    - args:
        cookiejar: null
        cookies: null
    requests.sessions.merge_hooks:
    - args:
        dict_class: null
        request_hooks: &id001 !!python/name:collections.OrderedDict ''
        session_hooks: null
    requests.sessions.merge_setting:
    - args:
        dict_class: null
        request_setting: *id001
        session_setting: null
    requests.sessions.requote_uri:
    - args:
        uri: null
...
```

to a much shorter and simpler:

```yaml
requests:
  filename: /home/vanessa/anaconda3/lib/python3.7/site-packages/requests/__init__.py
  tests:
    requests.api.get:
    - args:
        params: null
        url: https://google.com
      isinstance: Response
    requests.api.head:
    - args:
        url: https://google.com
      istrue: "self.result.status_code == 301"
    requests.api.options:
    - args:
        url: https://google.com
```

This recipe also shows a good example of how to check for an instance type.
The string "Response" should be given if I check the `type(result).__name__`
for the result. Also notice the `istrue` statement isn't targeting a `{% raw %}{{ result }}{% endraw %}`
template that can be converted to a string and evaluated, but rather the GridTest
instance directly (self.result) and more specifically, the status_code attribute.

You might next want to browse other [tutorials]({{ site.baseurl }}/tutorials/) available.
# Logo Colors

 - Blue: #1e7f9b
 - Yellow: #f8d901
 - https://vectr.com/tmp/brLUCYVEb/f2Arp8T1zH
# CHANGELOG

This is a manually generated log to track changes to the repository for each release.
Each section should include general headers such as **Implemented enhancements**
and **Merged pull requests**. Critical items to know are:

 - renamed commands
 - deprecated / removed commands
 - changed defaults
 - backward incompatible changes
 - migration guidance
 - changed behaviour

The versions coincide with releases on pip.

## [0.2.x](https://github.com/vsoch/gridtest/tree/master) (0.0.x)
 - code cleanup and JoSS review (0.0.15)
 - adding machine learning examples for grids (0.0.14)
 - refactoring to use Grid class, json-tricks export (0.0.13)
 - adding grid exports, and variables section (0.0.12)
 - adding "tests" level to config files (0.0.11)
 - adding grids section in test file, and function definition for args
 - first alpha release with basic grid test generation (0.0.1)
 - skeleton release (0.0.0)
