# Change Log

## 0.8.1

Released on August 11, 2019.

### Added

* Support even newer versions of cwltool (<=1.0.20190228155703)
* Support for Python 3.7

### Changed

* Ignore packed workflows when sorting for loading order

### Removed

* Support for Python 2
* Support for Python 3.4

## 0.8.0

### Added

* Add `MultipleInputFeatureRequirement` when a step gets a list of inputs (#105; see also #101)
* Convert input and output names with dashes (-) to underscores (\_), so they are valid Python names
* Allow printing of workflows (#86)
* Logging (for debugging)
* Support newer versions of cwltool (#108)

### Changed

* Using booleans to indicate how a workflow is saved is deprecated. Instead, a mode string should be used (e.g., `wf.save('wf.cwl', mode='rel')`) (#87)
* Inline saving of workflows is deprecated. When saving a workflow with `mode='inline'`, the workflow is saved as a packed workflow (#92)
* Make `scatter_method` optional when scattering over a single parameter (#103)

## 0.7.2

### Added

* Allow for list of step outputs/wf inputs as step input (#101)
* CFF files with citation metadata
* Link between a step's python name (i.e. how it is called on the WorkflowGenerator object) to step names (#100)
* Allow setting workflow labels
* Allow setting a label for a workflow input
* support for CommandInputEnumSchema as workflow input (#99)
* User manual in documentation

### Changed

* Ensure workflows without a requirements section are loaded into the steps library
* Raise real warning when duplicate cwl step (i.e. having the same file name as another step) is loaded

### Removed

* Method to convert a string to cwl file name

## 0.7.1

### Added

* Load tools before workflows when a working directory is used (#94)
* Make sure no duplicate workflow input names are used (#96)

### Changed

* Inputs with a default value are also recognized as optional

## 0.7.0

### Added

* Save packed workflows
* Save workflows using a working directory (a solution to the problem of dealing with paths to steps if steps are loaded from different local directories)

### Changed

* Prevent name clashes of embedded (sub)workflows (however, this doesn't work when a (sub)workflow is added multiple times)
* Use name of step in workflow to create unique ids when saving steps inline (#82)
* Allow saving workflows with inline steps for step files without shebang (#83)
* Document feature for adding documentation to a workflow (#81)
* Fix saving of relative paths for workflows with steps from urls
* By default, workflows are saved with absolute paths

## 0.6.0

### Added

* Make `WorkflowGenerator` into a context manager (#24)
* Type checking of input and output types (#22)
* Allow saving workflow with inline steps (#38)
* Allow saving workflow with relative paths (#25)
* Documentation on Read the Docs (#35)
* Allow loading of multiple CWL steps (from file, http url, and directory) at the same time

### Changed

* Rename `wf.add_inputs()` to `wf.add_input()` (#11)

### Removed

* Python 3.3 support (Sphinx needs Python 3.4)

## 0.5.1

### Added

* Allow addition of default values to workflow inputs (#32)
* List of steps and workflows in steps library is ordered alphabetically

## 0.5.0

### Added

* Python 3 compatibility (and testing with [tox](https://tox.readthedocs.io/en/latest/))

## 0.4.0

### Added

* Generate unique names for steps that are added to the workflow more than once (#31)
* Pass all outputs from a step, instead of just one (#27)
* Improve listing of workflow steps

## 0.3.1

### Added

* Load ExpressionTools as steps

### Changed

* Preserve the order in which steps were added to the workflow

## 0.3.0

### Changed

* Replace pyyaml by ruamel (fixes compatibility with cwltool)

## 0.2.0

### Added

* Documentation for WorkflowGenerator and Step (#15).
* Allow step to be scattered (#17)
* Tests (#9)
* Shebang to saved CWL file (#14)
* Preprocess shortcuts in CWL steps (#12)
* Allow workflows to be used as steps (subworkflows) (#4)
* Take into account optional arguments (#6)

### Removed

* Step.get_input() because it was not used (#21)

## 0.1.0

### Added

* WorkflowGenerator object that allows users to create CWL workflows. The WorkflowGenerator has functionality to
  * load CWL steps from a directory,
  * list available CWL steps
  * connect the inputs and outputs of CWL steps,
  * determine the types of a step's inputs
  * specify a workflow's inputs and outputs, and
  * add workflow documentation.
scriptcwl
=========

|codacy_grade| |codacy_coverage| |travis| |documentation| |pypi_version| |pypi_supported| |zenodo|

scriptcwl is a Python package for creating workflows in
`Common Workflow Language (CWL) <http://www.commonwl.org/>`_. If you give it a number of CWL
``CommandLineTools``, you can create a workflow by writing a Python script. This can
be done interactively using `Jupyter Notebooks <http://jupyter.org/>`_. The full
documentation can be found on `Read the Docs <http://scriptcwl.readthedocs.io/en/latest/>`_.

.. image:: docs/images/add-multiply-example-workflow.png
   :alt: add multiply example workflow
   :align: center

Given CWL ``CommandLineTools`` for ``add`` and ``multiply`` (these are available
in `scriptcwl <https://github.com/NLeSC/scriptcwl/tree/master/scriptcwl/examples>`_),
a CWL specification of this workflow can be written as:

.. code-block:: python

  from scriptcwl import WorkflowGenerator

  with WorkflowGenerator() as wf:
    wf.load(steps_dir='/path_to_scriptcwl/scriptcwl/examples/')

    num1 = wf.add_input(num1='int')
    num2 = wf.add_input(num2='int')

    answer1 = wf.add(x=num1, y=num2)
    answer2 = wf.multiply(x=answer1, y=num2)

    wf.add_outputs(final_answer=answer2)

    wf.save('add_multiply_example_workflow.cwl')

The workflow has two integers as inputs (``num1`` and ``num2``), and first adds
these two numbers (``wf.add(x=num1, y=num2)``), and then multiplies the answer
with the second input (``num2``). The result of that processing step is the output
of the workflow. Finally, the workflow is saved to a file. The result looks like:

.. code-block:: sh

  #!/usr/bin/env cwl-runner
  cwlVersion: v1.0
  class: Workflow
  inputs:
    num1: int
    num2: int
  outputs:
    final_answer:
      type: int
      outputSource: multiply/answer
  steps:
    add:
      run: add.cwl
      in:
        y: num2
        x: num1
      out:
      - answer
    multiply:
      run: multiply.cwl
      in:
        y: num2
        x: add/answer
      out:
      - answer

The Python and CWL files used in the example can be found in the `examples folder <https://github.com/NLeSC/scriptcwl/tree/master/scriptcwl/examples>`_.

Installation
############

Install using pip:

.. code-block:: sh

  pip install scriptcwl


For development:

.. code-block:: sh

  git clone git@github.com:NLeSC/scriptcwl.git
  cd scriptcwl
  python setup.py develop

Run tests (including coverage) with:

.. code-block:: sh

  python setup.py test

Useful tools
############

To use scriptcwl for creating CWL workflows, you need CWL ``CommandLineTools``.
There are some software packages that help with generating those
for existing command line tools written in Python:

* `argparse2tool <https://github.com/erasche/argparse2tool#cwl-specific-functionality>`_: Generate CWL CommandLineTool wrappers (and/or Galaxy tool descriptions) from Python programs that use argparse. Also supports the `click <http://click.pocoo.org>`_ argument parser.
* `pypi2cwl <https://github.com/common-workflow-language/pypi2cwl>`_: Automatically run argparse2cwl on any package in PyPi.
* `python-cwlgen <https://github.com/common-workflow-language/python-cwlgen>`_: Generate CommandLineTool and DockerRequirement programmatically

License
#######

Copyright (c) 2016-2018, Netherlands eScience Center, University of Twente

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

.. |codacy_grade| image:: https://api.codacy.com/project/badge/Grade/8f383bca18384d8187c10c27affa9d53
                     :target: https://www.codacy.com/app/jvdzwaan/scriptcwl?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/scriptcwl&amp;utm_campaign=Badge_Grade

.. |codacy_coverage| image:: https://api.codacy.com/project/badge/Coverage/8f383bca18384d8187c10c27affa9d53
                       :target: https://www.codacy.com/app/jvdzwaan/scriptcwl?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/scriptcwl&amp;utm_campaign=Badge_Coverage

.. |travis| image:: https://travis-ci.org/NLeSC/scriptcwl.svg?branch=master
              :target: https://travis-ci.org/NLeSC/scriptcwl

.. |documentation| image:: https://readthedocs.org/projects/scriptcwl/badge/?version=latest
                    :target: http://scriptcwl.readthedocs.io/en/latest/?badge=latest

.. |pypi_version| image:: https://badge.fury.io/py/scriptcwl.svg
                    :target: https://badge.fury.io/py/scriptcwl

.. |pypi_supported| image:: https://img.shields.io/pypi/pyversions/scriptcwl.svg
                      :target: https://pypi.python.org/pypi/scriptcwl

.. |zenodo| image:: https://zenodo.org/badge/70679474.svg
                      :target: https://zenodo.org/badge/latestdoi/70679474
                      :alt: DOI
User Manual
===========

.. toctree::
	:maxdepth: 2

	loading_steps
	listing_steps
	workflow_inputs
	setting_documentation
	adding_workflow_steps
	adding_outputs
	printing_workflows
	saving_workflows
	enable_logging
.. _examples:

Example workflows
=================

.. toctree::
	:maxdepth: 1

	nlppln_anonymize
CWL Tips and Tricks
===================

Have a look at the `CWL User Guide: Recommended Practices
<http://www.commonwl.org/user_guide/rec-practices/>`_.

Generate yaml file with workflow inputs
#######################################

You can use ``cwltool --make-template`` to generate a yaml file with all the workflow inputs:
::

	cwltool --make-template add_multiply_example.cwl > inputs.yml

``inputs.yml`` contains:
::

	num1: 0
	num2: 0

Use your favorite text editor to set the inputs to appropriate values. Save the
file, and use it as input for your workflow:
::

	cwltool add_multiply_example.cwl inputs.yml

Using cwl-runner
################

Install the ``cwlref-runner`` package to set ``cwl-runner`` to ``cwltool``:
::

 	pip install cwlref-runner

If ``cwl-runner`` is set, you can run workflows by typing:
::

	chmod +x workflow.cwl
	./workflow.cwl <arguments>

If you have other CWL implementations installed and want ``cwl-runner`` to use one
of these implementations, you should define a symlink that points to the implementation
you want to use; e.g., by manually creating a symlink and adding it to your ``$PATH``
variable, or by using the linux `alternatives <https://linux.die.net/man/8/update-alternatives>`_ system.
Adding workflow steps
=====================

After loading steps and adding workflow inputs, the steps of the workflow should
be specified. To add a step to a workflow, its method must
be called on the ``WorkflowGenerator`` object. For example, to add a step
called ``add`` [#]_ to the workflow, the following method must be called:
::

  answer1 = wf.add(x=num1, y=num2)

The method expects a list of ``key=value`` pairs as input parameters. (To find
out what inputs a step needs call ``wf.inputs(<step name>)``. This method prints
all inputs and their types.) ``wf.<step name>()`` returns a string if the step has
a single output and a tuple of strings if the step has multiple output parameters:
::

  output1, output2 = wf.someStep(input=input)

The order of the outputs is the same as in the step specification, and can be
determined by printing the step signatures using ``print(wf.list_steps())``.

The strings returned by ``wf.<step name>()`` contain output
names that can be used as input for later steps, or that can be connected
to workflow outputs. For example, in a later step, ``answer1`` can be used as input:
::

  answer2 = wf.multiply(x=answer1, y=num2)

Scattering steps
################

Scriptcwl supports `scattering steps <http://www.commonwl.org/v1.0/Workflow.html#WorkflowStep>`_.
To scatter a step, keyword arguments
``scatter`` and ``scatter_method`` must be provided when a step is added to the
workflow. To scatter a step called ``echo``, which has a single input argument
``message``, this would look like:
::

	output = wf.echo(message=input1, scatter='message', scatter_method='dotproduct')

The type of ``message``, should be array (e.g., an array of strings).

To scatter over multiple variables, ``scatter`` also accepts a list of input names:
::

	output = wf.echo(message1=input1, message2=input2, scatter=['message1', 'message2'], scatter_method='dotproduct')

.. [#] Scriptcwl contains two example command line tools, ``add`` and ``multiply``. The Python and CWL files can be found in the `examples folder <https://github.com/NLeSC/scriptcwl/tree/master/scriptcwl/examples>`_.
Useful tools
============

To use scriptcwl for creating CWL workflows, you need CWL ``CommandLineTools``.
There are some software packages that help with generating those.

* `argparse2tool <https://github.com/erasche/argparse2tool#cwl-specific-functionality>`_: Generate CWL ``CommandLineTool`` wrappers (and/or Galaxy tool descriptions) from Python programs that use argparse. Also supports the `click <http://click.pocoo.org>`_ argument parser
* `pypi2cwl <https://github.com/common-workflow-language/pypi2cwl>`_: Automatically run argparse2cwl on any package in PyPi
* `python-cwlgen <https://github.com/common-workflow-language/python-cwlgen>`_: Generate CommandLineTool and DockerRequirement programmatically
Specifying workflow outputs
===========================

When all steps of the workflow have been added, the user can specify
workflow outputs by calling ``wf.add_outputs()``:
::

  wf.add_outputs(final_answer=answer2)
Remove named entities from a directory of text files
====================================================

In this example, we create a pipeline that replaces named entities in a collection
of (Dutch) text documents.
Named entities are objects in text referred to by proper names, such as persons,
organizations, and locations. In the workflow, named entities will be
replaced with their named entity type (i.e., PER (person), ORG (organization),
LOC (location), or UNSP (unspecified)).
The workflow can be used as part of a data anonimization procedure.

The workflow consists of the following steps:

* Extract named entities from text documents using an existing tool called `frog <https://languagemachines.github.io/frog/>`_
* Convert frog output to `SAF, a generic representation for text data <https://github.com/vanatteveldt/saf>`_
* Aggregate data about named entities that occur in the text files
* Replace named entities with their named entity type in the SAF documents
* Convert SAF documents to text

All steps required for this workflow are available through `nlppln <https://github.com/nlppln/nlppln>`_.

Workflow
########

.. image:: images/nlppln-anonymize-workflow.png
  :alt: add multiply example workflow
  :align: center

Scriptcwl script
################

::

	from scriptcwl import WorkflowGenerator

	with WorkflowGenerator() as wf:
	  wf.load(steps_dir='/path/to/dir/with/cwl/steps/')

	  doc = """Workflow that replaces named entities in text files.

	Input:
	  txt_dir: directory containing text files

	Output:
	  ner_stats: csv-file containing statistics about named entities in the text files
	  txt: text files with named enities replaced
	"""
	  wf.set_documentation(doc)

	  txt_dir = wf.add_inputs(txt_dir='Directory')

	  frogout = wf.frog_dir(in_files=txt_dir)
	  saf = wf.frog_to_saf(in_files=frogout)
	  ner_stats = wf.save_ner_data(in_files=saf)
	  new_saf = wf.replace_ner(metadata=ner_stats, in_files=saf)
	  txt = wf.saf_to_txt(in_files=new_saf)

	  wf.add_outputs(ner_stats=ner_stats, txt=txt)

	  wf.save('anonymize.cwl')


CWL workflow
############

::

	cwlVersion: v1.0
	class: Workflow
	inputs:
	  txt-dir: Directory
	  mode: string?

	outputs:
	  ner_stats:
	    type: File
	    outputSource: save-ner-data/ner_statistics

	  out_files:
	    type:
	      type: array
	      items: File
	    outputSource: saf-to-txt/out_files

	steps:
	  frog-ner:
	    run: frog-dir.cwl
	    in:
	      dir_in: txt-dir
	    out: [frogout]

	  frog-to-saf:
	    run: frog-to-saf.cwl
	    in:
	      in_files: frog-ner/frogout
	    out: [saf]

	  save-ner-data:
	    run: save-ner-data.cwl
	    in:
	      in_files: frog-to-saf/saf
	    out: [ner_statistics]

	  replace-ner:
	    run: replace-ner.cwl
	    in:
	      metadata: save-ner-data/ner_statistics
	      in_files: frog-to-saf/saf
	      mode: mode
	    out: [out_files]

	  saf-to-txt:
	    run: saf-to-txt.cwl
	    in:
	      in_files: replace-ner/out_files
	out: [out_files]
Enable logging for debugging
============================

If you get errors while creating workflows, and scriptcwl doesn't give you a
proper error message, you might want to enable logging to try and figure out
what goes wrong.

To enable logging, do:

::

  import logging
  logging.basicConfig(format="%(asctime)s [%(process)d] %(levelname)-8s "
                      "%(name)s,%(lineno)s\t%(message)s")
  logging.getLogger().setLevel('DEBUG')
Installation
============

* pip

    Install using pip:

    .. code-block:: sh

        pip install scriptcwl

    For development:

    .. code-block:: sh

        git clone git@github.com:NLeSC/scriptcwl.git
        cd scriptcwl
        python setup.py develop

    Run tests (including coverage) with:

   .. code-block:: sh

        python setup.py test

* conda
* Windows issues
* for development
Loading steps
=============

Before you can create workflows with scriptcwl, you need to load processing steps
(i.e., ``CommandLineTools``, ``ExpressionTools`` and/or (sub) ``Workflows``).
To load a directory of .cwl files, type:
::

	from scriptcwl import WorkflowGenerator

	with WorkflowGenerator() as wf:
		wf.load(steps_dir='/path/to/dir/with/cwl/steps/')

To load a single cwl file, do:
::

	with WorkflowGenerator() as wf:
		wf.load(step_file='/path/to/workflow.cwl')

The path to the ``step_file`` can be a local file path or a url.

You can also load a list of step files and directories:
::

	al_my_steps = ['step.cwl', 'url.cwl', '/path/to/directory/']
	with WorkflowGenerator() as wf:
		wf.load(step_list=all_my_steps)

``wf.load()`` can be called multiple times. Step files are added to the
steps library one after the other. For every step that is added to the
steps library, a method with the same name is added to the
WorkflowGenerator object. To add a step to the workflow, this method must
be called (examples below).
Listing steps
=============

Steps that are loaded into the WorkflowGenerator's steps library can be listed by running:
::

  print(wf.list_steps())

For the example workflow, the output would be:
::

  Steps
  add...................... answer = wf.add(x, y)
  multiply................. answer = wf.multiply(x, y)

  Workflows

This means that there are two processing steps and no (sub)workflows loaded into the
steps library. The listing contains the complete command to add the step to the workflow
(e.g., ``answer = wf.add(x, y)``). The command is supplied for convenient copy/pasting.
.. scriptcwl documentation master file, created by
   sphinx-quickstart on Mon Nov 13 15:12:14 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the Scriptcwl Documentation!
=======================================

Scriptcwl is a Python package for creating `Common Workflow Language (CWL) <http://www.commonwl.org/>`_ workflows.

.. image:: images/add-multiply-example-workflow.png
     :alt: add multiply example workflow
     :align: center

Given CWL ``CommandLineTools`` for ``add`` and ``multiply`` (these are available
in `scriptcwl <https://github.com/NLeSC/scriptcwl/tree/master/scriptcwl/examples>`_),
a CWL specification of this workflow can be written as:
::

  from scriptcwl import WorkflowGenerator

  with WorkflowGenerator() as wf:
      wf.load(steps_dir='/path_to_scriptcwl/scriptcwl/examples/')

      num1 = wf.add_input(num1='int')
      num2 = wf.add_input(num2='int')

      answer1 = wf.add(x=num1, y=num2)
      answer2 = wf.multiply(x=answer1, y=num2)

      wf.add_outputs(final_answer=answer2)

      wf.save('add_multiply_example_workflow.cwl')

The workflow has two integers as inputs (``num1`` and ``num2``), and first adds
these two numbers (``wf.add(x=num1, y=num2)``), and then multiplies the answer
with the second input (``num2``). The result of that processing step is the output
of the workflow. Finally, the workflow is saved to a file. The result looks like:

.. code-block:: none

  #!/usr/bin/env cwl-runner
  cwlVersion: v1.0
  class: Workflow
  inputs:
    num1: int
    num2: int
  outputs:
    final_answer:
      type: int
      outputSource: multiply/answer
  steps:
    add:
      run: add.cwl
      in:
        y: num2
        x: num1
      out:
      - answer
    multiply:
      run: multiply.cwl
      in:
        y: num2
        x: add/answer
      out:
      - answer

More examples of workflows created using scriptcwl can be found under :ref:`examples`.

Contents
========

.. toctree::
   :maxdepth: 3

   user_manual
   installation
   examples
   useful_tools
   cwl_tips_tricks

API Reference
=============

.. toctree::
   :maxdepth: 2

   scriptcwl <apidocs/modules.rst>
Printing workflows
==================

To view its contents, a workflow can be printed at any time:

.. code-block:: python

  with scriptcwl.WorkflowGenerator() as wf:
    print(wf)

For an empty workflow, this looks like:

.. code-block:: none

  #!/usr/bin/env cwl-runner
  cwlVersion: v1.0
  class: Workflow
  inputs: {}
  outputs: {}
  steps: {}

In a printed workflow, steps are referred to by their absolute paths.
**Therefore, do not use this method for saving workflows.
The absolute paths make them unportable.**
Adding workflow documentation
==============================

To add documentation to your workflow, use the ``set_documentation()`` method:
::

	doc = """Workflow that performs a special calculation with two numbers

	The two numbers are added and the answer is multiplied by the second number.

	Input:
		num1: int
		num2: int

	Output:
		answer: int
	"""
	wf.set_documentation(doc)

Setting labels
##############

Instead of or in addition to documentation, it is also possible to set a label
for a workflow:
::

	wf.set_label('Calculate special number')
Workflow inputs
===============

Wokflow inputs can be added by calling ``add_input()``:
::

	num1 = wf.add_input(num1='int')
	num2 = wf.add_input(num2='int')

The ``add_input()`` method expects a ``name=type`` pair as input parameter.
The pair connects an input name (``num1`` in the example) to a CWL type
(``'int'``). An overview of CWL types can be found in the
`specification <http://www.commonwl.org/v1.0/Workflow.html#CWLType>`_.

Optional inputs
###############

Workflow inputs can be made optional by adding a questionmark to the type:
::

	num1 = wf.add_input(num1='int?')

Default values
##############

When adding an input parameter to a workflow, you can set a default value:
::

	num1 = wf.add_input(num1='int', default=5)

As a consequence, ``default`` cannot be used as a name for a workflow input parameter.

Labels
######

You can also add a label to a workflow input:
::

	num1 = wf.add_input(num1='int', label='The first number that is processed.')

Again, this means ``label`` cannot be used as a name for a workflow input parameter.

Arrays and other complex input types
####################################

Arrays of workflow inputs can be specified with ``[]``:
::

  numbers = wf.add_input(numbers='int[]')

You can also specify the input using a dictionary with two keys: ``{'type':
'array', 'items': 'int'}``.
::

  numbers = wf.add_input(numbers=dict(type='array', items='int'))

This way you also can specify more complex inputs. For example, to create an
array of arrays of strings, do:
::

  inp = dict(type='array', items=dict(type='array', items='string'))
  strings = wf.add_input(my_array_of_array_of_strings=inp)

Use ``print(wf)`` and ``wf.validate()`` to make sure your inputs are correct.

Enums
#####

To use an enum as a workflow input, do:
::

	mode = wf.add_input(mode='enum', symbols=['one', 'two', 'three'])

The ``symbols`` should be a list of strings (lists containing other types are
converted lists of to strings).
Again, ``symbols`` cannot be used as a name for a workflow input parameter.
Saving workflows
================

To save a workflow call the ``WorkflowGenerator.save()`` method:
::

  wf.save('workflow.cwl')

By default, the paths in the ``run`` field of workflow steps are absolute. This means
that a workflow created on one machine cannot be run on another machine. However,
there are multiple options for creating portable workflows.

Saving workflows with relative paths
####################################

To get relative paths in the ``run`` field of workflow steps, use ``mode='rel'``:
::

  wf.save('workflow.cwl', mode='rel')

The paths in the ``run`` field are relative to where the workflow is saved. This
option is convenient when you are creating workflows using a single directory
with possible workflow steps.

Using a working directory
#########################

If you have multiple directories containing workflow steps and the locations of
these directories may differ depending on where software is installed (for example,
if you want to use the generic NLP steps from nlppln, but also need project specific
data processing steps), it is possible to specify a working directory when creating
the ``WorkflowGenerator`` object. If you this, all steps are copied to the working
directory. When you save the workflow using ``mode='wd'``, the paths in the ``run``
fields are set to the basename of the step (because all steps are in the same
directory).
::

  from scriptcwl import WorkflowGenerator

  with WorkflowGenerator(working_dir='path/to/working_dir') as wf:
    wf.load(steps_dir='some/path/')
    wf.load(steps_dir='some/other/path/')

    # add inputs, steps and outputs

    wf.save('workflow', mode='wd')

The workflow is saved in the working directory and then copied to
the specified location. To be able to run the workflow, use the copy in the
working directory (please note that the working directory is not deleted automatically).

Also, steps from urls are not copied to the working directory.

Pack workflows
##############

Another way to create workflows with all steps in one file is to save it with ``mode='pack'``:
::

  wf.save('workflow.cwl', mode='pack')

Please note that packed workflows cannot be used as a building block in ``scriptcwl``.
If you try to load a packed workflow, you will get a warning.

Saved With ``mode='pack'``, the example workflow looks like:
::

  {
      "cwlVersion": "v1.0",
      "$graph": [
          {
              "class": "CommandLineTool",
              "baseCommand": [
                  "python",
                  "-m",
                  "scriptcwl.examples.add"
              ],
              "inputs": [
                  {
                      "type": "int",
                      "inputBinding": {
                          "position": 1
                      },
                      "id": "#add.cwl/x"
                  },
                  {
                      "type": "int",
                      "inputBinding": {
                          "position": 2
                      },
                      "id": "#add.cwl/y"
                  }
              ],
              "stdout": "cwl.output.json",
              "outputs": [
                  {
                      "type": "int",
                      "id": "#add.cwl/answer"
                  }
              ],
              "id": "#add.cwl"
          },
          {
              "class": "CommandLineTool",
              "baseCommand": [
                  "python",
                  "-m",
                  "scriptcwl.examples.multiply"
              ],
              "inputs": [
                  {
                      "type": "int",
                      "inputBinding": {
                          "position": 1
                      },
                      "id": "#multiply.cwl/x"
                  },
                  {
                      "type": "int",
                      "inputBinding": {
                          "position": 2
                      },
                      "id": "#multiply.cwl/y"
                  }
              ],
              "stdout": "cwl.output.json",
              "outputs": [
                  {
                      "type": "int",
                      "id": "#multiply.cwl/answer"
                  }
              ],
              "id": "#multiply.cwl"
          },
          {
              "class": "Workflow",
              "inputs": [
                  {
                      "type": "int",
                      "id": "#main/num1"
                  },
                  {
                      "type": "int",
                      "id": "#main/num2"
                  }
              ],
              "outputs": [
                  {
                      "type": "int",
                      "outputSource": "#main/multiply-1/answer",
                      "id": "#main/final_answer"
                  }
              ],
              "steps": [
                  {
                      "run": "#add.cwl",
                      "in": [
                          {
                              "source": "#main/num1",
                              "id": "#main/add-1/x"
                          },
                          {
                              "source": "#main/num2",
                              "id": "#main/add-1/y"
                          }
                      ],
                      "out": [
                          "#main/add-1/answer"
                      ],
                      "id": "#main/add-1"
                  },
                  {
                      "run": "#multiply.cwl",
                      "in": [
                          {
                              "source": "#main/add-1/answer",
                              "id": "#main/multiply-1/x"
                          },
                          {
                              "source": "#main/num2",
                              "id": "#main/multiply-1/y"
                          }
                      ],
                      "out": [
                          "#main/multiply-1/answer"
                      ],
                      "id": "#main/multiply-1"
                  }
              ],
              "id": "#main"
          }
      ]
  }

Workflow validation
###################

Before the workflow is saved, it is validated using ``cwltool``. Validation can also be
triggered manually:
::

	wf.validate()

It is also possible to disable workflow validation on save:
::

  wf.save('workflow.cwl', validate=False)

File encoding
#############

By default, the encoding used to save workflows is ``utf-8``. If necessary,
a different encoding can be specified:
::

  wf.save('workflow.cwl', encoding='utf-16')
