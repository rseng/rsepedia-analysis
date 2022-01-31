# Change Log

## 0.3.3

### Added

* Command to flatten a list of lists of files (`flatten-list.cwl`)
* Command to merge yaml files using [yamlreader](https://github.com/ImmobilienScout24/yamlreader) (`merge-yaml.cwl`)

### Changed

* By default, keep directory structure of extracted archive (`archive2dir.cwl`)
* Change input for `check-utf8.cwl` from single file to directory of files
* Change input for `delete-empty-files.cwl` from directory to list of files
* Update method for saving workflows (requires scriptcwl >= 0.8.0)

## 0.3.2

### Added

- Command to ls a chunk of files (`ls_chunk.cwl`) This can be used to run workflows for a number of files in a directory (e.g., if running it for all files takes too long)
- Command to create a chunked file list (`create_chunked_list.cwl`)
- Command to return a directory containing all the files from a list of directories (`flatten-dirs.cwl`)
- Command to return a directory containing as subdirectories the list of input directories (`gather-dirs.cwl`)

### Changed

- Expression tool `save-files-to-dir.cwl` extracts the name of the directory to save to from the input files if the `dir_name` is omitted

## 0.3.1

### Added

- Command to extract tar.gz files (`tar.cwl`)
- Command to create directories (`mkdir.cwl`)
- Command to remove XML elements from XML files (`remove-xml-elements.cwl`)
- Command to align strings using edlib (in separate repository)
- Command to parse text using pattern (in separate repository)
- Command to remove newlines from a text (`remove_newlines.cwl`)

### Changed

- Remove erroneous arguments section from `freqs.cwl`
- Add option to filter on filename to `ls.cwl`

### Removed

- Dependency on `pyjq`.

## 0.3.0

### Added

- CWL files for NLP functionality (so they don't have to be downloaded separately)
- Dockerfile to run nlppln on Windows
- First tests
- Python 3 support
- By default, save workflows using working directory
- Documentation on [Read the Docs](http://nlppln.readthedocs.io/en/latest/)
- Command to copy and rename files (`copy-and-rename-files.cwl`)
- Command to generate data for TextDNA visualization (`textDNA-generate.cwl`)
- Command to normalize whitespace and punctuation in text files (`normalize-whitespace-punctuation.cwl`)
- Command to run LIWC on tokenized text (`liwc_tokenized.cwl`)
- Command to save a directory to a subdirectory (`save-dir-to-subdir.cwl`)
- Command to save a list of files to a directory (`save-files-to-dir.cwl`)
- Command to merge csv files (`merge-csv.cwl`)
- Command to filter Named Entities extracted with frog (`frog-filter-nes.cwl`)
- Command to list al files in a directory (`ls.cwl`)
- Command to lowercase a text (`lowercase.cwl`)
- Command to clear xml elements (`clear-xml-elements.cwl`)

### Changed

- Command `saf_to_text.py` (`saf-to-text.cwl`) outputs space separated sentences
- Give Python commands a meaningful name (#5)
- Use `InitialWorkDirRequirement` instead of list of input files (#1 #16)

### Removed

- GUI (users are recommended to use Juyter notebooks for inspection and analysis tasks)
- Functionality to generate boilerplate Python commands and CWL files (moved to https://github.com/nlppln/nlppln-gen)
- Command to create a file list (use cwltool option `--make-template` instead) (#14)

## 0.1.0

### Added

#### Commands

- command to convert xml to text (`xml-to-txt.cwl`)
- command to extract word frequencies from saf files (`saf-to-freqs.cwl`)
- command to convert word frequencies to vocabulary ranked by frequency (`freqs.cwl`)
- command to list all files in directory (`ls.cwl`)
- command that takes as input a list of files and returns a list of lists of files (`chunk-list-of-files.cwl`)
- command that tokenizes text using ixa-pipe-tok (`ixa-pipe-tok.cwl`)
- command to download files (`download.cwl`)
- command to rename and copy files (rename-and-copy-files) -> must be used together with cwl-runner option `--relax-path-checks`
- command to convert Word documents (.doc and .docx) to text files using Apache Tika (`apachetika.cwl`)
- command to convert Word documents (.docx) to text files using `docx2txt` (`docx2txt.cwl`)
- command to convert frog output to saf (`frog-to-saf.cwl`)
- commands to frog a single text file and directory of files (`frog-single-text.cwl` and `frog-dir.cwl`)
- command line script to guess the language of all (txt) files in a directory (`language.cwl`)

#### Utils

- command to generate boilerplate python commands and CWL command line tools
- `WorkflowGenerator` class (from `scriptcwl`)
- generate `out_file` name from `in_file` name (using proper file extension)
- copy cwl files to default location (`$XDG_DATA_HOME/commonwl/`) (*CWL files are not yet copied automatically!*)

#### GUI

- functionality for inspecting named entities in text files
NLP Pipeline
============

|codacy_grade| |travis| |documentation| |pypi_version| |pypi_supported| |zenodo|

nlppln is a python package for creating NLP pipelines using `Common Workflow Language <http://www.commonwl.org/>`_ (CWL).
It provides steps for (generic) NLP functionality, such as tokenization,
lemmatization, and part of speech tagging, and helps users to construct workflows
from these steps.

A text processing step consist of a (Python) command line tool and a CWL
specification to use this tool.
Most tools provided by nppln wrap existing NLP functionality.
The command line tools are made with `Click <http://click.pocoo.org>`_, a Python
package for creating command line interfaces.

To create a workflow, you have to write a Python script:
::

  from nlppln import WorkflowGenerator

  with WorkflowGenerator() as wf:
    txt_dir = wf.add_input(txt_dir='Directory')

    frogout = wf.frog_dir(in_dir=txt_dir)
    saf = wf.frog_to_saf(in_files=frogout)
    ner_stats = wf.save_ner_data(in_files=saf)
    new_saf = wf.replace_ner(metadata=ner_stats, in_files=saf)
    txt = wf.saf_to_txt(in_files=new_saf)

    wf.add_outputs(ner_stats=ner_stats, txt=txt)

    wf.save('anonymize.cwl')

The resulting workflow can be run using a CWL runner, such as `cwltool <https://github.com/common-workflow-language/cwltool/>`_:

.. code-block:: sh

  cwltool anonymize.cwl --txt_dir /path/to/directory/with/txt/files/

For creating new (e.g., project specific) NLP functionality, you can use `nlppln-gen <https://github.com/nlppln/nlppln-gen>`_
to generate boilerplate (i.e., empty) command line tools and CWL specifications.

The full documentation can be found on `Read the Docs <http://nlppln.readthedocs.io/en/latest/>`_.

Installation
############

Install nlppln using pip:

.. code-block:: sh

  pip install nlppln

Please check the `installation guidelines <http://nlppln.readthedocs.io/en/latest/installation.html>`_ for additional required software.

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

.. |codacy_grade| image:: https://api.codacy.com/project/badge/Grade/24cd15fe1d9e4a51ab4be8c247e95c47
                     :target: https://www.codacy.com/app/jvdzwaan/nlppln?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=nlppln/nlppln&amp;utm_campaign=Badge_Grade
                     :alt: Codacy Badge

.. |travis| image:: https://travis-ci.org/nlppln/nlppln.svg?branch=master
              :target: https://travis-ci.org/nlppln/nlppln
              :alt: Build Status

.. |documentation| image:: https://readthedocs.org/projects/nlppln/badge/?version=latest
                     :target: http://nlppln.readthedocs.io/en/latest/?badge=latest
                     :alt: Documentation Status

.. |pypi_version| image:: https://badge.fury.io/py/nlppln.svg
                    :target: https://badge.fury.io/py/nlppln
                    :alt: PyPI version

.. |pypi_supported| image:: https://img.shields.io/pypi/pyversions/nlppln.svg
                      :target: https://pypi.python.org/pypi/nlppln
                      :alt: PyPI

.. |zenodo| image:: https://zenodo.org/badge/65198876.svg
              :target: https://zenodo.org/badge/latestdoi/65198876
              :alt: DOI
Running workflows
=================

To run a workflow created with nlppln, you need to use a CWL runner. The default
CWL runner cwltool is installed when you install nlppln. To run a tool:
::

	cwltool <workflow> <inputs>

The scriptcwl documentation contains some `tips and tricks for working with CWL files <http://scriptcwl.readthedocs.io/en/latest/cwl_tips_tricks.html>`_.

Windows
#######

To run a workflow created with nlppln on Windows, use the nlppln Docker container:
::

  cwltool --default-container nlppln:nlppln <workflow> <inputs>

Also, if you have to refer to file paths with spaces in them, use the
``--relax-path-checks`` option.

For more details about CWL on Windows, see the `Windows documentation <https://github.com/common-workflow-language/cwltool/blob/master/windowsdoc.md>`_.
Installation
============

Requirements
############

Before installing nlppln, please install:

* `Python 2 or 3 <https://www.python.org/downloads/>`_
* `node.js <https://nodejs.org/en/download/>`_
* `Docker <https://docs.docker.com/engine/installation/>`_

For more details about CWL on Windows, see the `Windows documentation <https://github.com/common-workflow-language/cwltool/blob/master/windowsdoc.md>`_.

Installing nlppln
#################

Install nlppln using pip:

.. code-block:: sh

  pip install nlppln

For development:

.. code-block:: sh

  git clone git@github.com:nlppln/nlppln.git
  cd nlppln
  pip install -r requirements.txt
  python setup.py develop

Run tests (including coverage) with:

.. code-block:: sh

  python setup.py test
.. nlppln documentation master file, created by
   sphinx-quickstart on Tue Nov 28 16:27:45 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the nlppln documentation!
====================================

nlppln is a python package for creating NLP pipelines using `Common Workflow Language <http://www.commonwl.org/>`_ (CWL).
It provides steps for (generic) NLP functionality, such as tokenization,
lemmatization, and part of speech tagging, and helps users to construct workflows
from these steps.

A text processing step consist of a (Python) command line tool and a CWL
specification to use this tool.
Most tools provided by nppln wrap existing NLP functionality.
The command line tools are made with `Click <http://click.pocoo.org>`_, a Python
package for creating command line interfaces.

To create a workflow, you have to write a Python script:
::

  from nlppln import WorkflowGenerator

  with WorkflowGenerator() as wf:
    txt_dir = wf.add_input(txt_dir='Directory')

    frogout = wf.frog_dir(in_dir=txt_dir)
    saf = wf.frog_to_saf(in_files=frogout)
    ner_stats = wf.save_ner_data(in_files=saf)
    new_saf = wf.replace_ner(metadata=ner_stats, in_files=saf)
    txt = wf.saf_to_txt(in_files=new_saf)

    wf.add_outputs(ner_stats=ner_stats, txt=txt)

    wf.save('anonymize.cwl')

The resulting workflow can be run using a CWL runner, such as `cwltool <https://github.com/common-workflow-language/cwltool/>`_:

.. code-block:: none

  cwltool anonymize.cwl --txt_dir /path/to/directory/with/txt/files/

For creating new (e.g., project specific) NLP functionality, you can use `nlppln-gen <https://github.com/nlppln/nlppln-gen>`_
to generate boilerplate (i.e., empty) command line tools and CWL specifications.

Flexible and reproducible text processing workflows
###################################################

One of the problems with existing NLP software is that to combine functionality
from different software packages, researchers have to write custom scripts.
Generally, these scripts duplicate at least some text processing tasks (e.g.,
tokenization), and need to be adapted when used for new datasets or in other software
or hardware environments. This has a negative impact on research reproducibility and
reuse of existing software.

The main advantages of nlppln are:

* Flexibility: it is easy to combine NLP tools written in different programming languages
* Reproducibility and portability: the resulting workflows can be run on a wide variety of hardware environments without changing them

Contents
########
.. toctree::
   :maxdepth: 2

   installation
   creating_workflows
   running_workflows
   tools
Creating workflows
==================

Pipelines or workflows can be created by writing a Python script:
::

	from nlppln import WorkflowGenerator

	with WorkflowGenerator() as wf:
		txt_dir = wf.add_input(txt_dir='Directory')

		frogout = wf.frog_dir(dir_in=txt_dir)
		saf = wf.frog_to_saf(in_files=frogout)
		ner_stats = wf.save_ner_data(in_files=saf)
		new_saf = wf.replace_ner(metadata=ner_stats, in_files=saf)
		txt = wf.saf_to_txt(in_files=new_saf)

		wf.add_outputs(ner_stats=ner_stats, txt=txt)

		wf.save('anonymize.cwl')

This workflow finds named entities in all Dutch text files in a directory. Named
entities are replaced with their type (PER, LOC, ORG). The output consists of
text files and a csv file that contains the named entities that have been replaced.

The workflow creation functionality in ``nlppln`` is provided by a library called
`scriptcwl <https://github.com/NLeSC/scriptcwl>`_. For a more detailed explanation
of how to create workflows, have a look at the
`scriptcwl documentation <http://scriptcwl.readthedocs.io/en/latest/>`_.

Setting workflow inputs
#######################

Wokflow inputs can be added by calling ``add_input()``:
::

	txt_dir = wf.add_input(txt_dir='Directory')

The ``add_input()`` method expects a ``name=type`` pair as input parameter.
The pair connects an input name (``txt_dir`` in the example) to a CWL type
(``'Directory'``). An overview of CWL types can be found in the
`specification <http://www.commonwl.org/v1.0/Workflow.html#CWLType>`_.

Check the `scriptcwl documentation <http://scriptcwl.readthedocs.io/en/latest/workflow_inputs.html>`_
to find out how to add optional workflow inputs and default values.

Adding processing steps
#######################

To add a processing step to the workflow, you have to call its method on the ``WorkflowGenerator`` object.
The method expects a list of (key, value) pairs as input parameters. (To find out what inputs a step
needs call ``wf.inputs(<step name>)``. This method prints all the inputs
and their types.) The method returns a list of strings containing output
names that can be used as input for later steps, or that can be connected
to workflow outputs.

For example, to add a step called ``frog-dir`` to the workflow, the
following method must be called:
::

    frogout = wf.frog_dir(dir_in=txt_dir)

In a next step, ``frogout`` can be used as input:
::
    saf = wf.frog_to_saf(in_files=frogout)
    txt = wf.saf_to_txt(in_files=saf)

Etcetera.

Listing steps
#############

To find out what steps are available in ``nlppln`` and to get copy/paste-able
specification of what needs to be typed to add a step to the workflow you can
type:
::

	print(wf.list_steps())

The result is:

.. code-block:: none

	Steps
  	apachetika............... out_files = wf.apachetika(in_files[, tika_server])
  	basic-text-statistics.... metadata_out = wf.basic_text_statistics(in_files, out_file)
  	chunk-list-of-files...... file_list = wf.chunk_list_of_files(chunk_size, in_files)
  	clear-xml-elements....... out_file = wf.clear_xml_elements(element, xml_file)
  	copy-and-rename.......... copy = wf.copy_and_rename(in_file[, rename])
  	docx2txt................. out_files = wf.docx2txt(in_files)
  	download................. out_files = wf.download(urls)
  	freqs.................... freqs = wf.freqs(in_files)
  	frog-dir................. frogout = wf.frog_dir(in_files[, skip])
  	frog-filter-nes.......... filtered_nerstats = wf.frog_filter_nes(nerstats[, name])
  	frog-single-text......... frogout = wf.frog_single_text(in_file)
  	frog-to-saf.............. saf = wf.frog_to_saf(in_files)
  	ixa-pipe-tok............. out_file = wf.ixa_pipe_tok(language, in_file)
  	language................. language_csv = wf.language(dir_in)
  	liwc-tokenized........... liwc = wf.liwc_tokenized(in_dir, liwc_dict[, encoding])
  	lowercase................ out_files = wf.lowercase(in_file)
  	ls....................... out_files = wf.ls(in_dir[, recursive])
  	merge-csv................ merged = wf.merge_csv(in_files[, name])
  	normalize-whitespace-punctuation metadata_out = wf.normalize_whitespace_punctuation(meta_in)
  	pattern-nl............... out_files = wf.pattern_nl(in_files)
  	rename-and-copy-files.... out_files = wf.rename_and_copy_files(in_files)
  	replace-ner.............. out_files = wf.replace_ner(metadata, in_files[, mode])
  	saf-to-freqs............. freqs = wf.saf_to_freqs(in_files[, mode])
  	saf-to-txt............... out_files = wf.saf_to_txt(in_files)
  	save-dir-to-subdir....... out = wf.save_dir_to_subdir(inner_dir, outer_dir)
  	save-files-to-dir........ out = wf.save_files_to_dir(dir_name, in_files)
  	save-ner-data............ ner_statistics = wf.save_ner_data(in_files)
  	textDNA-generate......... json = wf.textDNA_generate(dir_in, mode[, folder_sequences, name_prefix, output_dir])
  	xml-to-text.............. out_files = wf.xml_to_text(in_files[, tag])

	Workflows
  	anonymize................ ner_stats, out_files = wf.anonymize(in_files[, mode])

Setting workflow outputs
########################

When all steps of the workflow have been added, you can specify
workflow outputs by calling ``wf.add_outputs()``:
::

  wf.add_outputs(ner_stats=ner_stats, txt=txt)

In this case the workflow has two outputs, one called ``ner_stats``, which is a
csv file and one called ``txt``, which is a list of text files.

Saving workflows
################

To save a workflow call the ``WorkflowGenerator.save()`` method:
::

  wf.save('anonymize.cwl')

Other options when saving workflows are described in the `scriptcwl
documentation <http://scriptcwl.readthedocs.io/en/latest/saving_workflows.html>`_.
By default, ``nlppln`` saves workflows with embedded steps (``inline=True``).

Adding documentation
####################

To add documentation to your workflow, use the ``set_documentation()`` method:
::

	doc = """Workflow that replaces named entities in text files.

	Input:
		txt_dir: directory containing text files

	Output:
		ner_stats: csv-file containing statistics about named entities in the text files
		txt: text files with named enities replaced
	"""
	wf.set_documentation(doc)

Loading processing steps
########################

``nlppln`` comes with nlp functionality pre-loaded. If you need custom processing
steps, you can create them using `nlppln-gen <https://github.com/nlppln/nlppln-gen>`_.
To be able to add these custom processing steps to you workflow,
you have to load them into the ``WorkflowGenerator``. To load a single CWL file, do:
::

	wf.load(step_file='/path/to/step_or_workflow.cwl')

The ``step_file`` can also be a url.

To load all CWL files in a directory, do:
::

	wf.load(steps_dir='/path/to/dir/with/cwl/steps/')


Using a working directory
#########################

Once you need more functionality than nlppln provides, and start creating your
own processing steps, we recommend using a CWL working directory.
A CWL working directory is a directory containing all available CWL specifications.
To specify a working directory, do:
::

	from nlppln import WorkflowGenerator

  with WorkflowGenerator(working_dir='path/to/working_dir') as wf:
    wf.load(steps_dir='some/path/')
    wf.load(steps_dir='some/other/path/')

    # add inputs, steps and outputs

If you use a working directory when creating pipelines, nlppln copies all CWL files
to the working directory.

To copy these files manually, you can also use the ``nlppln_copy_cwl`` command on the
command line:
::

  nlppln_copy_cwl /path/to/cwl/working/dir

To copy CWL files from a different directory than the one containing the nlppln
CWL files, do:
::

  nlppln_copy_cwl --from_dir /path/to/your/dir/with/cwl/files /path/to/cwl/working/dir


If you use a working directory, please save your workflow using the ``wd=True`` option:
::

  wf.save('workflow.cwl', wd=True)

The workflow is saved in the working directory and then copied to you specified location.
Subsequently, the workflow should be run from the working directory.

Tips and tricks
###############

Create workflows you can run for a single file
----------------------------------------------

If you want to create a workflow that should be applied to each (text) file in a directory,
create a workflow that performs all the steps to a single file. Then, use this workflow
as a subworkflow that is scattered over a list of input files:
::

  from nlppln import WorkflowGenerator

  with WorkflowGenerator(working_dir='path/to/working_dir') as wf:
    wf.load(steps_dir='some/path/')

    in_dir = wf.add_input(in_dir='Directory')

    in_files = wf.ls(in_dir=in_dir)
    processed_files = wf.some_subworkflow(in_file=in_files, scatter='in_file', scatter_method='dotproduct' [, ...])

    wf.add_outputs(out_files=processed_files)


Having a workflow you can run for a single file makes it easier to test the workflow.

Test your workflow by running it for the largest or otherwise most complex file
-------------------------------------------------------------------------------

By running your workflow for the largest or otherwise most complex file, you can
identify problems, such as excessive memory usage, early and/or before running it
for all files in your dataset. There may, of course, still be problems with other
files, but starting analysis with the largest file is easy to do.


Use ``create_chunked_list`` and ``ls_chunk`` to run a workflow for a subset of files
------------------------------------------------------------------------------------

Sometimes running a workflow for all files in a directory takes too long, and you'd like
to run it for subsets of files. Using ``create_chunked_list``, you can create a JSON file
containing a division of the files in a directory in chunks. You can then create a workflow
that, instead of using ``ls`` to list all files in a directory, uses ``ls_chunk`` that
runs the workflow for a single chunk of files.

To create a division of the input files, do:
::

  python -m nlppln.commands.create_chunked_list [--size 500 --out_name output.json] /path/to/directory/with/input/files

The result is a JSON file named ``output.json`` that contains numbered chunks
containing 500 files each.

To run a workflow for a chunk of files, instead of all files in a directory, do:
::

  from nlppln import WorkflowGenerator

  with WorkflowGenerator(working_dir='path/to/working_dir') as wf:
    wf.load(steps_dir='some/path/')

    in_dir = wf.add_input(in_dir='Directory')
    chunks = wf.add_input(chunks='File')
    chunk_name = wf.add_input(name='string')

    in_files = wf.ls_chunk(in_dir=in_dir, chunks=chunks, name=chunk_name)
    processed_files = wf.some_subworkflow(in_file=in_files, scatter='in_file', scatter_method='dotproduct' [, ...])

    wf.add_outputs(out_files=processed_files)
