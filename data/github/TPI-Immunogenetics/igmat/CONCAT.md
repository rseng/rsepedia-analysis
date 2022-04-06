# IgMAT: Antibody Multispecies Annotation Tool

IgMAT is a tool for the automatic discrimination and annotation of antibody sequences, specifically designed to be integrated into analysis pipelines or being used as a cli tool. IgMAT is highly customizable, allowing the addition of custom antibody sequences datasets and generating a range of output formats including a bed file of FR and CDR coordinates allowing simple downstream analysis on individual regions.

## Requirements ##
* Python version 3.8 or greater
* [Hmmer](http://hmmer.org/download.html)

## Installing ##
IgMAT can be installed locally or in a python environment: 

    git clone git@github.com:TPI-Immunogenetics/igmat.git
    cd igmat
    python3 -m venv env
    source env/bin/activate
    pip install ./

The default HMM model generated from IMGT data needs to be build in order to run IgMAT:

    igmat build

Once done, to exit the environment, type:

    deactivate
## Configuration ##
IgMAT will automatically generate a configuration file in the user home directory: 

    ~/.igmat/config.yaml

All generated hmmer models will be stored in this folder and the configuration file will be automatically initialised with the path containing the hmmer executables. For each one of the cli tools, the option `--hmmerpath` is available to temporarily override the hmmer path.

## CLI tools ##
IgMAT comes with a set of cli tools for handling custom HMM models and processing data:

    igmat <tool> --help

Where `tool` is one of the following:
 
 - **run**: process the input sequence/file and run the annotation script.
 - **build**: build custom HMM models 
 - **list**: handles the available HMM models

## Examples ##
IgMAT can be used as a stand alone tool, or embedded into a custom script. Please check the [examples](/docs/examples.md) and a tutorial about [embedding](/docs/embedding.md) IgMAT in your scripts.# Examples #

IgMAT can be used to annotate sequences with the default IMGT dataset, or specific datasets can be created:

## Annotating a sequence ##
To annotate an input sequence with the default IMGT database:

	igmat run -i DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH -v

It will print a representation of the annotated sequence:

```
Query 'Input sequence':

 HMMR match: Vicugna+pacos_H [eValue: 2.000e-55]

   0 -----------FR1------------|----CDR1----|-------FR2-------|---CDR2---|-  70
   0 DVQLVESGG-GSVQAGGSLRLSCAVS|GSTY----SPCT|TGWYRQAPGKEREWVSS|ISSP---GTI|Y  70

  70 -----------------FR3------------------|-----CDR3----|----FR4----| 135
  70 YQDSVK-GRFTISRDNAKNTVYLQMNSLQREDTGMYYC|QI-QCGVRSIREY|WGQGTQVTVSS| 135
```

## Creating a custom dataset ##
To create a custom dataset, input files with sequences from V and J regions are needed, where the name of the files must adhere the following format: 
	
		<name>_<chain><type>
Where `name` is a species or sample name, and `type` is either *V* or *J*.
An example of input sequences for generating custom datasets is distributed in the `test/build` folder:

	igmat build -i ./test/build -n cattle
The command will generate a `cattle` dataset that will be available to be used with IgMAT with the command:
	
	igmat run -i DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH -m cattle# Embedding IgMAT #

This example shows how to embed IgMAT in your code. 

```python

from igmat import igmat

# The sequence that needs to be annotated
sequence = "DVQLVESGGGSVQAGGSLRLSCAVSGSTYSPCTTGWYRQAPGKEREWVSSISSPGTIYYQDSVKGRFTISRDNAKNTVYLQMNSLQREDTGMYYCQIQCGVRSIREYWGQGTQVTVSSHHHHHH"

# Generate the list of results
resultList = igmat.annotate(sequence, 'IMGT')
if not resultList:
  raise Exception('No result found')

# Iterate the list of results
for result in resultList:

  # Print the resulting sequence
  print(result)

  # Print some details
  print(result.sequence)
  print(result.start)
  print(result.end)

  # Print the annotations
  for feature in result.annotations():
    print('Annotation {type}: {start}-{stop}'.format(
      type=feature['type'],
      start=feature['start'],
      stop=feature['stop']
      ))
```