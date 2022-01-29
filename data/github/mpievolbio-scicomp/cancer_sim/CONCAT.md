---
title: 'CancerSim: A Cancer Simulation Package for Python 3'
tags:
  - stochastic simulation
  - tumour growth
  - tumour sampling
  - cancer biology
authors:
  - name: Luka Opasic
    orcid: 0000-0001-7595-1722
    affiliation: 1
  - name: Jacob G. Scott
    orcid: 0000-0003-2971-7673
    affiliation: 2
  - name: Arne Traulsen
    orcid: 0000-0002-0669-5267
    affiliation: 1
  - name: Carsten Fortmann-Grote
    orcid: 0000-0002-2579-5546
    affiliation: 1
affiliations:
  - name: Max-Planck-Institute for Evolutionary Biology, Plön, Germany
    index: 1
  - name: Cleveland Clinic, Cleveland, OH, US
    index: 2
date: 17 March 2020
bibliography: references.bib
---


Summary
----------

Cancer is a group of complex diseases characterized by excessive cell
proliferation, invasion, and destruction of the surrounding tissue
[@Kumar2017]. Its high division and mutation rates
lead to excessive genetic diversity among tumour cells (intra tumour genetic
heterogeneity). As a consequence, tumours can adapt very efficiently to
environmental pressures, in particular to cancer therapy [@Turajlic2019].
This process is known as somatic evolution of cancer.  

Throughout most of its existence a tumour is inaccessible to direct 
observation and experimental evaluation through genetic sequencing of tumour
samples. Therefore, computational modelling
 can be useful to study many aspects of cancer. Examples where theoretical 
 models can be of great use include (i) early carcinogenesis, as lesions are clinically
observable when they already contain millions of cells, (ii) seeding of metastases, 
and (iii) cancer cell dormancy [@Altrock2015].

### Statement of Need
Advanced cancer simulation software [@waclaw2015] often exhibit a prohibitively steep
learning curve especially for new students in the field of somatic evolution of cancer. A software package that is 
 accessible, simple to use, and yet covers the essential biological processes of cancer growth is needed to provide an entry
  point for students and newcomers to mathematical oncology. 

### Cancer growth model
Here we present `CancerSim`, a software that simulates somatic evolution of
tumours. The software produces virtual spatial tumours with variable extent of
intra tumour genetic heterogeneity and realistic mutational profiles (i.e., the order
of appearance of mutations and their distribution among tumour cells).
Simulated tumours can be virtually sampled at random or in specified positions.
By binning the frequency of mutations in
each sample and convolving the resulting histogram with a beta-binomial
distribution to add realistic sequencing noise [@Williams2016], the simulation produces mutation profiles from different tumour regions. 
This makes the software useful for studying various sampling strategies in clinical cancer
diagnostics such as needle biopsy sampling or liquid biopsy sampling. An early version of this 
cancer evolution model was used to simulate tumours subjected to sampling for 
classification of mutations based on their abundance [@Opasic2019].

Our model is abstract, not specific to any neoplasm type, and does not
consider a variety of biological features commonly found in the neoplasm,
such as blood vessels, immune cells, nutrient availability, and
architecture of the tumour surroundings. It resembles the most
superficially spreading tumours like carcinoma *in-situ*, skin cancers, or
gastric cancers, but it can be used to model any tumour on this abstract
level.

The tumour is simulated using a two-dimensional, on-lattice, agent-based
model. The tumour lattice structure is established by a sparse matrix
whose non-zero elements correspond to the individual cells. Each cell is
surrounded by eight neighbouring cells (Moore neighbourhood). The value
of the matrix element is an index pointing to the last mutation the cell
acquired in the list of mutations which is updated in each simulation
step.

### Simulation parameters
The simulation advances in discrete time-steps. In each simulation step,
every tumour cell on the lattice that has an unoccupied neighbour can
divide with a certain probability (controlled through the parameter `division_probability`). The
daughter cell resulting from a cell division inherits all mutations from
the parent cell and acquires a new mutation with a given probability
(`mutation_probability`). 
A new mutation that changes death and/or birth probabilities can be introduced
into a random cell at a specific time step defined by `adv_mutation_wait_time`.
By changing the fitness parameters of a mutant cell `adv_mutant_division_probability`
and `adv_mutant_death_probability` one can model various evolutionary processes
like emergence of a faster dividing sub-clone or selective effects of drug treatment.

The simulation allows the acquisition of more than one mutational event per cell
(`number_of_mutations_per_division`). In that case, variable amounts of
sequencing noise [@Williams2016] can be added to make
the output data biologically more realistic. The key parameters
`number_of_generations`, 
`division_probability` and `death_probability`
determine the final size of the tumour, while the degree of intra tumour heterogeneity can 
be varied by changing the `mutation_probability` parameter. 
For neutral tumour evolution, parameter `adv_mutant_division_probability`
and `adv_mutant_death_probability` must be the same as `division_probability`
and `death_probability`. 

Throughout the cancer growth phase, `CancerSim` stores information about
the parent cell and a designation of newly acquired mutations for every
cell. Complete mutational profiles of cells are reconstructed a
posteriori based on the stored lineage information.

The division rules which allow only cells with empty neighbouring nodes
to divide, cause exclusively peripheral growth and complete absence of
dynamics in the tumour centre. To allow for variable degrees of growth
inside the tumour, we introduced a death process. At every time step,
after all cells attempt their division, a number of random cells die according
to `death_probability` and `adv_mutant_death_probability` and
yield their position to host a new cancer cell in a subsequent time
step.

### Simulation results
After the simulation, the tumour matrix, and the lists of lineages and
frequencies of each mutation in the tumour are exported to files.
Furthermore, the virtual tumour can be sampled and a histogram over the
frequency of mutations will be visualised. Alternatively, a saved tumour
can be loaded from file and then subjected to the sampling process.

### Outlook
Future investigations and expansions of `CancerSim` will focus on  sampling of tumour specimens
in a specific spatial pattern. Additionally, the effects of chemotherapy can be
modelled by the introduction of different modes of cell death. Our model only
considers sequencing noise for mutated cells. Future improvements of our model
will also consider other sources of sequencing noise, e.g., read errors from
non-mutated sites which would lead to false-positive cancer detections.


Download and Installation
-------------------------
`CancerSim` is written in Python (version \>3.5). We recommend to install
it directly from the source code hosted at github <https://github.com/mpievolbio-scicomp/cancer_sim>.

 Detailed instructions including creation of a
`conda` environment are given in the online documentation at
<https://cancer-sim.readthedocs.io/en/master/include/README.html#installation>.
After installation, the software is available as a python module `casim`.

High-level functionality
-------------------------
### Setting up the cancer simulation parameters
The parameters of the cancer simulation are specified in a python module or
programmatically via the `CancerSimulationParameters` class. 
The table below lists all parameters, their function and acceptable values.
[.,.] denotes a range of permissive parameter values.

Parameter name | function | valid options
---------------|----------|--------------
`matrix_size` | Number of mesh points in each dimension  | >0
`number_of_generations` | Number of generations to simulate  | >0
`division_probability` | Probability of cell division per generation  | [0,1]
`adv_mutant_division_probability` | Probability of division for cells with advantageous mutation  | [0,1]
`death_probability` | Fraction of cells that die per generation  | [0,1]
`adv_mutant_death_probability` | Fraction of cells with advantageous mutation that die per generation  | [0,1]
`mutation_probability` | Probability of mutations  | [0,1]
`adv_mutant_mutation_probability` | Mutation probability for the adv. cells  | [0,1]
`number_of_mutations_per_division` | Number of mutations per cell division  | >=0
`adv_mutation_wait_time` | Number of generations after which adv. mutation occurs  | >0
`number_of_initial_mutations` | Number of mutations present in first cancer cell | >=0
`tumour_multiplicity` | Tumour multiplicity  | "single", "double"
`read_depth` | Sequencing read depth  | read length * number of reads / genome length
`sampling_fraction` | Fraction of cells to be sampled  | [0,1]
`sampling_positions` | A list of [(x,y)] coordinates specifying the grid points where to take samples from the tumour. If unset or None, a single random position is chosen. | [0, matrix_size-1]
`plot_tumour_growth` | Plot the tumour growth curve  | True, False
`export_tumour` | Export the tumour growth data to file  | True, False


The file [`params.py`](https://github.com/mpievolbio-scicomp/cancer_sim/blob/master/params.py) can serve as a 
template to setup a simulation with all above parameters. Users should start with the template
and adjust the parameters as needed for
their application by setting experimentally or theoretically known values or by
calibrating the simulation output against experiments or other models.

### Run the simulation
The simulation is started
either from the command line, through a python script or from within an
interactive python session. The command line interface accepts as arguments a
random seed, a path for the output directory and a log level. More detailed instructions
are given in the [online documentation](https://cancer-sim.readthedocs.io/en/latest/include/README.html#high-level-functionality) and in the [reference manual](https://cancer-sim.readthedocs.io/en/latest/refman.html#casim.casim.CancerSimulator.run).

For the impatient, we also provide a jupyter notebook with a more condensed
version of the above example (grid size 20x20) at
`docs/source/include/notebooks/quickstart_example.ipynb`. An interactive version
can be launched on the [Binder service](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim.git/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb).

A further notebook demonstrates the possibility to dump a simulation to disk,
reload it and continue the simulation with optionally changed parameters
(`docs/source/include/notebooks/run_dump_reload_continue.ipynb`). This
feature could be exploited to simulate cancer dormancy with `CancerSim`.

### Output
After the run has finished, all output is found in the specified output
directory.
Simulation output consists of the mutation list and a sample file representing
the sampled portion of the tumour.
Two pdf files are generated showing the mutation histograms of
the whole tumour and the sampled portion, respectively. If requested by the
parameter `plot_tumour_growth`, a third pdf is generated showing the
tumour growth curve (number of tumour cells over time).
Furthermore, the output contains serialized versions ("pickles") of the tumour
geometry as a 2D matrix, the death list, and the
mutation list. Another directory contains all logged information from the run.
The amount of logged information depends on the chosen log level.

As an example and possible starting point for further investigations beyond the
produced plots, the quick start example notebook demonstrates how to load the
tumour matrix into memory and how to visualize the mutational profile as a
heat map as shown in Figure 1.

![Tumour matrix as a 2D grid with colored cells indicating the order in which mutations appear.](img/mutation_profile.png)

The figure shows the central part of the tumour matrix that results from the
template parameter file `params.py` on the X-Y grid with
colors indicating the order in which mutations appear. 0 (deep purple) indicates
no mutation in this cell.

Documentation and support
-------------------------

The API reference manual and community guidelines including directions
for contributors and bug reports are given in the online documentation at
<https://cancer-sim.readthedocs.io>. 

Testing and continuous integration
----------------------------------
Each commit to the code base triggers a build and test run on
https://travis-ci.org. The build and test status can be observed at [https://travis-ci.org/mpievolbio-scicomp/cancer_sim](https://travis-ci.org/mpievolbio-scicomp/cancer_sim).


References
----------

[![Build Status](https://travis-ci.org/mpievolbio-scicomp/cancer_sim.svg?branch=master)](https://travis-ci.org/mpievolbio-scicomp/cancer_sim)
[![Documentation Status](https://readthedocs.org/projects/cancer-sim/badge/?version=latest)](https://cancer-sim.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02436/status.svg)](https://doi.org/10.21105/joss.02436)

## Documentation
Documentation for `CancerSim`, including this README and the API reference manual
 is hosted on [readthedocs](https://cancer-sim.readthedocs.io).

Background
----------
Cancer is a group of complex diseases characterized by excessive cell
proliferation, invasion, and destruction of the surrounding tissue
 \[[1](#ref-kumar2017)\]
. Its high division and mutation rates lead to excessive
intratumour genetic heterogeneity which makes cancer highly adaptable to
environmental pressures such as therapy 
\[[2](#ref-turajlic2019)\]
. 
This process is known as somatic evolution of cancer. Throughout most
of its existence a tumour is inaccessible to direct observation and
experimental evaluation. Therefore, computational modelling can be
useful to study many aspects of cancer. Some examples where theoretical
models can be of great use include early carcinogenesis, as lesions are
clinically observable when they already contain millions of cells,
seeding of metastases, and cancer cell dormancy 
\[[3](#ref-altrock2015)\]
.

### Statement of Need
Advanced cancer simulation software [@waclaw2015] often exhibit a prohibitively steep
learning curve especially for new students in the field of somatic evolution of cancer. A software package that is 
 accessible, simple to use, and yet covers the essential biological processes of cancer growth is needed to provide an entry
  point for students and newcomers to mathematical oncology. 
  
Here, we present `CancerSim`, a software that simulates somatic evolution
of tumours. The software produces virtual spatial tumours with variable
extent of intratumour genetic heterogeneity and realistic mutational
profiles. Simulated tumours can be subjected to multi-region sampling to
obtain mutation profiles that are realistic representation of the
sequencing data. This makes the software useful for studying various
sampling strategies in clinical cancer diagnostics. An early version of
this cancer evolution model was used to simulate tumours subjected to
sampling for classification of mutations based on their abundance
\[[4](#ref-opasic2019)\].

### Cancer growth model
Our model is abstract, not specific to any neoplasm type, and does not
consider a variety of biological features commonly found in neoplasm
such as vasculature, immune contexture, availability of nutrients, and
architecture of the tumour surroundings. It most closely resembles the 
superficially spreading tumours like carcinoma in situ, skin cancers, or
gastric cancers, but it can be used to model any tumour on this abstract
level. 

The tumour is simulated using a two-dimensional, on-lattice, agent-based
model. The tumour lattice structure is established by a sparse matrix
whose non-zero elements correspond to the individual cells. Each cell is
surrounded by eight neighbouring cells (Moore neighbourhood). The value
of the matrix element is an index pointing to the last mutation the cell
acquired in the list of mutations which is updated in each simulation
step.

### Simulation parameters
The simulation advances in discrete time-steps. In each simulation step,
every tumour cell in the tumour that has an unoccupied neighbour can
divide with a certain probability (set by the parameter `division_probability`). The
daughter cell resulting from a cell division inherits all mutations from
the parent cell and acquires a new mutation with a given probability
`mutation_probability`). A new mutation that changes death and birth probability of cell can be introduced
at into random cell at the specific time step defined by `adv_mutation_wait_time`.
By changing fitness parameters of a mutant cell `adv_mutant_division_probability`
and `adv_mutant_death_probability` one can model various evolutionary processes
like emergence of a faster dividing sub-clone or selective effects of a drug treatment.

The simulation allows the acquisition of more than one mutational event per cell
(`number_of_mutations_per_division`). In that case, variable amounts of
sequencing noise \[[6](#ref-williams2016)\] can be added to make
the output data more biologically realistic. The key parameters
`number_of_generations`, `division_probability` and `death_probability`
determine the final size of the tumour, while the degree of intratumour heterogeneity can 
be varied by changing the `mutation_probability` parameter. 
For neutral tumour evolution, parameter `adv_mutant_division_probability`
and `adv_mutant_death_probability` must be the same as `division_probability`
and `death_probability`.

Throughout the cancer growth phase, `CancerSim` stores information about
the parent cell and a designation of newly acquired mutations for every
cell. Complete mutational profiles of cells are reconstructed a
posteriori based on the stored lineage information.

The division rules which allow only cells with empty neighbouring nodes
to divide, cause exclusively peripheral growth and complete absence of
dynamics in the tumour centre. To allow for variable degree of growth
inside the tumour, we introduced a death process. At every time step,
after all cells attempt their division, a number of random cells die according
to `death_probability` and `adv_mutant_death_probability` and
yield their position to host a new cancer cell in a subsequent time
step.

### Simulation results
After the simulation, the tumour matrix, and the lists of lineages and
frequencies of each mutation in the tumour are exported to files.
Furthermore, the virtual tumour can be sampled and a histogram over the
frequency of mutations will be visualised. Alternatively, a saved tumour
can be loaded from file and then be subjected to the sampling process.


Installation 
------------
`CancerSim` is written in Python (version \>3.5). We recommend to install
it directly from the source code. To download the code:

**EITHER** clone the repository:

    $> git clone https://github.com/mpievolbio-scicomp/cancer_sim.git

**OR** download the source code archive:

    $> wget https://github.com/mpievolbio-scicomp/cancer_sim/archive/master.zip
    $> unzip master.zip -d cancer_sim

Change into the source code directory

    $> cd cancer_sim

We provide for two alternatives to install the software after it was
downloaded:

### Alternative 1: `conda`
#### New conda environment
We [provide](provide) an `environment.yml` to be consumed by `conda`. To create a
fully self-contained conda environment (named `casim`):

    $> conda env create -n casim --file environment.yml

This will also install the cancer simulation code into the new
environment.

To activate the new conda environment:

    $> source activate casim

or

    $> conda activate casim

if you have set up conda appropriately.

#### Install into existing and activated conda environment
To install the software into an already existing environment:

    $> conda activate <name_of_existing_conda_environment>
    $> conda env update --file environment.yml

### Alternative 2: `pip`
The file `requirements.txt` is meant to be consumed by `pip`:

    $> pip install -r requirements.txt [--user]

The option `--user` is needed to install without admin privileges.

### Installed module
After installation, `CancerSim` is available in python as the `casim` module.
E.g. in a python script, one would import the module as:

    >>> from casim import casim

Testing
-------
Although not strictly required, we recommend to run the test suite after
installation. Simply execute the `run_tests.sh` shell script:

    $> ./run_tests.sh

This will generate a test log named `casim_test@<timestamp>.log` with
`<timestamp>` being the date and time when the test was run. You should
see an `OK` at the bottom of the log. If instead errors or failures are
reported, something is wrong with the installation or the code itself.
Feel free to open a github issue at
<https://github.com/mpievolbio-scicomp/cancer_sim/issues> and attach the
test log plus any information that may be useful to reproduce the error
(version hash, computer hardware, operating system, python version, a
dump of `conda env export` if applicable).

The test suite is automatically run after each commit to the code base.
Results are published on
[travis-ci.org](https://travis-ci.org/mpievolbio-scicomp/cancer_sim).

High-level functionality
-------------------------
### Setting the simulation parameters
The parameters of the cancer simulation are specified in a python module or
programmatically via the `CancerSimulationParameters` class. A
documented example `params.py` is included in the source code (under
`casim/params.py`) and reproduced here:

    $> cat casim/params.py
    ################################################################################
    #                                                                              #
    # Commented casim parameter input file.                                        #
    # Valid settings are indicated in parentheses at the end of each comment line. #
    # [0,1] stands for the closed interval from 0 to 1, including the limits; ||   #
    # means "or".                                                                  #
    #                                                                              #
    ################################################################################
    
    # Number of mesh points in each dimension (>0)
    matrix_size = 1000

    # Number of generations to simulate (>0).
    number_of_generations = 20

    # Probability of cell division per generation ([0,1]).
    division_probability = 1

    # Probability of division for cells with advantageous mutation ([0,1]).
    adv_mutant_division_probability = 1

    # Fraction of cells that die per generation ([0,1]).
    death_probability = 0.1

    # Fraction of cells with advantageous mutation that die per generation ([0,1]).
    adv_mutant_death_probability = 0.0

    # Probability of mutations ([0,1]).
    mutation_probability = 1

    # Mutation probability for the adv. cells ([0,1]).
    adv_mutant_mutation_probability = 1

    # Number of mutations per cell division (>=0).
    number_of_mutations_per_division = 10

    # Number of generations after which the beneficial mutation occurs (>=1).
    adv_mutation_wait_time = 10

    # Number of mutations present in first cancer cell (>=0).
    number_of_initial_mutations = 150

    # Tumour multiplicity ("single" || "double").
    tumour_multiplicity = "single"

    # Sequencing read depth (read length * number of reads / genome length).
    read_depth = 100

    # Fraction of cells to be sampled ([0,1]).
    sampling_fraction = 0.1
        
	# Sampling position (list of (x,y) coordinates in the range [0,matrix_size-1]).
	# If left blank or None, random position will be chosen.
	# sampling_positions = None # This will randomly set a single sampling position.
    sampling_positions = [(500,500),(490,490)]
    
    # Plot the tumour growth curve (True || False).
    plot_tumour_growth = True
        
    # Export the tumour growth data to file (True || False).
    export_tumour = True

Here, we simulate a single 2D tumour on a 1000x1000 grid (`matrix_size=1000`) for a total of
20 generations (`number_of_generations=20`).
On average, both healthy and mutant cells divide once per generation
(`division_probability`). The first cancer cell carries 150 mutations
(`number_of_initial_mutations=150`); both healthy and mutant cells aquire 10 new
mutations (`number_of_mutations_per_division=10`) in
each generation with a certainty of 100% (`mutation_probability=1`). The advantageous mutation happens in
the 10th generation (`adv_mutation_wait_time=10`). Mutant cells with advantageous mutations live on forever
(`adv_mutant_death_probability=0`) while healthy cells die with a rate of 0.1
per generation (`death_probability=0.1`).

After completion of the last generation, two spatial samples are taken, one from the tumour center and one from a
slightly more lateral
position (`sampling_positions = [(500,500),(490,490)]`). Each sample contains 10% (with respect to the whole tumour size)
 closely positioned tumour cells
(`sampling_fraction=0.1`). The samples are subject to genetic sequencing with a read depth
of 100 (`read_depth=100`). The data is written to disk (`export_tumour=True`)
and plots showing the mutation histograms for the whole tumour as well as for the sampled part of the
tumour are generated. Furthermore, a plot showing the tumour growth over time is
saved (`plot_tumour_growth=True`).
 
Users should start with the template
and adjust the parameters as needed for
their application by setting experimentally or theoretically known values or by
calibrating the simulation output against experiments or other models.

### Run the example
The simulation is started from the command line. The syntax is

    $> python -m casim.casim [-h] [-s SEED] [-p PARAMS] [-o DIR]

`SEED` is the random seed. Using the
same seed in two simulation runs with identical parameters results in
identical results. If not given, `SEED` defaults to 1. `PARAMS` should point to
a python parameter file. If not given, it defaults to `params.py` in the current
working directory. If that file does not exist, default parameters are assumed.
`DIR` specifies the directory where to store the
simulation log and output data. If not given, output will be stored in
the directory `casim_out` in the current directory. For each seed, a
subdirectory `cancer_SEED` will be created. If that subdirectory already
exists because an earlier run used the same seed, the run will abort.
This is a safety catch to avoid overwriting data from previous runs.

### Output
After the run has finished, you should find the results in the specified output
directory:

    $> ls out/cancer_1/simOutput  
    growthCurve.pdf  mut_container.p              sample_out_490_490.txt
    mtx.p            sampleHistogram_490_490.pdf  sample_out_500_500.txt
    mtx_VAF.txt      sampleHistogram_500_500.pdf  wholeTumourVAFHistogram.pdf   
    
Let's take a look at the `.txt`  files. They contain the simulation output:
`mtx_VAF.txt` is a datafile with three columns: `mutation_id` lists the index of
each primary mutation, `additional_mut_id` indexes the subsequent mutations that occur in a cell of
a given `mutation_id`; `frequency` is the frequency which at a given mutation occurs.

Corresponding to the given sample positions, there is one
`sample_out_XXX_YYY.txt` and one `sampleHistogram_XXX_YYY.pdf` for each
position. The `.txt` filel lists all mutations of the artificial sample
taken from the whole tumour. Columns are identical to `mtx_VAF.txt`.
 
The `.pdf` files are plots of sampled
tumour histogram. `wholeTumourVAFHistogram.pdf` is the histogram for the
complete tumour. You should see figures similar to these:

Whole tumour histogram:  
![Whole tumour histogram](img/example_whole_tumour.png "Whole tumour histogram")  

Central sample histogram:   
![Sampled tumour histogram](img/example_sampled_tumour_pos1.png "Center sample histogram")  

Lateral sample histogram:   
![Sampled tumour histogram](img/example_sampled_tumour_pos2.png "Lateral sample histogram")  

The remaining output files are serialized versions ("pickles") of the tumour
geometry as a 2D matrix (`mtx.p`) and the
mutation list (list of tuples listing the cancer cell index and the mutation ID of each
tumour cell, `mut_container.p`).

### Example notebooks
An example demonstrating how to parametrize the simulation through the
`CancerSimulationParameters` API is provided in the accompanying jupyter notebook at
[quickstart_example.ipynb`](https://github.com/mpievolbio-scicomp/cancer_sim/blob/master/docs/source/include/notebooks/quickstart.ipynb). Launch it in
[binder](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim.git/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb).

In
[`run_dump_reload_continue.ipynb`](https://github.com/mpievolbio-scicomp/cancer_sim/blob/master/docs/source/include/notebooks/run_dump_reload_continue.ipynb), we
demonstrate how to use the restart capability to modifiy tumour growth parameters in the
middle of a run. In this way, one can model different phases of tumour growth,
e.g. tumour dormancy or onset of cancer therapy.

Community Guidelines
--------------------
As an Open Source project, we welcome all contributions to `CancerSim`. We
recommend the usual github workflow: Fork this repository, commit
your changes and additions to the fork and create a pull request back to the
master branch on this repository. If uncertain about anything, please create an
issue at [https://github.com/mpievolbio-scicomp/cancer_sim/issues](https://github.com/mpievolbio-scicomp/cancer_sim/issues).

Comments, bug reports, or other issues as well as requests
for support should be submitted as a [github
issue](https://github.com/mpievolbio-scicomp/cancer_sim/issues). Please check
the [list of issues](https://github.com/mpievolbio-scicomp/cancer_sim/issues?q=is%3Aissue) if your problem has already been addressed. We will do our best to respond in a timely manner.

References
----------
<a name=ref-kumar2017></a> \[1\] J. C. A. Vinay Kumar Abul K. Abbas,
*Robbins Basic Pathology*, 10th ed. (Elsevier, 2017). ISBN: 9780323353175. 

<a name=ref-turajlic2019></a> \[2\] S. Turajlic, A. Sottoriva, T. Graham,
and C. Swanton, Nat Rev Genet (2019). DOI:
[10.1038/s41576-019-0114-6](https://dx.doi.org/10.1038/s41576-019-0114-6)


<a name=ref-altrock2015></a> \[3\] P. M. Altrock, L. L. Liu, and
F. Michor, Nat Rev Cancer **15**, 730 (2015).
 DOI:
[10.1038/nrc4029](https://dx.doi.org/10.1038/nrc4029)

<a name=ref-opasic2019></a> \[4\] L. Opasic, D. Zhou, B. Werner, D.
Dingli, and A. Traulsen, BMC Cancer **19**, 403 (2019).
 DOI:
[10.1186/s12885-019-5597-1](https://dx.doi.org/10.1186/s12885-019-5597-1)

<a name=ref-waclaw2015></a> \[5\] B. Waclaw, I. Bozic, M. E. Pittman, R. H.  Hruban, B. Vogelstein, and M. A. Nowak, Nature **525**, 261 (2015). DOI: [10.1038/nature14971](https://dx.doi.org/10.1038/nature14971)

<a name=ref-williams2016></a> \[6\] M. J. Williams, B. Werner, C. P. Barnes,
T. A. Graham, and A. Sottoriva, Nature Genetics **48**, 238 (2016).  DOI:
[10.1038/ng.3489](https://dx.doi.org/10.1038/ng.3489)

[![Build Status](https://travis-ci.org/mpievolbio-scicomp/cancer_sim.svg?branch=master)](https://travis-ci.org/mpievolbio-scicomp/cancer_sim)
[![Documentation Status](https://readthedocs.org/projects/cancer-sim/badge/?version=latest)](https://cancer-sim.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02436/status.svg)](https://doi.org/10.21105/joss.02436)

## Documentation
Documentation for `CancerSim`, including this README and the API reference manual
 is hosted on [readthedocs](https://cancer-sim.readthedocs.io).

Background
----------
Cancer is a group of complex diseases characterized by excessive cell
proliferation, invasion, and destruction of the surrounding tissue
 \[[1](#ref-kumar2017)\]
. Its high division and mutation rates lead to excessive
intratumour genetic heterogeneity which makes cancer highly adaptable to
environmental pressures such as therapy 
\[[2](#ref-turajlic2019)\]
. 
This process is known as somatic evolution of cancer. Throughout most
of its existence a tumour is inaccessible to direct observation and
experimental evaluation. Therefore, computational modelling can be
useful to study many aspects of cancer. Some examples where theoretical
models can be of great use include early carcinogenesis, as lesions are
clinically observable when they already contain millions of cells,
seeding of metastases, and cancer cell dormancy 
\[[3](#ref-altrock2015)\]
.

### Statement of Need
Advanced cancer simulation software [@waclaw2015] often exhibit a prohibitively steep
learning curve especially for new students in the field of somatic evolution of cancer. A software package that is 
 accessible, simple to use, and yet covers the essential biological processes of cancer growth is needed to provide an entry
  point for students and newcomers to mathematical oncology. 
  
Here, we present `CancerSim`, a software that simulates somatic evolution
of tumours. The software produces virtual spatial tumours with variable
extent of intratumour genetic heterogeneity and realistic mutational
profiles. Simulated tumours can be subjected to multi-region sampling to
obtain mutation profiles that are realistic representation of the
sequencing data. This makes the software useful for studying various
sampling strategies in clinical cancer diagnostics. An early version of
this cancer evolution model was used to simulate tumours subjected to
sampling for classification of mutations based on their abundance
\[[4](#ref-opasic2019)\].

### Cancer growth model
Our model is abstract, not specific to any neoplasm type, and does not
consider a variety of biological features commonly found in neoplasm
such as vasculature, immune contexture, availability of nutrients, and
architecture of the tumour surroundings. It most closely resembles the 
superficially spreading tumours like carcinoma in situ, skin cancers, or
gastric cancers, but it can be used to model any tumour on this abstract
level. 

The tumour is simulated using a two-dimensional, on-lattice, agent-based
model. The tumour lattice structure is established by a sparse matrix
whose non-zero elements correspond to the individual cells. Each cell is
surrounded by eight neighbouring cells (Moore neighbourhood). The value
of the matrix element is an index pointing to the last mutation the cell
acquired in the list of mutations which is updated in each simulation
step.

### Simulation parameters
The simulation advances in discrete time-steps. In each simulation step,
every tumour cell in the tumour that has an unoccupied neighbour can
divide with a certain probability (set by the parameter `division_probability`). The
daughter cell resulting from a cell division inherits all mutations from
the parent cell and acquires a new mutation with a given probability
`mutation_probability`). A new mutation that changes death and birth probability of cell can be introduced
at into random cell at the specific time step defined by `adv_mutation_wait_time`.
By changing fitness parameters of a mutant cell `adv_mutant_division_probability`
and `adv_mutant_death_probability` one can model various evolutionary processes
like emergence of a faster dividing sub-clone or selective effects of a drug treatment.

The simulation allows the acquisition of more than one mutational event per cell
(`number_of_mutations_per_division`). In that case, variable amounts of
sequencing noise \[[6](#ref-williams2016)\] can be added to make
the output data more biologically realistic. The key parameters
`number_of_generations`, `division_probability` and `death_probability`
determine the final size of the tumour, while the degree of intratumour heterogeneity can 
be varied by changing the `mutation_probability` parameter. 
For neutral tumour evolution, parameter `adv_mutant_division_probability`
and `adv_mutant_death_probability` must be the same as `division_probability`
and `death_probability`.

Throughout the cancer growth phase, `CancerSim` stores information about
the parent cell and a designation of newly acquired mutations for every
cell. Complete mutational profiles of cells are reconstructed a
posteriori based on the stored lineage information.

The division rules which allow only cells with empty neighbouring nodes
to divide, cause exclusively peripheral growth and complete absence of
dynamics in the tumour centre. To allow for variable degree of growth
inside the tumour, we introduced a death process. At every time step,
after all cells attempt their division, a number of random cells die according
to `death_probability` and `adv_mutant_death_probability` and
yield their position to host a new cancer cell in a subsequent time
step.

### Simulation results
After the simulation, the tumour matrix, and the lists of lineages and
frequencies of each mutation in the tumour are exported to files.
Furthermore, the virtual tumour can be sampled and a histogram over the
frequency of mutations will be visualised. Alternatively, a saved tumour
can be loaded from file and then be subjected to the sampling process.


Installation 
------------
`CancerSim` is written in Python (version \>3.5). We recommend to install
it directly from the source code. To download the code:

**EITHER** clone the repository:

    $> git clone https://github.com/mpievolbio-scicomp/cancer_sim.git

**OR** download the source code archive:

    $> wget https://github.com/mpievolbio-scicomp/cancer_sim/archive/master.zip
    $> unzip master.zip -d cancer_sim

Change into the source code directory

    $> cd cancer_sim

We provide for two alternatives to install the software after it was
downloaded:

### Alternative 1: `conda`
#### New conda environment
We [provide](provide) an `environment.yml` to be consumed by `conda`. To create a
fully self-contained conda environment (named `casim`):

    $> conda env create -n casim --file environment.yml

This will also install the cancer simulation code into the new
environment.

To activate the new conda environment:

    $> source activate casim

or

    $> conda activate casim

if you have set up conda appropriately.

#### Install into existing and activated conda environment
To install the software into an already existing environment:

    $> conda activate <name_of_existing_conda_environment>
    $> conda env update --file environment.yml

### Alternative 2: `pip`
The file `requirements.txt` is meant to be consumed by `pip`:

    $> pip install -r requirements.txt [--user]

The option `--user` is needed to install without admin privileges.

### Installed module
After installation, `CancerSim` is available in python as the `casim` module.
E.g. in a python script, one would import the module as:

    >>> from casim import casim

Testing
-------
Although not strictly required, we recommend to run the test suite after
installation. Simply execute the `run_tests.sh` shell script:

    $> ./run_tests.sh

This will generate a test log named `casim_test@<timestamp>.log` with
`<timestamp>` being the date and time when the test was run. You should
see an `OK` at the bottom of the log. If instead errors or failures are
reported, something is wrong with the installation or the code itself.
Feel free to open a github issue at
<https://github.com/mpievolbio-scicomp/cancer_sim/issues> and attach the
test log plus any information that may be useful to reproduce the error
(version hash, computer hardware, operating system, python version, a
dump of `conda env export` if applicable).

The test suite is automatically run after each commit to the code base.
Results are published on
[travis-ci.org](https://travis-ci.org/mpievolbio-scicomp/cancer_sim).

High-level functionality
-------------------------
### Setting the simulation parameters
The parameters of the cancer simulation are specified in a python module or
programmatically via the `CancerSimulationParameters` class. A
documented example `params.py` is included in the source code (under
`casim/params.py`) and reproduced here:

    $> cat casim/params.py
    ################################################################################
    #                                                                              #
    # Commented casim parameter input file.                                        #
    # Valid settings are indicated in parentheses at the end of each comment line. #
    # [0,1] stands for the closed interval from 0 to 1, including the limits; ||   #
    # means "or".                                                                  #
    #                                                                              #
    ################################################################################
    
    # Number of mesh points in each dimension (>0)
    matrix_size = 1000

    # Number of generations to simulate (>0).
    number_of_generations = 20

    # Probability of cell division per generation ([0,1]).
    division_probability = 1

    # Probability of division for cells with advantageous mutation ([0,1]).
    adv_mutant_division_probability = 1

    # Fraction of cells that die per generation ([0,1]).
    death_probability = 0.1

    # Fraction of cells with advantageous mutation that die per generation ([0,1]).
    adv_mutant_death_probability = 0.0

    # Probability of mutations ([0,1]).
    mutation_probability = 1

    # Mutation probability for the adv. cells ([0,1]).
    adv_mutant_mutation_probability = 1

    # Number of mutations per cell division (>=0).
    number_of_mutations_per_division = 10

    # Number of generations after which the beneficial mutation occurs (>=1).
    adv_mutation_wait_time = 10

    # Number of mutations present in first cancer cell (>=0).
    number_of_initial_mutations = 150

    # Tumour multiplicity ("single" || "double").
    tumour_multiplicity = "single"

    # Sequencing read depth (read length * number of reads / genome length).
    read_depth = 100

    # Fraction of cells to be sampled ([0,1]).
    sampling_fraction = 0.1
        
	# Sampling position (list of (x,y) coordinates in the range [0,matrix_size-1]).
	# If left blank or None, random position will be chosen.
	# sampling_positions = None # This will randomly set a single sampling position.
    sampling_positions = [(500,500),(490,490)]
    
    # Plot the tumour growth curve (True || False).
    plot_tumour_growth = True
        
    # Export the tumour growth data to file (True || False).
    export_tumour = True

Here, we simulate a single 2D tumour on a 1000x1000 grid (`matrix_size=1000`) for a total of
20 generations (`number_of_generations=20`).
On average, both healthy and mutant cells divide once per generation
(`division_probability`). The first cancer cell carries 150 mutations
(`number_of_initial_mutations=150`); both healthy and mutant cells aquire 10 new
mutations (`number_of_mutations_per_division=10`) in
each generation with a certainty of 100% (`mutation_probability=1`). The advantageous mutation happens in
the 10th generation (`adv_mutation_wait_time=10`). Mutant cells with advantageous mutations live on forever
(`adv_mutant_death_probability=0`) while healthy cells die with a rate of 0.1
per generation (`death_probability=0.1`).

After completion of the last generation, two spatial samples are taken, one from the tumour center and one from a
slightly more lateral
position (`sampling_positions = [(500,500),(490,490)]`). Each sample contains 10% (with respect to the whole tumour size)
 closely positioned tumour cells
(`sampling_fraction=0.1`). The samples are subject to genetic sequencing with a read depth
of 100 (`read_depth=100`). The data is written to disk (`export_tumour=True`)
and plots showing the mutation histograms for the whole tumour as well as for the sampled part of the
tumour are generated. Furthermore, a plot showing the tumour growth over time is
saved (`plot_tumour_growth=True`).
 
Users should start with the template
and adjust the parameters as needed for
their application by setting experimentally or theoretically known values or by
calibrating the simulation output against experiments or other models.

### Run the example
The simulation is started from the command line. The syntax is

    $> python -m casim.casim [-h] [-s SEED] [-p PARAMS] [-o DIR]

`SEED` is the random seed. Using the
same seed in two simulation runs with identical parameters results in
identical results. If not given, `SEED` defaults to 1. `PARAMS` should point to
a python parameter file. If not given, it defaults to `params.py` in the current
working directory. If that file does not exist, default parameters are assumed.
`DIR` specifies the directory where to store the
simulation log and output data. If not given, output will be stored in
the directory `casim_out` in the current directory. For each seed, a
subdirectory `cancer_SEED` will be created. If that subdirectory already
exists because an earlier run used the same seed, the run will abort.
This is a safety catch to avoid overwriting data from previous runs.

### Output
After the run has finished, you should find the results in the specified output
directory:

    $> ls out/cancer_1/simOutput  
    growthCurve.pdf  mut_container.p              sample_out_490_490.txt
    mtx.p            sampleHistogram_490_490.pdf  sample_out_500_500.txt
    mtx_VAF.txt      sampleHistogram_500_500.pdf  wholeTumourVAFHistogram.pdf   
    
Let's take a look at the `.txt`  files. They contain the simulation output:
`mtx_VAF.txt` is a datafile with three columns: `mutation_id` lists the index of
each primary mutation, `additional_mut_id` indexes the subsequent mutations that occur in a cell of
a given `mutation_id`; `frequency` is the frequency which at a given mutation occurs.

Corresponding to the given sample positions, there is one
`sample_out_XXX_YYY.txt` and one `sampleHistogram_XXX_YYY.pdf` for each
position. The `.txt` filel lists all mutations of the artificial sample
taken from the whole tumour. Columns are identical to `mtx_VAF.txt`.
 
The `.pdf` files are plots of sampled
tumour histogram. `wholeTumourVAFHistogram.pdf` is the histogram for the
complete tumour. You should see figures similar to these:

Whole tumour histogram:  
![Whole tumour histogram](img/example_whole_tumour.png "Whole tumour histogram")  

Central sample histogram:   
![Sampled tumour histogram](img/example_sampled_tumour_pos1.png "Center sample histogram")  

Lateral sample histogram:   
![Sampled tumour histogram](img/example_sampled_tumour_pos2.png "Lateral sample histogram")  

The remaining output files are serialized versions ("pickles") of the tumour
geometry as a 2D matrix (`mtx.p`) and the
mutation list (list of tuples listing the cancer cell index and the mutation ID of each
tumour cell, `mut_container.p`).

### Example notebooks
An example demonstrating how to parametrize the simulation through the
`CancerSimulationParameters` API is provided in the accompanying jupyter notebook at
[quickstart_example.ipynb`](https://github.com/mpievolbio-scicomp/cancer_sim/blob/master/docs/source/include/notebooks/quickstart.ipynb). Launch it in
[binder](https://mybinder.org/v2/gh/mpievolbio-scicomp/cancer_sim.git/master?filepath=docs%2Fsource%2Finclude%2Fnotebooks%2Fquickstart_example.ipynb).

In
[`run_dump_reload_continue.ipynb`](https://github.com/mpievolbio-scicomp/cancer_sim/blob/master/docs/source/include/notebooks/run_dump_reload_continue.ipynb), we
demonstrate how to use the restart capability to modifiy tumour growth parameters in the
middle of a run. In this way, one can model different phases of tumour growth,
e.g. tumour dormancy or onset of cancer therapy.

Community Guidelines
--------------------
As an Open Source project, we welcome all contributions to `CancerSim`. We
recommend the usual github workflow: Fork this repository, commit
your changes and additions to the fork and create a pull request back to the
master branch on this repository. If uncertain about anything, please create an
issue at [https://github.com/mpievolbio-scicomp/cancer_sim/issues](https://github.com/mpievolbio-scicomp/cancer_sim/issues).

Comments, bug reports, or other issues as well as requests
for support should be submitted as a [github
issue](https://github.com/mpievolbio-scicomp/cancer_sim/issues). Please check
the [list of issues](https://github.com/mpievolbio-scicomp/cancer_sim/issues?q=is%3Aissue) if your problem has already been addressed. We will do our best to respond in a timely manner.

References
----------
<a name=ref-kumar2017></a> \[1\] J. C. A. Vinay Kumar Abul K. Abbas,
*Robbins Basic Pathology*, 10th ed. (Elsevier, 2017). ISBN: 9780323353175. 

<a name=ref-turajlic2019></a> \[2\] S. Turajlic, A. Sottoriva, T. Graham,
and C. Swanton, Nat Rev Genet (2019). DOI:
[10.1038/s41576-019-0114-6](https://dx.doi.org/10.1038/s41576-019-0114-6)


<a name=ref-altrock2015></a> \[3\] P. M. Altrock, L. L. Liu, and
F. Michor, Nat Rev Cancer **15**, 730 (2015).
 DOI:
[10.1038/nrc4029](https://dx.doi.org/10.1038/nrc4029)

<a name=ref-opasic2019></a> \[4\] L. Opasic, D. Zhou, B. Werner, D.
Dingli, and A. Traulsen, BMC Cancer **19**, 403 (2019).
 DOI:
[10.1186/s12885-019-5597-1](https://dx.doi.org/10.1186/s12885-019-5597-1)

<a name=ref-waclaw2015></a> \[5\] B. Waclaw, I. Bozic, M. E. Pittman, R. H.  Hruban, B. Vogelstein, and M. A. Nowak, Nature **525**, 261 (2015). DOI: [10.1038/nature14971](https://dx.doi.org/10.1038/nature14971)

<a name=ref-williams2016></a> \[6\] M. J. Williams, B. Werner, C. P. Barnes,
T. A. Graham, and A. Sottoriva, Nature Genetics **48**, 238 (2016).  DOI:
[10.1038/ng.3489](https://dx.doi.org/10.1038/ng.3489)

