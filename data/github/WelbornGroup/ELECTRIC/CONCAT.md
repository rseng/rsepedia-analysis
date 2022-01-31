---
title: 'ELECTRIC: Electric fields Leveraged from multipole Expansion Calculations in Tinker Rapid Interface Code'
tags:
  - Tinker
  - molecular dynamics
  - polarizable force fields
  - electric fields
authors:
  - name: Jessica Nash
    orcid: 0000-0003-1967-5094
    affiliation: "1, 2" 
  - name: Taylor Barnes
    orcid: 0000-0001-9396-1094
    affiliation: "1, 2" 
  - name: Valerie Vaissier Welborn^[Corresponding author]
    orcid: 0000-0003-0834-4441
    affiliation: 2 
affiliations:
 - name: The Molecular Sciences Software Institute, 1880 Pratt Drive, Suite 1100, Blacksburg, VA 24060
   index: 1
 - name: Department of Chemistry, Davidson Hall, Virginia Tech, 1040 Drillfield Drive, Blacksburg, VA 24061
   index: 2
date: July 1, 2020
bibliography: paper.bib
---


# Summary

Polarizable force fields have changed the landscape of biomolecular simulation, mostly by featuring improved electrostatic potential energy terms [@doi:10.1146/annurev-biophys-070317-033349]. These novel energy functions allow environment-driven changes in charge distribution, which yield simulations with improved geometries and molecular properties. In particular, the AMOEBA polarizable force field exhibits two fundamental changes compared to more traditional fixed charge force fields [@doi:10.1021/ct4003702; @doi:10.1021/acs.jctc.7b01169]. The first one relates to permanent electrostatics, expressed in AMOEBA in terms of atomic multipoles (truncated at quadrupoles) that account for anisotropy in the computed charge distributions. The second one represents polarizability through an induced dipole term that can respond to the chemical environment. These modified terms make the AMOEBA force field more physically grounded than other force fields and is the basis for more realistic simulations of biomolecular systems.

Improved electrostatics with AMOEBA give a unique opportunity to accurately compute electric fields, powerful metrics of catalytic activity in enzymes and other systems [@doi:10.1021/jacs.6b12265; @doi:10.1021/acscatal.7b03151; @natcat]. Electric fields projected onto specific bonds report on the effect of the surroundings (interacting via coulombic interactions, solvent effects, hydrogen bonding or other forces, all mostly electrostatic in nature) on the flow of electrons along these bonds. Therefore, projected electric fields are correlated to the probability of breaking these bonds, making them a useful probe of chemical reactivity.

`ELECTRIC` [@electric] is a MolSSI Driver Interface (MDI) [@mdi_repo; @barnes_taylor_arnold_2020_3659285] driver that utilizes Tinker [@doi:10.1021/acs.jctc.8b00529] to analyze specific components of electric fields that are modeled using the AMOEBA force field.  `ELECTRIC` parses Tinker trajectories and orchestrates additional Tinker calculations in order to project components of the electric fields onto user-defined bonds (specified by two atoms). It outputs the field in MV/cm, which is the sum of the direct field (from permanent electrostatics) and the induced field (from the induce dipole term), projected onto the bond unit vector (i.e., normalized by the bond length). `ELECTRIC` enables splitting of the total field into contributions from different components of the system, by molecules or by residues as specified in a reference PDB file. In summary, `ELECTRIC` was designed to expand quantitative system characterization via the computation of electric fields with user-friendly processing tools of Tinker-AMOEBA simulations.

In practice, the user needs a Tinker trajectory file (i.e., `filename.arc`), a Tinker input file (i.e., `tinker.key`) stripped from any keywords corresponding to periodic boundary conditions and a Tinker reference snapshot (`filename.xyz` - this can be the first frame in `filename.arc`) where the dimensions of the box, usually printed on line 2, has been deleted. Tinker is then launched as an MDI engine before the `ELECTRIC` driver. Options regarding the electric fields calculations, such as probe index number, number of frames at the beginning of the trajectory file to skip (equilibration procedure), reference `PDB` file, etc. are specified as command line arguments when launching the driver. A complete usage procedure is provided in the `README` file. 

# Statement of Need
Since electric fields assist the motion of charges such as ions or electrons, they can link structure to function in molecular dynamics simulations. By taking advantage of their additivity property, we can decompose the total electric field into contributions from each system component (protein residues, solvent molecules, etc.). While this approach has already been used to probe electron flow along the bonds that break and form during a catalyzed reaction [@doi:10.1021/jacs.6b12265; @doi:10.1021/acscatal.7b03151], its applicability reaches many other research areas as it can also be used to probe ion transport or electron reorganization upon molecular binding, for example.

# Mathematics
The electric field at atom $i$, $\vec{E}^i$, has components defined as:
\begin{equation}
E^i_x=\sum_jE^{j\to i}_x =\sum_j\left( E^{j\to i}_{x,\text{perm}} +E^{j\to i}_{x,\text{ind}} \right)
\end{equation}
where $E^i_x$ is the $x$-component of the electric field on atom $i$, $E^{j\to i}_x$ is the $x$-component of the electric field on atom $i$ due to atom $j$ and "perm" and "ind" refer to permanent and induced fields, respectively. 

Each atom $i$ is characterized by permanent atomic multipoles, including a monopole (charge) $q^i$, a dipole $\{\mu^i_x,\mu^i_y,\mu^i_z\}$ and a quadrupole $\{Q^i_{xx}, Q^i_{xx},Q^i_{xy},Q^i_{xz}...Q^i_{zz}\}$, such that
\begin{equation}
E^{j\to i}_{x,\text{perm}}=-T_xq^j+\sum_{m=y,z}T_{xm}\mu^j_m-\frac{1}{3}\sum_{m=y,z}\sum_{n=y,z}T_{xmn}Q^j_{mn},
\end{equation}
and
\begin{equation}
E^{j\to i}_{x,\text{ind}}=\sum_{m=y,z}=T_{xm}\mu^j_{\text{ind},m},
\end{equation}
where
\begin{equation}
T_{xy...}=\frac{1}{4\pi\epsilon_0}\nabla_x\nabla_y...\frac{1}{r_{ij}},
\end{equation}
as also defined in @doi:10.1021/jacs.6b12265.

Similar equations can be written for the $y$- and $z$-components of the field. 
The field projected onto a specific bond, say bond $ij$ between atoms $i$ and $j$, is then calculated as:  
\begin{equation}
E^{ij}_\text{proj}=\left( \frac{\vec{E}^i+\vec{E}^j}{2}\right).\vec{u}^{ij},
\end{equation}
where $E^{ij}_\text{proj}$ is the electric field projected onto bond $ij$ and $\vec{u}^{ij}$ the unitary vector defining bond $ij$.

# Acknowledgements

The authors thank Yi Zheng for producing AMOEBA Tinker trajectories to test `ELECTRIC` and the Virginia Tech Department Faculty Start-up Funds for financial support.

# References
[![Build Status](https://travis-ci.com/WelbornGroup/ELECTRIC.svg?branch=master)](https://travis-ci.com/WelbornGroup/ELECTRIC)
[![codecov](https://codecov.io/gh/WelbornGroup/ELECTRIC/branch/master/graph/badge.svg)](https://codecov.io/gh/WelbornGroup/ELECTRIC)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/WelbornGroup/ELECTRIC.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/WelbornGroup/ELECTRIC/context:python)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02576/status.svg)](https://doi.org/10.21105/joss.02576)

# ELECTRIC

The `ELECTRIC` driver has been written by Jessica Nash and Taylor Barnes, from the [MolSSI](https://molssi.org/).  Read the documentation [here](https://welborngroup.github.io/ELECTRIC/).


## Overview

This repository contains a driver that uses the [MolSSI Driver Interface](https://github.com/MolSSI-MDI/MDI_Library) to perform electric field analysis of [Tinker](https://dasher.wustl.edu/tinker/) trajectories which use the AMOEBA forcefield. This currently works as a post-processing tool, meaning that you run simulations as normal using Tinker, then analyze the trajectories using
MDI-enabled Tinker and this driver.

Using this tool, you can calculate the electric field along a bond or between atoms due to molecules or residues in the system.

### Compiling MDI-Tinker and ELECTRIC

Installation of ELECTRIC and MDI-enabled Tinker are bundled in one convenient build script. 

To install ELECTRIC and MDI-enabled Tinker, you should have cmake and a fortran compiler installed. Then, you can download and build ELECTRIC and MDI-enabled Tinker using the following command in your terminal. Make sure you are in the directory where you want your ELECTRIC driver to be. You should note this location, because you will need to specify the path to some files built during this process in order to perform analysis.

```
git clone --recurse-submodules https://github.com/WelbornGroup/ELECTRIC.git
cd ELECTRIC
./build.sh
```

This will download and build ELECTRIC and MDI-enabled Tinker.

In certain environments, it may be necessary to manually set the compilers.
This can be done when calling the `build.sh` script.
For example, to compile on NERSC's Cori system, you can do:
```
module load PrgEnv-gnu
git clone --recurse-submodules https://github.com/WelbornGroup/ELECTRIC.git
cd ELECTRIC
CC=cc FC=ftn ./build.sh
```

Upon successfull building, you will have the ELECTRIC driver in ELECTRIC/ELECTRIC/ELECTRIC.py, and the needed Tinker executable (dynamic.x) in ELECTRIC/modules/Tinker/build/tinker/source/dynamic.x . The location of these files can be found in text files in ELECTRIC/test/locations/ELECTRIC and ELECTRIC/test/locations/Tinker_ELECTRIC. You will need these for using ELECTRIC.

### Python Dependencies

In order to run ELECTRIC, you will need to be in a python environment which has numpy and pandas installed. We recommend installing these packages in a conda environment created for ELECTRIC analysis.

``` 
conda install -c conda-forge numpy pandas
```

## Testing

You can now run a quick test of the driver by changing directory to the `ELECTRIC/test/bench5` directory and running the `tcp.sh` script:

    ./tcp.sh

This script will run a short Tinker dynamics simulation that includes periodic boundary conditions. This command is on line 20 of the provided file. This is a standard Tinker call, as you would normally run a simulation. If you are performing post processing on a simulation, you will not use this line.

```
${TINKER_LOC} bench5 -k bench5.key 10 1.0 0.001999 2 300.00 > Dynamics.log
```

The script then launches an instance of Tinker as an MDI engine, which will request a connection to the driver and then listen for commands from the driver. This command is similar to running a simulation with Tinker, except that it uses a modified Tinker input file (more on this below), and adds an additional command line argument which passes information to MDI (`-mdi "role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost"`):

```
{TINKER_LOC} bench5 -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost" 10 1.0 0.001999 2 300.00 > no_ewald.log &
```

The script will then launch an instance of the driver in the background, which will listen for connections from an MDI engine:

```
python ${DRIVER_LOC} -probes "1 40" -snap bench5.arc -mdi "-role DRIVER -name driver -method TCP -port 8022" --bymol &
```

The driver's output should match the reference output file (`proj_totfield.csv`) in the `sample_analysis` directory.

## Usage

In general, running a calculation with the driver requires the following steps:

1. **Run a dynamics simulation with Tinker.**  
This simulation should be run with periodic boundary conditions (if desired), and should print snapshots of its results to a single file (i.e., `coordinates.arc`).
If each snapshot was instead written to a different file (i.e., `coordinates.001`, `coordinates.002`, etc.) then you may concatenate them into a single file.

2. **Create a new Tinker keyfile.**   
This keyfile should be identical to the one used in Step 1, except that it **must not** include periodic boundary conditions and **must not** use an Ewald summation. This means that in the `.key` file for running the driver, you should not have an `a-axis` keyword, or keywords related to Ewald.

3. **Launch one (or more; see the `--nengines` option below) instance(s) of Tinker as an MDI engine, using the keyfile created in Step 2.**  
This is done in the same way you launch a normal Tinker simulation (by launching the `dynamic.x` executable) except that the `-mdi` command-line option is added. However, it is **very important** that the reference coordinates you use do not have periodic boundary information. So, if when you originally ran the simulation you started it with a snapshot from a previous simulation run, make sure to create a new snapshot to launch the simulation from which does not include box information on line 2.

  The argument to the `-mdi` command-line option details how Tinker should connect to the driver; its possible arguments are described in the [MDI documentation](https://molssi.github.io/MDI_Library/html/library_page.html#library_launching_sec).
  When in doubt, we recommend doing `-mdi "-role ENGINE -name NO_EWALD -method TCP -port 8021 -hostname localhost"`
  When run as an engine, Tinker should be launched in the background; this is done by adding an ampersand (`&`) at the end of the launch line.

4. **Launch the driver.**
The driver accepts a variety of command-line options, which are described in detail below.
One possible launch command would be:

    `python ${DRIVER_LOC} -probes "1 2 10" -snap coordinates.arc -mdi "-role DRIVER -name driver -method TCP -port 8021" --byres ke15.pdb --equil 51 --nengines 15 &`

where `DRIVER_LOC` is the path to ELECTRIC.py which you set during the configuration step.
The output will be written to `proj_totfield.csv`.

It is useful to write a script that performs Steps 3 and 4, especially if the calculations are intended to be run on a shared cluster.
Such a script might look like:

    # location of required codes
    DRIVER_LOC=$(cat ../locations/ELECTRIC)
    TINKER_LOC=$(cat ../locations/Tinker_ELECTRIC)

    # number of instances of Tinker to run as an engine
    nengines=18

    # set the number of threads used by each code
    export OMP_NUM_THREADS=1

    # launch Tinker as an engine
    for i in $( eval echo {1..$nengines} )
    do
    ${TINKER_LOC} coordinates.in -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8021 -hostname localhost" 10 1.0 1.0 2 300 > no_ewald${i}.log &
    done

    # launch the driver
    python ${DRIVER_LOC} -probes "32 33 59 60" -snap coordinates.arc -mdi "-role DRIVER -name driver -method TCP -port 8021" --byres ke15.pdb --equil 51 --nengines ${nengines} &

    wait

### Using MPI for communication

In addition to the above examples of using TCP/IP sockets for communication between ELECTRIC and Tinker, it is also possible to establish communication between the codes using the Message Passing Interface (MPI).
Information about launching codes using MPI communication is available [here](https://molssi-mdi.github.io/MDI_Library/html/library_page.html#library_launching_sec).
This approach is likely to be preferable when running on large supercomputing clusters.

Note that on systems that manage MPI jobs using SLURM, it is necessary to use `srun` to launch jobs rather than direct calls to `mpiexec`.
Example scripts for launching on NERSC's Cori system are provided at `ELECTRIC/test/bench5/cori.sh` (for running with a single instance of Tinker) and `ELECTRIC/test/bench5/cori5.sh` (for running with multiple instances of Tinker).
Before running one of these scripts, you will need to change the `--account` SBATCH option to your own account.




## Command-Line Options

You can see command line arguments for this driver using the following command from the top level of this repositry:

    python ELECTRIC/ELECTRIC.py --help

Here is the help information for the command line arguments:

    usage: ELECTRIC.py [-h] -mdi MDI -snap SNAP -probes PROBES
                              [--nengines NENGINES] [--equil EQUIL]
                              [--stride STRIDE] [--byres BYRES] [--bymol]

    required arguments:
      -mdi MDI            flags for mdi (default: None)
      -snap SNAP          The file name of the trajectory to analyze. (default:
                          None)
      -probes PROBES      Atom indices which are probes for the electric field
                          calculations. For example, if you would like to
                          calculate the electric field along the bond between
                          atoms 1 and 2, you would use -probes "1 2". (default:
                          None)

    optional arguments:
      -h, --help          show this help message and exit
      --nengines NENGINES This option allows the driver to farm tasks out to
                          multiple Tinker engines simultaneously, enabling
                          parallelization of the electric field analysis
                          computation. The argument to this option **must** be
                          equal to the number of Tinker engines that are launched
                          along with the driver. (default: 1)
      --equil EQUIL       The number of frames to skip performing analysis on at
                          the beginning of the trajectory file (given by the -snap
                          argument) For example, using --equil 50 will result in
                          the first 50 frames of the trajectory being skipped.
                          (default: 0)
      --stride STRIDE     The number of frames to skip between analysis
                          calculations. For example, using --stride 2 would result
                          in analysis of every other frame in the trajectory.
                          (default: 1)
      --byres BYRES       Flag which indicates electric field at the probe atoms
                          should be calculated with electric field contributions
                          given per residue. If --byres is indicated, the argument
                          should be followed by the filename for a pdb file which
                          gives residues. (default: None)
      --bymol             Flag which indicates electric field at the probe atoms
                          should be calculated with electric field contributions
                          given per molecule. (default: False)

## Output

The driver will output a file called `proj_totfield.csv`. This is a CSV file which contains data on the projected electric field at the point between each probe atom due to each fragment , depending on input (`--byres` for by residue, `--bymol` for by molecule, or by atom if neither argument is given.). Each column will contain a header which indicates which probe atoms the measurement is between, followed by the frame number, while the rows will be the electric field at the mean location between the probe atoms due to a particular fragment

Consider the example (`bench5`), which was run with the following command:

    python ${DRIVER_LOC} -probes "1 40" -snap bench5.arc -mdi "-role DRIVER -name driver -method TCP -port 8022" --bymol

Here, we have set the probe atoms to be atoms 1 and 40, and we have indicated we want the the electric field between the probe atoms based on contributions by molecule. Headers will be "`i and j - frame n`. Where `i` and `j` are the atom indices of the probes, and `n` is the frame number.

For the example, headers are:

    "1 and 40 - frame 0"
    "1 and 40 - frame 1"
    "1 and 40 - frame 2"
    "1 and 40 - frame 3"
    "1 and 40 - frame 4"

Since this calculation was run using `--bymol`, there are 216 rows, one for each molecule in the system.

The first entry, column `1 and 40 - frame 0`, header `molecule 1`, gives the projected total electric field at the midway point between `atom 1` and `atom 40` due to `molecule 1`. The electric field has been projected along the vector which points from `atom 1` to `atom 40`. The projection will always be along the vector from atom 1 to atom 2. You can reverse the sign of the number if you would like the vector to point the opposite way.

A sample script which calculates the time average for each probe pair is given in the directory `sample_analysis`.
# How to contribute

We welcome contributions from external contributors, and this document
describes how to merge code changes into ELECTRIC. 

## Getting Started

* Make sure you have a [GitHub account](https://github.com/signup/free).
* [Fork](https://help.github.com/articles/fork-a-repo/) this repository on GitHub.
* On your local machine,
  [clone](https://help.github.com/articles/cloning-a-repository/) your fork of
  the repository.

## Making Changes

* Add some really awesome code to your local fork.  It's usually a [good
  idea](http://blog.jasonmeridth.com/posts/do-not-issue-pull-requests-from-your-master-branch/)
  to make changes on a
  [branch](https://help.github.com/articles/creating-and-deleting-branches-within-your-repository/)
  with the branch name relating to the feature you are going to add.
* When you are ready for others to examine and comment on your new feature,
  navigate to your fork of ELECTRIC on GitHub and open a [pull
  request](https://help.github.com/articles/using-pull-requests/) (PR). Note that
  after you launch a PR from one of your fork's branches, all
  subsequent commits to that branch will be added to the open pull request
  automatically.  Each commit added to the PR will be validated for
  mergability, compilation and test suite compliance; the results of these tests
  will be visible on the PR page.
* If you're providing a new feature, you must add test cases and documentation.
* Please also make sure to run the code formattter [black](https://black.readthedocs.io/en/stable/) 
  on your code before submitting your pull request.
* When the code is ready to go, make sure you run the test suite using pytest.
* When you're ready to be considered for merging, check the "Ready to go"
  box on the PR page to let the ELECTRIC devs know that the changes are complete.
  The code will not be merged until this box is checked, the continuous
  integration returns checkmarks,
  and multiple core developers give "Approved" reviews.

# Additional Resources

* [General GitHub documentation](https://help.github.com/)
* [PR best practices](http://codeinthehole.com/writing/pull-requests-and-other-good-practices-for-teams-using-github/)
* [A guide to contributing to software packages](http://www.contribution-guide.org)
* [Thinkful PR example](http://www.thinkful.com/learn/github-pull-request-tutorial/#Time-to-Submit-Your-First-PR)# Sample Analysis

This directory contains a sample script to process the output of the driver. The driver will write a file called `proj_totfield.csv`. For an in depth explanation of this file, please see the `README.md` file in the top level of this repository.

The file `calculate_average.py` will calculate the time average projected field per fragment. This script should work to calculate the time average for any `proj_totfield.csv` output by the driver.

This script, will automatically look for a file called `proj_totfield.csv` in the current directory. To run it:

    python calculate_average.py

If you would like to change the file name, add the file name after the script name when you run.

    python calculate_average.py -filename FILENAME

The output will be a file for each pairwise probe interaction. The file gives the time average value of the projected field and the standard deviation.
MolSSI Driver Interface (MDI) Library
=====================================

[![Build Status](https://travis-ci.org/MolSSI-MDI/MDI_Library.svg?branch=master)](https://travis-ci.org/MolSSI-MDI/MDI_Library)
[![Build Status](https://dev.azure.com/taylorabarnes/MDI_Library/_apis/build/status/MolSSI.MDI_Library?branchName=master)](https://dev.azure.com/taylorabarnes/MDI_Library/_build/latest?definitionId=1&branchName=master)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/MolSSI-MDI/MDI_Library.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MolSSI-MDI/MDI_Library/context:python)
[![Language grade: C/C++](https://img.shields.io/lgtm/grade/cpp/g/MolSSI-MDI/MDI_Library.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/MolSSI-MDI/MDI_Library/context:cpp)
[![codecov](https://codecov.io/gh/MolSSI-MDI/MDI_Library/branch/master/graph/badge.svg)](https://codecov.io/gh/MolSSI-MDI/MDI_Library)

## Overview

The MolSSI Driver Interface (MDI) project provides a standardized API for fast, on-the-fly communication between computational chemistry codes.  This greatly simplifies the process of implementing methods that require the cooperation of multiple software packages and enables developers to write a single implementation that works across many different codes.  The API is sufficiently general to support a wide variety of techniques, including QM/MM, ab initio MD, machine learning, advanced sampling, and path integral MD, while also being straightforwardly extensible.  Communication between codes is handled by the MDI Library, which enables tight coupling between codes using either the MPI or TCP/IP methods.

## Documentation

Complete documentation can be found at https://molssi-mdi.github.io/MDI_Library

## License

The MDI Library is released under the BSD 3-clause license. See LICENSE for details.

## Acknowledgements

This work was supported by the Molecular Sciences Software Institute under U.S. National Science Foundation grant ACI-1547580.

MDI builds upon the work of numerous development groups, without whom it could not exist.
The syntactical structure of the MDI Standard, including the fundamental command-response communication pattern, is modelled after that used by the <a href="http://ipi-code.org/">i-PI</a> project, as is the string-based representation of commands.
The Node System draws inspiration from the techniques used by several molecular mechanics packages, especially <a href="https://lammps.sandia.gov/">LAMMPS</a> and <a href="http://openmm.org/">OpenMM</a>, to enable modular code additions.
The unit conversions available through the MDI Library were provided by the <a href="https://github.com/MolSSI/QCElemental">QCElemental</a> project.
Certain details of the communication protocols implemented by the MDI Library, especially pertaining to MPI-based communication in the MPMD regime, were informed by the accomplishments of the <a href="https://cslib.sandia.gov/">CSlib</a> library.
The library-based communication protocol was developed in response to discussions with the <a href="https://gitlab.com/exaalt">EXAALT</a> team.
The interface, error handling, data types, and numerous other elements of the MDI Library are modelled after the <a href="https://www.mpi-forum.org/">MPI Standard</a>.
A distribution of the MDI Library for Python is provided by <a href="https://conda-forge.org/">Conda Forge</a>.
Usage
=====

Information Flow
----------------

ELECTRIC is a post-processing tool for simulations running with the AMOEBA polarizable force field using the Tinker software package.

.. image:: images/inputs_and_outputs.svg
   :width: 600

Procedure
---------

In general, running a calculation with the driver requires the following steps:

1. **Run a dynamics simulation with Tinker.**  
This simulation should be run with periodic boundary conditions (if desired), and should print snapshots of its results to a single file (i.e., :code:`coordinates.arc`).
If each snapshot was instead written to a different file (i.e., :code:`coordinates.001`, :code:`coordinates.002`, etc.) then you may concatenate them into a single file.

2. **Create a new Tinker keyfile.**   
This keyfile should be identical to the one used in Step 1, except that it **must not** include periodic boundary conditions and **must not** use an Ewald summation. This means that in the :code:`.key` file for running the driver, you should not have an :code:`a-axis` keyword, or keywords related to Ewald.

3. **Launch one (or more; see the `-nengines` option below) instance(s) of Tinker as an MDI engine, using the keyfile created in Step 2.**  
This is done in the same way you launch a normal Tinker simulation (by launching the :code:`dynamic.x` executable) except that the :code:`-mdi` command-line option is added. However, it is **very important** that the reference coordinates you use do not have periodic boundary information. So, if when you originally ran the simulation you started it with a snapshot from a previous simulation run, make sure to create a new snapshot to launch the simulation from which does not include box information on line 2.

The argument to the :code:`-mdi` command-line option details how Tinker should connect to the driver; its possible arguments are described in the `MDI documentation`_ .
When in doubt, we recommend doing :code:`-mdi "-role ENGINE -name NO_EWALD -method TCP -port 8021 -hostname localhost"`
When run as an engine, Tinker should be launched in the background; this is done by adding an ampersand (:code:`&`) at the end of the launch line.

4. **Launch the driver**
The driver accepts a variety of command-line options, which are described in detail below.
One possible launch command would be:

.. code-block:: bash

    `python ${DRIVER_LOC} -probes "1 2 10" -snap coordinates.arc -mdi "-role DRIVER -name driver -method TCP -port 8021" --byres ke15.pdb --equil 51 --nengines 15 &`

where `DRIVER_LOC` is the path to ELECTRIC.py which you set during the configuration step. See the section :ref:`electric settings` for a detailed explanation of command line arguments for ELECTRIC.

The output will be written to `proj_totfield.csv`.

It is useful to write a script that performs Steps 3 and 4, especially if the calculations are intended to be run on a shared cluster.
Such a script might look like:

.. _example:

Example Script
^^^^^^^^^^^^^^

.. code-block:: bash
    :linenos:

    # location of required codes
    DRIVER_LOC=$(cat ../locations/ELECTRIC)
    TINKER_LOC=$(cat ../locations/Tinker_ELECTRIC)

    # number of instances of Tinker to run as an engine
    nengines=18

    # set the number of threads used by each code
    export OMP_NUM_THREADS=1

    # launch Tinker as an engine
    for i in $( eval echo {1..$nengines} )
    do
    ${TINKER_LOC} coordinates.in -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8021 -hostname localhost" 10 1.0 1.0 2 300 > no_ewald${i}.log &
    done

    # launch the driver
    python ${DRIVER_LOC} -probes "32 33 59 60" -snap coordinates.arc -mdi "-role DRIVER -name driver -method TCP -port 8021" --byres ke15.pdb --equil 51 --nengines ${nengines} &

    wait

You can read more below, or you can try out the tutorial_ to run a calculation yourself.


.. _electric settings:

ELECTRIC Calculation Settings
-----------------------------

You can change the options for your electric calculation through command line arguments. 

.. argparse::
   :filename: ../ELECTRIC/util.py
   :func: create_parser
   :prog: python ELECTRIC.py


Output
------

The driver will output a file called :code:`proj_totfield.csv`. This is a CSV file which contains data on the projected electric field at the point between each probe atom due to each fragment , depending on input (`--byres` for by residue, `--bymol` for by molecule, or by atom if neither argument is given.). Each column will contain a header which indicates which probe atoms the measurement is between, followed by the frame number, while the rows will be the electric field at the mean location between the probe atoms due to a particular fragment

Consider the example (:code:`bench5`), which was run with the following command:

.. code-block:: bash

    python ${DRIVER_LOC} -probes "1 40" -snap bench5.arc -mdi "-role DRIVER -name driver -method TCP -port 8022" --bymol

Here, we have set the probe atoms to be atoms 1 and 40, and we have indicated we want the the electric field between the probe atoms based on contributions by molecule. Headers will be "`i and j - frame n`. Where `i` and `j` are the atom indices of the probes, and `n` is the frame number.

For the example, headers are:

.. code-block:: text

    "1 and 40 - frame 1"
    "1 and 40 - frame 2"
    "1 and 40 - frame 3"
    "1 and 40 - frame 4"
    "1 and 40 - frame 5"

Since this calculation was run using :code:`--bymol`, there are 216 rows, one for each molecule in the system.

The first entry, column :code:`1 and 40 - frame 1`, header :code:`molecule 1`, gives the projected total electric field at the midway point between :code:`atom 1` and :code:`atom 40` due to :code:`molecule 1`. The electric field has been projected along the vector which points from :code:`atom 1` to :code:`atom 40`. The projection will always be along the vector from atom 1 to atom 2. You can reverse the sign of the number if you would like the vector to point the opposite way.


Running ELECTRIC in Parallel
-----------------------------

.. note::

    You must have mpi4py installed to run ELECTRIC in parallel. You can install it from conda
    
    .. code-block:: bash

        conda install -c anaconda mpi4py

ELECTRIC is parallelized using MPI4Py. You can take advantage of this parallelization by making sure MPI4Py is installed and using more than one ELECTRIC engine using the :code:`-nengines` command. Note that if you are using the :code:`-nengines` argument with a number greater than one, you must launch the equivalent number of Tinker instances. In the :ref:`example`, this is acheived by setting a variable :code:`nengines` and using this number to launch Tinker instances in a loop (:code:`lines 12-15`) and inputting the same variable into the ELECTRIC launch on :code:`line 18`.

.. warning::

    Launching an unmatching number of MDI-Tinker and ELECTRIC instances will result in your calculation hanging. Make sure that you launch an equivalent number of MDI-Tinker instances to your :code:`-nengines` argument.

.. _tutorial: tutorial.html
.. _`MDI documentation`: https://molssi.github.io/MDI_Library/html/library_page.html#library_launching_sec
Installation
============

Compiling MDI-Tinker and ELECTRIC
----------------------------------

Installation of ELECTRIC and MDI-enabled Tinker are bundled in one convenient build script. 

To install ELECTRIC and MDI-enabled Tinker, you should have cmake and a fortran compiler installed. Then, you can download and build ELECTRIC and MDI-enabled Tinker using the following command in your terminal. Make sure you are in the directory where you want your ELECTRIC driver to be. You should note this location, because you will need to specify the path to some files built during this process in order to perform analysis.

.. code-block:: bash

    git clone --recurse-submodules https://github.com/WelbornGroup/ELECTRIC.git
    cd ELECTRIC
    ./build.sh

This will download and build ELECTRIC and MDI-enabled Tinker. 

Upon successful building, you will have the ELECTRIC driver in ELECTRIC/ELECTRIC/ELECTRIC.py, and the needed Tinker executable (dynamic.x) in ELECTRIC/modules/Tinker/build/tinker/source/dynamic.x . The location of these files can be found in text files in ELECTRIC/test/locations/ELECTRIC and ELECTRIC/test/locations/Tinker_ELECTRIC. You will need these for using ELECTRIC.

Python Dependencies
-------------------

In order to run ELECTRIC, you will need to be in a python environment which has numpy and pandas installed. If you want to run ELECTRIC with more than one engine, you should also install MPI4Py. We recommend installing these packages in a conda environment created for ELECTRIC analysis.

.. code-block:: bash   

    conda install -c anaconda mpi4py
    conda install -c conda-forge numpy pandas

Testing Your Installation
--------------------------

You can now run a quick test of the driver by changing directory to the `ELECTRIC/test/bench5` directory and running the `tcp.sh` script:

.. code-block:: bash

    ./tcp.sh

This script will run a short Tinker dynamics simulation that includes periodic boundary conditions. This command is on line 20 of the provided file. This is a standard Tinker call, as you would normally run a simulation. If you are performing post processing on a simulation, you will not use this line.

.. code-block:: bash

    ${TINKER_LOC} bench5 -k bench5.key 10 1.0 0.001999 2 300.00 > Dynamics.log

The script then launches an instance of Tinker as an MDI engine, which will request a connection to the driver and then listen for commands from the driver. This command is similar to running a simulation with Tinker, except that it uses a modified Tinker input file (more on this below), and adds an additional command line argument which passes information to MDI (`-mdi "role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost"`):

.. code-block:: bash

    ${TINKER_LOC} bench5 -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost" 10 1.0 0.001999 2 300.00 > no_ewald.log &

The script will then launch an instance of the driver in the background, which will listen for connections from an MDI engine:

.. code-block:: bash

    python ${DRIVER_LOC} -probes "1 40" -snap bench5.arc -mdi "-role DRIVER -name driver -method TCP -port 8022" --bymol &

The driver's output should match the reference output file (`proj_totfield.csv`) in the `sample_analysis` directory.




Tutorial
========

This tutorial will walk you through using ELECTRIC to analyze the electric field in a small protein. This tutorial assumes you have ELECTRIC and MDI-enabled Tinker installed. If you don't, navigate to the installation_ instructions.

.. note::
    This tutorial will assume the following:
        - You are able to run a molecular dynamics simulation using the Tinker software, or are familiar enough with molecular dynamics to follow along.
        - You have installed ELECTRIC and MDI-enabled Tinker. If you have not, see the installation_ instructions.
        - You are able to download or clone a directory from git.
        - You are familiar with bash scripts.


The pdb code for this protein is 1l2y_, and you can see the structure below. We have chosen a small protein for demonstrative purposes. The image below shows only the protein, but our simulation is solvated.

.. moleculeView:: 
    
    data-pdb: 1l2y
    data-backgroundcolor: white
    width: 300px
    height: 300px
    data-style: cartoon:color=spectrum

Running an ELECTRIC calculation
###############################

Preparing Files
----------------
To follow along with this tutorial, clone the `tutorial repository`_. Included in this repository is folder called :code:`data`. The :code:`data` directory has all of the data you will need for this tutorial.

ELECTRIC is a post-processing analysis, meaning that you should first run your simulations using the AMOEBA forcefield and save the trajectory. After you have a trajectory from Tinker simulation, use ELECTRIC to perform electric field analysis on that trajectory. We will be analyzing the trajectory included in the tutorial repository, :code:`1l2y_npt.arc`. This trajectory is a text file containing coordinates for a simulation that has already been run.

To get started with ELECTRIC, we will need to prepare our input files. We will need:
    - a molecular dynamics trajectory
    - a Tinker input file (usually called a key file) which does not have settings for periodic boundaries or Ewald summation
    - the forcefield parameter file
    - a bash script file 

We already have the molecular dynamics trajectory, so let's look at each of these additional files.

Simulation Parameter File
^^^^^^^^^^^^^^^^^^^^^^^^^
This file contains the force field parameters for the simulation you have run. If you are using this software for analysis, use the same force field you used to run the molecular dynamics simulation. For this tutorial, the parameter file is :code:`amoebabio18.prm`. We will need this for our Tinker input file (next section).

Tinker input file
^^^^^^^^^^^^^^^^^
Next, we must prepare an input file which tells Tinker settings for our calculation. This input file should be a modified version of the one which you used to run your initial simulation. Consider the input file, :code:`tinker.key` used to obtain this trajectory. The parameter file in the previous step is given on :code:`line 1`.

.. code-block:: text
    :linenos:

    parameters amoebabio18.prm 
    openmp-threads 16

    a-axis 50.00 
    b-axis 50.00
    c-axis 50.00

    polar-eps 0.000010
    polar-predict
    polarization mutual


    cutoff 10.0
    ewald
    neighbor-list
    integrator beeman

    thermostat nose-hoover
    barostat nose-hoover

    maxiter 8000
    printout 1000


The input file used for this simulation uses periodic boundaries and an Ewald summation for electrostatics. During a Tinker simulation using AMOEBA, electric fields are evaluated in order to calculate the induced dipoles at each step. In order to get electric field contributions from specific residues, we must calculate the electric field using the real space interactions only (no periodic boundaries or Ewald). 

Remove settings related to cutoffs (:code:`cutoff` keyword), periodic boundaries (:code:`a-axis`, :code:`b-axis`,:code:`c-axis`) and Ewald summation (:code:`ewald`). You can also remove settings having to do with neighbor lists (:code:`neighbor-list`), as they are not needed and can cause an error for this calculation if included.

The modifed input file for ELECTRIC is given below. This file is saved in the data directory with the name :code:`noewald.key`.

.. code-block:: text

    parameters amoebabio18.prm
    openmp-threads 16

    polar-eps 0.000010
    polar-predict
    polarization mutual

    integrator beeman

    thermostat nose-hoover
    barostat nose-hoover

    maxiter 8000
    printout 100


Bash script - run_analysis.sh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
When you run analysis uisng ELECTRIC, ELECTRIC parses your given trajectory sends snapshots to Tinker for electric field calculation. The MDI-enabled version of Tinker then calculates the electric field information for that snapshot. 

You use ELECTRIC from the command line. Consider the following bash script provided for analysis, :code:`run_analysis.sh`. We will explain this script in detail.

.. code-block:: bash
    :linenos:

    #location of required codes
    DRIVER_LOC=LOCATION/TO/ELECTRIC/ELECTRIC.py
    TINKER_LOC=LOCATION/TO/DYNAMIC/dynamic.x

    #remove old files
    if [ -d work ]; then
    rm -r work
    fi

    #create work directory
    cp -r data work
    cd work

    #set the number of threads
    export OMP_NUM_THREADS=2

    #launch MDI enabled Tinker
    ${TINKER_LOC} 1l2y -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost"  10 1.0 0.002 2 300.00 > no_ewald.log &

    #launch driver
    python ${DRIVER_LOC} -snap 1l2y_npt.arc -probes "93 94" -mdi "-role DRIVER -name driver -method TCP -port 8022" --byres 1l2y_solvated.pdb  --equil 120 --stride 2 &

    wait

.. note:: 

    For this tutorial, we use the approach of having all data needed for analysis in a directory called `data`. During analysis, we copy everything from :code:`data` into a folder :code:`work`. This part of the tutorial is stylistic. The authors prefer this method to keep files separated, and original files unaltered.

In lines :code:`2` and :code:`3`, you should change the location to your installed ELECTRIC.py file and MDI-enabled :code:`dynamic.x`. Recall from the installation instructions that you can find these in the ELECTRIC directory in the files :code:`ELECTRIC/test/locations/ELECTRIC` and :code:`ELECTRIC/test/locations/Tinker_ELECTRIC`. 

The next section removes the folder called :code:`work` if it exists. This bash script is written to put all analysis files into a folder called :code:`work` to keep our original files clean. 

MDI-enabled Tinker is launched on :code:`line 18` with the command

.. code-block:: bash

    ${TINKER_LOC} 1l2y -k no_ewald.key -mdi "-role ENGINE -name NO_EWALD -method TCP -port 8022 -hostname localhost"  10 1.0 0.002 2 300.00 > no_ewald.log &

The first thing on this line, :code:`${TINKER_LOC}` fills in the location for :code:`dynamic.x` which you put in line 2. Next, `1l2y` is the file name (without an extension) of the xyz file for this calculation (provided vile :code:`12ly.xyz`). You should have this from your original simulation. However, make sure that there is no box information on line two of this :code:`xyz` file, as this could cause Tinker to use periodic boundaries. Next, we give the input file (key file) we have prepared in the previous step using :code:`-k noewald.key`. Then, we give our MDI options. The given options should work for most analysis. After the MDI options are some Tinker input options. For our analysis, it will not really matter what we put here since we are running calculations on one snapshot at a time. However, you must have these present for Tinker to run. Very importantly, note the ampersand (:code:`&`) at the end of this line. This will launch Tinker in the background, where it will be waiting for commands from ELECTRIC.

.. warning::
    
    Make sure that there is no box information on line two of the :code:`xyz` file used to launch MDI-enabled Tinker. This could cause Tinker to use periodic boundaries.

In the next command (:code:`line 21`), we launch ELECTRIC.

.. code-block:: bash   

    python ${DRIVER_LOC} -snap 1l2y_npt.arc -probes "78 93 94"  -mdi "-role DRIVER -name driver -method TCP -port 8022" --byres 1l2y_solvated.pdb  --equil 120 --stride 2 &

Here, we first give the location of our ELECTRIC driver. We indicate our trajectory file using the `-snap` argument with the filename to analyze, followed by MDI options.

Probe Atoms 
++++++++++++

To run an ELECTRIC calculation, you must give the indices of your probe atoms. The probe atoms are the atoms which are used as 'probes' for the electric field. ELECTRIC reports the projected total electric field at the midpoint between all probe atom pairs. This allows you to calculate electric fields along bonds `as reported in literature <https://pubs.acs.org/doi/10.1021/jacs.9b05323>`_.

You should obtain the number of the probe atoms from the :code:`xyz` file you use to launch MDI-enabled Tinker. Note that the index you use here should match the number given in the first column of your xyz file. The projection of the electric field at the midpoint of these two atoms will be reported for each analyzed frame. If you indicate more than two probes, all pairwise fields will be reported (ie, if using "78 93 94", you will get "78 and 93", "78 and 94" and "93 and 94"). You can see the atoms we have chosen as probes highlighted below:

.. moleculeView:: 
    
    data-pdb: 1l2y
    data-backgroundcolor: 0xffffff
    width: 300px
    height: 300px
    data-style: cartoon:color=spectrum
    data-select1: serial:78,93,94
    data-style1: sphere

The argument `--byres` gives information to ELECTRIC about how we would like the electric field reported. When we use the :code:`--byres` argument, it should be followed by a pdb which contains residue information for the system you are studying. When using this argument, electric field contributions from each residue will be reported. Other options are :code:`--byatom` top report electric field contributions from each atom, and :code:`--bymol` to report electric field contributions from each molecule. 

When using :code:`--byres`, solvent should be at the end of the pdb and xyz files. Solvent (ions and water) will be grouped together into a single residue.

.. warning::

    When using the :code:`byres` option, you should verify that the residues in your pdb file match what you expect for your xyz file. You can do this with the utility function :code:`residue_report.py`. ELECTRIC will check that the :code:`xyz` and :code:`pdb` have the same number of atoms. However, all residue information will come from the PDB, so make sure the residue information in your provided PDB is as you expect.

.. note::

    The utility script :code:`residue_report.py` is provided in the same directory as :code:`ELECTRIC.py`. To use it,

    .. code-block:: bash

        python residue_report.py PDB_FILENAME

    This will output a report which gives the residue number, the atom index on which the residue starts and the residue name. When using :code:`--byres`, you should first verify that your pdb file has residues defined as you want and matches your xyz file and trajectory. ELECTRIC only checks that the pdb and xyz file have the same number of atoms, it does not check atom identity or order. For this tutorial, our output is

    .. code-block:: text

        Found 12199 atoms and 21 residues.
        Residue Number       Starting atom        Residue Name        
                1                    1                   ASN         
                2                    17                  LEU         
                3                    36                  TYR         
                4                    57                  ILE         
                5                    76                  GLN         
                6                    93                  TRP         
                7                   117                  LEU         
                8                   136                  LYS         
                9                   158                  ASP         
                10                  170                  GLY         
                11                  177                  GLY         
                12                  184                  PRO         
                13                  198                  SER         
                14                  209                  SER         
                15                  220                  GLY         
                16                  227                  ARG         
                17                  251                  PRO         
                18                  265                  PRO         
                19                  279                  PRO         
                20                  293                  SER         
                21                  305                solvent       


Finally, we give arguments which gives information about the frame we want to analyze. Using `--equil 120` tells ELECTRIC to skip the first 120 frames for analysis, and :code:`--stride 2` tells ELECTRIC to analyze every other frame after 120.

Running the calculation
-----------------------

After you have prepared your files, you can run analysis using the command

.. code-block:: bash

    ./run_analysis.sh > analysis.out &

This will launch ELECTRIC. Again, using the ampersand :code:`&` will run this in the background. Now, you just have to wait for your analysis to finish running.

Analyzing Results from ELECTRIC
###############################

ELECTRIC will output a csv file with the electric field information :code:`proj_totfield.csv` in the :code:`work` folder. Below, we show results (numbers rounded for clarity) for probes 78 and 93 from :code:`proj_totfield.csv`. When these numbers are reported, they are the electric field in Mv/cm projected along the vector pointing from atom 1 to atom 2 due to each residue.

.. datatable::

    csv_file: data/proj_totfield.csv


You are free to analyze this as you like, but we recommend using `pandas`_ to process the csv file. A script to perform averaging of probe pairs across frames is provided in :code:`ELECTRIC/sample_analysis/calculate_average.py`. For example, you can run this script

.. code-block :: bash

    python PATH/TO/calculate_average.py -filename work/proj_totfield.csv

This will output a file with the average projected field for each residue pair. In our case, three files should be output: :code:`78 _and_93.csv`, :code:`78_and_94.csv`, and :code:`93_and_94.csv`. The output for the :code:`78_and_93.csv` is shown in the table below:

.. datatable::

    csv_file: data/78_and_93.csv

.. _1l2y: https://www.rcsb.org/structure/1l2y
.. _installation: installation.html
.. _`tutorial repository`: http://www.github.com/janash/ELECTRIC_tutorial
.. _pandas: https://pandas.pydata.org/
.. ELECTRIC documentation master file, created by
   sphinx-quickstart on Fri Sep 25 17:00:58 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

ELECTRIC
====================================

ELECTRIC: Electric fields Leveraged from multipole Expansion Calculations in Tinker Rapid Interface Code.

ELECTRIC uses the MolSSI Driver Interface to perform electric field analysis of Tinker trajectories which use the AMOEBA forcefield. This currently works as a post-processing tool, meaning that you run simulations as normal using Tinker, then analyze the trajectories using MDI-enabled Tinker and this driver. 

ELECTRIC is written by Jessica A. Nash and Taylor A. Barnes from `The Molecular Sciences Software Institute <https://molssi.org/>`_, in collaboration with `Prof. Valerie Vassier Welborn <https://www.valeriewelborn.com/>`_.

Using this tool, you can calculate the electric field along a bond or between atoms due to molecules or residues in the system.

This method has been reported in the following publications:

- `Computational optimization of electric fields for better catalysis design <https://www.nature.com/articles/s41929-018-0109-2>`_, Nature Catalysis

- `Fluctuations of Electric Fields in the Active Site of the Enzyme Ketosteroid Isomerase <https://pubs.acs.org/doi/10.1021/jacs.9b05323>`_, Journal of the American Chemical Society

- `Computational Optimization of Electric Fields for Improving Catalysis of a Designed Kemp Eliminase <https://pubs.acs.org/doi/10.1021/acscatal.7b03151>`_, ACS Catalysis

You can read about the underlying principles of this analysis in 

- `Computational Design of Synthetic Enzymes <https://pubs.acs.org/doi/10.1021/acs.chemrev.8b00399>`_, Chemical Reviews

ELECTRIC is now available as an open source software package. To get started, head to the installation_ instructions, see the usage_, or try out the tutorial_.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   tutorial


.. _installation: installation.html
.. _usage: usage.html
.. _tutorial: tutorial.html
