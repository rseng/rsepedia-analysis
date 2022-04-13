This repository contains DSOs for Linux in htmd/lib/Linux
Source to these is in the restricted access repo Acellera/htmdlib, only available to Acellera-approved developers.
If you have access to that repo and want to build the DSOs, do the following:

cd [root of htmd checkout]
git clone ssh://git@github.com/acellera/htmdlib --depth 1
htmdlib/C/build.sh $PWD/htmd/lib/Linux/

[![Build Status](https://dev.azure.com/stefdoerr/htmd/_apis/build/status/Acellera.htmd?branchName=master)](https://dev.azure.com/stefdoerr/htmd/_build/latest?definitionId=3&branchName=master)
[![Language Grade: Python](https://img.shields.io/lgtm/grade/python/g/Acellera/htmd.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Acellera/htmd/context:python) 
[![Conda](https://anaconda.org/acellera/htmd/badges/version.svg)](https://anaconda.org/acellera/HTMD)
<!---[![Build status](https://ci.appveyor.com/api/projects/status/m1bxrop34b2qw68x/branch/master?svg=true)](https://ci.appveyor.com/project/acelleraci/htmd/branch/master)--->


# HTMD: Programming Environment for Molecular Discovery
[HTMD](https://www.htmd.org) (acronym for High-Throughput Molecular Dynamics) is a programmable, extensible platform 
written in Python. It provides a complete workspace for simulation-based discovery through molecular simulations while 
aiming to solve the data generation and analysis problem as well as increase reproducibility.

## Licensing
HTMD Community Edition is free to use for non-profit work. Contact Acellera 
[www.acellera.com/contact](https://www.acellera.com/contact/) for information on the full version HTMD Pro or if you need a different license.

## Download HTMD

### Using released versions
HTMD is distributed through conda package manager. The instructions for downloading HTMD can be found in
[https://software.acellera.com/download.html](https://software.acellera.com/download.html). 

### Using this repository 
If you want to use this repository, we recommend to still download a released version of HTMD to have all dependencies 
and then set PYTHONPATH to the git directory.

## HTMD Documentation and User Guide
For HTMD Documentation, please visit: 
[https://software.acellera.com/docs/latest/htmd/api.html](https://software.acellera.com/docs/latest/htmd/api.html).

For a User Guide (easy to start examples), please visit: 
[https://software.acellera.com/docs/latest/htmd/tutorials.html](https://software.acellera.com/docs/latest/htmd/tutorials.html)

## Support and Development

Please report bugs via [GitHub Issues](https://github.org/acellera/htmd/issues).

HTMD is an open-source software and we welcome contributions from the community. For more information on how to 
contribute to HTMD, please visit:
[https://software.acellera.com/docs/latest/htmd/developers/howto.html](https://software.acellera.com/docs/latest/htmd/developers/howto.html)

## Citing HTMD

If you use HTMD in your publication please cite:

Stefan Doerr, Matthew J. Harvey, Frank Noé, and Gianni De Fabritiis. 
**HTMD: High-throughput molecular dynamics for molecular discovery.** 
*Journal of Chemical Theory and Computation*, **2016**, *12* (4), pp 1845–1852.
[doi:10.1021/acs.jctc.6b00049](http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049)


# Notes

When doing a `make diff`, the following differences are expected:

```
Only in str/lipid: toppar_all36_lipid_cholesterol_model_1.str
Only in top: top_all22star_prot.rtf
Only in par: par_all22star_prot.prm
```

The `toppar_all36_lipid_cholesterol_model_1.str` is a processed file created by us to use the right cholesterol. 
The `top_all22star_prot.rtf` and `par_all22star_prot.prm` are files added by us that are not present in the standard 
CHARMM force-field.# Acellera Contributor License ("CL") V1.0

Thank you for your interest in contributing to Acellera, Ltd. ("Acellera"). Acellera is very interested in receiving Your Contribution (defined below). In order to participate, we need to confirm how the rights in Your Contribution will be handled. Following the practices of other open source communities, Acellera requests that you grant Acellera a license, as indicated below, to the intellectual property rights in Your Contributions. Acellera requires that you have a Contributor License ("CL") on file prior to using any of Your Contributions. This helps us ensure that the intellectual property embodied within Acellera Products remains unencumbered for use by the whole of the community. This license is for your protection as a Contributor as well as the protection of Acellera and its users; it does not change your rights to use Your Contributions for any other purpose.

Please read the following document carefully before signing and keep a copy for your records.

Full name: 

Github login: 

E-Mail: 

(optional) Organisation:

## Terms and Conditions

You accept and agree to the following terms and conditions for Your present and future Contributions submitted to Acellera, in consideration for the potential inclusion of Your Contributions in Acellera Products. Except for the license and rights granted herein to Acellera and recipients of software distributed or otherwise made available by Acellera, You reserve all right, title, and interest in and to Your Contributions.

### 1. Definitions

1.1 "You" (or "Your") shall mean the copyright owner or legal entity authorized by the copyright owner that is making this Agreement with Acellera. For legal entities, the entity making a Contribution and all other entities that control, are controlled by, or are under common control with that entity are considered to be a single Contributor. For the purposes of this definition, "control" means (a) the power, direct or indirect, to cause the direction or management of such entity, whether by contract or otherwise, or (b) ownership of fifty percent (50%) or more of the outstanding shares, or (c) beneficial ownership of such entity.

1.2 "Contribution" means any original work of authorship (including software, documentation, or other material), including any modifications or additions to an existing work, that is intentionally submitted by You to Acellera for inclusion in, or documentation of, any of the products owned or managed by Acellera (the "Work"). For the purposes of this definition, "submitted" means any form of electronic, verbal, or written communication sent to Acellera or its representatives, including but not limited to communication on electronic mailing lists, source code control systems, and issue tracking systems that are managed by, or on behalf of, Acellera for the purpose of discussing and improving the Work, but excluding communication that is conspicuously marked or otherwise designated in writing by You as "Not a Contribution."

### 2. Grant of Copyright License.

Subject to the terms and conditions of this Agreement, You hereby grant to Acellera and to recipients of software distributed by Acellera a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable copyright license to reproduce, prepare derivative works of, publicly display, publicly perform, sublicense, and distribute Your Contributions and such derivative works.

### 3. Grant of Patent License.

Subject to the terms and conditions of this Agreement, You hereby grant to Acellera and to recipients of software distributed by Acellera a perpetual, worldwide, non-exclusive, no-charge, royalty-free, irrevocable (except as stated in this section) patent license to make, have made, use, offer to sell, sell, import, and otherwise transfer the Work, where such license applies only to those patent claims existing as of the effective date of this CL and licensable by You that are necessarily infringed by Your Contribution(s) alone or by combination of Your Contribution(s) with the Work to which such Contribution(s) was submitted. If any entity institutes patent litigation against You or any other entity (including a cross-claim or counterclaim in a lawsuit) alleging that your Contribution, or the Work to which you have contributed, constitutes direct or contributory patent infringement, then any patent licenses granted to that entity under this Agreement for that Contribution or Work shall terminate as of the date such litigation is filed.

### 4. Representations

4.1 You represent that you are legally entitled to grant the above license under this CL, whether on behalf of Yourself (if you are an individual person) or on behalf of the entity that You represent (if You are an entity). If You are an individual and Your employer(s) has rights to intellectual property that You create that includes Your Contributions, You represent that You have received permission to make Contributions on behalf of that employer, that Your employer has waived such rights for Your Contributions to Acellera, or that Your employer has executed a separate CL with Acellera.

4.2 You represent that each of Your Contributions is Your original creation (see section 6 for submissions on behalf of others). You represent that Your Contribution submissions include complete details of any third-party license or other restriction (including, but not limited to, related patents and trademarks) of which you are personally aware and which are associated with any part of Your Contributions.

### 5. Support

You are not expected to provide support for Your Contributions, except to the extent You desire to provide support. You may provide support for free, for a fee, or not at all. Unless required by applicable law or agreed to in writing, You provide Your Contributions on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied, including, without limitation, any warranties or conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR A PARTICULAR PURPOSE.

### 6. Work from Others

Should You wish to submit work that is not Your original creation, You may submit it to Acellera separately from any Contribution, identifying the complete details of its source and of any license or other restriction (including, but not limited to, related patents, trademarks, and license agreements) of which you are personally aware, and conspicuously marking the work as "Submitted on behalf of a third-party: [named here]".

### 7. Changes

You agree to notify Acellera of any facts or circumstances of which you become aware that would make these representations inaccurate in any respect.

Signed: 

Date: 

Projections
===========
Simulation trajectories have very large dimensionality not only in terms of number of frames but also for the
intrinsic dimensionality of each frame made by all three-dimensional coordinates. It is therefore necessary to reduce
this space by projecting these coordinates into a simpler space. HTMD provides many projection types,
e.g. MetricCoordinate to only keep the coordinate of few atoms, or MetricDistance to keep the matrix distance between
two sets of atoms.

Contents:

.. toctree::
    :maxdepth: 1

    MetricData - Storage for projected data <htmd.metricdata>
    Metric - Helper class for combining Metrics for projection <htmd.projections.metric>
    MetricCoordinate - coordinates of an atom selection <../moleculekit/moleculekit.projections.metriccoordinate.rst>
    MetricDistance - (Self-)distance-based metrics between atoms selections <../moleculekit/moleculekit.projections.metricdistance.rst>
    MetricDihedral - Dihedral-based metrics <../moleculekit/moleculekit.projections.metricdihedral.rst>
    MetricRmsd - RMSD-based metric <../moleculekit/moleculekit.projections.metricrmsd.rst>
    MetricShell - occupancy of an atom selection (e.g. water) around another selection <../moleculekit/moleculekit.projections.metricshell.rst>
    MetricSecondaryStructure - secondary structure-based metric <../moleculekit/moleculekit.projections.metricsecondarystructure.rst>
    MetricPlumed2 - access all Plumed2 metrics (CVs) <../moleculekit/moleculekit.projections.metricplumed2.rst>
    MetricSasa - Solvent accessible surface area <../moleculekit/moleculekit.projections.metricsasa.rst>
    MetricTMscore - TMscore-based metric <../moleculekit/moleculekit.projections.metrictmscore.rst>
    MetricFluctuation - RMSF-based metric <../moleculekit/moleculekit.projections.metricfluctuation.rst>Molecule
========

The Molecule class is a central object in HTMD. Most of HTMD functionalities are implemented via Molecules.
A Molecule can actually contain a molecular System (e.g. read from a PDB file), which can of course be composed of
several independent "chemical" molecules (e.g. protein + solvent + ions). This nomenclature is drawn from VMD molecules
and, in order to define "chemical" molecules, one uses segment IDs.

Molecule can read many input formats (PDB, PSF, PRMTOP, etc.) and some trajectory files (XTC, etc). Molecules can be
viewed (VMD or WebGL), aligned, selected, rotated, truncated, appended, and so on.

A very important feature is atom selection. This is identical to `VMD's Atom Selection Language`_, so that it is
possible to verify an atom selection visually and then apply it programmatically.

.. _VMD's Atom Selection Language: http://www.ks.uiuc.edu/Research/vmd/current/ug/node89.html

Contents:

.. toctree::
    :maxdepth: 1

    Molecule module <moleculekit.molecule>
Building
========

Building is the most important module for preparing Molecules for simulation. HTMD includes tools to solvate molecular
systems in water, prepare proteins for simulation by assigning appropriate protonation states, and can build all the
necessary files to simulate the system in either two of the most widely used force-fields, CHARMM and AMBER.

The workflow of the building process is normally:

#. Obtain structures
#. Clean structures
#. Define segments
#. Combine structures
#. Solvate
#. Choose force-field
#. Build and ionize

Contents:

.. toctree::
    :maxdepth: 1

    Solvating <htmd.builder.solvate>
    Protein preparation <htmd.builder.preparation>
    CHARMM builder <htmd.builder.charmm>
    AMBER builder <htmd.builder.amber>
HTMD
====

**What is it?**

HTMD is a molecular-specific programmable environment to prepare, handle, simulate, visualize and analyze molecular systems.
HTMD is based on Python, so scientists can easily extend it to their needs. With HTMD is possible to do very
complex protocols in just few lines.

In a single script, it is possible to plan an entire computational experiment, from manipulating PDBs, building,
executing and analyzing simulations, computing Markov state models, kinetic rates, affinities and pathways.

**Citing HTMD:**

If you are using HTMD in your publications please cite following papers:

| **HTMD: High-Throughput Molecular Dynamics for Molecular Discovery**
| S. Doerr, M. J. Harvey, Frank Noé, and G. De Fabritiis
| Journal of Chemical Theory and Computation 2016 12 (4), 1845-1852
| DOI: `10.1021/acs.jctc.6b00049 <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049>`_.
 


**Where to start**

The user guide is a good place to start playing around.

.. toctree::
    :maxdepth: 3

    Installation <https://software.acellera.com/install-htmd.html>
    tutorials
    documentation
    API <htmd.rst>
    developers



**Indices and tables**

* :ref:`genindex`
* :ref:`modindex`


MD simulations
==============

HTMD allows to prepare molecular simulations and run them with minimal knowledge of the technical details 
required to do so, thanks to prepared protocols. 

All functionalities of HTMD are general, it is possible to use HTMD with any MD engine. However, ACEMD the MD engine
embedded in HTMD is easier to use and is supported.   In particular, there are functions to setup the configuration of
a simulation into a specific directory and way to perform adaptive sampling methods which is the primary method of
sampling in HTMD.

ACEMD, a powerful and simple MD engine which has pioneered GPU computing since 2009, is distributed together with HTMD
in a standard version. A professional version is available by writing at `info@acellera.com`_ (for more information
visit `www.acellera.com/acemd`_).

.. _info@acellera.com: info@acellera.com
.. _www.acellera.com/acemd: http://www.acellera.com/acemd

References

*  S. Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular
   simulations, J. Chem. Theory Comput., 2014, 10 (5), pp 2064–2069. doi: `10.1021/ct400919u`_

*  M. J. Harvey, G. Giupponi and G. De Fabritiis, ACEMD: Accelerated molecular dynamics simulations in the microseconds
   timescale, J. Chem. Theory and Comput., 2009, 5 (6), 1632. doi: `10.1021/ct9000685`_

*  M. J. Harvey and G. De Fabritiis, AceCloud: Molecular Dynamics Simulations in the Cloud, J. Chem. Inf. Model.,
   2015, 55 (5), pp 909–914. doi: `10.1021/acs.jcim.5b00086`_

.. _10.1021/ct400919u: http://dx.doi.org/10.1021/ct400919u
.. _10.1021/ct9000685: http://dx.doi.org/10.1021/ct9000685
.. _10.1021/acs.jcim.5b00086: http://dx.doi.org/10.1021/acs.jcim.5b00086

Contents:

.. toctree::
    :maxdepth: 1

    MD Engines <htmd.mdengines>
    Queues <../jobqueues/jobqueues.rst>
    Adaptive sampling <adaptive>
    State-of-the-art protocols for molecular simulations <htmd.protocols>
Adaptive sampling
=================

HTMD is build around adaptive sampling, i.e. use on-the-fly information from the current data to decide where to restart new simulations. This in practice works really well with speed-up of over an order of magnitude compared to standard sampling methods. Adaptive sampling replaces to a large extent any need for biased sampling in a much safer way than biasing. Also it does not require to have very good reaction coordinates between it is building its own reaction space on-the-fly while sampling. It also integrate very well with Markov state models.

Relevant papers to read:

*  S.Doerr and G. De Fabritiis, On-the-fly learning and sampling of ligand binding by high-throughput molecular simulations, J. Chem. Theory Comput. 10 (5), pp 2064–2069(2014).

*  S.Doerr , M.J. Harvey, F. Noé ,G. De Fabritiis, HTMD: High-throughput molecular dynamics for molecular discovery, J. Chem. Theory Comput., 2016, 12 (4), pp 1845–1852


Contents:

.. toctree::
   :maxdepth: 2

   Adaptive sampling <htmd.adaptive.adaptiverun>
Tutorials
=========

The tutorials are divided in a set of sections, each containing several guides into several
functionalities of HTMD:

* The first section focuses on teaching the basics of HTMD, ideal for those starting with the software
* The second section introduces the HTMD tools available to prepare and build systems for MD simulations
* The third section explains how to setup and run MD simulations using HTMD
* The fourth section is about analysis of MD simulations using Markov state models
* The fifth section showcases some advanced uses of HTMD
* The sixth section contains examples on how to use useful tools included in HTMD
* The last section is a changelog of the modifications to HTMD along versions

Contents:

.. toctree::
    :maxdepth: 2

    Introduction to HTMD <userguide/introduction>
    Building with HTMD <userguide/building>
    Simulations in HTMD <userguide/running>
    Analysis in HTMD <userguide/analysing>
    Advanced uses <userguide/applications>
    Tools <userguide/tools>
    Changelog <userguide/changelog>

htmd
====

.. toctree::
   :maxdepth: 4

   htmd
Developer's Guide
=================

These are several tips and guidelines for developers of HTMD.


Contents:

.. toctree::
    :maxdepth: 2

    How to contribute <developers/howto>
    Guidelines for developers <developers/guidelines>
    Versioning <developers/versioning>Dimensionality Reduction
========================

On top of projection methods it is also highly recommended to use dimensionality reduction methods to further reduce
the space. HTMD provides, for instance, time independent component analysis (TICA) and Kmeans with triangle inequality.
These can be used on MetricData objects. TICA is recommended for Markov Model construction.

Contents:

.. toctree::
    :maxdepth: 1

    TICA - Time independent component analysis <htmd.projections.tica>
    KMeansTri - kmeans triangle inequality <htmd.projections.kmeanstri>
    GWPCA Principal Component Analysis <htmd.projections.gwpca>Clustering
==========

Clustering is done using the `scikit-learn clustering library`_. Other clustering classes can be used as long as they
adhere to the same interface (Methods: fit; Attributes: cluster_centers\_, labels\_).

For example, `MiniBatchKMeans`_ can be directly passed to the cluster command of MetricData::

    metricdata.cluster(MiniBatchKMeans(n_clusters=1000), mergesmall=3)


.. _scikit-learn clustering library: http://scikit-learn.org/stable/modules/clustering.html#clustering
.. _MiniBatchKMeans: http://scikit-learn.org/stable/modules/generated/sklearn.cluster.MiniBatchKMeans.html

Contents:

.. toctree::
    :maxdepth: 1

    KCenters clustering method <htmd.clustering.kcenters>
    RegCluster regular sized clustering  <htmd.clustering.regular>
Documentation
=============

You can find here documentation for most of the modules that HTMD provides. One can either explore the tree of contents
or, if looking for a particular module/class/function, use the search box on the sidebar.

Contents:

.. toctree::
    :maxdepth: 2

    Molecule <../moleculekit/moleculekit.molecule.rst>
    Building <building>
    MD Simulations <simulation>
    Simulation List <htmd.simlist>
    Projections <projections>
    Dimensionality Reduction <dimreduction>
    Clustering <clustering>
    Markov state models <htmd.model>
    Kinetics <htmd.kinetics>
    Molecular viewer <htmd.vmdviewer>


#########
Changelog
#########

1.12.0 (stable) - 2018/03/22
============================

**Main new features (relative to 1.10.0):**

- Improved User Guide: tutorials updated, new refactored projections tutorial available
- Compability with NGLview 1.0 visualization
- Support for multiple queues in Slurm and LSF
- Improvements in MOL2 reading and writing
- Introduced Membrane Builder
- Parameterization tooling gradual improvements

1.10.0 (stable) - 2017/10/25
============================

**Main new features (relative to 1.8.0):**

- Parameterize reloaded
   - More robust Psi4 settings
   - New dihedral angle fitting procedure
   - Improvement on ESP partial charges calculation
   - Refactoring of the back-end
   - Parallelized local runs and improved configurations to run on clusters
- New cofactors and non-canonical residues added for Amber
- Improvements on MOL2 file reading and writing
- Update on AceCloudQueue to launch ACEMD jobs on AceCloud from HTMD

1.8.0 (stable) - 2017/07/17
===========================

Concomitant release of 1.9.0 (latest).

**Main new features (relative to 1.6.0):**

- Adds new metric for spherical coordinates (``MetricSphericalCoordinate``).
- ``Model.markovModel``: automatic reduction of macrostates number
- Improves stability of builders, membrane tiling and equilibration protocol.
- Faster workflow from build to ``proteinprepare``.
- Psi4 version updated to 1.1
- Improves treatment of atom alternative locations
- Adds ``MetricData.sampleRegion``, which can return conformations from a specific region of data-space
- Improves file reading (CHARMM CRD CARD file format , ENT PDB files and gunzipped files).

1.6.0 (stable) - 2017/03/01
===========================

Concomitant release of 1.7.0 (latest).

**Main new features (relative to 1.4.0):**

- AdaptiveGoal implemented
- Parameterize tool released
- Apps changed for Queues
- New Metrics implemented: MetricSasa, MetricTMscore, MetricFluctuation
- Equilibration and Production protocols updated
- Documentation improved

Adaptive sampling
=================

Concept
-------

-  Exploration of the conformational space using MD simulations can
   waste lots of simulation time sampling the same conformational
   regions which does not provide any new knowledge.

-  Instead it would be desirable to explore more under-sampled regions
   of the conformational space, to overcome energetic barriers and
   eventually reach the desired conformation (folded protein / bound
   ligand etc.)

-  In adaptive sampling, instead of launching thousands of simulations
   at once from a small set of initial structures as in naive
   high-throughput MD, simulations are launched in sequential batches
   called epochs utilizing knowledge of the conformational space
   obtained from all previous epochs.

-  The starting points of the simulations in each epoch are chosen based
   on some criteria; in this case, it selects conformations from the
   most under-sampled conformational regions detected in all previous
   simulations. This is done by using Markov state models which
   discretize the conformational space into a set of most metastable
   states. Then the starting conformations are sampled based on a
   distribution related to the population of each state.

S. Doerr and G. De Fabritiis, `On-the-fly learning and sampling of
ligand binding by high-throughput molecular
simulations <http://pubs.acs.org/doi/abs/10.1021/ct400919u>`__, J. Chem.
Theory Comput. 10 (5), pp 2064–2069(2014).

S. Doerr, M. J. Harvey, Frank Noé, and G. De Fabritiis, `HTMD:
High-Throughput Molecular Dynamics for Molecular
Discovery <http://pubs.acs.org/doi/abs/10.1021/acs.jctc.6b00049>`__ J.
Chem. Theory Comput. 2016 12 (4), 1845-1852

Unit of execution
-----------------

Each simulation in adaptive is associated to a single directory which
contains all files to run it. To run a project it is therefore necessary
to provide one or more initial simulation directories, called
generators.

How to start
------------

Adaptive creates multiple directories which it uses to organize the
simulations. The user only needs to provide a generators folder
containing one sub-folder for each starting conformation containing all
files needed to execute that simulation. For example:

::

    └── generators/
        ├── gen1/
        │   ├── structure.pdb
        │   ├── input
        │   └── ...
        ├── gen2/
        │   ├── structure.pdb
        │   ├── input
        │   └── ...

Then the adaptive will generate an ``input``, ``data`` and later a
``filtered`` folder as well, looking like this:

::

    ├── data/          # Contains the completed simulations (automatically created)
    ├── filtered/      # Contains the completed simulations without water atoms (automatically created)
    ├── generators/    # Contains the initial generators provided by the user
    └── input/         # Contains the files needed to start all simulations of all epochs (automatically created)

Adaptive uses a naming scheme for simulations which follows the pattern:
``e4s3_e2s1p0f45``. This name tells us that this simulation was
generated in epoch 4 as the 3rd simulation of the batch. The starting
conformation was taken from simulation 1 of epoch 2 from the first piece
of the simulation [*]_ and from frame 45 of that simulation piece.

.. [*] some MD software might fragment simulations into pieces. Usually
       though this number will be 0 and can be ignored.

Simulation length
-----------------

1. The length of each simulation is really system dependent.
2. It could be anything like tens of nanoseconds to hundred of
   nanoseconds.
3. As a rule of thumb use twice the expected lag time for your molecular
   process (e.g. for binding anything between 30 and 100 ns).

Simulation details
------------------

As only the coordinates files are seeded for every new epoch,
simulations cannot use a velocity file. Velocities are therefore
reinitialized to the Maxwell Boltzmann distribution at the given
temperature.

E.g. if setting up the simulations with the :class:`~htmd.protocols.production_v6.Production` class:

.. code:: python

    from htmd.protocols.production_v6 import Production
    md = Production()
    md.adaptive = True
    [...]

or directly modifying the ACEMD ``input`` file of the simulations and
removing the binvelocities line.

Adaptive script example
-----------------------

The power of adaptive sampling is accessible on HTMD through the :class:`~htmd.adaptive.adaptiverun.AdaptiveMD` class:

.. code:: python

    from htmd.ui import *
    app = LocalGPUQueue()
    app.datadir = './data'
    md = AdaptiveMD()
    md.nmin=5
    md.nmax=10
    md.nepochs = 30
    md.app = app
    md.projection = MetricDistance('name CA', '(resname BEN) and ((name C7) or (name C6))', periodic='selections', metric='contacts')
    md.ticadim = 3
    md.updateperiod = 14400 # execute every 4 hours
    md.run()

Execution in a notebook
-----------------------

1. It is possible to run the adaptive scheme syncronosly or
   asyncrounsly.
2. The option ``updateperiod`` controls this behaviour.
3. The default is to run and exit, so ``updateperiod`` needs to be specified
   if adaptive should be run synchronously

Setting a simple cron job
-------------------------

1. This is useful for having the script execute automatically every x
   hours.
2. Do not set ``updateperiod`` then, or set it to zero such that the
   program will execute and exit

.. code:: bash

    #!/bin/bash -login
    # cron.sh file
    # use crontab -e to add this line:
    # 0 */4 * * * cd /pathtomydir/; ./cron.sh
    #
    python conf.py

Visualizing the starting conformations
--------------------------------------

If we want to look at what structures were chosen as the starting
conformations of a given epoch we can use a code snippet like the
following:

.. code:: python

    for s in glob('input/e28s*'):  # Visualize all starting conf of epoch 28
       mol = Molecule(s+'/structure.pdb')
       mol.read(s+'/input.coor')
       mol.view()
Advanced System Building in HTMD
================================

PDB format considerations
-------------------------

The PDB format is very old. In an effort to handle its legacy shortcomings, several versions have been made over the
years, they are not all readily interchangeable, and not all software can handle each version perfectly. The most
important things to watch out for are:

- Columns: the PDB format has very rigid rules about what values can go in each space. Keep in mind that it is not a
  space/tab/comma delimited format, but rather has rigid definitions of what should be in each space/column.

- The PDB format as originally designed cannot handle more than 9,999 resids or 99,999 atoms (due to the column format
  issue). Several workarounds have been devised, such as using hexadecimal numbers or other compact number formats. VMD
  has no trouble saving more atoms/residues.

Physical and chemical considerations
------------------------------------

One needs to know well the working system, thus:

- Always review your PDB file: inspect the REMARK sections of the PDB file. You can often find key specific information
  regarding the structure (e.g. disulphide bonds, missing atoms, etc.).

- Disulphide bonds present in the system must be identified. This is automatically done by HTMD in both CHARMM and Amber.

- Metalloproteins: if the metal ion is not an active part of an interaction it may be acceptable to just allow it to act
  as a cation perhaps restraining it with some harmonic constraints if necessary.

- Duplicate atoms in the PDB file: typically simply delete one of the duplicated groups. However, if both conformations
  are potentially important (e.g. such loops involved in molecular recognition) it might be necessary to simulate both
  conformations separately.

Protonation/pH
--------------

The protonation state of the system is critical. Since MD simulations typically don't allow for bond breaking, the
initial protonation of the system must be accurate. Knowing what pH you are trying to reproduce is therefore important
to obtain the correct results. If you suspect changing protonation is important to your system and you still want to use
classical mechanics, consider simulating both states (protonated and not protonated).

Histidine residues can have three different protonations states even at pH 7, therefore, a correct protonation of this
residue is particularly critical. This residue can be protonated at either delta (most common; HSD/HID), epsilon (very
common also; HSE/HIE) or at both nitrogens (special situations and low pH; HSP/HIP).

.. image:: http://pub.htmd.org/tutorials/advancedbuilding/histidines.png

The best way to determine how histidine should be protonated is to look at the the structure. Typically, a histidine
residue is protonated if it is close enough to an electron donor (e.g. a glutamic acid), thus creating a hydrogen bond.
Since histidines are frequently present at protein active sites, a correct protonation state is particularly important
in ligand binding simulations.

In HTMD, one can use :class:`~moleculekit.tools.preparation.systemPrepare` to help with protonation.

List of useful tools
--------------------

====================================================== ==================
:class:`~moleculekit.molecule.Molecule` class        Building functions
====================================================== ==================
:func:`~moleculekit.molecule.Molecule.append`        :func:`~moleculekit.util.maxDistance`
:func:`~moleculekit.molecule.Molecule.center`        :func:`~moleculekit.util.uniformRandomRotation`
:func:`~moleculekit.molecule.Molecule.mutateResidue` :func:`~moleculekit.util.boundingBox`
:func:`~moleculekit.molecule.Molecule.moveBy`        :func:`~moleculekit.util.sequenceID`
:func:`~moleculekit.molecule.Molecule.rotateBy`      :func:`~htmd.builder.builder.embed`
---                                                    :func:`~moleculekit.tools.autosegment.autoSegment`
====================================================== ==================

List of common CHARMM patches
-----------------------------

- C-terminal patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
CTER -1     standard C-terminus
CT1  0      methylated C-terminus from methyl acetate
CT2  0      amidated C-terminus
CT3  0      N-Methylamide C-terminus
==== ====== ===========

- N-terminal patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
NTER +1     standard N-terminus
ACE  0      acetylated N-terminus (to create dipeptide)
ACP  0      acetylated N-terminus (for proline dipeptide)
PROP +1     Proline N-Terminal
GLYP +1     Glycine N-terminus
==== ====== ===========

- Side-chain patches

==== ====== ===========
Name Charge Description
==== ====== ===========
ASPP 0      patch for protonated aspartic acid, proton on OD2
GLUP 0      patch for protonated glutamic acid, proton on OE2
CYSD -1     patch for deprotonated CYS
DISU +1     patch for disulfides. Patch must be 1-CYS and 2-CYS
HS2  +1     patch for neutral His, move proton from ND1 to NE2
TP1  -1     patch to convert tyrosine to monoanionic phosphotyrosine
TP1A -1     patch to convert tyrosine to monoanionic phenol-phosphate model compound when generating tyr, use first none last none for terminal patches
TP2  -2     patch to convert tyrosine to dianionic phosphotyrosine
TP2A -2     patch to convert tyrosine to dianionic phosphotyrosine when generating tyr, use first none last none for terminal patches this converts a single tyrosine to a phenol phosphate
TMP1 -1     patch to convert tyrosine to monoanionic phosphonate ester O -> methylene (see RESI BMPH)
TMP2 -2     patch to convert tyrosine to dianionic phosphonate ester O -> methylene (see RESI BMPD)
TDF1 -1     patch to convert tyrosine to monoanionic difluoro phosphonate ester O -> methylene (see RESI BDFH)
==== ====== ===========

- Circular protein chain patches:

==== ====== ===========
Name Charge Description
==== ====== ===========
LIG1 0      linkage for cyclic peptide, 1 refers to the C terminus which is a glycine , 2 refers to the N terminus
LIG2 0      linkage for cyclic peptide, 1 refers to the C terminus, 2 refers to the N terminus which is a glycine
LIG3 0      linkage for cyclic peptide, 1 refers to the C terminus which is a glycine, 2 refers to the N terminus which is a glycine
==== ====== ===========
Analysis in HTMD
================

Contents:

.. toctree::
    :maxdepth: 1

    Getting Started with Projections <../tutorials/projections>
    Ligand binding analysis <../tutorials/ligand-binding-analysis>
    Protein folding analysis <../tutorials/protein-folding-analysis>
    CXCL12 conformational analysis <../tutorials/cxcl12-conformational-analysis>Building with HTMD
==================

Contents:

.. toctree::
    :maxdepth: 1

    Protein preparation <../tutorials/protein-preparation>
    Advanced building <advancedbuilding>
    System building protein-ligand <../tutorials/system-building-protein-ligand>
    System building protein in Membrane <../tutorials/system-building-protein-in-membrane>Handy tools in HTMD
===================

Contents:

.. toctree::
    :maxdepth: 1

    Maximal Substructural Alignment <../tutorials/MaxSubstructuralAlignment>
    Sequence Based Alignment <../tutorials/SequenceAlignment>
    FFEvaluate <../tutorials/FFEvaluate>
    Membrane Builder <../tutorials/MembraneBuilder>
Introduction to HTMD
====================

Contents:

.. toctree::
    :maxdepth: 1

    Python primer <../tutorials/python-primer>
    Getting started with HTMD <../tutorials/getting-started-with-htmd>
    HTMD Molecules <../tutorials/htmd-molecules>
    Visualization in HTMD <../tutorials/visualization>Simulations in HTMD
===================

Contents:

.. toctree::
    :maxdepth: 1

    Molecular Dynamics in HTMD <../tutorials/md-protocols>
    Wrapping Simulations <../tutorials/WrappingTutorial>
    Adaptive Sampling Explained <adaptive-sampling-explained>
    Adaptive Sampling Tutorial <../tutorials/adaptive-sampling>
    Adaptive Bandit Tutorial <../tutorials/adaptive-bandit>
Advanced Uses of HTMD
=====================

Contents:

.. toctree::
    :maxdepth: 1

    System building Protein-Protein <../tutorials/system-building-protein-protein>
    Using docking to initialize positions <../tutorials/docking-simulation-generators>
    Equilibration protocol for the mu opioid GPCR <../tutorials/mu-opioid-receptor-gpcr-equilibration>##############################
Guidelines for HTMD developers
##############################

Coding Style Tips
=================

Note: some style is at variance with Python PIP recommendations.
 
* Class names start with capital letters
* Method names start with lower-case, then ``camelCase``
* Same for generic function, e.g. ``testMe()``
* Modules should be nouns
* Methods and functions should be verbs
* Make the main a test case, where possible.
* Use namespaces instead of composite name, e.g. ``charmm.build()`` instead of ``charmmBuild()`` when possible
* Try to keep single names when possible, so there is no need for camelCase


Using modules inside HTMD
=========================

Creating efficient code is different from generating scripts for personal use. As such, when developing code for HTMD,
be as minimal as possible when importing modules. As examples, please do _not_ use any of these types of imports inside
HTMD code:

.. code:: python

    import htmd
    from htmd.ui import *

Furthermore, from ``from <module> import *`` should _never_ be used, as it pollutes the namespace and can shadow same-name
functionalities. Try as much as possible to only import the function/class you specifically need instead of importing an
entire module, unless one wants to use that module heavily on that implementation. Keep a simple
``htmd.modulename.submodulename`` structure. So file names if they are not meant to be modules (e.g. util.py) should be
imported in the upper module namespace.

Do not pollute the module and submodule names. Changes to the modules structure requires consensus and approval, as well
as documentation creation.


Using Docstrings and Doctests
=============================

See this template:

.. code:: python

    def home(dataDir=None, libDir=False):
        """Return the pathname of the HTMD root directory (or a data subdirectory).

        Parameters
        ----------
        dataDir : str
            If not None, return the path to a specific data directory
        libDir : bool
            If True, return path to the lib directory

        Returns
        -------
        dir : str
            The directory

        Example
        -------
        >>> htmd.home()                                 # doctest: +ELLIPSIS
        '.../htmd'
        >>> htmd.home(dataDir="dhfr")                   # doctest: +ELLIPSIS
        '.../data/dhfr'
        >>> os.path.join(htmd.home(dataDir="dhfr"),"dhfr.pdb")  # doctest: +ELLIPSIS
        '.../data/dhfr/dhfr.pdb'
        """

Docstrings can be test cases (as above). This is convenient because you have four things in one place:

#. the test case
#. the expected result
#. an example
#. the rest of the documentation

It's sufficient to add this in the main:

.. code:: python

    if __name__ == "__main__":
        import doctest

        failure_count, _ = doctest.testmod()
        if failure_count != 0:
            raise Exception('Doctests failed')


The ``doctest: +ELLIPSIS`` comment on the docstring indicates that match with ``...`` is flexible.
Other possibly useful directives are ``SKIP`` and ``NORMALIZE_WHITESPACE``.

One can also:

- run tests placed in external files with ``doctest.testfile('doctest_in_help.rst')``
- test a different module with ``doctest.testmod(doctest_simple)``
#########################
How to Contribute to HTMD
#########################

Introduction
============

Acellera Ltd. (`<https://github.com/Acellera>`_) is the owner of the HTMD repository (`<https://github.com/Acellera/htmd>`_).
The HTMD repository is composed of a ``master`` branch, where all development is carried on.

How to Contribute Fixes
=======================

If you have found a bug in HTMD, the correct way to fix it is to fork the Acellera/htmd repository
(`<https://github.com/Acellera/htmd/fork>`_), switch to the branch corresponding to the latest stable
release, fix the bug, test the fix, and open a pull request (`<https://github.com/Acellera/htmd/compare>`_).

How to Contribute New Features
==============================

If you'd like to develop for HTMD and contribute new features, it would be best to contact us through the issues
(`<https://github.com/Acellera/htmd/issues>`_), so we can better advice you on your contribution. The process of
contributing is also to fork the Acellera/htmd repository (`<https://github.com/Acellera/htmd/fork>`_),
develop your feature (rebase often), document it, test it, and when it's ready, open a pull request to the master branch
(`<https://github.com/Acellera/htmd/compare>`_).

Contributor Agreement
=====================

To find the contributor agreement for signing, please follow this link
(`<https://github.com/Acellera/htmd/blob/master/doc/source/developers/contributor_agreement.md>`_). Please add the
signed document to your PR (inside ``doc/source/developers/signed_agreements/``). We can only accept your PR after you
have added this file.

List of Features Interesting for HTMD
=====================================

TBD##########
Versioning
##########

HTMD Versions Explained
=======================

HTMD versions are of the following format:

```
<big release>.<major release (stable and develment tagging)>.<minor release (bug fixes)>
```

Big release
-----------

Only changed at the developers' discretion. This may only happen when there are ground-breaking changes to the software.

Major release - Definition of _stable_ and _development_ versions
-----------------------------------------------------------------

Major releases are periodic (every 3 or 6 months, check the [milestones](https://github.com/Acellera/htmd/milestones)
for when the next one comes out). When a release is done, actually two versions come out:
a __stable__ one and a __development__ one.

When the second/middle number (major release) is even (0,2,...,2n), it means that this major release is a __stable__
release. Examples: 1.0.0; 1.2.0; 1.24.0.
In this context, the third/right number will be 0 (zero) when the major release comes out and will increase until a
new major release is done.

For each stable release, there is a corresponding __development__ release, whose second/middle number is odd (1,3,...,2n+1) and always one more than the stable. Examples: 1.1.0; 1.3.0; 1.25.0.
The third/right number is also a 0 (zero) when the major release comes out and will increase until a new major release is done.

Minor release
-------------

Minor releases correspond to the third/right number and are, in general, periodic (every two weeks). In the case of being a _stable_ version, they correspond to critical bug-fixes and may be released as soon as the fix is made. 

This number is always 0 (zero) at major releases and the corresponding _stable_ and _development_ versions (example: 1.0.0 and 1.1.0 correspond) are "in sync" when this happens.

When this number is different from 0 (zero), nothing should be assumed in the relationship between corresponding _stable_ and _development_ versions.

GitHub Tagging
--------------

Stable releases are tagged on a specific branch created for each of these (named ``rel-<big release>.<major release>.x``). Development releases are tagged on the master branch.

Useful tag listing commands:
```
git tag -n
git describe --tags
```

How to release a new HTMD version?
==================================

These are examples of how to release HTMD versions.

Big and/or major releases
-------------------------

Imagine one wants to do a big/major release (in this case, let's assume it's the first major release of big release 1).

1. Make sure you are working on `Acellera/htmd:master` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git checkout master
   ```

1. Make sure the ``master`` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. On ``master``, create the new stable branch, check it out, and tag it:

   ```
   git branch rel-1.0.x
   git checkout rel-1.0.x
   git tag -a 1.0.0 -m "new stable release"
   ```

1. Push the new branch and tag to the remote (``origin``):

   `git push --tags origin rel-1.0.x`

   This will trigger two Travis builds: one due to the branch and another due to the tag. A conda release will be made.

1. Check out ``master``, tag the development version, and push the tag:

   ```
   git checkout master
   git tag -a 1.1.0 -m "new development release"
   git push --tags
   ```

   This will trigger another Travis build. A conda release will be made.

These two tags will point to the same commit (in sync).

Minor releases - stable
-----------------------

Imagine one wants to do a minor release (bug-fix) on release 1.0.0.

1. Make sure you are working on `Acellera/htmd:rel-1.0.x` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git fetch
   git checkout rel-1.0.x
   ```

1. Make sure the ``rel-1.0.x`` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. Do the fix, add the files, and commit it.
1. Tag the new bug-fix:

   `git tag -a 1.0.1 -m "new bug-fix"`

1. Push the fix and the tag to the remote (``origin``):

   `git push --tags origin rel-1.0.x`

   This will trigger a new Travis build. A conda release will be made.

In alternative, push the commit but do not tag it, if it is not critical to release right now.

Minor releases - development
----------------------------

Imagine one wants to do a minor release on release 1.1.0.

1. Make sure you are working on `Acellera/htmd:master` (https://github.com/Acellera/htmd.git):

   ```
   git remote -v
   git checkout master
   ```

1. Make sure the ``master`` branch is up-to-date:

   ```
   git fetch
   git pull
   ```

1. Tag the new minor release:

   `git tag -a 1.1.1 -m "new minor release"`

1. Push the tag to the remote (``origin``):

   `git push --tags origin master`

   This will trigger a new Travis build. A conda release will be made.
