# Missing wedge restoration (MWR)

The code in this repository is described in [this pre-print](https://hal.inria.fr/hal-01966821/document). This paper has been submitted to Journal of Structural Biology and is currently under revision.

## Contents
- [System requirements](##System requirements)
- [Installation guide](##Installation guide)
- [Instructions for use](##Instructions for use)

## System requirements
__MWR__ has been implemented using Matlab R2015a. It has been tested on Linux (Debian 8.6) and Mac OSX (10.12.3) and should also work on Windows.

### Dependencies
Our code needs following conditions to run:
- Matlab Image Processing Toolbox to display images (the processing itself does not depend on this toolbox)
- BM4D denoising algorithm implemented by M. Maggioni

## Installation guide
Download the BM4D implementation [here](http://www.cs.tut.fi/~foi/GCF-BM3D/BM4D_v3p2.zip), and unzip the folder in utils/. That's it.

## Instructions for use
Instructions for using MWR are contained in the Matlab scripts example_*.m. 
- example_proteasome.m applies the method on synthetic data.
- example_gold_particle.m applies the method on experimental data.=====================
ContinuousFlex plugin
=====================

This plugin provides the latest Scipion protocols for cryo-EM continuous conformational flexibility/heterogeneity analysis of biomolecular complexes.


Installation
------------

You will need to use `3.0 <https://github.com/I2PC/scipion/releases>`_ version of Scipion to be able to run these protocols. To install the plugin, you have two options:
We you need help installing Scipion3, please refer to the Scipion Documentation `here <https://scipion-em.github.io/docs/docs/scipion-modes/how-to-install.html>`__

Make sure that you have cmake installed on your Linux system. For example, if you are using Ubuntu
 .. code-block::

    sudo apt install cmake


a) Stable version

	Install Scipion3 and use the plugin manager to install the plugin.

b) Developer's version

   * download repository

   .. code-block::

      git clone https://github.com/scipion-em/scipion-em-continuousflex.git
      git checkout devel

   * install

   .. code-block::

      scipion3 installp -p path_to_scipion-em-continuousflex --devel

continuousflex sources will be downloaded automatically with the plugin.


Note: Xmipp and Chimerax plugins should be installed (from Scipion3 plugin manager) to run continuousflex protocols.
You should also consider having VMD on your system for visualization.
We assume that VMD is installed on your system in "/usr/local/lib/vmd".
If VMD is installed but does not work, you may run the command "scipion3 config" and look for VMD_HOME in the config file (the config file is usually at ~/scipion3/config/scipion.conf)

Supported versions
------------------

versions > 3.0.15

Protocols
---------

* HEMNMA: Hybrid Electron Microscopy Normal Mode Analysis method to interpret heterogeneity of a set of single particle cryo-EM images in terms of continuous macromolecular conformational transitions [1-3]
* StructMap: Structural Mapping method to interpret heterogeneity of a set of single particle cryo-EM maps in terms of continuous conformational transitions [4]
* HEMNMA-3D: Extension of HEMNMA to continuous conformational variability analysis of macromolecules from in situ cryo-ET subtomograms [5]
* TomoFlow: Method for analyzing continuous conformational variability of macromolecules in in vitro and in situ cryogenic subtomograms based on 3D dense optical flow [7]

Notes:

* The plugin additionally provides the test data and automated tests of the protocols in Scipion 3. The following two types of tests of HEMNMA and HEMNMA-3D can be produced by running, in the terminal, "scipion3 tests continuousflex.tests.test_workflow_HEMNMA" and “scipion3 tests continuousflex.tests.test_workflow_HEMNMA3D”, respectively: (1) tests of the entire protocol with the flexible references coming from an atomic structure and from an EM map; and (2) test of the alignment module (test run using 5 MPI threads). The automated tests of the TomoFlow method are also available and can be run using scipion3 tests continuousflex.tests.test_workflow_TomoFlow.
* HEMNMA additionally provides tools for synthesizing noisy and CTF-affected single particle cryo-EM images with flexible or rigid biomolecular conformations, for several types of conformational distributions, from a given atomic structure or an EM map. One part of the noise is applied on the ideal projections before and the other after the CTF, as described in [6].
* HEMNMA-3D additionally provides tools for synthesizing noisy, CTF and missing wedge affected cryo-ET tomograms and single particle subtomograms with flexible or rigid biomolecular conformations, for several types of conformational distributions, from a given atomic structure or an EM map. One part of the noise is applied on the ideal projections before and the other after the CTF, as described in [6].
* A reproduction of some utility codes with their corresponding licenses are contained in this plugin for subtomogram averaging, missing wedge correction, denoising and data reading. These codes are not used in the methods above, but they are made optional for data preprocessing and visualization.


References
----------
[1] Jin Q, Sorzano CO, de la Rosa-Trevin JM, Bilbao-Castro JR, Nunez-Ramirez R, Llorca O, Tama F, Jonic S: Iterative elastic 3D-to-2D alignment method using normal modes for studying structural dynamics of large macromolecular complexes. Structure 2014, 22:496-506. `[Open-access] <http://www-ext.impmc.upmc.fr/~jonic/Papers/HEMNMA.pdf>`__

[2] Jonic S: Computational methods for analyzing conformational variability of macromolecular complexes from cryo-electron microscopy images. Curr Opin Struct Biol 2017, 43:114-121. `[Link] <http://dx.doi.org/10.1016/j.sbi.2016.12.011>`__ `[Author’s version] <http://www-ext.impmc.upmc.fr/~jonic/Papers/CurrentOpinionStructBiol_Jonic_2017.pdf>`__

[3] Harastani M, Sorzano CO, Jonic S: Hybrid Electron Microscopy Normal Mode Analysis with Scipion. Protein Sci 2020, 29:223-36. `[Open-access] <https://onlinelibrary.wiley.com/doi/epdf/10.1002/pro.3772>`__

[4] Sanchez Sorzano CO, Alvarez-Cabrera AL, Kazemi M, Carazo JM, Jonic S: StructMap: Elastic Distance Analysis of Electron Microscopy Maps for Studying Conformational Changes. Biophys J 2016, 110:1753-1765. `[Open-access] <http://www-ext.impmc.upmc.fr/~jonic/Papers/StructMap.pdf>`__

[5] Harastani M, Eltsov M, Leforestier A, Jonic S: HEMNMA-3D: Cryo Electron Tomography Method Based on Normal Mode Analysis to Study Continuous Conformational Variability of Macromolecular Complexes. Front Mol Biosci 2021, 8:663121. `[Open-access] <https://www.frontiersin.org/articles/10.3389/fmolb.2021.663121/abstract>`__

[6] Jonic S, Sorzano CO, Thevenaz P, El-Bez C, De Carlo S, Unser M: Spline-based image-to-volume registration for three-dimensional electron microscopy. Ultramicroscopy 2005, 103:303-317. `[Author’s version] <http://www-ext.impmc.upmc.fr/~jonic/Papers/Ultramicroscopy_2005_v103_p303.pdf>`__

[7] Harastani M, Eltsov M, Leforestier A, Jonic S: TomoFlow: Analysis of continuous conformational variability of macromolecules in cryogenic subtomograms based on 3D dense optical flow. J Mol Biol 2021,167381. `[Author’s version] <https://hal.archives-ouvertes.fr/hal-03452809>`__ `[Journal] <https://doi.org/10.1016/j.jmb.2021.167381>`__

# scipion-em-continuousflex
