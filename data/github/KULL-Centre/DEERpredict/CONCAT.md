[![Build Status](https://app.travis-ci.com/KULL-Centre/DEERpredict.svg?branch=main)](https://app.travis-ci.com/KULL-Centre/DEERpredict)
[![Documentation Status](https://readthedocs.org/projects/deerpredict/badge/?version=latest)](https://deerpredict.readthedocs.io)
[![DOI](https://zenodo.org/badge/217526987.svg)](https://zenodo.org/badge/latestdoi/217526987)
[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/KULL-Centre/DEERpredict/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/KULL-Centre/DEERpredict)

DEER-PREdict
===========

Overview
--------

A package for double electron-electron resonance (DEER) and paramagnetic relaxation enhancement (PRE) predictions from molecular dynamics ensembles.

Installation
------------

To install DEER-PREdict, use the [PyPI package](https://pypi.org/project/DEERPREdict):

```bash
  pip install DEERPREdict
```

or clone the repo:

```bash
  git clone https://github.com/KULL-Centre/DEERpredict.git
  cd DEERpredict

  pip install -e . 
```

The software requires Python 3.6+.
    
Documentation
-------------

[![Documentation Status](https://readthedocs.org/projects/deerpredict/badge/?version=latest&style=for-the-badge)](https://deerpredict.readthedocs.io)


Testing
-------

Run all the tests in one go

```bash
  cd DEERpredict

  python -m pytest
```
or run single tests, e.g.

```bash
  cd DEERpredict

  python -m pytest tests/test_PRE.py::test_ACBP
  python -m pytest tests/test_DEER.py::test_T4L
```


Authors
-------------

[Giulio Tesei (@gitesei)](https://github.com/gitesei)

[João M Martins (@joaommartins)](https://github.com/joaommartins)

[Micha BA Kunze (@mbakunze)](https://github.com/mbakunze)

[Ramon Crehuet (@rcrehuet)](https://github.com/rcrehuet)

[Kresten Lindorff-Larsen (@lindorff-larsen)](https://github.com/lindorff-larsen)


Article
-------------

Tesei G, Martins JM, Kunze MBA, Wang Y, Crehuet R, et al. (2021) 
DEER-PREdict: Software for efficient calculation of spin-labeling EPR and NMR data from conformational ensembles. 
PLOS Computational Biology 17(1): e1008551. [https://doi.org/10.1371/journal.pcbi.1008551](https://doi.org/10.1371/journal.pcbi.1008551)
# Libraries

New libraries should be included in the lib folder and listed in the yaml file `DEERPREdict/lib/libraries.yml`.

The currently available libraries are the 175 K and 298 K 1-Oxyl-2,2,5,5-tetramethylpyrroline-3-methyl methanethiosulfonate (MTSSL) 
libraries developed by Yevhen Polyhach, Enrica Bordignon and Gunnar Jeschke (DOI 10.1039/C0CP01865A) and implemented in MMM (DOI 10.1002/pro.3269).

| name                      | coordinates                  | reference            |
|:-------------------------:|:----------------------------:|:--------------------:|
| MTSSL 175K X1X2           | MTSSL_175K_X1X2_46.txt       |DOI 10.1039/C0CP01865A|   
| MTSSL 175K CASD           | MTSSL_175K_CaSd_216.txt      |DOI 10.1039/C0CP01865A|
| MTSSL 298K UFF r1 CASD    | R1A_298K_UFF_216_r1_CASD.txt |DOI 10.1002/pro.3269  |
| MTSSL 298K UFF r2 CASD    | R1A_298K_UFF_216_r2_CASD.txt |DOI 10.1002/pro.3269  |
| MTSSL 298K UFF r3 CASD    | R1A_298K_UFF_216_r3_CASD.txt |DOI 10.1002/pro.3269  |
| MTSSL 298K UFF r4 CASD    | R1A_298K_UFF_216_r4_CASD.txt |DOI 10.1002/pro.3269  |
| MTSSL 298K UFF r5 CASD    | R1A_298K_UFF_216_r5_CASD.txt |DOI 10.1002/pro.3269  |

Copyright © 2009-2013, Yevhen Polyhach, Stefan Stoll & Gunnar Jeschke

All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
- Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer 
  in the documentation and/or other materials provided with the distribution

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS 
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE 
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# Running Calculations

## PREpredict

The class for running the PRE calculations is `PREpredict`.

Here is how to calculate intensity ratios and PRE rates for a mutant of PDB 1NTI (20 conformations) using the rotamer-library approach

~~~ python
from DEERPREdict.PRE import PREpredict

u = MDAnalysis.Universe('1nti.pdb')
PRE = PREpredict(u, residue = 36, log_file = 'log', temperature = 298, atom_selection = 'H')
PRE.run(output_prefix = 'res', tau_c = 2*1e-09, tau_t = .5*1e-9, delay = 10e-3, r_2 = 10, wh = 750)
~~~

The program generates a data file called `res-36.dat`. The first column contains the residue numbers while the second 
and third contain the intensity ratios and the PRE rates (in Hz), respectively.
Additionally, the Boltzmann weighted averages of $r^{-3}$, $r^{-6}$ and $( 3 \cos^{2} \Omega - 1 ) / 2$ over the rotamer states are saved to the pickle file `res-36.pkl`.
The sum over the Boltzmann weights for the Lennard-Jones probe-protein interaction energies are saved to `res-Z-36.pkl`.

### Reweighting

Per-frame distances and angles are saved in the pickle file `res-36.pkl`. These intermediate quantities can be used to reweight the trajectory
by statistical weights obtaine e.g. from BME reweighting:

~~~ python
PRE.run(output_prefix = 'calcPREs/res', tau_c = 2*1e-09, tau_t = .5*1e-9, delay = 10e-3, r_2 = 10, wh = 750, weights = weights, load_file = 'res-36.pkl')
~~~

### Intermolecular PREs

To calculate intermolecular PREs, a list of two strings can be set to the `chains` option indicating the segment id of the labeled chain and the NMR-active chain.

~~~ python
u = MDAnalysis.Universe('3BVB.pdb')
PRE = PREpredict(u, residue = 55, chains = ['A', 'B'], log_file = 'calcPREs/log', temperature = 298, atom_selection = 'N')
PRE.run(output_prefix = 'calcPREs/res', tau_c = 2*1e-09, tau_t = .5*1e-9, delay = 10e-3, r_2 = 10, wh = 750)
~~~

### Approximate electron positions to C$\beta$ coordinates

Instead of using the rotamer library approach, the position of the unpaired electron can be approximated to the position of the C$\beta$ atom of the spin-labeled residue
setting `Cbeta=True`.

~~~ python
from DEERPREdict.PRE import PREpredict

u = MDAnalysis.Universe('1nti.pdb')
PRE = PREpredict(u, residue = 36, log_file = 'log', temperature = 298, atom_selection = 'H', Cbeta = True)
PRE.run(output_prefix = 'calcPREs/res', tau_c = 2*1e-09, tau_t = .5*1e-9, delay = 10e-3, r_2 = 10, wh = 750)
~~~


## DEERpredict

Here is an example of how to run DEERpredict to calculate the DEER distribution for HIV-1 protease (PDB ID 3BVB) labeled with nitroxide groups at residue 55.

~~~ python
from DEERPREdict.DEER import DEERpredict

DEER = DEERpredict(MDAnalysis.Universe('3BVB.pdb'), residues = [55, 55], chains=['A', 'B'], log_file = 'log', temperature = 298 )
DEER.run(output_prefix = 'res')
~~~

DEERpredict generates `res-55-55.dat` containing the smoothed distance distribution and `res-55-55_time_domain.dat` containing the time-domain 
DEER data (Eq. 3 in DOI: 10.1126/sciadv.aat5218).
Per-frame distance distributions are saved to the hdf5 file 'res-55-55.hdf5' making it possible to quickly reweight the data, as shown above for PREpredict.
The function to back-calculate the time-domain data from a distance distribution can also be accessed externally from the `Operations` class.

~~~ python
from DEERPREdict.utils import Operations

r, p = np.loadtxt('res-55-55.dat', unpack=True)
t = np.linspace(0.01, 5.5, 512)
dt = Operations.calcTimeDomain(t,r,p)
~~~

The sums over the Boltzmann weights for the Lennard-Jones probe-protein interaction energies of positions K55 and K55' are saved to `res-Z-55-55.pkl`.
The upper bound for the inter-probe distances can be set using the `rmax` option (default 12 nm). The interval of the time variable can be set with the options `tmin`, `tmax` and `dt`, whose default values are 0.01, 5.5 and 0.01 microseconds, respectively.
The standard deviation of the Gaussian low pass filter (default 0.05 nm), which is applied to the distance distribution for noise reduction, can be set via the option `filter_stdev` of the `run()` function, as shown below.

~~~ python
from DEERPREdict.DEER import DEERpredict

DEER = DEERpredict(MDAnalysis.Universe('3BVB.pdb'), residues = [55, 55], chains=['A', 'B'], log_file = 'log', temperature = 298, rmax = 11)
DEER.run(output_prefix = 'res', filter_stdev = 0.07)
~~~

# Installation

DEER-PREdict requires Python 3.6+. An environment with Python 3.6+ can be readily set up using conda:

~~~ bash
conda create --name DEER-PRE python>=3.6
conda activate DEER-PRE
~~~

Within an environment with Python 3.6+, DEER-PREdict can be installed using the PyPI package:

~~~ bash
pip install DEERPREdict
~~~

or from a local copy of the repository:

~~~ bash
git clone https://github.com/KULL-Centre/DEERpredict.git
cd DEERpredict
pip install -e . 
~~~
# DEER-PREdict API Reference

DEER-PREdict has three main classes:
- `DEERpredict` performs Double Electron-Electron Resonance predictions. <br>
   Functions: `trajectoryAnalysis`, `run` and `save`.
- `PREpredict` performs Paramagnetic Relaxation Enhancement calculations. <br>
   Functions: `trajectoryAnalysis`, `trajectoryAnalysisCbeta`, `run` and `save`.
- `Operations` is the base class containing attributes and methods inherited and used by the calculation classes. <br> 
   Functions: `precalculate_rotamer`, `rotamer_placement`, `lj_calculation`, `rotamerWeights`, `rotamerPREanalysis`, `calc_gamma_2`, `calc_gamma_2_Cbeta` and `calcTimeDomain`.

New rotamer libraries should be defined in the `DEERpredict.libraries.LIBRARIES` dictionary.

## RotamerLibrary class

`DEERpredict.libraries.LIBRARIES`: Loaded from `DEERPREdict/lib/libraries.yml`, rotamers libraries consist of a PDB file, a DCD files and a text file for the weights. These files are included in the `DEERPREdict/lib` folder.

`DEERpredict.libraries.RotamerLibrary`: Makes available the attributes `top` (rotamer topology), `coord` (rotamer coordinates) and `weights` (intrinsic probability of each rotamer).

## Lennard-Jones parameters

`DEERpredict.lennardjones`: Lennard-Jones parameters of the CHARMM36 force field used to calculate the external 
energy contribution to the Boltzmann weight of each conformer.

~~~python 
DEERpredict.lennardjones.lj_parameters = {
    'C': {
        'vdw': 1.70,
        'p_q': 0,
        'p_Rmin2': 2.275,
        'eps': -0.020
    }, 
    ...
}
~~~
<script type="text/x-mathjax-config">
MathJax.Hub.Config({
  tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
});
</script>
<script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.0/MathJax.js?config=TeX-AMS-MML_HTMLorMML" type="text/javascript"></script>

# Introduction

DEER-PREdict is a software for the prediction of Double Electron-Electron Resonance (DEER) distance distributions and Paramagnetic Relaxation Enhancement (PRE) rates from conformational ensembles.

The development is a team effort with contributions from João Martins, Micha BA Kunze, Ramon Crehuet and Giulio Tesei.
# Testing

Tests can be run using the test running tool pytest:

~~~ bash
git clone https://github.com/KULL-Centre/DEERpredict.git
cd DEERpredict
python -m pytest
~~~

The tests reproduce reference data for four protein systems:
- HIV-1 Protease
- T4 Lysozyme
- Acyl-CoA-Binding Protein
- A discoidal complex of Apolipoprotein A-I

Test systems can be further explored through the Jupyter Notebooks in the `tests/data` folder:
- `HIV-1PR/HIV-1PR.ipynb`
- `nanodisc/nanodisc.ipynb`
- `ACBP/ACBP.ipynb`
- `article.ipynb`
.. DEERpredict documentation master file, created by
   sphinx-quickstart on Wed Apr 22 11:08:18 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to DEER-PREdict!
=======================

.. toctree::
   :maxdepth: 2

   _docs/docs
   _docs/installation
   _docs/running
   _docs/api
   _docs/libraries
   _docs/tests

