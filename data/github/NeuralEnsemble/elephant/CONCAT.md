# Elephant - Electrophysiology Analysis Toolkit

[![Build Status](https://travis-ci.org/NeuralEnsemble/elephant.svg?branch=master)](https://travis-ci.org/NeuralEnsemble/elephant)
[![Coverage Status](https://coveralls.io/repos/github/NeuralEnsemble/elephant/badge.svg?branch=master)](https://coveralls.io/github/NeuralEnsemble/elephant?branch=master)
[![Documentation Status](https://readthedocs.org/projects/elephant/badge/?version=latest)](https://elephant.readthedocs.io/en/latest/?badge=latest)
[![![PyPi]](https://img.shields.io/pypi/v/elephant)](https://pypi.org/project/elephant/)
[![Statistics](https://img.shields.io/pypi/dm/elephant)](https://seladb.github.io/StarTrack-js/#/preload?r=neuralensemble,elephant)
[![Gitter](https://badges.gitter.im/python-elephant/community.svg)](https://gitter.im/python-elephant/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge)

*Elephant* package analyses all sorts of neurophysiological data:
spike trains, LFP, analog signals. The input-output data format is either
[Neo](https://github.com/NeuralEnsemble/python-neo), Quantity or Numpy array.


### More Information

* Documentation: https://elephant.readthedocs.io/en/latest/
* Mailing list: https://groups.google.com/group/neuralensemble


#### Visualization of Elephant analysis objects

Viziphant package (https://github.com/INM-6/viziphant) is developed by Elephant
team for plotting and visualization of the output of Elephant functions in a
few lines of code.


#### License
 
Modified BSD License, see [LICENSE.txt](LICENSE.txt) for details.


#### Copyright

:copyright: 2014-2021 by the [Elephant team](doc/authors.rst).


#### Acknowledgments

See [acknowledgments](doc/acknowledgments.rst).


#### Citation

See [citations](doc/citation.rst).

Here, are CSD methods for different electrode configurations.

Keywords: Local field potentials; Current-source density; CSD;
Multielectrode; Laminar electrode; Barrel cortex

1D - laminar probe like electrodes. 
2D - Microelectrode Array like
3D - UtahArray or multiple laminar probes.

The following methods have been implemented here, for, 

1D - StandardCSD, DeltaiCSD, SplineiCSD, StepiCSD, KCSD1D
2D - KCSD2D, MoIKCSD (Saline layer on top of slice)
3D - KCSD3D

Each of these methods listed have some advantages - except StandardCSD which is
not recommended. The KCSD methods can handle broken or irregular electrode
configurations electrode

iCSD
----
Python-implementation of the inverse current source density (iCSD) methods from
http://software.incf.org/software/csdplotter

The Python iCSD toolbox lives on GitHub as well:
https://github.com/espenhgn/iCSD

The methods were originally developed by Klas H. Pettersen, as described in:
Klas H. Pettersen, Anna Devor, Istvan Ulbert, Anders M. Dale, Gaute
T. Einevoll, Current-source density estimation based on inversion of
electrostatic forward solution: Effects of finite extent of neuronal activity
and conductivity discontinuities, Journal of Neuroscience Methods, Volume 154,
Issues 1Ð2, 30 June 2006, Pages 116-133, ISSN 0165-0270,
http://dx.doi.org/10.1016/j.jneumeth.2005.12.005.
(http://www.sciencedirect.com/science/article/pii/S0165027005004541)

To see an example of usage of the methods, see
[demo_icsd.py](https://github.com/espenhgn/iCSD/blob/master/demo_icsd.py)

KCSD 
---- 
This is 1.0 version of kCSD inverse method proposed in

J. Potworowski, W. Jakuczun, S. Łęski, D. K. Wójcik
"Kernel Current Source Density Method"
Neural Computation 24 (2012), 541–575

Some key advantages for KCSD methods are
-- irregular grid of electrodes - accepts arbitrary electrode placement.
-- crossvalidation to ensure no over fitting
-- CSD is not limited to electrode positions - it can obtained at any location

For citation purposes, 
If you use this software in published research please cite the following work
- kCSD1D - [1, 2]
- kCSD2D - [1, 3]
- kCSD3D - [1, 4]
- MoIkCSD - [1, 3, 5]

[1] Potworowski, J., Jakuczun, W., Łęski, S. & Wójcik, D. (2012) 'Kernel
current source density method.' Neural Comput 24(2), 541-575.

[2] Pettersen, K. H., Devor, A., Ulbert, I., Dale, A. M. & Einevoll,
G. T. (2006) 'Current-source density estimation based on inversion of
electrostatic forward solution: effects of finite extent of neuronal activity
and conductivity discontinuities.' J Neurosci Methods 154(1-2), 116-133.

[3] Łęski, S., Pettersen, K. H., Tunstall, B., Einevoll, G. T., Gigg, J. &
Wójcik, D. K. (2011) 'Inverse Current Source Density method in two dimensions:
Inferring neural activation from multielectrode recordings.' Neuroinformatics
9(4), 401-425.

[4] Łęski, S., Wójcik, D. K., Tereszczuk, J., Świejkowski, D. A., Kublik, E. &
Wróbel, A. (2007) 'Inverse current-source density method in 3D: reconstruction
fidelity, boundary effects, and influence of distant sources.' Neuroinformatics
5(4), 207-222.

[5] Ness, T. V., Chintaluri, C., Potworowski, J., Łeski, S., Głabska, H.,
Wójcik, D. K. & Einevoll, G. T. (2015) 'Modelling and Analysis of Electrical
Potentials Recorded in Microelectrode Arrays (MEAs).' Neuroinformatics 13(4),
403-426.

For your research interests of Kernel methods of CSD please see,
https://github.com/Neuroinflab/kCSD-python 

Contact: Prof. Daniel K. Wojcik

Here (https://github.com/Neuroinflab/kCSD-python/tree/master/tests), are
scripts to compare different KCSD methods with different CSD sources. You can
play around with the different parameters of the methods.

The implentation is based on the Matlab version at INCF
(http://software.incf.org/software/kcsd), which is now out-dated. A python
version based on this was developed by Grzegorz Parka
(https://github.com/INCF/pykCSD), which is also not supported at this
point. This current version of KCSD methods in elephant is a mirror of
https://github.com/Neuroinflab/kCSD-python/commit/8e2ae26b00da7b96884f2192ec9ea612b195ec30
---
name: Bug report
about: Report python-elephant bugs/errors here
title: "[Bug] "
labels: ''
assignees: ''

---

**Describe the bug**
<!-- A clear and concise description of what the bug is. -->

**To Reproduce**
1. 
2.

<!-- If you have a code sample, error messages, or stack traces, please provide it here as well. -->


**Expected behavior**
<!-- A clear and concise description of what you expected to happen. -->

**Environment**
 - OS (e.g., Linux):
 - How you installed elephant (`conda`, `pip`, source):
 - Python version:
 - `neo` python package version: 
 - `elephant` python package version:
 - (Any additional python package you want to include here)
---
name: Questions/Help/Support
about: Do you have a question or need help to start running elephant?
title: ''
labels: ''
assignees: ''

---

**Questions and Help**

_Please note that this issue tracker is not a help form and this issue may be closed and redirected to a more appropriate communication channel. Before asking a new question, please check our [NeuralEnsembe forum](https://groups.google.com/forum/#!forum/neuralensemble), where you can discuss any NeuralEnsemble resource, including Elephant._
---
name: Feature request
about: Suggest an idea for this project
title: "[Feature] "
labels: ''
assignees: ''

---

**Feature**
<!-- A clear and concise description of the feature proposal. -->

**Motivation**
<!-- Please outline the motivation for the proposal. Is your feature request related to a problem? E.g., I'm always frustrated when [...]. If this problem is related to another GitHub issue, please link here too. -->

**Alternatives**
<!-- A clear and concise description of any alternative solutions or features you've considered, if any. -->

**Additional context**
<!-- Add any other context about the feature request here. This may be screenshots, online resources, potential experts to contact, scientific publications, demos, etc.-->
