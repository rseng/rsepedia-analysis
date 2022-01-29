[![CI general checks](https://github.com/tylunel/pvpumpingsystem/workflows/CI%20general%20checks/badge.svg)](https://github.com/tylunel/pvpumpingsystem/actions)
[![Coverage](https://codecov.io/gh/tylunel/pvpumpingsystem/branch/master/graph/badge.svg)](https://codecov.io/gh/tylunel/pvpumpingsystem)
[![Documentation Status](https://readthedocs.org/projects/pvpumpingsystem/badge/?version=latest)](https://pvpumpingsystem.readthedocs.io/en/latest/?badge=latest)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tylunel/pvpumpingsystem/master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02637/status.svg)](https://doi.org/10.21105/joss.02637)

![Logo](/docs/images/logo_pvpumpingsystem.jpg)
# pvpumpingsystem
*pvpumpingsystem* is a package providing tools for modeling and sizing
photovoltaic water pumping systems.

![Schema of a PV pumping system](/docs/images/schema_pvps.jpg)

It can model the whole functioning of such pumping system on an hourly basis
and eventually provide key financial and technical findings on a year.
Conversely it can help choose some elements of the pumping station
depending on output values wanted (like daily water consumption and
acceptable risk of water shortage). Further details are provided on the [scope page of the documentation](https://pvpumpingsystem.readthedocs.io/en/latest/package_overview.html).


# Documentation
The full package documentation is available on readthedocs:

[pvpumpingsystem docs](https://pvpumpingsystem.readthedocs.io/en/latest/?badge=latest)


# Installation
*pvpumpingsystem* works with Python 3.5 and superior only.

## With pip

For a rapid installation with pip, type in a command line interface:
```
python -m pip install pvpumpingsystem
```

Consult the docs for more information on installation:
https://pvpumpingsystem.readthedocs.io/en/latest/installation.html


# Hands-on start

Three examples of how the software can be used are in the folder
'docs/examples'.

For a given system, the first two show how to obtain the outflows,
probability of water shortage, life cycle cost and many other results:

[Basic usage example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_basic.ipynb)

[More advanced usage example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/simulation_tunis_advanced.ipynb)

The third shows how to optimize the selection of one or more component
on the pumping station based on user requirements:

[Sizing example](https://nbviewer.jupyter.org/github/tylunel/pvpumpingsystem/blob/master/docs/examples/sizing_example.ipynb)


# Contributions

All kind of contributions (documentation, testing, bug reports,
new features, suggestions...) are highly appreciated.
They can be reported as issues, pull requests, or simple message via
Github (prefered) or via mail of the maintainer.
﻿---
title: 'pvpumpingsystem: A Python package for modeling and sizing photovoltaic water pumping systems'
tags:
  - Python
  - sizing
  - modeling
  - water pumping
  - photovoltaics
  - solar energy
authors:
  - name: Tanguy R. Lunel
    orcid: 0000-0003-3008-1422
    affiliation: "1, 2"
  - name: Daniel R. Rousse
    orcid: 0000-0002-7247-5705
    affiliation: 1
affiliations:
 - name: Industrial Research Group In Technologies of Energy and Energy Efficiency (T3E), Department of Mechanical Sciences, Ecole de Technologie Supérieure Montreal
   index: 1
 - name: Department of Material Science and Engineering, Institut National des Sciences Appliquées Rennes
   index: 2
date: 14 April 2020
bibliography: paper.bib
---

# Summary

According to the World Health Organization, one tenth of the world population still lacks access to
basic water supply. One of the reasons for this is the remoteness of these populations from modern
water collection and distribution technologies, often coupled with an unfavorable socio-economic
situation. Photovoltaic (PV) pumping technology makes it possible to respond both to this problem
and to the criteria of sustainable development. However, these pumping systems must be carefully
modeled and sized in order to make the water supply cost efficient and reliable.

Pvpumpingsystem was conceived in order to tackle this issue. It is an open source package
providing various tools aimed at facilitating the modeling and the sizing of photovoltaic
powered water pumping systems. Even though the package is originally targeted at researchers
and engineers, three practical examples are provided in order to help anyone to use pvpumpingsystem.

Python is the programming language used in the software, and the code is structured with an
object-oriented approach. Continuous integration services allow for lint checking
and to test automation. Each class and function are documented with reference to the
literature when applicable. Pvpumpingsystem is released under a GPL-v3 license.

Pvpumpingsystem relies on already existing packages for photovoltaic and fluid mechanics modeling,
namely “pvlib-python” [@pvlib-python] and “fluids” [@fluids]. pvpumpingsystem's originality lies
in the implementation of various motor-pump models for finite power sources and in the coupling
of the distinct component models. In order to increase the understandability of the code,
each physical component of the PV pumping system corresponds to a class, like for example
the classes `Pump()`, `MPPT()`, `PipeNetwork()`, `Reservoir()`, and `PVGeneration()`. The previous objects
are then gathered in the class `PVPumpSystem()` which allows running a comprehensive model of
the pumping system.

The main inputs to the simulation are an hourly weather file, water source characteristics, expected water
consumption profile, and specifications of the photovoltaic array, motor-pump and water reservoir.
Typical outputs are hourly flow rates, unused electric power, efficiency of components, life
cycle cost and loss of load probability. The sizing module then builds on the
modeling tools, and uses them to provides functions to help choose
the best combination of components in order to minimize the total life cycle cost. Nevertheless,
sizing such complex systems is still an active field of research, and this module is subsequently
expected to be expanded with time.

Two software packages with similar scope already exist: PVsyst and online tool SISIFO, developed by
the MASLOWATEN consortium. Nevertheless, both are closed source, with restricted information
on models used internally, and no API is made accessible. Pvpumpingsystem also has the advantage
of providing ways to size PV pumping systems thanks to automation of pump and PV array choices.

Pvpumpingsystem is the second academic contribution of a broader research program on photovoltaic
water pumping launched in the Technologies of Energy and Energy Efficiency research group at Ecole de Technologie Supérieure Montreal, and is expected to grow with new
features and accuracy assessment provided by experimental studies. The authors also want to give
full access and help to anyone else interested in the use of the software.


# Acknowledgements

The authors would like to acknowledge Mr. Michel Trottier for his generous support to the T3E research group, as well as the NSERC and the FRQNT for their grants and subsidies.
The first author acknowledges the contributions and fruitful discussions with Louis Lamarche and Sergio Gualteros that inspired and helped with the current work.

# Reference
