Contributing to SWIFT-Emulator
==============================

Contributions for the SWIFT-Emulator should come through our GitHub repository,
available at https://github.com/swiftsim/emulator.

Contributions are always welcome, but you should make sure of the following:

+ Your contributions pass all unit tests (you can check this with `pytest`)
+ Your contributions add unit tests for new functionality
+ Your contributions are formatted with the `black` formatter (see `format.sh`)
+ Your contributions are documented fully under `/docs`.

You should also abide by the following code of conduct:

### Code of Conduct

The community of participants in open source Astronomy projects is made up of
members from around the globe with a diverse set of skills, personalities,
and experiences. It is through these differences that our community
experiences success and continued growth. We expect everyone in our community
to follow these guidelines when interacting with others both inside and
outside of our community. Our goal is to keep ours a positive, inclusive,
successful, and growing community.

As members of the community,

+ We pledge to treat all people with respect and provide a harassment- and
  bullying-free environment, regardless of sex, sexual orientation and/or
  gender identity, disability, physical appearance, body size, race,
  nationality, ethnicity, and religion. In particular, sexual language and
  imagery, sexist, racist, or otherwise exclusionary jokes are not appropriate.
+ We pledge to respect the work of others by recognizing
  acknowledgement/citation requests of original authors. As authors, we pledge
  to be explicit about how we want our own work to be cited or acknowledged.
+ We pledge to welcome those interested in joining the community, and realize
  that including people with a variety of opinions and backgrounds will only
  serve to enrich our community. In particular, discussions relating to
  pros/cons of various technologies, programming languages, and so on are
  welcome, but these should be done with respect, taking proactive measure to
  ensure that all participants are heard and feel confident that they can
  freely express their opinions.
+ We pledge to welcome questions and answer them respectfully, paying
  particular attention to those new to the community. We pledge to provide
  respectful criticisms and feedback in forums, especially in discussion
  threads resulting from code contributions.
+ We pledge to be conscientious of the perceptions of the wider community and
  to respond to criticism respectfully. We will strive to model behaviours that
  encourage productive debate and disagreement, both within our community and
  where we are criticized. We will treat those outside our community with the
  same respect as people within our community.
+ We pledge to help the entire community follow the code of conduct, and to
  not remain silent when we see violations of the code of conduct. We will
  take action when members of our community violate this code such as
  contacting josh@joshborrow.com with the subject line SWIFT-emulator Code
  of Conduct (all emails sent in this fashion will be treated with the
  strictest confidence) or talking privately with the person.
+ This code of conduct applies to all community situations online and
  offline, including mailing lists, forums, social media, conferences,
  meetings, associated social events, and one-to-one interactions.

Any related activity or project organized by members of the SWIFT
community, including affiliated packages, are welcome to have their own codes
of conduct, but agree to also abide by the present code of conduct.

Parts of this code of conduct have been adapted from the PSF code of conduct and
the Astropy code of conduct: https://www.astropy.org/code_of_conduct.html.SWIFT-Emulator
==============

[![Documentation Status](https://readthedocs.org/projects/swiftemulator/badge/?version=latest)](https://swiftemulator.readthedocs.io/en/latest/?badge=latest)
![Test Status](https://github.com/swiftsim/emulator/actions/workflows/pytest.yml/badge.svg)
[![PyPI version](https://badge.fury.io/py/swiftemulator.svg)](https://badge.fury.io/py/swiftemulator)
[![status](https://joss.theoj.org/papers/61d082196ef861cc0b612486c1fa6d40/status.svg)](https://joss.theoj.org/papers/61d082196ef861cc0b612486c1fa6d40)

The SWIFT-emulator (henceforth 'the emulator') was initially designed for [SWIFT](http://swift.dur.ac.uk)
outputs, and includes utilities to read and write SWIFT data.

The emulator can be used used to predict 
outputs of simulations without having to run them, by employing Gaussian Process
Regression with `george` and sensitivity analysis with `SALib`.

Dcumentation is available at [ReadTheDocs](https://swiftemulator.readthedocs.io/).

Predicting Simulations
----------------------

The emulator can predict a given scaling relation
(the relationship between two output variables in a simulation, for instance the
masses of galaxies and their size) when varying the underlying physical model
simulated in a continuous way.

As an example from cosmological simulations, imagine varing the energy that supernovae
release when they explode as a parameter `x`. This affects both the sizes and masses of galaxies.
The emulator, using a few 'base' simulations, performed with the real code,
at various values of `x` spaced evenly throughout the viable region, can predict
what the shape of the relationship between mass and size would be at other values
of `x`, given that it has been trained on the base simulation outputs.

Why SWIFT Emulator?
-------------------

The emulator works at a much higher level than other Gaussian Process Regression
libraries, such as `george` (which it is built on).

Working with base machine learning libraries can be tricky, as it typically
requries knowledge of how to structure input arrays in a specific way (both for
training and prediction). They also rarely also include model design routines.
Additionally, validation and visualisation routines
are typically missing.

The SWIFT Emulator package provides a one-stop solution, with a consistent API,
for developing a model design, running it (if using SWIFT), reading in data (
again if using SWIFT), building an emulation model, validating said model,
comparing the model against ground-truth data _across parameter space_ 
(e.g. observations), and visualising the results.

Installation
------------

The package can be installed easily from PyPI under the name `swiftemulator`,
so:

```
pip3 install swiftemulator
```

This will install all necessary dependencies.

The package can be installed from source, by cloning the repository and
then using `pip install -e .` for development purposes.


Requirements
------------

The package requires a number of numerical and experimental design packages.
These have been tested (and are continuously tested) using GitHub actions CI
to use the latest versions available on PyPI. See `requirements.txt` for
details for the packages required to develop SWIFT-Emulator. The packages
will be installed automatically by `pip` when installing from PyPI.


Authors
-------

+ Roi Kugel (@Moyoxkit)
+ Josh Borrow (@jborrow)
---
title: 'swift-emulator: A Python package for emulation of simulated scaling relations'
tags:
  - Python
  - astronomy
  - Simulations
  - Cosmology
  - Machine Learning
authors:
  - name: Roi Kugel^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: "1" # (Multiple affiliations must be quoted)
    orcid: 0000-0003-0862-8639
  - name: Josh Borrow^[co-first author] # note this makes a footnote saying 'co-first author'
    affiliation: 2
    orcid: 0000-0002-1327-1921
affiliations:
 - name: Leiden Observatory, Leiden University, PO Box 9513, NL-2300 RA Leiden, The Netherlands
   index: 1
 - name: Department of Physics, Kavli Institute for Astrophysics and Space Research, Massachusetts Institute of Technology, Cambridge, MA 02139, USA
   index: 2
date: 28 January 2022
bibliography: paper.bib

---

# Summary

`swift-emulator` is a Python toolkit for using Gaussian processes machine
learning to emulate scaling relations from cosmological simulations. 
`swift-emulator` focusses on implementing a clear, easy to use design and API to
remove the barrier to entry for using emulator techniques. `swift-emulator`
provides tools for every step: the design of the parameter sampling, the
training of the Gaussian process model, and validating and anaylsing the trained
emulators. By making these techniques easier to use, in particular in
combination with the SWIFT code [@Schaller2018; @Borrow2020], it will be
possible use fitting methods (like MCMC) to calibrate and better understand
theoretical simulation models.

# Statement of need

One of the limits of doing cosmological (hydrodynamical) simulations is
that any simulation is limited to only a single set of parameters, be these
choices of cosmology, or the implemented physics (e.g., stellar feedback).
These parameters need to be tuned to calibrate against observational data.
At odds with this, cosmological simulations are computationally expensive,
with the cheapest viable runs costing thousands of CPU hours, and running up to
tens of millions for the largest volumes at the highest resolutions.
This makes the use of cosmological simulations in state-of-the-art
fitting pipelines (e.g., MCMC), where tens of thousands to millions of
evaluations of the model are required to explore the parameter space,
computationally unfeasable. In order to get a statistical grip on the models
of cosmology and galaxy formation, a better solution is needed.

This problem is a major limiting factor in "calibration" of the sub-resolution
(subgrid) models that are often used. Works like Illustris [@Vogelsberger2014], 
EAGLE [@Crain2015], BAHAMAS [@McCarthy2017], and Illustris-TNG [@Pillepich2018] are
able to "match" observed relations by eye, but a statistical ground for the
chosen parameters is missing. This poses a signifcant problem for cosmology,
where a deeper understanding of our subgrid models will be required to
interpret results from upcoming surveys like LSST and EUCLID.

A solution here comes through the use of machine learning techniques. Training
'emulators' on a limited amount of simulations enables the evaluation of a
fully continuous model based on changes in the underlying parameters. Instead
of performing a new simulation for each required datapoint, the emulator can
predict the results a simulation would give for that set of parameters. This
makes it feasable to use methods like MCMC based purely on simulation results.

# Emulator Requirements

For emulation in hydro simulations we want to use Gaussian processes to
emulate scaling relations in the following form:

$$GP(y,x,\vec\theta).$$

We want to emulate scaling relations between a dependent variable $y$,
as a function of the independent variable $x$ and the model parameters
$\vec\theta$. For each simulation many of these individual scaling relations can be
calculated, for example the sizes of galaxies relative to their stellar mass,
or the mass fraction of gas in galaxy clusters as a function of their mass. The
individual object properties used in scaling realtions can be measured
from each individual simulation using a tool like VELOCIraptor [@Elahi2019].

Between simulations, the underlying parameters $\vec\theta$ can change,
for instance the energy injected by each supernovae.
Using an emulator, we want to be able to see how many scaling relations
change as a function of these parameters like the supernova strength.

Emulators do not make a distinction between the independent $x$
and the model parameters $\vec\theta$. An emulator will model $y$ as a
function of the combined vector $\vec\theta'=(x,\vec\theta)$. Getting the training
data in the correct format can pose a significant challenge.

In order to save computational time, it is important
to have an efficient sampling of the parameter space represented by $\vec\theta$. 
It may be more efficient to search the parameter space in a transformed
coordinate space, like logarithmic space, if the expected viable range
is over several orders of magnitude.

Once the emulator is working it can be challenging to perform
standard tests to validate it.
Things like cross-checks or parameter sweeps have to be implemented
by hand, making proper use of emulators more difficult.

# Why `swift-emulator`?

Many packages exist for Gausian process emulation, like
`george` (@Ambikasaran2015; this provides the basis for `swift-emulator`),
`gpytorch` [@Gardner2018] and `GPy` [@gpy2014]. Additionally, a package like
`pyDOE` [@pyDOE2012] can be used to set up efficient parameter samplings.
However, most of these packages operate close to theory, and create
a significant barrier for entry.

With `swift-emulator` we aim to provide a single `python` package
that interfaces with available tools at a high level. Additionaly
we aim to streamline the processes by providing i/o tools for the
SWIFT simulation code [@Schaller2018; @Borrow2020]. This is done in a modular
fashion, giving the users the freedom to change any steps along the way.
`swift-emulator` provides many methods that work out of the box,
removing the barrier to entry, and aim at making emulator methods easy to
use. The more wide-spread use of emulators will boost the potential of 
future simulation projects.

`swift-emulator` combines these tools to streamline the complete emulation
process. There are tools for experimental design, such as producing latin
hypercubes or uniform samplings of $n$-dimensional spaces. For simulations
performed with SWIFT, parameter files can be created and simulation outputs can
be loaded in through helper methods in the library. The results can then be used
to train an emulator that can make predictions for the scaling relations in the
simulation. There are also methods to perform cross-checks to find the accuracy
of the emulator. In addition, for investigating the impact of individual
parameters on a given scaling relation, there is a simple method to do a
parameter sweep implemented. Finally, there are tools for comparing the emulated
relations with other data, from a simple $\chi^2$ method to complex model
discrepancy structures.

`swift-emulator` is currently being used for two of the flagship simulation
projects using the SWIFT simulation code, ranging across five orders of 
magnitude in mass resolution. The package is being used to allow modern
simulations to reporduce key observations with high accuracy.

Finally `swift-emulator` has many options to optimise the methods for
specific emulation problems. While the focus so far has been on integration
with SWIFT, the underlying API is structured in a simple enough way that
using the emulator with a different simulation code is easy. `swift-emulator`
is currently being used for simulation projects outside of the SWIFT
project for the calibration of postprocessing models.

# Acknowledgements

We acknowledge support from the SWIFT collaboration whilst developing this
project, with notable involvement from Richard Bower, Ian Vernon, Joop Schaye, 
and Matthieu Schaller. This work is partly funded by Vici grant 639.043.409 from 
the Dutch Research Council (NWO). This work used the DiRAC@Durham facility managed 
by the Institute for Computational Cosmology on behalf of the STFC DiRAC HPC Facility
(www.dirac.ac.uk). The equipment was funded by BEIS capital funding via STFC
capital grants ST/K00042X/1, ST/P002293/1, ST/R002371/1 and ST/S002502/1, Durham
University and STFC operations grant ST/R000832/1. DiRAC is part of the National
e-Infrastructure.

# References
Welcome to SWIFT-Emulator's Documentation
=========================================

The SWIFT-Emulator is a python toolkit for using Gaussian
Process machine learning to produce synthetic simulation data
by interpolating between base outputs. It excels at creating
synthetic scaling relations across large swathes of model
parameter space, as it was created to model galaxy scaling
relations as a function of galaxy formation model parameters
for calibration purposes.

It includes functionality to:

- Generate parameters to perform ground truth runs with in an
  efficient way as a latin hypercube.
- Train machine learning models, including linear models and
  Gaussian Process Regression models (with mean models), on
  this data in a very clean way.
- Generate densly populated synthetic data across the original
  parameter space, and tools to generate complex model
  discrepancy descriptions (known here as penalty functions).
- Generate sweeps across model parameter space for the emulated
  scaling relations to assist in physical insight, as well as
  sensitivity analysis tools based upon raw and synthetic data.
- Validate predictions through cross validation.
- Visualise the resulting penalty data to assist in model choice
  decisions.
- Produce inputs and read outputs from the cosmological code
  SWIFT that processed by VELOCIraptor and the swift-pipeline.

Information about `SWIFT` can be found 
`here <http://swift.dur.ac.uk/>`_, Information about 
`VELOCIraptor` can be found 
`here <https://velociraptor-stf.readthedocs.io/en/latest/>`_
and tnformation about the `SWIFT-pipeline` can be found 
`here <https://github.com/SWIFTSIM/pipeline>`_.

By combining a selection of SWIFT-io and GP analysis 
tools, the SWIFT-Emulator serves to make emulation of 
SWIFT outputs very easy, while staying flexible enough 
to emulate anything, given a good set of training data.

.. toctree::
   :maxdepth: 2

   getting_started/index
   
   emulator_analysis/index

   emulator_options/index
   
   swift_io/index

   comparisons/index

   modules/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
.. API Documentation

API Documentation
=================

.. toctree::
   :maxdepth: 3

   swiftemulator




Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Emulator Options
================

For any smooth function with low dynamic range
the default gaussian process will provide all 
the accuracy that is necesarry for most purposes.
However, it is not hard to imagine certain
situations where taking some extra care before
invoking the GP can lead to more accuracy.
Here, some of the extra options that are
available for the emulation are highlighted.
The example data will be the Schecter function 
example:

.. code-block:: python

    import swiftemulator as se
    from swiftemulator.emulators import gaussian_process
    import numpy as np

    def log_schecter_function(log_M, log_M_star, alpha):
        M = 10 ** log_M
        M_star = 10 ** log_M_star
        return np.log10( (1 / M_star) * (M / M_star) ** alpha * np.exp(- M / M_star ))

    model_specification = se.ModelSpecification(
        number_of_parameters=2,
        parameter_names=["log_M_star","alpha"],
        parameter_limits=[[11.,12.],[-1.,-3.]],
        parameter_printable_names=["Mass at knee","Low mass slope"],
    )

    log_M_star = np.random.uniform(11., 12., 100)
    alpha      = np.random.uniform(-1., -3., 100)

    modelparameters = {}
    for unique_identifier in range(100):
        modelparameters[unique_identifier] = {"log_M_star": log_M_star[unique_identifier],
                                            "alpha": alpha[unique_identifier]}

    model_parameters = se.ModelParameters(model_parameters=modelparameters)

    modelvalues = {}
    for unique_identifier in range(100):
        independent = np.linspace(10,12,10)
        dependent = log_schecter_function(independent,
                                        log_M_star[unique_identifier],
                                        alpha[unique_identifier])
        dependent_error = 0.02 * dependent
        modelvalues[unique_identifier] = {"independent": independent,
                                        "dependent": dependent,
                                        "dependent_error": dependent_error}

    model_values = se.ModelValues(model_values=modelvalues)

    schecter_emulator = gaussian_process.GaussianProcessEmulator()
    schecter_emulator.fit_model(model_specification=model_specification,
                                model_parameters=model_parameters,
                                model_values=model_values)


Mean models
-----------

The most basic addition to a GP is to exchange
the constant (or zero) mean for a more complete
model. For the SWIFT-Emulator these can be found
under :meth:`swiftemulator.mean\_models`. All 
currently implement models come in the form of
different order polynomials. Much like the GP
you have to define your model first. The model can
then be passed to the GP which will fit the
coefficients and use that as a mean model.
The GP will then be used to model the residuals
between the polynomial fit and the data.

.. code-block:: python

    from swiftemulator.mean_models.polynomial import PolynomialMeanModel

    polynomial_model = PolynomialMeanModel(degree=1)

    schecter_emulator = gaussian_process.GaussianProcessEmulator(mean_model=polynomial_model)
    schecter_emulator.fit_model(model_specification=model_specification,
                            model_parameters=model_parameters,
                            model_values=model_values)

The polynomial model fits a polynomial surface
to all parameters. This includes not just the 
polynomial coefficients for each parameter but
also the linear combinations up to the degree
of the model. Be carefull picking a degree that
is very large, as it can quickly lead to 
over-fitting.

Emulating Bin-by-Bin
--------------------

The scaling relations obtained from simulations
are often binned relations. For the all-purpose
emulator we use the `independent` as an additional
parameter which the emulator uses for prediction,
but there are situation where you would prefer
modeling the response at each `independent` bin
seperately, instead of modeling it all at once.

When emulating bin-to-bin the main difference
is that your GP now comes from
:meth:`swiftemulator.emulators.gaussian_process_bins`.
Each bin will have a unique emulator, that is
trained on all data available for that bin.
A bin is created for each unique value of
the `independent` found in the `ModelValues`
container. It is extremely important that
each bin of the original data-set has exactly 
the same value for the `independent`. However,
the individual models do not need the same
sample of bins. If some models are missing
values for some of the bins, this is not
a problem.

Using the binned emulator is as simple as

.. code-block:: python

    from swiftemulator.emulators import gaussian_process_bins

    schecter_emulator_binned = gaussian_process_bins.GaussianProcessEmulatorBins()
    schecter_emulator_binned.fit_model(model_specification=model_specification,
                                       model_parameters=model_parameters,
                                       model_values=model_values)

Which has the same predicion functionality
as the standard `gaussian_process`.
Note that there is also a binned version
of the cross checks, 
:meth:`swiftemulator.sensitivity.cross\_check\_bins`,
which acts the same as the normal `cross_check`
but instead uses the binned emulator, making
it easy to compare the two methods.

1D Emulation
------------

Sometimes the emulation problem is better solved as

.. math::
    f(\vec\theta)

In this case we only have the model parameters.
the emulator won't be a function of an additional x
parameter stored in the model values. In this case the
use can use :meth:`swiftemulator.emulators.gaussian_process_one_dim`.
This method has similar functionality as the other
emulator types. It will still need a ModelValues
container. Here is an example of how such a container
should look like:

.. code-block:: python

    modelvalues = {}
    for unique_identifier in range(100):
        dependent = func(a_arr[unique_identifier], b_arr[unique_identifier])
        dependent_error = 0.02 * dependent
        modelvalues[unique_identifier] = {"independent": [None],
                                          "dependent": [dependent],
                                        "dependent_error": [dependent_error]}

In order to make use of the general emulator
containers, it is still required to provide the values
as list. In this case the lists will only contain a single
value. The independent value will not be read. When your
data is in the correct format the emulator can be trained
like all the other methods.

.. code-block:: python

    from swiftemulator.emulators import gaussian_process_one_dim

    schecter_emulator_one_dim = gaussian_process_one_dim.GaussianProcessEmulator1D()
    schecter_emulator_one_dim.fit_model(model_specification=model_specification,
                                       model_parameters=model_parameters,
                                       model_values=model_values)

The only other thing of note is that while
`predict_values` retains the same functionality,
you are no longer required to specify any independent
values. The prediction is now based purely of the
given values of the model parameters.Comparing With Data
===================

To unlock the full power of emulation it
is often usefull to compare your results
with observational data. With
SWIFT-Emulator you can directly compare
the emulated outputs with the observational
data stored in `velociraptor-comparison-data`
which can be found 
`here <https://github.com/SWIFTSIM/velociraptor-comparison-data>`_.
You can also download just the required ``Vernon.hdf5`` file
`instead <http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/Vernon.hdf5>`_.

In this case we will again use the data from
`http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/emulator_output.zip`
First we have to set up the emulator

.. code-block:: python

    from swiftemulator.io.swift import load_parameter_files, load_pipeline_outputs
    from swiftemulator.emulators.gaussian_process import GaussianProcessEmulator
    from swiftemulator.mean_models import LinearMeanModel
    from velociraptor.observations import load_observations

    from glob import glob
    from pathlib import Path
    from tqdm import tqdm
    from matplotlib.colors import Normalize

    import matplotlib.pyplot as plt
    import numpy as np
    import corner

    import os

    files = [Path(x) for x in glob("./emulator_output/input_data/*.yml")]

    filenames = {filename.stem: filename for filename in files}

    spec, parameters = load_parameter_files(
        filenames=filenames,
        parameters=[
            "EAGLEFeedback:SNII_energy_fraction_min",
            "EAGLEFeedback:SNII_energy_fraction_max",
            "EAGLEFeedback:SNII_energy_fraction_n_Z",
            "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",
            "EAGLEFeedback:SNII_energy_fraction_n_n",
            "EAGLEAGN:coupling_efficiency",
            "EAGLEAGN:viscous_alpha",
            "EAGLEAGN:AGN_delta_T_K",
        ],
        log_parameters=[
            "EAGLEAGN:AGN_delta_T_K",
            "EAGLEAGN:viscous_alpha",
            "EAGLEAGN:coupling_efficiency",
            "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",
        ],
        parameter_printable_names=[
            "$f_{\\rm E, min}$",
            "$f_{\\rm E, max}$",
            "$n_{Z}$",
            "$\\log_{10}$ $n_{\\rm H, 0}$",
            "$n_{n}$",
            "$\\log_{10}$ $C_{\\rm eff}$",
            "$\\log_{10}$ $\\alpha_{\\rm V}$",
            "AGN $\\log_{10}$ $\\Delta T$",
        ],
    )

    value_files = [Path(x) for x in glob("./emulator_output/output_data/*.yml")]

    filenames = {filename.stem: filename for filename in value_files}

    values, units = load_pipeline_outputs(
        filenames=filenames,
        scaling_relations=["stellar_mass_function_100"],
        log_independent=["stellar_mass_function_100"],
        log_dependent=["stellar_mass_function_100"],
    )

    # Train an emulator for the space.
    scaling_relation = values["stellar_mass_function_100"]
    scaling_relation_units = units["stellar_mass_function_100"]

    emulator = GaussianProcessEmulator()
    emulator.fit_model(model_specification=spec,
        model_parameters=parameters,
        model_values=scaling_relation,
    )


In this case we are gonna look at the stellar
mass function. To compare we load the calibration
SMF for EAGLE-XL.

.. code-block:: python

    observation = load_observations(
        "../velociraptor-comparison-data/data/GalaxyStellarMassFunction/Vernon.hdf5"
    )[0]

Penalty Functions
-----------------

There is a large selection of "Penalty" functions
available. We define a penalty function as an
analogous to a likelihood.

.. math::
    \mathcal{L} = 1 -  P(x,\theta),

where :math:`\mathcal{L}` is the likelihood and
:math:`P(x,\theta)` is the accompanying penalty
function.

As an example we will use an L2 norm. This will
calculate the mean squared distance between the
emulator and the data. 

.. code-block:: python

    from swiftemulator.comparison.penalty import L2PenaltyCalculator
    from unyt import Msun, Mpc

    L2_penalty = L2PenaltyCalculator(offset = 0.5, lower=9,upper=12)
    L2_penalty.register_observation(observation,log_independent=True
                                ,log_dependent=True
                                ,independent_units=Msun
                                ,dependent_units=Mpc**-3)

    L2_penalty.plot_penalty(9,12,-6,-1,"penalty_example",x_label="Stellar mass",y_label="dn/dlogM")

.. image:: penalty_example.png

Now we can combine this with the emulator to compare models
in terms of how good they fit the data. Without using the
emulator we can use interpolation to be able to quickly check
which node of the parameter space best fits the data via
:meth:`swiftemulator.comparison.penalty.L2PenaltyCalculator.penalties`

.. code-block:: python

    all_penalties = L2_penalty.penalties(emulator.model_values,np.mean)

    all_penalties_array = []
    node_number = []
    for key in all_penalties.keys():
        all_penalties_array.append(all_penalties[key])
        node_number.append(int(key))
        
    print("Best fit node = ",node_number[np.argmin(all_penalties_array)])

.. code-block:: python

    Best fit node =  107

If we want to check the simulation that is best without rerunning
anything we can use node 107. In general we can use this to check
not just models at the nodes, but use the emulator to check the
complete parameter range. Starting with node 107, let's see if
we can improve the fit by chaning one of the parameters.

.. code-block:: python

    predictparams = emulator.model_parameters["107"].copy()
    x_to_predict = np.log10(L2_penalty.observation.x.value)
    pred, pred_var = emulator.predict_values(x_to_predict, predictparams)

    print("Mean Penalty of node 107 = ",np.mean(L2_penalty.penalty(x_to_predict,pred)))

    #Let's change one of the parameters and see if it improves the fit
    predictparams["EAGLEFeedback:SNII_energy_fraction_max"] = 1
    x_to_predict = np.log10(L2_penalty.observation.x.value)
    pred, pred_var = emulator.predict_values(x_to_predict, predictparams)

    print("Mean after change = ",np.mean(L2_penalty.penalty(x_to_predict,pred)))

.. code-block:: python

    Mean Penalty of node 107 =  0.21988119507121354
    Mean after change =  0.3344361855742612

This change makes the fit worse, so no luck. In general you would
not do this by hand, but use for example MCMC to sample all the
parameters.

Defining New Penalty Functions
------------------------------

What you want out of these penalty functions can vary wildy,
but it is very easy to define your own. There is a large set 
of functions available within
:meth:`swiftemulator.comparison.penalty`. It is also possible
to add your own functions. The base class
:meth:`swiftemulator.comparison.penalty.PenaltyCalculator`
covers the most important part, which is loading and
interpolating the data. You can then add whichever calculattion
of the penalties you want. In the example below we create a
function that is Gaussian weighted, with a constent error
term.

.. code-block:: python

    from swiftemulator.comparison.penalty import PenaltyCalculator
    import unyt

    class ExamplePenaltyCalculator(PenaltyCalculator):
        
        def penalty(self,independent, dependent, dependent_error):
            #We can use the observational data from the base class.
            #We calculate the observational y-values to compare with
            #from the interpolated observations.
            obs_dependent = self.interpolator_values(independent)
            
            penalties = np.exp(-np.abs(dependent - obs_dependent)**2/0.1)
            return penalties
        
    my_penalty = ExamplePenaltyCalculator()
    my_penalty.register_observation(observation,log_independent=True,log_dependent=True
                                ,independent_units=Msun,dependent_units=Mpc**-3)

    my_penalty.plot_penalty(9,12,-6,-1,"my_penalty",x_label="Stellar mass",y_label="dn/dlogM")

.. image:: example_penalty_example.png

For the simplest models you can also still use the `plot_penalty`
functionality. There are also PF's available that use the
errors on the data, for example
:meth:`swiftemulator.comparison.penalty.GaussianDataErrorsPenaltyCalculator`.
When creating new penalty functions you can use different parts 
of already existing ones to make the process very easy.
.. _getting_started:

Getting Started
===============

In this section the basics of using the SWIFT-Emulator will
be explained, with examples of how to make your first GP 
predictions.

Installation
------------

The package can be installed easily from PyPI under the name `swiftemulator`,
so:

``pip3 install swiftemulator``

This will install all necessary dependencies.

The package can be installed from source, by cloning the repository and
then using `pip install -e .` for development purposes.

Requirements
------------

The package requires a number of numerical and experimental design packages.
These have been tested (and are continuously tested) using GitHub actions CI
to use the latest versions available on PyPI. See `requirements.txt` for
details for the packages required to develop SWIFT-Emulator. The packages
will be installed automatically by `pip` when installing from PyPI.

Loading data
------------

In the way we set up the emulator, loading the data is the
most cumbersome part of emulation. Once everything is in the
right format the emulation itself will be very easy.

At the basis of the SWIFT-Emulator lies the ability to train
a Gaussian process (GP) based on a set of training data. As 
the main goal is emulating scaling relations on the back of 
hydro simulations you should think of the emulation being in 
the following form

.. math::
    GP(y,x,\theta),

where we want to predict the dependent :math:`y` as a function 
of the independent :math:`x` and model parameters :math:`\theta`.
The distinction between :math:`x` and :math:`\theta` is made
to distinquish the relation that can be obtained from a single
simulation output (Like the number density of galaxies as a 
function of their stellar mass, where y is the number density
and x the stellar mass) from the parameters that span different
outputs (Like redshift, or AGN feedback strength). For this
example we will predict the stellar mass function using some
data generated with a Schecter function

.. code-block:: python

    import swiftemulator as se
    import numpy as np

    def log_schecter_function(log_M, log_M_star, alpha):
        M = 10 ** log_M
        M_star = 10 ** log_M_star
        return np.log10( (1 / M_star) * (M / M_star) ** alpha * np.exp(- M / M_star ))

where we set the normalisation to unity. In this case we will 
use `log_M` as the independent, while `M_star` and `alpha` are
our model parameters. The choice to emulate in log space is 
important, this massively decreases the dynamic range which
makes it a lot easier for a Gaussian process to emulate
accurately.

In order to get the data in the correct form we need to
define three containers. We start by specifying some of the 
basic information of our model. This is done via :meth:`swiftemulator.backend.model\_specification`.

.. code-block:: python

    model_specification = se.ModelSpecification(
        number_of_parameters=2,
        parameter_names=["log_M_star","alpha"],
        parameter_limits=[[11.,12.],[-1.,-3.]],
        parameter_printable_names=["Mass at knee","Low mass slope"],
    )

The mode specification is used to store some of the metadata
of the training set.

Lets assume our training set consists of 100 simulations where
our model parameters are randomly sampled

.. code-block:: python

    log_M_star = np.random.uniform(11., 12., 100)
    alpha      = np.random.uniform(-1., -3., 100)

This can be used to set up the second container, which 
contains the model parameters. For each unique model, named
by `unique_identifier`, we store the values in a dictionary

.. code-block:: python

    modelparameters = {}
    for unique_identifier in range(100):
        modelparameters[unique_identifier] = {"log_M_star": log_M_star[unique_identifier],
                                            "alpha": alpha[unique_identifier]}
        
    model_parameters = se.ModelParameters(model_parameters=modelparameters)

The `unique_identifier` is really important, as this will be used to
link the model parameters to the model values. There are some major
advantages to splitting this up. By splitting the model from the 
scaling relation we only have to define the model once, allowing
us to attach as many relations to it as we want.

Final thing is adding the values of the function we want to emulate.
This is done in a similar way to how we add the model parameters,
except that we now attach a complete array for each model.

.. code-block:: python

    modelvalues = {}
    for unique_identifier in range(100):
        independent = np.linspace(10,12,10)
        dependent = log_schecter_function(independent, 
                                          log_M_star[unique_identifier], 
                                          alpha[unique_identifier])
        dependent_error = 0.02 * dependent
        modelvalues[unique_identifier] = {"independent": independent, 
                                        "dependent": dependent, 
                                        "dependent_error": dependent_error}
        
    model_values = se.ModelValues(model_values=modelvalues)

For the model values it is important that you use the names
`independent`, `dependent` and `dependent_error` for `x`, `y`
and `y_err` respectively. These specific names are used when 
setting up the emulator

Training the emulator
---------------------

After setting up de model containers, training the emulator becomes
very simple. First we create an empty GP, which we can then train
on the data we have just loaded.

.. code-block:: python

    from swiftemulator.emulators import gaussian_process
    schecter_emulator = gaussian_process.GaussianProcessEmulator()
    schecter_emulator.fit_model(model_specification=model_specification,
                                model_parameters=model_parameters,
                                model_values=model_values)

This might take a little bit of time. At this point the GP is
fully trained and can be used to make predictions. There are 
a lot more options when setting up the GP, like indlucing a 
model for mean, but if your input is smooth this is likely
all you will need.

Making predictions
------------------

The real reason to use an emulator is to eventually predict
the shape of the scaling relation continuously over the
parameterspace. Just like training the emulator, making
predictions is extremely simple

.. code-block:: python

    predictparams = {"log_M_star": 11.5, "alpha": -2}
    predict_x = np.linspace(10,12,100)

    pred, pred_var = schecter_emulator.predict_values(predict_x, predictparams)


The main thing to keep in mind is that you give the
model parameters as a dictionary again, with the same 
names as how they are defined in the `model_parameters`.
In this case we can directly compare with the original
model.

.. code-block:: python

    import matplotlib.pyplot as plt
    plt.plot(predict_x,pred,label="Emulator")
    plt.plot(predict_x,log_schecter_function(predict_x,
                                            predictparams["log_M_star"]
                                            ,predictparams["alpha"])
            ,color="black",ls=":",label="Model")
    plt.xlabel("Stellar mass")
    plt.ylabel("dn/dlogM")
    plt.legend()

Which shows that the emulator can predict the model with
high accuracy.

.. image:: predict_vs_model.png

This covers the most basic way to use SWIFT-Emulator and
should give a good baseline for using some of the
additional features it offers.Experimental Design
-------------------

One part of SWIFT-Emulator's i/o features is
the ability to generate Latin Hypercube (LH)
designs and save them to SWIFT parameter files,
one file for each set of parameters. For this
example you will need to download some data
from `http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/emulator_output.zip`.

We do this by combining the :meth:`swiftemulator.design`
with :meth:`swiftemulator.io.swift`. First we have 
to specify what parameter we want to vary.

.. code-block:: python

    from swiftemulator.design import latin
    from swiftemulator.io.swift import write_parameter_files
    from swiftemulator import ModelSpecification

    spec = ModelSpecification(
        number_of_parameters=5,
        parameter_names=[
            "EAGLEFeedback:SNII_energy_fraction_min",
            "EAGLEFeedback:SNII_energy_fraction_max",
            "EAGLEFeedback:SNII_energy_fraction_n_Z",
            "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",
            "EAGLEFeedback:SNII_energy_fraction_n_n",
        ],
        parameter_printable_names=[
            "$f_{\\rm E, min}$",
            "$f_{\\rm E, max}$",
            "$n_{Z}$",
            "$\\log_{10}$ $n_{\\rm H, 0}$",
            "$n_{n}$",
        ],
        parameter_limits=[
            [0.0, 1.0],
            [1.0, 7.0],
            [-0.5, 5.0],
            [-1.0, 1.5],
            [-0.5, 5.0],
        ],
    )

    parameter_transforms = {"SNII_energy_fraction_n_0_H_p_cm3": lambda x: 10.0 ** x}

In this case it is important that your
`parameter_names` are identical to the
names in the SWIFT parameter file. The parameter
file is a `.yml` file so the individual parameters
should be named in that format. 

In this case the fourth parameter, 
`EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3`
is sampled in log-space. The `parameter_limits`
are given in log-space as well in this case, but
you need to define the transformation needed when
going from the design space, to the value you
want to put in the parameter file.

Generating the LH can then be done with
:meth:`swiftemulator.design.latin.create\_hypercube`.

.. code-block:: python

    number_of_simulations = 30

    model_parameters = latin.create_hypercube(
        model_specification=spec,
        number_of_samples=number_of_simulations,
    )

Now we can use the SWIFT i/o to write these
to a set of parameter files. You will have
noticed that we only need to provide the
parameters that we want to vary. This is 
because we provide `write_parameter_files`
with a base parameter file. This file
should hold the base values for all
other parameters. Example parameter files
can be found on the main `SWIFT` repository.
For this example we will use one of the
parameter files from the example hypercube.

.. code-block:: python

    base_parameter_file = "emulator_output/input_data/1.yml"
    output_path = "."

    write_parameter_files(
    filenames={
        key: f"{output_path} / {key}.yml"
        for key in model_parameters.model_parameters.keys()
    },
    model_parameters=model_parameters,
    parameter_transforms=parameter_transforms,
    base_parameter_file=base_parameter_file,
    )

This writes 30 files to the current
directory. These files can then be used to run
SWIFT for each of the models.Loading SWIFT data
------------------

In order for this example to work you will need to
download some data from
`http://virgodb.cosma.dur.ac.uk/swift-webstorage/IOExamples/emulator_output.zip`.
This will contain a set of parameter files and a 
set of data files. The parameters files will be
used to retrieve the parameters of the Latin 
Hypercube, while the data files contain the
results for the scaling relations for each model.

It is adviced to use the SWIFT-io options if you
want to compare directly with observational data.
The main advantage being that loading the data in
this way will ensure that you use the correct units.

What will be required to load the data is a list 
of all the parameter files, and a list for all
the data files. This can easily be obtained using
:mod:`glob` and :mod:`Path`.

.. code-block:: python

    from glob import glob
    from pathlib import Path

    parameter_files = [Path(x) for x in glob("./emulator_output/input_data/*.yml")]
    parameter_filenames = {filename.stem: filename for filename in parameter_files}

    data_files = [Path(x) for x in glob("./emulator_output/output_data/*.yml")]
    data_filenames = {filename.stem: filename for filename in data_files}

For the parameter we use
:meth:`swiftemulator.io.swift.load_parameter_files`.
This reads in the parameters and returns both a `ModelSpecification`
and a `ModelParameters` container to pass to the emulator.

.. code-block:: python

    from swiftemulator.io.swift import load_parameter_files

    spec, parameters = load_parameter_files(
    filenames=parameter_filenames,
    parameters=[
        "EAGLEFeedback:SNII_energy_fraction_min",
        "EAGLEFeedback:SNII_energy_fraction_max",
        "EAGLEFeedback:SNII_energy_fraction_n_Z",
        "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",
        "EAGLEFeedback:SNII_energy_fraction_n_n",
        "EAGLEAGN:coupling_efficiency",
        "EAGLEAGN:viscous_alpha",
        "EAGLEAGN:AGN_delta_T_K",
    ],
    log_parameters=[
        "EAGLEAGN:AGN_delta_T_K",
        "EAGLEAGN:viscous_alpha",
        "EAGLEAGN:coupling_efficiency",
        "EAGLEFeedback:SNII_energy_fraction_n_0_H_p_cm3",
    ],
    parameter_printable_names=[
        "$f_{\\rm E, min}$",
        "$f_{\\rm E, max}$",
        "$n_{Z}$",
        "$\\log_{10}$ $n_{\\rm H, 0}$",
        "$n_{n}$",
        "$\\log_{10}$ $C_{\\rm eff}$",
        "$\\log_{10}$ $\\alpha_{\\rm V}$",
        "AGN $\\log_{10}$ $\\Delta T$",
    ],
    )

Just like for the experimental design, it is
important that the name used for `parameters`
is the same as the one used in the parameter 
file. Note also that you have to supply a list
of parameters that where sampled in log-space.
These are then trasformed to log-space before
being stored in a `ModelParameters` container.

To read the `ModelValues` the function 
:meth:`swiftemulator.io.swift.load_pipeline_outputs`
is used. In this case you have to supply the
filenames, and the name(s) of the scaling relation(s).
These names can be easily found in the data file 
and are set by the config used for the pipeline.
`log_independent` and `log_dependent` will cause
the x or y to be loaded in log-space.

.. code-block:: python

    from swiftemulator.io.swift import load_pipeline_outputs

    values, units = load_pipeline_outputs(
        filenames=data_filenames,
        scaling_relations=["stellar_mass_function_100"],
        log_independent=["stellar_mass_function_100"],
        log_dependent=["stellar_mass_function_100"],
    )

    scaling_relation = values["stellar_mass_function_100"]
    scaling_relation_units = units["stellar_mass_function_100"]

`load_pipeline_outputs` can return as many scaling
relations as required. `values` is dictionary that
contains a `ModelValues` container for each requested
scaling relation. A `ModelValues` container for a
single relation can be obtained by parsing it with
the correct name.

At this point the data is loaded and you can build
and train your emulator.

.. code-block:: python

    from swiftemulator.emulators import gaussian_process

    emulator = gaussian_process.GaussianProcessEmulator()
    emulator.fit_model(model_specification=spec,
        model_parameters=parameters,
        model_values=scaling_relation,
    )Interaction with SWIFT
======================

The SWIFT-Emulator is designed to easily generate
:mod:`SWIFT` parameter files, and read
:mod:`swift-pipeline` outputs. Here you can find
a description of how to use the `SWIFT-Emulator`'s
tools to set up an experimental design, and then
load the resulting scaling relations.

.. toctree::
    :maxdepth: 2
 
    design
    data_loadingAnalysis tools
==============

Here we will outline some of the available
tools that can help inspect the performance
of the emulator. The example data will be
the Schecter function example:

.. code-block:: python

    import swiftemulator as se
    from swiftemulator.emulators import gaussian_process
    import numpy as np

    def log_schecter_function(log_M, log_M_star, alpha):
        M = 10 ** log_M
        M_star = 10 ** log_M_star
        return np.log10( (1 / M_star) * (M / M_star) ** alpha * np.exp(- M / M_star ))

    model_specification = se.ModelSpecification(
        number_of_parameters=2,
        parameter_names=["log_M_star","alpha"],
        parameter_limits=[[11.,12.],[-1.,-3.]],
        parameter_printable_names=["Mass at knee","Low mass slope"],
    )

    log_M_star = np.random.uniform(11., 12., 100)
    alpha      = np.random.uniform(-1., -3., 100)

    modelparameters = {}
    for unique_identifier in range(100):
        modelparameters[unique_identifier] = {"log_M_star": log_M_star[unique_identifier],
                                            "alpha": alpha[unique_identifier]}

    model_parameters = se.ModelParameters(model_parameters=modelparameters)

    modelvalues = {}
    for unique_identifier in range(100):
        independent = np.linspace(10,12,10)
        dependent = log_schecter_function(independent,
                                        log_M_star[unique_identifier],
                                        alpha[unique_identifier])
        dependent_error = 0.02 * dependent
        modelvalues[unique_identifier] = {"independent": independent,
                                        "dependent": dependent,
                                        "dependent_error": dependent_error}

    model_values = se.ModelValues(model_values=modelvalues)

    schecter_emulator = gaussian_process.GaussianProcessEmulator()
    schecter_emulator.fit_model(model_specification=model_specification,
                                model_parameters=model_parameters,
                                model_values=model_values)

Cross checks
------------

To set up a cross check, the emulator is
trained on all but one of the input data-sets.
The resulting emulator can then be compared
against the model that was left out.

Cross checks are the main way of quantifying
emulator performance in the absence of validation
data. When emulating via cosmological simulations
it is likely to be very expensive to generate a 
validation dataset of sufficient size. for cases
like this SWIFT-Emulator has an easy way of setting up
cross-checks.

The :meth:`swiftemulator.sensitivity.cross\_check`
object acts identically to :meth:`swiftemulator.emulators.gaussian\_process`
and takes the same inputs. By setting the cross-checks
up in this way you can directly compare the results
with the main GP that you use for predictions.

.. code-block:: python

    from swiftemulator.sensitivity import cross_check

    schecter_ccheck = cross_check.CrossCheck()
    schecter_ccheck.build_emulators(model_specification=model_specification,
                            model_parameters=model_parameters,
                            model_values=model_values)

In this case `build_emulators` takes the place of `fit_model`.
Note that build_emulators now creates N independent trained
emulators, where N is the number of models, so this can take
quite a long time. For this example the amount of models was
reduced from 100 to 20.

Once the emulators have been build there are some inherent
tools to have a look at the result (see :meth:`swiftemulator.sensitivity.cross\_check`).
We will use `build_mocked_model_values_original_independent()`
to compare the cross-check predictions with the original
data.

.. code-block:: python

    import matplotlib.pyplot as plt

    data_by_cc = schecter_ccheck.build_mocked_model_values_original_independent()

    for unique_identifier in range(20):
        cc_over_og = data_by_cc[unique_identifier]["dependent"] / \
                    model_values[unique_identifier]["dependent"]
        plt.plot(data_by_cc[unique_identifier]["independent"],cc_over_og)
        plt.xlabel("Mass")
        plt.ylabel("Cross-check / Truth")
        
    plt.savefig("Cross_check_accuracy.png",dpi=200)

.. image:: Cross_check_accuracy.png

Just with a few line we are able to quantify how accurate
the emulator is. Also note that any `ModelValues` container
can be parsed as if it is a dictionary.

Sweeps Of Parameter Space
-------------------------

One of the advantages of using emulators is that it supplies
you with a fully continuous model of the given function.
Besides fitting the parameters it is often interesting to see
the effect of changing a single parameter, by doing a sweep.

This is implemented into the SWIFT-Emulator with 
:meth:`swiftemulator.mocking.mock\_sweep`.

.. code-block:: python

    from swiftemulator.mocking import mock_sweep

    center = {"log_M_star": 11.5, "alpha": -2.0}

    Mock_values, Mock_parameters = mock_sweep(schecter_emulator
                                          ,model_specification
                                          ,6,"alpha",center)

    for mock_name in Mock_values.keys():
        plt.plot(Mock_values[mock_name]["independent"],
                Mock_values[mock_name]["dependent"],
                label = "Alpha = " +str(Mock_parameters[mock_name]["alpha"])[:4])

    plt.xlabel("Stellar mass")
    plt.ylabel("dn/dlogM")    
    plt.legend()
    plt.savefig("parameter_sweep.png",dpi=200)

.. image:: parameter_sweep.png

`mock_sweep` returns the values and parameter of the 
sweep as `ModelValues` and `ModelParameters`
containers, that are easy to parse. 

Model Parameters Features
-------------------------

This highlights two small functions that are attached to
the :meth:`swiftemulator.backend.model\_parameters`
object. The first is the ability to generate a quick plot
of the experimental design using :mod:`corner`.

.. code-block:: python

    model_parameters.plot_model(model_specification)

.. image:: experimental_design.png

Note that the axis label used here are the one passed to
the model specification. This can be used to have a quick
look at whether your space is well sampled.

After finding a set of best fit model parameters it is
sometimes usefull to see if there are any individual model
that has similar values. `find_closest_model` takes a
dictionary of input values and finds the training model
that is closets to those values. 

.. code-block:: python

    best_model = {"log_M_star": 11.3, "alpha": -2.1}

    model_parameters.find_closest_model(best_model,number_of_close_models=5)

which outputs

.. code-block:: python

    ([2, 12, 18, 19, 3],
    [{'log_M_star': 11.26347510702813, 'alpha': -1.9614226414699145},
    {'log_M_star': 11.507944778215956, 'alpha': -1.9818583963792449},
    {'log_M_star': 11.19527147203741, 'alpha': -1.8330160108907092},
    {'log_M_star': 11.033961506507945, 'alpha': -2.275313906753826},
    {'log_M_star': 11.67912812994198, 'alpha': -2.0664526312834353}])

It returns a list with the `unique_identifier` of each close
model, and the model parameters belonging to that model. This
can be used to explore the models close to you best fit model,
for example to check how well sampled that part of parameter
space is.

Checking Hyperparameters
------------------------

In general one should not look at the hyperparameters. They
should only be used as a diagnostic when the emulator is
giving strange results. The SWIFT-Emulator provides an
easy way to check the parameterspace of the hyperparameters.
The hyperparameters are optimised to using the
marginalised likelihood, so we can inspect how well converged
they are by looking at the probability distribution of each
individual hyperparameter. This is done via
:meth:`swiftemulator.emulators.gaussian\_process\_mcmc`.
In this case MCMC implies the use of Markov chian
Monte Carlo (via :mod:`emcee`) to find the best
hyperparameters, allowing us to look at the complete
parameter space.

.. code-block:: python

    from swiftemulator.emulators import gaussian_process_mcmc
    schecter_emulator_mcmc = gaussian_process_mcmc.GaussianProcessEmulatorMCMC(burn_in_steps=1
                                                                              ,mcmc_steps=1000)
    schecter_emulator_mcmc.fit_model(model_specification=model_specification,
                            model_parameters=model_parameters,
                            model_values=model_values)

    schecter_emulator_mcmc.plot_hyperparameter_distribution()

.. image:: hyperparameters.png

This method is a lot slower than the default hyperparameter
optimisation, and may take some time to compute. The main
take away from plots like this is to see whether the
hyperparameters are converged, and whether they are 
consistent with the faster optimisation method.