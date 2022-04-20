## Table of Contents
- [Introduction](https://github.com/UNSW-CEEM/nempy#introduction)
- [Motivation](https://github.com/UNSW-CEEM/nempy#motivation)
- [Documentation](https://github.com/UNSW-CEEM/nempy#documentation)
- [Community](https://github.com/UNSW-CEEM/nempy#community)
- [Installation](https://github.com/UNSW-CEEM/nempy#installation)
- [A simple example](https://github.com/UNSW-CEEM/nempy#a-simple-example)
- [A detailed example](https://github.com/UNSW-CEEM/nempy#a-detailed-example)

## Introduction
Nempy is a python package for modelling the dispatch procedure of the Australian National Electricity Market (NEM). The idea is 
that you can start simple, like in the [example below](https://github.com/UNSW-CEEM/nempy#a-simple-example), and grow the complexity of your model by adding features such as 
ramping constraints, interconnectors, FCAS markets and more.

A brief introduction to the NEM can be found at the link below:

https://aemo.com.au/-/media/Files/Electricity/NEM/National-Electricity-Market-Fact-Sheet.pdf

## Documentation
A more detailed introduction to Nempy, examples, and reference documentation can be found on the 
[readthedocs](https://nempy.readthedocs.io/en/latest/) page.

## Community
Nempy is open-source and we welcome all forms of community engagement, some more info is provided below.

### Support
You can seek suport for using Nempy using the discussion tab on GitHub (https://github.com/UNSW-CEEM/nempy/discussions), checking the issues register (https://github.com/UNSW-CEEM/nempy/issues), or by contacting Nick directly n.gorman at unsw.edu.au.

### Reporting issues
Issues with Nempy can be reported via the issues register (https://github.com/UNSW-CEEM/nempy/issues), issues submissions do not need to adhere to any particular format.

### Contributing
Contributions via pull requests are welcome. Contributions should; follow the PEP8 style guide (with exception of line length up to 120 rather than 80), ensure that all existing automated tests continue to pass (unless you are explicitly changing intended behavour, please highlight this in your pull request description), implement automated tests for new features, and provided doc strings for public interfaces. 

### Get in touch
This project is being lead by Nick Gorman, a PhD candidate at the University for New South Wales and the Collaboration 
on Energy and Environmental Markets (CEEM). As part the project we hope to engage with and support prospective users of
the software. Feel welcome to get in touch if you have any questions, want to provide feed back, have a feature request,
are interested in collaborating or just want to discuss this project. You can contact Nick via n.gorman at unsw.edu.au.      

## Installation
Installing nempy to use in your project is easy.

`pip install nempy`

To install for development purposes, such as adding new features. Download the source code, unzip, cd into the directory, then install.

`pip install e .[dev]`

Then the test suite can be run using.

`python -m pytest`

## A simple example
```python
import pandas as pd
from nempy import markets

# Volume of each bid, number of bands must equal number of bands in price_bids.
volume_bids = pd.DataFrame({
    'unit': ['A', 'B'],
    '1': [20.0, 50.0],  # MW
    '2': [20.0, 30.0],  # MW
    '3': [5.0, 10.0]  # More bid bands could be added.
})

# Price of each bid, bids must be monotonically increasing.
price_bids = pd.DataFrame({
    'unit': ['A', 'B'],
    '1': [50.0, 50.0],  # $/MW
    '2': [60.0, 55.0],  # $/MW
    '3': [100.0, 80.0]  # . . .
})

# Other unit properties
unit_info = pd.DataFrame({
    'unit': ['A', 'B'],
    'region': ['NSW', 'NSW'],  # MW
})

# The demand in the region\s being dispatched
demand = pd.DataFrame({
    'region': ['NSW'],
    'demand': [120.0]  # MW
})

# Create the market model
market = markets.SpotMarket(unit_info=unit_info, 
                            market_regions=['NSW'])
market.set_unit_volume_bids(volume_bids)
market.set_unit_price_bids(price_bids)
market.set_demand_constraints(demand)

# Calculate dispatch and pricing
market.dispatch()

# Return the total dispatch of each unit in MW.
print(market.get_unit_dispatch())
#   unit service  dispatch
# 0    A  energy      40.0
# 1    B  energy      80.0

# Return the price of energy in each region.
print(market.get_energy_prices())
#   region  price
# 0    NSW   60.0
```

## A detailed example
The example demonstrates the broad range of market features that can be implemented with nempy and the use of auxiliary 
modelling tools for accessing historical market data published by AEMO and preprocessing it for compatibility with nempy.
    
Warning: this example downloads approximately 8.5 GB of data from AEMO.
```python
# Notice: this script downloads large volumes of historical market data from AEMO's nemweb portal.

import sqlite3
import pandas as pd
from nempy import markets
from nempy.historical_inputs import loaders, mms_db, \
    xml_cache, units, demand, interconnectors, \
    constraints

con = sqlite3.connect('market_management_system.db')
mms_db_manager = mms_db.DBManager(connection=con)

xml_cache_manager = xml_cache.XMLCacheManager('cache_directory')

# The second time this example is run on a machine this flag can
# be set to false to save downloading the data again.
download_inputs = True

if download_inputs:
    # This requires approximately 5 GB of storage.
    mms_db_manager.populate(start_year=2019, start_month=1,
                            end_year=2019, end_month=1)

    # This requires approximately 3.5 GB of storage.
    xml_cache_manager.populate_by_day(start_year=2019, start_month=1, start_day=1,
                                      end_year=2019, end_month=1, end_day=1)

raw_inputs_loader = loaders.RawInputsLoader(
    nemde_xml_cache_manager=xml_cache_manager,
    market_management_system_database=mms_db_manager)

# A list of intervals we want to recreate historical dispatch for.
dispatch_intervals = ['2019/01/01 12:00:00',
                      '2019/01/01 12:05:00',
                      '2019/01/01 12:10:00',
                      '2019/01/01 12:15:00',
                      '2019/01/01 12:20:00',
                      '2019/01/01 12:25:00',
                      '2019/01/01 12:30:00']

# List for saving outputs to.
outputs = []

# Create and dispatch the spot market for each dispatch interval.
for interval in dispatch_intervals:
    raw_inputs_loader.set_interval(interval)
    unit_inputs = units.UnitData(raw_inputs_loader)
    interconnector_inputs = interconnectors.InterconnectorData(raw_inputs_loader)
    constraint_inputs = constraints.ConstraintData(raw_inputs_loader)
    demand_inputs = demand.DemandData(raw_inputs_loader)

    unit_info = unit_inputs.get_unit_info()
    market = markets.SpotMarket(market_regions=['QLD1', 'NSW1', 'VIC1',
                                                'SA1', 'TAS1'],
                                unit_info=unit_info)

    # Set bids
    volume_bids, price_bids = unit_inputs.get_processed_bids()
    market.set_unit_volume_bids(volume_bids)
    market.set_unit_price_bids(price_bids)

    # Set bid in capacity limits
    unit_bid_limit = unit_inputs.get_unit_bid_availability()
    market.set_unit_bid_capacity_constraints(unit_bid_limit)
    cost = constraint_inputs.get_constraint_violation_prices()['unit_capacity']
    market.make_constraints_elastic('unit_bid_capacity', violation_cost=cost)

    # Set limits provided by the unconstrained intermittent generation
    # forecasts. Primarily for wind and solar.
    unit_uigf_limit = unit_inputs.get_unit_uigf_limits()
    market.set_unconstrained_intermitent_generation_forecast_constraint(
        unit_uigf_limit)
    cost = constraint_inputs.get_constraint_violation_prices()['uigf']
    market.make_constraints_elastic('uigf_capacity', violation_cost=cost)

    # Set unit ramp rates.
    ramp_rates = unit_inputs.get_ramp_rates_used_for_energy_dispatch()
    market.set_unit_ramp_up_constraints(
        ramp_rates.loc[:, ['unit', 'initial_output', 'ramp_up_rate']])
    market.set_unit_ramp_down_constraints(
        ramp_rates.loc[:, ['unit', 'initial_output', 'ramp_down_rate']])
    cost = constraint_inputs.get_constraint_violation_prices()['ramp_rate']
    market.make_constraints_elastic('ramp_up', violation_cost=cost)
    market.make_constraints_elastic('ramp_down', violation_cost=cost)

    # Set unit FCAS trapezium constraints.
    unit_inputs.add_fcas_trapezium_constraints()
    cost = constraint_inputs.get_constraint_violation_prices()['fcas_max_avail']
    fcas_availability = unit_inputs.get_fcas_max_availability()
    market.set_fcas_max_availability(fcas_availability)
    market.make_constraints_elastic('fcas_max_availability', cost)
    cost = constraint_inputs.get_constraint_violation_prices()['fcas_profile']
    regulation_trapeziums = unit_inputs.get_fcas_regulation_trapeziums()
    market.set_energy_and_regulation_capacity_constraints(regulation_trapeziums)
    market.make_constraints_elastic('energy_and_regulation_capacity', cost)
    scada_ramp_down_rates = unit_inputs.get_scada_ramp_down_rates_of_lower_reg_units()
    market.set_joint_ramping_constraints_lower_reg(scada_ramp_down_rates)
    market.make_constraints_elastic('joint_ramping_lower_reg', cost)
    scada_ramp_up_rates = unit_inputs.get_scada_ramp_up_rates_of_raise_reg_units()
    market.set_joint_ramping_constraints_raise_reg(scada_ramp_up_rates)
    market.make_constraints_elastic('joint_ramping_raise_reg', cost)
    contingency_trapeziums = unit_inputs.get_contingency_services()
    market.set_joint_capacity_constraints(contingency_trapeziums)
    market.make_constraints_elastic('joint_capacity', cost)

    # Set interconnector definitions, limits and loss models.
    interconnectors_definitions = \
        interconnector_inputs.get_interconnector_definitions()
    loss_functions, interpolation_break_points = \
        interconnector_inputs.get_interconnector_loss_model()
    market.set_interconnectors(interconnectors_definitions)
    market.set_interconnector_losses(loss_functions,
                                     interpolation_break_points)

    # Add generic constraints and FCAS market constraints.
    fcas_requirements = constraint_inputs.get_fcas_requirements()
    market.set_fcas_requirements_constraints(fcas_requirements)
    violation_costs = constraint_inputs.get_violation_costs()
    market.make_constraints_elastic('fcas', violation_cost=violation_costs)
    generic_rhs = constraint_inputs.get_rhs_and_type_excluding_regional_fcas_constraints()
    market.set_generic_constraints(generic_rhs)
    market.make_constraints_elastic('generic', violation_cost=violation_costs)
    unit_generic_lhs = constraint_inputs.get_unit_lhs()
    market.link_units_to_generic_constraints(unit_generic_lhs)
    interconnector_generic_lhs = constraint_inputs.get_interconnector_lhs()
    market.link_interconnectors_to_generic_constraints(
        interconnector_generic_lhs)

    # Set the operational demand to be met by dispatch.
    regional_demand = demand_inputs.get_operational_demand()
    market.set_demand_constraints(regional_demand)
    
    # Get unit dispatch without fast start constraints and use it to
    # make fast start unit commitment decisions.
    market.dispatch()
    dispatch = market.get_unit_dispatch()
    fast_start_profiles = unit_inputs.get_fast_start_profiles_for_dispatch(dispatch)
    market.set_fast_start_constraints(fast_start_profiles)
    if 'fast_start' in market.get_constraint_set_names():
        cost = constraint_inputs.get_constraint_violation_prices()['fast_start']
        market.make_constraints_elastic('fast_start', violation_cost=cost)

    # If AEMO historical used the over constrained dispatch rerun
    # process then allow it to be used in dispatch. This is needed
    # because sometimes the conditions for over constrained dispatch
    # are present but the rerun process isn't used.
    if constraint_inputs.is_over_constrained_dispatch_rerun():
        market.dispatch(allow_over_constrained_dispatch_re_run=True,
                        energy_market_floor_price=-1000.0,
                        energy_market_ceiling_price=14500.0,
                        fcas_market_ceiling_price=1000.0)
    else:
        # The market price ceiling and floor are not needed here
        # because they are only used for the over constrained
        # dispatch rerun process.
        market.dispatch(allow_over_constrained_dispatch_re_run=False)

    # Save prices from this interval
    prices = market.get_energy_prices()
    prices['time'] = interval
    outputs.append(prices.loc[:, ['time', 'region', 'price']])

con.close()
print(pd.concat(outputs))
#                   time region      price
# 0  2019/01/01 12:00:00   NSW1  91.870167
# 1  2019/01/01 12:00:00   QLD1  76.190796
# 2  2019/01/01 12:00:00    SA1  86.899534
# 3  2019/01/01 12:00:00   TAS1  89.805037
# 4  2019/01/01 12:00:00   VIC1  84.984255
# 0  2019/01/01 12:05:00   NSW1  91.870496
# 1  2019/01/01 12:05:00   QLD1  64.991736
# 2  2019/01/01 12:05:00    SA1  87.462599
# 3  2019/01/01 12:05:00   TAS1  90.178036
# 4  2019/01/01 12:05:00   VIC1  85.556009
# 0  2019/01/01 12:10:00   NSW1  91.870496
# 1  2019/01/01 12:10:00   QLD1  64.991736
# 2  2019/01/01 12:10:00    SA1  86.868556
# 3  2019/01/01 12:10:00   TAS1  89.983716
# 4  2019/01/01 12:10:00   VIC1  84.936150
# 0  2019/01/01 12:15:00   NSW1  91.870496
# 1  2019/01/01 12:15:00   QLD1  64.776456
# 2  2019/01/01 12:15:00    SA1  86.844540
# 3  2019/01/01 12:15:00   TAS1  89.582288
# 4  2019/01/01 12:15:00   VIC1  84.990796
# 0  2019/01/01 12:20:00   NSW1  91.870496
# 1  2019/01/01 12:20:00   QLD1  64.776456
# 2  2019/01/01 12:20:00    SA1  87.496112
# 3  2019/01/01 12:20:00   TAS1  90.291144
# 4  2019/01/01 12:20:00   VIC1  85.594840
# 0  2019/01/01 12:25:00   NSW1  91.870167
# 1  2019/01/01 12:25:00   QLD1  64.991736
# 2  2019/01/01 12:25:00    SA1  87.519993
# 3  2019/01/01 12:25:00   TAS1  90.488064
# 4  2019/01/01 12:25:00   VIC1  85.630617
# 0  2019/01/01 12:30:00   NSW1  91.870496
# 1  2019/01/01 12:30:00   QLD1  64.991736
# 2  2019/01/01 12:30:00    SA1  87.462000
# 3  2019/01/01 12:30:00   TAS1  90.196284
# 4  2019/01/01 12:30:00   VIC1  85.573321
```
---
title: 'Nempy: A Python package for modelling the Australian National Electricity Market dispatch procedure'

tags:
  - Python
  - electricity markets
  - economic dispatch
  - Australian National Electricity Market
  - NEM
  - dispatch
authors:
  - name: Nicholas Gorman
    affiliation: "1, 3"
  - name: Anna Bruce
    affiliation: "1, 3"
  - name: Iain MacGill
    affiliation: "2, 3"
affiliations:
 - name: School of Photovoltaics and Renewable Energy Engineering, University of New South Wales, Australia
   index: 1
 - name: School of Electrical Engineering and Telecommunications, University of New South Wales, Australia
   index: 2
 - name: Collaboration on Energy and Environmental Markets (CEEM), University of New South Wales, Australia
   index: 3
date: 16 August 2021
bibliography: paper.bib
---

# Summary

Nempy is a python package for modelling the dispatch procedure of the Australian National Electricity Market (NEM).
Electricity markets are a way of coordinating the supply of electricity by private firms. The NEM is a gross pool spot 
market that operates on 5 min dispatch basis [@nemfactsheet]. Described simply, this means all generators wishing to sell electricity 
must bid into the market every 5 minutes, market clearing proceeds by calculating the cheapest combination of generator 
operating levels to meet forecast demand at the end of the 5 minute dispatch interval. The price of electricity is set as the 
marginal cost of generation, which, under a simple market formulation, would be the cost of the next generation bid to be 
dispatched if demand for electricity were to increase. Real-world formulations require significant additional adjustment 
in order to manage the technical complexity of securely and reliably operating an electricity grid. For example, in the 
case of the NEM additional markets for ancillary services have been introduced. One set of ancillary markets that have 
been integrated into the market dispatch procedure are the Frequency Control Ancillary Services (FCAS) markets. In these 
markets generators compete to provide the ability to rapidly change generation levels in order to control the grid frequency. 
Nempy is flexible in that it allows for the formulation of very simple market models, for the formulation of market models 
of near real-world complexity, and at the various levels of intermediate complexity. Simple models can be constructed 
using just generator bids and electricity demand, so called bid stack models. More complete models can be constructed by 
using the inbuilt features to create multiple market regions, ramp rate limits, loss factors, FCAS markets, FCAS trapezium 
constraints, dynamic interconnector loss models, generic constraints and fast start dispatch inflexibility profiles. 
Outputs include market clearing prices, generator and scheduled load dispatch targets, FCAS enablement levels, unit FCAS 
availability levels, interconnector flows, interconnector losses and region net inflows. Nempy is written in Python 3, 
and uses a relatively small number of first-order dependencies; pandas [@reback2020pandas; @mckinney-proc-scipy-2010], 
Numpy [@harris2020array], MIP-Python [@coin-orpython-mip], xmltodict [@xmltodict], and Requests [@psf].

# Statement of need

In modern industrialised economies, the electricity sector plays a key role in societal welfare and progress, yet 
commonly is also associated with major environmental harms, particularly where primary energy is sourced mainly through 
the burning of fossil fuels. As such, all of us are stakeholders in the continued successful operation of the 
electricity industry, while it transitions to cleaner energy sources, and beyond. Computer models are often used to 
study the operation, interactions and potential future direction of the electricity sector, review papers highlight the 
large body of work in this space [@ringkjob2018review; @chang2021trends; @fattahi2020systemic]. Such tools are, 
invariably, simplifications of the underlying processes of electricity industry operation and investment for reasons 
including the underlying complexity of the processes and the difficulty of gathering representative data. Commonly they 
tackle only a subset of the decision making that must operate from milliseconds (for example, under frequency relay 
trips) through to decades (investment in large generation units with long lead times). A particularly challenging and 
key task is that of operational dispatch – setting generator outputs and controllable network elements to meet expected 
demand over the next five to thirty minutes, and minimising costs while ensuring secure and reliable operation. To the 
best of the author's knowledge `Nempy` is the only open-source software that provides a detailed model of the NEM's 
dispatch procedure. Other more generalised models of the NEM [@grozev2005nemsim; @mcconnell2013retrospective; 
@ANEMWorkingReport; @mountain; @wood] or commercial tools such as PLEXOS [@energy_exemplar_plexos_2021] and Prophet 
[@IES] are used to model NEM dispatch, at various levels of complexity, but are not open-source. More recent work by 
Xenophon and Hill provides open-source code and data for modelling the NEM, but the dispatch functionality does not 
include many of the NEM wholesale market features [@xenophon2018open].

Nempy has been designed as a flexible model of the NEM's dispatch procedure and to be re-usable in a number of 
contexts. The software is aimed at analysts and modellers studying the NEM, either in industry or academia. It can be 
used as is, or as a building block in a larger modelling tool. Some potential use cases are outlined below:

1. As a tool for studying the dispatch process itself. Models of any energy system or electricity market are necessarily 
simplifications, however, to improve model performance it is often desirable to add additional detail. Nempy can be used 
to study the impact of different simplifications on modelling outcomes, and thus provide guidance on how model 
performance could be improved by adding additional detail. Figure 1 shows a simple example of such an analysis. The price
results from the New South Wales region for 1000 randomly selected intervals in the 2019 calendar year are shown. When
Nempy is configured with a full set of market features price results closely match actual prices. When the FCAS 
markets and generic constraints (network and security) are removed from the model, results differ significantly. Resorting
the results of the simpler market model, we can see that both models produce a similar number of medianly priced 
intervals. However, the highest and lowest priced intervals of the simpler model are significantly lower. The average
historical price is 81.4 $/MWh, the average price of the full featured model is 81.3 $/MWh, and the average price of the 
simpler model is 75.0 $/MWh. The close match between the results of the full featured model and actual prices allows 
for the attribution of the deviation of the simpler model explicitly to the simplifications that have been made.  

![Dispatch price results from the New South Wales region for 1000 randomly selected intervals in the 2019 calendar year.
  The actual prices, prior to scaling or capping, are also shown for comparison. Results from two Nempy models are
  shown, one with a full set of dispatch features, and one without FCAS markets or generic constraints (network and 
  security constraints). Actual prices, results from the full featured model, and the simpler model are shown in 
  descending order for actual prices, results from the simpler model are also shown resorted.\label{fig:example}](plot.png)

2. As a building block in agent based market models. Agent based models can be used to study electricity market 
operation, and are particularly useful in modelling both the competitive nature of electricity markets and their complex 
operational constraints [@ventosa]. In such models, agents must interact with a modelled environment, and a key part of that 
environment is the market dispatch process. Thus, Nempy could be useful as a building block to create agent based models 
of the NEM, and play a role in answering various questions about market operational outcomes. Such questions could 
include: 

    * How does changing the demand for electricity affect market outcomes? 
    * How does the entry of new generating technologies affect market outcomes? 
    * How do patterns of generator ownership affect market outcomes? 

    Of course, another necessary component of agent based models are the behavioural models of the agents, a prototype 
    behavioural model of NEM participants is being developed as part of the NEMPRO project [@nempro].

3. To answer counter factual questions about historical dispatch outcomes. For example:

    * What would have been the impact on market dispatch if a particular network constraint had not been present? 
    * How would have dispatch outcomes differed if a unit had offered a different bid into the market? 

    The answers to such questions have direct, and  potentially large, financial implications for market participants. 
    AEMO offers access to a production version of the market dispatch engine to allow participants to answer such questions 
    [@nemde]. However, access is restricted to registered participants and is provided at a cost of $15,000 per year. 
    Additionally, users of this service are not provided with a copy of the dispatch engine, but access it by submitting 
    input files to AEMO. This prevents the use of this service to answer questions about how changes to the dispatch 
    process, rather than the inputs, would affect dispatch outcomes. In contrast, access to Nempy is not restricted, it is 
    free to use, and is open to modification.

4. As a reference implementation of the NEM's dispatch procedure. While the Australian Energy Market Operator (AEMO) 
has published several documents that describe aspects of the dispatch process [@fcasmodel; @faststart; @lossfactors; 
@constraintviolation; @treatmentlossfactors], our experience developing Nempy has indicated that key 
implementation details are often missing from the publicly available documentation. Through a process of testing various 
implementation options, where the documentation was not explicit, Nempy has been refined in an attempt to better reflect 
the actual dispatch procedure. As a result, Nempy is a useful additional reference for analysts and modellers 
looking to understand the NEM's dispatch procedure.

# ReferencesPublications
============
Links to publications and associate source code.

Nempy Technical Brief
---------------------
The nempy technical brief is available :download:`here  <../../docs/pdfs/Nempy Technical Brief v1.0.0.pdf>`.
as pdf, and is also used as the introduction for readthedocs page.

Source code for Figure 1
************************

.. literalinclude:: ../../publications/all_features_example.py
    :linenos:
    :language: python

Source code for Figure 2
************************

.. literalinclude:: ../../publications/energy_only_market.py
    :linenos:
    :language: python
.. _spotmarket:

markets module
===============================
A model of the NEM spot market dispatch process.

Overview
--------
The market, both in real life and in this model, is implemented as a linear program. Linear programs consist of three
elements:

1.  **Decision variables**: the quantities being optimised for. In an electricity market these will be things like the
    outputs of generators, the consumption of dispatchable loads and interconnector flows.
2.  An **objective function**: the linear function being optimised. In this model of the spot market the cost of production
    is being minimised, and is defined as the sum of each bids dispatch level multiplied by the bid price.
3.  A set of **linear constraints**: used to implement market features such as network constraints and interconnectors.

The class :class:`nempy.SpotMarket` is used to construct these elements and then solve the linear program to calculate
dispatch and pricing. The examples below give an overview of how method calls build the linear program.

*   Initialising the market instance, doesn't create any part of the linear program, just saves general information for
    later use.

.. code-block:: python

    market = markets.SpotMarket(unit_info=unit_info, market_regions=['NSW'])

*   Providing volume bids creates a set of n decision variables, where n is the number of bids with a volume greater
    than zero.

.. code-block:: python

    market.set_unit_volume_bids(volume_bids)

*   Providing price bids creates the objective function, i.e. units will be dispatch to minimise cost, as determined
    by the bid prices.

.. code-block:: python

    market.set_unit_price_bids(price_bids)

*   Providing unit capacities creates a constraint for each unit that caps its total dispatch at a set capacity

.. code-block:: python

    market.set_unit_bid_capacity_constraints(unit_limits)

*   Providing regional energy demand creates a constraint for each region that forces supply from units and
    interconnectors to equal demand

.. code-block:: python

    market.set_demand_constraints(demand)

Specific examples for using this class are provided on the `examples1`_ page, detailed documentation of the class
:class:`nempy.markets.SpotMarket` is provided in the `Reference`_ material below.

.. _reference:

Reference
---------------------

.. automodule:: nempy.markets
    :autosummary:
    :members:




.. _examples1:

Examples
====================
A number of examples of how to use Nempy are provided below. Examples 1 to 5 are simple and aim introduce various
market features that can be modelled with Nempy in an easy to understand way, the dispatch and pricing outcomes are
explained in inline comments where the results are printed. Examples 6 and 7 show how to use the historical data input
preparation tools provided with Nempy to recreate historical dispatch intervals. Historical dispatch and pricing
outcomes can be difficult to interpret as they are usually the result of complex interactions between the many features
of the dispatch process, for these example the results are plotted in comparison to historical price outcomes.
Example 8 demonstrates how the outputs of one dispatch interval can be used as the initial conditions of the
next dispatch interval to create a time sequential model, additionally the current limitations with the approach are
briefly discussed.

1. Bid stack equivalent market
---------------------------
This example implements a one region bid stack model of an electricity market. Under the bid stack model, generators are
dispatched according to their bid prices, from cheapest to most expensive, until all demand is satisfied. No loss factors,
ramping constraints or other factors are considered.

.. literalinclude:: ../../examples/bidstack.py
    :linenos:
    :language: python


2. Unit loss factors, capacities and ramp rates
-----------------------------------------------
A simple example with two units in a one region market, units are given loss factors, capacity values and ramp rates.
The effects of loss factors on dispatch and market prices are explained.

.. literalinclude:: ../../examples/ramp_rates_and_loss_factors.py
    :linenos:
    :language: python


3. Interconnector with losses
-----------------------------
A simple example demonstrating how to implement a two region market with an interconnector. The interconnector is
modelled simply, with a fixed percentage of losses. To make the interconnector flow and loss calculation easy to
understand a single unit is modelled in the NSW region, NSW demand is set zero, and VIC region demand is set to 90 MW,
thus all the power to meet VIC demand must flow across the interconnetcor.

.. literalinclude:: ../../examples/interconnector_constant_loss_percentage.py
    :linenos:
    :language: python


4. Dynamic non-linear interconnector losses
----------------------------------------
This example demonstrates how to model regional demand dependant interconnector loss functions as decribed in the AEMO
:download:`Marginal Loss Factors documentation section 3 to 5  <../../docs/pdfs/Marginal Loss Factors for the 2020-21 Financial year.pdf>`.
To make the interconnector flow and loss calculation easy to understand a single unit is modelled in the NSW region,
NSW demand is set zero, and VIC region demand is set to 800 MW, thus all the power to meet VIC demand must flow across
the interconnetcor.


.. literalinclude:: ../../examples/interconnector_dynamic_losses.py
    :linenos:
    :language: python


5. Simple FCAS markets
----------------------------------------
This example implements a market for energy, regulation raise and contingency 6 sec raise, with
co-optimisation constraints as described in section 6.2 and 6.3 of
:download:`FCAS Model in NEMDE <../../docs/pdfs/FCAS Model in NEMDE.pdf>`.

.. literalinclude:: ../../examples/simple_FCAS_markets.py
    :linenos:
    :language: python


6. Simple recreation of historical dispatch
----------------------------------------
Demonstrates using Nempy to recreate historical dispatch intervals by implementing a simple energy market with unit bids,
unit maximum capacity constraints and interconnector models, all sourced from historical data published by AEMO.

.. image:: ../../examples/charts/energy_market_only_qld_prices.png
  :width: 600

*Results from example: for the QLD region a reasonable fit between modelled prices and historical prices is obtained.*

.. warning:: Warning this script downloads approximately 8.5 GB of data from AEMO. The download_inputs flag can be set
             to false to stop the script re-downloading data for subsequent runs.

.. note:: This example also requires plotly >= 5.3.1, < 6.0.0 and kaleido == 0.2.1. Run pip install plotly==5.3.1 and pip
          install kaleido==0.2.1

.. literalinclude:: ../../examples/recreating_historical_dispatch.py
    :linenos:
    :language: python

7. Detailed recreation of historical dispatch
------------------------------------------
This example demonstrates using Nempy to recreate historical dispatch intervals by implementing a  energy market using all the
features of the Nempy market model, all inputs sourced from historical data published by AEMO. Note each interval is
dispatched as a standalone simulation and the results from one dispatch interval are not carried over to be the initial
conditions of the next interval, rather the historical initial conditions are always used.

.. image:: ../../examples/charts/full_featured_market_qld_prices.png
  :width: 600

*Results from example: for the QLD region a very close fit between modelled prices and historical prices is obtained.*

.. warning:: Warning this script downloads approximately 8.5 GB of data from AEMO. The download_inputs flag can be set
             to false to stop the script re-downloading data for subsequent runs.

.. note:: This example also requires plotly >= 5.3.1, < 6.0.0 and kaleido == 0.2.1. Run pip install plotly==5.3.1 and pip
          install kaleido==0.2.1

.. literalinclude:: ../../examples/all_features_example.py
    :linenos:
    :language: python

8. Time sequential recreation of historical dispatch
-------------------------------------------------
This example demonstrates using Nempy to recreate historical dispatch in a dynamic or time sequential manner, this means the outputs
of one interval become the initial conditions for the next dispatch interval. Note, currently there is not the infrastructure
in place to include features such as generic constraints in the time sequential model as the rhs values of many constraints
would need to be re-calculated based on the dynamic system state. Similarly, using historical bids in this example is
some what problematic as participants also dynamically change their bids based on market conditions. However, for the sake
of demonstrating how Nempy can be used to create time sequential models, historical bids are used in this example.

.. warning:: Warning this script downloads approximately 8.5 GB of data from AEMO. The download_inputs flag can be set
             to false to stop the script re-downloading data for subsequent runs.

.. note:: This example also requires plotly >= 5.3.1, < 6.0.0 and kaleido == 0.2.1. Run pip install plotly==5.3.1 and pip
          install kaleido==0.2.1

.. literalinclude:: ../../examples/time_sequential.py
    :linenos:
    :language: python

8. Nempy performance on recent data (Jan 2022)
----------------------------------------------
This example demonstrates using Nempy to recreate historical dispatch intervals by implementing a energy market using all the
features of the Nempy market model, all inputs sourced from historical data published by AEMO. A set of 100 random dispatch
intervals from a recent month are dispatched and compared to historical results to see if Nempy is keeping up with any
recent changes to the NEM's dispatch procedure. Comparison is against ROP, the region price prior to any post dispatch
adjustments, scaling, capping etc.

Summary of results:

| Mean price error: -0.255
| Median price error: 0.00
| 5% percentile price error: -0.050
| 95% percentile price error: 1.084

.. warning:: Warning this script downloads approximately 84 GB of data from AEMO. The download_inputs flag can be set
             to false to stop the script re-downloading data for subsequent runs.

.. literalinclude:: ../../examples/recent_performance.py
    :linenos:
    :language: python


8. Nempy performance on older data (Jan 2015)
---------------------------------------------
This example demonstrates using Nempy to recreate historical dispatch intervals by implementing a energy market using all the
features of the Nempy market model, all inputs sourced from historical data published by AEMO. A set of 100 random dispatch
intervals from January 2015 are dispatched and compared to historical results to see how well Nempy performs for
replicating older versions of the NEM's dispatch procedure. Comparison is against ROP, the region price prior to any post
dispatch adjustments, scaling, capping etc.

Summary of results:

| Mean price error: -0.240
| Median price error: 0.000
| 5% percentile price error: 0.000
| 95% percentile price error: 0.051

.. warning:: Warning this script downloads approximately 54 GB of data from AEMO. The download_inputs flag can be set
             to false to stop the script re-downloading data for subsequent runs.

.. literalinclude:: ../../examples/performance_on_older_data.py
    :linenos:
    :language: python


.. _historical:

historical_inputs modules
===============================
The module provides tools for accessing historical market data and preprocessing for compatibility with the SpotMarket
class.

xml_cache
---------------

.. automodule:: nempy.historical_inputs.xml_cache
    :autosummary:
    :members:
    :inherited-members:
    :exclude-members: timedelta,datetime,time,Path,xmltodict

mms_db
-----------------------------------------

.. automodule:: nempy.historical_inputs.mms_db
    :autosummary:
    :members:
    :inherited-members:
    :exclude-members: timedelta,datetime

loaders
-----------------------------------------

.. automodule:: nempy.historical_inputs.loaders
    :autosummary:
    :members:

units
--------

.. automodule:: nempy.historical_inputs.units
    :autosummary:
    :members:

interconnectors
------------------

.. automodule:: nempy.historical_inputs.interconnectors
    :autosummary:
    :members:

demand
------------------

.. automodule:: nempy.historical_inputs.demand
    :autosummary:
    :members:

constraints
------------------

.. automodule:: nempy.historical_inputs.constraints
    :autosummary:
    :members:


.. _time_sequential:

time_sequential modules
=======================
The module provides tools constructing time sequential models using nempy.

.. automodule:: nempy.time_sequential
    :autosummary:
    :members:.. nempy documentation master file, created by
   sphinx-quickstart on Tue Apr 14 21:20:52 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to nempy's documentation!
=================================

.. toctree::
   :maxdepth: 4
   :caption: Contents:

   intro
   installation
   examples
   markets
   historical
   time_sequential
   publications

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
Installation
============
Installing nempy to use in your project is easy.

`pip install nempy`

To install for development purposes, such as adding new features. Download the source code, unzip, cd into the directory, then install.

`pip install e .[dev]`

Then the test suite can be run using.

`python -m pytest`Introduction
============
Nempy is an open-source python package that can be used to model the dispatch procedure of the Australian National
Electricity Market (NEM). The dispatch process is at the core of many market modelling projects. As the
NEM evolves, constraints and ancillary service markets are becoming increasingly important in shaping dispatch outcomes.
As part of the ongoing reform process significant changes to the dispatch process have also been proposed, for example,
the introduction of an operating reserve market is being considered, as is the adoption of a dynamic intra-regional loss
model. Nempy allows users to easily configure a dispatch model to fit the relevant research question. Furthermore, if
extra functionality is needed, the python implementation, open-source licencing and planned ongoing support from developers
make it possible to adapt Nempy to your needs. Version 1 will be a stable release and ongoing minor updates or patches
will remain backwards compatible. Currently the latest minor release is v1.1.0. However, we are interested in receiving
feedback to inform the ongoing maintenance, development and any future major updates.

Nempy is feature rich, flexible, can recreate historical dispatch with a high degree of accuracy, runs fast, has detailed
documentation and has planned support until mid-2023.

The Nempy source code is on GitHub: https://github.com/UNSW-CEEM/nempy.

A brief introduction to the NEM can be found here: https://aemo.com.au/-/media/Files/Electricity/NEM/National-Electricity-Market-Fact-Sheet.pdf

Example use cases
-----------------
Nempy is intended for analysts and modellers studying the NEM either in industry or academic. It can be
used either as is, or as building block in a large modelling tool. Some potential use case are:

#. As a tool for studying the dispatch process itself. The example shown in the:ref:`section on model accuracy <Accuracy>`
   below demonstrates how model simplifications effects accuracy, this is potentially useful information for other
   NEM modeller either using Nempy or other modelling tools.
#. As a building block in agent based market models, as part of the environment for agents to interact with.
#. To answer counter factual questions about historical dispatch outcomes. For example, how removing a network
   constraint would have effected dispatch and pricing outcomes?
#. As a reference implementation of the NEM's dispatch procedure. Published documentation can lack detail, studying the
   source code of Nempy may be useful for some NEM analysts to gain a better understanding of the dispatch procedure.


Dispatch Procedure Outline
--------------------------
The main task of the dispatch procedure is the construction and solving of a mixed integer linear problem (MIP) to find the
least cost set of dispatch levels for generators and scheduled loads. Note, in this optimisation the dispatch of
scheduled loads is treated as a negative cost, this makes the least cost optimisation equivalent to maximising the value of
market trade. The construction of the MIP as implemented by Nempy proceeds roughly as follows:

#. Bids from generators and loads are preprocessed, some FCAS bids are excluded if they do not meet a set of inclusion
   criteria set out by AEMO (:download:`FCAS Model in NEMDE <../../docs/pdfs/FCAS Model in NEMDE.pdf>`).
#. For each bid a decision variable in the MIP is created, the cost of the variable in the objective function is the bid
   price, and the price is adjusted by a loss factor if one is provided.
#. For each market region a constraint forcing generation to equal demand is created.
#. The rest of the market features are implemented as additional variables and/or constraints in the MIP, for example:

   - unit ramp rates are converted to a set MW ramp that units can achieve over the dispatch interval, and the sum of a
     unit's dispatch is limited by this MW value
   - interconnectors are formulated as additional decision variables that link the supply equals demand constraints
     of the interconnected regions, and are combined with constraints sets that enforce interconnector losses as a
     function of power flow

#. The MIP is solved to determined interconnector flows and dispatch targets, the MIP is then converted to a linear
   problem, and re-solved, such that market prices can be determined from constraint shadow prices.

Other steps in the dispatch procedure that are not implemented by Nempy are:

#. The preprocessing step to calculate of FCAS market and network constraint right hand side values (right hand side
   values need to be provided as inputs to Nempy)
#. Multiple re-runs of the optimisation to the operational settings for DC link between mainland synchronous region and
   the Tasmainian synchronous region


Features
--------
- **Energy bids**: between one and ten price quantity bid pairs can be provided for each generator or load bidding in the energy market
- **Loss factors**: loss factors can be provided for each generator and load
- **FCAS bids**: between one and ten price quantity bid pairs can be provided for each generator or load bidding in each of the eight FCAS markets
- **Ramp rates**: unit ramp rates can be set
- **FCAS trapezium constraints**: a set of trapezium constraints can be provided for each FCAS bid, these ensure FCAS is co-optimised with energy dispatch and would be technically deliverable
- **Fast start dispatch inflexibility profiles**: dispatch inflexibility profiles can be provided  for unit commitment of fast-start plants
- **Interconnectors and losses**: interconnectors between each market region can be defined, non-linear loss functions and interpolation breakpoints for their linearisation can be provided
- **Generic constraints**: generic constraints that link across unit output, FCAS enablement and interconnector flows can be defined
- **Elastic constraints**: constraints can be made elastic, i.e. a violation cost can be set for constraints
- **Tie-break constraints**: constraints that minimise the difference in dispatch between energy bids for the same price can be enabled
- **Market clearing prices**: market prices are returned for both energy and FCAS markets, based on market constraint shadow prices
- **Historical inputs**: tools for downloading dispatch inputs from AEMO's NEMWeb portal and preprocessing them for compatibility with the nempy SpotMarket class are available
- **Input validation**: optionally check user inputs and raise descriptive errors when they do not meet the expected criteria
- **Adjustable dispatch interval**: a dispatch interval of any length can be used

Flexibility
-----------
Nempy is designed to have a high degree of flexibility, it can be used to implement very simple merit order dispatch models,
highly detailed models that seek to re-create the real world dispatch procedure, or a model at the many levels of intermediate
complexity. A set of :ref:`examples, <examples1>` demonstrating this flexibility are available. Most inputs are passed to nempy as pandas DataFrame
objects, which means Nempy can easily source inputs from other python code, SQL databases, CSVs and other formats supported by
the pandas' interface.

Accuracy
--------
The accuracy with which Nempy represents the NEM's dispatch process can be measured by re-creating historical dispatch results.
This is done for a given dispatch interval by downloading the relevant historical inputs such as unit initial operating levels,
bids and generic constraints, processing these inputs so they are compatible with the Nempy SpotMarket class, and finally
dispatching the spot market. The results can then be compared to historical results to gauge the model's accuracy.
Figure 1 shows the results of this process for 1000 randomly selected dispatch intervals in 2019, comparing the modelled
NSW energy price with historical prices. Here the model is configured to maximally reflect the NEM's dispatch procedure.
The code to produce the results shown in this figure is available `here <https://nempy.readthedocs.io/en/latest/publications.html#source-code-for-figure-1>`_.
Figure 2 shows a similar comparison, but without FCAS markets or generic constraints. The code to produce the results
shown in Figure 2 is available `here <https://nempy.readthedocs.io/en/latest/publications.html#source-code-for-figure-2>`_.
The simpler model produces a similar number of medianly priced intervals, however, outcomes for extreme ends of the price
duration curve differ significantly from historical values.

.. image:: nempy_vs_historical.svg
  :width: 600

*Figure 1: A comparison of the historical NSW reference node price, prior to scaling or capping, with the price calculated using nempy.
The nempy model was configured to maximally replicated the NEM dispatch process and 1000 randomly selected intervals were used.*

.. image:: nempy_vs_historical_simple.svg
  :width: 600

*Figure 2: A comparison of the historical NSW reference node price, prior to scaling or capping, with the price calculated
using Nempy. The Nempy model was configured without FCAS markets or generic constraints and 1000 randomly selected intervals were used.*

Run-time
--------
The run-time for Nempy to calculate dispatch depends on several factors, the complexity of the model implemented, time
taken to load inputs, the mixed-integer linear solver used and of course the hardware. Run-times reported here used an
Intel® Xeon(R) W-2145 CPU @ 3.70 GHz. For the model results shown in Figure 1, including time taken to load inputs from
the disk and using the open-source solver CBC, the average run-time per dispatch interval was 2.54 s. When the proprietary
solver Gurobi was used, a run-time of 1.84 s was achieved. For the results shown in Figure 2, the run-times with CBC and
Gurobi were 1.02 s and 0.98 s respectively, indicating that for simpler models the solver used has a smaller impact on
run-time. For the simpler model, the time to load inputs is increased significantly by the loading of historical NEMDE
input/output XML files which takes approximately 0.4 s. Importantly, this means it will be possible to speed up simpler
models by sourcing inputs from different data storage formats.

Notes:

- Information on solvers is provided is provided in the `reference documentation <https://nempy.readthedocs.io/en/latest/markets.html#nempy.markets.SpotMarket.solver_name>`_
  of the SpotMarket class.
- The total runtime was calculated using the python time module and measuring the time taken from the loading of inputs
  to the extraction of results from the model. The runtime of different sub-process, i.e. loading of the XML file, was
  measured by inserting timing code into the Nempy source code where required.

Documentation
-------------
Nempy has a detailed set of documentation, mainly comprising of two types: examples and reference documentation. The
examples aim to show how Nempy can be used and how it works in a practical manner. A number of simple examples focus on
demonstrating the use of subsets of the package's features in isolation in order to make them easier to understand. The
more complex examples show how features can be combined to build models more suitable for analysis. The reference
documentation aims to cover all the package's public APIs (the classes, methods and functions accessible to the user),
describing their use, inputs, outputs and any side effects.

Support
-------
Nempy's development is being led by Nick Gorman as part of his PhD candidature at the Collaboration on Energy and Environmental
Markets at the University of New South Wales' School of Photovoltaics and Renewable Energy Engineering. As part of this
project we plan to engage with and support software users, this can be facilitated through the PhD until mid-2023. If
Nempy is used sufficiently broadly we would look to continue support beyond this timeframe.


Ongoing work
------------
Maintenance:

1. Retest Nempy on 2020 and 2021 historical data, previous testing has been against 2019 data.

Enhancements:

* No enhancements are currently planned for Nempy. However, development is active on a market participant behavioural
  modelling package that would strongly complement the functionality of Nempy, https://github.com/UNSW-CEEM/NEMPRO .

Dependencies
------------
* pandas >=1.0.0, <2.0.0
* mip>=1.11.0, <2.0.0: https://github.com/coin-or/python-mip)
* xmltodict==0.12.0:  https://github.com/martinblech/xmltodict)
* requests>=2.0.0, <3.0.0

