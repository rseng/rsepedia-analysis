<h1 align="center">epimargin</h1>

<div align="center"> <img src="./docs/logo.svg" height="150"> </div>

<hr>

a public health policy analysis toolkit consisting of: 
1. Bayesian simulated annealing estimator for the reproductive rate (<i>R<sub>t</sub></i>)
2. a stochastic compartmental model class supporting multiple compartment schemes (SIR, SEIR, SIRV, subpopulations, etc.)
3. policy impact evaluation that calculates longevity benefits and economic activity disruption/resumption under scenarios including: 
    - targeted lockdowns
    - urban-rural migration and reverse-migration flows
    - multiple vaccine allocation prioritizations 

For examples of how the package is used, see the docs folder (featuring a [toy example](https://github.com/COVID-IWG/epimargin/blob/master/docs/toy_example.py) and a [full tutorial](https://github.com/COVID-IWG/epimargin/blob/master/docs/tutorial.py)), or the [epimargin-studies](https://github.com/COVID-IWG/epimargin-studies) repository. 

# installation  
The `epimargin` package is available on PyPI and can be installed via `pip`: 
    
    pip install epimargin

# support and issues 
Please [file an issue on Github](https://github.com/COVID-IWG/epimargin/issues/new/choose) if you run into any problems with the software.

# contributions and development 
Contributions are always welcome! Please fork the repository and open a new pull request for any new features.

For development, we recommending installing the development dependencies, and then installing the package in editable mode: 

    git clone https://github.com/COVID-IWG/epimargin
    cd epimargin

    pip install -r requirements.txt
    pip install -e . 

We also recommend using a virtual environment for development.

# tutorial 
In this tutorial, we will download a timeseries of daily confirmed COVID-19 cases in Mumbai from COVID19India.org, estimate the reproductive rate for the city over time, plug these estimates into a compartmental model, and compare two policy scenarios by running the compartmental model forward. The entire tutorial can be found in the `docs/tutorial` directory.

## 1. setup
After installing the package, import some commonly-used tools and set up convenience functions/variables:

    from itertools import cycle

    import epimargin.plots as plt
    import numpy as np
    import pandas as pd
    from epimargin.utils import setup

    (data, figs) = setup() # creates convenient directories
    plt.set_theme("minimal")


## 2. download and clean data 
Next, download data on daily COVID19 cases in the city of Mumbai from COVID19India.org. The data are noisy, so we apply a notch filter to remove weekly reporting artifacts and smooth using a convolution:
    
    from epimargin.etl import download_data
    from epimargin.smoothing import notched_smoothing

    download_data(data, "districts.csv", "https://api.covid19india.org/csv/latest/") 

    # load raw data
    daily_reports = pd.read_csv(data / "districts.csv", parse_dates = ["Date"])\
        .rename(str.lower, axis = 1)\
        .set_index(["state", "district", "date"])\
        .sort_index()\
        .loc["Maharashtra", "Mumbai"]
    daily_cases = daily_reports["confirmed"]\
        .diff()\
        .clip(lower = 0)\
        .dropna()\

    # smooth/notch-filter timeseries
    smoother = notched_smoothing(window = 5)
    smoothed_cases = pd.Series(
        data  = smoother(daily_cases),
        index = daily_cases.index
    )

    # plot raw and cleaned data 
    beg = "December 15, 2020"
    end = "March 1, 2021"
    training_cases = smoothed_cases[beg:end]

    plt.scatter(daily_cases[beg:end].index, daily_cases[beg:end].values, color = "black", s = 5, alpha = 0.5, label = "raw case count data")
    plt.plot(training_cases.index, training_cases.values, color = "black", linewidth = 2, label = "notch-filtered, smoothed case count data")
    plt.PlotDevice()\
        .l_title("case timeseries for Mumbai")\
        .axis_labels(x = "date", y = "daily cases")\
        .legend()\
        .adjust(bottom = 0.15, left = 0.15)\
        .format_xaxis()\
        .size(9.5, 6)\
        .save(figs / "fig_1.svg")\
        .show()

![raw and cleaned case count timeseries data](./docs/fig_1.svg)

NOTE: a copy of the reference timeseries for all districts available through the API is checked into the `docs/data` folder in case you run into download issues or if the upstream API changes.

## 3. estimate the reproductive rate, <i>R<sub>t</sub></i>
From these data, we can estimate the reproductive rate, or the number of secondary infections caused by a single active infection. A pandemic is under control if the reproductive rate stays below 1. A number of estimation procedures are provided; we show the Bettencourt/Soman estimator as an example:

    from epimargin.estimators import analytical_MPVS

    (dates, Rt, Rt_CI_upper, Rt_CI_lower, *_) = analytical_MPVS(training_cases, smoother, infectious_period = 10, totals = False)
    plt.Rt(dates[1:], Rt[1:], Rt_CI_upper[1:], Rt_CI_lower[1:], 0.95, legend_loc = "upper left")\
        .l_title("$R_t$ over time for Mumbai")\
        .axis_labels(x = "date", y = "reproductive rate")\
        .adjust(bottom = 0.15, left = 0.15)\
        .size(9.5, 6)\
        .save(figs / "fig_2.svg")\
        .show()

![estimated reproductive rate over time](./docs/fig_2.svg)

## 4. set up a model and run it forward to compare policy scenarios

Finally, we can use the case count data and estimated reproductive rate to project forward cases. We also show how the input data can be modified to test hypotheses about specific policies. For example, you might expect a lockdown policy to reduce the reproductive rate by 25% given historical mobility data or lockdown stringency indices. Assuming a successful reduction in <i>R<sub>t</sub></i>, what does the trajectory of daily cases look like?

    from epimargin.models import SIR

    num_sims = 100
    N0 = 12.48e6
    R0, D0 = daily_reports.loc[end][["recovered", "deceased"]]
    I0  = smoothed_cases[:end].sum()
    dT0 = smoothed_cases[end]
    S0  = N0 - I0 - R0 - D0
    Rt0 = Rt[-1] * N0 / S0
    no_lockdown = SIR(
        name = "no lockdown", 
        population = N0, 
        dT0 = np.ones(num_sims) * dT0, Rt0 = np.ones(num_sims) * Rt0, I0 = np.ones(num_sims) * I0, R0 = np.ones(num_sims) * R0, D0 = np.ones(num_sims) * D0, S0 = np.ones(num_sims) * S0, infectious_period = 10
    )
    lockdown = SIR(
        name = "partial lockdown", 
        population = N0, 
        dT0 = np.ones(num_sims) * dT0, Rt0 = np.ones(num_sims) * 0.75 * Rt0, I0 = np.ones(num_sims) * I0, R0 = np.ones(num_sims) * R0, D0 = np.ones(num_sims) * D0, S0 = np.ones(num_sims) * S0, infectious_period = 10
    )

    # run models forward 
    simulation_range = 7
    for _ in range(simulation_range):
        lockdown   .parallel_forward_epi_step(num_sims = num_sims)
        no_lockdown.parallel_forward_epi_step(num_sims = num_sims)

    # compare policies 
    test_cases = smoothed_cases["February 15, 2021":pd.Timestamp(end) + pd.Timedelta(days = simulation_range)]
    date_range = pd.date_range(start = end, periods = simulation_range + 1, freq = "D")
    legend_entries = [plt.predictions(date_range, model, color) for (model, color) in zip([lockdown, no_lockdown], cycle(plt.SIM_PALETTE))]
    train_marker, = plt.plot(test_cases[:end].index, test_cases[:end].values, color = "black")
    test_marker,  = plt.plot(test_cases[end:].index, test_cases[end:].values, color = "black", linestyle = "dotted")
    markers, _ = zip(*legend_entries)
    plt.PlotDevice()\
        .l_title("projected case counts")\
        .axis_labels(x = "date", y = "daily cases")\
        .legend(
            [train_marker, test_marker] + list(markers),
            ["case counts (training)", "case counts (actual)", "case counts (partial lockdown; 95% simulation range)", "case counts (no lockdown; 95% simulation range)"],
            loc = "upper left"
        )\
        .adjust(bottom = 0.15, left = 0.15)\
        .size(9.5, 6)\
        .format_xaxis()\
        .save(figs / "fig_3.svg")\
        .show()

![policy projections over time](./docs/fig_3.svg)
The median projections from the no-lockdown model (crimson) mirror the observed timeseries (dotted black) fairly well, and the model predicts that even an imperfect lockdown would have changed the trajectory of the pandemic at the time period we looked at. As the model is stochastic, we show a range of outcomes (shaded) and note the model accuracy decreases as the projection period goes on. In real-time settings, we encourage daily updating of projections to handle new data.
---
title: '`epimargin`: A Toolkit for Epidemiological Estimation, Prediction, and Policy Evaluation'
tags:
  - Python
  - epidemiology
  - stochastic processes
  - economics
  - COVID-19
  - Bayesian inference
authors:
  - name: Satej Soman 
    orcid: 0000-0001-8450-7025
    affiliation: 1
  - name: Caitlin Loftus
    affiliation: 1
  - name: Steven Buschbach
    affiliation: 1
  - name: Manasi Phadnis
    affiliation: 1
  - name: Luís M. A. Bettencourt
    orcid: 0000-0001-6176-5160
    affiliation: "1, 2, 3, 4"
affiliations:
  - name: Mansueto Institute for Urban Innovation, University of Chicago
    index: 1
  - name: Department of Ecology & Evolution, University of Chicago
    index: 2
  - name: Department of Sociology, University of Chicago
    index: 3
  - name: Santa Fe Institute
    index: 4
date: 22 May 2021
bibliography: paper.bib

---

# Summary
As pandemics (including the COVID-19 crisis) pose threats to societies, public health officials, epidemiologists, and policymakers need improved tools to assess the impact of disease, as well as a framework for understanding the effects and tradeoffs of health policy decisions. The `epimargin` package provides functionality to answer these questions in a way that incorporates and quantifies irreducible uncertainty in both the input data and complex dynamics of disease propagation.  

The `epimargin` software package primarily consists of: 

1. a set of Bayesian estimation procedures for epidemiological metrics such as the reproductive rate ($R_t$), which is the average number of secondary infections caused by an active infection

2. a flexible, stochastic epidemiological model informed by estimated metrics and reflecting real-world epidemic and geographic structure, and 

3. a set of tools to evaluate different public health policy choices simulated by the model.

The software is implemented in the Python 3 programming language and is built using commonly-used elements of the Python data science ecosystem, including NumPy [@harris2020array], Scipy [@virtanen2020scipy], and Pandas [@mckinney2011pandas].

# Statement of need

The `epimargin` software package is designed for the data-driven analysis of policy choices related to the spread of disease. It consists primarily of a set of estimators for key epidemiological metrics, a stochastic model for predicting near-future disease dynamics, and evaluation tools for various policy scenarios. 

Included with the package are connectors and download utilities for common sources of disease data for the COVID-19 pandemic (the pressing concern at the time of writing), as well as a set of tools to prepare and clean data in a format amenable to analysis. It is widely understood that preprocessing epidemiological data is necessary to make inferences about disease progression [@gostic2020practical]. To that end, `epimargin` provides commonly-used preprocessing routines to encourage explicit documentation of data preparation, but is agnostic to which procedures are used due to the fact that all metadata required for certain preparations may not be uniformly available across geographies. 

This same modularity extends to both the estimation procedures and epidemiological models provided by `epimargin`. While the package includes a novel Bayesian estimator for key metrics, classical approaches based on rolling linear regressions and Markov chain Monte Carlo sampling are also included. The core model class in `epimargin` in which these estimates are used is known as a <i>compartmental</i> model: a modeled population is split into a number of mutually-exclusive compartments (uninfected, infected, recovered, vaccinated, etc) and flows between these compartments are estimated from empirical data. The exact choice of compartments and interactions is left to the modeler, but the package includes several commonly-used models, as well as variations customized for specific policy questions (such as large-scale migration during pandemics, or the effects of various vaccine distribution policies).

For similar data downloading tools, see `covidregionaldata` [@Palmer2021]; for similar estimation tools, see `EpiEstim` [@EpiEstim] and `EpiNow2` [@EpiNow]. While many of these tools are used in conjunction with each other, `epimargin` aims to offer tools for an end-to-end epidemiological workflow in one package, while offering the flexibility in estimator choice and data preparation methods.

Attempts to use a compartmental model to drive policy decisions often treat the systems under study as deterministic and vary parameters such as the reproductive rate across a range deemed appropriate by the study authors [@bubar2021model]. This methodology complicates incorporation of recent disease data and the development of theories for why the reproductive rate changes due to socioeconomic factors external to the model. The incorporation of stochasticity into the models from the outset allows for the quantification of uncertainty and the illustration of a range of outcomes for a given public health policy under consideration.

The `epimargin` package has been used to drive a number of research projects and inform policy decisions in a number of countries:

1. lockdown, quarantine planning, migrant return policies, and vaccine distribution in India and Indonesia (at the behest of national governments, regional authorities, and various NGOs)

2. an illustration of a novel Bayesian estimator for the reproductive rate as well as general architectural principles for real-time epidemiological systems [@bettencourt2020systems]

3. a trigger-based algorithmic policy for determining when administrative units of a country should exit or return to a pandemic lockdown based on projected reproductive rates and case counts [@malani2020adaptive]

4. a World Bank study of vaccination policies in South Asia [@southasiavaccinates]

5. a general framework for quantifying the health and economic benefits to guide vaccine prioritization and distribution [@vaccineallocation]


# Figures
Sample output for common workflows are illustrated in the following figures:

## downloaded and cleaned time series
![Raw and cleaned case count timeseries for Mumbai downloaded from COVID19India.org.\label{fig:fig1}](fig_1.png){ width=80% }

## estimated reproductive rate
![Estimated reproductive rate over time for Mumbai](fig_2.png){ width=80% }

## forward projection/policy comparison
![Projected case counts using a stochastic compartmental model and reproductive rate estimates](fig_3.png){ width=80% }

# Acknowledgements

We acknowledge code review and comments from Gyanendra Badgaiyan (IDFC Institute), ongoing conversations with Anup Malani (University of Chicago) and Nico Marchio (Mansueto Institute) and helpful discussions with Katelyn Gostic (University of Chicago) and Sarah Cobey (University of Chicago).

# References The `adaptive-control` project is a Python3 package, and can be installed by following these steps:

0. Install Python3 on your development machine, if you do not already have it.

1. Clone or download the `adaptive-control` repository.

2. Within the repository folder, create a virtual environment.
```
python3 -mvenv venv
```

3. Activate your virtual environment.
```
source venv/bin/activate
```

4. Install the requirements.
```
pip install -r requirements.txt
```

5. Install the project. 
```
pip install -e . 
```
