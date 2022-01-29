# PyEI

[![DOI](https://joss.theoj.org/papers/10.21105/joss.03397/status.svg)](https://doi.org/10.21105/joss.03397)

PyEI is a Python library for ecological inference. The target audience is the analyst with an interest in the phenomenon called Racially Polarized Voting.

Racially Polarized Voting is a legal concept developed through case law under the Voting Rights Act of 1965; its genesis is in the majority opinion of ***Thornburg v. Gingles (1982)***. Considered the “evidentiary linchpin” for vote dilution cases, RPV is a necessary, but not sufficient, condition that plaintiffs must satisfy for a valid claim. 

Toward that end, ecological inference uses observed data (historical election results), pairing voting outcomes with demographic information
for each precinct in a given polity, to infer voting patterns for each demographic group.

PyEI brings together a variety of ecological inference methods in one place and facilitates reporting and plotting results; quantifying the uncertainty associated with results under a given model; making comparisons between methods; and bringing relevant diagnostic tools to bear on ecological inference methods.

PyEI is relatively new and under active development, so expect rough edges and bugs -- and for additional features and documentation to be coming quickly!

## Want to use PyEI? Start here.

### Installation
You can install the latest release from `PyPi` with:

```
pip install pyei
```

Or, install directly from GitHub for the most up-to-date (but potentially less stable) version:

```
pip install git+git://github.com/mggg/ecological-inference.git
 ```
 
If you would like to explore PyEI without installation, you can explore this [interactive Colab notebook](https://colab.research.google.com/drive/1w4vWJBMEY_ULal9LWTOa_TrimPVWfum0#scrollTo=tLPaJ279zsG_) (just note that inference might be slow!)


### Example notebooks

Check out the [intro notebooks](https://github.com/mggg/ecological-inference/tree/main/pyei/intro_notebooks) and [example notebooks](https://github.com/mggg/ecological-inference/tree/main/pyei/examples) for sample code
that shows how to run and adjust the various models in PyEI on datesets.  

If you are new to ecological inference generally, start with `intro_notebooks/Introduction_toEI.ipynb`.

If you are familiar with ecological inference and want an overview of PyEI and how to use it, with examples start with `intro_notebooks/PyEI_overview.ipynb`.

To explore EI's plotting functionality, check out `intro_notebooks/Plotting_with_PyEI.ipynb`.

For more work with two-by-two examples, see in `examples/santa_clara_demo.ipynb`.

For more work with r-by-c examples, see `examples/santa_clara_demo_r_by_c.ipynb`.

For examples of depth model comparison and checking steps with PyEI, see `examples/model_eval_and_comparison_demo.ipynb`.

### Issues

Feel free to file an issue if you are running into trouble or if there is a feature you'd particularly like to see, and we will do our best to get to it!


## Want to contribute to PyEI? Start here.

Contributions are welcome! 

Uses Python>=3.7. After cloning the environment, you should be able to use either `virtualenv` or `conda` to run the code. The second (`conda`) is probably easier for development, but `virtualenv` is used for the project's CI.

Here is how to create and activate each environment. See the docs for more elaborate details:

### Install with virtualenv

```bash
virtualenv pyei_venv           # create virtualenv
source pyei_venv/bin/activate  # activate virtualenv
python -m pip install -U pip   # upgrade pip
python -m pip install -e .     # install project locally
python -m pip install -r requirements-dev.txt  # install dev requirements
```

### Install with conda

```bash
conda create --name pyei --channel conda-forge python=3.8 --file requirements.txt --file requirements-dev.txt # create conda environment and install requirements
pip install -e . #install project locally
```

### Testing

After making changes, make sure everything works by running

```bash
./scripts/lint_and_test.sh
```

This will also run automatically when you make a pull request, so if you have trouble getting that to run, just open the PR, and we can help!


## Citation

If you are using PyEI, please cite it as: 

Knudson et al., (2021). PyEI: A Python package for ecological inference. Journal of Open Source Software, 6(64), 3397, https://doi.org/10.21105/joss.03397

BibTeX:

```
@article{Knudson2021,
  doi = {10.21105/joss.03397},
  url = {https://doi.org/10.21105/joss.03397},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {64},
  pages = {3397},
  author = {Karin C. Knudson and Gabe Schoenbach and Amariah Becker},
  title = {PyEI: A Python package for ecological inference},
  journal = {Journal of Open Source Software}
}
```


---
title: 'PyEI: A Python package for ecological inference'
tags:
  - Python
  - voting rights
  - ecological inference
  - Bayesian inference
authors:
  - name: Karin C. Knudson
    orcid: 0000-0003-1312-6473
    affiliation: 1
  - name: Gabe Schoenbach
    orcid: 0000-0002-4300-0356
    affiliation: 2
  - name: Amariah Becker
    orcid: 0000-0002-7392-8307
    affiliation: 2
affiliations:
 - name: Data Intensive Studies Center, Tufts University
   index: 1
 - name: MGGG Redistricting Lab, Tufts University
   index: 2
date: 10 May 2021
bibliography: paper.bib
---


# Summary

An important question in some voting rights and redistricting litigation in the U.S. is whether and to what degree voting is racially polarized.
In the setting of voting rights cases, there is a family of methods called "ecological inference" [see especially @king1997solution] that uses
observed data, pairing voting outcomes with demographic information
for each precinct in a given polity, to infer voting patterns for each demographic group.

More generally, we can think of ecological inference as seeking to use knowledge about the margins of a set of tables (\autoref{fig:table_ex}) to infer associations between the row and column variables, by making (typically probabilistic) assumptions. In the context of assessing racially polarized voting, a table like the one in \autoref{fig:table_ex} will correspond to a precinct, where each column corresponds to a candidate or voting outcome and each row to a racial group. Ecological inference methods then use the vote counts and demographic data for each precinct to make inferences about the overall voting preferences by demographic group, thus addressing questions like: "What percentage of East Asian voters voted for Hardy?". This example is an instance of what is referred to in the literature as "R by C" ecological inference, where here we have R $=$ 2 groups and C $=$ 3 voting outcomes.
`PyEI` was created to support performing ecological inference with voting data; however, ecological inference methods also applicable in other fields, such as epidemiology [@elliot2000spatial] and sociology [@goodman1953ecological].

\newpage
|                 | Hardy         | Kolstad        | Nadeem |   |
| :-------------  | :-----------: | :------------: | :-------------: | :------------- |
| East Asian      | ?             |  ?                    | ? | Total East Asian |
| non- East Asian | ?     | ?              | ? | Total non- East Asian |
|                 | Total for    | Total for   | Total for |  |
|                 |  Hardy     | Kolstad   |  Nadeem |  |

Table: In ecological inference we have information about the marginal counts for a set of tables like the one here and would like to make inferences about, for example, the number or proportion of East Asian voters who voted for Hardy. The system is underdetermined and ecological inference methods proceed by making statistical assumptions. \label{fig:table_ex}

# Statement of need

The results of ecological inference for inferring racially polarized voting are routinely used in
US voting rights cases [@king1997solution]; therefore, easy to use and high quality tools for performing ecological inference are of practical interest. There is a need for an ecological inference library that 
brings together a variety of ecological inference methods in one place to facilitate
crucial tasks such as: quantifying the uncertainty associated with ecological inference
results under a given model; making comparisons between methods; and bringing relevant 
diagnostic tools to bear on ecological inference methods. To address this need, 
we introduce `PyEI`, a Python package for ecological inference. 

`PyEI` is meant to be useful to two main groups of researchers. First, it serves application-oriented researchers and practitioners who seek to run ecological inference on domain data (e.g., voting data), report the results, and understand the uncertainty related to those results.
Second, it facilitates exploration and benchmarking for researchers who are seeking to understand properties of existing
ecological inference methods in different settings and/or develop new statistical methods for ecological inference.

`PyEI` brings together the following ecological inference methods in a common framework alongside plotting, reporting, and diagnostic tools:

- Goodman's ecological regression [@goodman1953ecological] and a Bayesian linear regression variant
- A truncated-normal based approach [@king1997solution]
- Binomial-Beta hierarchical models [@king1999binomial]
- Dirichlet-Multinomial hierarchical models [@rosen2001bayesian]
- A Bayesian hierarchical method for ${2 \times 2}$ EI following the approach of @wakefield2004ecological

In several of these cases, `PyEI` includes modifications to the models as originally proposed in the cited literature, such as reparametrizations or other changes to upper levels of the hierarchical models in order to ease sampling difficulties.

`PyEI` is intended to be easily extensible, so that additional methods from the literature can continue to be incorporated (for example, work is underway to add the method of @greiner2009r, currently implemented in the R package `RxCEcolInf` [@RxCEcolInf]). Newly developed statistical methods for ecological inference can be included and conveniently compared with existing methods.

Several R libraries implementing different ecological inference methods exist, such as `eiPack` [@eiPack], `RxCEcolInf` [@RxCEcolInf], `ei` [@ei], and `eiCompare` [@eiCompare]. In addition to presenting a Python-based option that researchers who primarily use Python may appreciate, `PyEI` 
incorporates the following key features and characteristics.

First, the Bayesian hierarchical methods implemented in `PyEI` rest on modern probabilistic programming tooling [@salvatier2016probabilistic] and gradient-based MCMC methods such as the No U-Turn Sampler (NUTS) [@hoffman2014no; @betancourt2018conceptual]. Using NUTS where possible should allow for faster convergence than existing implementations that rest primarily on Metropolis-Hastings and Gibbs sampling steps. Consider effective sample size, which is a measure of how the variance of a Monte Carlo estimate of a posterior expectation computed from dependent samples compares to the variance of the corresponding estimate computed from independent samples from the posterior distribution (or, very roughly, how “effective” the samples are for estimating a posterior expectation, compared to independent samples) [@BDA3]. Under certain assumptions on the target posterior distribution, in Metropolis-Hastings the number of evaluations of the log-posterior required for a given effective sample size scales linearly with the dimensionality of the parameter space, while in Hamiltonian Monte Carlo approaches such as NUTS, the number of required evaluations of the gradient of the log-posterior scales only as the fourth root of the dimension [@neal2011mcmc]. Reasonable scaling with the dimensionality of the parameter space is important in ecological inference, as that dimensionality is large when there are many precincts.

Second, integration with the existing tools `PyMC3` [@salvatier2016probabilistic] and `ArviZ` [@arviz_2019] makes the results amenable to state of the art diagnostics (e.g. convergence diagnostics) and some reasonable checks are automatically performed. 
 
Third, summary and plotting utilities for reporting, visualizing, and comparing results are included (e.g. \autoref{fig:kdes}, \autoref{fig:polarization}), with an emphasis on visualizations and reports that clarify the uncertainty of estimates under a model.

Lastly, clear documentation is provided, including a set of introductory and example notebooks.

![Kernel density estimation plots for visualizing uncertainty of support for candidates within each group.\label{fig:kdes}](figs/figure2.png){ width=100% } 

![Visualizing and quantifying degree of polarization.\label{fig:polarization}](figs/figure4.png){ width=100% }

# Acknowledgments

This software development is part of a research project comparing methods, joint with Moon Duchin and Thomas Weighill. We thank Colin Carroll, JN Matthews, and Matthew Sun for their helpful contributions to `PyEI`. 


<!-- ![Bayesian credible intervals for support of candidates within groups.\label{fig:credible_interval}](figs/figure3.png){ width=100% } -->

<!-- ![Visualizing estimates and uncertainty for precinct-level estimates.\label{fig:precinct_level}](figs/figure5.png){ width=50% }

![Tomography plots for two-by-two ecological inference.\label{fig:tomography}](figs/figure6.png){ width=40% } -->

# References

