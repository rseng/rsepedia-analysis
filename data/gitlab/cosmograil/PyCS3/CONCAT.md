# PyCS3


[![Documentation Status](https://cosmograil.gitlab.io/PyCS3/doc_status.svg)](https://cosmograil.gitlab.io/PyCS3/)
[![pipeline status](https://gitlab.com/cosmograil/PyCS3/badges/master/pipeline.svg)](https://gitlab.com/cosmograil/PyCS3/commits/master)
[![coverage report](https://gitlab.com/cosmograil/PyCS3/badges/master/coverage.svg)](https://cosmograil.gitlab.io/PyCS3/coverage/)
[![Python 3.7](https://img.shields.io/badge/python-3.7-blue.svg)](https://www.python.org/downloads/release/python-370/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02654/status.svg)](https://doi.org/10.21105/joss.02654)


PyCS3 is a software toolbox to estimate time delays between multiple images of gravitationally lensed quasars, developed within the [COSMOGRAIL](http://www.cosmograil.org) collaboration. This is an update of [PyCS](https://github.com/COSMOGRAIL/PyCS), which is no longer maintained. 


Proceed to the [documentation](https://cosmograil.gitlab.io/PyCS3/) to get further information. In case of any questions, feel free to open an issue here on GitLab.

## Installation 

    git clone https://gitlab.com/cosmograil/PyCS3
    cd PyCS3 
    python setup.py install

or if you prefer to install it locally : 

    python setup.py install --user 
    
## Requirements 

PyCS3 requires the following standard python packages : 
* `numpy`
* `scipy`
* `matplotlib`
* `multiprocess`

If you want to use the regdiff optimiser, you will also need : 
* `scikit-learn`
    
## Example Notebooks and Documentation
The full documentation can be found [here](https://cosmograil.gitlab.io/PyCS3/). 

Example notebooks are located in the [notebook](https://gitlab.com/cosmograil/PyCS3/-/tree/master/notebook) folder : 
* [Importing, exporting and displaying light curves](https://gitlab.com/cosmograil/PyCS3/-/blob/master/notebook/Import_export_and_display.ipynb)
* [Measuring time delays with regdiff and the splines](https://gitlab.com/cosmograil/PyCS3/-/blob/master/notebook/Measuring%20time%20delays%20with%20spline%20and%20regdiff.ipynb)
* [Estimating uncertainties with PyCS3](https://gitlab.com/cosmograil/PyCS3/-/blob/master/notebook/Uncertainties%20estimation.ipynb)

## Attribution

If you use this code, please cite [the papers](https://cosmograil.gitlab.io/PyCS3/citing.html) indicated in the documentation.

## License
PyCS3 is a free software ; you can redistribute it and/or modify it under the terms of the 
GNU General Public License as published by the Free Software Foundation ; either version 3 
of the License, or (at your option) any later version.

PyCS3 is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY ; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
General Public License for more details ([LICENSE.txt](LICENSE)).

# Contributing


You can contribute to PyCS3 in many ways and this is very appreciated. 

## Report bugs 

Please use GitLab issues to report bug.If you report a bug, please include:

* Your operating system name and version.
* Detailed steps to reproduce the bug.

## Implement new features 

Please use pull requests if you want to modify the master branch. 

* The pull request should include tests.
* If the pull request adds functionality, the docs should be updated. Put
  your new functionality into a function with a docstring. 

# PyCS3 test pipeline


This folder contains all the script to process multiple light curves from the TDC at the same time. It aims at providing a test framework to check the precision and accuracy of PyCS3. It is based on the original pipeline which is in the `script` folder. A list of light curves from the TDC1 ([Liao et al. (2015)](https://arxiv.org/abs/1409.1254)) that can be used with this sub-package are available [here](https://lsstdesc.org/TimeDelayChallenge/downloads.html).

Four metrics were selected in the Time-Delay challenge to evaluate the performance of the curve-shifting techniques : 
* the accuracy 
* the precision 
* the $`\chi^2`$ 
* X, the fraction of estimates without outliers (with a $`\chi^2`$ < 10). 

We recently ran this test pipeline on the first 200 curves of the rung 3 of the TDC1. Those data closely mimics the real Euler light curves in terms of cadence, photometric noise and microlensing. We reproduce here Figure 8 of [Liao et al. (2015)](https://arxiv.org/abs/1409.1254), which summarizes the 4 metrics in one plot : 

![](figure/SS_Final_Plot.png)*TDC1 rung 3 results. The dots are blind submissions from [Liao et al. (2015)](https://arxiv.org/abs/1409.1254). Our new results with PyCS3 and both the spline and regression difference estimators are marked with triangles. The different strategies to marginalise between estimator parameters are parametrized with the `sigma` threshold parameter (see. [Millon et al. (2020)](https://arxiv.org/abs/2002.05736) for details.)*

The results presented here are for our *Silver Sample* (SS), which contains more than 60% of the curves. We excluded the curves with a precision >40% and ones with a time delay > 100 days, which does not have enough overlap to measure a robust time delay. Of course, the results were not *blinded* here, but this demonstrates that `PyCS3` was able to measure the time-delay in most of curves of the TDC1 rung3 with an automated procedure with an accuracy better than 2%. 


## Using the pipeline 

You should first define a working directory that must contain your light curves in the sub-folder `data`. 
Your data should have the naming convention `name_double_pairX_ECAM.rdb` with X the number of the light curve and must contain one column for the days of the observation and two columns per lens image for the measured magnitudes and its associated uncertainties. 
You also need to have the true time delay for your light curves in the sub-folder `data/truth` with the naming convention `truth_name.txt` where the second column are the true time delay like in the following example :
        test_double_pair1     -70.89
        test_double_pair2     104.22
        test_double_pair3     56.59


## 0. Set up
 
Run the command : 

    python3 multiple_0.py name double numberofcurves --dir='path_to_working_directory'

This will set up the directories for this multiple pipeline and generate the multiple config file in the `./config/multiple` directory under the name `config_multiple_name.py` and can be modified before running the next script. Check the [README](../scripts/README.md) from the main pipeline for more information on the different parameters.

## 1. Create data set 

Run the command : 

    python3 multiple_1.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `1_create_dataset.py` from the PyCS3 main pipeline for all the curves. The specific config file for each of the curve is created in `config` directory under the name `config_name_double_pairX_EXAM.py`. The single config files are then updated with the parameters from the multiple config file previously mentioned.
The initial guess will be randomly taken around the true time delay from the `./data/truth/truth_name.txt`.
The initial guesses are saved in `./Simulation/multiple/name_double/post_gaussian_guess.txt`. 

## 2. Fit spline 

Run the command : 

    python3 multiple_2.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `2_fit_spline.py` from the PyCS3 main pipeline for each curve.

## 3. Launch spline 

Modify the `./cluster/multiple_launcher.sh` with the parameters from your cluster. Then make sure the line 18 execute the `start_3all.slurm` and update the number of curves you will run on the line 28. Then run the command while in the `./cluster` folder :

    ./multiple_launcher.sh

This will run the script `3a_generate_tweakml.py`, `3b_draw_copy_mocks.py`, `3c_optimise_copy_mocks.py` and `3d_check_statistics.py` from the PyCS3 main pipeline for each curves.

## 4a. Get the time delay for the spline 

Run the command : 

    python3 multiple_4a.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `4a_plot_results.py` from the PyCS3 main pipeline for each curve.

##  3bis. Launch regdiff 

Run the command : 

    python3 multiple_update_config.py name double numberofcurves --dir='path_to_working_directory'

Update every config file in order to use the regdiff instead of the spline.
Then modify the `./cluster/multiple_launcher.sh` with the parameters from your cluster. Then make sure the line 18 execute the `start_3c.slurm` and update the number of curves you will run on the line 28. Then run the command while in the `./cluster` folder :

    ./multiple_launcher.sh

This will run the script `3c_optimise_copy_mocks.py` from the PyCS3 main pipeline for each curves, but this time with regdiff.

## 4a. Get the time delay for the regdiff

Run the command : 

    python3 multiple_4a.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `4a_plot_results.py` from the PyCS3 main pipeline for each curve, but this time with regdiff.

## 4b. Margininalise the spline Estimates 

Run the command : 

    python3 multiple_4b.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `4b_marginalise_spline.py` from the PyCS3 main pipeline for each curve for three different `sigmathresh = 0, 0.5, 1000`. You can find information on the marginalisation in the README of the main pipeline.

## 4b. Margininalise the regdiff Estimates 

Run the command : 

    python3 multiple_4c.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `4c_marginalise_regdiff.py` from the PyCS3 main pipeline for each curve for three different `sigmathresh = 0, 0.5, 1000`. You can find information on the marginalisation in the README of the main pipeline.

## 4d. Marginalise spline and regression difference Estimates 

Run the command : 

    python3 multiple_4d.py name double numberofcurves --dir='path_to_working_directory'

This will run the script `4d_marginalise_rall.py` from the PyCS3 main pipeline for each curve with both the `sigmathresh = 0.5`. 

## 6. Compute TDC1 metrics
You need to update the `config_multiple_name.py` with the simulation that failed in the `failed_sim` list. You can then choose if you want to display the gold and silver sample and add the pair you want in the silver/golden sample.
You can also decide what plot you want to have displayed. Then, run the command :

    python3 multiple_6.py name double numberofcurves --dir='path_to_working_directory'

This will compute the statistics for each sample, and for each estimate. Summary files are created in the `./Simulation/multiple/name` folder. A sub-folder `figure` is also created and contains all the figure produced.
This folder contains the copies and mock light curves used to measure time-delay uncertainties. This is where all the products of the example notebook are saved. ---
title: 'PyCS3: A Python toolbox for time-delay measurements in lensed quasars'
tags:
  - Python
  - astronomy
authors:
  - name: Martin Millon
    orcid: 0000-0001-7051-497X
    affiliation: 1
  - name: Malte Tewes
    orcid: 0000-0002-1155-8689
    affiliation: 2
  - name: Vivien Bonvin
    orcid: 0000-0003-1471-3952
    affiliation: 1
  - name: Bastian Lengen
    affiliation: 1
  - name: Frederic Courbin
    orcid: 0000-0003-0758-6510
    affiliation: 1
affiliations:
  - name: Institute of Physics, Laboratory of Astrophysique, Ecole Polytechnique Fédérale de Lausanne (EPFL), Switzerland
    index: 1
  - name: Argelander-Institut für Astronomie, Bonn, Germany
    index: 2 
date: 24 June 2020
bibliography: paper.bib
---

# Summary
Time-delay cosmography is a competitive technique for measuring the current expansion rate of the Universe, that is, the Hubble Constant [see, e.g., @Riess2019 for a review; @Wong:2019 for recent results]. It relies on the strong gravitational lensing effect that happens when a massive foreground galaxy deviates the light from a background object, producing 2 or 4 mirage images of the same background source. In this configuration, the optical path length is slightly different in each multiple image and thus, the travel time of the photons along those paths is also slightly different. If the background source is varying, the same variations are visible in all multiple images with different delays. Lensed quasars or lensed supernovae are ideal targets to measure such (relative) time delays, because they are variable on short timescale, and are sufficiently bright to be observed at cosmological distance. These measured delays can be used to infer the so-called time-delay distance, $D_{\Delta t}$, which is directly inversely proportional to the Hubble Constant. The method relies on three main ingredients : 

 - a precise and accurate determination of the time delays
 - a model of the mass distribution of the lensing galaxy 
 - an estimate of the mass of all the galaxies along the line of sight that also deviates the light rays, and thus perturbs the time delays. 
 
 Obtaining the time delays of lensed quasars typically requires a decade of continuous observation to produce long light curves which contain several variations of the quasars that can unambiguously matched between the multiple images. An example of the light curves of the multiple images of the lensed quasar RXJ1131-1231 is presented in \autoref{fig:lcs}. The aim of the PyCS3 software is to measure time delays between such curves.

![Light curves of the lensed quasar RXJ1131-1231 presented in @Millon1:2020 (left panel). The same quasar variations can be seen in image D 92 days after in image A, whereas images A, B, and C arrive approximately at the same time. The right panel shows an Hubble Space Telescope image of RXJ1131-1231 [@Suyu2017].\label{fig:lcs}](RXJ1131.png)
 
 
# Statement of need
 The "simple" problem of measuring time delays between irregularly sampled light curves has received attention for almost three decades [e.g., @Press1992]. A complicating factor is the microlensing of the multiple images which happens when stars in the lens galaxy are passing in front of the quasar images, also acting as gravitational lenses and affecting the shape of the light curves. Microlensing effects perturb the light curve of each quasar image individually. Ignoring these perturbations would often result in significant biases on the measured time delays, that directly propagate to the Hubble Constant.
 
 ``PyCS3`` is a python package developed by the [COSMOGRAIL](www.cosmograil.org) collaboration and designed to address this problem. It allows us to measure time delays in lensed quasars in the presence of microlensing by providing a flexible and data-driven model of the "extrinsic" variations with splines to account for microlensing and recover an accurate estimate of the time delays. A realistic estimation of the uncertainties is also extremely important for cosmography. The approach followed by ``PyCS3`` is based on best-fit point estimation and a framework to faithfully simulate the input data and assess the uncertainties following a Monte Carlo approach. While a Bayesian inference of the delay, given a model for quasar variability and microlensing, would be very attractive, it is hampered by the difficulty to accurately model microlensing and the high number of nuisance parameters in the problem.
 
  The previous version of the package (``PyCS``) was first presented in @Tewes1:2013 and successfully applied to real data in @Tewes2:2013; @Bonvin:2017; @Bonvin:2018 and @Bonvin:2019. The method was also tested on simulated light curves of the Time-Delay Challenge [@Liao:2016; @Bonvin:2016] and was empirically shown to provide excellent results, even if the light curves are strongly affected by microlensing.
   
We have now developed an automated pipeline based on ``PyCS3`` to measure time delays in a large sample of lensed quasars [@Millon1:2020; @Millon2:2020]. Such improvements toward automation of the procedure is necessary with the hundreds of new lensed quasars expected to be discovered in the near future. 

# Functionality
 The basic functionality of ``PyCS3`` is built around a LightCurve class to manipulate photometric monitoring data. It has methods to import, shift, fit and export light curves. These are located in the `pycs3.gen` subpackage.

``PyCS3`` contains two time-delay estimators, namely the free-knot splines and the regression difference, that are in the `pycs3.spl` and `pycs3.regdiff` subpackages. These two estimators are fundamentally different and allows us to check the robustness of the measured time delays. The subpackage `pycs3.sim` is used to generate simulated light curves in order to estimate the uncertainties of the time-delay measurements. ``PyCS3`` ensures that the simulated curves have the same constraining power than the original data which is crucial for a correct estimation of the uncertainties. These simulated curves can then be shifted with either the free-knot splines or the regression difference estimator. 

The ``script`` folder contains a pipeline to automate the measurement of the time delays. The functions that are used by this pipeline are located in the `pycs3.pipe` subpackage. It automatically explores several set of estimator parameters, generates simulated light curves, returns the best fit value and the associated uncertainties before selecting and combining the different sets of estimator parameters. It also makes use of the subpackage `pycs3.tdcomb` to display and combine the final time-delay estimates.  The details of the method can be found in @Millon1:2020. 

Finally, the ``tdlmc_test`` folder contains an ensemble of python scripts to apply ``PyCS3`` on the simulated data generated for the Time-Delay Challenge [@Liao:2016]. It provides a benchmark test framework to assess the precision and accuracy of the ``PyCS3``'s curve-shifting algorithms. 

# Acknowledgement

We acknowledge the support of the Swiss National Science Foundation (SNSF) and the European Research Council (ERC) under the European Union’s Horizon 2020 research and innovation programme (COSMICLENS: grant agreement No 787886).

# References
# PyCS3 pipeline

This folder contains all the script to process the light curves and measure the time delays. The complete description of the pipeline can be found in [Millon et al. (2020)](https://arxiv.org/abs/2002.05736).

You should first define a working directory that must contain your light curves in the sub-folder `data`. 
Your data should have the naming convention `lensname_dataset.rdb` and must contain one column for the days of the observation and two columns per lens image for the measured magnitudes and its associated uncertainties. A dataset typically corresponds to a monitoring campaign conducted by a single telescope.  

## 1. Create data set 

Run the command : 

    python3 1_create_dataset.py lensname dataset --dir='path_to_working_directory'

to setup the config file. It will create a new folder for your lens and dataset under the subfolder `Simulation/lensname_dataset` from your working directory. If you are not specifying a working directory, it will simply use the current folder `'./'`. The config file is created in the `config` directory under the name `config_lensname_dataset.py` and can be modified. 

## 2. Fit spline 
The first step consists in fitting splines to your original data. You can change the mean knotstep and the microlensing model from the config file. The key parameters are : 

   - `knotstep` : mean initial spacing between the knots of the intrinsic spline. Give a list with the value that you want to test. 
   - `mlknotsteps` or alternatively `nmlspl` : if you use the option `forcen`, you will use the value given in the  `nmlspl` list, which then contains the number of knots of the extrinsic splines, equally distributed over the monitoring period. These knots have a fixed position. This option is recommended for short light curves to ensure that the microlensing model do not have too much freedom. Note that using `nmlspl`=1, will place 0 internal knot and the microlensing is modelled with a polynomial of degree 3. If you have long curves (more than one season), we recommend to turn off the `forcen` option and use `mlknotsteps` instead of `nmlspl`. You need to provide the mean initial spacing between the knots of the microlensing spline that are then free to move. If you insert 0 in the `nmlspl` or `mlknotstep` lists, this will not introduce any microlensing model.
   - `preselection_file` if `use_preselected_regdiff` : provide directly the regdiff set of parameters in dictionaries saved in a `json` readable file. Alternatively, you can directly provide the regdiff parameter that you want to test in the `covkernel`, `pointdensity`, `pow` and `errscale` lists. 
  
   
Once you have chosen the estimator parameters, you can run the command : 

    python3 2_fit_spline.py lensname dataset --dir='path_to_working_directory
    
This will produce the original fit with all the combination of intrinsic and microlensing models that you provided. You can check visually the quality of the fit with the figures saved in the figure/spline_and_residuals_plots. 

## 3. Generate mock light curves
#### 3a. Fit the parameter of the generative noise model 
This is the tricky part of the process. We will need to adjust two parameters to ensure that the generative noise model creates simulated curves with the same constraining power than the original data. The fit of these two parameters is done with a simple dichotomous search algorithm. The key parameters to control the behaviour of the algorithm are : 

 - `n_curve_stat` : number of curves to compute statistics. A large number of curves will reduce the random fluctuations and accelerate the convergence but this comes at the price of heavy computation. Curves are computed in parallel, so you might want to choose a multiple of the number of cores available. For example, 16 is usually good number. 
 - `max_iter` : maximum number of iteration before stopping. If the algorithm finds a good solution before that, it will stop automatically. 
 
 The algorithm will try to find parameters that produces curves with the  mean *rms residuals* &sigma; and the mean *number of runs* z<sub>run</sub> that falls within 0.75-&sigma; from the real curves. This script will save the functions and optimal parameters in a python file named `tweakml_ps.py`. Now, you can run :
 
    python3 3a_generate_tweakml.py lensname dataset --dir='path_to_working_directory
 
 __Warning__ : For quad, it happens regularly that the algorithm do not match the 0.75-&sigma; criterion in 15 iteration. It will then choose the closest value. Matching &sigma; and z<sub>run</sub> within 1 or 2&sigma; is sufficient in 99% of the cases. You can still continue until script 3d and check a posteriori that the simulated curves are indeed not too different. 

#### 3b. Draw mock lights from the generative noise model 
The next step is to generate the copies and mocks. You can adjust the number of copies and mocks generated by changing `ncopy`, `ncopypkls`, `nsim` and `nsimpkls`. The total number of copies generated is `ncopy`*`ncopypkls`, for every combination of `knotstep` and `mlknotstep`, so the computation time increases very quickly. We recommend to generate 500 copies and 800 mocks for reliable estimates of the uncertainties. 

     python3 3b_draw_copy_mocks.py lensname dataset --dir='path_to_working_directory
     
If you turn on the `ask_question` option, you can run this script many times to generate more mocks. By default, it will erase the already existing files. 

#### 3c. Optimise the mock light curves 
We will now fit the copies and mocks. If you want to run only on the mocks, turn off the `run_on_copies` option. By default, you need to run the optimiser on both the mocks and copies. You can choose which optimiser to use with the `simoptfctkw`. You can now run : 

    python3 3c_optimise_copy_mocks.py lensname dataset --dir='path_to_working_directory
    
This script will optimise batches of curves saved in pickle files. You can run on many cores in parallel by using the `max_core` option. 
You can run this script twice, once with the spline optimiser, once with the regdiff optimiser. 

#### 3d. Check statistics 
You can now check *a posteriori* that the generated curves roughly match the original data in term of z<sub>run</sub> and &sigma;. To do so, run : 

    python3 3d_check_statistics.py lensname dataset --dir='path_to_working_directory
    
This will create a series of plots in the figure/check_stat_plots folder. You can now check that the residuals of the mock curves looks similar to the residuals of the data. 

## 4. Marginalise over the estimator parameters

#### 4a. Get the time delay Estimates 
The first step is to collect the results of the optimisation of the mocks and copies. This is done by running : 

    python3 4a_plot_results.py lensname dataset --dir='path_to_working_directory
    
You should run that script once for each estimator. You can switch from one estimator to the other with the `simoptfctkw` parameter. 
This will creates plots in the figure/final_results folder. You can verify that the distribution of the delay when running on the copies is peaked. This indicates that your data indeed constrain the time delays. You can also visualize the random and systematic errors made on the mocks.

#### 4b. Margininalise the spline Estimates 
To obtain the final spline estimate, you can either select the most precise model, marginalize over all available models or follow the hybrid approach described in [Millon et al. (2020)](https://arxiv.org/abs/2002.05736). This can be controlled with the `sigmathresh` parameter (0 is a true marginalisation, 0.5 is the hybrid approach, and some very high number will pick only the most precise estimate). You can ensure a high numerical precision by turning off the `testmode` option but this will take a bit longer. 
You can name your combination with the `name_marg_spline` option. This is useful if you want to combine this estimate with the something else in the future. You can also select manually the model to be included in the combination. To do so, you need to provide the list of parameters in `knotstep_marg` and `mlknotsteps_marg` and run : 

    python3 4b_marginalise_spline.py lensname dataset --dir='path_to_working_directory

By default, the code will take the models provided in your `knotstep` and `mlknotsteps` or `nmlspl` lists. 

#### 4c. Marginalise the regression difference Estimates
Similarly to the spline estimate, you can choose which generative noise model to use with the `knotstep_marg_regdiff` and `mlknotsteps_marg_regdiff` lists. By default we choose the generative noise model that gives the most precise estimates and we use `sigmathresh` to combine the different set of regdiff parameters. 

    python3 4c_marginalise_regdiff.py lensname dataset --dir='path_to_working_directory

#### 4d. Marginalise spline and regression difference Estimates 

The last step is to combine regdiff and the spline together. You can choose which estimate to combine in the `name_marg_list` list and the threshold that was used at the previous steps in `sigmathresh_list`. We also recommend to marginalise over the two estimators, which correspond to use `sigmathresh_final = 0`. Once you decided which estimator to include in your final measurement, you can run : 
    
    python3 4d_marginalise_all.py lensname dataset --dir='path_to_working_directory

## 5. Combine Estimates from different data set. 
If you have several data sets for the same object, you might be interested in combining the data set together. You can then run 

     python3 5_combine_dataset.py lensname --dir='path_to_working_directory
     
It will then create the folder Combination/lensname. You can edit the file config_combination_lensname.py located in this folder and rerun the script 5. 
