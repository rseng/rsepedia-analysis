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
Warning & Disclaimer
====================


By *itself*, the use of PyCS3 does not guarantee realistic time-delay uncertainty estimates. All we have (blindly) demonstrated regarding the quality of PyCS3 uncertainty estimates assumes a careful generation of the simulated light curves. Keep in mind that the PyCS3 approach involves setting estimator parameters or ranges of estimator parameters. In principle, the impact of these settings has to be assessed on each light curve set, as described in the first PyCS paper (`Tewes et al. 2013 <http://dx.doi.org/10.1051/0004-6361/201220123>`_). The estimator parameters can also be marginalised over as described in `Millon et al. 2020 <https://arxiv.org/abs/2002.05736>`_.

.. warning:: In particular, obtaining uncertainty estimates by running on simulated light curves that all have a same true delay (i.e., setting ``truetsr = 0``, in PyCS terminology) can lead to **vastly underestimated error bars**.

It is of **high importance** to use simulations with a range of plausible true time delays to tune and/or verify the accuracy and precision of time-delay estimators. Tests on simulations with only a single true time delay do not probe reliably the quality of a time-delay estimation, as many time-delay estimators are prone to responding unsteadily to the true delay. Do not hesitate to contact us (see front page of this documentation) in case of doubts!

.. note:: Please try to avoid publishing time delays with overly optimistic uncertainty estimates. Seemingly accurate measurements might be propagated into unrealistic cosmological inferences, and harm the community. Be critical, and aim for the most comprehensive instead of the smallest error bar.

Thanks!Citing PyCS3 in a publication
=============================

If you want to acknowledge PyCS3 in a publication, we suggest you to cite both the original PyCS paper and the new version of the software presented in the Journal of Open Source Software (JOSS) :

* Tewes et al. 2013

  - `COSMOGRAIL XI: Techniques for time delay measurement in presence of microlensing, A&A 553 A120 <http://dx.doi.org/10.1051/0004-6361/201220123>`_.
  - This paper describes the original curve-shifting algorithms of PyCS.


* Millon et al. 2020b

  - `PyCS3: A Python toolbox for time-delay measurements in lensed quasars, JOSS 5 53 2654 <https://doi.org/10.21105/joss.02654>`_.
  - This paper briefly describes the functionality of the package


If you make use of the automated time-delay measurement pipeline, please also cite :

* Millon et al. 2020a

  - `COSMOGRAIL XIX: Time delays in 18 strongly lensed quasars from 15 years of optical monitoring <https://arxiv.org/abs/2002.05736>`_.
  - This paper presents the automated version of PyCS3 and apply this technique on large monitoring data set.



.. note:: Please use the URL ``http://www.cosmograil.org`` when indicating a website for PyCS3. While the hosting of the code might change, we will try hard to keep this URL valid for the years to come.Intra-dependency chart
======================

This chart lists the intra-dependencies of the various modules and files in PyCS.


.. image:: _static/pycs3.svg
    :align: center

This plot was generated with ``pydeps`` (`documentation <https://pydeps.readthedocs.io/en/latest/>`_).
Last updated on 2020-06-10. Open the image in a new tab for a better experience.
Published work with PyCS
========================

You can find the published articles making use of PyCS below. Please let us know if you have published a paper that make use of PyCS. We will include it in this list.

PyCS methodology and software publication
-----------------------------------------
* Tewes et al. 2013a

  - `COSMOGRAIL - XI. Techniques for time delay measurement in presence of microlensing, A&A 553 A120 <http://dx.doi.org/10.1051/0004-6361/201220123>`_.
  - This paper describes the original curve-shifting algorithms of PyCS.

* Millon et al. 2020c

  - `PyCS3: A Python toolbox for time-delay measurements in lensed quasars, JOSS 5 53 2654 <https://doi.org/10.21105/joss.02654>`_.
  - This paper briefly describes the functionality of the package.


Measuring time delays in lensed quasars
---------------------------------------

* Tewes et al. 2013b

  - `COSMOGRAIL - XIII. Time delays and 9-yr optical monitoring of the lensed quasar RX J1131−1231, A&A 556 A22 <https://doi.org/10.1051/0004-6361/201220352>`_
  - This work presents the application of PyCS to real data. It reports the measurement of time delays in the lensed quasar RXJ1131-1231 from the COSMOGRAIL data.

* Eulaers et al. 2013

  - `COSMOGRAIL - XII. Time delays of the doubly lensed quasars SDSS J1206+4332 and HS 2209+1914, A&A 553 A121 <https://doi.org/10.1051/0004-6361/201321140>`_
  - Time-delay measurements of lensed quasars SDSS J1206+4332 and HS 2209+1914 with data from the Swiss Euler telescope, Himalayan Chandra Telescope (HCT) and Mercator telescope.

* Rathna Kumar et al. 2013

  - `COSMOGRAIL - XIV. Time delay of the doubly lensed quasar SDSS J1001+5027, A&A 557 A44 <https://doi.org/10.1051/0004-6361/201322116>`_
  - Time-delay measurements of lensed quasars SDSS J1001+5027 with data from the Himalayan Chandra Telescope (HCT), Mercator and 1-5m Maidanak Observatory telescopes.


* Goiecoechea and Shalyapin 2016

  - `Gravitational lens system SDSS J1339+1310: microlensing factory and time delay, A&A 596 A77 <https://doi.org/10.1051/0004-6361/201628790>`_
  - Presentation of the delay of SDSS J1339+1310 and of the microlensing variability affecting this system.

* Giannini et al. 2017

  - `MiNDSTEp differential photometry of the gravitationally lensed quasars WFI 2033-4723 and HE 0047-1756: microlensing and a new time delay, A&A 597 A49 <https://doi.org/10.1051/0004-6361/201527422>`_
  - Time-delay measurement of lensed quasars HE 0047-1756 and WFI 2033-4723.

* Bonvin et al. 2017

  - `H0LiCOW – V. New COSMOGRAIL time delays of HE 0435−1223: H0 to 3.8 per cent precision from strong lensing in a flat ΛCDM model, MNRAS 465 4 <https://doi.org/10.1093/mnras/stw3006>`_
  - This paper presents the new measurement of Hubble constant by the `H0LiCOW <https://shsuyu.github.io/H0LiCOW/site/>`_ collaboration.

* Courbin et al. 2017

  - `COSMOGRAIL - XVI. Time delays for the quadruply imaged quasar DES J0408−5354 with high-cadence photometric monitoring, A&A 609 A71 <https://doi.org/10.1051/0004-6361/201731461>`_
  - This paper presents the time delays measurement of DES0408 with only one season of monitoring at high-cadence and high signal-to-noise.

* Bonvin et al. 2018

  - `COSMOGRAIL - XVII. Time delays for the quadruply imaged quasar PG 1115+080, A&A 616 A183 <https://doi.org/10.1051/0004-6361/201833287>`_
  - Time-delay measurement of lensed quasar PG1115+080. This papers also introduce a new framework for combining the PyCS estimator parameters.

* Birrer et al. 2019

  - `H0LiCOW - IX. Cosmographic analysis of the doubly imaged quasar SDSS 1206+4332 and a new measurement of the Hubble constant, MNRAS 484 4 <https://doi.org/10.1093/mnras/stz200>`_
  - New measurement of the Hubble Constant by the `H0LiCOW <https://shsuyu.github.io/H0LiCOW/site/>`_ collaboration, using the lens system SDSS J1206+4332.

* Bonvin et al. 2019

  - `COSMOGRAIL - XVIII. time delays of the quadruply lensed quasar WFI2033−4723, A&A 629 A97 <https://doi.org/10.1051/0004-6361/201935921>`_
  - This paper present the time delays of lensed quasar WFI 2033-4723 from 14 years of data taken at the Euler Swiss telescope, 13 years from the SMARTS telescope and from one season of high-cadence data at the MPIA 2m2 telescope.

* Millon et al. 2020a

  - `COSMOGRAIL - XIX. Time delays in 18 strongly lensed quasars from 15 years of optical monitoring, A&A 640 A105 <https://doi.org/10.1051/0004-6361/202037740>`_
  - This paper reports the measurement of time-delay in 18 lensed quasars. It also presents the new PyCS automated time-delay measurement pipeline.

* Millon et al. 2020b

  - `TDCOSMO - II. 6 new time delays in lensed quasars from high-cadence monitoring at the MPIA 2.2m telescope, accepted in A&A <https://arxiv.org/abs/2006.10066>`_
  - Data release and time-delay measurements of 6 lensed quasars monitored at high cadence and high signal-to-noise from the MPIA 2m2 telescope.


Time-Delay Challenge (TDC)
--------------------------

* Liao et al. 2015

  - `Strong Lens Time Delay Challenge: II. Results of TDC1, ApJ 800 1 <https://doi.org/10.1088/0004-637X/800/1/11>`_
  - This paper presents the results of the TDC and compares the performance of several curve-shifting algorithms

* Bonvin et al. 2016

  - `COSMOGRAI - XV. Assessing the achievability and precision of time-delay measurements, A&A 585 A88 <https://doi.org/10.1051/0004-6361/201526704>`_
  - This paper shows the performance of PyCS on the TDC simulated light curves.



Reverberation Mapping
---------------------

* Chan et al. 2020

  - `Twisted quasar light curves: implications for continuum reverberation mapping of accretion disks, A&A 636 A52 <https://doi.org/10.1051/0004-6361/201935423>`_
  - This paper explores the implication of the deformation of the light curves between different bands when measuring time delays.

Microlensing studies
--------------------

* Sluse & Tewes 2014

  - `Imprints of the quasar structure in time-delay light curves: Microlensing-aided reverberation mapping, A&A 571 A60 <https://doi.org/10.1051/0004-6361/201424776>`_
  - This work demonstrates that microlensing can help to disentangle the light coming from the accretion and light being reverberated in the Broad Line Region (BLR), allowing the measurement of the size of the BLR.

* Cornachione et al. 2020

  - `A Microlensing Accretion Disk Size Measurement in the Lensed Quasar WFI 2026-4536, ApJ 895 2 <https://doi.org/10.3847/1538-4357/ab557a>`_
  - This paper presents the accretion disk size measurement of the lensed quasar WFI 2026-4536 using microlensing.


Measuring time delays in lensed Supernovae
------------------------------------------

* Rodney et al. 2016

  - `SN Refsdal: Photometry and time delay measurements of th efirst Einstein cross supernovae, ApJ 820 1 <https://doi.org/10.3847/0004-637X/820/1/50>`_
  - This paper presents the measurement of the time delays of supernovae Refsdal in the Hubble Frontier Field Cluster MACS J1149

* Huber et al. 2019

  - `Strongly lensed SNe Ia in the era of LSST: observing cadence for lens discoveries and time-delay measurements, A&A 631 A161 <https://doi.org/10.1051/0004-6361/201935370>`_
  - This papers presents forecast of the precision on time delays of future lensed supernovae that can be reached with the LSST.
Download & Installation
=======================


Dependencies
------------

PyCS3 is developed using python 3.7 and might be fine with older versions.

It requires ``numpy``, ``scipy``, ``matplotlib`` and ``multiprocess``.
Those are all you need to run the free-knot splines.

The regression difference technique have further dependencies:

* `scikit-learn <http://scikit-learn.org>`_


Download
--------

Get the latest PyCS by cloning it from `GitLab <https://gitlab.com/cosmograil/PyCS3>`_::

	git clone https://gitlab.com/cosmograil/PyCS3


Installation
------------

If you  want to update or tweak the sources, we suggest to just add your cloned repository to your ``PYTHONPATH`` or type :

::

    python setup.py develop

If you don't plan to tweak the code, you can also simply

::

	python setup.py install

or maybe

::

	python setup.py install --user

... if you don't have write access to the global site-packages directory of your machine.

Tests
-----

PyCS3 now have automatic tests to verify that everything works correctly. If you want to check your installation, you will first need to install `PyTest` with the command :

::

    pip install pytest --user

Then, you can simply go in the PyCS3 repository and run the command :

::

    pytest

Running all the tests should take between 5 and 10 minutes... PyCS3 documentation master file, created by
   sphinx-quickstart on Fri Aug 10 11:19:31 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PyCS3's documentation!
=================================

.. image:: _static/cover_image_1.png
	:align: center


About
-----

PyCS3 is a software toolbox to estimate time delays between multiple images of strongly lensed quasars, from resolved light curves such as obtained by the `COSMOGRAIL <http://www.cosmograil.org/>`_ monitoring program. It comes in the form of a python package, and heavily depends on ``numpy``, ``scipy``, and ``matplotlib`` for its core functionality. The `repository is on GitLab <https://gitlab.com/cosmograil/PyCS3>`_.

To measure time delays with ``pycs3``, you'll typically write a script calling some high-level functions provided by the package. PyCS3 allows you to compare different point estimators (including your own), without much code integration. You can follow the `example notebooks <https://gitlab.com/cosmograil/PyCS3/-/tree/master/notebook>`_ to learn how to use the core functionnality of PyCS3.

If you have already read our :doc:`papers<citing>`, you might want to proceed with :doc:`installation`, or the :doc:`Tutorial <tutorial/index>`.

.. warning:: Please read this :doc:`important warning about using PyCS3<warning>`.

Questions ?
-----------

Feel free to post an `issue on GitLab <https://gitlab.com/cosmograil/PyCS3/-/issues>`_, or to contact the code authors `Martin Millon <http://people.epfl.ch/martin.millon>`_, `Malte Tewes <https://astro.uni-bonn.de/~mtewes>`_ and `Vivien Bonvin <http://people.epfl.ch/vivien.bonvin>`_.


Contents
--------


.. toctree::
   :maxdepth: 2

	A word of warning <warning>
	Download & Installation <installation>
	Tutorial<tutorial/index>
    Intra-dependency chart <chart>
	Citing <citing>
    Published work with PyCS <publication>
	Autogenerated Full API <apidoc/pycs3>


Last build of this documentation : |today|.


Getting delay histograms and results
====================================


Generalities
------------

To obtain time-delay histograms and point and uncertainty estimates, several steps need to be done; they are presented in the following subsections.
We assume that you have generated some mock curves as described in the previous section.



Run curve shifting algorithms on the mock curves
------------------------------------------------


The wrapper function to run any method on any set of mock curves has the funny name :py:func:`pycs3.sim.run.multirun`.
As input, it takes those pkl files containing simulated curves, made in the previous section. The output of this step are pkl files of runresult objects.

One important argument of ``multirun`` is ``tsrand``. It is the radius (in days) of a uniform randomization of the input time shifts.

Here is an example script running the two methods. We assume that you've defined your favorite optimizers in ``myopt.py``.


::
	
	lcs = pycs3.gen.util.readpickle("merged.pkl")

	
	# Spline method
	"""
	
	# Initial conditions for the analysis :
	# (You might want to set some timeshifts, fluxshifts, as well)
	# In this case we'll just add some spline ML :
	pycs3.gen.splml.addtolc(lcs[0], knotstep=100)
	pycs3.gen.splml.addtolc(lcs[2], knotstep=300)
	pycs3.gen.splml.addtolc(lcs[3], knotstep=300)
	
	#pycs3.sim.run.multirun("copies", lcs, myopt.spl, optset="spl1", tsrand=10.0)
	#pycs3.sim.run.multirun("sim1tsr10", lcs, myopt.spl, optset="spl1", tsrand=10.0)
	"""
	
	# Regdiff method
	"""
	# No need for ml, but time delays should be set about right.
	
	#pycs3.sim.run.multirun("copies", lcs, myopt.regdiff, optset="regdiff2", tsrand=10.0)
	#pycs3.sim.run.multirun("sim1tsr10", lcs, myopt.regdiff, optset="regdiff2", tsrand=10.0)
	"""
	
Parallel computing
------------------
You will probably need more than 500 mock curves to have reliable estimates of the uncertainties. To speed up the computation, you can launch several identical "calls" to this ``multirun`` function (for instance simply by launching your script on several CPUs), and this will indeed process the pkl files in parallel. For this to work, the ``multirun`` function stores a temporary file in its results directory as soon as it starts working on a pickle file, so that other scripts know that they should not run on this same pkl file as well. You can see those temporary files, of course. If something goes wrong and they don't get deleted automatically as the script crashed, you might have to remove the ``.workingon`` files by hand, otherwise ``multirun`` will just skip those pkl files.

Here is an example on how you could perform parallel execution of the ``multirun`` function :

::

    from multiprocess import Pool
    import time

    def exec_worker_mocks_aux(args):
        return exec_worker_mocks(*args)


    def exec_worker_mocks(i, simset_mock, lcs, simoptfct, kwargs_optim, optset, tsrand, destpath):
        print("worker %i starting..." % i)
        time.sleep(i)
        sucess_dic = pycs3.sim.run.multirun(simset_mock, lcs, simoptfct, kwargs_optim=kwargs_optim,
                                           optset=optset, tsrand=tsrand, keepopt=True, destpath=destpath)
        return sucess_dic


    nworkers = 8
    kwargs = {} # if your optimiser, i.e. myopt.spl, takes argument you can pass them here.
    job_args = [(j, "sim1tsr10", lcs, myopt.spl, kwargs, "spl1", 10.0, "./") for j in range(nworkers)]
    p = Pool(nworkers)
    success_list_copies = p.map(exec_worker_copie_aux, job_args)

For a detailed example, you can check this `script <https://gitlab.com/cosmograil/PyCS3/-/blob/master/scripts/3c_optimise_copy_mocks.py>`_.

Analysing the measurement results
---------------------------------


We read the "runresults" pickle files created at the previous step, and turn them into plots.
This is very flexible, as you might want to plot and analyse many things.

To start, we have the function :py:func:`pycs3.sim.run.collect` that collects all the results from one directory::

	results = pycs3.sim.run.collect(directory="./for/example/sims_copies_opt_spl")

The resulting object ``results`` is an instance of the class :py:class:`pycs3.sim.run.RunResults`. If you want to perform your own analysis of the results, you could directly access the following attributes::

	print(results.labels) # A list of the QSO image names (defines the order of QSO images with which the following results are given)
	print(results.tsarray) # A 2D array with the measured time shifts. Shape is (number of sets, number of QSO images)
	print(results.truetsarray) # Idem, for the TRUE time shifts, in case of simulated data
	print(results.qs) # A 1D array with the "chi2" or dispersion values. Shape is (number of sets).

Note that these "tsarrays" contain time shifts, not time delays. To get time delays between images "A" and "B" (i.e., ``results.labels[0]`` and ``results.labels[1]``), you would have to compute the differences yourself::

	measured_delays = results.tsarray[:,1] - results.tsarray[:,0]
	print(measured_delays)


If you want to go straight to some more or less automatic plots showing the results, here is a typical example:

::

		
	copiesres = [
		pycs3.sim.run.collect("sims_copies_opt_spl1", "blue", "Spline"),
		pycs3.sim.run.collect("sims_copies_opt_regdiff1", "green", "Regdiff")
	]
	
	pycs3.sim.plot.hists(copiesres, r=30.0, nbins=100, dataout =True)
	
	
	simres = [
		pycs3.sim.run.collect("sims_sim1tsr10_opt_spl1", "blue", "Splines"),
		pycs3.sim.run.collect("sims_sim1tsr10_opt_regdiff1", "green", "Regdiff")
	]
	
	
	pycs3.sim.plot.hists(simres, r=30.0, nbins=100, dataout =True)
	
	pycs3.sim.plot.measvstrue(simres, r=10.0, nbins = 1, plotpoints=True, ploterrorbars=True, sidebyside=True, errorrange=8, binclip=False, binclipr=20.0, dataout =True)


The measured time delays and their associated uncertainties are stored in pkl files that can be later processed with the :py:mod:`pycs3.tdcomb` :

::


    group_list = [pycs3.tdcomb.comb.getresults(pycs3.tdcomb.comb.CScontainer("Free-knot Spline",
                                                result_file_delays='sims_copies_opt_regdiff1_delays.pkl',
                                                result_file_errorbars='sims_sim1tsr10_opt_regdiff1_errorbars.pkl',
                                                colour = 'blue')),

                  pycs3.tdcomb.comb.getresults(pycs3.tdcomb.comb.CScontainer("Regression Difference",
                                                result_file_delays='sims_copies_opt_spl1_delays.pkl',
                                                result_file_errorbars='sims_sim1tsr10_opt_spl1_errorbars.pkl',
                                                colour = 'red'))
                    ]

    pycs3.tdcomb.plot.delayplot(group_list, rplot=10, hidedetails=True,
                                  showbias=False, showran=False, showlegend=True, figsize=(15, 10), auto_radius=True,
                                  tick_step_auto=True)


.. image:: ../_static/tutorial/delays.png
	:align: center
	:width: 800
	

	
	
	
	


	
Tips and Tricks, code snippets
==============================


.. contents::


Which PyCS am I using ?
-----------------------

If you have several copies, from SVN or installed...
::

	import pycs3
	print(pycs3.__file__)
	


"Faking" PyCS3 delay or error measurements
------------------------------------------

This is useful if you want to include measurements from non-pycs techniques on the typical pycs3 plots
::

	# To fake a PyCS3 delay measurement :
	data = [{"label":"AB", "mean":-118.6, "med":0.0, "std":0.0}]

	# To fake a PyCS3 errorbar :
	"""
	data = [{
	"label":"AB",
	"sys":3.0,
	"ran":4.0,
	"tot":5.0,
	"bias":-3.0
	}]
	"""
	outname = "sims_copies_runresults_rathna_May07_delays.pkl"
	dc = pycs3.sim.plot.delaycontainer(data = data, name = "Difference-smoothing technique", plotcolour = "darkorange", objects=["A", "B"])
	pycs3.gen.util.writepickle(dc, outname)


Building scrolling plots for long curves
----------------------------------------

::

	# Animated plot for talk :

	startjd = 52900.0
	width = 1000.0
	endjd = 55800.0
	n = 1000
	
	for i in range(n):
		
		a = startjd + i* (endjd - width - startjd)/(n-1)
		b = a + width
		
		filename = "mov/%i.png" % (i)
		pycs3.gen.lc_func.display(lcs, nicefont=True, showdelays=False, showlegend=False, showdates=True, showgrid=True, magrange=(4.3, 0), jdrange=(a, b), filename=filename)
	

And then use ffmpeg (or any other similar tool) to turn this into a movie.

Tweaking magnitudes for individual seasons
------------------------------------------

For a lightcurve ``l``, ``l.mags`` is just a numpy array.
To *lower* the third season by 0.03 mags :
::
	
	seasons = pycs3.gen.sea.autofactory(l)
	l.mags[seasons[2].indices] += 0.03
	



Playing with custom properties
------------------------------

You can perfectly create your own properties. It's just a list of dicts ...
::
	
	for i in range(len(l)):
		l.properties[i]["my_new_prop"] = "brocoli"
		
	# To see what properties a curve has :
	print(l.longinfo())

"Common" properties are properties that all points of the curve have (this is usually the case). Only those "common" properties can be exported as columns in rdb files, for instance.


Splitting a curve by properties
-------------------------------

::
	
	def splitbyprop(l, prop = "telescope"):
		"""
		kills mask ...
		"""
		
		vals = sorted(list(set([l.properties[i][prop] for i in range(len(l))])))
		
		out = []
		for val in vals:
			lcp = l.copy()
			lcp.mask = np.array([l.properties[i][prop] == val for i in range(len(l))])
			lcp.cutmask()
			lcp.telescopename = val
			out.append(lcp)
			
		#pycs3.gen.mrg.colourise(out)
		return out




Correcting for flux sharing
---------------------------

March 2012, only implemented for the spline method. Simple code works well, but quick tests on simulated data (HE2149) show degeneracies.
Need complete tests on simulated data with a little flux sharing, to see if it reduces systematic error.

::

	# draw fake curves :
	flcs = pycs3.sim.draw.draw(lcs, spline, shotnoise="none", keepshifts=False)
	pycs3.sim.draw.shareflux(flcs[0], flcs[1], frac=0.02)
	pycs3.gen.lc_func.display(flcs)

	# then run pycs3.spl.topopt.opt_fine, it has the option "redistribfluxes"
Drawing mock curves
===================


The splines introduced in the previous section are used to draw simulated curves with known time delays. These mock curves can be made very similar to the real observations.

At the bottom of this page a single function to build an entire population of simulated curves is presented. This will be extensively used to empirically evaluate the accuracy of curve shifting methods. But before looking at this wrapper, we have to introduce the details of how these curves are drawn, step by step. If you are not interested in the details, you can skip this part and simply follow this `notebook <https://gitlab.com/cosmograil/PyCS3/-/blob/master/notebook/Uncertainties%20estimation.ipynb>`_


Drawing individual curves
-------------------------

The idea is the following : you provide some real light curves, and a spline usually obtained from a spline optimizer. Your curves might have some (optimized) microlensing representations, and they will typically be shifted so to match to the spline (i.e., your curves have some well defined time/mag/flux shifts set).

The functions of the module :py:mod:`pycs3.sim.draw` will use your curves' time/mag/flux shifts, as well as their sampling, errorbars, and microlensing, to draw mock curves *from the spline you give them*. This means that if we would **not** add noise to these mock curves, they would *perfectly match to the spline*, given all the shifts that your real curves had.

Here comes an example. The function to "manually" draw some mock curves has the funny name :py:func:`pycs3.sim.draw.draw`. It takes a list of light curves as well as a spline as input, and returns a corresponding list of mock curves. To use this option, and for practical reasons that will become clear later, we need to "save" the residuals (:py:func:`pycs3.sim.draw.saveresiduals`) of your lcs before drawing the mock curves.

::
	
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	# Some lightcurves, and a spline to play with.
	
	pycs3.sim.draw.saveresiduals(lcs, spline)
	mocklcs = pycs3.sim.draw.draw(lcs, spline, shotnoise=None)
	# We draw our first mock curves, here without adding any noise !
		
	pycs3.gen.lc_func.display(mocklcs, [spline])
	# As you see, by default those mocklcs "come with" the same shifts and microlensing as your lcs.

	# Let's manually remove these shifts from all our curves :
	for l in mocklcs:
		l.resetshifts() # also removes microlensing
		l.plotcolour = "black"
	for l in lcs:
		l.resetshifts()
	
	# And overplot the original and mock curves :
	pycs3.gen.lc_func.display(lcs + mocklcs, showdelays=False)


.. note:: The mock curves just created are brand new independent :py:class:`pycs3.gen.lc.LightCurve` objects. Everything you did so far with lightcurve objects applies also to these mock curves.

The last plot looks like this (the black points are the mock curves, as expected without noise in this case) :

.. image:: ../_static/tutorial/mock_no_noise.png
	:align: center
	:width: 800



Adding some noise
-----------------

Of course, for about any purpose, we want our mock curves to be noisy. The first trivial way to do this is to add some random "white" (i.e., independent) noise to each magnitude measurement.
This could be done by drawing random gaussian errors according to the errorbars of each point (option shotnoise="magerrs" below), or, to avoid explicitly using the errorbars, we could use the actual observed mismatch between your shifted lcs and the spline, which are the residuals we saved just before.

::
	
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	
	# So these lightcurves match to the spline
	pycs3.sim.draw.saveresiduals(lcs, spline)

	mocklcs = pycs3.sim.draw.draw(lcs, spline, shotnoise="mcres")
	# "mcres" adds some random gaussian noise to the mock curves,
	# using gaussian distributions whose sigma are the previously saved residuals.
	
	pycs3.gen.lc_func.display(mocklcs, [spline])


These new mock curves will now already look rather similar to your observed data.

But the whole point is that we *know* the "true" delays of these mock curves. In fact, the mock curves have an extra "secret" attribute (no need to remember, later functions will do all the calculations for you) :

::

	for l in mocklcs:
		print(l.truetimeshift)

... that stores what shifts where used to obtain those curves, and hence what the true delays between them are.



Choosing your own shifts
------------------------

Simply shift the curves (or modify their microlensing) *before* calling :py:func:`pycs3.sim.draw.draw` (but after having saved the residuals if you want to use them (``shotnoise = mcres`` or ``res``), otherwise these residuals will be crap or course, as the curves won't match to the spline anymore).




Randomizing the microlensing
----------------------------

The aim here is to randomly add some "fast" extrinsic variability on top of the existing microlensing splines.

For illustration purposes, let's start by doing this manually with the high level function :py:func:`pycs3.sim.twk.tweakml`. It takes as argument some lightcurve objects, and adds power-law "noise" to their microlensing, using under the hood the algorithm by Timmer and Koening 1995.
For this to work, the lightcurve objects must have spline microlensing (otherwise they simply won't be tweaked).
Once the function has run on them, they will still have spline microlensing objects, but with many many knots. So these microlensing objects are not meant to be be optimized -- they are just meant to be used as models to draw light curves from ! Of course you can display these lightcurves with tweaked ML.

To illustrate this, we can just tweak the ML of the "observed" data::

	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	
	# I assume here that at least one of your lcs has some spline ML.
	
	pycs3.sim.twk.tweakml(lcs, spline, beta=-2.0, sigma=0.05, fmin=1/500.0)

	# And plot this, to see the tweaked ML :
	pycs3.gen.lc_func.display(lcs)
	


.. image:: ../_static/tutorial/tweakml.png
	:align: center
	:width: 800

.. note:: In fact, the microlensing curves are noisier on small scales then suggested by the above image, but the display function does not sample the microlensing objects finely enough. This is especially true if you interactively zoom in.


You can experiment a little with different beta, sigma, fmin, fmax, that control the power law noise that will be added to the microlensing.
Also you can try setting the option psplot=True of tweakml, it will show you power spectra.

``beta = -2.0`` corresponds to a random walk !


As you guess, you could use :py:func:`pycs3.sim.draw.draw` to draw light curves from these tweaked ones.

 
So this was a nice example to get the idea, but in fact, you don't want to tweak the ML of your lcs *once*, but you want to draw mock curves with always newly tweaked ML.

That's why instead of explicitly calling your mytweakml function, we will just pass this function as an argument to :py:func:`pycs3.sim.draw.draw`, and the latter will take care of tweaking the ML itself.


Here is a (new) example :

::
	
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	
	# Maybe you need to add some spline ML to curves that don't have it yet :
	pycs3.gen.splml.addtolc(lcs[0])
	
	# We define our own tweakml function (you can also do this in myopt.py ...)
	def mytweakml(lcs, spline):
		return pycs3.sim.twk.tweakml(lcs, spline, beta=-2.0, sigma=0.05, fmin=1/500.0, fmax=None, psplot=False)

	# And directly draw mock curves :
	mocklcs = pycs3.sim.draw.draw(lcs, spline, shotnoise="none", tweakml = mytweakml)

	pycs3.gen.lc_func.display(mocklcs, [spline])
	
	# These mocklcs are drawn without any "shotnoise", all the noise comes from tweakml.

Alternatively, you can use the new functionality of PyCS3 that is using the directly the power spectrum of the residuals to inject correlated noise at the same frequencies than measured in the real data. You can define this function as follow :

::

    def tweakml_PS(lcs, spline):
        return pycs3.sim.twk.tweakml_PS(lcs, spline, B=1.0, f_min=1 / 300.0, psplot=False, verbose=False,
                          interpolation='linear', A_correction=1.0)

and use it similarly to :py:func:`pycs3.sim.twk.tweakml`. Note that this do not garantee that your mock light curves will be similar as the real data. You might need to adjust the spectral window (controlled by the parameter `B`) and the amplitude of the power spectrum with the parameter `A_correction`. This can be done automatically with this `script <https://gitlab.com/cosmograil/PyCS3/-/blob/master/scripts/3a_generate_tweakml.py>`_ and the :py:mod:`pycs3.pipe.optimiser` module.


.. note:: Instead of providing a single "mytweakml" function to draw, you can also provide a *list* of mytweakml-like functions, each item of this list corresponding to a light curve in your lcs. This way you can individually adapt the tweakml to the noise properties in each curve.
	Same is true for :py:func:`pycs3.sim.draw.multidraw` described below !
	
	::
		
		# Define different tweakml functions, and then (example) : 
		mocklcs = pycs3.sim.draw.draw(lcs, spline, tweakml=[Atweakml, othertweakml, othertweakml, othertweakml], shotnoise="none")

		


To generate adequate simulations, we now want to adjust tweakml (and shotnoise) so to get the same kind of residuals between the spline and the real lcs and between the spline and the mocklcs.
We compute those residuals in the next section.

Checking spline residuals
-------------------------

Here are some functions to take a curve, take a spline, "subtract" the spline from the curve, and analyse/look at the scatter of the residuals :

* :py:func:`pycs3.gen.stat.subtract`
* :py:func:`pycs3.gen.stat.mapresistats`
* :py:func:`pycs3.gen.stat.plotresiduals`


Here is how to get a plot of the residuals :

::
	
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	
	rls = pycs3.gen.stat.subtract(lcs, spline) # This simply subtracts the spline from the datapoints.
	# rls is a list of new lightcurve objects, corresponding to "lcs - spline".
	# You could display it as usual.
	
	# Stats about the residuals :
	print pycs3.gen.stat.mapresistats(rls)
	
	# A special function to plot residuals :
	pycs3.gen.stat.plotresiduals([rls])


Putting this together with some mocklcs:

::
	
	pycs3.gen.splml.addtolc(lcs[1]) # So that all curves have some SplineML !

	def mytweakml(lcs):
		return pycs3.sim.twk.tweakml(lcs, beta=-0.5, sigma=1.5, fmin=1/500.0, fmax=None, psplot=False)
	
	mocklcs = pycs3.sim.draw.draw(lcs, spline, tweakml=mytweakml, shotnoise="none", keeptweakedml=False)
	
	for l in mocklcs:
		l.plotcolour = "black"
	
	rmocklcs = pycs3.gen.stat.subtract(mocklcs, spline) # Same as for the real data.
	# Note that it would be better to fit a new spline to the mocklcs, using the old one is a shortcut ...


	pycs3.gen.stat.plotresiduals([rlcs, rmocklcs])
	# Yes, this function takes lists of corresponding lightcurve-lists, exactly for this purpose.


The resulting plot (coloured points are the real curve, black points are a mock curves) :

.. image:: ../_static/tutorial/resi_tweakml.png
	:align: center
	:width: 800



	


Building sets of mock curves
----------------------------

This is done with one single function, the topmost wrapper, called :py:func:`pycs3.sim.draw.multidraw`. It uses :py:func:`pycs3.sim.draw.draw`, and stores the drawn curves in pickle files. The same function is also used to simply make a set that contains plain copies of your original curves (I agree, this seems stupid, but hey its flexible).

.. note:: In any case, the curves returned by :py:func:`pycs3.sim.draw.multidraw` are **raw observations** : they have no shifts, no ML. Just datapoints !

These mock curves will later be analysed by :py:func:`pycs3.sim.run.multirun`.

.. note:: The files I save are just pickles of lists of "lcs". You are welcome to read such a pickle and display it.
	

Define a function to tweak the ml, as above (for instance in ``myopt.py``) :

::

	def tweakml(lcs):
	    return pycs3.sim.twk.tweakml(lcs, beta=-2.0, sigma=0.03, fmin=1/300.0, fmax=None, psplot=False)

.. warning:: You will probably want to add some spline microlensing to **all** your lcs before calling ``multidraw`` or ``draw``, as they will tweak the microlensing only of those curves that have microlensing !

::
	
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	pycs3.sim.draw.saveresiduals(lcs, spline)
	
	pycs3.gen.splml.addtolc(lcs[0]) # So that all curves have some SplineML !
	
	#pycs3.gen.lc_func.display(lcs, [spline])
	
	#pycs3.sim.draw.multidraw(lcs, onlycopy=True, n=20, npkl=10, simset="copies")
	
	#pycs3.sim.draw.multidraw(lcs, spline, onlycopy=False, n=20, npkl=30, simset="sim1tsr5", shotnoise="mcres", shotnoisefrac=1.0, truetsr=5.0, tweakml=myopt.tweakml, tweakspl=None)




Displaying some curves drawn with multidraw
-------------------------------------------

Just to show that the structure of those pkl files is very easy

::

	# We read in the original data, to overplot :
	lcs = pycs3.gen.util.readpickle("merged.pkl")
	for l in lcs:
		l.resetshifts()
		l.plotcolour = "black"
	
	# Reading in a random pickle file :
	mocklcslist = pycs3.gen.util.readpickle("sims_sim1tsr5/2_1334738572.78151.pkl")
	pycs3.gen.lc_func.display(mocklcslist[0] + lcs, showdelays=False)


.. image:: ../_static/tutorial/mockJ1001.png
	:align: center
	:width: 800







First steps with microlensing
=============================

PyCS3 offers two families of microlensing representations : polynoms, and splines. Ok of course splines are piecewise polynoms as well, but we mean free-knot splines that are flexible and stable enough to cover decade long curves.

Physically, a microlensing representation can only be "relative", i.e., represent the mismatch *between* curves. Indeed we do not have access to the absolute microlensing magnification that affects each QSO image.

Nevertheless, in ``pycs3`` a microlensing representation is "attached" to a lightcurve, as if it would be absolute. This allows you to control individually for each curve what kind of "flexibility" you want to allow.
If you analyse a quad, you might attach such microlensings to 3 of the 4 curves. This would mean that the microlensings are relative with respect to the curve without microlensing. But you can as well put microlensing on all 4 curves. The microlensing representations will then be degenerate, but this will not prevent you from getting the right delays.


Adding ML to your lightcurve
----------------------------

Before running any optimization, you will have to choose and add (aka "attach") the microlensing representation you want for each curve.

Two very easy functions let you add microlensing to a curve with one single line :py:func:`pycs3.gen.splml.addtolc` (splines) and :py:func:`pycs3.gen.polyml.addtolc` (polynoms):
::
		
	# Adding a spline, with one knot every 200 days (initial position) :	
	pycs3.gen.splml.addtolc(l, knotstep=200)
	
	# Or, if you prefer polynoms, this adds a straight line to each season of l :
	pycs3.gen.polyml.addtolc(l, nparams=2)
	
	# The latter function "automatically" recognizes seasons. But with this same function,
	# you can also add one single polynom over your full multi-season curve.
	# To do this, simply increase the autoseasongap.
	pycs3.gen.polyml.addtolc(l, nparams=3, autoseasonsgap=1000.0)
	# So this would be a parabola over all seasons.
	
	
Note that microlensing can also be configured without these functions, "by hand", if you have special wishes.
This was just the easy and automatic way.


So now that the microlensing is attached, it will be represented on the plots by continuous coloured lines.

::
	
	pycs3.gen.lc_func.display([l])
	
	
As you see, the microlensig is there, but completely "flat". Its initial coefficients are set so to not affect the curve.
	

Fun with ML
-----------

Mix it !

.. image:: ../_static/tutorial/fun_with_ml.png
	:align: center



Curve shifting, finally
=======================


PyCS3 makes it very easy to compare different point estimators. You can also add your own, without much code integration or inheritance. You would simply have to provide a python function that shifts instances of the lightcurve class.
High level functions are provided to run these point estimators on both real observations and the synthetic curves. This can be done in parallel on several CPUs, again using a very low tech approach.

This section describes how to do a first optimization of the time shifts (and of course the microlensing) so that the lightcurves of different QSO images "match".

Getting accurate time delay estimations (including the error analysis) is described in the next section.

We assume here that you have some curves (e.g. two, or four), nicely prepared according to the previous steps of this tutorial. Save these nicely merged curves into a pkl, as described earlier.

PyCS3 now support two different optimizers. *Optimizers* are functions that take a list of lightcurves as only argument, and shift these curves (potentially also adjusting the microlensing) so that they match.

.. warning:: The dispersion optimizer is no longer included in PyCS3. This is because the spline and regdiff optimizer outperforms the dispersion optimizer in term of accuracy and are, in addition, much faster. For these reasons, we decided to stop the development of the dispersion technique.

The spline optimizer
--------------------

The idea here is to fit one **single** spline (representing the "intrinsic" variation of the QSO) to all your curves, shifting the latter so that the chi2 between all the datapoints and this single spline is minimal.

This is a very parametric problem, and not a trivial task. The optimizer has to "simultaneously" adjust at least :

* the time/magnitude/flux -shifts of the curves
* the microlensing representation of the curves (polynom coefficients, or spline coeffs and knots)
* the intrinsic spline (both coefficients and knot positions)

In this first approach, we won't describe the internal details. In fact, the spline optimizer works in a very similar way than the dispersion optimizer described above : it shifts your curves, and adjusts their microlensing representations. But of course, it also involves one new object : the spline (that you will want to display on top of your lightcurves) !
Try it :

::

	lcs = pycs3.gen.util.readpickle("merged.pkl")
	
	# We might need some microlensing representations for some of the curves :
	pycs3.gen.splml.addtolc(lcs[1], knotstep=300) # Yep, let's try spline microlensing !
	
	# And now the optimizer. Note that it returns the spline object !
	spline = pycs3.spl.topopt.opt_rough(lcs, nit=5, knotstep=150)
	

To show this spline, we make use of :py:func:`pycs3.gen.lc_func.display` that you've met many times before. Indeed, this function can also display an arbitrary number of splines !

As usual, it's first argument is simply a list of lightcurve objects (or an empty list, if you don't want to show any lightcurves at all). But you can also specify as optional second argument a **list of spline objects** (hence the ``[ ]``) :

::

	pycs3.gen.lc_func.display(lcs, [spline])


Don't expect a perfect fit at this point -- this was just a first demo of the principles of the optimizer.
This particular spline is shown in black, with vertical ticks indicating its knot positions.

The spline optimizer seen above takes a few optional arguments. Some words about the arguments seen at this point :

* ``nit`` is a number of iterations, it's fine to leave it at 5 unless you dig into the details of these optimizers.
* ``knotstep`` sets the initial "spacing" (in days) of the knots. These knots will then move around, so the spacing will change... but the number of knots will not !

.. warning:: A *lower* knotstep corresponds to *more* knots !

When adding the spline microlensing to the curve, we specified ``knotstep=300`` : this is the same parameter, but for the knots of the microlensing. So choose a lower microlensing-``knotstep`` to get a more flexible microlensing.

The above example used the "rough" optimizer. This one is not made to get accurate time delays, but to roughly (and quickly) adjust the shifts and the microlensing so that the curves match somehow. Hence, for this *rough* part, leave a relatively high ``knotstep``.

Directly after this rough optimization, add a call to a finer optimizer :

::

	spline = pycs3.spl.topopt.opt_fine(lcs, nit = 5, knotstep=100)
	

This optimizer will build a new spline from scratch (and return it), using a (usually finer) ``knotstep`` of your choice. Add this line just after the call to opt_rough, and play with the knotstep (e.g. 50) to see the effect. Also note that the knots are now effectively moving (the opt_rough didn't move them).


.. image:: ../_static/tutorial/spline.png
	:align: center
	:width: 800


It's now a good idea to add these optimizers to your ``myopt.py`` file, directly concatenating them ! This allows you to build a custom optimizer for your particular light curve. Here is an example (you should probably update the knotsteps, depending on the curves you want to process) :

::

	def spl(lcs):
		spline = pycs3.spl.topopt.opt_rough(lcs, nit=5, knotstep=40)
		spline = pycs3.spl.topopt.opt_fine(lcs, nit=10, knotstep=30)

		return spline # Do not forget to return the spline !


You can now use ``myopt.spl(lcs)`` to return a spline, that you might want to "catch" by writing

::

	spline = myopt.spl(lcs)
	
As usual, after such an optimization, it might be convenient to save the shifted curves and in this case also the spline into a pickle file, so that you can work on them without rerunning the optimization. Tip : save both the curves and the spline into the same pickle file ! 

::
	
	pycs3.gen.util.writepickle((lcs, spline), "optspl.pkl")
		
	# ...
		
	(lcs, spline) = pycs3.gen.util.readpickle("optspl.pkl")
	pycs3.gen.lc_func.display(lcs, [spline])
		



To learn more about the optional arguments of the spline optimizers, see the doc of :py:func:`pycs3.spl.topopt.opt_rough` and :py:func:`pycs3.spl.topopt.opt_fine`.

These spline optimizers also work with polynomial microlensing. You can mix the microlensing representations at will.

.. note:: Formally, the linear optimization of splines requires data points that are not only sorted, but also *strictly* increasing in jds : it cannot deal with lightcurves that have several data points taken at exactly the same epoch (which may happen as we shift the curves in time). This issue is automatically adressed by the class :py:class:`pycs3.gen.datapoints.DataPoints`. As a user you don't have to worry about this in principle.


The regdiff optimizer
---------------------
This is de facto the easiest method to use, as it does not involve any explicit microlensing representation.

The idea is to shift the light curves so to minimize the variability of their "differences". To compute these difference curves, we need a regression, and in particular we use Gaussian Process Regression as provided by the ``scikit-learn`` module.

.. note:: Therefore, to use the regdiff optimizer, you will have to **install** ``scikit-learn`` first.
	
	Here is the website : `https://scikit-learn.org/stable/ <https://scikit-learn.org/stable/>`_



See `Tewes et al. 2013 <http://dx.doi.org/10.1051/0004-6361/201220123>`_ for a more detailed description of the idea. In practice, as for the spline optimizer, there is a simple top-level wrapper function, that you can add to your ``myopt.py`` :

::

	# we need to call the regdiff module explicitely
	import pycs3.regdiff

	def regdiff(lcs):
		return pycs3.regdiff.multiopt.opt_ts(lcs, pd=5, verbose=True)


But before blindly using the above optimizer, it is a good idea to test by yourself if the the Gaussian Process Regression (GPR) performs well on your lightcurve. The regressions are represented by "regularly sampled light curve" objects, implemented by the class :py:class:`pycs3.regdiff.rslc.Rslc`.
It is easy to perform a regression "manually", i.e. to obtain such a regularly sampled light curve starting from a usual light curve. The function that performs this GPR is :py:func:`pycs3.regdiff.rslc.factory`, and you could for instance directly apply this directly to all your light curves :

::

	myrslcs = [pycs3.regdiff.rslc.factory(l, pd=2) for l in lcs]
	# As this can take a minute, you might want to save the results :
	pycs3.gen.util.writepickle(myrslcs, "myrslcs.pkl")


The parameter ``pd`` is a point density of the regression. Usually this is set to 2 (corresponding to one point every 0.5 days). Less points will give you a faster regression.

You can display these ``Rslc`` objects with the usual display function, simply by putting them in the second argument list, as you would do for spline objects.

::
	
	myrslcs = pycs3.gen.util.readpickle("myrslcs.pkl")
	pycs3.gen.lc_func.display(lcs, myrslcs)


.. note:: These ``Rslc`` have some attributes very similar to the usual ``LightCurve`` objects, like ``jds``, ``mags``, ``magerrs``, ``plotcolour``. To shift an ``Rslc`` in time, use ``myrslc.shifttime(12.3)``. To perform other operations, directly modify the attributes, for instance : ``myrslc.mags += 3.0``.


The reason why we want these finely sampled light curves is that we can easily subtract them from each other to get difference curves. This operation is implemented by :py:func:`pycs3.regdiff.rslc.subtract`.

::
	
	diffrslc = pycs3.regdiff.rslc.subtract(myrslcs[0], myrslcs[1])
	# This diffrslc is the difference curve, and its again a rslc object.
	
	# Hence you can display the difference easily by putting it in the list, for instance :
	pycs3.gen.lc_func.display([], myrslcs + [diffrslc])



Finally, the *WAV* of any ``Rslc`` can be computed by calling the method :py:meth:`pycs3.regdiff.rslc.rslc.wtv`.



.. _matchtels:

Matching telescopes
===================

Lightcurves from different telescopes usually need to be empirically "matched" in terms of (at least) magnitude offset. Even if you try to perform relative photometry with respect to the same reference stars, a small shift in magnitude is typically needed, as the detectors + filters of the telescopes are different. Often a shift in **flux** is also needed.

A note about flux shifts
------------------------


In PyCS3, lightcurves are encoded in terms of magnitudes, and thus also shown in this way on plots. Thus, multiplicative macrolensing amplification (i.e., the flux ratio between quasar images), but also stationary microlensing amplifications correspond to a simple "vertical" shift of the lightcurves.

But aside of the attribute magshift, lightcurve objects also have an attribute fluxshift. **Note that a "shift" in flux does not correspond to a vertical shift of the curve plotted on a magnitude scale !** A shift in flux *deforms* the curves. That's why the flux shift represents a welcome parameter when matching lightcurves from different telescopes. The physical origin of a shift in flux can be any additive contamination of the QSO image (like a constant foreground star or some structure of the lens galaxy). As this contamination do not necessarily have the same colour as the QSO images, they can yield flux shifts not only between the QSO images, but also between curves of the same QSO image observed with different telescopes !

If your lightcurves comes from cosmouline/lcmanip, the "magnitudes" are simply related to the instrumental flux (in electrons) via :math:`m = -2.5 \log(f)`, i.e. :math:`f = 10^{(-0.4 \cdot m)}` : there is no zeropoint. Thus you can easily calculate yourself the fluxes of your curves.

Play a bit around with these fluxshifts (take a red curve, make a copy of it and set its plotcolour to blue, then set some fluxshifts, and plot both ...).


.. note:: Keep an eye on the values of the flux shifts, when optimizing them with the functions described below. There is probably no physical reason for a fluxshift to be larger than the actual measured flux of the source ...


Matching lightcurves from different telescopes
----------------------------------------------

Here we describe how to empirically adjust the magshift and fluxshift of different curves of the *same* QSO image, so that they match.

There's a module for this : :py:mod:`pycs3.gen.mrg`.

Of course, this makes only sense if your curves show some overlap !

First of all, we need to define a *matching criterion*. A simple yet effective choice is a dispersion measure. This can be done with the :py:func:`pycs3.gen.lc_func.linintnp` function

::

	pycs3.gen.lc_func.linintnp(lc1, lc2, interpdist = 30.0)


The function :py:func:`pycs3.gen.mrg.matchtels` will now optimize the magnitude shifts and flux shifts of some curves so that the match, given this criterion. It takes 3 mandatory arguments : a list of reference lightcurves (they will not be modified), a corresponding list of lightcurves to be adjusted so that they match to the reference curves, and the dispersionmethod defined above. Note that :py:func:`pycs3.gen.mrg.matchtels` will attribute **one same magshift** and (if asked) **different individual fluxshifts** to all the curves ! If you want to adjust an individual magshift for each curve, simply call ``matchtels`` on dedicated lightcurve lists.

::

	# Let's assume we have a quad observed with telescopes Euler and Mercator.
	# This means we probably already have two corresponding lightcurves lists like e.g.
	# eulerlcs = [a_euler, b_euler, c_euler, d_euler]
	# mercatorlcs = [a_mercator, b_mercator, c_mercator, d_mercator]

	pycs3.gen.mrg.matchtels(eulerlcs, mercatorlcs, pycs3.gen.lc_func.linintnp, fluxshifts=True)

	# Let's see if it worked :
	
	pycs3.gen.lc_func.display(eulerlcs + mercatorlcs, showdates=True, showdelays=False)


Of course you can also set the fluxshifts and magshifts by hand ...
If you are happy with the match, merge the lightcurves. Instead of doing this one by one, there is a function to do this operating directly on lightcurve lists :py:func:`pycs3.gen.mrg.merge`

::
	
	lcs = pycs3.gen.mrg.merge([eulerlcs, mercatorlcs])
	pycs3.gen.mrg.colourise(lcs)
	
	pycs3.gen.lc_func.display(lcs, showdates=True, showdelays=False)
	
If you have another telescope, you can at this point call ``matchtels`` again to match this third telescope on the merged curves ``lcs`` that you have just obtained.




	
	
Generalities
============


The documentation
-----------------

... consists of several parts:

* Runnable `Example notebooks <https://gitlab.com/cosmograil/PyCS3/-/tree/master/notebook>`_, for a quick introduction on how to measure time delays with PyCS3. You might want to start from this if you just want a quick overlook of PyCS3's key functionality.
* A tutorial that covers the main objects/functions in details. (You are here !)
* A `README <https://gitlab.com/cosmograil/PyCS3/-/blob/master/scripts/README.md>`_, for the documentation of the automated pipeline described in `Millon et al. (2020) <https://arxiv.org/abs/2002.05736>`_.
* The autogenerated full API documentation: :doc:`../apidoc/pycs3`


Importing PyCS3
---------------

If PyCS3 is installed, simply use ``import pycs3`` at the beginning of your script.
A minimal script (just to check that it works)::
	
	import pycs3.gen.lc
	lca = pycs3.gen.lc.LightCurve(object="A", plotcolour="red")
	lcb = pycs3.gen.lc.LightCurve(object="B", plotcolour="blue")
	lcb.shifttime(1.0)
	pycs3.gen.lc_func.display([lca, lcb])

Interacting with PyCS3
----------------------

PyCS3 is now using a logger for its output message. Make sure you include these lines at the beginning of all your scripts if you want to print the info messages in the terminal ::

    import logging
    loggerformat='%(message)s'
    logging.basicConfig(format=loggerformat,level=logging.INFO)


Nomenclature
------------

* For simplicity we call "microlensing" (abreviated ML) any form of extrinsic variability.


About dates : JD, MJD, HJD, ...
-------------------------------

To represent observing epochs, we use the same convention as ``cosmouline``, that is a "Modified Heliocentric Julian Date".

The "modification" is defined by : **mjd = jd - 2400000.5** . In the same way, mhjd = hjd - 2400000.5 .


.. note:: Despite this, in the code we simply use ``jd`` or ``jds`` (it's a plural s) to refer to mhjd. For instance, the numpy array that stores the observing epochs of a light curve ``l`` is ``l.jds``.

Note that in principle, it is perfectly OK to use PyCS3 with any time axis, as long as the latter is linear, and can be stored as floats. Any advanced interpretations of the observing epochs are only used for uncritical stuff like plots and exports.


Tutorial
========

.. toctree::
	:maxdepth: 5
	:numbered:
	
	generalities
	basic_lc
	merging_telescopes
	microlensing
	shifting
	drawing
	makehists
	tips_and_tricks
	
The lightcurve object
=====================

The module :py:mod:`pycs3.gen.lc` defines the class :py:class:`~pycs3.gen.lc.LightCurve` that constitutes the backbone of PyCS3.
First step is to learn about this fundamental object, and the functions to manipulate lightcurve instances. The module :py:mod:`pycs3.gen.lc_func` is now gathering all the higher level functions to manipulate lists of :py:class:`~pycs3.gen.lc.LightCurve`



Importing lightcurves
---------------------

You probably have some light curves in the form of simple text/csv/rdb files. If not, get the trial curves from :download:`here <trialcurves.txt>` or get some real COSMOGRAIL curves from `here <https://obswww.unige.ch/~millon/d3cs/COSMOGRAIL_public/code.php>`_


Importing from cosmouline/lcmanip or similar rdb files
""""""""""""""""""""""""""""""""""""""""""""""""""""""

Cosmouline and lcmanip write their lightcurves as rdb files (i.e., lots of tab-separated collumns, and a short header).
By *importing* such data into PyCS3, I mean *reading such rdb files* and turning them into the mentionned lightcurve objects.
Here are some lines of a cosmouline/lcmanip output file (for a double lens) :

::

	mhjd	datetime	telescope	setname	nbimg	fwhm	ellipticity	airmass	relskylevel	normcoeff	flag	mag_A	magerr_A_1	magerr_A_2	magerr_A_3	magerr_A_4	mag_B	magerr_B_1	magerr_B_2	magerr_B_3	magerr_B_4
	====	========	=========	=======	=====	====	===========	=======	===========	=========	====	=====	==========	==========	==========	==========	=====	==========	==========	==========	==========
	54061.171611	2006-11-22T04:07:07	EulerC2	1	5	1.782	0.236	1.081	 507.4	1.018	True	-12.594276	0.0041	0.0100	0.0376	0.0174	-11.372701	0.0093	0.0130	0.0897	0.0266
	54065.165098	2006-11-26T03:57:44	EulerC2	1	5	1.637	0.207	1.086	 363.0	0.830	True	-12.613698	0.0040	0.0111	0.0210	0.0066	-11.333895	0.0092	0.0139	0.0706	0.0130
	54140.086307	2007-02-09T02:04:17	EulerC2	1	5	1.792	0.116	1.882	 558.4	0.873	True	-12.623750	0.0042	0.0073	0.0220	0.0056	-11.298380	0.0109	0.0124	0.1088	0.0387
	54161.021235	2007-03-02T00:30:35	EulerC2	1	4	2.117	0.077	1.792	2845.2	0.863	True	-12.617584	0.0075	0.0112	0.0210	0.0062	-11.202763	0.0254	0.0267	0.1433	0.0747
	54169.004465	2007-03-10T00:06:26	EulerC2	1	2	1.701	0.197	1.879	2147.2	0.868	False	-12.560599	0.0059	0.0064	0.0708	0.0259	-11.266728	0.0169	0.0171	0.2567	0.0902
	54316.327864	2007-08-04T07:52:07	EulerC2	1	5	1.600	0.186	1.405	1637.6	1.242	True	-12.533420	0.0052	0.0224	0.0518	0.0229	-11.134323	0.0153	0.0266	0.0779	0.0435


Let's assume we have such an rdb file called "Mercator.rdb", made with lcmanip, and containing data from a quad, i.e., 4 images : A, B, C and D.
This will give us 4 lightcurve objects; here is how :

::

	filename = "data/Mercator.rdb"
	lcs = [
		pycs3.gen.lc_func.rdbimport(filename, object="A", magcolname="mag_A", magerrcolname="magerr_A_4", telescopename="Mercator", plotcolour="red", mhjdcolname="mhjd", flagcolname = "flag", propertycolnames = ["fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff"], verbose = True),
		pycs3.gen.lc_func.rdbimport(filename, object="B", magcolname="mag_B", magerrcolname="magerr_B_4", telescopename="Mercator", plotcolour="red", mhjdcolname="mhjd", flagcolname = "flag", propertycolnames = ["fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff"], verbose = True),
		pycs3.gen.lc_func.rdbimport(filename, object="C", magcolname="mag_C", magerrcolname="magerr_C_4", telescopename="Mercator", plotcolour="red", mhjdcolname="mhjd", flagcolname = "flag", propertycolnames = ["fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff"], verbose = True),
		pycs3.gen.lc_func.rdbimport(filename, object="D", magcolname="mag_D", magerrcolname="magerr_D_4", telescopename="Mercator", plotcolour="red", mhjdcolname="mhjd", flagcolname = "flag", propertycolnames = ["fwhm", "ellipticity", "airmass", "relskylevel", "normcoeff"], verbose = True)
	]
	

So this uses the function :py:func:`pycs3.gen.lc_func.rdbimport` (click on such links to jump to detailed documentation of the function).
In the example above I've put all the options to illustrate them. Not all these options are required !
The option ``propertycolnames`` allows to read in extra columns of data, that can for instance be used to make colour-coded plots, or custom stuff.
If you omit this option, the function will by default import all the fields specified in the above example, if these are available in the rdb file. 
Note that you are free to choose whatever "object" or "telescopename" you want, these attributes are not related to the content of the rdb file.

Here is a more concise example that you could use as "default" to import lightcurves from lcmanip:

::
	
	mercatorlcs = [pycs3.gen.lc_func.rdbimport("data/2013-04-02_config_J1001_mercator_lcmanip.rdb",
		magcolname="mag_%s"%e, magerrcolname="magerr_%s_5"%e, flagcolname = "flag",
		telescopename="Mercator", object=e, plotcolour="blue") for e in ["A", "B"]]





A note about object instance names
""""""""""""""""""""""""""""""""""

In the following tutorial (and in the code, actually), I often use the short name ``l`` to designate a given lightcurve object. If we have several lightcurve objects in play, and if we need to give them explicit names for the purpose of the tutorial, we will often call them ``lca``, ``lcb``, ... as they typically represent curves from QSO image A, image B, etc.
A list of lightcurves (as we have just imported above) is usually named ``lcs``. It will be very common to do manipulations like ::

	for l in lcs:
		l.do_something()
		

.. warning:: Note that if possible you should **avoid** to give explicit names to your lightcurve variables  (like ``lca``, ``lcb``) and refer to them using such names all across your code. That's what lists are for ! Just read on :-)


Importing from simple tab-separated files
"""""""""""""""""""""""""""""""""""""""""

The function :py:func:`pycs3.gen.lc_func.flexibleimport` does essentially the same as :py:func:`pycs3.gen.lc_func.rdbimport`, but for simpler input files that do not have headers :

::
	
	l = pycs3.gen.lc_func.flexibleimport("curve.dat", , jdcol=1, magcol=2, errcol=3, startline=1, flagcol=None, propertycols=None,
                   telescopename="Unknown", object="Unknown", plotcolour="crimson", verbose=True, absmagerrs=False)
	

The shown arguments values are the defaults.


Set some plotcolour !
"""""""""""""""""""""

After such an import, and before plotting the curves to see if it went well, you might want to

::
	
	pycs3.gen.mrg.colourise(lcs)
	
	
them. This function sets the attribute ``plotcolour`` of each lightcurve to a different colour.
It is important to understand that you can do such operations by hand at any time :

::
	
	l.plotcolour = "brown"

Colours are simply matplotlib colours, so you can use whatever matplotlib accepts.




Plotting lightcurves
--------------------

It's time to see these colours.
Many functions of pycs3 work with lists of lightcurves, instead of individual lightcurves. These lists usually contain simply one curve for every QSO image. But be aware that all this works in the same way if you use lists that contain curves from different telescopes, or "identical" curves with different settings, or simulated curves, etc.

As mentionned, in the tutorials we will usually call such lists ``lcs``; *s* is a plural s.
The single most important function that uses such a list of curves as argument is the function that displays them :

::

	lcs = [lca, lcb, lcc, lcd] # So that's a list 

	pycs3.gen.lc_func.display(lcs)
	

This function has lots of options, it can be used for many tasks. As we will see in this tutorial, the same function is used to plot microlensing, splines, etc. Here is a link to the full documentation for this function: :py:func:`pycs3.gen.lc_func.display` (have a look).
For now, just as an example, try these options ::

	pycs3.gen.lc.display(lcs, title=r"$\mathrm{SDSS\,J1234-5678}$", nicefont=True, showlogo=True)
	# The option nicefont is your friend if you like serif fonts
	# (I don't, except for the title in LaTeX which is always in serif)



.. image:: ../_static/tutorial/display.png
	:align: center



Manually shifting lightcurves in time, magnitude, and flux
----------------------------------------------------------


We have 3 elementary methods to do this :

* :py:meth:`pycs3.gen.lc.lightcurve.shifttime`
* :py:meth:`pycs3.gen.lc.lightcurve.shiftmag`
* :py:meth:`pycs3.gen.lc.lightcurve.shiftflux`

::
	
	l.shifttime(5.0) # Shifts the curve by +5 days with respect to its current shift.
	l.shiftmag(-0.2) # Shifts the curve by -0.2 mags (i.e., it gets brighter) with respect to its current shift.
	
	l.shiftflux(2000.0) # "Shifts" the curve by +2000.0 electrons with respect to its current shift.
	# Note that on a magnitude plot, that's actually not a shift, it deforms the curve !


A lightcurve object is always "aware" of its shifts. These shifts don't get *applied* to the data (as long as you don't ask for it). They just set attributes of the lightcurve, telling them by how much they are shifted. The actual data is not modified. It is also perfectly ok to directly tweak the attributes :

::
	
	l.timeshift = 0.0 # "Resets" the curve
	

In a nutshell, we could now see a curve shifting method as a python function that sets these shifts for you, so to minimize a given criteria (for instance a dispersion measure) between curves. More on this later.


Displaying info about lightcurves
---------------------------------

::
	
	print(l) # Short oneliner; corresponds to str(l), that is also used in plot legends, etc.

The ouput might come with a paranthesis containing 3 numbers, like for instance ``[Mercator/A](10.000,-0.500,1200)``. This would mean that the curve is shifted by 10 days in time, -0.5 mag in magnitude, and 1200 counts in flux.

::
	
	print(l.longinfo())


Gives you a wider picture. Try it !

To display time delays between some curves, try this :

::
	
	print(pycs3.gen.lc.getnicetimedelays(lcs, separator = " | "))
	print(pycs3.gen.lc.getnicetimedelays(lcs, separator = " | ", sorted = True)) # Sorts according to object names
	

About "properties"
------------------

You saw how to import them, you saw how to use them in plots. 
Properties are very flexible. You can access/modify them from within your scripts, to store just about anything you want.
Properties are stored as entries of dictionnaries in a list as long as your curve (i.e., one dict per data point).

::
	
	print (l.properties)   # That's a long list of dicts.
	
	l.properties[0]["fwhm"] = "10.0"   # Tweak fwhm of first point
	
	for point in l.properties:
		point["w"] =  ... # Add your own properties !


.. note:: To keep all the import/export functionality, store your custom properties as strings. Indeed all the stuff like "fwhm" and "ellipticity" is stored as strings as well.

.. warning:: Some functions of ``pycs3`` might get significantly slower when you use properties. For instance stuff that requires merging of curves.




Cutting seasons
---------------

The module :py:mod:`pycs3.gen.sea` contains a class and functions to handle seasons.
You can define seasons "by hand", but usually for cosmograil curves the default automatic season detection is fine.

The concept of seaons can be important when defining microlensing representations.
Seasons are also handy to cut curves. There is a very easy function to do just this. In the following example we want to keep only the first and second seasons of some long lightcurves. 

::

	lcs = [lca, lcb] # That's a list of long lightcurves...
	pycs3.gen.lc_func.display(lcs)
	
	pycs3.gen.sea.easycut(lcs, keep=[1, 2])
	
	# Each lightcurve is processed individually. Check your results :

	pycs3.gen.lc_func.display(lcs)
	

If you are not happy with how the seasons where identified, try to add the option ``seasongap = 100`` to your call of :py:func:`pycs3.gen.sea.easycut`.
This is the number of days without points that start a new season. Default is 60.



Copying lightcurves
-------------------

... can be useful for instance to try out or compare things, and is very easy :

::

	testl = l.copy() # Makes a full deep copy of the entire lightcurve objects, with all properties, labels, mask, etc.
	
	testl.plotcolour = "blue"
	testl.shiftflux(5000)
	
	pycs3.gen.lc_func.display([l, testl])
	




Masking points
--------------

Each lightcurve object has a mask. This is simply a boolean numpy array of the same length as the curve. That's convenient, as such boolean arrays can be used to index normal numpy arrays. In the mask array, ``True`` means that the point is ok, ``False`` means that the point is masked.
Some demo of the flexibility :

::

	l.mask[17] = False # Manual way of masking a point
	l.mask[17:22] = False # Yes, it's a numpy array after all
	
	print(l.jds[l.mask]) # This gives you only the non-masked raw jds
 
 	l.mask = l.magerrs < 0.1 # Sets the mask to be False for all points with large errorbars.
	# Note that this would also set the mask of all other points to True.
	
.. note:: Masked points are shown with black circles on plots.


Some methods of lightcurve objects related to masks :

* :py:meth:`pycs3.gen.lc.lightcurve.hasmask`
* :py:meth:`pycs3.gen.lc.lightcurve.clearmask`
* :py:meth:`pycs3.gen.lc.lightcurve.cutmask`
* :py:meth:`pycs3.gen.lc.lightcurve.maskskiplist`
* :py:meth:`pycs3.gen.lc.lightcurve.maskinfo`


Buiding a mask "by hand"
""""""""""""""""""""""""

The best way to do this is to write a "skiplist" of the dates that you want to mask (this is much better than just specifying array indexes, as your skiplist will stay valid even if you merge/cut/tweak your curves). To help you writing such a list, use the function :py:meth:`pycs3.gen.lc.LightCurve.setjdlabels`. What are labels ? Labels are a bit like properties (see below), you can use them to attach any string to data points, and show them on plots. This particular functions puts the observation epochs as label to each point.
::

	for l in lcs:
		l.setjdlabels() # Sets the approximate epoch as label of each point.
		l.showlabels = True # Show the labels on plots
	pycs3.gen.lc_func.display(lcs)


Now you can write your skiplist; it's just a plain textfile with one line per data point to mask.
Any text following the MHJD is considered as a comment. One decimal is sufficient.
::

	# Some comment
	55111.3		Bad night

To apply this list to mask points of a curve, use the method :py:meth:`pycs3.gen.lc.LightCurve.maskskiplist` (click for details).
Of course you can use one file to set the same mask on A and B, or define separate masks.
::
	
	#l.clearmask() # Maybe you want to clear the mask first ?
	l.maskskiplist("myskiplist.txt")

This will mask the point within 0.2 days of the dates specified in the skiplist. You will be warned if there's anything fishy (like two separate points within 0.2 days or so).

Once you are happy with your masking, you could :
::
	
	for l in lcs:
		l.cutmask() # Removes all the masked points from your curve.
		l.clearlabels()
		l.showlabels = False

.. note:: It's a good idea to use cutmask to get "definitively" rid of points that you don't want to use *before* feeding the curves into a curve shifting algorithm. Some curve shifting methods might not accept curves with masked points.


Merging lightcurves
-------------------

When you import lightcurves from several telescopes, you might want to *merge* them, i.e. transform them into one single lightcurve object per quasar image. For instance to pass the resulting merged curves to some curve shifting algorithms.

.. note:: The operation described here is about merging any two lightcurve objects *as they are*. It does not involve optimizing any shifts between the curves so that they *match*. This is described later, in section :ref:`matchtels`. For now let's assume that you have for instance shifted your curves by hand (in magnitude and flux, not in time, usually...) so that they match.

There is a low-level method to merge one lightcurve into another one : :py:meth:`pycs3.gen.lc.LightCurve.merge` :
::

	# l and otherl are 2 lightcurve objects.
	
	otherl.shiftmag(0.23)
	
	pycs3.gen.lc_func.display([l, otherl])
	
	l.merge(otherl)
	
	print(l.longinfo())
	
	pycs3.gen.lc_func.display([l])
	

.. note:: Any lightcurve, at any time, has to be sorted according to its mhjds. We require that the jds
	are either increasing or (flat). This method thus takes care of this sorting ! Furthermore the properties, masks, labels etc are merged as well, as expected.

.. warning:: Any shifts of ``l`` or ``otherl`` will be *applied* to the data arrays, i.e. the resulting curve is no longer aware of previous shifts.

Often we want to merge a *list* of lightcurve from telescope 1 with a corresponding list of curves from telescope 2 and so on. :py:func:`pycs3.gen.mrg.merge` is a wrapper to do exactly this.
::
	
	# You have imported two lists of lightcurves : eulerlcs and mercatorlcs
	# Both lists contain n corresponding lightcurve objects, in the same order (image A, B, C and D).
	
	lcs = pycs3.gen.mrg.merge([eulerlcs, mercatorlcs])




Writing and reading  pickles
----------------------------

You will do this all the time, mostly with lightcurve objects. It allows to split up your workflow into different parts, making it a lot more effective and user-friendly. For instance, a first script imports your curves from various sources, masks some outliers and merges telescopes (i.e. all the stuff seen so far in this tutorial), and other scripts use these processed curves to measure the time delays. Writing and reading pickles is **the** easy-to-use connection between these scripts.

::
	
	# Say you have some lightcurves (perhaps just imported, or already heavily processed) :
	lcs = [lca, lcb, lcc, lcd]
	
	pycs3.gen.util.writepickle(lcs, "data/lcs_v2_merged.pkl") # Choose your own file name ...
	
	# And "later", in the next script :
	lcs = pycs3.gen.util.readpickle("data/lcs_v2_merged.pkl")


If you don't like these "``lcs``" lists, you are free to use other "containers" of your choice, like for instance dicts. Or just directly store one single lightcurve object into your pkl file.
You can of course also store other stuff using these same functions. If working with splines, this is typical (as we will see later in the tutorial) :

::
	
	pycs3.gen.util.writepickle((lcs, spline), "opt_test4.pkl")

	# And later ...
	
	(lcs, spline) = pycs3.gen.util.readpickle("opt_test4.pkl")


.. note:: Avoid relying on such pickle files to store actual data for eternity. Indeed the definitions of e.g. the LightCurve class might change, and this would make your pickles incompatible.


Writing lightcurves into rdb/ascii files
----------------------------------------

We come to the last point of this first chapter : what to do if your colleague doesn't accept pickle files ?
It is easy to write lightcurve objects into plain rdb files, using :py:meth:`pycs3.gen.lc.LightCurve.rdbexport` (click for details). This method nicely works together with :py:func:`pycs3.gen.lc_func.rdbimport`, in the sense that "written" lightcurves can then be "read" again :
::
	
	l.rdbexport(filename="test.txt", properties=["fwhm", "ellipticity"]) # l is a lightcurve object.
	
	imported_l = pycs3.gen.lc_func.rdbimport(filepath="test.txt", telescopename="Test", object="A", plotcolour="blue")

	pycs3.gen.lc_func.display([l, imported_l])


Both of these functions can handle properties. If you want to store properties in your exported file, you will have to specify them as optional arguments, as shown.
To see what properties are available, remember that you can use
::
	
	print (l.longinfo())
	

.. note:: As suggested by these functions, you should always write one file per lightcurve, when working with PyCS3. This is indeed natural, as you might have deleted or masked different points of a lightcurve. PyCS3 can perfectly process lightcurves of different lengths ! But on some occasions, you may want to to write several lightcurves into one single flat ascii file. For instance to submit to CDS... See function :py:func:`pycs3.gen.util.multilcsexport`.

As you probably expect, when writing a lightcurve object into an ASCII file, all "shifts" (and also microlensing models which we will see later) get applied to the datapoints before these are written to disk. Of course, when you then read the lightcurve again from this ASCII file, PyCS3 will no longer be aware that your lightcurve has previously been shifted.


pycs3\.gen package
==================

Submodules
----------

pycs3.gen.lc module
-------------------

.. automodule:: pycs3.gen.lc
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.lc_func module
------------------------

.. automodule:: pycs3.gen.lc_func
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.mrg module
--------------------

.. automodule:: pycs3.gen.mrg
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.polyml module
-----------------------

.. automodule:: pycs3.gen.polyml
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.sea module
--------------------

.. automodule:: pycs3.gen.sea
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.spl module
--------------------

.. automodule:: pycs3.gen.spl
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.spl_func module
-------------------------

.. automodule:: pycs3.gen.spl_func
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.splml module
----------------------

.. automodule:: pycs3.gen.splml
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.stat module
---------------------

.. automodule:: pycs3.gen.stat
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.gen.util module
---------------------

.. automodule:: pycs3.gen.util
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: pycs3.gen
    :members:
    :undoc-members:
    :show-inheritance:pycs3\.regdiff package
======================

Submodules
----------

pycs3.regdiff.multiopt module
-----------------------------

.. automodule:: pycs3.regdiff.multiopt
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.regdiff.rslc module
-------------------------

.. automodule:: pycs3.regdiff.rslc
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.regdiff.scikitgp module
-----------------------------

.. automodule:: pycs3.regdiff.scikitgp
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: pycs3.regdiff
    :members:
    :undoc-members:
    :show-inheritance:pycs3\.pipe package
===================

Submodules
----------

pycs3.pipe.optimiser module
---------------------------
.. automodule:: pycs3.pipe.optimiser
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.pipe.pipe_utils module
----------------------------

.. automodule:: pycs3.pipe.pipe_utils
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: pycs3.pipe
    :members:
    :undoc-members:
    :show-inheritance:pycs3 package
=============

Subpackages
-----------

.. toctree::

    pycs3.gen
    pycs3.regdiff
    pycs3.sim
    pycs3.spl
    pycs3.pipe
    pycs3.tdcomb


Module contents
---------------

.. automodule:: pycs3
    :members:
    :undoc-members:
    :show-inheritance:pycs3\.sim package
==================

Submodules
----------

pycs3.sim.draw module
---------------------

.. automodule:: pycs3.sim.draw
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.sim.plot module
---------------------

.. automodule:: pycs3.sim.plot
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.sim.power_spec module
---------------------------

.. automodule:: pycs3.sim.power_spec
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.sim.run module
--------------------

.. automodule:: pycs3.sim.run
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.sim.src module
--------------------

.. automodule:: pycs3.sim.src
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.sim.twk module
--------------------

.. automodule:: pycs3.sim.twk
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: pycs3.sim
    :members:
    :undoc-members:
    :show-inheritance:pycs3\.spl package
==================

Submodules
----------

pycs3.spl.multiopt module
-------------------------

.. automodule:: pycs3.spl.multiopt
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.spl.topopt module
-----------------------

.. automodule:: pycs3.spl.topopt
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------

.. automodule:: pycs3.spl
    :members:
    :undoc-members:
    :show-inheritance:pycs3\.tdcomb package
=====================

Submodules
----------

pycs3.tdcomb.comb module
------------------------

.. automodule:: pycs3.tdcomb.comb
    :members:
    :undoc-members:
    :show-inheritance:

pycs3.tdcomb.plot module
------------------------

.. automodule:: pycs3.tdcomb.plot
    :members:
    :undoc-members:
    :show-inheritance:

Module contents
---------------
.. automodule:: pycs3.tdcomb
    :members:
    :undoc-members:
    :show-inheritance: