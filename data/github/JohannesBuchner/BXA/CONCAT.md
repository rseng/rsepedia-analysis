.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/JohannesBuchner/BXA/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Sherpa or xspec version
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

BXA could always use more documentation, whether as part of the
official BXA docs, in docstrings, or even on the web in blog posts,
articles, tutorials, and such.

Notebooks demonstrating how to use BXA are also appreciated.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/JohannesBuchner/BXA/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `BXA` for local development.

1. Fork the `BXA` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:JohannesBuchner/BXA.git

3. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

4. Check that your changes still allow the examples to run::

    $ cd examples/xspec/
    $ PYTHONPATH=../../:$PYTHONPATH bash runall.sh

5. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

6. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests or be covered by an example.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for up-to-date Python versions. Check
   https://travis-ci.org/JohannesBuchner/BXA/pull_requests
   and make sure that the tests pass for all supported Python versions.

Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ bump2version patch # possible: major / minor / patch
$ git push
$ git push --tags

Travis will then deploy to PyPI if tests pass.
==============
Release Notes
==============

4.0.0 (2021-01-06)
------------------

* Make ultranest default, remove multinest.

3.4 (2019-05-26)
------------------

3.3 (2020-01-28)
------------------

3.2 (2020-01-12)
------------------

* circumvent xspec numerical interpolation issues (Lepsilon parameter)

3.1 (2019-10-21)
------------------

* Make multinest optional, allow ultranest

3.0 (2019-10-15)
------------------

2.10 (2019-06-12)
------------------

* Added PCA background models

2.9 (2019-06-11)
-----------------

2.8 (2019-05-31)
-----------------

2.7 (2019-05-31)
-----------------

2.6 (2019-05-21)
-----------------

* replace outdated progressbar with tqdm

2.5 (2019-04-08)
-----------------

* added PCA background model
* added galactic absorption fetcher
* added generic AGN fitting script

2.4 (2015-06-17)
-----------------

2.3 (2015-06-09)
-----------------

* added  Chandra background model
* added  Swift background model

2.0.0 (2015-06-05)
------------------

* Python 3 support
* added corner plots
* added more convenience functions for plotting
* added Chandra background model
* added acceleration of models

1.0.0 (2013-12-01)
------------------

* simple multinest interfaces for xspec and sherpa based on pyblocxs code
About Bayesian X-ray Analysis (BXA)
------------------------------------

BXA connects the X-ray spectral analysis environments Xspec/Sherpa
to the nested sampling algorithm UltraNest 
for **Bayesian Parameter Estimation** and **Model comparison**.

BXA provides the following features:

* parameter estimation in arbitrary dimensions, which involves:
   * finding the best fit
   * computing error bars
   * computing marginal probability distributions
   * parallelisation with MPI
* plotting of spectral model vs. the data:
   * for the best fit
   * for each of the solutions (posterior samples)
   * for each component
* model selection:
   * computing the evidence for the considered model, 
     ready for use in Bayes factors
   * unlike likelihood-ratios, not limited to nested models 
* model discovery:
   * visualize deviations between model and data with Quantile-Quantile (QQ) plots.
     QQ-plots do not require binning and are more comprehensive than residuals.
     This will give you ideas on when to introduce more complex models, which 
     may again be tested with model selection

BXA shines especially

* when systematically analysing a large data-set, or
* when comparing multiple models
* when analysing low counts data-set with realistic models

because its robust and unsupervised fitting algorithm explores
even complicated parameter spaces in an automated fashion.
The user does not need to initialise to good starting points.
The `algorithm <https://johannesbuchner.github.io/UltraNest/method.html>`_ automatically runs until convergence, and slows down to sample
carefully if complicated parameter spaces are encountered. This allows building automated analysis pipelines.

.. image:: https://img.shields.io/pypi/v/BXA.svg
        :target: https://pypi.python.org/pypi/BXA

.. image:: https://coveralls.io/repos/github/JohannesBuchner/BXA/badge.svg
        :target: https://coveralls.io/github/JohannesBuchner/BXA

.. image:: https://img.shields.io/badge/docs-published-ok.svg
        :target: https://johannesbuchner.github.io/BXA/
        :alt: Documentation Status

.. image:: https://img.shields.io/badge/GitHub-JohannesBuchner%2FBXA-blue.svg?style=flat
        :target: https://github.com/JohannesBuchner/BXA/
        :alt: Github repository

Who is using BXA?
-------------------------------

* Dr. Antonis Georgakakis, Dr. Angel Ruiz (NOA, Athens)
* Dr. Mike Anderson (MPA, Munich)
* Dr. Franz Bauer, Charlotte Simmonds (PUC, Jonathan Quirola Vásquez, Santiago)
* Dr. Stéphane Paltani, Dr. Carlo Ferrigno (ISDC, Geneva)
* Dr. Zhu Liu (NAO, Beijing)
* Dr. Georgios Vasilopoulos (Yale, New Haven)
* Dr. Francesca Civano, Dr. Aneta Siemiginowska (CfA/SAO, Cambridge)
* Dr. Teng Liu, Adam Malyali, Riccardo Arcodia, Sophia Waddell, Torben Simm, ... (MPE, Garching)
* Dr. Sibasish Laha, Dr. Alex Markowitz (UCSD, San Diego)
* Dr. Arash Bahramian (Curtin University, Perth)
* Dr. Peter Boorman (U of Southampton, Southampton; ASU, Prague)
* and `you <https://ui.adsabs.harvard.edu/search/q=citations(bibcode%3A2014A%26A...564A.125B)%20full%3A%22BXA%22&sort=date%20desc%2C%20bibcode%20desc&p_=0>`_?

Documentation
----------------

BXA's `documentation <http://johannesbuchner.github.io/BXA/>`_ is hosted at http://johannesbuchner.github.io/BXA/

Installation
-------------

First, you need to have either `Sherpa`_ or `Xspec`_ installed and its environment loaded.

BXA itself can installed easily using pip or conda::

	$ pip install bxa

If you want to install in your home directory, install with::

	$ pip install bxa --user

The following commands should not yield any error message::

	$ python -c 'import ultranest'
	$ python -c 'import xspec'
	$ sherpa

You may need to install python and some basic packages through your package manager. For example::

	$ yum install ipython python-matplotlib scipy numpy matplotlib
	$ apt-get install python-numpy python-scipy python-matplotlib ipython

BXA requires the following python packages: requests corner astropy h5py cython scipy tqdm.
They should be downloaded automatically. If they are not, install them
also with pip/conda.

The source code is available from https://github.com/JohannesBuchner/BXA,
so alternatively you can download and install it::
	
	$ git clone https://github.com/JohannesBuchner/BXA
	$ cd BXA
	$ python setup.py install

Or if you only want to install it for the current user::

	$ python setup.py install --user

**Supported operating systems**: 
BXA runs on all operating systems supported by 
`ciao/sherpa <https://cxc.cfa.harvard.edu/ciao/watchout.html#install>`_ or 
`heasoft/xspec <https://heasarc.gsfc.nasa.gov/lheasoft/issues.html>`_.
The support is systematically tested for every BXA release by 
`Travis CI <https://travis-ci.com/github/JohannesBuchner/BXA>`_, but only for Ubuntu Linux.


Running
--------------

In *Sherpa*, load the package::

	jbuchner@ds42 ~ $ sherpa
	-----------------------------------------------------
	Welcome to Sherpa: CXC's Modeling and Fitting Package
	-----------------------------------------------------
	CIAO 4.4 Sherpa version 2 Tuesday, June 5, 2012

	sherpa-1> import bxa.sherpa as bxa
	sherpa-2> bxa.BXASolver?

For *Xspec*, start python or ipython::
	
	jbuchner@ds42 ~ $ ipython
	In [1]: import xspec
	
	In [2]: import bxa.xspec as bxa
	
	In [3]:	bxa.BXASolver?

Now you can use BXA. See the documentation pages for how
to perform analyses. Several examples are included.

.. _ultranest: http://johannesbuchner.github.io/UltraNest/

.. _Sherpa: http://cxc.harvard.edu/sherpa/

.. _Xspec: http://heasarc.gsfc.nasa.gov/docs/xanadu/xspec/

Code
-------------------------------

See the `code repository page <https://github.com/JohannesBuchner/BXA>`_ 

.. _cite:

Citing BXA correctly
---------------------

Refer to the `accompaning paper Buchner et al. (2014) <http://www.aanda.org/articles/aa/abs/2014/04/aa22971-13/aa22971-13.html>`_ which gives introduction and 
detailed discussion on the methodology and its statistical footing.

We suggest giving credit to the developers of Sherpa/Xspec, UltraNest and of this software.
As an example::

	For analysing X-ray spectra, we use the analysis software BXA (\ref{Buchner2014}),
	which connects the nested sampling algorithm UltraNest (\ref{ultranest})
	with the fitting environment CIAO/Sherpa (\ref{Fruscione2006}).

Where the BibTex entries are:

* for BXA and the contributions to X-ray spectral analysis methodology (model comparison, model discovery, Experiment design, Model discovery through QQ-plots):

	- Buchner et al. (2014) A&A
	- The paper is available at `arXiv:1402.0004 <http://arxiv.org/abs/arXiv:1402.0004>`_
	- `bibtex entry <https://ui.adsabs.harvard.edu/abs/2014A%26A...564A.125B/exportcitation>`_

* for UltraNest: see https://johannesbuchner.github.io/UltraNest/issues.html#how-should-i-cite-ultranest
* for Sherpa: see `Sherpa`_
* for Xspec: see `Xspec`_
================================
Fitting with BXA & XAGNFitter
================================

With this tool you can automatically fit X-ray spectra of AGN.

What it does
---------------

* Fit the source spectrum with 

  * Obscured AGN component (powerlaw reprocessed by a torus structure)
  * Warm mirror component (powerlaw, Thompson scattering within the opening angle)
  * Stellar process component (apec whose luminosity remains below 1e42 erg/s)
  * All three are absorbed by Milky Way absorption.

* Additional features:

  * Uses the MultiNest global parameter space exploration algorithm, which can deal with many parameters and does not get stuck in local minima.
  * If the redshift is uncertain, these uncertainties are propagated through. Or you can give a fixed value.
  * The background spectrum is approximated automatically with an empirical model. This extracts more information compared to default xspec wstat, even for very few source counts (see e.g., Simmonds+18).

* Produces parameter posterior distributions, which describe the uncertainties in:

  * AGN Luminosity           (log-uniform uninformative prior)
  * Column density N_H       (log-uniform uninformative prior)
  * Powerlaw photon index    (1.95+-0.15 informative prior)
  * Redshift (if not fixed)  (user-supplied)
  * apec temperature         (log-uniform uninformative prior)
  * apec normalisation       (log-uniform uninformative prior)
  * torus opening angle parameters (CTKcover, TORtheta, uniform uninformative priors) 


Preparing your data
---------------------

* You only need your spectrum as a .pi/.pha file and the associated RMF, ARF and background files.
* Check that the keywords are set correctly:

  * You can use the fixkeywords.py script (see BXA repository) to make sure the src points to the ARF, RMF and background files
  * python fixkeywords.py pn_src.fits pn_bgd.fits pn_rmf.fits pn_arf.fits

* Set the galactic N_H

  * Create a file next to your spectrum file with the extention .nh and put in the galactic column density to the source. For example:
  * $ echo 1e21 > combined_src.pi.nh
  * You can also use the gal.py script (see BXA repository). It automatically fetches the nhtot value from http://www.swift.ac.uk/analysis/nhtot/ if RA/DEC are set:
  * python gal.py pn_src.fits

* Give redshift information

  * Create a file next to your spectrum file with the extention .z and put in the source redshift. For example:
  * $ echo 6.31 > combined_src.pi.z
  * If you do not know the redshift, don't create this file.
  * If you have redshift uncertainties as a probability distribution, store here two columns: cumulative probability (from 0 to 1, uniformly sampled) and corresponding redshift. For example, for a redshift 1.1+-0.2, you could do::

	cdfsteps = numpy.linspace(0, 1, 100)
	z = scipy.stats.norm(1.1, 0.2).cdf(cdfsteps)
	numpy.savetxt("combined_src.pi.z", numpy.transpose([cdfsteps, z]))


Finally, you should have files like in the testsrc/ folder:

* combined_src.pi
* combined_src.pi.nh (galactic NH)
* combined_src.pi.z  (redshift information)
* combined_src.arf
* combined_src.rmf
* combined_bkg.arf  (optional)
* combined_bkg.pi 
* combined_bkg.rmf  (optional)


How to use
---------------

You only need to **install docker** on your computer. 

All the necessary software (ciao, sherpa, multinest, BXA) is in a container image that will be automatically downloaded for you.

First, go to the directory where your spectrum is::

	$ cd testsrc

Optionally, allow docker containers to show stuff on your screen::

	$ xhost +local:`docker inspect --format='{{ .Config.Hostname }}' johannesbuchner/bxa_absorbed` 
	$ XDOCKERARGS="-v /tmp/.X11-unix:/tmp/.X11-unix -e DISPLAY=$DISPLAY"

Then, run docker::

	$ docker run $XDOCKERARGS -v $PWD:/opt/example/ -e FILENAME=combined_src.pi -e ELO=0.5 -e EHI=8 -ti johannesbuchner/bxa_absorbed 

The arguments mean the following:

* -v arguments mount local directories into the virtual machine:

  * -v $PWD:/opt/example/ -- This connects the /opt/example folder inside the container, to your host machines folder $PWD (the folder you are currently in). Otherwise the container cannot read your files.

* -e arguments set environment variables

  * -e FILENAME=combined_src.pi sets the file to analyse
  * -e ELO=0.5 -e EHI=8 sets the energy range (in keV) to use

* -ti -- Get a interactive terminal
* johannesbuchner/bxa_absorbed -- name of the container image to run


What it will do internally
----------------------------

1. Initialises the container. It is similar to a virtual machine, with a completely separated Ubuntu Linux installation inside.
2. Inside, it run the command ". /opt/ciao-4.8/bin/ciao.sh; sherpa /opt/scripts/fitagn.py". This loads ciao and runs sherpa. 

  * If you want to replace or edit the fitting script, add  in the docker run command "-v mydirectory/scripts:/opt/scripts/" before "-ti" to replace the /opt/scripts folder with your own scripts folder.

3. The sherpa script sets up the source model, parameters and priors. The background is fitted. Finally, multinest is run to constrain the parameters. 

4. Output files are combined_src.pi_out_*. The most important ones are

  * params.json -- contains the parameter names 
  * post_equal_weights.dat -- contains the posterior samples, each column is a parameter

What you can do with the results
-------------------------------------

* Plot the parameter distribution. The multinest_marginals.py tool in the pymultinest repository can help::

  $ python pymultinest/multinest_marginals.py combined_src.piout_withapec_

Alternatively, you can also use corner.py or any other plotting tool.

Modify the behaviour
-------------------------

* to change redshift, alter the .z file (see above)
* to disable apec component, set the environment variable "-e WITHAPEC=0"
* to change the fitagn.py script altogether, edit it in the scripts/ folder and pass "-v mydirectory/scripts:/opt/scripts/". You have to give the absolute path to your scripts/ directory.






BXA/Sherpa example scripts
==========================

This folder contains simple and complex examples
how BXA can be invoked for spectral analysis.

Please refer to https://johannesbuchner.github.io/BXA/
for full documentation, including how to install BXA.


Test data
-------------------

As test data, this includes the spectral file example-file.fak
representing a ATHENA observation of an absorbed AGN, and the corresponding 
response.

The spectrum looks something like this:

.. image:: ../xspec/reference-output/data.gif

Simple analysis
-----------------

Have a look at the file example_simplest.py. It contains:

* Loading data
* setting up a model and its parameter ranges
* running a BXA fit with specified priors

See https://johannesbuchner.github.io/BXA/ to understand the code.
See https://johannesbuchner.github.io/UltraNest/ to understand the output of the
fitting engine (for example, its `FAQ page <https://johannesbuchner.github.io/UltraNest/issues.html>`_).

Expected output::

	$ python3 example_simplest.py
	read RMF file athenapp_ir_b4c_wfi_withfilter_fov40.0arcmin_avg.rsp
	[ultranest] Sampling 400 live points from prior ...


	Mono-modal Volume: ~exp(-4.69) * Expected Volume: exp(0.00) Quality: ok

	mypow.gamma:      +1.0|************************************************************************************|     +3.0
	mypow.ampl :  +1.0e-10|************ ** * * *  ***    * * * * *  *      * * *    * *        *     **  *   * | +1.0e+01

	Z=-1436853.7(0.00%) | Like=-1425796.15..-4079.21 [-1.09e+08..-4333] | it/evals=80/481 eff=98.7654% N=400 0 

	Mono-modal Volume: ~exp(-4.69)   Expected Volume: exp(-0.23) Quality: ok

	mypow.gamma:      +1.0|************************************************************************************|     +3.0
	mypow.ampl :  +1.0e-10|*********** *** ***  * *** ***   **        **  *      *      * ** *  +7.9e-02       | +1.0e-01

	Z=-49990.0(0.00%) | Like=-48482.70..-4079.21 [-1.09e+08..-4333] | it/evals=160/566 eff=96.3855% N=400 0 

	Mono-modal Volume: ~exp(-4.69)   Expected Volume: exp(-0.45) Quality: ok

	mypow.gamma:      +1.0|************************************************************************************|     +3.0
	mypow.ampl :  +1.0e-10|************* ***   *   * *  +3.2e-03                                               | +1.0e-02

	Z=-28599.1(0.00%) | Like=-28584.11..-4079.21 [-1.09e+08..-4333] | it/evals=240/650 eff=96.0000% N=400 

	Mono-modal Volume: ~exp(-5.23) * Expected Volume: exp(-0.67) Quality: ok

	mypow.gamma:      +1.0|**************************** *******************************************************|     +3.0
	mypow.ampl :  +1.0e-10|************* ****  *  +2.4e-03                                                     | +1.0e-02

	Z=-23814.6(0.00%) | Like=-23785.78..-4079.21 [-1.09e+08..-4333] | it/evals=320/748 eff=91.9540% N=400 

	Mono-modal Volume: ~exp(-5.23)   Expected Volume: exp(-0.90) Quality: ok

	mypow.gamma:      +1.0|******************************************************************** ***************|     +3.0
	mypow.ampl :  +1.0e-10|************* * ** *  +2.3e-03                                                      | +1.0e-02

	Z=-18339.2(0.00%) | Like=-18303.01..-4079.21 [-1.09e+08..-4333] | it/evals=440/880 eff=91.6667% N=400 

	Mono-modal Volume: ~exp(-5.23)   Expected Volume: exp(-1.12) Quality: ok

	mypow.gamma:      +1.0|********************************************************************** *************|     +3.0
	mypow.ampl :  +1.0e-10|*************   *  +2.0e-03                                                         | +1.0e-02

	Z=-15475.0(0.00%) | Like=-15444.04..-4067.92 [-1.09e+08..-4333] | it/evals=520/977 eff=90.1213% N=400 

	Mono-modal Volume: ~exp(-5.52) * Expected Volume: exp(-1.35) Quality: ok

	mypow.gamma:      +1.0|************************************ ******************************************* ***|     +3.0
	mypow.ampl :  +1.0e-10|************  +1.4e-03                                                              | +1.0e-02

	Z=-13226.8(0.00%) | Like=-13216.88..-4067.92 [-1.09e+08..-4333] | it/evals=600/1070 eff=89.5522% N=400 

	Mono-modal Volume: ~exp(-5.52)   Expected Volume: exp(-1.57) Quality: ok

	mypow.gamma:      +1.0|******************************************* ********************************** **** |     +3.0
	mypow.ampl :  +1.0e-10|***********  +1.3e-03                                                               | +1.0e-02

	Z=-11727.7(0.00%) | Like=-11707.66..-4042.40 [-1.09e+08..-4333] | it/evals=680/1175 eff=87.7419% N=400 

	Mono-modal Volume: ~exp(-5.82) * Expected Volume: exp(-1.80) Quality: ok

	mypow.gamma:      +1.0|******************************************* ********************************** **** |     +3.0
	mypow.ampl :  +1.0e-10|**********  +1.1e-03                                                                | +1.0e-02

	Z=-9850.1(0.00%) | Like=-9833.59..-4042.40 [-1.09e+08..-4333] | it/evals=800/1321 eff=86.8621% N=400 0 

	Mono-modal Volume: ~exp(-6.36) * Expected Volume: exp(-2.02) Quality: ok

	mypow.gamma:      +1.0|**********************************************************************  +2.7        |     +3.0
	mypow.ampl :  +1.0e-10|*********  +1.0e-03                                                                 | +1.0e-02

	Z=-9111.3(0.00%) | Like=-9088.80..-4014.82 [-1.09e+08..-4333] | it/evals=880/1422 eff=86.1057% N=400 

	Mono-modal Volume: ~exp(-6.36)   Expected Volume: exp(-2.25) Quality: ok

	mypow.gamma:      +1.0|************************************************************** *  +2.5              |     +3.0
	mypow.ampl :  +1.0e-10|*****************************************  ** ****    *** ****  *** **   *  +8.8e-04| +1.0e-03

	Z=-8390.6(0.00%) | Like=-8379.03..-4014.82 [-1.09e+08..-4333] | it/evals=960/1541 eff=84.1367% N=400 

	Mono-modal Volume: ~exp(-6.66) * Expected Volume: exp(-2.47) Quality: ok

	mypow.gamma:      +1.0|*******************************************************  +2.3                       |     +3.0
	mypow.ampl :  +1.0e-10|********************************************* ****  ***** ****  *** **   *  +8.8e-04| +1.0e-03

	Z=-7775.5(0.00%) | Like=-7762.79..-4014.82 [-1.09e+08..-4333] | it/evals=1040/1650 eff=83.2000% N=400 

	Mono-modal Volume: ~exp(-6.97) * Expected Volume: exp(-2.70) Quality: ok

	mypow.gamma:      +1.0|*************************************************  +2.1                             |     +3.0
	mypow.ampl :  +1.0e-10| ************************************* ****** ****  ***** ****  +7.3e-04            | +1.0e-03

	Z=-7136.5(0.00%) | Like=-7127.33..-4014.82 [-1.09e+08..-4333] | it/evals=1160/1803 eff=82.6800% N=400 

	Mono-modal Volume: ~exp(-7.29) * Expected Volume: exp(-2.92) Quality: ok

	mypow.gamma:      +1.0|*******************************************  +2.0                                   |     +3.0
	mypow.ampl :  +1.0e-10| ************************************* ****** ****** **** * *  +7.2e-04             | +1.0e-03

	Z=-6781.9(0.00%) | Like=-6772.22..-4014.82 [-1.09e+08..-4333] | it/evals=1240/1895 eff=82.9431% N=400 

	Mono-modal Volume: ~exp(-7.59) * Expected Volume: exp(-3.15) Quality: ok

	mypow.gamma:      +1.0|**************************************  +1.9                                        |     +3.0
	mypow.ampl :  +1.0e-10| ************************************* ************* ****  +6.8e-04                 | +1.0e-03

	Z=-6439.0(0.00%) | Like=-6429.44..-4014.82 [-1.09e+08..-4333] | it/evals=1320/2008 eff=82.0896% N=400 

	Mono-modal Volume: ~exp(-7.59)   Expected Volume: exp(-3.37) Quality: ok

	mypow.gamma:      +1.0|*********************************  +1.8                                             |     +3.0
	mypow.ampl :  +1.0e-10|  **************************************** ********* *  +6.4e-04                    | +1.0e-03

	Z=-6156.2(0.00%) | Like=-6135.97..-4014.82 [-1.09e+08..-4333] | it/evals=1400/2118 eff=81.4901% N=400 

	Mono-modal Volume: ~exp(-8.05) * Expected Volume: exp(-3.60) Quality: ok

	mypow.gamma:      +1.0|************************** *  +1.7                                                  |     +3.0
	mypow.ampl :  +1.0e-10|   *************************************** *** **  *  +6.1e-04                      | +1.0e-03

	Z=-5693.3(0.00%) | Like=-5673.11..-3982.43 [-1.09e+08..-4333] | it/evals=1520/2263 eff=81.5888% N=400 

	Mono-modal Volume: ~exp(-8.16) * Expected Volume: exp(-3.82) Quality: ok

	mypow.gamma:      +1.0|*************************  +1.6                                                     |     +3.0
	mypow.ampl :  +1.0e-10|   ***************************************  **  +5.4e-04                            | +1.0e-03

	Z=-5437.6(0.00%) | Like=-5420.08..-3975.18 [-1.09e+08..-4333] | it/evals=1600/2362 eff=81.5494% N=400 

	Mono-modal Volume: ~exp(-8.16)   Expected Volume: exp(-4.05) Quality: ok

	mypow.gamma:      +1.0|**********************  +1.5                                                        |     +3.0
	mypow.ampl :  +1.0e-10|   *************************************** *  +5.2e-04                              | +1.0e-03

	Z=-5275.0(0.00%) | Like=-5263.43..-3975.18 [-1.09e+08..-4333] | it/evals=1680/2469 eff=81.1986% N=400 

	Mono-modal Volume: ~exp(-8.40) * Expected Volume: exp(-4.27) Quality: ok

	mypow.gamma:      +1.0|******************  +1.4                                                            |     +3.0
	mypow.ampl :  +1.0e-10|    ******************************** *****  +5.0e-04                                | +1.0e-03

	Z=-5089.4(0.00%) | Like=-5078.33..-3975.18 [-1.09e+08..-4333] | it/evals=1760/2567 eff=81.2183% N=400 

	Mono-modal Volume: ~exp(-8.83) * Expected Volume: exp(-4.50) Quality: ok

	mypow.gamma:      +1.0|****************  +1.4                                                              |     +3.0
	mypow.ampl :  +1.0e-10|     *********************************  +4.5e-04                                    | +1.0e-03

	Z=-4863.3(0.00%) | Like=-4852.01..-3975.18 [-1.09e+08..-4333] | it/evals=1880/2729 eff=80.7213% N=400 

	Mono-modal Volume: ~exp(-8.83)   Expected Volume: exp(-4.73) Quality: ok

	mypow.gamma:      +1.0|**************  +1.3                                                                |     +3.0
	mypow.ampl :  +1.0e-10|     ********************************  +4.3e-04                                     | +1.0e-03

	Z=-4752.0(0.00%) | Like=-4740.26..-3975.18 [-1.09e+08..-4333] | it/evals=1960/2829 eff=80.6916% N=400 

	Mono-modal Volume: ~exp(-9.23) * Expected Volume: exp(-4.95) Quality: ok

	mypow.gamma:      +1.0|************  +1.3                                                                  |     +3.0
	mypow.ampl :  +1.0e-10|      ************************** *  +4.0e-04                                        | +1.0e-03

	Z=-4640.0(0.00%) | Like=-4628.87..-3975.18 [-1.09e+08..-4333] | it/evals=2040/2932 eff=80.5687% N=400 

	Mono-modal Volume: ~exp(-9.40) * Expected Volume: exp(-5.18) Quality: ok

	mypow.gamma:      +1.0|**********  +1.2                                                                    |     +3.0
	mypow.ampl :  +1.0e-10|       *************************  +3.8e-04                                          | +1.0e-03

	Z=-4552.1(0.00%) | Like=-4540.22..-3975.18 [-1.09e+08..-4333] | it/evals=2120/3036 eff=80.4249% N=400 

	Mono-modal Volume: ~exp(-9.40)   Expected Volume: exp(-5.40) Quality: ok

	mypow.gamma:      +1.0|*********  +1.2                                                                     |     +3.0
	mypow.ampl :  +1.0e-10|       ************************  +3.7e-04                                           | +1.0e-03

	Z=-4456.3(0.00%) | Like=-4444.47..-3952.75 [-1.09e+08..-4333] | it/evals=2240/3192 eff=80.2292% N=400 

	Mono-modal Volume: ~exp(-9.76) * Expected Volume: exp(-5.63) Quality: ok

	mypow.gamma:      +1.0|********  +1.2                                                                      |     +3.0
	mypow.ampl :  +1.0e-10|       ************************  +3.6e-04                                           | +1.0e-03

	Z=-4399.4(0.00%) | Like=-4388.16..-3952.75 [-1.09e+08..-4333] | it/evals=2320/3305 eff=79.8623% N=400 

	Mono-modal Volume: ~exp(-10.34) * Expected Volume: exp(-5.85) Quality: ok

	mypow.gamma:      +1.0|*******  +1.2                                                                       |     +3.0
	mypow.ampl :  +1.0e-10|        *********************  +3.4e-04                                             | +1.0e-03

	Z=-4349.4(0.00%) | Like=-4338.36..-3952.75 [-1.09e+08..-4333] | it/evals=2400/3402 eff=79.9467% N=400 

	Mono-modal Volume: ~exp(-10.34)   Expected Volume: exp(-6.08) Quality: ok

	mypow.gamma:      +1.0|******  +1.1                                                                        |     +3.0
	mypow.ampl :   +0.0000|        ********************  +0.0003                                               |  +0.0010

	Z=-4300.8(0.00%) | Like=-4289.26..-3952.75 [-4332.2367..-3960.1176] | it/evals=2480/3509 eff=79.7684% N=400 

	Mono-modal Volume: ~exp(-10.59) * Expected Volume: exp(-6.30) Quality: ok

	mypow.gamma:      +1.0|******  +1.1                                                                        |     +3.0
	mypow.ampl :   +0.0000|         ******************  +0.0003                                                |  +0.0010

	Z=-4243.0(0.00%) | Like=-4230.88..-3952.75 [-4332.2367..-3960.1176] | it/evals=2600/3665 eff=79.6325% N=400 

	Mono-modal Volume: ~exp(-10.77) * Expected Volume: exp(-6.53) Quality: ok

	mypow.gamma:      +1.0|*****  +1.1                                                                         |     +3.0
	mypow.ampl :   +0.0000|         ******************  +0.0003                                                |  +0.0010

	Z=-4208.7(0.00%) | Like=-4197.23..-3952.75 [-4332.2367..-3960.1176] | it/evals=2680/3756 eff=79.8570% N=400 

	Mono-modal Volume: ~exp(-10.77)   Expected Volume: exp(-6.75) Quality: ok

	mypow.gamma:      +1.0|****  +1.1                                                                          |     +3.0
	mypow.ampl :   +0.0000|          ****************  +0.0003                                                 |  +0.0010

	Z=-4180.7(0.00%) | Like=-4168.58..-3952.75 [-4332.2367..-3960.1176] | it/evals=2760/3872 eff=79.4931% N=400 

	Mono-modal Volume: ~exp(-11.38) * Expected Volume: exp(-6.98) Quality: ok

	mypow.gamma:      +1.0|****  +1.1                                                                          |     +3.0
	mypow.ampl :   +0.0000|          ***************  +0.0003                                                  |  +0.0010

	Z=-4148.6(0.00%) | Like=-4135.69..-3952.75 [-4332.2367..-3960.1176] | it/evals=2840/3967 eff=79.6187% N=400 

	Mono-modal Volume: ~exp(-11.38)   Expected Volume: exp(-7.20) Quality: ok

	mypow.gamma:      +1.0|***  +1.1                                                                           |     +3.0
	mypow.ampl :   +0.0000|          **************  +0.0003                                                   |  +0.0010

	Z=-4109.3(0.00%) | Like=-4097.52..-3952.75 [-4332.2367..-3960.1176] | it/evals=2960/4125 eff=79.4631% N=400 

	Mono-modal Volume: ~exp(-11.56) * Expected Volume: exp(-7.43) Quality: ok

	mypow.gamma:      +1.0|***  +1.1                                                                           |     +3.0
	mypow.ampl :   +0.0000|  +0.0001  *************  +0.0003                                                   |  +0.0010

	Z=-4091.7(0.00%) | Like=-4079.63..-3952.75 [-4332.2367..-3960.1176] | it/evals=3040/4236 eff=79.2492% N=400 

	Mono-modal Volume: ~exp(-11.56)   Expected Volume: exp(-7.65) Quality: ok

	mypow.gamma:     +1.00|***  +1.05                                                                          |    +3.00
	mypow.ampl :   +0.0000|  +0.0001  ************  +0.0003                                                    |  +0.0010

	Z=-4077.2(0.00%) | Like=-4064.74..-3952.75 [-4332.2367..-3960.1176] | it/evals=3120/4330 eff=79.3893% N=400 

	Mono-modal Volume: ~exp(-12.19) * Expected Volume: exp(-7.88) Quality: ok

	mypow.gamma:     +1.00|**  +1.04                                                                           |    +3.00
	mypow.ampl :   +0.0000|  +0.0001  ************  +0.0003                                                    |  +0.0010

	Z=-4059.8(0.00%) | Like=-4046.99..-3952.75 [-4332.2367..-3960.1176] | it/evals=3200/4432 eff=79.3651% N=400 

	Mono-modal Volume: ~exp(-12.23) * Expected Volume: exp(-8.10) Quality: ok

	mypow.gamma:     +1.00|**  +1.03                                                                           |    +3.00
	mypow.ampl :   +0.0000|   +0.0001  **********  +0.0003                                                     |  +0.0010

	Z=-4042.9(0.00%) | Like=-4030.36..-3952.75 [-4332.2367..-3960.1176] | it/evals=3320/4600 eff=79.0476% N=400 

	Mono-modal Volume: ~exp(-12.67) * Expected Volume: exp(-8.33) Quality: ok

	mypow.gamma:     +1.00|**  +1.03                                                                           |    +3.00
	mypow.ampl :   +0.0000|   +0.0001  **********  +0.0003                                                     |  +0.0010

	Z=-4032.5(0.00%) | Like=-4020.32..-3952.75 [-4332.2367..-3960.1176] | it/evals=3400/4706 eff=78.9596% N=400 

	Mono-modal Volume: ~exp(-12.67)   Expected Volume: exp(-8.55) Quality: ok

	mypow.gamma:     +1.00|**  +1.03                                                                           |    +3.00
	mypow.ampl :  +0.00000|  +0.00015  *********  +0.00025                                                     | +0.00100

	Z=-4025.4(0.00%) | Like=-4012.93..-3952.75 [-4332.2367..-3960.1176] | it/evals=3480/4810 eff=78.9116% N=400 

	Mono-modal Volume: ~exp(-12.67)   Expected Volume: exp(-8.78) Quality: ok

	mypow.gamma:     +1.00|*  +1.02                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00015  ********  +0.00025                                                     | +0.00100

	Z=-4018.9(0.00%) | Like=-4007.04..-3952.58 [-4332.2367..-3960.1176] | it/evals=3560/4926 eff=78.6567% N=400 

	Mono-modal Volume: ~exp(-13.22) * Expected Volume: exp(-9.00) Quality: ok

	mypow.gamma:     +1.00|*  +1.02                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00016  ********  +0.00024                                                     | +0.00100

	Z=-4008.1(0.00%) | Like=-3995.35..-3952.58 [-4332.2367..-3960.1176] | it/evals=3680/5093 eff=78.4147% N=400 

	Mono-modal Volume: ~exp(-13.22)   Expected Volume: exp(-9.23) Quality: ok

	mypow.gamma:     +1.00|*  +1.02                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00016  *******  +0.00024                                                      | +0.00100

	Z=-4002.5(0.00%) | Like=-3989.71..-3952.58 [-4332.2367..-3960.1176] | it/evals=3760/5198 eff=78.3660% N=400 

	Mono-modal Volume: ~exp(-13.87) * Expected Volume: exp(-9.45) Quality: ok

	mypow.gamma:     +1.00|*  +1.01                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00016  *******  +0.00023                                                      | +0.00100

	Z=-3996.8(0.00%) | Like=-3984.08..-3952.58 [-4332.2367..-3960.1176] | it/evals=3840/5301 eff=78.3514% N=400 

	Mono-modal Volume: ~exp(-13.87)   Expected Volume: exp(-9.68) Quality: ok

	mypow.gamma:     +1.00|*  +1.01                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00016  *******  +0.00023                                                      | +0.00100

	Z=-3992.8(0.00%) | Like=-3980.03..-3952.58 [-4332.2367..-3960.1176] | it/evals=3920/5399 eff=78.4157% N=400 

	Mono-modal Volume: ~exp(-14.11) * Expected Volume: exp(-9.90) Quality: ok

	mypow.gamma:     +1.00|*  +1.01                                                                            |    +3.00
	mypow.ampl :  +0.00000|   +0.00017  *******  +0.00023                                                      | +0.00100

	Z=-3987.3(0.00%) | Like=-3974.22..-3952.51 [-4332.2367..-3960.1176] | it/evals=4040/5556 eff=78.3553% N=400 

	Mono-modal Volume: ~exp(-14.41) * Expected Volume: exp(-10.13) Quality: ok

	mypow.gamma:    +1.000|*  +1.009                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00017  *****  +0.00023                                                       | +0.00100

	Z=-3984.3(0.00%) | Like=-3971.47..-3952.37 [-4332.2367..-3960.1176] | it/evals=4120/5648 eff=78.5061% N=400 

	Mono-modal Volume: ~exp(-14.89) * Expected Volume: exp(-10.35) Quality: ok

	mypow.gamma:    +1.000|*  +1.008                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00017  *****  +0.00022                                                       | +0.00100

	Z=-3982.3(0.00%) | Like=-3969.36..-3951.35 [-4332.2367..-3960.1176] | it/evals=4200/5743 eff=78.6075% N=400 

	Mono-modal Volume: ~exp(-15.11) * Expected Volume: exp(-10.58) Quality: ok

	mypow.gamma:    +1.000|*  +1.007                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00017  *****  +0.00022                                                       | +0.00100

	Z=-3979.9(0.00%) | Like=-3966.89..-3951.35 [-4332.2367..-3960.1176] | it/evals=4280/5843 eff=78.6331% N=400 

	Mono-modal Volume: ~exp(-15.24) * Expected Volume: exp(-10.80) Quality: ok

	mypow.gamma:    +1.000|*  +1.006                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00017  *****  +0.00022                                                       | +0.00100

	Z=-3977.2(0.00%) | Like=-3964.12..-3951.35 [-4332.2367..-3960.1176] | it/evals=4400/5995 eff=78.6416% N=400 

	Mono-modal Volume: ~exp(-15.60) * Expected Volume: exp(-11.02) Quality: ok

	mypow.gamma:    +1.000|*  +1.005                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00018  *****  +0.00022                                                       | +0.00100

	Z=-3975.7(0.00%) | Like=-3962.51..-3951.35 [-4332.2367..-3960.1176] | it/evals=4480/6095 eff=78.6655% N=400 

	Mono-modal Volume: ~exp(-15.60)   Expected Volume: exp(-11.25) Quality: ok

	mypow.gamma:    +1.000|*  +1.004                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00018  *****  +0.00021                                                       | +0.00100

	Z=-3974.4(0.02%) | Like=-3961.31..-3951.35 [-4332.2367..-3960.1176] | it/evals=4560/6190 eff=78.7565% N=400 

	Mono-modal Volume: ~exp(-15.84) * Expected Volume: exp(-11.47) Quality: ok

	mypow.gamma:    +1.000|*  +1.004                                                                           |   +3.000
	mypow.ampl :  +0.00000|    +0.00018  *****  +0.00021                                                       | +0.00100

	Z=-3973.3(0.05%) | Like=-3959.89..-3951.35 [-3960.1160..-3955.4960] | it/evals=4640/6292 eff=78.7508% N=400 

	Mono-modal Volume: ~exp(-16.06) * Expected Volume: exp(-11.70) Quality: ok

	mypow.gamma:    +1.000|*  +1.003                                                                           |   +3.000
	mypow.ampl :  +0.00000|     +0.00018  ***  +0.00021                                                        | +0.00100

	Z=-3971.7(0.27%) | Like=-3958.30..-3951.35 [-3960.1160..-3955.4960] | it/evals=4760/6446 eff=78.7297% N=400 

	Mono-modal Volume: ~exp(-16.14) * Expected Volume: exp(-11.92) Quality: ok

	mypow.gamma:    +1.000|*  +1.003                                                                           |   +3.000
	mypow.ampl :  +0.00000|     +0.00018  ***  +0.00021                                                        | +0.00100

	Z=-3970.9(0.62%) | Like=-3957.53..-3951.35 [-3960.1160..-3955.4960] | it/evals=4840/6550 eff=78.6992% N=400 

	Mono-modal Volume: ~exp(-16.38) * Expected Volume: exp(-12.15) Quality: ok

	mypow.gamma:    +1.000|*  +1.002                                                                           |   +3.000
	mypow.ampl :  +0.00000|     +0.00018  ***  +0.00021                                                        | +0.00100

	Z=-3970.2(1.22%) | Like=-3956.76..-3951.30 [-3960.1160..-3955.4960] | it/evals=4920/6660 eff=78.5942% N=400 

	Mono-modal Volume: ~exp(-16.72) * Expected Volume: exp(-12.37) Quality: ok

	mypow.gamma:    +0.000|                    +1.000  *  +1.002                                               |   +3.000
	mypow.ampl :  +0.00000|     +0.00018  ***  +0.00021                                                        | +0.00100

	Z=-3969.6(2.35%) | Like=-3955.93..-3951.30 [-3960.1160..-3955.4960] | it/evals=5000/6775 eff=78.4314% N=400 

	Mono-modal Volume: ~exp(-17.25) * Expected Volume: exp(-12.60) Quality: ok

	mypow.gamma:    +0.000|                    +1.000  *  +1.002                                               |   +3.000
	mypow.ampl :  +0.00000|     +0.00018  ***  +0.00021                                                        | +0.00100

	Z=-3968.8(5.38%) | Like=-3955.08..-3951.30 [-3955.4918..-3954.6140] | it/evals=5120/6924 eff=78.4795% N=400 

	Mono-modal Volume: ~exp(-17.25)   Expected Volume: exp(-12.82) Quality: ok

	mypow.gamma:    +0.000|                    +1.000  *  +1.002                                               |   +3.000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00021                                                        | +0.00100

	Z=-3968.3(8.30%) | Like=-3954.58..-3951.30 [-3954.6009..-3954.3077] | it/evals=5200/7020 eff=78.5498% N=400 

	Mono-modal Volume: ~exp(-17.51) * Expected Volume: exp(-13.05) Quality: ok

	mypow.gamma:    +0.000|                    +1.000  *  +1.001                                               |   +3.000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00021                                                        | +0.00100

	Z=-3967.9(11.79%) | Like=-3954.15..-3951.22 [-3954.1493..-3954.1122]*| it/evals=5280/7115 eff=78.6299% N=400 

	Mono-modal Volume: ~exp(-17.51)   Expected Volume: exp(-13.27) Quality: ok

	mypow.gamma:    +0.000|                    +1.000  *  +1.001                                               |   +3.000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00021                                                        | +0.00100

	Z=-3967.6(15.79%) | Like=-3953.76..-3951.22 [-3953.7594..-3953.7505]*| it/evals=5360/7222 eff=78.5693% N=400 

	Mono-modal Volume: ~exp(-17.51)   Expected Volume: exp(-13.50) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0009                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00021                                                        | +0.00100

	Z=-3967.2(23.68%) | Like=-3953.35..-3951.20 [-3953.3535..-3953.3517]*| it/evals=5480/7383 eff=78.4763% N=400 

	Mono-modal Volume: ~exp(-17.72) * Expected Volume: exp(-13.72) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0008                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00021                                                        | +0.00100

	Z=-3967.0(29.66%) | Like=-3953.09..-3951.20 [-3953.0935..-3953.0875]*| it/evals=5560/7488 eff=78.4424% N=400 

	Mono-modal Volume: ~exp(-17.77) * Expected Volume: exp(-13.95) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0007                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00020                                                        | +0.00100

	Z=-3966.8(35.86%) | Like=-3952.84..-3951.20 [-3952.8395..-3952.8389]*| it/evals=5640/7605 eff=78.2790% N=400 

	Mono-modal Volume: ~exp(-18.67) * Expected Volume: exp(-14.17) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0006                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00020                                                        | +0.00100

	Z=-3966.6(41.40%) | Like=-3952.65..-3951.20 [-3952.6483..-3952.6472]*| it/evals=5720/7707 eff=78.2811% N=400 

	Mono-modal Volume: ~exp(-18.74) * Expected Volume: exp(-14.40) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0005                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00020                                                        | +0.00100

	Z=-3966.4(50.27%) | Like=-3952.37..-3951.20 [-3952.3749..-3952.3708]*| it/evals=5840/7863 eff=78.2527% N=400 

	Mono-modal Volume: ~exp(-19.02) * Expected Volume: exp(-14.62) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0004                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  ***  +0.00020                                                        | +0.00100

	Z=-3966.3(55.97%) | Like=-3952.23..-3951.20 [-3952.2308..-3952.2306]*| it/evals=5920/7964 eff=78.2655% N=400 

	Mono-modal Volume: ~exp(-19.43) * Expected Volume: exp(-14.85) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0004                                              |  +3.0000
	mypow.ampl :  +0.00000|     +0.00019  **  +0.00020                                                         | +0.00100

	Z=-3966.2(61.58%) | Like=-3952.10..-3951.20 [-3952.0962..-3952.0958]*| it/evals=6000/8072 eff=78.2065% N=400 

	Mono-modal Volume: ~exp(-19.43)   Expected Volume: exp(-15.07) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0003                                              |  +3.0000
	mypow.ampl :  +0.00000|      +0.00019  *  +0.00020                                                         | +0.00100

	Z=-3966.2(66.61%) | Like=-3952.01..-3951.20 [-3952.0116..-3952.0102]*| it/evals=6080/8165 eff=78.3001% N=400 

	Mono-modal Volume: ~exp(-19.43)   Expected Volume: exp(-15.30) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0003                                              |  +3.0000
	mypow.ampl :  +0.00000|      +0.00019  *  +0.00020                                                         | +0.00100

	Z=-3966.1(73.10%) | Like=-3951.87..-3951.19 [-3951.8731..-3951.8724]*| it/evals=6200/8321 eff=78.2729% N=400 

	Mono-modal Volume: ~exp(-19.90) * Expected Volume: exp(-15.52) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0003                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000192  *  +0.000201                                                        |+0.001000

	Z=-3966.0(76.93%) | Like=-3951.78..-3951.19 [-3951.7795..-3951.7786]*| it/evals=6280/8429 eff=78.2165% N=400 

	Mono-modal Volume: ~exp(-19.92) * Expected Volume: exp(-15.75) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0002                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000192  *  +0.000201                                                        |+0.001000

	Z=-3966.0(80.28%) | Like=-3951.70..-3951.19 [-3951.7006..-3951.6982]*| it/evals=6360/8528 eff=78.2480% N=400 

	Mono-modal Volume: ~exp(-19.92)   Expected Volume: exp(-15.97) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0002                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000192  *  +0.000200                                                        |+0.001000

	Z=-3965.9(83.26%) | Like=-3951.64..-3951.19 [-3951.6351..-3951.6350]*| it/evals=6440/8641 eff=78.1459% N=400 

	Mono-modal Volume: ~exp(-20.54) * Expected Volume: exp(-16.20) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0002                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000193  *  +0.000200                                                        |+0.001000

	Z=-3965.9(86.99%) | Like=-3951.55..-3951.19 [-3951.5503..-3951.5501]*| it/evals=6560/8798 eff=78.1138% N=400 

	Mono-modal Volume: ~exp(-20.90) * Expected Volume: exp(-16.42) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0001                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000193  *  +0.000200                                                        |+0.001000

	Z=-3965.9(89.05%) | Like=-3951.50..-3951.19 [-3951.4965..-3951.4961]*| it/evals=6640/8897 eff=78.1452% N=400 

	Mono-modal Volume: ~exp(-21.13) * Expected Volume: exp(-16.65) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0001                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000193  *  +0.000200                                                        |+0.001000

	Z=-3965.9(90.83%) | Like=-3951.46..-3951.19 [-3951.4617..-3951.4610]*| it/evals=6720/9004 eff=78.1032% N=400 

	Mono-modal Volume: ~exp(-21.13)   Expected Volume: exp(-16.87) Quality: ok

	mypow.gamma:   +0.0000|                   +1.0000  *  +1.0001                                              |  +3.0000
	mypow.ampl : +0.000000|     +0.000193  *  +0.000200                                                        |+0.001000

	Z=-3965.8(92.35%) | Like=-3951.43..-3951.19 [-3951.4309..-3951.4307]*| it/evals=6800/9108 eff=78.0891% N=400 

	Mono-modal Volume: ~exp(-21.61) * Expected Volume: exp(-17.10) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00009                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000199                                                        |+0.001000

	Z=-3965.8(94.18%) | Like=-3951.38..-3951.19 [-3951.3804..-3951.3803]*| it/evals=6920/9262 eff=78.0862% N=400 

	Mono-modal Volume: ~exp(-21.81) * Expected Volume: exp(-17.32) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00007                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000199                                                        |+0.001000

	Z=-3965.8(95.16%) | Like=-3951.36..-3951.19 [-3951.3589..-3951.3586]*| it/evals=7000/9359 eff=78.1337% N=400 

	Mono-modal Volume: ~exp(-21.81)   Expected Volume: exp(-17.55) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00006                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000199                                                        |+0.001000

	Z=-3965.8(95.98%) | Like=-3951.34..-3951.19 [-3951.3359..-3951.3357]*| it/evals=7080/9462 eff=78.1284% N=400 

	Mono-modal Volume: ~exp(-21.81)   Expected Volume: exp(-17.77) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00005                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000199                                                        |+0.001000

	Z=-3965.8(96.68%) | Like=-3951.32..-3951.19 [-3951.3161..-3951.3161]*| it/evals=7160/9562 eff=78.1489% N=400 

	Mono-modal Volume: ~exp(-22.78) * Expected Volume: exp(-18.00) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00005                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000198                                                        |+0.001000

	Z=-3965.8(97.50%) | Like=-3951.29..-3951.19 [-3951.2918..-3951.2912]*| it/evals=7280/9727 eff=78.0530% N=400 

	Mono-modal Volume: ~exp(-22.78)   Expected Volume: exp(-18.23) Quality: ok

	mypow.gamma:  +0.00000|                  +1.00000  *  +1.00004                                             | +3.00000
	mypow.ampl : +0.000000|     +0.000194  *  +0.000198                                                        |+0.001000

	[ultranest] Explored until L=-4e+03  
	[ultranest] Likelihood function evaluations: 9743
	[ultranest] Writing samples and results to disk ...
	[ultranest] Writing samples and results to disk ... done
	[ultranest]   logZ = -3966 +- 0.1345
	[ultranest] Posterior uncertainty strategy is satisfied (KL: 0.46+-0.08 nat, need <0.50 nat)
	[ultranest] Evidency uncertainty strategy is satisfied (dlogz=0.27, need <0.5)
	[ultranest]   logZ error budget: single: 0.18 bs:0.13 tail:0.02 total:0.14 required:<0.50
	[ultranest] done iterating.

	logZ = -3965.760 +- 0.273
	  single instance: logZ = -3965.760 +- 0.181
	  bootstrapped   : logZ = -3965.761 +- 0.272
	  tail           : logZ = +- 0.024
	insert order U test : converged: True correlation: inf iterations

		mypow.gamma         1.00039 +- 0.00038
		mypow.ampl          0.0001964 +- 0.0000045


Output files::

	$ find simplest-/
	simplest-/
	simplest-/debug.log
	simplest-/plots
	simplest-/plots/corner.pdf
	simplest-/plots/trace.pdf
	simplest-/plots/run.pdf
	simplest-/info
	simplest-/info/post_summary.csv
	simplest-/info/results.json
	simplest-/results
	simplest-/results/points.hdf5
	simplest-/extra
	simplest-/chains
	simplest-/chains/run.txt
	simplest-/chains/weighted_post_untransformed.txt
	simplest-/chains/equal_weighted_post.txt
	simplest-/chains/weighted_post.txt


"simplest-" is the `outputfiles_basename` defined in the script.

The most important files are:

* plots/corner.pdf:

	.. image:: reference-output/corner.png
	
	Plot of the parameter constraints and uncertainties and their correlations.
	The photon index parameter is hitting the edge of the parameter space,
	and its uncertainties are tiny. This can be a hint that it is a poor model.

* info/results.json: summary of all parameters, their uncertainties and estimated lnZ
* info/post_summary.csv: summary of all parameters and their uncertainties as CSV
* chains/equal_weighted_post.txt: contains posterior samples: each row is a model parameter vector. You can iterate through these, set up the model in pyxspec, and then do something with it (compute fluxes and luminosities, for example).

You probably want to plot the fit as well (after setting to the best fit).

Try modifying the model.

For more information, see https://johannesbuchner.github.io/BXA/sherpa-analysis.html

Other examples
---------------

* Example of PCA background:

  This uses the Swift data file swift/interval0pc.pi.

  First, store the galactic NH value (1.68e+20)
  into the text file swift/interval0pc.pi.nh.

  Then run with::

	$ python3 example_pcabackground.py
	
	....

	loading nH from swift/interval0pc.pi.nh (expecting something like 1e21 in there)
	setting galactic nH to 0.0168 [units of 1e22/cm²]
	[bxa.Fitter INFO]: PCAFitter(for ID=2)
	[bxa.Fitter INFO]: loading PCA information from /home/user/bin/ciao-4.13/ots/lib/python3.7/site-packages/bxa/sherpa/background/swift_xrt_1024.json
	[bxa.Fitter INFO]: fitting background of ID=2 using PCA method
	[bxa.Fitter INFO]: have 2751 background counts for deconvolution
	[bxa.Fitter INFO]: fit: initial PCA decomposition: [ 3.43950638e+00 -2.18629410e-02  7.52750306e-03 -4.07883039e-03
	 -3.49918117e-03 -3.20861431e-03  3.52942831e-03 -5.05089198e-03
	 -9.34656625e-04 -4.86905140e-03  2.29800943e-03]
	[bxa.Fitter INFO]: fit: first full fit done
	[bxa.Fitter INFO]: fit: parameters: [-0.8139987117963805, 0.42489117817206506, 0.03088268390136437, 0.19696313135650556, 0.09137494506325541, -0.17493295963368954, -0.09507225292526847, 0.16435598097773643, -0.058544963240419884, 0.25546836854960586, 0.08241814841520864]
	[bxa.Fitter INFO]: fit: stat: 551.3592848191211
	[bxa.Fitter INFO]: fit: second full fit from zero
	[bxa.Fitter INFO]: fit: parameters: [-0.8139987117963805, 0.42489117817206506, 0.03088268390136437, 0.19696313135650556, 0.09137494506325541, -0.17493295963368954, -0.09507225292526847, 0.16435598097773643, -0.058544963240419884, 0.25546836854960586, 0.08241814841520864]
	[bxa.Fitter INFO]: fit: stat: 551.3592848191096
	[bxa.Fitter INFO]: fit: using zero-fit
	11 parameters, stat=551.36
	--> 10 parameters, stat=552.44
	--> 9 parameters, stat=582.41
	--> 8 parameters, stat=583.58
	--> 7 parameters, stat=682.99
	--> 6 parameters, stat=696.18
	--> 5 parameters, stat=698.64
	--> 4 parameters, stat=707.46
	--> 3 parameters, stat=716.11
	--> 2 parameters, stat=716.63
	--> 1 parameters, stat=1145.24

	Background PCA fitting AIC results:
	-----------------------------------

	stat Ncomp AIC
	1145.2  1 1147.2
	716.6  2 720.6
	716.1  3 722.1
	707.5  4 715.5
	698.6  5 708.6
	696.2  6 708.2
	683.0  7 697.0
	583.6  8 599.6
	582.4  9 600.4
	552.4 10 572.4
	551.4 11 573.4

	Increasing parameters again...
	11 parameters, aic=573.36
	Final choice: 10 parameters, aic=572.44

	Adding Gaussian#1
	largest remaining discrepancy at 1.855keV[185], need 5959 counts
	placing gaussian at 1.86keV, with power 0.6227582993302021
	with Gaussian: 579.3593637901457 ; change: 6.9 (negative is good)
	not significant, rejecting
	creating prior functions...
	running BXA ...
	[ultranest] Sampling 400 live points from prior ...
	[ultranest INFO]: Sampling 400 live points from prior ...


	Mono-modal Volume: ~exp(-2.94) * Expected Volume: exp(0.00) Quality: ok

	src.level   :      -8.0|*************************************************** ******************************************|     +2.0
	src.PhoIndex:      +1.0|      *                *   * *********************************                               *|     +3.0
	src.nh      :     +19.0|**********************************************************************************************|    +24.0
	src.redshift:      +0.0|                +0.2  * ********************************* *  +0.4                             |     +0.7
	pca2.lognorm:      -5.0|******************************************************************************************** *|    +20.0

	Z=-1e+18(0.00%) | Like=-1.3e+18..-7.9e+02 [-1.142e+23..-1100] | it/evals=88/505 eff=83.8095% N=400 

	....

	logZ = -356.466 +- 0.259
	  single instance: logZ = -356.466 +- 0.183
	  bootstrapped   : logZ = -356.488 +- 0.259
	  tail           : logZ = +- 0.011

		src.level           -2.404 +- 0.037
		src.PhoIndex        2.029 +- 0.045
		src.nh              19.75 +- 0.46
		src.redshift        0.302 +- 0.050
		pca2.lognorm        -0.769 +- 0.031



* Example of empirical background model (and different priors). Redshift is a free parameter here:
  
  Run with::

    $ python3 example_automatic_background_model.py

	calling singlefitter...
	[bxa.Fitter INFO]: SingleFitter(for ID=2, storing to "swift/interval0pc")
	[bxa.Fitter INFO]: prepare_stage 2 of ID=2
	[bxa.Fitter INFO]: prepare_stage 2 of ID=2 done
	[bxa.Fitter INFO]: prepare_stage 2 of ID=2
	[bxa.Fitter INFO]: prepare_stage 2 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 2 of ID=2
	[bxa.Fitter INFO]: fit_stage 2 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 2 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 2 of ID=2.  stage done
	[bxa.Fitter INFO]: prepare_stage 3 of ID=2
	[bxa.Fitter INFO]: prepare_stage 3 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 3 of ID=2
	[bxa.Fitter INFO]: fit_stage 3 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 3 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 3 of ID=2.  stage done
	[bxa.Fitter INFO]: prepare_stage 4 of ID=2
	[bxa.Fitter INFO]: prepare_stage 4 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 4 of ID=2
	[bxa.Fitter INFO]: fit_stage 4 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 4 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 4 of ID=2.  stage done
	[bxa.Fitter INFO]: prepare_stage 5 of ID=2
	[bxa.Fitter INFO]: prepare_stage 5 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 5 of ID=2
	[bxa.Fitter INFO]: fit_stage 5 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 5 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 5 of ID=2.  stage done
	[bxa.Fitter INFO]: prepare_stage 6 of ID=2
	[bxa.Fitter INFO]: prepare_stage 6 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 6 of ID=2
	[bxa.Fitter INFO]: fit_stage 6 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 6 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 6 of ID=2.  stage done
	[bxa.Fitter INFO]: prepare_stage 7 of ID=2
	[bxa.Fitter INFO]: prepare_stage 7 of ID=2 done
	[bxa.Fitter INFO]: fit_stage 7 of ID=2
	[bxa.Fitter INFO]: fit_stage 7 of ID=2.  fine fit ...
	[bxa.Fitter INFO]: fit_stage 7 of ID=2.  fitted
	[bxa.Fitter INFO]: fit_stage 7 of ID=2.  stage done
	[bxa.Fitter INFO]: Background fit complete.

	freezing background params
	loading nH from swift/interval0pc.pi.nh (expecting something like 1e21 in there)
	setting galactic nH to 0.0168 [units of 1e22/cm²]
	apply_rmf(apply_arf((9504.67 * (((xszpowerlw.src * xszwabs.abso) * xswabs.galabso) + (0.05553358353932973 * (((1.0 - box1d.dip_2) * (((xsbknpower.pbknpl_2 + gauss1d.gauss1_2) + gauss1d.gauss2_2) + gauss1d.gauss3_2)) + gauss1d.gauss4_2))))))
	   Param        Type          Value          Min          Max      Units
	   -----        ----          -----          ---          ---      -----
	   src.PhoIndex thawed            1           -2            9           
	   src.redshift frozen            0       -0.999           10           
	   src.norm     thawed            1            0        1e+24           
	   abso.nH      thawed            1            0       100000 10^22 atoms / cm^2
	   abso.redshift frozen            0       -0.999           10           
	   galabso.nH   thawed            1            0       100000 10^22 atoms / cm^2
	   dip_2.xlow   frozen      1.76028         1.75         2.25           
	   dip_2.xhi    frozen      3.21915         2.75         3.25           
	   dip_2.ampl   frozen     0.937426        0.001        0.999           
	   pbknpl_2.PhoIndx1 frozen      1.53607          0.8            4           
	   pbknpl_2.BreakE frozen      4.10858          0.2            5        keV
	   pbknpl_2.PhoIndx2 frozen      2.71323          0.8            4           
	   pbknpl_2.norm frozen    0.0140262        1e-10            1           
	   gauss1_2.fwhm frozen     0.817724         0.01            1           
	   gauss1_2.pos frozen     0.607138          0.1          1.1           
	   gauss1_2.ampl frozen      0.01218        1e-06            1           
	   gauss2_2.fwhm frozen    0.0199524         0.01            1           
	   gauss2_2.pos frozen      2.19733            2          2.5           
	   gauss2_2.ampl frozen   0.00376109        1e-06            1           
	   gauss3_2.fwhm frozen    0.0303601         0.01            1           
	   gauss3_2.pos frozen      1.37532            1          1.4           
	   gauss3_2.ampl frozen   0.00415471        1e-06            1           
	   gauss4_2.fwhm frozen      0.93892         0.01            1           
	   gauss4_2.pos frozen        0.125            0          0.5           
	   gauss4_2.ampl frozen  0.000819252        1e-06            1           
	creating prior functions...
	running BXA ...
	[ultranest] Sampling 400 live points from prior ...
	[ultranest INFO]: Sampling 400 live points from prior ...


	Mono-modal Volume: ~exp(-3.96) * Expected Volume: exp(0.00) Quality: ok

	src.level   :      -8.0|*************************** **** ***** *************** **************************************** *************************************** |     -1.0
	src.PhoIndex:      +1.0|***** *************************** **************************************** ** ******* ********  ********** ************************* ***|     +3.0
	src.nh      :     +19.0|**************************************************** **************************** **** *********************************************** *|    +24.0
	gal.nh      :     +20.8|                                        +21.6  ************************************** * ***  +22.4                                      |    +23.2

	Z=-10118.6(0.00%) | Like=-10114.95..-778.63 [-34687.3038..-1560.0087] | it/evals=80/491 eff=87.9121% N=400 
	...

	logZ = -362.501 +- 0.288
	  single instance: logZ = -362.501 +- 0.221
	  bootstrapped   : logZ = -362.528 +- 0.288
	  tail           : logZ = +- 0.010

		src.level           -2.394 +- 0.016
		src.PhoIndex        2.200 +- 0.046
		src.nh              19.51 +- 0.35
		gal.nh              20.807 +- 0.035


Compare the models with::

	$ python3 model_compare.py superfit/ wabs_noz/

	Model comparison
	****************

	model superfit/ : log10(Z) =    -2.6  XXX ruled out
	model wabs_noz/ : log10(Z) =     0.0    <-- GOOD

	The last, most likely model was used as normalization.
	Uniform model priors are assumed, with a cut of log10(30) to rule out models.


Beware of the caveats of these log10(Z) differences (log-Bayes factors),
and derive thresholds with simulated data.

For the full documentation, see https://johannesbuchner.github.io/BXA/sherpa-analysis.html

Please explore this folder for other demo scripts.

For example, go into the chandra folder, and run the `xagnfitter.py <https://johannesbuchner.github.io/BXA/xagnfitter.html>`_ in this folder against it.
BXA/Xspec example scripts
==========================

This folder contains simple and complex examples
how BXA can be invoked for spectral analysis.

Please refer to https://johannesbuchner.github.io/BXA/
for full documentation, including how to install BXA.


Generate test data
-------------------

Generate test data with xspec using the commands in gen.xspec
This will produce a spectral file example-file.fak
representing a ATHENA observation of an absorbed AGN.

Expected output::

	$ < gen.xspec xspec

The spectrum looks something like this:

.. image:: reference-output/data.gif

Simple analysis
-----------------

Have a look at the file example_simplest.py. It contains:

* Loading data
* setting up a model and its parameter ranges
* running a BXA fit with specified priors
* plotting the posterior predictions (convolved with the response)
* plotting the model (posterior predictions, not convolved)
* making a Q-Q plot

See https://johannesbuchner.github.io/BXA/ to understand the code.
See https://johannesbuchner.github.io/UltraNest/ to understand the output of the
fitting engine (for example, its `FAQ page <https://johannesbuchner.github.io/UltraNest/issues.html>`_).

Expected output::

	$ python3 example_simplest.py
	Default fit statistic is set to: C-Statistic
	   This will apply to all current and newly loaded spectra.

	1 spectrum  in use
	 
	Spectral Data File: example-file.fak  Spectrum 1
	Net count rate (cts/s) for Spectrum:1  4.224e+00 +/- 9.191e-02
	 Assigned to Data Group 1 and Plot Group 1
	  Noticed Channels:  1-4096
	  Telescope: ATHENA+ Instrument: WFI  Channel Type: PI
	  Exposure Time: 500 sec
	 Using fit statistic: cstat
	 Using Response (RMF) File            athenapp_ir_b4c_wfi_withfilter_fov40.0arcmin_avg.rsp for Source 1

	  4096 channels (1,4096) ignored in spectrum #     1

	   801 channels (11-811) noticed in spectrum #     1


	========================================================================
	Model powerlaw<1> Source No.: 1   Active/On
	Model Model Component  Parameter  Unit     Value
	 par  comp
	   1    1   powerlaw   PhoIndex            1.00000      +/-  0.0          
	   2    1   powerlaw   norm                1.00000      +/-  0.0          
	________________________________________________________________________


	Fit statistic  : C-Statistic              2.056191e+07     using 801 bins.

	Test statistic : Chi-Squared              4.197601e+11     using 801 bins.

	***Warning: Chi-square may not be valid due to bins with zero variance
				in spectrum number(s): 1 

	 Null hypothesis probability of 0.000000e+00 with 799 degrees of freedom
	 Current data and model not fit yet.

	Fit statistic  : C-Statistic              2.056191e+07     using 801 bins.

	Test statistic : Chi-Squared              4.197601e+11     using 801 bins.

	***Warning: Chi-square may not be valid due to bins with zero variance
				in spectrum number(s): 1 

	 Null hypothesis probability of 0.000000e+00 with 799 degrees of freedom
	 Current data and model not fit yet.

	Fit statistic  : C-Statistic              2.056191e+07     using 801 bins.

	Test statistic : Chi-Squared              4.197601e+11     using 801 bins.

	***Warning: Chi-square may not be valid due to bins with zero variance
				in spectrum number(s): 1 

	 Null hypothesis probability of 0.000000e+00 with 799 degrees of freedom
	 Current data and model not fit yet.
	  uniform prior for PhoIndex between 1.000000 and 3.000000 
	  jeffreys prior for norm between 1.000000e-10 and 1.000000e+01 
	   note: this parameter spans *many* dex. Double-check the limits are reasonable.
	running analysis ...
	[ultranest] Resuming from 7774 stored points


	Mono-modal Volume: ~exp(-4.24) * Expected Volume: exp(0.00) Quality: ok

	PhoIndex :      +1.0|*** ****************************** ****************************************************** **************|     +3.0
	log(norm):     -10.0|********************************************************************************************************|     +1.0

	Z=-1199206.7(0.00%) | Like=-1089578.51..-4277.72 [-1.045e+08..-4464] | it/evals=80/9998 eff=inf% N=400 

	Mono-modal Volume: ~exp(-4.24)   Expected Volume: exp(-0.23) Quality: correlation length: 3 (+)

	PhoIndex :      +1.0|********************************** ************* *******************************************************|     +3.0
	log(norm):     -10.0|************************************************************************************  -1.2              |     +1.0

	Z=-24616.1(0.00%) | Like=-24611.37..-4277.72 [-1.045e+08..-4464] | it/evals=160/9998 eff=inf% N=400 0 

	...
	...
	...

	Mono-modal Volume: ~exp(-22.05) * Expected Volume: exp(-18.00) Quality: correlation length: 1913 (+)

	PhoIndex :  +0.00000|                        +1.00000  *  +1.00005                                                           | +3.00000
	log(norm):   -10.000|                                                   -3.709  *  -3.699                                    |   +1.000

	Z=-3996.5(96.93%) | Like=-3981.81..-3981.70 [-3981.8147..-3981.8146]*| it/evals=7280/9998 eff=inf% N=400 

	Mono-modal Volume: ~exp(-22.35) * Expected Volume: exp(-18.23) Quality: correlation length: 1913 (+)

	PhoIndex :  +0.00000|                        +1.00000  *  +1.00004                                                           | +3.00000
	log(norm):   -10.000|                                                   -3.709  *  -3.700                                    |   +1.000

	[ultranest] Explored until L=-4e+03  981.70 [-3981.7988..-3981.7987]*| it/evals=7360/9998 eff=inf% N=400 
	[ultranest] Likelihood function evaluations: 9998
	[ultranest] Writing samples and results to disk ...
	[ultranest] Writing samples and results to disk ... done
	[ultranest]   logZ = -3996 +- 0.1528
	[ultranest] Posterior uncertainty strategy is satisfied (KL: 0.46+-0.08 nat, need <0.50 nat)
	[ultranest] Evidency uncertainty strategy is satisfied (dlogz=0.39, need <0.5)
	[ultranest]   logZ error budget: single: 0.18 bs:0.15 tail:0.02 total:0.15 required:<0.50
	[ultranest] done iterating.

	logZ = -3996.484 +- 0.389
	  single instance: logZ = -3996.484 +- 0.183
	  bootstrapped   : logZ = -3996.490 +- 0.389
	  tail           : logZ = +- 0.024
	insert order U test : converged: False correlation: 3.0 iterations

		PhoIndex            1.00038 +- 0.00038
		log(norm)           -3.7043 +- 0.0094
	running analysis ... done!
	creating plot of posterior predictions against data ...
	100%|████████████████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:00<00:00, 107.90it/s]
	binning for plot...
	100%|█████████████████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:01<00:00, 85.53it/s]
	saving plot...
	creating plot of posterior predictions ...
	100%|████████████████████████████████████████████████████████████████████████████████████████████████| 100/100 [00:00<00:00, 117.24it/s]
	saving plot...
	creating quantile-quantile plot ...
	saving plot...


Output files::

	$ find simplest/
	simplest/
	simplest/debug.log
	simplest/convolved_posterior.pdf
	simplest/chain.fits
	simplest/plots
	simplest/plots/corner.pdf
	simplest/plots/trace.pdf
	simplest/plots/run.pdf
	simplest/unconvolved_posterior.pdf
	simplest/info
	simplest/info/post_summary.csv
	simplest/info/results.json
	simplest/qq_model_deviations.pdf
	simplest/results
	simplest/results/points.hdf5
	simplest/extra
	simplest/chains
	simplest/chains/run.txt
	simplest/chains/weighted_post_untransformed.txt
	simplest/chains/equal_weighted_post.txt
	simplest/chains/weighted_post.txt

"simplest/" is the `outputfiles_basename` defined in the script.

The most important files are:

* unconvolved_posterior.pdf : 

	.. image:: reference-output/unconvolved_posterior.png
	
	The model itself is a powerlaw, and the uncertainties are too narrow to see.

	For further explanation of this plot, see https://johannesbuchner.github.io/BXA/xspec-analysis.html

* convolved_posterior.pdf : 

	.. image:: reference-output/convolved_posterior.png
	
	The model and the data convolved through the response. 
	Red means the data are poorly fitted by this model.
	The model is clearly off -- For example, the lower energy X-rays are overpredicted.

	For further explanation of this plot, see https://johannesbuchner.github.io/BXA/xspec-analysis.html

* plots/corner.pdf:

	.. image:: reference-output/corner.png
	
	Plot of the parameter constraints and uncertainties and their correlations.
	The photon index parameter is hitting the edge of the parameter space,
	and its uncertainties are tiny. Another hint of a poor model.

	For further explanation of this plot, see https://johannesbuchner.github.io/BXA/xspec-analysis.html

* qq_model_deviations.pdf : 
	
	.. image:: reference-output/qq_model_deviations.png
	
	`Q-Q plot <https://en.wikipedia.org/wiki/Q%E2%80%93Q_plot>`_:
	The red curve is far from the 1:1 line. That it is on the bottom right
	indicates the model produces many more counts than the data.
	The tickmarks indicate that the problem is accumulating below 2keV.

	For further explanation of this plot, see https://johannesbuchner.github.io/BXA/xspec-analysis.html

* info/results.json: summary of all parameters, their uncertainties and estimated lnZ
* info/post_summary.csv: summary of all parameters and their uncertainties as CSV
* chains/equal_weighted_post.txt: contains posterior samples: each row is a model parameter vector. You can iterate through these, set up the model in pyxspec, and then do something with it (compute fluxes and luminosities, for example).

Other examples
---------------

* example_advanced_priors.py shows a absorbed powerlaw fit, which is better. It 
  also demonstrates how to specify custom prior functions.

  Run with::
	
	$ python3 example_advanced_priors.py example-file.fak absorbed/
	
  Here the spectral file and output folder are command line arguments,
  which is convenient for analysing many sources.

* example_custom_run.py finally adds a emission line. Run with::

	$ python3 example_custom_run.py example-file.fak line/

Compare the models with::

	$ python3 model_compare.py absorbed simplest line

	Model comparison
	****************

	model simplest  : log10(Z) = -1519.1  XXX ruled out
	model absorbed  : log10(Z) =    -5.6  XXX ruled out
	model line      : log10(Z) =     0.0    <-- GOOD

	The last, most likely model was used as normalization.
	Uniform model priors are assumed, with a cut of log10(30) to rule out models.

Beware of the caveats of these log10(Z) differences (log-Bayes factors),
and derive thresholds with simulated data. 

For the full documentation, see https://johannesbuchner.github.io/BXA/xspec-analysis.html

Please explore this folder for other demo scripts.
.. include:: ../CONTRIBUTING.rst
.. include:: ../HISTORY.rst
Obscured Active Galactic Nuclei
=======================================

A script for fitting Active Galactic Nuclei is provided at
`xagnfitter.py <https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xagnfitter.py>`_.
This is the method used in Buchner+14, Buchner+15, Simmonds+17.

Features:

* Maximum information extraction in the low count regime, by Bayesian inference and background models.
* Provides robust uncertainty estimation of all parameters, including the

  * obscuring column density NH
  * Photon index Gamma
  * rest-frame, intrinsic accretion luminosity
  * etc.

* Redshift can be fixed, unknown or come from a probability distribution (photo-z)
* Realistic nuclear obscurer model (UXCLUMPY) that 

  * extends to the highest, Compton-thick column densities
  * is clumpy; does not confuse the line-of-sight inclination and viewing angle parameters
  * fits objects in the local Universe well

* Corrects for galactic absorption
* Optional: add an apec :math:`L<10^{42}` erg/s contamination component (set WITHAPEC=1)
* Can fit multiple observations simultaneously

I strongly recommend using the `xagnfitter <https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xagnfitter.py>`_
script instead of hardness ratios.

xagnfitter.py script
--------------------

It is included verbatim below:

.. literalinclude:: ../examples/sherpa/xagnfitter.py

Example run
--------------------

You can try running this against the AGN spectrum provided in the examples/sherpa/chandra/ folder.

The output should look something like this::

    $ WITHAPEC=0 MODELDIR=$HOME/Downloads/specmodels/ python3 ../xagnfitter.py 
    read ARF file cdfs4Ms_179.arf
    [sherpa.astro.io INFO]: read ARF file cdfs4Ms_179.arf
    read RMF file cdfs4Ms_179.rmf
    [sherpa.astro.io INFO]: read RMF file cdfs4Ms_179.rmf
    read ARF (background) file cdfs4Ms_179.arf
    [sherpa.astro.io INFO]: read ARF (background) file cdfs4Ms_179.arf
    read RMF (background) file cdfs4Ms_179.rmf
    [sherpa.astro.io INFO]: read RMF (background) file cdfs4Ms_179.rmf
    read background file cdfs4Ms_179_bkg.pi
    [sherpa.astro.io INFO]: read background file cdfs4Ms_179_bkg.pi
     Solar Abundance Vector set to wilm:  Wilms, J., Allen, A. & McCray, R. ApJ 542 914 (2000) (abundances are set to zero for those elements not included in the paper).
     Cross Section Table set to vern:  Verner, Ferland, Korista, and Yakovlev 1996
    loading nH from 179.pi.nh (expecting something like 1e21 in there)
    setting galactic nH to 0.0088 [units of 1e22/cm²]
    combining components
    linking parameters
    setting redshift
    creating priors
    setting source and background model ...
    [bxa.Fitter INFO]: PCAFitter(for ID=1)
    [bxa.Fitter INFO]: loading PCA information from /home/user/bin/ciao-4.13/ots/lib/python3.7/site-packages/bxa/sherpa/background/chandra_1024.json
    [bxa.Fitter INFO]: fitting background of ID=1 using PCA method
    [bxa.Fitter INFO]: have 2521 background counts for deconvolution
    [bxa.Fitter INFO]: fit: initial PCA decomposition: [ 3.40159007e+00 -6.29801294e-04  1.15351627e-02  8.72323308e-03
      7.12865031e-03 -8.44044155e-03 -4.73346905e-03  5.83916607e-04
     -4.79908423e-04 -9.18585670e-04 -1.87723042e-04]
    [bxa.Fitter INFO]: fit: first full fit done
    [bxa.Fitter INFO]: fit: parameters: [-0.24972366790501033, 0.1237617574156954, 1.619517445974407, -0.3466035071997995, -0.7134659455645368, 1.364537126802302, 0.10510218498377895, -0.41467820359992413, 0.9144121142234537, -1.144703704590311, -2.4736204436412694]
    tbvabs Version 2.3
    Cosmic absorption with grains and H2, modified from
    Wilms, Allen, & McCray, 2000, ApJ 542, 914-924
    Questions: Joern Wilms
    joern.wilms@sternwarte.uni-erlangen.de
    joern.wilms@fau.de

    http://pulsar.sternwarte.uni-erlangen.de/wilms/research/tbabs/

    PLEASE NOTICE:
    To get the model described by the above paper
    you will also have to set the abundances:
       abund wilm

    Note that this routine ignores the current cross section setting
    as it always HAS to use the Verner cross sections as a baseline.
    [bxa.Fitter INFO]: fit: stat: 2512.626888813398
    [bxa.Fitter INFO]: fit: second full fit from zero
    [bxa.Fitter INFO]: fit: parameters: [-0.24972366790501033, 0.1237617574156954, 1.619517445974407, -0.3466035071997995, -0.7134659455645368, 1.364537126802302, 0.10510218498377895, -0.41467820359992413, 0.9144121142234537, -1.144703704590311, -2.4736204436412694]
    [bxa.Fitter INFO]: fit: stat: 2512.421580816705
    [bxa.Fitter INFO]: fit: using zero-fit
    11 parameters, stat=2512.42
    --> 10 parameters, stat=2521.66
    --> 9 parameters, stat=2535.61
    --> 8 parameters, stat=2541.22
    --> 7 parameters, stat=2544.23
    --> 6 parameters, stat=2547.09
    --> 5 parameters, stat=2656.04
    --> 4 parameters, stat=2723.56
    --> 3 parameters, stat=2796.09
    --> 2 parameters, stat=2798.28
    --> 1 parameters, stat=5724.59

    Background PCA fitting AIC results:
    -----------------------------------

    stat Ncomp AIC
    5724.6  1 5726.6
    2798.3  2 2802.3
    2796.1  3 2802.1
    2723.6  4 2731.6
    2656.0  5 2666.0
    2547.1  6 2559.1
    2544.2  7 2558.2
    2541.2  8 2557.2
    2535.6  9 2553.6
    2521.7 10 2541.7
    2512.4 11 2534.4

    Increasing parameters again...
    Final choice: 11 parameters, aic=2534.42

    Adding Gaussian#1
    largest remaining discrepancy at 7.643keV[489], need 5 counts
    placing gaussian at 7.64keV, with power 9.18976317187143e-05
    with Gaussian: 2540.1578140930155 ; change: 5.7 (negative is good)
    not significant, rejecting
    running BXA ...
    [ultranest] Sampling 400 live points from prior ...
    [ultranest INFO]: Sampling 400 live points from prior ...


    Mono-modal Volume: ~exp(-3.81) * Expected Volume: exp(0.00) Quality: ok

    src.level       :      -8.0|************************************** ******************************************|     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ************************************  +2.3                  |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -2.2|********** ******************************************* **************************|     +1.8

    Z=-4e+08(0.00%) | Like=-3.8e+08..-2.2e+03 [-2.967e+11..-2604] | it/evals=80/490 eff=88.8889% N=400 

    Mono-modal Volume: ~exp(-3.81)   Expected Volume: exp(-0.23) Quality: ok

    src.level       :      -8.0|**************************************************************** ******* **      |     +3.0
    torus.phoindex  :      +1.2|           +1.5  *** ****************************** *****  +2.3                  |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|******************************************************** ************************|     -1.0
    pca1.lognorm    :      -2.2|******************** ********************************* ******************** *****|     +1.8

    Z=-9978548.5(0.00%) | Like=-9923310.86..-2028.74 [-2.967e+11..-2604] | it/evals=160/596 eff=81.6327% N=400 0 

    Mono-modal Volume: ~exp(-3.81)   Expected Volume: exp(-0.45) Quality: ok

    src.level       :      -8.0|*********************************************************** ****  +0.6           |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ****************************** *****  +2.3                  |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0| ******************************************************* ************************|     -1.0
    pca1.lognorm    :      -2.2|**************************************************** ****************************|     +1.8

    Z=-594817.6(0.00%) | Like=-586915.40..-2028.74 [-2.967e+11..-2604] | it/evals=265/765 eff=72.6027% N=400 0 

    Mono-modal Volume: ~exp(-4.47) * Expected Volume: exp(-0.67) Quality: ok

    src.level       :      -8.0|************************************************* *****  -0.6                    |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ****************************** *****  +2.3                  |     +2.8
    src.nH          :     +20.0|********* ***********************************************************************|    +26.0
    src.softscatnorm:      -7.0|**************************************** ****************************************|     -1.0
    pca1.lognorm    :      -2.2|**************************************************** ************************    |     +1.8

    Z=-173190.1(0.00%) | Like=-160428.37..-2028.74 [-2.967e+11..-2604] | it/evals=353/918 eff=68.1467% N=400 

    Mono-modal Volume: ~exp(-4.72) * Expected Volume: exp(-0.90) Quality: ok

    src.level       :      -8.0|********************************************** ** *  -1.1                        |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** * **************************** *****  +2.3                  |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|**************************************** ** *********************** *************|     -1.0
    pca1.lognorm    :      -2.2|********************************************************************  +1.1       |     +1.8

    Z=-59053.5(0.00%) | Like=-58881.01..-2028.74 [-2.967e+11..-2604] | it/evals=440/1072 eff=65.4762% N=400  

    Mono-modal Volume: ~exp(-4.99) * Expected Volume: exp(-1.12) Quality: ok

    src.level       :      -8.0|**********************************************  -1.8                             |     +3.0
    torus.phoindex  :      +1.2|            +1.6  *  ****************************** *****  +2.3                  |     +2.8
    src.nH          :     +20.0|**** *********************************************** ****************************|    +26.0
    src.softscatnorm:      -7.0|******************************************************************* *************|     -1.0
    pca1.lognorm    :      -2.2|************************************************************  +0.7               |     +1.8

    Z=-35508.0(0.00%) | Like=-35436.48..-2028.74 [-2.967e+11..-2604] | it/evals=520/1218 eff=63.5697% N=400 

    Mono-modal Volume: ~exp(-4.99)   Expected Volume: exp(-1.35) Quality: ok

    src.level       :      -8.0|*********************************************  -1.9                              |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ****************************** ** **  +2.3                  |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -2.2|*********************************************************  +0.6                  |     +1.8

    Z=-25985.5(0.00%) | Like=-25963.87..-2028.74 [-2.967e+11..-2604] | it/evals=626/1443 eff=60.0192% N=400 

    Mono-modal Volume: ~exp(-4.99)   Expected Volume: exp(-1.57) Quality: ok

    src.level       :      -8.0|******************************************** *  -1.8                             |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ***********************************  +2.3                   |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|***************************************************************** ***************|     -1.0
    pca1.lognorm    :      -2.2|        ** *********************************************  +0.5                   |     +1.8

    Z=-19777.5(0.00%) | Like=-19690.19..-2028.74 [-2.967e+11..-2604] | it/evals=714/1689 eff=55.3918% N=400 

    Mono-modal Volume: ~exp(-4.99)   Expected Volume: exp(-1.80) Quality: ok

    src.level       :      -8.0|********************************************  -2.1                               |     +3.0
    torus.phoindex  :      +1.2|            +1.6  *  ******************************* ***  +2.3                   |     +2.8
    src.nH          :     +20.0|*************************************************************************** *****|    +26.0
    src.softscatnorm:      -7.0|********************************************************* ***********************|     -1.0
    pca1.lognorm    :      -2.2|         -1.4  ***************************************  +0.4                     |     +1.8

    Z=-15401.1(0.00%) | Like=-15359.21..-1674.24 [-2.967e+11..-2604] | it/evals=800/1903 eff=53.2269% N=400 

    Mono-modal Volume: ~exp(-5.44) * Expected Volume: exp(-2.02) Quality: ok

    src.level       :      -8.0|***************************************** **  -2.1                               |     +3.0
    torus.phoindex  :      +1.2|            +1.6  *  ******************************* ***  +2.3                   |     +2.8
    src.nH          :     +20.0|*************************************************************************** *****|    +26.0
    src.softscatnorm:      -7.0|********************************************************* ******* ***************|     -1.0
    pca1.lognorm    :      -2.2|            -1.3  * *********************************  +0.4                      |     +1.8

    Z=-12325.4(0.00%) | Like=-12258.25..-1674.24 [-2.967e+11..-2604] | it/evals=880/2085 eff=52.2255% N=400 

    Mono-modal Volume: ~exp(-5.59) * Expected Volume: exp(-2.25) Quality: ok

    src.level       :      -8.0|***************************************** *  -2.2                                |     +3.0
    torus.phoindex  :      +1.2|            +1.6  *  ******************************* ***  +2.3                   |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -2.2|                -1.1  *****************************  +0.3                        |     +1.8

    Z=-9408.8(0.00%) | Like=-9399.06..-1674.24 [-2.967e+11..-2604] | it/evals=988/2343 eff=50.8492% N=400 0 

    Mono-modal Volume: ~exp(-6.22) * Expected Volume: exp(-2.47) Quality: ok

    src.level       :      -8.0|*****************************************  -2.4                                  |     +3.0
    torus.phoindex  :      +1.2|               +1.6  ***********************************  +2.3                   |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|   ** ********************************************  +0.2                         |     +1.0

    Z=-7834.1(0.00%) | Like=-7815.01..-1674.24 [-2.967e+11..-2604] | it/evals=1069/2510 eff=50.6635% N=400 

    Mono-modal Volume: ~exp(-6.25) * Expected Volume: exp(-2.70) Quality: ok

    src.level       :      -8.0|*****************************************  -2.4                                  |     +3.0
    torus.phoindex  :      +1.2|            +1.6  *  ***********************************  *  +2.4                |     +2.8
    src.nH          :     +20.0|************************** ******************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|        * ***************************************  +0.2                          |     +1.0

    Z=-6716.3(0.00%) | Like=-6695.76..-1674.24 [-2.967e+11..-2604] | it/evals=1160/2702 eff=50.3910% N=400 

    Mono-modal Volume: ~exp(-6.25)   Expected Volume: exp(-2.92) Quality: ok

    src.level       :      -8.0|*****************************************  -2.4                                  |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ** ************************************ **  +2.4               |     +2.8
    src.nH          :     +20.0|* ************************ ******************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|          * ************************************  +0.2                           |     +1.0

    Z=-6084.1(0.00%) | Like=-6075.22..-1666.27 [-2.967e+11..-2604] | it/evals=1240/2898 eff=49.6397% N=400 

    Mono-modal Volume: ~exp(-6.62) * Expected Volume: exp(-3.15) Quality: ok

    src.level       :      -8.0|*****************************************  -2.4                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  **************************************** **  +2.4               |     +2.8
    src.nH          :     +20.0|* *******************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|       -0.7  **********************************  +0.1                            |     +1.0

    Z=-5625.1(0.00%) | Like=-5603.01..-1666.27 [-2.967e+11..-2604] | it/evals=1338/3177 eff=48.1815% N=400 

    Mono-modal Volume: ~exp(-6.62)   Expected Volume: exp(-3.37) Quality: ok

    src.level       :      -8.0|*****************************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  *** ************************************ **  +2.4               |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|******************** ******** ********************************** ****************|     -1.0
    pca1.lognorm    :      -1.0|      -0.7  * *********************************  +0.1                            |     +1.0

    Z=-5348.2(0.00%) | Like=-5337.89..-1666.27 [-2.967e+11..-2604] | it/evals=1435/3538 eff=45.7298% N=400 

    Mono-modal Volume: ~exp(-6.62)   Expected Volume: exp(-3.60) Quality: ok

    src.level       :      -8.0|  **** **********************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  *** ************************************ **  +2.4               |     +2.8
    src.nH          :     +20.0|****************************************** **************************************|    +26.0
    src.softscatnorm:      -7.0|**************************************************************** ****************|     -1.0
    pca1.lognorm    :      -1.0|        -0.6  ********************************  +0.1                             |     +1.0

    Z=-5057.7(0.00%) | Like=-5041.99..-1664.55 [-2.967e+11..-2604] | it/evals=1520/3827 eff=44.3537% N=400 

    Mono-modal Volume: ~exp(-6.80) * Expected Volume: exp(-3.82) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.79
    src.level       :      -8.0|       **********************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  *** ************************************ **  +2.4               |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|         -0.6  *******************************  +0.1                             |     +1.0

    Z=-4648.4(0.00%) | Like=-4631.08..-1664.55 [-2.967e+11..-2604] | it/evals=1600/4128 eff=42.9185% N=400 

    Mono-modal Volume: ~exp(-6.80)   Expected Volume: exp(-4.05) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|         *******************************  -2.6                                   |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************ ** **  +2.4                |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|        -0.61  * ****************************  +0.10                             |    +1.00

    Z=-4286.3(0.00%) | Like=-4268.66..-1664.55 [-2.967e+11..-2604] | it/evals=1700/4431 eff=42.1732% N=400 

    Mono-modal Volume: ~exp(-6.95) * Expected Volume: exp(-4.27) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|     -6.4  ******************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************ ** **  +2.4                |     +2.8
    src.nH          :     +20.0|*********************************************************************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|         -0.58  * **************************  +0.08                              |    +1.00

    Z=-3849.9(0.00%) | Like=-3839.24..-1664.55 [-2.967e+11..-2604] | it/evals=1796/4792 eff=40.8925% N=400 

    Mono-modal Volume: ~exp(-6.95)   Expected Volume: exp(-4.50) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|       -6.2  ****************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************ ** **  +2.4                |     +2.8
    src.nH          :     +20.0|*********************************************************** *********************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|           -0.53  **************************  +0.07                              |    +1.00

    Z=-3599.9(0.00%) | Like=-3588.20..-1664.55 [-2.967e+11..-2604] | it/evals=1880/5224 eff=38.9718% N=400 

    Mono-modal Volume: ~exp(-6.95)   Expected Volume: exp(-4.73) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.83
    src.level       :      -8.0|       -6.2  ****************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************ ** **  +2.4                |     +2.8
    src.nH          :     +20.0|**************************************************** ****************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|            -0.52  ************************  +0.06                               |    +1.00

    Z=-3331.9(0.00%) | Like=-3317.44..-1664.55 [-2.967e+11..-2604] | it/evals=1971/5680 eff=37.3295% N=400 

    Mono-modal Volume: ~exp(-6.95)   Expected Volume: exp(-4.95) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|      -6.2  ****************************  -2.7                                   |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************ **  +2.3                   |     +2.8
    src.nH          :     +20.0|*********************************************  ****  ****************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|              -0.48  *********************  +0.03                                |    +1.00

    Z=-3126.5(0.00%) | Like=-3115.87..-1587.90 [-2.967e+11..-2604] | it/evals=2063/6274 eff=35.1209% N=400 

    Mono-modal Volume: ~exp(-8.08) * Expected Volume: exp(-5.18) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|        -6.1  **************************  -2.7                                   |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ***************************************  +2.3                   |     +2.8
    src.nH          :     +20.0|******************************************* * *** *  ** ****************** ******|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|              -0.47  *********************  +0.02                                |    +1.00

    Z=-2950.6(0.00%) | Like=-2939.53..-1587.90 [-2.967e+11..-2604] | it/evals=2147/6645 eff=34.3795% N=400 

    Mono-modal Volume: ~exp(-8.08)   Expected Volume: exp(-5.40) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.83
    src.level       :      -8.0|         -5.9  *************************  -2.7                                   |     +3.0
    torus.phoindex  :      +1.2|            +1.6  **************************************    *  +2.4              |     +2.8
    src.nH          :     +20.0|******************************************  * *** *   * *** *********************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :    -1.000|             -0.461  ********************  +0.004                                |   +1.000

    Z=-2743.1(0.00%) | Like=-2730.97..-1587.90 [-2.967e+11..-2604] | it/evals=2240/7157 eff=33.1508% N=400 

    Mono-modal Volume: ~exp(-8.08)   Expected Volume: exp(-5.63) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|         -5.8  ************* ************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ***************************************   *  +2.4              |     +2.8
    src.nH          :     +20.0|******************************************     *      * *** *********************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                -0.43  *****************  -0.01                                  |    +1.00

    Z=-2613.0(0.00%) | Like=-2600.30..-1587.90 [-2603.9145..-1664.5521] | it/evals=2331/7782 eff=31.5768% N=400 

    Mono-modal Volume: ~exp(-8.89) * Expected Volume: exp(-5.85) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.80
    src.level       :      -8.0|          -5.8  *************************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|            +1.6  ***************************************   *  +2.4              |     +2.8
    src.nH          :     +20.0|**************************************** *            ***** *********** *********|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                -0.43  *****************  -0.01                                  |    +1.00

    Z=-2470.2(0.00%) | Like=-2457.99..-1587.90 [-2603.9145..-1664.5521] | it/evals=2427/8467 eff=30.0855% N=400 

    Mono-modal Volume: ~exp(-8.89)   Expected Volume: exp(-6.08) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.80
    src.level       :      -8.0|          -5.8  ************ ************  -2.5                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ** ***************************************   *  +2.4              |     +2.8
    src.nH          :     +20.0|******************************************            * *************** *** *****|    +26.0
    src.softscatnorm:      -7.0|******************* ** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                 -0.39  ****************  -0.03                                  |    +1.00

    Z=-2374.8(0.00%) | Like=-2362.86..-1587.90 [-2603.9145..-1664.5521] | it/evals=2512/9047 eff=29.0505% N=400 

    Mono-modal Volume: ~exp(-8.89)   Expected Volume: exp(-6.30) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|          -5.8  ************ **********  -2.7                                    |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ** **************************************    *  +2.4              |     +2.8
    src.nH          :     +20.0|**************************************** *              *************** *********|    +26.0
    src.softscatnorm:      -7.0|************* *************************************************** ***************|     -1.0
    pca1.lognorm    :     -1.00|                 -0.39  ***************  -0.04                                   |    +1.00

    Z=-2292.2(0.00%) | Like=-2280.42..-1587.90 [-2603.9145..-1664.5521] | it/evals=2609/9756 eff=27.8858% N=400 

    Mono-modal Volume: ~exp(-8.89)   Expected Volume: exp(-6.53) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|          -5.8  ************ *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *** ********************************  ****  +2.3                   |     +2.8
    src.nH          :     +20.0|****************************************                *************************|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                  -0.38  **************  -0.04                                   |    +1.00

    Z=-2213.4(0.00%) | Like=-2200.91..-1587.90 [-2603.9145..-1664.5521] | it/evals=2699/10432 eff=26.9039% N=400 

    Mono-modal Volume: ~exp(-8.89)   Expected Volume: exp(-6.75) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.85
    src.level       :      -8.0|          -5.8  ************ *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *** ******************************** *** *  +2.3                   |     +2.8
    src.nH          :     +20.0|****************************************                * ********* *** *********|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                  -0.37  **************  -0.05                                   |    +1.00

    Z=-2154.3(0.00%) | Like=-2142.24..-1587.90 [-2603.9145..-1664.5521] | it/evals=2777/11170 eff=25.7846% N=400 

    Mono-modal Volume: ~exp(-8.89)   Expected Volume: exp(-6.98) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.87
    src.level       :      -8.0|          -5.8  ************ *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ****** ******************************** **  +2.3                  |     +2.8
    src.nH          :     +20.0|****************************************                ************ ** *********|    +26.0
    src.softscatnorm:      -7.0|********************** ****************************************** ***************|     -1.0
    pca1.lognorm    :     -1.00|                   -0.36  *************  -0.06                                   |    +1.00

    Z=-2087.4(0.00%) | Like=-2074.70..-1587.90 [-2603.9145..-1664.5521] | it/evals=2876/12296 eff=24.1762% N=400 

    Mono-modal Volume: ~exp(-9.84) * Expected Volume: exp(-7.20) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.87
    src.level       :      -8.0|           -5.7  **********  *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ************************************* * **  +2.3                  |     +2.8
    src.nH          :     +20.0| ********* ****************************                 ******* **** ************|    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                   -0.35  ************  -0.06                                    |    +1.00

    Z=-2028.9(0.00%) | Like=-2016.01..-1587.90 [-2603.9145..-1664.5521] | it/evals=2960/13307 eff=22.9333% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-7.43) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.86
    src.level       :      -8.0|           -5.6  **********  *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ****** ****************************** * **  +2.3                  |     +2.8
    src.nH          :     +20.0| ***** *** ****************************               *** ********** ************|    +26.0
    src.softscatnorm:      -7.0|********************** ********************************** ***********************|     -1.0
    pca1.lognorm    :     -1.00|                   -0.34  ************  -0.08                                    |    +1.00

    Z=-1985.8(0.00%) | Like=-1973.22..-1575.73 [-2603.9145..-1664.5521] | it/evals=3054/14321 eff=21.9381% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-7.65) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|           -5.6  *********   *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *************************************** ***  +2.3                 |     +2.8
    src.nH          :     +20.0|  ** * * ** **** *********************                 ** ************ **********|    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                   -0.34  ************  -0.08                                    |    +1.00

    Z=-1950.5(0.00%) | Like=-1938.08..-1575.73 [-2603.9145..-1664.5521] | it/evals=3139/15245 eff=21.1452% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-7.88) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|           -5.6  *********   *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *************************************** ***  +2.3                 |     +2.8
    src.nH          :     +20.0| *** *        ************************                 ** ******** *** ******  **|    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                    -0.32  ***********  -0.08                                    |    +1.00

    Z=-1909.4(0.00%) | Like=-1896.02..-1575.73 [-2603.9145..-1664.5521] | it/evals=3237/16462 eff=20.1532% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-8.10) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.80
    src.level       :      -8.0|            -5.5  ********   *********  -2.9                                     |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ************************************** ***  +2.3                 |     +2.8
    src.nH          :     +20.0|   *       ** ************************                  * ******** *** ******  **|    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                    -0.32  **********  -0.09                                     |    +1.00

    Z=-1873.0(0.00%) | Like=-1859.18..-1565.16 [-2603.9145..-1664.5521] | it/evals=3329/18036 eff=18.8762% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-8.33) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.79
    src.level       :      -8.0|            -5.5  *********  ********  -3.1                                      |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ************************************** ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|   *    *   *   **********************                  * ********  ***** ** * **|    +26.0
    src.softscatnorm:      -7.0|****************************************************************** * ************|     -1.0
    pca1.lognorm    :     -1.00|                     -0.31  *********  -0.10                                     |    +1.00

    Z=-1840.8(0.00%) | Like=-1826.57..-1565.16 [-2603.9145..-1664.5521] | it/evals=3419/20101 eff=17.3544% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-8.55) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.77
    src.level       :      -8.0|            -5.5  ********    ******  -3.1                                       |     +3.0
    torus.phoindex  :      +1.2|           +1.5  * ******************************** **  **       *  +2.5         |     +2.8
    src.nH          :     +20.0|     +20.9  * * *** ******************                    ** ****  ****** *  *  *|    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :      -1.0|                      -0.3  *********  -0.1                                      |     +1.0

    Z=-1815.0(0.00%) | Like=-1801.26..-1565.16 [-2603.9145..-1664.5521] | it/evals=3505/21751 eff=16.4161% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-8.78) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.78
    src.level       :      -8.0|            -5.5  *********    *****  -3.2                                       |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************* ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|     +20.9  * * * *  ******************                   ** *  *  *****        *|    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :      -1.0|                      -0.3  *********  -0.1                                      |     +1.0

    Z=-1785.7(0.00%) | Like=-1771.13..-1565.16 [-2603.9145..-1664.5521] | it/evals=3596/23869 eff=15.3223% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-9.00) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.77
    src.level       :      -8.0|             -5.4  *******     ***  -3.4                                         |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ************************************* ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|            +21.5  * *****************                    **       *   *         |    +26.0
    src.softscatnorm:      -7.0|****************************************** ****** *******************************|     -1.0
    pca1.lognorm    :      -1.0|                      -0.3  *********  -0.1                                      |     +1.0

    Z=-1760.5(0.00%) | Like=-1746.68..-1564.50 [-2603.9145..-1664.5521] | it/evals=3688/25859 eff=14.4860% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-9.23) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******    **  -3.7                                           |     +3.0
    torus.phoindex  :      +1.2|*                * *********************************** ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|              +21.6  *****************                    *       *  +25.0       |    +26.0
    src.softscatnorm:      -7.0|****************************************** **************************************|     -1.0
    pca1.lognorm    :      -1.0|                      -0.3  ********  -0.1                                       |     +1.0

    Z=-1737.9(0.00%) | Like=-1724.21..-1564.50 [-2603.9145..-1664.5521] | it/evals=3774/28135 eff=13.6074% N=400 

    Mono-modal Volume: ~exp(-9.84)   Expected Volume: exp(-9.45) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******    *  -3.9                                            |     +3.0
    torus.phoindex  :      +1.2|*                ************************************* ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|               +21.6  ****************                            *  +25.0       |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                      -0.3  ********  -0.1                                       |     +1.0

    Z=-1717.0(0.00%) | Like=-1703.15..-1564.50 [-2603.9145..-1664.5521] | it/evals=3869/32066 eff=12.2182% N=400 

    Mono-modal Volume: ~exp(-12.71) * Expected Volume: exp(-9.68) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|           +1.5  ********************************* *** ***       *  +2.5         |     +2.8
    src.nH          :     +20.0|               +21.6  * **************  +22.8                                    |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  *******  -0.1                                       |     +1.0

    Z=-1700.2(0.00%) | Like=-1686.12..-1564.50 [-2603.9145..-1664.5521] | it/evals=3952/32430 eff=12.3384% N=400 

    Mono-modal Volume: ~exp(-12.71)   Expected Volume: exp(-9.90) Quality: ok

       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|           +1.6  ********************************* *******       *  +2.5         |     +2.8
    src.nH          :     +20.0|               +21.6  * **************  +22.8                                    |    +26.0
    src.softscatnorm:      -7.0|**************** ****************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  *******  -0.1                                       |     +1.0

    Z=-1685.9(0.00%) | Like=-1671.48..-1563.17 [-2603.9145..-1664.5521] | it/evals=4040/32889 eff=12.4350% N=400 

    Mono-modal Volume: ~exp(-13.07) * Expected Volume: exp(-10.13) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.75
       positive degeneracy between src.nH and src.level: rho=0.77
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *      ***************************************** **    *  +2.5         |     +2.8
    src.nH          :     +20.0|                 +21.8  **************  +22.8                                    |    +26.0
    src.softscatnorm:      -7.0|**************** ****************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  ******  -0.1                                        |     +1.0

    Z=-1672.4(0.00%) | Like=-1657.80..-1563.17 [-1664.5027..-1569.9864] | it/evals=4135/33285 eff=12.5741% N=400 

    Mono-modal Volume: ~exp(-13.32) * Expected Volume: exp(-10.35) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.77
       positive degeneracy between src.nH and src.level: rho=0.77
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *      ***************************************** *     *  +2.5         |     +2.8
    src.nH          :     +20.0|                 +21.8  *************  +22.7                                     |    +26.0
    src.softscatnorm:      -7.0|* *******************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  ******  -0.1                                        |     +1.0

    Z=-1662.5(0.00%) | Like=-1648.08..-1563.17 [-1664.5027..-1569.9864] | it/evals=4214/33643 eff=12.6764% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-10.58) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.80
    src.level       :      -8.0|             -5.3  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *    * ************************************** ** *  +2.4               |     +2.8
    src.nH          :     +20.0|                 +21.8  *************  +22.7                                     |    +26.0
    src.softscatnorm:      -7.0|* ******** **********************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  ******  -0.1                                        |     +1.0

    Z=-1651.0(0.00%) | Like=-1636.39..-1563.17 [-1664.5027..-1569.9864] | it/evals=4315/34202 eff=12.7655% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-10.80) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.82
       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.3  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *    ************************************** * ** *  +2.4               |     +2.8
    src.nH          :     +20.0|                 +21.8  * **********  +22.7                                      |    +26.0
    src.softscatnorm:      -7.0|* *******************************************************************************|     -1.0
    pca1.lognorm    :      -1.0|                       -0.3  ******  -0.1                                        |     +1.0

    Z=-1642.1(0.00%) | Like=-1627.18..-1561.61 [-1664.5027..-1569.9864] | it/evals=4400/34699 eff=12.8284% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-11.02) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.83
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *   *************************************** *  * *  +2.4               |     +2.8
    src.nH          :     +20.0|                 +21.8  * **********  +22.7                                      |    +26.0
    src.softscatnorm:      -7.0|************************************ ********************************************|     -1.0
    pca1.lognorm    :      -1.0|                        -0.3  *****  -0.1                                        |     +1.0

    Z=-1634.0(0.00%) | Like=-1619.31..-1561.61 [-1664.5027..-1569.9864] | it/evals=4489/35443 eff=12.8100% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-11.25) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.82
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *   *********************************** *** *  * *  +2.4               |     +2.8
    src.nH          :     +20.0|                   +22.0  **********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|************************** ********* ********************************************|     -1.0
    pca1.lognorm    :      -1.0|                        -0.3  *****  -0.2                                        |     +1.0

    Z=-1627.5(0.00%) | Like=-1612.63..-1557.21 [-1664.5027..-1569.9864] | it/evals=4584/36151 eff=12.8220% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-11.47) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.85
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *   *************************************** *    *  +2.4               |     +2.8
    src.nH          :     +20.0|                   +22.0  **********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********** **********************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.25  *****  -0.15                                       |    +1.00

    Z=-1621.8(0.00%) | Like=-1607.06..-1557.21 [-1664.5027..-1569.9864] | it/evals=4676/36891 eff=12.8141% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-11.70) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.87
       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *    ************************************** *    *  +2.4               |     +2.8
    src.nH          :     +20.0|                   +22.0  **********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********** **********************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.25  *****  -0.15                                       |    +1.00

    Z=-1617.8(0.00%) | Like=-1602.92..-1557.21 [-1664.5027..-1569.9864] | it/evals=4760/37607 eff=12.7933% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-11.92) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.88
       positive degeneracy between src.nH and src.level: rho=0.76
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *    ************************************** *    *  +2.4               |     +2.8
    src.nH          :     +20.0|                   +22.0  **********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|***************************************************************** ***************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.24  ****  -0.16                                        |    +1.00

    Z=-1611.6(0.00%) | Like=-1596.33..-1557.21 [-1664.5027..-1569.9864] | it/evals=4859/38531 eff=12.7429% N=400 

    Mono-modal Volume: ~exp(-13.32)   Expected Volume: exp(-12.15) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.89
       positive degeneracy between src.nH and src.level: rho=0.78
    src.level       :      -8.0|             -5.4  *******  -4.5                                                 |     +3.0
    torus.phoindex  :      +1.2|          *    ************************************        *  +2.4               |     +2.8
    src.nH          :     +20.0|                   +22.0  **********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.24  ****  -0.16                                        |    +1.00

    Z=-1606.9(0.00%) | Like=-1591.79..-1557.21 [-1664.5027..-1569.9864] | it/evals=4947/39540 eff=12.6392% N=400 

    Mono-modal Volume: ~exp(-13.40) * Expected Volume: exp(-12.37) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.91
       positive degeneracy between src.nH and src.level: rho=0.78
    src.level       :      -8.0|             -5.3  ******  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|          *     *********************************  +2.2                          |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|*********** *********************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.24  ****  -0.17                                        |    +1.00

    Z=-1603.2(0.00%) | Like=-1588.00..-1557.21 [-1664.5027..-1569.9864] | it/evals=5039/40505 eff=12.5645% N=400 

    Mono-modal Volume: ~exp(-14.16) * Expected Volume: exp(-12.60) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.91
       positive degeneracy between src.nH and src.level: rho=0.81
    src.level       :      -8.0|             -5.3  ******  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|          *     ********************************  +2.1                           |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                       -0.24  ****  -0.17                                        |    +1.00

    Z=-1600.5(0.00%) | Like=-1585.37..-1557.21 [-1664.5027..-1569.9864] | it/evals=5126/41358 eff=12.5153% N=400 

    Mono-modal Volume: ~exp(-14.16)   Expected Volume: exp(-12.82) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.91
       positive degeneracy between src.nH and src.level: rho=0.81
    src.level       :      -8.0|             -5.3  ******  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  *********************************  +2.2                          |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1597.8(0.00%) | Like=-1582.26..-1557.21 [-1664.5027..-1569.9864] | it/evals=5217/42394 eff=12.4232% N=400 

    Mono-modal Volume: ~exp(-14.16)   Expected Volume: exp(-13.05) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.92
       positive degeneracy between src.nH and src.level: rho=0.83
    src.level       :      -8.0|             -5.3  ******  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ****************************** * *  +2.2                         |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********************** **********************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1595.3(0.00%) | Like=-1579.82..-1556.41 [-1664.5027..-1569.9864] | it/evals=5308/43204 eff=12.4007% N=400 

    Mono-modal Volume: ~exp(-14.16)   Expected Volume: exp(-13.27) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.93
       positive degeneracy between src.nH and src.level: rho=0.81
    src.level       :      -8.0|             -5.3  ******  -4.6                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ** *************************** *  *  *  +2.3                     |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|******************************************************* *************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1593.4(0.00%) | Like=-1577.83..-1556.41 [-1664.5027..-1569.9864] | it/evals=5397/44719 eff=12.1776% N=400 

    Mono-modal Volume: ~exp(-14.16)   Expected Volume: exp(-13.50) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.94
       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|             -5.3  ******  -4.6                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ** ****************************   ** *  +2.3                     |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|******************************************************* *************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1591.8(0.00%) | Like=-1576.21..-1556.41 [-1664.5027..-1569.9864] | it/evals=5489/46123 eff=12.0049% N=400 

    Mono-modal Volume: ~exp(-14.16)   Expected Volume: exp(-13.72) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.94
       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|              -5.3  *****  -4.6                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ** ****************************   ** *  +2.3                     |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1590.2(0.00%) | Like=-1574.28..-1556.41 [-1664.5027..-1569.9864] | it/evals=5579/47164 eff=11.9301% N=400 

    Mono-modal Volume: ~exp(-15.03) * Expected Volume: exp(-13.95) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.94
       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|              -5.3  *****  -4.6                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  *******************************   *  *  +2.3                     |     +2.8
    src.nH          :     +20.0|                    +22.0  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|************** ******************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1588.8(0.00%) | Like=-1572.79..-1556.41 [-1664.5027..-1569.9864] | it/evals=5664/48269 eff=11.8323% N=400 

    Mono-modal Volume: ~exp(-15.36) * Expected Volume: exp(-14.17) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.82
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|          +1.5  ******************************** *   *  +2.3                     |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.6                                       |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.17                                        |    +1.00

    Z=-1587.3(0.00%) | Like=-1571.04..-1556.41 [-1664.5027..-1569.9864] | it/evals=5756/49283 eff=11.7751% N=400 

    Mono-modal Volume: ~exp(-15.36)   Expected Volume: exp(-14.40) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *********************************     *  +2.3                     |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.6                                       |    +26.0
    src.softscatnorm:      -7.0|******************************************************* *************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.18                                        |    +1.00

    Z=-1585.9(0.00%) | Like=-1569.60..-1556.41 [-1569.9552..-1562.1643] | it/evals=5846/50379 eff=11.6969% N=400 

    Mono-modal Volume: ~exp(-15.92) * Expected Volume: exp(-14.62) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.85
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *****************************  **     *  +2.3                     |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.6                                       |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.23  ***  -0.18                                        |    +1.00

    Z=-1584.8(0.01%) | Like=-1568.21..-1555.91 [-1569.9552..-1562.1643] | it/evals=5938/51466 eff=11.6281% N=400 

    Mono-modal Volume: ~exp(-15.92)   Expected Volume: exp(-14.85) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.84
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *****************************   *     *  +2.3                     |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.6                                       |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1583.6(0.03%) | Like=-1566.95..-1555.91 [-1569.9552..-1562.1643] | it/evals=6029/53159 eff=11.4274% N=400 

    Mono-modal Volume: ~exp(-15.92)   Expected Volume: exp(-15.07) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.85
       positive degeneracy between src.nH and torus.phoindex: rho=0.75
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *****************************   *    **  +2.3                     |     +2.8
    src.nH          :     +20.0|                     +22.1  ********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1582.6(0.07%) | Like=-1566.02..-1555.91 [-1569.9552..-1562.1643] | it/evals=6117/54814 eff=11.2416% N=400 

    Mono-modal Volume: ~exp(-15.92)   Expected Volume: exp(-15.30) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.86
       positive degeneracy between src.nH and torus.phoindex: rho=0.77
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************************   *    *  +2.2                      |     +2.8
    src.nH          :     +20.0|                    +22.1  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|*********************************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1581.9(0.11%) | Like=-1565.27..-1555.35 [-1569.9552..-1562.1643] | it/evals=6206/56873 eff=10.9893% N=400 

    Mono-modal Volume: ~exp(-16.11) * Expected Volume: exp(-15.52) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.88
       positive degeneracy between src.nH and torus.phoindex: rho=0.79
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************************   *    *  +2.2                      |     +2.8
    src.nH          :     +20.0|                    +22.1  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|********************************************************************** **********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1581.2(0.24%) | Like=-1564.32..-1555.35 [-1569.9552..-1562.1643] | it/evals=6298/58887 eff=10.7682% N=400 

    Mono-modal Volume: ~exp(-16.11)   Expected Volume: exp(-15.75) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.95
       positive degeneracy between src.nH and src.level: rho=0.89
       positive degeneracy between src.nH and torus.phoindex: rho=0.80
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|         +1.5  *****************************   *    *  +2.2                      |     +2.8
    src.nH          :     +20.0|                    +22.1  *********  +22.6                                      |    +26.0
    src.softscatnorm:      -7.0|************** ******************************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1580.6(0.48%) | Like=-1563.67..-1555.35 [-1569.9552..-1562.1643] | it/evals=6389/61590 eff=10.4412% N=400 

    Mono-modal Volume: ~exp(-16.11)   Expected Volume: exp(-15.97) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.90
       positive degeneracy between src.nH and torus.phoindex: rho=0.83
    src.level       :      -8.0|              -5.3  *****  -4.7                                                  |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *******************************  *  +2.1                           |     +2.8
    src.nH          :     +20.0|                    +22.1  ********  +22.6                                       |    +26.0
    src.softscatnorm:      -7.0|***************************** ******************************** ******************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1580.0(0.60%) | Like=-1562.94..-1555.22 [-1569.9552..-1562.1643] | it/evals=6478/63946 eff=10.1942% N=400 

    Mono-modal Volume: ~exp(-16.72) * Expected Volume: exp(-16.20) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.89
       positive degeneracy between src.nH and torus.phoindex: rho=0.82
    src.level       :      -8.0|              -5.3  ****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *************************** ***  +2.1                              |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|******* ********************* ***************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1579.5(0.98%) | Like=-1562.37..-1555.07 [-1569.9552..-1562.1643] | it/evals=6568/65759 eff=10.0491% N=400 

    Mono-modal Volume: ~exp(-16.72)   Expected Volume: exp(-16.42) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.90
       positive degeneracy between src.nH and torus.phoindex: rho=0.84
    src.level       :      -8.0|              -5.3  ****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|         +1.5  ************************** ***  +2.1                              |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.5                                       |    +26.0
    src.softscatnorm:      -7.0|***************************** ***************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1579.1(1.50%) | Like=-1561.74..-1555.07 [-1562.1619..-1559.8406] | it/evals=6654/67999 eff=9.8433% N=400  

    Mono-modal Volume: ~exp(-16.88) * Expected Volume: exp(-16.65) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.91
       positive degeneracy between src.nH and torus.phoindex: rho=0.85
    src.level       :      -8.0|             -5.3  *****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|      +1.5  *  **************************  **  +2.1                              |     +2.8
    src.nH          :     +20.0|                    +22.1  ********  +22.5                                       |    +26.0
    src.softscatnorm:      -7.0|******** ******************** ***************************************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1578.7(1.71%) | Like=-1561.22..-1553.53 [-1562.1619..-1559.8406] | it/evals=6747/70297 eff=9.6528% N=400 

    Mono-modal Volume: ~exp(-17.05) * Expected Volume: exp(-16.87) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.91
       positive degeneracy between src.nH and torus.phoindex: rho=0.86
    src.level       :      -8.0|              -5.3  ****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|         +1.5  **************************  **  +2.1                              |     +2.8
    src.nH          :     +20.0|                     +22.1  *******  +22.5                                       |    +26.0
    src.softscatnorm:      -7.0|******** **************************************** ********** ********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1578.3(2.07%) | Like=-1560.64..-1553.53 [-1562.1619..-1559.8406] | it/evals=6836/72814 eff=9.4402% N=400 

    Mono-modal Volume: ~exp(-17.05)   Expected Volume: exp(-17.10) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.91
       positive degeneracy between src.nH and torus.phoindex: rho=0.85
    src.level       :      -8.0|             -5.3  *****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  *************************** *  +2.0                                 |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|************************************************************ ****** *************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1577.9(3.29%) | Like=-1560.14..-1553.53 [-1562.1619..-1559.8406] | it/evals=6929/76201 eff=9.1410% N=400 

    Mono-modal Volume: ~exp(-17.05)   Expected Volume: exp(-17.32) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.96
       positive degeneracy between src.nH and src.level: rho=0.92
       positive degeneracy between src.nH and torus.phoindex: rho=0.87
    src.level       :      -8.0|             -5.3  *****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  *************************** *  +2.0                                 |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|*********************************** ******  **************** ********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  ***  -0.18                                        |    +1.00

    Z=-1577.6(4.67%) | Like=-1559.76..-1553.53 [-1559.8397..-1558.4448] | it/evals=7019/79336 eff=8.8920% N=400 

    Mono-modal Volume: ~exp(-17.23) * Expected Volume: exp(-17.55) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.97
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  *****  -4.8                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  *************************** *  +2.0                                 |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|***** ** ******************** ************* **** ***** **************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1577.3(5.72%) | Like=-1559.38..-1553.36 [-1559.8397..-1558.4448] | it/evals=7109/81518 eff=8.7638% N=400 

    Mono-modal Volume: ~exp(-17.27) * Expected Volume: exp(-17.77) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.97
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  *****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ***************************  +2.0                                   |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|************************** ** * **********  **** ***** **************************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1577.0(7.21%) | Like=-1558.87..-1553.36 [-1559.8397..-1558.4448] | it/evals=7198/84998 eff=8.5085% N=400 

    Mono-modal Volume: ~exp(-17.84) * Expected Volume: exp(-18.00) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|              -5.3  ****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ***************************  +2.0                                   |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|*** * ******** * ******** *** ************* **** ***** ***** ********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1576.8(8.11%) | Like=-1558.52..-1552.95 [-1559.8397..-1558.4448] | it/evals=7287/87763 eff=8.3411% N=400 

    Mono-modal Volume: ~exp(-18.04) * Expected Volume: exp(-18.23) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
       positive degeneracy between src.softscatnorm and src.nH: rho=0.77
    src.level       :      -8.0|              -5.3  ****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ***************************  +2.0                                   |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|*** * ********** ******** *** ************* **** ***** *** **********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1576.6(10.54%) | Like=-1558.11..-1552.95 [-1558.4392..-1557.5018] | it/evals=7379/90194 eff=8.2177% N=400 

    Mono-modal Volume: ~exp(-18.57) * Expected Volume: exp(-18.45) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.92
       positive degeneracy between src.softscatnorm and src.nH: rho=0.77
    src.level       :      -8.0|             -5.3  *****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ***************************  +2.0                                   |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|****** ********* ** ***** *** ******* *** * **** ***** *** **********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1576.3(12.76%) | Like=-1557.78..-1552.95 [-1558.4392..-1557.5018] | it/evals=7468/93284 eff=8.0401% N=400 

    Mono-modal Volume: ~exp(-18.61) * Expected Volume: exp(-18.68) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
       positive degeneracy between src.softscatnorm and src.nH: rho=0.76
    src.level       :      -8.0|             -5.3  *****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **************************  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|**  ** * ******* ** ** ** *** **** ****** ****** ** ** *** **********************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1576.1(14.63%) | Like=-1557.39..-1552.83 [-1557.4972..-1556.9304] | it/evals=7559/96755 eff=7.8449% N=400 

    Mono-modal Volume: ~exp(-18.61)   Expected Volume: exp(-18.90) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.92
       positive degeneracy between src.softscatnorm and src.nH: rho=0.76
    src.level       :      -8.0|             -5.3  *****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **************************  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|**  *    *******  * *****  ** ********************* **  ** ***** ****************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1576.0(16.77%) | Like=-1557.09..-1551.92 [-1557.4972..-1556.9304] | it/evals=7647/100684 eff=7.6253% N=400 

    Mono-modal Volume: ~exp(-18.61)   Expected Volume: exp(-19.13) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
       positive degeneracy between src.softscatnorm and src.nH: rho=0.75
    src.level       :      -8.0|             -5.3  *****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **************************  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|*** * *   *****   *  ****  ** ********************* ****** ** **  * *************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1575.8(19.94%) | Like=-1556.76..-1551.92 [-1556.9302..-1556.5382] | it/evals=7738/105087 eff=7.3916% N=400 

    Mono-modal Volume: ~exp(-18.61)   Expected Volume: exp(-19.35) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
    src.level       :      -8.0|              -5.3  ****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **************************  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|*** * *** *****      ** *  ** ** ** *********** *** ****** ** **  *  ************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1575.6(22.45%) | Like=-1556.46..-1551.92 [-1556.5371..-1556.2426] | it/evals=7829/111277 eff=7.0610% N=400 

    Mono-modal Volume: ~exp(-18.61)   Expected Volume: exp(-19.58) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
    src.level       :      -8.0|              -5.3  ****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **************************  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0| **   *** *****     *** ** *   **** * ** *** **  ** * ******* **  ** ************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1575.5(25.17%) | Like=-1556.14..-1551.92 [-1556.2404..-1556.0691] | it/evals=7918/117786 eff=6.7453% N=400 

    Mono-modal Volume: ~exp(-18.61)   Expected Volume: exp(-19.80) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  ****  -4.9                                                   |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ************************ *  +2.0                                    |     +2.8
    src.nH          :     +20.0|                     +22.1  ******  +22.5                                        |    +26.0
    src.softscatnorm:      -7.0|**    **  *****     **  **     ***  * *  ***  *  ** *   *  *  **   **************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.22  **  -0.19                                         |    +1.00

    Z=-1575.3(28.17%) | Like=-1555.84..-1551.87 [-1555.8381..-1555.7736]*| it/evals=8009/129235 eff=6.2165% N=400 

    Mono-modal Volume: ~exp(-18.96) * Expected Volume: exp(-20.03) Quality: correlation length: 4 (+)

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  ***  -4.9                                                    |     +3.0
    torus.phoindex  :      +1.2|       +1.5  ***********************  +1.9                                       |     +2.8
    src.nH          :     +20.0|                     +22.1  *****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|**     *  * * *     *   **     ***  * *   *   *   * *   *  *  **   **************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1575.2(31.36%) | Like=-1555.54..-1551.87 [-1555.5437..-1555.5419]*| it/evals=8099/142743 eff=5.6898% N=400 

    Mono-modal Volume: ~exp(-19.01) * Expected Volume: exp(-20.25) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  ***  -4.9                                                    |     +3.0
    torus.phoindex  :      +1.2|       +1.5  **********************  +1.9                                        |     +2.8
    src.nH          :     +20.0|                     +22.1  *****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0| *      *  **       *   **     *      *   *   *   **       *  **   **************|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1575.1(35.17%) | Like=-1555.29..-1551.87 [-1555.2946..-1555.2931]*| it/evals=8189/157441 eff=5.2146% N=400 

    Mono-modal Volume: ~exp(-19.78) * Expected Volume: exp(-20.48) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.97
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *********************  +1.9                                        |     +2.8
    src.nH          :     +20.0|                     +22.1  *****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|         *  *       *    ** *  *      *   *     * **       *       **** *********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1575.0(39.34%) | Like=-1555.06..-1551.87 [-1555.0551..-1555.0517]*| it/evals=8279/170357 eff=4.8712% N=400 

    Mono-modal Volume: ~exp(-19.78)   Expected Volume: exp(-20.70) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.87
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  *********************  +1.9                                        |     +2.8
    src.nH          :     +20.0|                     +22.1  *****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|         *                *    *      *   *     *          *       **** *********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.9(43.32%) | Like=-1554.82..-1551.87 [-1554.8210..-1554.8210]*| it/evals=8369/183052 eff=4.5819% N=400 

    Mono-modal Volume: ~exp(-19.78)   Expected Volume: exp(-20.93) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************* *  +1.9                                        |     +2.8
    src.nH          :     +20.0|                     +22.1  *****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|         *                *          **                    *     * *  * *********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.8(47.39%) | Like=-1554.58..-1551.59 [-1554.5798..-1554.5769]*| it/evals=8459/192463 eff=4.4043% N=400 

    Mono-modal Volume: ~exp(-19.78)   Expected Volume: exp(-21.15) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.92
       positive degeneracy between src.nH and torus.phoindex: rho=0.87
    src.level       :      -8.0|             -5.3  ****  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|          *   *******************  +1.9                                          |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                           -2.1  *       ********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.7(51.01%) | Like=-1554.37..-1551.59 [-1554.3701..-1554.3676]*| it/evals=8549/202006 eff=4.2404% N=400 

    Mono-modal Volume: ~exp(-20.87) * Expected Volume: exp(-21.38) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|             -5.3  ****  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|          *   ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                 -1.7  * ********|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.6(54.68%) | Like=-1554.15..-1551.59 [-1554.1479..-1554.1422]*| it/evals=8638/205659 eff=4.2083% N=400 

    Mono-modal Volume: ~exp(-20.87)   Expected Volume: exp(-21.60) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  ****  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|          *   ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                 -1.7  *  *******|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.6(58.99%) | Like=-1553.94..-1551.59 [-1553.9423..-1553.9392]*| it/evals=8727/210178 eff=4.1601% N=400 

    Mono-modal Volume: ~exp(-22.50) * Expected Volume: exp(-21.83) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                     -1.4  ******|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.5(63.37%) | Like=-1553.75..-1551.59 [-1553.7459..-1553.7455]*| it/evals=8819/211888 eff=4.1700% N=400 

    Mono-modal Volume: ~exp(-22.50)   Expected Volume: exp(-22.05) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                     -1.4  ******|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.4(67.79%) | Like=-1553.55..-1551.59 [-1553.5544..-1553.5517]*| it/evals=8909/213905 eff=4.1727% N=400 

    Mono-modal Volume: ~exp(-22.63) * Expected Volume: exp(-22.28) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  ***  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                     -1.4  ******|     -1.0
    pca1.lognorm    :     -1.00|                        -0.21  **  -0.19                                         |    +1.00

    Z=-1574.4(71.84%) | Like=-1553.39..-1551.59 [-1553.3907..-1553.3904]*| it/evals=8999/216150 eff=4.1710% N=400 

    Mono-modal Volume: ~exp(-22.63)   Expected Volume: exp(-22.50) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|             -5.3  ****  -5.0                                                    |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ******************  +1.8                                           |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                      -1.3  *****|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.3(74.76%) | Like=-1553.25..-1551.59 [-1553.2478..-1553.2465]*| it/evals=9089/218373 eff=4.1698% N=400 

    Mono-modal Volume: ~exp(-22.63)   Expected Volume: exp(-22.73) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ****************  +1.8                                             |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                      -1.3  *****|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.3(77.99%) | Like=-1553.11..-1551.59 [-1553.1059..-1553.1031]*| it/evals=9179/220699 eff=4.1666% N=400 

    Mono-modal Volume: ~exp(-22.90) * Expected Volume: exp(-22.95) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ****************  +1.8                                             |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                      -1.3  *****|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.3(80.72%) | Like=-1552.98..-1551.53 [-1552.9756..-1552.9731]*| it/evals=9266/223137 eff=4.1601% N=400 

    Mono-modal Volume: ~exp(-22.90)   Expected Volume: exp(-23.18) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  *  ***************  +1.8                                             |     +2.8
    src.nH          :     +20.0|                      +22.2  ****  +22.4                                         |    +26.0
    src.softscatnorm:      -7.0|                                                                        -1.2  ***|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.2(83.31%) | Like=-1552.88..-1551.53 [-1552.8772..-1552.8767]*| it/evals=9354/224790 eff=4.1686% N=400 

    Mono-modal Volume: ~exp(-23.18) * Expected Volume: exp(-23.40) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ****************  +1.8                                             |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.4                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                        -1.2  ***|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.2(85.47%) | Like=-1552.76..-1551.41 [-1552.7619..-1552.7619]*| it/evals=9448/226892 eff=4.1714% N=400 

    Mono-modal Volume: ~exp(-23.31) * Expected Volume: exp(-23.63) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.87
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ****************  +1.8                                             |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.4                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                        -1.2  ***|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.19                                         |    +1.00

    Z=-1574.2(87.39%) | Like=-1552.65..-1551.41 [-1552.6486..-1552.6479]*| it/evals=9538/229325 eff=4.1664% N=400 

    Mono-modal Volume: ~exp(-23.32) * Expected Volume: exp(-23.85) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.4                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                        -1.2  ***|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.2(89.23%) | Like=-1552.57..-1551.41 [-1552.5707..-1552.5678]*| it/evals=9626/232205 eff=4.1526% N=400 

    Mono-modal Volume: ~exp(-24.74) * Expected Volume: exp(-24.08) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.4                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                        -1.2  ***|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(90.80%) | Like=-1552.45..-1551.41 [-1552.4544..-1552.4532]*| it/evals=9717/234184 eff=4.1564% N=400 

    Mono-modal Volume: ~exp(-24.94) * Expected Volume: exp(-24.30) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(92.08%) | Like=-1552.37..-1551.40 [-1552.3716..-1552.3709]*| it/evals=9804/235548 eff=4.1693% N=400 

    Mono-modal Volume: ~exp(-25.28) * Expected Volume: exp(-24.53) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|             -5.3  ***  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|      +1.5  * ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(93.30%) | Like=-1552.28..-1551.40 [-1552.2766..-1552.2766]*| it/evals=9895/236970 eff=4.1827% N=400 

    Mono-modal Volume: ~exp(-25.52) * Expected Volume: exp(-24.75) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  **  -5.0                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(94.40%) | Like=-1552.21..-1551.40 [-1552.2052..-1552.2048]*| it/evals=9987/238449 eff=4.1954% N=400 

    Mono-modal Volume: ~exp(-25.52)   Expected Volume: exp(-24.98) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(95.31%) | Like=-1552.15..-1551.39 [-1552.1456..-1552.1439]*| it/evals=10077/240066 eff=4.2046% N=400 

    Mono-modal Volume: ~exp(-25.52)   Expected Volume: exp(-25.20) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(96.12%) | Like=-1552.08..-1551.39 [-1552.0780..-1552.0775]*| it/evals=10169/241859 eff=4.2115% N=400 

    Mono-modal Volume: ~exp(-25.52)   Expected Volume: exp(-25.43) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  ***************  +1.8                                              |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :     -1.00|                         -0.21  *  -0.20                                         |    +1.00

    Z=-1574.1(96.75%) | Like=-1552.02..-1551.39 [-1552.0210..-1552.0194]*| it/evals=10258/244611 eff=4.2005% N=400 

    Mono-modal Volume: ~exp(-26.59) * Expected Volume: exp(-25.65) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.88
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.8                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :    -1.000|                        -0.207  *  -0.197                                        |   +1.000

    Z=-1574.1(97.29%) | Like=-1551.97..-1551.39 [-1551.9659..-1551.9658]*| it/evals=10345/246068 eff=4.2110% N=400 

    Mono-modal Volume: ~exp(-26.59)   Expected Volume: exp(-25.88) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.92
       positive degeneracy between src.nH and torus.phoindex: rho=0.87
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :    -1.000|                        -0.206  *  -0.197                                        |   +1.000

    Z=-1574.1(97.76%) | Like=-1551.91..-1551.35 [-1551.9109..-1551.9102]*| it/evals=10436/247800 eff=4.2183% N=400 

    Mono-modal Volume: ~exp(-26.59)   Expected Volume: exp(-26.10) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.91
       positive degeneracy between src.nH and torus.phoindex: rho=0.85
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                         -1.1  **|     -1.0
    pca1.lognorm    :    -1.000|                        -0.206  *  -0.197                                        |   +1.000

    Z=-1574.1(98.17%) | Like=-1551.86..-1551.35 [-1551.8639..-1551.8639]*| it/evals=10528/249800 eff=4.2213% N=400 

    Mono-modal Volume: ~exp(-26.63) * Expected Volume: exp(-26.33) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.93
       positive degeneracy between src.nH and torus.phoindex: rho=0.89
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                          -1.1  *|     -1.0
    pca1.lognorm    :    -1.000|                        -0.206  *  -0.197                                        |   +1.000

    Z=-1574.1(98.48%) | Like=-1551.83..-1551.35 [-1551.8282..-1551.8280]*| it/evals=10617/251253 eff=4.2324% N=400 

    Mono-modal Volume: ~exp(-27.74) * Expected Volume: exp(-26.55) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.98
       positive degeneracy between src.nH and src.level: rho=0.94
       positive degeneracy between src.nH and torus.phoindex: rho=0.90
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                          -1.1  *|     -1.0
    pca1.lognorm    :    -1.000|                        -0.206  *  -0.197                                        |   +1.000

    Z=-1574.1(98.76%) | Like=-1551.78..-1551.35 [-1551.7806..-1551.7802]*| it/evals=10708/252264 eff=4.2515% N=400 

    Mono-modal Volume: ~exp(-27.74)   Expected Volume: exp(-26.78) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.99
       positive degeneracy between src.nH and src.level: rho=0.95
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                          -1.1  *|     -1.0
    pca1.lognorm    :    -1.000|                        -0.205  *  -0.197                                        |   +1.000

    Z=-1574.1(98.98%) | Like=-1551.75..-1551.35 [-1551.7514..-1551.7510]*| it/evals=10798/253662 eff=4.2636% N=400 

    Mono-modal Volume: ~exp(-27.74)   Expected Volume: exp(-27.00) Quality: ok

       positive degeneracy between torus.phoindex and src.level: rho=0.99
       positive degeneracy between src.nH and src.level: rho=0.95
       positive degeneracy between src.nH and torus.phoindex: rho=0.91
    src.level       :      -8.0|              -5.3  **  -5.1                                                     |     +3.0
    torus.phoindex  :      +1.2|        +1.5  **************  +1.7                                               |     +2.8
    src.nH          :     +20.0|                      +22.2  ***  +22.3                                          |    +26.0
    src.softscatnorm:      -7.0|                                                                          -1.1  *|     -1.0
    pca1.lognorm    :    -1.000|                        -0.205  *  -0.197                                        |   +1.000

    [ultranest] Explored until L=-2e+03  551.35 [-1551.7492..-1551.7487]*| it/evals=10806/253825 eff=4.2640% N=400 
    [ultranest INFO]: Explored until L=-2e+03  
    [ultranest] Likelihood function evaluations: 253825
    [ultranest INFO]: Likelihood function evaluations: 253825
    [ultranest] Writing samples and results to disk ...
    [ultranest INFO]: Writing samples and results to disk ...
    [ultranest] Writing samples and results to disk ... done
    [ultranest INFO]: Writing samples and results to disk ... done
    [ultranest]   logZ = -1574 +- 0.1586
    [ultranest INFO]:   logZ = -1574 +- 0.1586
    [ultranest] Posterior uncertainty strategy is satisfied (KL: 0.46+-0.06 nat, need <0.50 nat)
    [ultranest INFO]: Posterior uncertainty strategy is satisfied (KL: 0.46+-0.06 nat, need <0.50 nat)
    [ultranest] Evidency uncertainty strategy is satisfied (dlogz=0.36, need <0.5)
    [ultranest INFO]: Evidency uncertainty strategy is satisfied (dlogz=0.36, need <0.5)
    [ultranest]   logZ error budget: single: 0.22 bs:0.16 tail:0.01 total:0.16 required:<0.50
    [ultranest INFO]:   logZ error budget: single: 0.22 bs:0.16 tail:0.01 total:0.16 required:<0.50
    [ultranest] done iterating.
    [ultranest INFO]: done iterating.

    logZ = -1574.044 +- 0.363
      single instance: logZ = -1574.044 +- 0.219
      bootstrapped   : logZ = -1574.033 +- 0.363
      tail           : logZ = +- 0.010
    insert order U test : converged: False correlation: 10667.0 iterations

        src.level           -5.106 +- 0.071
        torus.phoindex      1.684 +- 0.091
        src.nH              22.292 +- 0.060
        src.softscatnorm    -1.5 +- 1.2
        pca1.lognorm        -0.2012 +- 0.0048
    plotting spectrum ...
    9221/11207 (82.29%)
    calculating intrinsic fluxes and distribution of model spectra
    WARNING: Clearing convolved model
    '(apply_rmf(apply_arf((3746510.0 * (xstablemodel.torus + xstablemodel.scat)))) + (apply_rmf(apply_arf(pca1)) * 0.03082443021040756))'
    for dataset 1
    [sherpa.ui.utils WARNING]: Clearing convolved model
    '(apply_rmf(apply_arf((3746510.0 * (xstablemodel.torus + xstablemodel.scat)))) + (apply_rmf(apply_arf(pca1)) * 0.03082443021040756))'
    for dataset 1
    saving distribution plot data
    

This produces the following output files::

    $ find multiple_out_zspec_*
    multiple_out_zspec_
    multiple_out_zspec_/debug.log
    multiple_out_zspec_/plots
    multiple_out_zspec_/plots/corner.pdf
    multiple_out_zspec_/plots/trace.pdf
    multiple_out_zspec_/plots/run.pdf
    multiple_out_zspec_/info
    multiple_out_zspec_/info/post_summary.csv
    multiple_out_zspec_/info/results.json
    multiple_out_zspec_/results
    multiple_out_zspec_/results/points.hdf5
    multiple_out_zspec_/extra
    multiple_out_zspec_/chains
    multiple_out_zspec_/chains/run.txt
    multiple_out_zspec_/chains/weighted_post_untransformed.txt
    multiple_out_zspec_/chains/equal_weighted_post.txt
    multiple_out_zspec_/chains/weighted_post.txt
    multiple_out_zspec_bkg_1.txt.gz
    multiple_out_zspec_intrinsic_photonflux.dist.gz
    multiple_out_zspec_src_1.txt.gz

"multiple_out_zspec_" is the prefix used when a .z redshift file is provided and it contains a single number.

The most important files are:

* plots/corner.pdf: plot of the parameter constraints and uncertainties and their correlations
* info/results.json: summary of all parameters, their uncertainties and estimated lnZ
* info/post_summary.csv: summary of all parameters and their uncertainties as CSV
* chains/equal_weighted_post.txt: contains posterior samples: each row is a model parameter vector. You can iterate through these, set up the model in pyxspec, and then do something with it (compute fluxes and luminosities, for example).
* multiple_out_zspec_intrinsic_photonflux.dist.gz: contains the posterior distribution of redshifts, photon and energy fluxes (see step 9 in the script code)
Convenience features
=======================================

.. contents:: :local:

Fixing FITS keywords
----------------------------

Source, background, ARF and RMF files should reference each other
by FITS header keywords.

The `fixkeywords.py <https://github.com/JohannesBuchner/BXA/blob/master/fixkeywords.py>`_
script adjusts the keywords::

	fixkeywords.py src.pi bkg.pi rmf.rmf arf.arf

This also fixes ARF/RMF that start the energy bounds at zero (which is invalid)
instead of a small number.

Galactic absorption
----------------------------

The galactic column density in the direction of the source is often needed.
https://www.swift.ac.uk/analysis/nhtot/donhtot.php provides a look-up service.

The `gal.py <https://github.com/JohannesBuchner/BXA/blob/master/gal.py>`_
script fetches the value from there::

	gal.py src.pi

and stores it in src.pi.nh.

You can also give this script multiple spectral files and it avoids duplicate requests.


Accelerating slow models
----------------------------

Some models are very slow (such as convolutions).
It is worthwhile to cache them or produce interpolation grids.

For Sherpa, caching of arbitrary models is provided by CachedModel, which you can use as a wrapper:

.. autoclass:: bxa.sherpa.cachedmodel.CachedModel
	:noindex:


Automatic production of an interpolation model is possible with the RebinnedModel:

.. autoclass:: bxa.sherpa.rebinnedmodel.RebinnedModel
	:noindex:


Xspec chain files
----------------------------

BXA, when run from pyxspec, also provides chain fits file compatible with the
mcmc feature in xspec/pyxspec. xspec error propagation tools can thus be used
after a BXA fit. In xspec, one can load it with::

	XSPEC12> chain load path/to/mychain.fits


Parallelisation
---------------------------

BXA supports parallelisation with MPI (Message Passing Interface).
This allows scaling the inference from laptops all the way to computing clusters.

To use it, install mpi4py and run your python script with mpiexec::

	$ mpiexec -np 4 python3 myscript.py

No modifications of your scripts are needed. However, you may want to 
run plotting and other post-analysis only on rank 0.

Analysis of many data sets, or of many models are trivial to parallelise.
If your script accepts a command line argument, unix tools such as
"make -j 10" and "xargs --max-args 1 --max-procs 10" 
can help run your code in parallel.


Verbosity
--------------------------

If you want to make the fitting more quiet, set `verbose=False` when calling `run()`.

You can find more instructions 
`how to reduce the output of the UltraNest fitting engine here <https://johannesbuchner.github.io/UltraNest/issues.html#how-do-i-suppress-the-output>`_.

Code inside a XSilence container disables Xspec chatter::

	from bxa.xspec.solver import XSilence
	
	with XSilence():
		# do something here

Model discovery
---------------------

Is the model the right one? Is there more in the data? These questions can not
be answered in a statistical way, **but** what we can do is 

1. generate ideas on what models could fit better
2. test those models for significance with model selection

For the first point, **Quantile-Quantile plots** provide a unbinned, less noisy alternative to 
residual plots.

.. image:: absorbed-qq_model_deviations.*

.. image:: absorbed-convolved_posterior.*

QQ plot example (left), with the corresponding spectrum for comparison (right).

In these plots, for each energy the number of counts observed with lower energy
are plotted on one axis, while the predicted are on the other axis.
If model and data agree perfectly, this would be a straight line. 
Deviances are indications of possible mis-fits.

This example is almost a perfect fit!
You can see a offset growing at 6-7 keV, which remains at higher energies.
This indicates that the data has more counts than the model there.

As the growth is in a S-shape, it is probably a Gaussian (see its `cumulative density function <https://en.wikipedia.org/wiki/Normal_distribution>`_).

Refer to the appendix of the `accompaning paper <cite>`_ for more examples.

Sherpa with BXA
=======================================

Begin by loading bxa in a session with sherpa loaded::

   import bxa.sherpa as bxa

.. _sherpa-priors:

Defining priors
---------------------

Define your background model and source model as usual in sherpa.
Then define the priors over the free parameters, for example::

   # three parameters we want to vary
   param1 = xsapec.myapec.norm
   param2 = xspowerlaw.mypowerlaw.norm
   param3 = xspowerlaw.mypowerlaw.PhoIndex

   # list of parameters
   parameters = [param1, param2, param3]
   # list of prior transforms
   priors = [
      bxa.create_uniform_prior_for(param1),
      bxa.create_loguniform_prior_for(param2),
      bxa.create_gaussian_prior_for(param3, 1.95, 0.15),
      # and more priors
   ]
   
   # make a single function:
   priorfunction = bxa.create_prior_function(priors)

Make sure you set the parameter minimum and maximum values to appropriate (a priori reasonable) values.
The limits are used to define the uniform and loguniform priors.

You can freeze the parameters you do not want to investigate, but BXA only modifies the parameters specified.
As a hint, you can find all thawed parameters of a model with::

   parameters = [for p in get_model().pars if not p.frozen and p.link is None]

You can also define your own prior functions, which transform 
unit variables unto the values needed for each parameter.
See the `UltraNest documentation on priors <https://johannesbuchner.github.io/UltraNest/priors.html>`_ 
for more details about this concept.
The script `examples/sherpa/example_automatic_background_model.py <https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/example_automatic_background_model.py>`_ 
gives an example of such a custom prior function (`limited_19_24`).

API information:

.. autofunction:: bxa.sherpa.create_jeffreys_prior_for
.. autofunction:: bxa.sherpa.create_uniform_prior_for
.. autofunction:: bxa.sherpa.create_gaussian_prior_for
.. autofunction:: bxa.sherpa.create_prior_function

.. _sherpa-prior-predictive-checks:

Prior Predictive Checks
------------------------

To check that your priors and model is okay and working,
create a flipbook of prior samples.

1) Pick a random sample from the prior::

   for parameter, prior_function in zip(parameters, priors):
       parameter.val = prior_function(numpy.random.uniform())

2) make a plot (plot_model, plot_source, etc.)

Repeat this 20 times and look at the plots.

Do the shapes and number of counts expected
look like a reasonable representation of your prior expectation?


.. _sherpa-run:

Running the analysis
---------------------

You need to specify a prefix, called *outputfiles_basename* where the files are stored.

::

   # see the pymultinest documentation for all options
   solver = bxa.BXASolver(prior=priorfunction, parameters=parameters,
		outputfiles_basename = "myoutputs/")
   results = solver.run(resume=True)

API information:

.. autoclass:: bxa.sherpa.BXASolver
   :members: run

.. _sherpa-analyse:

Parameter posterior plots
--------------------------

Credible intervals of the model parameters, and histograms 
(1d and 2d) of the marginal parameter distributions are plotted
in 'myoutputs/plots/corner.pdf' for you.

.. figure:: absorbed-corner.*
	:scale: 50%

You can also plot them yourself using corner, triangle and getdist, by
passing `results['samples']` to them.

For more information on the corner library used here, 
see https://corner.readthedocs.io/en/latest/.

Error propagation
---------------------

`results['samples']` provides access to the posterior samples (similar to a Markov Chain).
Use these to propagate errors:

* For every row in the chain, compute the quantity of interest
* Then, make a histogram of the results, or compute mean and standard deviations.

This preserves the structure of the uncertainty (multiple modes, degeneracies, etc.)

BXA also allows you to compute the fluxes corresponding to the 
parameter estimation, giving the correct probability distribution on the flux.
With distance information (fixed value or distribution), you can later infer
the correct luminosity distribution.

::

     dist = solver.get_distribution_with_fluxes(lo=2, hi=10)
     numpy.savetxt(out + prefix + "dist.txt", dist)

API information:

.. automethod:: bxa.sherpa.BXASolver.get_distribution_with_fluxes

This does nothing more than::

    r = []
    for row in results['samples']:
	    # set the parameter values to the current sample
	    for p, v in zip(parameters, row):
		    p.val = v
	    r.append(list(row) + [calc_photon_flux(lo=elo, hi=ehi), 
		    calc_energy_flux(lo=elo, hi=ehi)])

Such loops can be useful for computing obscuration-corrected, rest-frame luminosities,
(modifying the nH parameter and the energy ranges before computing the fluxes).

.. _sherpa-models:

.. include:: model_comparison.rst

.. _sherpa-design:

.. include:: experiment_design.rst

.. _sherpa-qq:

.. include:: model_discovery.rst

The *qq* function in the *qq* module allows you to create such plots easily, by
exporting the cumulative functions into a file.

.. autofunction:: bxa.sherpa.qq.qq_export
   :noindex:

Refer to the :ref:`accompaning paper <cite>`, which gives an introduction and 
detailed discussion on the methodology.
Recommendations for X-ray spectral analysis
--------------------------------------------

For good and valid results, experienced users of XSpec or Sherpa already do these:

1. Using a continuous background model (parameteric, albeit not necessarily physical),
   instead of "subtracting" or using bin-wise backgrounds (XSpec default).
   Commonly, a extraction region near the source is used to estimate the background.
2. Use C-Stat/Cash (poisson likelihood) instead of Chi^2 (gaussian likelihood)
3. Use unbinned spectra (except perhaps for visualization, albeit you can use QQ-plots there without loss of resolution)

Beyond these already accepted practices, we recommend:

4. Estimating the values with uncertainties using Bayesian inference (this software, or MCMC methods)
   instead of Contour-search, Fisher matrix, stepping, or other approximations
   
   Instead of a local optimization, the benefit is that a global search can deal with multiples solutions.
   Error propagation is easy too when using the posterior samples (similar to a Markov chain),
   and it preserves the structure of the error (dependence between parameters, etc.)

5. Comparing models using the computed evidence (this software)
   instead of Likelihood ratio tests (which are invalid for non-nested models)
   
   Bayesian model selection based on the "evidence" Z resolves a number of limitations
   of current methods, and is easy to do with this software.

See the :ref:`accompaning paper <cite>` for a detailed discussion and comparison! 
Xspec with BXA
=======================================

This documentation shows how to use BXA with `PyXSpec`_.

Begin by importing bxa in a python session with heasoft loaded:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 4-5

Load your data, define your background model and source model as usual
with `PyXSpec`_. 

.. _PyXspec: https://heasarc.gsfc.nasa.gov/xanadu/xspec/python/html/

.. _xspec-priors:

Defining priors
---------------------

Define your background model and source model as usual in xspec.
Then define the priors over the free parameters, for example:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 13-27

The above is taken from 
`examples/xspec/example_simplest.py <https://github.com/JohannesBuchner/BXA/blob/master/examples/xspec/example_simplest.py>`_.

Make sure you set the parameter minimum and maximum values to appropriate (a priori reasonable) values.
The limits are used to define the uniform and loguniform priors.

You can freeze the parameters you do not want to investigate, but BXA only modifies the parameters specified.

You can also define your own prior functions, which transform 
unit variables unto the values needed for each parameter.
See the `UltraNest documentation on priors <https://johannesbuchner.github.io/UltraNest/priors.html>`_ 
for more details about this concept.
The script `examples/xspec/example_advanced_priors.py <https://github.com/JohannesBuchner/BXA/blob/master/examples/xspec/example_advanced_priors.py>`_ 
gives an example of such a custom prior function (`my_custom_prior`).

API information:

.. autofunction:: bxa.xspec.create_uniform_prior_for
.. autofunction:: bxa.xspec.create_loguniform_prior_for
.. autofunction:: bxa.xspec.create_gaussian_prior_for
.. autofunction:: bxa.xspec.create_custom_prior_for


Prior Predictive Checks
------------------------

To check that your priors and model is okay and working,
create a flipbook of prior sample predictions.

1) Pick a random sample from the prior::

   import numpy
   from bxa.xspec.solver import set_parameters

   prior_function = bxa.create_prior_function(transformations)
   values = prior_function(numpy.random.uniform(size=len(transformations)))
   set_parameters(transformations, values)
   print("set to parameters:", values)

2) make a plot

Repeat this 20 times and look at the plots.

Do the shapes and number of counts expected
look like a reasonable representation of your prior expectation?


.. _xspec-run:

Running the analysis
---------------------

This runs the fit and stores the result in the specified output folder:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 30-32

.. autoclass:: bxa.xspec.BXASolver
   :noindex:

The returned results contain posterior samples and the Bayesian evidence.
These are also reported on the screen for you.

.. _xspec-analyse:

Parameter posterior plots
--------------------------

Credible intervals of the model parameters, and histograms 
(1d and 2d) of the marginal parameter distributions are plotted
in 'myoutputs/plots/corner.pdf' for you.

.. figure:: absorbed-corner.*
	:scale: 50%

You can also plot them yourself using corner, triangle and getdist, by
passing `results['samples']` to them.

For more information on the corner library used here, 
see https://corner.readthedocs.io/en/latest/.

Model checking
-----------------------

For this functionality, you also need scipy installed.

The following code creates a plot of the convolved posterior model:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 37-60

.. figure:: absorbed-convolved_posterior.*
	
	Example of the convolved spectrum with data.
	For each posterior sample (solution), the parameters are taken and put
	through the model. All such lines are plotted. Where the region is darker,
	more lines ended up, and thus it is more likely.
	
	The data points are adaptively binned to contain at least 20 counts.
	The error bars are created by asking: which model count rate can produce
	this amount of counts. 
	In a Poisson process, the inverse incomplete gamma
	function provides this answer. The 10%-90% probability range is used.
	
	On the colors of the data points:
	
	For all intents and purposes, you can ignore the colors.
	
	The colors are intended to aid the discovery of discrepancies, by using
	a custom Goodness of Fit measure. In this procedure (gof module), 
	a tree of the bins is built, i.e. in the first layer, every 2 bins 
	are merged, in the second, every 4 bins are merged, etc.
	Then, the counts in the bins are compared against with 
	the poisson process of the model. The worst case, i.e. the least likely 
	probability over the whole tree is considered. That is, for each bin,
	the lowest probability of all its merges is kept. Finally, this is multiplied
	by the number of nodes in the tree (as more comparisons lead to more 
	random chances).
	
	Then, if the probability for the bin is below :math:`10^{-2}`, the point is marked orange,
	and if it reaches below :math:`10^{-6}`, it is marked red.
	
	It is ok to ignore the colors, this computation is not used otherwise.

The following code creates a plot of the unconvolved posterior:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 65-79

.. figure:: absorbed-unconvolved_posterior.*
	
	Example of the unconvolved spectrum.
	For each posterior sample (solution), the parameters are taken and put
	through the model. All such lines are plotted. Where the region is darker,
	more lines ended up, and thus it is more likely.

For plotting the model parameters found against the data, use these functions.

.. automethod:: bxa.xspec.BXASolver.posterior_predictions_unconvolved
.. automethod:: bxa.xspec.BXASolver.posterior_predictions_convolved
.. automethod:: bxa.xspec.sinning.binning
    :noindex:


Error propagation
---------------------

`results['samples']` provides access to the posterior samples (similar to a Markov Chain).
Use these to propagate errors:

* For every row in the chain, compute the quantity of interest
* Then, make a histogram of the results, or compute mean and standard deviations.

This preserves the structure of the uncertainty (multiple modes, degeneracies, etc.)

*Continuing in Xspec*: A chain file, compatible with Xspec chain commands is 
written for you into *<outputfiles_basename>chain.fits*. In Xspec, load it using `"chain load" <https://heasarc.gsfc.nasa.gov/xanadu/xspec/manual/XSchain.html>`_.
This should set parameters, and compute flux estimates.

.. _xspec-models:

.. include:: model_comparison.rst

.. _xspec-design:

.. include:: experiment_design.rst

.. _xspec-qq:

.. include:: model_discovery.rst

For Xspec, the *qq* function in the *qq* module allows you to create such plots easily:

.. literalinclude:: ../examples/xspec/example_simplest.py
   :lines: 83-90

.. autofunction:: bxa.xspec.qq.qq

Refer to the :ref:`accompaning paper <cite>`, which gives an introduction and 
detailed discussion on the methodology.

Experiment design
------------------

We want to to evaluate whether a planned experiment can detect features or constrain parameters, 
i.e. determine the discriminatory power of future configurations/surveys/missions.

For this, simulate a few spectra using the appropriate response.

* Case 1: Can the experiment constrain the parameters?
	* Analyse and check what fraction of the posterior samples lie inside/outside the region of interest.
* Case 2: Can the experiment distinguish between two models?
	* Model selection as above.
* Case 3: Which sources (redshift range, luminosity, etc) can be distinguished?
	* Compute a grid of spectra. Do model selection at each point in the grid.

PCA-based background models
=======================================

The advantages of using background spectral models is that more information 
can be extracted from low-count data, as correlations between bins and instrument behaviours are known.

This page describes the machine-learning approach to derive empirical background models
published in `Simmonds et al., 2018, A&A, 618A, 66 <https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..66S/abstract>`_.
For Chandra/ACIS, XMM/EPIC, Swift/XRT, NuSTAR/FPMA, RXTE, 
Large archives of background spectra were used to derive principal components (PCA)
that empirically describe the background and its variations.

BXA includes these PCA models, which can be fitted to a specific background spectrum.

These PCA models are trained in log10(counts + 1) space to avoid negative counts.
The PCA models operate on detector channels and thus should never pass through the response.

The PCA models are limited in how well they can describe additive components such as 
gaussian emission lines. For this reason, the fitters also try adding Gaussian lines at 
the location of strongest fit mismatch.

The fits keep increasing complexity (first, the number of PCA components and then, the gaussians) as long as the AIC (Akaike information criterion) improves.

In Sherpa
-------------------------

After setting your source model (with set_model), use::

	from bxa.sherpa.background.pca import auto_background
	convmodel = get_model(id)
	bkg_model = auto_background(id)
	set_full_model(id, get_response(id)(model) + bkg_model * get_bkg_scale(id))

.. autofunction:: bxa.sherpa.background.pca.auto_background

A full example for fitting obscured Active Galactic Nuclei is available:
https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xagnfitter.py

In Xspec
-------------------------

In xspec, there are two steps. 

First, fit the background spectrum outside xspec using the autobackgroundmodel/fitbkg.py script.

You can find it here: https://github.com/JohannesBuchner/BXA/tree/master/autobackgroundmodel

It will give you instructions how to load the PCA model in xspec.


Creating a new PCA model
---------------------------

See https://github.com/JohannesBuchner/BXA/tree/master/autobackgroundmodel#create-a-model


If the fit is bad
------------------

This should be judged from the q-q plot that the background fitting produces.

If it is off, it shows at what energies the problem is.

It may be that your spectrum is somehow different than typical spectra.

You can try another method (e.g., empirical background models, or fall back to Wstat statistics).
bxa package
===========

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   bxa.sherpa
   bxa.xspec

Module contents
---------------

.. automodule:: bxa
   :members:
   :undoc-members:
   :show-inheritance:
Empirical background models
=======================================

The advantages of using background spectral models is that more information 
can be extracted from low-count data, as correlations between bins and instrument behaviours are known.

The :py:mod:`bxa.sherpa.background` module includes hand-crafted models to 
empirically describe the background spectra
for Chandra/ACIS, XMM/EPIC, Swift/XRT. 
This requires that you extracted a background spectrum (they are not ab initio predictions).

A mixture of powerlaws, gaussian lines and mekals are fitted to the background.
The best-fit model can then be used to fit the source. Optionally,
background parameters (such as the overall normalisation) can be varied with the source fit.

You may also be interested in the `PCA models <pca-background-models>`_ for Chandra/ACIS, XMM/EPIC, Swift/XRT, NuSTAR/FPMA, RXTE observations.


Empirical models for XMM/EPIC
--------------------------------

An example for XMM is available at
https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xmm/fit.py

.. autofunction:: bxa.sherpa.background.xmm::get_pn_bkg_model
	:noindex:
.. autofunction:: bxa.sherpa.background.xmm::get_mos_bkg_model
	:noindex:

The XMM model was developed by Richard Sturm at MPE.
The citation is `Maggi P., et al., 2014, A&A, 561, AA76 <https://ui.adsabs.harvard.edu/abs/2012A&A...546A.109M/abstract>`_.

Empirical models for Swift/XRT and Chandra/ACIS
------------------------------------------------

An example for Swift/XRT is available at
https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/swift/fit.py

First define which background you want to use and where you want to store intermediate results:
::
	
	# where to store intermediary fit information
	# usually the name of the spectral file
	filename = 'mybackgroundspecfile'
	
	# create a fitter for the desired type of spectrum
	from bxa.sherpa.background.models import SwiftXRTBackground, ChandraBackground
	from bxa.sherpa.background.fitter import SingleFitter
	fitter = SingleFitter(id, filename, SwiftXRTBackground)
	# or 
	fitter = SingleFitter(id, filename, ChandraBackground)

Finally, run the background fit::

	fitter.fit(plot=True)

The Chandra model was developed by Johannes Buchner at MPE. 
The citation is `Buchner et al., 2014, A&A, 564A, 125 <https://ui.adsabs.harvard.edu/abs/2014A%26A...564A.125B/abstract>`_.

The Swift/XRT model was developed by Johannes Buchner at PUC. 
The citation is `Buchner et al., 2017, MNRAS, 464, 4545 <https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.4545B/abstract>`_.

.. autoclass:: bxa.sherpa.background.models.SwiftXRTBackground
	:noindex:
.. autoclass:: bxa.sherpa.background.models.SwiftXRTWTBackground
	:noindex:
.. autoclass:: bxa.sherpa.background.models.ChandraBackground
	:noindex:

If the fit is bad
------------------

Try another method (e.g., PCA models, or fall back to Wstat statistics).
API
===

.. toctree::
   :maxdepth: 4

   bxa

Model comparison
----------------------

*examples/xspec/model_compare.py* shows an example of model selection. Keep in mind what model prior you would like to use.

* Case 1: Multiple models, want to find one best one to use from there on:
	* follow *examples/model_compare.py*, and pick the model with the highest evidence
* Case 2: Simpler and more complex models, want to find out which complexity is justified:
	* follow *examples/model_compare.py*, and keep the models above a certain threshold
* Case 3: Multiple models which could be correct, only interested in a parameter
	* Marginalize over the models: Use the posterior samples from each model, and weigh them by the 
	  relative probability of the models (weight = exp(lnZ))

Example output::

	jbuchner@ds42 $ python model_compare.py absorbed/ line/ simplest/

	Model comparison
	****************

	model simplest : log10(Z) = -1632.7  XXX ruled out
	model absorbed : log10(Z) =    -7.5  XXX ruled out
	model line     : log10(Z) =     0.0    <-- GOOD

	The last, most likely model was used as normalization.
	Uniform model priors are assumed, with a cut of log10(30) to rule out models.
	
	jbuchner@ds42 $ 

Here, the probability of the second-best model, "absorbed", is :math:`10^7.5` times
less likely than the model "line". As this exceeds our threshold (by a lot!)
we can claim the detection of an iron line!

Monte Carlo simulated spectra are recommended to derive a 
Bayes factor threshold for a preferred false selection rate.
You can find an example in the `Appendix of Buchner+14 <https://ui.adsabs.harvard.edu/abs/2014A%26A...564A.125B/abstract>`_
and in `the ultranest tutorial <https://johannesbuchner.github.io/UltraNest/example-sine-modelcomparison.html>`_.
bxa.sherpa.background package
=============================

Submodules
----------

bxa.sherpa.background.fitters module
------------------------------------

.. automodule:: bxa.sherpa.background.fitters
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.background.models module
-----------------------------------

.. automodule:: bxa.sherpa.background.models
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.background.pca module
--------------------------------

.. automodule:: bxa.sherpa.background.pca
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.background.xmm module
--------------------------------

.. automodule:: bxa.sherpa.background.xmm
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: bxa.sherpa.background
   :members:
   :undoc-members:
   :show-inheritance:
Welcome to BXA's documentation!
=======================================

.. include:: ../README.rst

Usage
--------

* We start with some notes on :doc:`best practice <best-practise>` in currently common 
  Maximum Likelihood analysis (fitting) and Bayesian analysis.

The usage is similar in Sherpa and Xspec conceptionally:

* Define **priors** as transformations from the unit cube to the parameter space
	* :ref:`In Sherpa <sherpa-priors>`
	* :ref:`In Xspec <xspec-priors>`
* **Run** the analysis
	* :ref:`In Sherpa <sherpa-run>`
	* :ref:`In Xspec <xspec-run>`
* **Analyse** the results: Get marginal distributions and computed evidence for the model
	* :ref:`In Sherpa <sherpa-analyse>`
	* :ref:`In Xspec <xspec-analyse>`
* **Model comparison**
	* :ref:`In Sherpa <sherpa-models>`
	* :ref:`In Xspec <xspec-models>`
* **Experiment design**
	* :ref:`In Sherpa <sherpa-design>`
	* :ref:`In Xspec <xspec-design>`
* **Model discovery**: Quantile-Quantile (QQ) plots
	* :ref:`In Sherpa <sherpa-qq>`
	* :ref:`In Xspec <xspec-qq>`

.. toctree::
   :maxdepth: 2
   :caption: Contents

   sherpa-analysis
   xspec-analysis
   convenience
   background-models
   pca-background-models
   xagnfitter
   contributing
   history
   modules

Indices and tables
==================
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
bxa.sherpa package
==================

Subpackages
-----------

.. toctree::
   :maxdepth: 4

   bxa.sherpa.background

Submodules
----------

bxa.sherpa.cachedmodel module
-----------------------------

.. automodule:: bxa.sherpa.cachedmodel
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.galabs module
------------------------

.. automodule:: bxa.sherpa.galabs
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.invgauss module
--------------------------

.. automodule:: bxa.sherpa.invgauss
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.priors module
------------------------

.. automodule:: bxa.sherpa.priors
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.qq module
--------------------

.. automodule:: bxa.sherpa.qq
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.rebinnedmodel module
-------------------------------

.. automodule:: bxa.sherpa.rebinnedmodel
   :members:
   :undoc-members:
   :show-inheritance:

bxa.sherpa.solver module
------------------------

.. automodule:: bxa.sherpa.solver
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: bxa.sherpa
   :members:
   :undoc-members:
   :show-inheritance:
bxa.xspec package
=================

Submodules
----------

bxa.xspec.gof module
--------------------

.. automodule:: bxa.xspec.gof
   :members:
   :undoc-members:
   :show-inheritance:

bxa.xspec.priors module
-----------------------

.. automodule:: bxa.xspec.priors
   :members:
   :undoc-members:
   :show-inheritance:

bxa.xspec.qq module
-------------------

.. automodule:: bxa.xspec.qq
   :members:
   :undoc-members:
   :show-inheritance:

bxa.xspec.sinning module
------------------------

.. automodule:: bxa.xspec.sinning
   :members:
   :undoc-members:
   :show-inheritance:

bxa.xspec.solver module
-----------------------

.. automodule:: bxa.xspec.solver
   :members:
   :undoc-members:
   :show-inheritance:

Module contents
---------------

.. automodule:: bxa.xspec
   :members:
   :undoc-members:
   :show-inheritance:
=======================================
Machine-learned PCA Background models
=======================================

The advantages of using background spectral models is that more information 
can be extracted from low-count data, as correlations between bins and instrument behaviours are known.

This page describes the machine-learning approach to derive empirical background models
published in `Simmonds et al., 2018, A&A, 618A, 66 <https://ui.adsabs.harvard.edu/abs/2018A%26A...618A..66S/abstract>`_.
For Chandra/ACIS, XMM/EPIC, Swift/XRT, NuSTAR/FPMA, RXTE, 
Large archives of background spectra were used to derive principal components (PCA)
that empirically describe the background and its variations.

BXA includes these PCA models, which can be fitted to a specific background spectrum.

These PCA models are trained in log10(counts + 1) space to avoid negative counts.
The PCA models operate on detector channels and thus should never pass through the response.

The PCA models are limited in how well they can describe additive components such as 
gaussian emission lines. For this reason, the fitters also try adding Gaussian lines at 
the location of strongest fit mismatch.

The fits keep increasing complexity (first, the number of PCA components and then, the gaussians) as long as the AIC (Akaike information criterion) improves.

In Sherpa
-------------------------

After setting your source model (with set_model), use::

	from bxa.sherpa.background.pca import auto_background
	convmodel = get_model(id)
	bkg_model = auto_background(id)
	set_full_model(id, get_response(id)(model) + bkg_model * get_bkg_scale(id))

.. autofunction:: bxa.sherpa.background.pca.auto_background

A full example for fitting obscured Active Galactic Nuclei is available:
https://github.com/JohannesBuchner/BXA/blob/master/examples/sherpa/xagnfitter.py

In Xspec
-------------------------

In xspec, there are two steps. 

1) precompute the best-fit continuum model

   fit the background spectrum outside xspec using the autobackgroundmodel/fitbkg.py script.

   You can find it here: https://github.com/JohannesBuchner/BXA/tree/master/autobackgroundmodel
   
   It needs the json files in the same folder.
   
   It will give you instructions how to load the PCA model in xspec.

2) Use the resulting atable with a dummy unit response matrix.

   fitbkg.py in step 1 gives you instructions how to load it in xspec.

Creating a new PCA model
---------------------------

If you have a new instrument or survey, create a model with:

`python autobackgroundmodel/__init__.py bkg1.pha bkg2.pha bkg3.pha`

The resulting file is telescope.json or telescope_instrument.json.


Fitting a model
----------------

Once you want to fit a specific sources, you can obtain its background model
in two different ways:

**For use in xspec**, you can precompute a background model:

`python autobackgroundmodel/fitbkg.py bkg1.pha [src1.pha]`

This will make a file bkg1.pha.bstat.out with the estimated per-channel count rate.

**For use in sherpa**, use the autobackground feature.
see `examples/sherpa/example_automatic_background.py`

If you have a new model (json file), it needs to be installed in the BXA folder,
probably at `~/.local/lib/python*/site-packages//bxa/sherpa/background/`.
