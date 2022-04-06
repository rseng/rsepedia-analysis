[![CI tests](https://github.com/BjoernBiltzinger/pyspi/actions/workflows/publish_pypi.yml/badge.svg)](https://github.com/BjoernBiltzinger/pyspi/actions/workflows/publish_pypi.yml)
[![Docs](https://github.com/BjoernBiltzinger/pyspi/actions/workflows/docs.yml/badge.svg)](https://pyspi.readthedocs.io/en/latest/)
[![codecov](https://codecov.io/gh/BjoernBiltzinger/pyspi/branch/master/graph/badge.svg)](https://codecov.io/gh/BjoernBiltzinger/pyspi)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6376003.svg)](https://doi.org/10.5281/zenodo.6376003)
[![JOSS](https://joss.theoj.org/papers/10.21105/joss.04017/status.svg)](https://doi.org/10.21105/joss.04017)
# pyspi
![alt text](https://raw.githubusercontent.com/BjoernBiltzinger/pyspi/master/docs/media/pypsi_logo2.png)

A python analysis framework for INTEGRAL/SPI

```PySPI``` provides a plugin for [3ML](https://threeml.readthedocs.io/en/stable/) for INTEGRAL/SPI data, which allows to analyze GRB data at the moment. In the future we plan to also add support for non transient sources.

## Installation

### Pip
```PySPI``` can be installed via pip.
```bash
pip install py-spi
```

### Github

To install the latest release from Github run
```bash
git clone https://github.com/BjoernBiltzinger/pyspi.git
```
After that first install the packages from the requirement.txt file with
```bash
cd pyspi
pip install -r requirements.txt
```
Now you can install ```PySPI``` with
```bash
python setup.py install
```

### Additional Data Files

There are a few large data files for the background model and the response that are not included in the Github repository. To get these data files run and specify the path where this data folder should be stored on your local machine. Here you have to change the /path/to/internal/data with the path you want to use on your local computer.
```bash
wget https://grb.mpe.mpg.de/pyspi_datafolder && unzip pyspi_datafolder
mv data /path/to/internal/data && rm -f pyspi_datafolder
```

### Environment Variables

Next you have to set two environment variable. One to define the path to the folder of the external data like the different SPI data files that will be downloaded by ```PySPI``` and one to define the path to the internal data folder we downloaded earlier.
```bash
export PYSPI=/path/to/external/datafolder
export PYSPI_PACKAGE_DATA=/path/to/internal/data
```

You should add these two line to your bashrc (or similar) file to automatically set this variable in every new terminal.

Now we are ready to go.

## Features

Please have a look at the [documentation](https://pyspi.readthedocs.io/en/latest/) to check out the features ```PySPI``` provides. There is also a [full example](https://pyspi.readthedocs.io/en/latest/notebooks/grb_analysis/), how to perform a spectral fit for the data for GRB120711A, as well as how to localize the GRB with ```PySPI```.

# Contributing 

Contributions to ```PySPI``` are always welcome. They can come in the form of:

## Issues

Please use the [Github issue tracking system for any bugs](https://github.com/BjoernBiltzinger/pyspi/issues), for questions, bug reports and or feature requests.

## Add to Source Code

To directly contribute to the source code of ```PySPI```, please fork the Github repository, add the changes to one of the branches in your forked repository and then create a [pull request to the master of the main repository](https://github.com/BjoernBiltzinger/pyspi/pulls) from this branch. Code contribution is welcome for different topics:

### Add Functionality

If ```PySPI``` is missing some functionality that you need, you can either create an issue in the Github repository or add it to the code and create a pull request. Always make sure that the old tests do not break and adjust them if needed. Also please add tests and documentation for the new functionality in the pyspi/test folder. This ensures that the functionality will not get broken by future changes to the code and other people will know that this feature exists.

### Code Improvement

You can also contribute code improvements, like making calculations faster or improve the style of the code. Please make sure that the results of the software do not change in this case.

### Bug Fixes

Fixing bugs that you found or that are mentioned in one of the issues is also a good way to contribute to ```PySPI```. Please also make sure to add tests for your changes to check that the bug is gone and that the bug will not recur in future versions of the code.

### Documentation

Additions or examples, tutorials, or better explanations are always welcome. To ensure that the documentation builds with the current version of the software, we are using [jupytext](https://jupytext.readthedocs.io/en/latest/) to write the documentation in Markdown. These are automatically converted to and executed as jupyter notebooks when changes are pushed to Github. 

## Testing

If one wants to run the test suite, simply install `pytest` and `pytest-cov`, then run

```bash
pytest -v

```

in the top level directory. 


---
title: 'PySPI: A python analysis framework for INTEGRAL/SPI'
tags:
  - Python
  - Astronomy
  - Gamma-Ray Bursts
  - INTEGRAL/SPI
authors:
  - name: Björn Biltzinger
    orcid: 0000-0003-3381-0985
    affiliation: "1, 2"
  - name: J. Michael Burgess
    orcid: 0000-0003-3345-9515
    affiliation: "1"
  - name: Thomas Siegert
    orcid: 0000-0002-0552-3535
    affiliation: "1"
bibliography: paper.bib
affiliations:
 - name: Max Planck Institute for Extraterrestrial Physics, Giessenbachstrasse 1, 85748 Garching, Germany
   index: 1
 - name: Technical University of Munich, Boltzmannstrasse 2, 85748 Garching, Germany
   index: 2
date: "27 October 2021"
---

# Summary

PySPI is a newly developed pure python analysis framework for 
Gamma-Ray Burst (GRB) data from the spectrometer (SPI) onboard the 
International Gamma-Ray Astrophysics Laboratory (INTEGRAL). The INTEGRAL 
satellite is a gamma-ray observatory hosting four instruments that 
operate in the energy range between 3 keV and 10 MeV. It was launched in 
2002 and is still working today. The main goals of PySPI are to provide 
an easy to install and develop analysis software for SPI, which includes 
improvements on the statistical analysis of GRB data.

At the moment PySPI is designed for transient sources. One interesting example of 
transient sources are GRBs, which are extremely bright but short flashes 
of Gamma-Rays, with a typical duration between a few ms and a few hundred seconds. They are
believed to be produced by the collapse of massive stars and
mergers of compact objects, like for example neutron stars. In the future we plan to add support for other types of sources than transients, such as persistent point sources as well as extended emission.


# Statement of need

The main analysis tool to analyze SPI data up to now is the
"Off-line Scientific Analysis" (OSA) [@osa], which is maintained by
the INTEGRAL Science Data Centre (ISDC). While it is comprehensive
in its capabilities for manipulating data obtained from all
instrument on-board INTEGRAL, it exists as an IDL interface to a
variety of low-level C++ libraries and is very difficult to install on modern computers. 
While there are containerized 
versions of OSA now available, the modern workflow of simply installing 
the software from a package manager and running on a local workstation is 
not possible and often students rely on a centralized installation which 
must be maintained by a seasoned expert. Moreover, adding more sophisticated 
and/or correct data analysis methods to the software requires an expertise 
that is not immediately accessible to junior researchers or non-experts in 
the installation of OSA. Also due to the
increased computational power that is available today compared to that
of 20 years ago, many of the analysis methods can be improved.

PySPI addresses both these problems: It is providing an easy to install 
software, that can be developed further by everyone who wants to participate. 
It also allows Bayesian fits of the data with true forward folding of the physical spectra 
into the data space via the response. This improves the sensitivity 
and the scientific output of GRB analyses with INTEGRAL/SPI. 


# SPectrometer on INTEGRAL (SPI)

SPI is a coded mask instrument covering the energy range between 20 keV and 8 MeV. It consists of a detector plane with 19 Germanium detectors and a mask plane 1.7 meters above the detectors with 3 cm thick tungsten elements. The mask produces a shadow pattern on the detectors depending on the source position. This information can be used to construct an image from an observation. Also SPI has an excellent energy resolution of 2.5 keV at 1.3 MeV, which makes SPI an ideal instrument to analyze fine spectral features, such as lines from radioactive decays [@spi]. 


# Procedure

To analyze GRB data, PySPI accepts inputs such as the time of the GRB
and the spectral energy bins that will be used in an analysis. With this
information, it automatically downloads all the data files required
for a specific analysis and constructs a response 
as well as a time series for the observation that contains the GRB time. The time series 
can be used to select active time intervals for the source, and time intervals before 
and after the GRB signal for background estimation. After this has been done, a plugin 
for `3ML` [@3mlpaper;@3ML] can be constructed. This allows for all the benefits the 3ML 
framework offers like the modeling framework `astromodels` [@astromodels], joint 
fits with other instruments, many different Bayesian samplers and much more. 
In the [documentation](https://pyspi.readthedocs.io/en/latest/) there is an
[example](https://pyspi.readthedocs.io/en/latest/notebooks/grb_analysis/)
for this workflow procedure.


# Acknowledgments

B. Biltzinger acknowledges financial support from the `German Aerospaces Center (Deutsches Zentrum für Luft- und Raumfahrt, DLR)` under FKZ 50 0R 1913. 


# References

---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region -->

# Installation

## Pip
To install PySPI via pip just use
```bash
pip install py-spi
```

## Conda/Mamba

If you have problems installing PySPI within a Conda environment try to create your environment with this command

```bash
conda create -n pyspi -c conda-forge python=3.9 numpy scipy ipython numba astropy matplotlib h5py pandas pytables
```

or for Mamba

```bash
mamba create -n pyspi -c conda-forge python=3.9 numpy scipy ipython numba astropy matplotlib h5py pandas pytables
```

and then run 

```bash
pip install py-spi
```

with the environment activated.

## Github

To install the latest release from Github run
```bash
git clone https://github.com/BjoernBiltzinger/pyspi.git
```
After that first install the packages from the requirement.txt file with
```bash
cd pyspi
pip install -r requirements.txt
```
Now you can install PySPI with
```bash
python setup.py install
```

## Additional Data Files

There are a few large data files for the background model and the response that are not included in the Github repository. To get these data files run the following commands. Here the data folder is downloaded and is moved to a user defined path where this data folder should be stored on your local machine. Here you have to change the /path/to/internal/data to the path you want to use on your local computer. This only needs to be downloaded once and will not change afterwards.
```bash
wget https://grb.mpe.mpg.de/pyspi_datafolder && unzip pyspi_datafolder
mv data /path/to/internal/data && rm -f pyspi_datafolder
```

## Environment Variables

Next you have to set two environment variable. One to define the path to the folder of the external data like the different SPI data files that will be downloaded by PySPI and one to define the path to the internal data folder we downloaded earlier.
```bash
export PYSPI=/path/to/external/datafolder
export PYSPI_PACKAGE_DATA=/path/to/internal/data
```
Here /path/to/external/datafolder is the path to a folder on your local machine, where PySPI should save all the downloaded data needed for the analysis. The data that will be saved into this folder are the SPI data files as well as one housekeeping data file of SPI and one housekeeping data file of INTEGRAL per analyzed GRB. In total this adds up to roughly 30-70 MB per analyzed GRB. 
It is not recommended to use the same path for both environment variables.

You should also add these two line to your bashrc (or similar) file to automatically set these variables in every new terminal.

Now we are ready to go.

## Run Unit Test Locally

PySPI includes unit test to check that non of its functionality break in new versions. These run automatically for every push on GitHub via GitHub Actions. But you can also run the tests locally. To run the test you need to install pytest and pytest-cov.
```bash
pip install pytest pytest-cov
```
After this run
```bash
pytest -v
```
in the top level directory.

<!-- #endregion -->
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Active Detectors

Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

During the life of INTEGRAL/SPI several detectors stopped working correctly and were therefore disabled. In our analysis we need to take this into account, to not include a detector with 0 counts all the time and because the response for the surrounding detectors change when a detector is deactivated. 

With PySPI you can calculate for a given time, which detectors are active and which response version is valid at that time.

```python
time_string = "051212 205010" #"YYMMDD HHMMSS"; astropy time object also possible
```

To get the active single detectors for this time use:

```python
from pyspi.utils.livedets import get_live_dets
ld = get_live_dets(time_string, event_types="single")
print(f"Active detectors: {ld}")
```

It is also possible to plot the same information visually. This shows the detector plane, and all the inactive detectors at the given time are colored red.
```python
from pyspi.io.plotting.spi_display import SPI
s = SPI(time=time_string)
fig = s.plot_spi_working_dets()
```

Also the response version at that time can be calculated.
```python
from pyspi.utils.function_utils import find_response_version
v = find_response_version(time_string)
print(f"Response version number: {v}")
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Fit for the PSD Efficiency


Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

The first thing we need to do, is to specify the time of the GRB. We do this by specifying a astropy time object or a string in the format YYMMDD HHMMSS.
```python
from astropy.time import Time
grbtime = Time("2012-07-11T02:44:53", format='isot', scale='utc')
#grbtime = "120711 024453" # works also
```

Now we want to analyze in total the energy between 20 and 2000 keV. So we have to take into account the spurious events in the Non-PSD events (see electronic noise section). For the energy bins up to 500 keV we will use all the single events and from 500 to 2000 keV, we will only use the PSD events.
```python
import numpy as np
ein = np.geomspace(20,3000,300)
ebounds_sgl = np.geomspace(20,500,30)
ebounds_psd = np.geomspace(500,2000,30)
```

Due to detector failures there are several versions of the response for SPI. Therefore we have to find the version number for the time of the GRB and construct the base response object for this version.
```python
from pyspi.utils.function_utils import find_response_version
from pyspi.utils.response.spi_response_data import ResponseDataRMF
version = find_response_version(grbtime)
print(version)
rsp_base = ResponseDataRMF.from_version(version)
```

Now we can create the response object for detector 0 and set the position of the GRB, which we already know.
```python
from pyspi.utils.response.spi_response import ResponseRMFGenerator
from pyspi.utils.response.spi_drm import SPIDRM
det=0
ra = 94.6783
dec = -70.99905
drm_generator_sgl = ResponseRMFGenerator.from_time(grbtime, 
                                                    det,
                                                    ebounds_sgl, 
                                                    ein,
                                                    rsp_base)
sd_sgl = SPIDRM(drm_generator_sgl, ra, dec)
```

With this we can build a time series and we use all the single events in this case (PSD + non PSD; see section about electronic noise). To be able to convert the time series into 3ML plugins later, we need to assign them a response object.

```python
from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
tsb_sgl = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
                                            det, 
                                            grbtime, 
                                            response=sd_sgl,
                                            sgl_type="both",
                                            )
```

Now we can have a look at the light curves from -50 to 150 seconds around the specified GRB time.

```python
fig = tsb_sgl.view_lightcurve(-50,150)
```

With this we can select the active time and some background time intervals.

```python
active_time = "65-75"
bkg_time1 = "-500--10"
bkg_time2 = "150-1000"
tsb_sgl.set_active_time_interval(active_time)
tsb_sgl.set_background_interval(bkg_time1, bkg_time2)
```

We can check if the selection and background fitting worked by looking again at the light curve

```python
fig = tsb_sgl.view_lightcurve(-50,150)
```

In this example we use three detectors (IDs: 0, 3 and 4). For these three detectors we build the times series, fit the background and construct the SPILikeGRB plugins which we can use in 3ML.

```python
from pyspi.SPILike import SPILikeGRB
from threeML import DataList
spilikes_sgl = []
spilikes_psd = []
for d in [0,3,4]:
    drm_generator_sgl = ResponseRMFGenerator.from_time(grbtime, 
                                                        d,
                                                        ebounds_sgl, 
                                                        ein,
                                                        rsp_base)
    sd_sgl = SPIDRM(drm_generator_sgl, ra, dec)
    tsb_sgl = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{d}", 
                                                d,
                                                grbtime, 
                                                response=sd_sgl,
                                                sgl_type="both",
                                                )
    tsb_sgl.set_active_time_interval(active_time)
    tsb_sgl.set_background_interval(bkg_time1, bkg_time2)

    sl_sgl = tsb_sgl.to_spectrumlike()
    spilikes_sgl.append(SPILikeGRB.from_spectrumlike(sl_sgl,
                                                    free_position=False))
                                                    
    drm_generator_psd = ResponseRMFGenerator.from_time(grbtime, 
                                                        d,
                                                        ebounds_psd, 
                                                        ein,
                                                        rsp_base)
    sd_psd = SPIDRM(drm_generator_psd, ra, dec)
    tsb_psd = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDetPSD{d}", 
                                                d,
                                                grbtime, 
                                                response=sd_psd,
                                                sgl_type="both",
                                                )
    tsb_psd.set_active_time_interval(active_time)
    tsb_psd.set_background_interval(bkg_time1, bkg_time2)
    
    sl_psd = tsb_psd.to_spectrumlike()
    spilikes_psd.append(SPILikeGRB.from_spectrumlike(sl_psd,
                                                    free_position=False))
                                                    
datalist = DataList(*spilikes_sgl, *spilikes_psd)
```

Now we set a nuisance parameter for the 3ML fit. Nuisance parameter are parameters that only affect one plugin. In this case it is the PSD efficiency for every plugin that uses only PSD events. We do not link the PSD efficiencies in this case, so we determine the PSD efficiency per detector.

```python
for i, s in enumerate(spilikes_psd):
    s.use_effective_area_correction(0,1)
```

Now we have to specify a model for the fit. We use [astromodels](https://astromodels.readthedocs.io/en/latest/) for this.

```python
from astromodels import *
band = Band()
band.K.prior = Log_uniform_prior(lower_bound=1e-6, upper_bound=1e4)
band.K.bounds = (1e-6, 1e4)
band.alpha.set_uninformative_prior(Uniform_prior)
band.beta.set_uninformative_prior(Uniform_prior)
band.xp.prior = Uniform_prior(lower_bound=10,upper_bound=8000)
band.piv.value = 200
ps = PointSource('GRB',ra=ra, dec=dec, spectral_shape=band)

model = Model(ps)
```

Everything is ready to fit now! We make a Bayesian fit here with MultiNest. To use MultiNest you need to install [pymultinest](https://github.com/JohannesBuchner/PyMultiNest) according to its [documentation](https://johannesbuchner.github.io/PyMultiNest/install.html). 

```python
from threeML import BayesianAnalysis
import os
os.mkdir("./chains_psd_eff")
ba_spi = BayesianAnalysis(model, datalist)
for i, s in enumerate(spilikes_psd):
    s.use_effective_area_correction(0,1)
ba_spi.set_sampler("multinest")

ba_spi.sampler.setup(500,
                    chain_name='./chains_psd_eff/docsfit1_',
                    resume=False,
                    verbose=False)
ba_spi.sample()
```

We can use the 3ML features to create a corner plot for this fit:

```python tags=["nbsphinx-thumbnail"]
from threeML.config.config import threeML_config
threeML_config.bayesian.corner_style.show_titles = False
fig = ba_spi.results.corner_plot(components=["cons_SPIDetPSD0", "cons_SPIDetPSD3", "cons_SPIDetPSD4"])
```

So we see we have a PSD efficiency of ~60 +/- 10 % in this case. 
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Contributing 

Contributions to ```PySPI``` are always welcome. They can come in the form of:

## Issues

Please use the [Github issue tracking system for any bugs](https://github.com/BjoernBiltzinger/pyspi/issues), for questions, bug reports and or feature requests.

## Add to Source Code

To directly contribute to the source code of ```PySPI```, please fork the Github repository, add the changes to one of the branches in your forked repository and then create a [pull request to the master of the main repository](https://github.com/BjoernBiltzinger/pyspi/pulls) from this branch. Code contribution is welcome for different topics:

### Add Functionality

If ```PySPI``` is missing some functionality that you need, you can either create an issue in the Github repository or add it to the code and create a pull request. Always make sure that the old tests do not break and adjust them if needed. Also please add tests and documentation for the new functionality in the pyspi/test folder. This ensures that the functionality will not get broken by future changes to the code and other people will know that this feature exists.

### Code Improvement

You can also contribute code improvements, like making calculations faster or improve the style of the code. Please make sure that the results of the software do not change in this case.

### Bug Fixes

Fixing bugs that you found or that are mentioned in one of the issues is also a good way to contribute to ```PySPI```. Please also make sure to add tests for your changes to check that the bug is gone and that the bug will not recur in future versions of the code.

### Documentation

Additions or examples, tutorials, or better explanations are always welcome. To ensure that the documentation builds with the current version of the software, we are using [jupytext](https://jupytext.readthedocs.io/en/latest/) to write the documentation in Markdown. These are automatically converted to and executed as jupyter notebooks when changes are pushed to Github. 


---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Response

Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

For every analysis of SPI data we need the correct response for the observation, which is the connection between physical spectra and detected counts. Normally the response is a function of the position of the source in the satellite frame, the input energies of the physical spectrum and the output energy bins of the experiment. For SPI, there is also a time dependency, because a few detectors failed during the mission time and this changed the response of the surrounding detectors.

In PySPI we construct the response from the official IRF and RMF files, which we interpolate for a given source position and user chosen input and output energy bins.

We start by defining a time, for which we want to construct the response, to get the pointing information of the satellite at this time and the version number of the response. 

```python
from astropy.time import Time
rsp_time = Time("2012-07-11T02:42:00", format='isot', scale='utc')
```

Next we define the input and output energy bins for the response.

```python
import numpy as np
ein = np.geomspace(20,8000,1000)
ebounds = np.geomspace(20,8000,100)
```

Get the response version and construct the rsp base, which is an object holding all the information of the IRF and RMF for this response version. We use this, because if we want to combine many observations later, we don't want to read in this for every observation independently, because this would use a lot of memory. Therefore all the observations with the same response version can share this rsp_base object.

```python
from pyspi.utils.function_utils import find_response_version
from pyspi.utils.response.spi_response_data import ResponseDataRMF
version = find_response_version(rsp_time)
print(version)
rsp_base = ResponseDataRMF.from_version(version)
```

Now we can construct the response for a given detector and source position (in ICRS coordinates)

```python
from pyspi.utils.response.spi_response import ResponseRMFGenerator
from pyspi.utils.response.spi_drm import SPIDRM
det = 0
ra = 94.6783
dec = -70.99905
drm_generator = ResponseRMFGenerator.from_time(rsp_time,
                                                det,
                                                ebounds, 
                                                ein,
                                                rsp_base)
sd = SPIDRM(drm_generator, ra, dec)
```

SPIDRM is a child class of [InstrumentResponse](https://threeml.readthedocs.io/en/stable/api/threeML.utils.OGIP.response.html#threeML.utils.OGIP.response.InstrumentResponse) from threeML, therefore we can use the plotting functions from 3ML.

```python
fig = sd.plot_matrix()
```

---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.2'
      jupytext_version: 1.7.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Access the Underlying Data

Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

Sometime you maybe want to access the underlying data of the analysis to do your own analysis or tests with this data. This section shows how to access some basic quantities, like for example the detected counts per energy channel and the response matrix. First we have to initialize the usual objects in PySPI. 

```python
from astropy.time import Time
import numpy as np
from pyspi.utils.function_utils import find_response_version
from pyspi.utils.response.spi_response_data import ResponseDataRMF
from pyspi.utils.response.spi_response import ResponseRMFGenerator
from pyspi.utils.response.spi_drm import SPIDRM
from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
from pyspi.SPILike import SPILikeGRB

grbtime = Time("2012-07-11T02:44:53", format='isot', scale='utc')
ein = np.geomspace(20,800,300)
ebounds = np.geomspace(20,400,30)
version = find_response_version(grbtime)
rsp_base = ResponseDataRMF.from_version(version)
det=0
ra = 94.6783
dec = -70.99905
drm_generator = ResponseRMFGenerator.from_time(grbtime, 
                                                det,
                                                ebounds, 
                                                ein,
                                                rsp_base)
sd = SPIDRM(drm_generator, ra, dec)
tsb = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
    det, 
    grbtime, 
    response=sd,
    sgl_type="both",
    )
active_time = "65-75"
bkg_time1 = "-500--10"
bkg_time2 = "150-1000"
tsb.set_active_time_interval(active_time)
tsb.set_background_interval(bkg_time1, bkg_time2)
sl = tsb.to_spectrumlike()
plugin = SPILikeGRB.from_spectrumlike(sl,free_position=False)
```

In the following it is listed how you can access some of the basic underlying data.

## Response Matrix

Get response matrix and plot the response for one incoming energy.

```python
import matplotlib.pyplot as plt
ein_id = 200
matrix = sd.matrix

fig, ax = plt.subplots(1,1)
ax.step(ebounds[1:], matrix[:,ein_id])
ax.set_title(f"Response for Ein={round(ein[ein_id], 1)} keV")
ax.set_xlabel("Detected Energy [keV]")
ax.set_ylabel("Effective Area [$cm^2$]")
ax.set_yscale("log");
```

## Event Data

The data is saved as time tagged events. You can access the arrival time and reconstructed energy bin of every photons. It is important to keep in mind that the reconstructed energy is not the true energy, it is just the energy assigned to one of the energy channels.

```python
#arrival times (time in seconds relative to given trigger time)
arrival_times = tsb.time_series.arrival_times

#energy bin of the events
energy_bin = tsb.time_series.measurement
```

## Lightcurve Data

With the event data you can create the lightcurves manually

```python
# plot lightcurves for all echans summed together
bins = np.linspace(-100,200,300)
cnts, bins = np.histogram(arrival_times, bins=bins)

fig, ax = plt.subplots(1,1)
ax.step(bins[1:], cnts)
ax.set_xlabel("Time [s]")
ax.set_ylabel("Counts [cnts]")
ax.set_title("Lightcurve")
ax.legend();
```

## Observed Data Active Time

Get the observed data of the active time and background time selection

```python
# counts
active_time_counts = plugin.observed_counts
estimated_background_counts = plugin.background_counts

# exposure
exposure = plugin.exposure

fig, ax = plt.subplots(1,1)
ax.step(ebounds[1:], active_time_counts/exposure, label="Data")
ax.step(ebounds[1:], estimated_background_counts/exposure, label="Bkg Estimation")

ax.set_xlabel("Detected Energy [keV]")
ax.set_ylabel("Count Rate [cnts/s]")
ax.set_yscale("log")
ax.legend();
```



---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Electronic Noise Region

Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

Since shortly after the launch of INTEGRAL it is known that there are spurious events in the SPI data around ~1.5 MeV. A paper from [Roques & Jourdain](https://arxiv.org/abs/1811.06391) gives an explanation for this problem. Luckily this problem exists only in the events that only triggered the analog front-end electronics (AFEE). The events that trigger in addition the pulse shape discrimination electronics (PSD) do not show this problem. According to [Roques & Jourdain](https://arxiv.org/abs/1811.06391), one should therefore use the PSD events whenever this is possible, which is for events between ~500 and 2500 keV (the precise boundaries were changed during the mission a few times). In the following the events that trigger both the AFEE and PSD are called "PSD events" and the other normal "single events" or "Non-PSD events", even thought the PSD events are of course also single events.

To account for this problem in our analysis we can construct plugins for the "PSD events" and the for the "Non-PSD events" and use only the events with the correct flags, when we construct the time series.

Let's check the difference between the PSD and the Non-PSD events, to see the effect in real SPI data. 

First we define the time and the energy bins we want to use. Then we construct the time series for the three cases:

1. Only the events that trigger AFEE and not PSD
2. Only the events that trigger AFEE and PSD
3. All the single events

```python
from astropy.time import Time
import numpy as np
from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
grbtime = Time("2012-07-11T02:44:53", format='isot', scale='utc')
ebounds = np.geomspace(20,8000,300)
det = 0

from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
tsb_sgl = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
    det,  
    grbtime, 
    ebounds=ebounds,
    sgl_type="sgl",
    )
    
tsb_psd = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
    det,  
    grbtime, 
    ebounds=ebounds,
    sgl_type="psd",
    )

tsb_both = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
    det, 
    grbtime, 
    ebounds=ebounds,
    sgl_type="both",
    )
```

We can check the light curves for all three cases.

```python
print("Only AFEE:")
fig = tsb_sgl.view_lightcurve(-100,300)
```
```python
print("AFFE and PSD trigger:")
fig = tsb_psd.view_lightcurve(-100,300)
```
```python
print("Both Combined:")
fig = tsb_both.view_lightcurve(-100,300)
```

We can see that the PSD event light curve has way less counts. This is due to the fact, that the PSD trigger only starts detecting photons with energies >~ 400 keV.

Next we can get the time integrated counts per energy channel.

```python
tstart = -500
tstop = 1000
counts_sgl = tsb_sgl.time_series.count_per_channel_over_interval(tstart, tstop)
counts_psd = tsb_psd.time_series.count_per_channel_over_interval(tstart, tstop)
counts_both = tsb_both.time_series.count_per_channel_over_interval(tstart, tstop)
```

We can now plot the counts as a function of the energy channel energies
```python
import matplotlib.pyplot as plt
fig, ax = plt.subplots(1,1)
ax.step(ebounds[1:], counts_sgl, label="Only AFEE")
ax.step(ebounds[1:], counts_psd, label="AFEE and PSD")
ax.step(ebounds[1:], counts_both, label="All")
ax.set_xlabel("Detected Energy [keV]")
ax.set_ylabel("Counts")
ax.set_xlim(20,3500)
ax.set_yscale("log")
ax.legend();
```

Several features are visible. 

1. A sharp cutoff at small energies for the PSD events, which is due to the low energy threshold in the PSD electronics. 
2. For energies>~2700 keV the PSD events decrease again faster than the other events.
3. In the Non-PSD events we see a peak at ~ 1600 keV that is not visible in the PSD events. This is the so called electronic noise, which consists of spurious events.
4. The fraction of PSD events to all single events between ~500 and ~2700 keV is very stable and can be explained by an additional dead time for the PSD electronics.

---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Analyse GRB data


Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

The first thing we need to do, is to specify the time of the GRB. We do this by specifying a astropy time object or a string in the format YYMMDD HHMMSS.
```python
from astropy.time import Time
grbtime = Time("2012-07-11T02:44:53", format='isot', scale='utc')
#grbtime = "120711 024453" # works also
```

Next, we need to specify the output and input energy bins we want to use.
```python
import numpy as np
ein = np.geomspace(20,800,300)
ebounds = np.geomspace(20,400,30)
```

Due to detector failures there are several versions of the response for SPI. Therefore we have to find the version number for the time of the GRB and construct the base response object for this version.
```python
from pyspi.utils.function_utils import find_response_version
from pyspi.utils.response.spi_response_data import ResponseDataRMF
version = find_response_version(grbtime)
print(version)
rsp_base = ResponseDataRMF.from_version(version)
```

Now we can create the response object for detector 0 and set the position of the GRB, which we already know.
```python
from pyspi.utils.response.spi_response import ResponseRMFGenerator
from pyspi.utils.response.spi_drm import SPIDRM
det=0
ra = 94.6783
dec = -70.99905
drm_generator = ResponseRMFGenerator.from_time(grbtime, 
                                                det,
                                                ebounds, 
                                                ein,
                                                rsp_base)
sd = SPIDRM(drm_generator, ra, dec)
```

With this we can build a time series and we use all the single events in this case (PSD + non PSD; see section about electronic noise). To be able to convert the time series into 3ML plugins later, we need to assign them a response object.
```python
from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
tsb = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
    det, 
    grbtime, 
    response=sd,
    sgl_type="both",
    )
```

Now we can have a look at the light curves from -50 to 150 seconds around the specified GRB time.
```python
fig = tsb.view_lightcurve(-50,150)
```

With this we can select the active time and some background time intervals.
```python
active_time = "65-75"
bkg_time1 = "-500--10"
bkg_time2 = "150-1000"
tsb.set_active_time_interval(active_time)
tsb.set_background_interval(bkg_time1, bkg_time2)
```
We can check if the selection and background fitting worked by looking again at the light curve
```python tags=["nbsphinx-thumbnail"]
fig = tsb.view_lightcurve(-50,150)
```
For the fit we of course want to use all the available detectors. So we first check which detectors were still working at that time.
```python
from pyspi.utils.livedets import get_live_dets
active_dets = get_live_dets(time=grbtime, event_types=["single"])
print(active_dets)
```

Now we loop over these detectors, build the times series, fit the background and construct the SPILikeGRB plugins which we can use in 3ML.
```python
from pyspi.SPILike import SPILikeGRB
from threeML import DataList
spilikes = []
for d in active_dets:
    drm_generator = ResponseRMFGenerator.from_time(grbtime, 
                                                    d,
                                                    ebounds, 
                                                    ein,
                                                    rsp_base)
    sd = SPIDRM(drm_generator, ra, dec)
    tsb = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{d}", 
                                                d,
                                                grbtime, 
                                                response=sd,
                                                sgl_type="both",
                                                )
    tsb.set_active_time_interval(active_time)
    tsb.set_background_interval(bkg_time1, bkg_time2)

    sl = tsb.to_spectrumlike()
    spilikes.append(SPILikeGRB.from_spectrumlike(sl,
                                                free_position=False))
datalist = DataList(*spilikes)
```

Now we have to specify a model for the fit. We use [astromodels](https://astromodels.readthedocs.io/en/latest/) for this.
```python
from astromodels import *
pl = Powerlaw()
pl.K.prior = Log_uniform_prior(lower_bound=1e-6, upper_bound=1e4)
pl.K.bounds = (1e-6, 1e4)
pl.index.set_uninformative_prior(Uniform_prior)
pl.piv.value = 200
ps = PointSource('GRB',ra=ra, dec=dec, spectral_shape=pl)

model = Model(ps)
```

Everything is ready to fit now! We make a Bayesian fit here with emcee
```python
from threeML import BayesianAnalysis
ba_spi = BayesianAnalysis(model, datalist)
ba_spi.set_sampler("emcee", share_spectrum=True)
ba_spi.sampler.setup(n_walkers=20, n_iterations=500)
ba_spi.sample()
```

We can inspect the fits with residual plots

```python
from threeML import display_spectrum_model_counts
fig = display_spectrum_model_counts(ba_spi, 
                                data_per_plot=5, 
                                source_only=True,
                                show_background=False,
                                model_cmap="viridis", 
                                data_cmap="viridis",
                                background_cmap="viridis")
```

and have a look at the spectrum

```python
from threeML import plot_spectra
fig = plot_spectra(ba_spi.results, flux_unit="keV/(s cm2)", ene_min=20, ene_max=600)
```
We can also get a summary of the fit and write the results to disk (see 3ML documentation).


It is also possible to localize GRBs with PySPI, to this we simply free the position of point source with:

```python
for s in spilikes:
    s.set_free_position(True)
datalist = DataList(*spilikes)
```
Initialize the Bayesian Analysis and start the sampling with MultiNest. To use MultiNest you need to install [pymultinest](https://github.com/JohannesBuchner/PyMultiNest) according to its [documentation](https://johannesbuchner.github.io/PyMultiNest/install.html). 
```python
import os
os.mkdir("./chains_grb_example")
ba_spi = BayesianAnalysis(model, datalist)
ba_spi.set_sampler("multinest")
ba_spi.sampler.setup(500, 
                    chain_name='./chains_grb_example/docsfit1_',
                    resume=False,
                    verbose=False)
ba_spi.sample()
```

We can use the 3ML features to create a corner plot for this fit:

```python
from threeML.config.config import threeML_config
threeML_config.bayesian.corner_style.show_titles = False
fig = ba_spi.results.corner_plot(components=["GRB.position.ra", "GRB.position.dec"])
```

When we compare the results for ra and dec, we can see that this matches with the position from [Swift-XRT for the same GRB (RA, Dec = 94.67830, -70.99905)](https://gcn.gsfc.nasa.gov/gcn/other/120711A.gcn3)
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.13.0
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Light Curves

Setup to make the output clean for the docs:
```python
%%capture
from threeML import silence_logs
import warnings
warnings.filterwarnings("ignore")
silence_logs()
import matplotlib.pyplot as plt
%matplotlib inline
from jupyterthemes import jtplot
jtplot.style(context="talk", fscale=1, ticks=True, grid=False)
```

Gamma-Ray Bursts are transient sources with a typical duration between milliseconds and a few tens of seconds. Therefore they are nicely visible in light curves. In the following we will see how we can get the light curve of a real GRB as seen by an INTEGRAL/SPI detector.

First we have to define the rough time of the GRB.
```python
from astropy.time import Time
grbtime = Time("2012-07-11T02:44:53", format='isot', scale='utc')
```

Next we need to define the bounds of the energy bins we want to use.

```python
import numpy as np
ebounds = np.geomspace(20,8000,100)
```

Now we can construct the time series.

```python
from pyspi.utils.data_builder.time_series_builder import TimeSeriesBuilderSPI
det = 0
tsb = TimeSeriesBuilderSPI.from_spi_grb(f"SPIDet{det}", 
                                        det, 
                                        grbtime,
                                        ebounds=ebounds, 
                                        sgl_type="both",
                                        )
```

We can now plot the light curves for visualization, in which we can clearly see a transient source in this case.

```python
fig = tsb.view_lightcurve(-50,250)
```
.. pySPI documentation master file, created by
   sphinx-quickstart on Sun Feb  4 11:24:43 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PySpi's documentation!
====================================
.. image:: media/pypsi_logo2.png

PySPI is pure python interface to analyze Gamma-Ray Burst (GRB) data from the spectrometer (SPI) onboard the International Gamma-Ray Astrophysics Laboratory (INTEGRAL). The INTEGRAL satellite is a gamma-ray observatory hosting four instruments that operate in the energy range between 3 keV and 10 MeV. It was launched in 2002 and is still working today. The main goals of PySPI are to provide an easy to install and develop analysis software for SPI, which includes improvements on the statistical analysis of GRB data. At the moment PySPI is designed for transient sources, like Gamma Ray Bursts (GRBs). In the future we plan to add support for other types of sources, such as persistent point sources as well as extended emission.

Comparison to OSA
------------------------------------

The main analysis tool to analyze SPI data up to now is the “Off-line Scientific Analysis” (OSA) `Chernyakova et al., 2020 <https://www.isdc.unige.ch/integral/download/osa/doc/11.1/osa_um_intro/man.html>`__), which is maintained by the INTEGRAL Science Data Centre (ISDC). While it is comprehensive in its capabilities for manipulating data obtained from all instrument on-board INTEGRAL, it exists as an IDL interface to a variety of low-level C++ libraries and is very difficult to install on modern computers. While there are containerized versions of OSA now available, the modern workflow of simply installing the software from a package manager and running on a local workstation is not possible and often students rely on a centralized installation which must be maintained by a seasoned expert. Moreover, adding more sophisticated and/or correct data analysis methods to the software requires an expertise that is not immediately accessible to junior researchers or non-experts in the installation of OSA. Also due to the increased computational power that is available today compared to that of 20 years ago, many of the analysis methods can be improved. PySPI addresses both these problems: It is providing an easy to install software, that can be developed further by everyone who wants to contribute. It also allows Bayesian fits of the data with true forward folding of the physical spectra into the data space via the response. This improves the sensitivity and the scientific output of GRB analyses with INTEGRAL/SPI.

Multi Mission Analysis
------------------------------------

PySPI provides a plugin for `3ML <https://threeml.readthedocs.io/en/stable>`__. This makes multi missions analysis with other instruments possible. Also all the spectral models from `astromodels <https://astromodels.readthedocs.io/en/latest/>`__ are available for the fits. Check out these two software packages for more information.


.. toctree::
   :maxdepth: 5
   :hidden:

   notebooks/installation.ipynb
   notebooks/time_series.ipynb
   notebooks/response.ipynb
   notebooks/psd.ipynb
   notebooks/active_detectors.ipynb
   notebooks/access_data.ipynb
   notebooks/contributing.ipynb
   api/API

.. nbgallery::
   :caption: Features and examples:

   notebooks/grb_analysis.ipynb
   notebooks/psd_eff.ipynb
{% if referencefile %}
.. include:: {{ referencefile }}
{% endif %}

{{ objname }}
{{ underline }}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: Functions

   .. autosummary::
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: Classes

   .. autosummary::
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: Exceptions

   .. autosummary::
   {% for item in exceptions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
{% if referencefile %}
.. include:: {{ referencefile }}
{% endif %}

{{ objname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :show-inheritance:

   {% if '__init__' in methods %}
     {% set caught_result = methods.remove('__init__') %}
   {% endif %}

   {% block attributes_summary %}
   {% if attributes %}

   .. rubric:: Attributes Summary

   .. autosummary::
      :toctree:
      :template: base.rst
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   {% endif %}
   {% endblock %}

   {% block methods_summary %}
   {% if methods %}

   .. rubric:: Methods Summary

   .. autosummary::
      :toctree:
      :template: base.rst
   {% for item in methods %}
      ~{{ name }}.{{ item }}
   {%- endfor %}

   {% endif %}
   {% endblock %}
{% if referencefile %}
.. include:: {{ referencefile }}
{% endif %}

{{ objname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. auto{{ objtype }}:: {{ objname }}
pyspi
========

.. toctree::
   :maxdepth: 4

   pyspi
API
===

Here you can find the documentation of all classes and methods:

.. include:: modules.rst
