---
title: 'popsynth: A generic astrophysical population synthesis framework'
tags:
  - Python
  - astronomy
  - population synthesis
  - cosmology
authors:
  - name: J. Michael Burgess
    orcid: 0000-0003-3345-9515
    affiliation: "1"
  - name: Francesca Capel
    orcid: 0000-0002-1153-2139
    affiliation: "2"
bibliography: paper.bib
affiliations:
 - name: Max Planck Institute for Extraterrestrial Physics, Giessenbachstrasse, 85748 Garching, Germany
   index: 1
 - name: Technical University of Munich, Boltzmannstrasse 2, 85748 Garching, Germany
   index: 2
date: "07 April 2021"
---

# Summary

Simulating a survey of fluxes and redshifts (distances) from an
astrophysical population is a routine task. `popsynth` provides a
generic, object-oriented framework to produce synthetic surveys from
various distributions and luminosity functions, apply selection
functions to the observed variables and store them in a portable (HDF5)
format. Population synthesis routines can be constructed either using
classes or from a serializable YAML format allowing flexibility and
portability. Users can not only sample the luminosity and distance of
the populations, but they can create auxiliary distributions for
parameters which can have arbitrarily complex dependencies on one
another. Thus, users can simulate complex astrophysical populations
which can be used to calibrate analysis frameworks or quickly test
ideas.

# Statement of need

`popsynth` provides a generic framework for simulating astrophysical
populations with an easily extensible class inheritance scheme that
allows users to adapt the code to their own needs. As understanding
the rate functions of astrophysical populations (e.g., gravitational
wave sources, gamma-ray bursts, active galactic nuclei) becomes an
increasingly important field [@Loredo:2019], researchers develop
various ways to estimate these populations from real data. `popsynth`
provides a way to calibrate these analysis techniques by producing
synthetic data where the inputs are known
[e.g. @Mortlock:2019]. Moreover, selection effects are an important
part of population analysis and the ability to include this property
when generating a population is vital to the calibration of any survey
analysis method which operates on an incomplete sample.

Similar frameworks exist for simulating data from specific catalogs
such as `SkyPy` [@skypy] and `firesong` [@firesong], however, these
have much more focused applications and do not include the ability to
impose selection functions.

# Procedure

Once a rate function and all associated distributions are specified in
`popsynth`, a numeric integral over the rate function produces the
total rate of objects in the populations. A survey is created by
making a draw from a Poisson distribution with mean equal to the total
rate of objects multiplied by survey duration for the number of
objects in the survey. For each object, the properties such as
distance and luminosity are sampled from their associated
distributions. Selection functions are then applied to latent or
observed variables as specified by the user. Finally, all population
objects and variables are returned in an object that can be serialized
to disk for later examination. Further details on the mathematics,
procedure, and details on customization can be found in the extensive
[documentation](https://popsynth.readthedocs.io/).


# Acknowledgments

This project was inspired by conversations with Daniel J. Mortlock
wherein we tried to calibrate an analysis method we will eventually
get around to finishing. Inspiration also came from wanting to
generalize the examples from Will Farr's lecture note
[@selection]. J. Michael Burgess acknowledges support from the
Alexander von Humboldt Stiftung. Francesca Capel acknowledges
financial support from the Excellence Cluster ORIGINS, which is funded
by the Deutsche Forschungsgemeinschaft (DFG, German Research
Foundation) under Germany’s Excellence Strategy
- EXC-2094-390783311.

# References

![CI](https://github.com/grburgess/popsynth/workflows/CI/badge.svg?branch=master)
[![codecov](https://codecov.io/gh/grburgess/popsynth/branch/master/graph/badge.svg)](https://codecov.io/gh/grburgess/popsynth)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/5d02c9e6f5c540989a615eb1575863e3)](https://app.codacy.com/gh/grburgess/popsynth?utm_source=github.com&utm_medium=referral&utm_content=grburgess/popsynth&utm_campaign=Badge_Grade_Settings)
[![Documentation Status](https://readthedocs.org/projects/popsynth/badge/?version=latest)](https://popsynth.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5109590.svg)](https://doi.org/10.5281/zenodo.5109590)
![PyPI](https://img.shields.io/pypi/v/popsynth)
![PyPI - Downloads](https://img.shields.io/pypi/dm/popsynth)
 [![status](https://joss.theoj.org/papers/a52e4c2c355396e7946917996502aac0/status.svg)](https://joss.theoj.org/papers/a52e4c2c355396e7946917996502aac0)
# popsynth

![alt text](https://raw.githubusercontent.com/grburgess/popsynth/master/external/logo.png)

`popsynth` core function is to create **observed** surveys from **latent** population models. 

First, let's define what a population of objects is in terms of a
generative model. The two main ingredients are the objects' spatial
distribution (<img src="https://render.githubusercontent.com/render/math?math=\lambda(\vec{r},\vec{\psi})">) and the distribution of
their inherent properties (<img src="https://render.githubusercontent.com/render/math?math=\pi(\vec{\phi} | \vec{\psi})">). Here,
<img src="https://render.githubusercontent.com/render/math?math=\vec{\psi}"> are the latent population parameters, <img src="https://render.githubusercontent.com/render/math?math=\vec{r}"> are the
spatial locations of the objects, and <img src="https://render.githubusercontent.com/render/math?math=\vec{\phi}"> are the properties
of the individual objects (luminosity, spin, viewing angle, mass,
etc.). The spatial distribution is defined such that:

<img src="https://render.githubusercontent.com/render/math?math=\frac{d \Lambda}{dt}(\vec{\psi}) = \int d r \frac{dV}{dr} \lambda(\vec{r}, \vec{\psi}))">

is the intensity of objects for a given set of population
parameters. With these definitions we can define the probability for
an object to have position <img src="https://render.githubusercontent.com/render/math?math=\vec{r}"> and properties <img src="https://render.githubusercontent.com/render/math?math=\vec{\phi}"> as

<img src="https://render.githubusercontent.com/render/math?math=\pi(\vec{r}, \vec{\phi} | \vec{\psi}) = \frac{\lambda(\vec{r}, \vec{\psi})  \pi(\vec{\phi} | \vec{\psi})}{ \int d r \frac{dV}{dr} \lambda(\vec{r}, \vec{\psi})}">

`popsynth` allows you to specify these spatial and property
distributions in an object-oriented way to create surveys. The final
ingredient to creating a sample for a survey is knowing how many
objects to sample from the population (before any selection effects
are applied). Often, we see this number in simulation frameworks
presented as "we draw N objects to guarantee we have enough." This is
incorrect. A survey takes place over a given period of time (<img src="https://render.githubusercontent.com/render/math?math=\Delta t">) in which observed objects are counted. This is a description of a
Poisson process. Thus, the number of objects in a simulation of this
survey is a draw from a Poisson distribution:

<img src="https://render.githubusercontent.com/render/math?math=N \sim Poisson \left( \Delta t \frac{d\Lambda}{dt} \right)">

Thus, ```popsynth``` first numerically integrates the spatial
distribution to determine the Poisson rate parameter for the given
$\vec{\psi}$, then makes a Poisson draw for the number of objects in
the population survey. For each object, positions and properties are
drawn with arbitrary dependencies between them. Finally, selection
functions are applied to either latent or observed (with or without
measurement error) properties.


**Note:** If instead we draw a preset number of objects, as is done in
many astrophysical population simulation frameworks, it is equivalent
to running a survey up until that specific number of objects is
detected. This process is distributed as a negative binomial process,
i.e, wait for a number of successes and requires a different
statistical framework to compare models to data.


## Installation
```bash
pip install popsynth
```


Note: **This is not synth pop!** If you were looking for some hard driving beats out of a yamaha keyboard with bells... look elsewhere

![alt text](https://raw.githubusercontent.com/grburgess/popsynth/master/external/pop.gif)


## Contributing

Contributions to ```popsynth``` are always welcome. They can come in the form of:

### Bug reports

Please use the [Github issue tracking system for any
bugs](https://github.com/grburgess/popsynth/issues), for questions,
and or feature requests.

### Code and more distributions

While it is easy to create custom distributions in your local setup,
if you would like to add them to popsynth directly, go ahead. Please
include tests to ensure that your contributions are compatible with
the code and can be maintained in the long term.

### Documentation

Additions or examples, tutorials, or better explanations are always
welcome. To ensure that the documentation builds with the current
version of the software, I am using
[jupytext](https://jupytext.readthedocs.io/en/latest/) to write the
documentation in Markdown. These are automatically converted to and
executed as jupyter notebooks when changes are pushed to Github.


---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Core Concept
`popsynth` core function is to create **observed** surveys from **latent** population models. 

First, let's define what a population of objects is in terms of a
generative model. The two main ingredients are the objects' spatial
distribution ($\lambda(\vec{r}; \vec{\psi})$) and the distribution of
their inherent properties ($\pi(\vec{\phi} | \vec{\psi})$). Here,
$\vec{\psi}$ are the latent population parameters, $\vec{r}$ are the
spatial locations of the objects, and $\vec{\phi}$ are the properties
of the individual objects (luminosity, spin, viewing angle, mass,
etc.). The spatial distribution is defined such that:

$$\frac{d \Lambda}{dt}(\vec{\psi}) = \int d r \frac{dV}{dr} \lambda(\vec{r}; \vec{\psi})$$

is the intensity of objects for a given set of population
parameters. With these definitions we can define the probability for
an object to have position $\vec{r}$ and properties $\vec{\phi}$ as

$$\pi(\vec{r}, \vec{\phi} | \vec{\psi}) = \frac{\lambda(\vec{r}; \vec{\psi})  \pi(\vec{\phi} | \vec{\psi})}{ \int d r \frac{dV}{dr} \lambda(\vec{r}; \vec{\psi})} $$

`popsynth` allows you to specify these spatial and property
distributions in an object-oriented way to create surveys. The final
ingredient to creating a sample for a survey is knowing how many
objects to sample from the population (before any selection effects
are applied). Often, we see this number in simulation frameworks
presented as "we draw N objects to guarantee we have enough." This is
incorrect. A survey takes place over a given period of time ($\Delta
t$) in which observed objects are counted. This is a description of a
Poisson process. Thus, the number of objects in a simulation of this
survey is a draw from a Poisson distribution:

$$N \sim \mathrm{Poisson}\left(\Delta t \frac{d\Lambda}{dt}\right) \mathrm{.}$$

Thus, `popsynth` first numerically integrates the spatial
distribution to determine the Poisson rate parameter for the given
$\vec{\psi}$, then makes a Poisson draw for the number of objects in
the population survey. For each object, positions and properties are
drawn with arbitrary dependencies between them. Finally, selection
functions are applied to either latent or observed (with or without
measurement error) properties.


**Note:** If instead we draw a preset number of objects, as is done in
many astrophysical population simulation frameworks, it is equivalent
to running a survey up until that specific number of objects is
detected. This process is distributed as a negative binomial process,
i.e, wait for a number of successes and requires a different
statistical framework to compare models to data.

In the following, the process for constructing distributions and
populations is described.
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

<!-- #region -->
# Installation

```popsynth``` can be installed with

```bash
pip install popsynth
```
Alternatively, one can install via git

```bash
git clone https://github.com/grburgess/popsynth
cd popsynth
python setup.py install
```


In order to produce graphs of populations, one needs the optional [graphiz](https://graphviz.readthedocs.io/en/stable/) library.
The 3D plots seen in the documentation require properly setting up the [ipyvolume](https://ipyvolume.readthedocs.io/en/latest/) package and configuring it for your local setup. 

<!-- #endregion -->

```python

```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Stellar Mass-Luminosity Bias

Suppose that stars have a mass-luminosity relationship such that $L
\propto M^3$. If we have a flux-limited survey, it will bias us
towards observing more massive stars which are not representative of
the full mass distribution. Let's see how to set this up in `popsynth`. 

## Setup the problem
First, we will assume that we have some initial mass function (IMF)
for our stars that describes how there masses are distributed. For
simplicity, we will assume that this IMF is just a log-normal distribution:

```python

%matplotlib inline

import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import networkx as nx
import numpy as np
import warnings

warnings.simplefilter("ignore")

import popsynth

# create a sampler for mass
# we do not directly observe the mass as it is a latent quantity

initial_mass_function = popsynth.LogNormalAuxSampler(name="mass", observed = False)

```

We now assume the dependent variable is the luminosity, so we need a
`DerivedLumAuxSampler` that generates luminosities given a mass:

```python
class MassLuminosityRelation(popsynth.DerivedLumAuxSampler):
    _auxiliary_sampler_name = "MassLuminosityRelation"
  
    def __init__(self, mu=0.0, tau=1.0, sigma=1):
        # this time set observed=True
        super(MassLuminosityRelation, self).__init__("mass_lum_relation", uses_distance=False)

    def true_sampler(self, size):
        
		# the secondary quantity is mass 
		
        mass = self._secondary_samplers["mass"].true_values
        
		# we will store the log of mass cubed
        self._true_values = 3 * np.log10(mass)

    def compute_luminosity(self):
        # compute the luminosity
		# from the relation
        return np.power(10., self._true_values)


luminosity = MassLuminosityRelation()


```

Now we can put everything together. First, we need to assign
`mass` as a secondary quantity to the luminosity

```python
luminosity.set_secondary_sampler(initial_mass_function)
```

Finally, we will use a simple spherical geometry to hold our stars. We
will also put a hard flux limit on our survey to simulate a
flux-limited catalog.

```python
pop_gen = popsynth.populations.SphericalPopulation(1, r_max=10)

# create the flux selection

flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-2
pop_gen.set_flux_selection(flux_selector)

# now add the luminisity sampler

pop_gen.add_observed_quantity(luminosity)

```

Now let's draw our survey.

```python

pop = pop_gen.draw_survey(flux_sigma=0.5)

```

We can now look at the distribution of the masses:


```python tags=["nbsphinx-thumbnail"]
fig, ax = plt.subplots()

bins = np.linspace(0,20,50)

ax.hist(pop.mass,
        bins=bins,
        label='all', 
        color=purple,
        histtype="step",
        lw=2
       
       )

ax.hist(pop.mass[pop.selection],
        bins=bins,
        label='selected',
        color=yellow,
        histtype="step",
        lw=2
       )

ax.set_xlabel('stellar mass')
ax.legend()

```

We can see that indeed our selected masses are biased towards higher values. 

Let's look in the mass-luminostiy plane:

```python
fig, ax = plt.subplots()

bins = np.linspace(0,20,50)

ax.scatter(pop.mass[~pop.selection],
           pop.luminosities_latent[~pop.selection],        
           label='not selected', 
           color=purple,
           alpha=0.5,
           s=10
            
            
       
       )

ax.scatter(pop.mass[pop.selection],
           pop.luminosities_latent[pop.selection],        
           label='selected', 
           color=yellow,
           alpha=0.5,
           s=5
          )



ax.set_xlabel('stellar mass')
ax.set_ylabel('luminosity')
ax.set_yscale('log')


ax.legend()


```

---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Custom distributions and populations

Custom populations can be created either by piecing together existing populations (spatial and luminosity populations) or building them from scratch with distributions.

**popsynth** comes loaded with many combinations of typical population distributions. However, we demonstrate here how to create your own.


## Creating distributions

The population samplers rely on distributions. Each population has an internal spatial and luminosity distribution. For example, lets look at a simple spatial distribution:


```python
%matplotlib inline

import numpy as np
import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"


import popsynth

popsynth.update_logging_level("INFO")
import warnings

warnings.simplefilter("ignore")
```

```python
from popsynth.distribution import SpatialDistribution


class MySphericalDistribution(SpatialDistribution):

    # we need this property to register the class

    _distribution_name = "MySphericalDistribution"

    def __init__(
        self,
        seed=1234,
        form=None,
    ):

        # the latex formula for the ditribution
        form = r"4 \pi r2"

        # we do not need a "truth" dict here because
        # there are no parameters

        super(MySphericalDistribution, self).__init__(
            seed=seed,
            name="sphere",
            form=form,
        )

    def differential_volume(self, r):

        # the differential volume of a sphere
        return 4 * np.pi * r * r

    def transform(self, L, r):

        # luminosity to flux
        return L / (4.0 * np.pi * r * r)

    def dNdV(self, r):

        # define some crazy change in the number/volume for fun

        return 10.0 / (r + 1) ** 2
```

<!-- #region -->
We simply define the differential volume and how luminosity is transformed to flux in the metric. Here, we have a simple sphere out to some *r_max*. We can of course subclass this object and add a normalization.


Now we define a luminosity distribution.
<!-- #endregion -->

```python
from popsynth.distribution import LuminosityDistribution, DistributionParameter


class MyParetoDistribution(LuminosityDistribution):
    _distribution_name = "MyParetoDistribution"

    Lmin = DistributionParameter(default=1, vmin=0)
    alpha = DistributionParameter(default=2)

    def __init__(self, seed=1234, name="pareto"):

        # the latex formula for the ditribution
        lf_form = r"\frac{\alpha L_{\rm min}^{\alpha}}{L^{\alpha+1}}"

        super(MyParetoDistribution, self).__init__(
            seed=seed,
            name="pareto",
            form=lf_form,
        )

    def phi(self, L):

        # the actual function, only for plotting

        out = np.zeros_like(L)

        idx = L >= self.Lmin

        out[idx] = self.alpha * self.Lmin ** self.alpha / L[idx] ** (self.alpha + 1)

        return out

    def draw_luminosity(self, size=1):
        # how to sample the latent parameters
        return (np.random.pareto(self.alpha, size) + 1) * self.Lmin
```


<div class="alert alert-info">

**Note:** If you want to create a cosmological distribution, inherit from from ComologicalDistribution class!

</div>

## Creating a population synthesizer

Now that we have defined our distributions, we can create a population synthesizer that encapsulated them

```python
from popsynth.population_synth import PopulationSynth


class MyPopulation(PopulationSynth):
    def __init__(self, Lmin, alpha, r_max=5, seed=1234):

        # instantiate the distributions
        luminosity_distribution = MyParetoDistribution(seed=seed)

        luminosity_distribution.alpha = alpha
        luminosity_distribution.Lmin = Lmin

        spatial_distribution = MySphericalDistribution(seed=seed)
        spatial_distribution.r_max = r_max

        # pass to the super class
        super(MyPopulation, self).__init__(
            spatial_distribution=spatial_distribution,
            luminosity_distribution=luminosity_distribution,
            seed=seed,
        )
```

```python
my_pop_gen = MyPopulation(Lmin=1, alpha=1, r_max=10)

flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-2

my_pop_gen.set_flux_selection(flux_selector)

population = my_pop_gen.draw_survey()
```

```python
fig = population.display_obs_fluxes_sphere(cmap="magma", background_color="black" ,s=50)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Short GRBS 

In [Ghirlanda et al. 2016](https://arxiv.org/abs/1607.07875) a fitting algorithm was used to determine the redshift and luminosity of short GRBS. We can use the parameters to reproduce the population and the observed GBM survey.

```python
from popsynth import SFRDistribution, BPLDistribution, PopulationSynth, NormalAuxSampler, AuxiliarySampler, HardFluxSelection
from popsynth import update_logging_level
update_logging_level("INFO")
```

```python
%matplotlib inline

import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import networkx as nx
import numpy as np
import warnings

warnings.simplefilter("ignore")
```

In the work, the luminosity function of short GRBs is model as a broken power law.

```python
bpl = BPLDistribution()

bpl.alpha = -0.53
bpl.beta = -3.4
bpl.Lmin = 1e47 # erg/s
bpl.Lbreak = 2.8e52
bpl.Lmax = 1e55

```

To model the redshift distribution, an empirical form from [Cole et al 2001](https://academic.oup.com/mnras/article/326/1/255/1026734?login=true) is used. In ```popsynth``` we call this the ```SFRDistribution``` (but perhaps a better name is needed).

```python
sfr = SFRDistribution()
```

```python
sfr.r0 = 5.
sfr.a = 1
sfr.rise = 2.8
sfr.decay = 3.5
sfr.peak = 2.3
```

We can checkout how the rate changes with redshift

```python
fig, ax = plt.subplots()

z = np.linspace(0,5,100)

ax.plot(z, sfr.dNdV(z), color=purple)
ax.set_xlabel("z")
ax.set_ylabel(r"$\frac{\mathrm{d}N}{\mathrm{d}V}$")
```

<!-- #region -->
In their model, the authors also have some secondary parameters that are connected to the luminosity. These are the  parameters for the spectrum of the GRB. It is proposed that the spectra peak energy (Ep) is linked to the luminosity by a power law relation:


$$ \log E_{\mathrm{p}} \propto a + b \log L$$

We can build an auxiliary sample to simulate this as well. But we will also add a bit of scatter to the intercept of the relation.
<!-- #endregion -->

```python
intercept =NormalAuxSampler(name="intercept", observed=False)

intercept.mu = 0.034
intercept.sigma = .005



```

```python
class EpSampler(AuxiliarySampler):
    
    _auxiliary_sampler_name = "EpSampler"

    def __init__(self):

        # pass up to the super class
        super(EpSampler, self).__init__("Ep", observed=True, uses_luminosity = True)

    def true_sampler(self, size):

        # we will get the intercept's latent (true) value
        # from its sampler
        
        intercept = self._secondary_samplers["intercept"].true_values
        
        slope = 0.84

        self._true_values = np.power(10., intercept + slope * np.log10(self._luminosity/1e52) + np.log10(670.))
        
    def observation_sampler(self, size):
        
        # we will also add some measurement error to Ep
        self._obs_values = self._true_values + np.random.normal(0., 10, size=size)
        
```

Now we can put it all together.

```python
pop_synth = PopulationSynth(spatial_distribution=sfr, luminosity_distribution=bpl)
```

We will have a hard flux selection which is Fermi-GBM's fluz limit of ~ 1e-7 erg/s/cm2

```python
selection = HardFluxSelection()
selection.boundary = 1e-7
```

```python
pop_synth.set_flux_selection(selection)
```

We need to add the Ep sampler. Once we set the intercept sampler as a secondary it will automatically be added to the population synth.

```python
ep = EpSampler()
```

```python
ep.set_secondary_sampler(intercept)
```

```python
pop_synth.add_auxiliary_sampler(ep)
```

We are ready to sample our population. We will add some measurement uncertainty to the fluxes as well.

```python
population = pop_synth.draw_survey(flux_sigma=0.2)
```

```python tags=["nbsphinx-thumbnail"]
population.display_fluxes(true_color=purple, obs_color=yellow, with_arrows=False, s= 5);
```

Let's look at our distribution of Ep

```python
fig, ax = plt.subplots()

ax.hist(np.log10(population.Ep_obs[population.selection]), histtype="step", color=yellow, lw=3, label="Ep observed")
ax.hist(np.log10(population.Ep[~population.selection]), histtype="step", color=purple, lw=3,  label="Ep hidden")
ax.set_xlabel("log Ep")

ax.legend()
```

```python
fig, ax = plt.subplots()


ax.scatter(population.fluxes_observed[~population.selection],
           population.Ep_obs[~population.selection],c=purple, alpha=0.5)
ax.scatter(population.fluxes_observed[population.selection],
           population.Ep_obs[population.selection],c=yellow, alpha=0.5)

ax.set_xscale("log")
ax.set_yscale("log")

ax.set_xlabel("log Ep")
ax.set_xlabel("log Flux")
```

Does this look like the observed catalogs?
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Fun with the Milky Way

While not entirely useful at the moment. There is support for generating simplistic spiral galaxy distribtuions.

```python
import popsynth
import ipyvolume as ipv


from astropy.coordinates import SkyCoord

%matplotlib inline

import matplotlib.pyplot as plt
from jupyterthemes import jtplot

purple = "#B833FF"

popsynth.update_logging_level("INFO")
from popsynth.populations.spatial_populations import MWRadialPopulation
```

```python
ld = popsynth.distributions.pareto_distribution.ParetoDistribution()
ld.alpha = 3
ld.Lmin = 1
```

```python
synth = MWRadialPopulation(rho=1, luminosity_distribution=ld)
```

```python
population = synth.draw_survey()
```

```python
fig = population.display_obs_fluxes_sphere(
    cmap="magma", background_color="black", size=0.1
)
```

---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# BL Lac blazars
A model for the luminosity function and cosmic evolution of BL Lac type blazars is presented in [Ajello et al. 2014](https://arxiv.org/abs/1310.0006), based on observations in gamma-rays with the Fermi-LAT instrument.

We can use the results of this paper to build a BL Lac population that is able to reproduce the results reported in the recent [4FGL Fermi-LAT catalog](https://arxiv.org/abs/1902.10045) reasonably well.

```python
from scipy import special as sf
from astropy.coordinates import SkyCoord
from popsynth import (ZPowerCosmoDistribution, SoftFluxSelection,
                      GalacticPlaneSelection)
					  
from popsynth import SFRDistribution, BPLDistribution, PopulationSynth, NormalAuxSampler, AuxiliarySampler, HardFluxSelection

%matplotlib inline

import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import networkx as nx
import numpy as np
import warnings

warnings.simplefilter("ignore")

```

The work mentioned above presents 3 models for the BL Lac luminosity function. Here, we focus on the case of pure density evolution (referred to as PDE in the paper). In this case, the BL Lac population is parametrised as having a broken power law luminosity distribution, with an independent density evolution following a cosmological power law distribution.

We work with a luminosity range of $L_\mathrm{min} = 7\times 10^{43}$ erg $\mathrm{s}^{-1}$ and $L_\mathrm{max} = 10^{52}$ erg $\mathrm{s}^{-1}$ following Ajello et al. 2014. All luminosities are in units of erg $\mathrm{s}^{-1}$. Similarly, the maxium redshift considered in $z=6$.

We start by setting up the broken power law distribution (`BPLDistribution`).

```python
bpl = BPLDistribution()
bpl.alpha = -1.5
bpl.beta = -2.5
bpl.Lmin = 7e43
bpl.Lmax = 1e52
bpl.Lbreak = 1e47

fig, ax = plt.subplots()
L = np.geomspace(bpl.Lmin, bpl.Lmax)
ax.plot(L, bpl.phi(L), color=purple)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"L [erg $\mathrm{s}^{-1}$]")
```

We now move to the redshift distribution. Following the paper, we parametrize this as a negative power law in $z$. for the purpose of this example, we assume that Bl Lac blazars emit with a steady state. This means that we need to set the `is_rate` parameter to `False` when defining the `ZPowerCosmoDistribution` cosmological distribution. What we mean when we do this is that our local number density, `Lambda` is not per unit time. We also want to survey the whole sky, so we integrate over $4\pi$ sr in the value that we pass to the `Lambda`. 

```python
zpow = ZPowerCosmoDistribution(is_rate=False)
zpow.Lambda = 9000 # Gpc^-3 sr 
zpow.delta = -6

fig, ax = plt.subplots()
z = np.linspace(0.01, 6)
ax.plot(z, zpow.dNdV(z), color=purple)
ax.set_yscale("log")
ax.set_xlabel("z")
ax.set_ylabel(r"$\frac{\mathrm{d}N}{\mathrm{d}V}$")
```

Apart from their redshifts and luminosities, BL Lacs also have other properties. As a simple example, we can consider their spectral index, assuming that the gamma-ray emission is well modelled in the energy range of interest (0.1 to 100 GeV) by a simple power law. 

We assume these true values of these indices are normally distributed with mean, $\mu$, and standard deviation, $\tau$. Additionally, we recognise that these are reconstructed quantities in real surveys, with uncertain values. This is reflected in the error, $\sigma$, on the observed values.

```python
index = NormalAuxSampler(name="index")
index.mu = 2.1
index.tau = 0.25
index.sigma = 0.1
```

We know that the Fermi-LAT detector cannot detect all objects in the Universe, and it is necessary to model some kind of selection function. In general, brighter and spectrally harder objects are easier to detect. We take this into acount by selecting on the flux, $F=L/4\pi d_L^2(z)$, where $d_L$ is the luminosity distance in cm.

This selection effect will not really be a hard boundary, although we could approximate it as such. In reality, the probability to detect an object increases as a function of its flux. To capture this effect, we consider a `SoftFluxSelection` as follows.

```python
flux_selector = SoftFluxSelection()
flux_selector.boundary = 4e-12 # erg cm^-2 s^-1
flux_selector.strength = 2

# This is what is happening under the hood
fig, ax = plt.subplots()
F = np.geomspace(1e-15, 1e-8)
ax.plot(F, sf.expit(flux_selector.strength * (np.log10(F) - np.log10(flux_selector.boundary))), color=purple)
ax.set_xscale("log")
ax.set_xlabel("F [erg $\mathrm{cm}^{-2}$ $\mathrm{s}^{-1}$]")
ax.set_ylabel("Detection prob.")
```

Finally, sometimes it is harder to detect objects near the bright Galactic plane. We take this into account by excluding $10^\circ$ either side of the plane in Galactic longitude using the `GalacticPlaneSelector`.

```python
gp = GalacticPlaneSelection()
gp.b_limit = 10
```

Now, lets finally bring all this together to make a simulated population. Here, we defined our luminosity and spatial distributions already, so we can use them directly in `PopulationSynth`, but there is also the `BPLZPowerCosmoPopulation` available as a quick interface.

```python
# Main pop synth
pop_synth = PopulationSynth(spatial_distribution=zpow, luminosity_distribution=bpl)

# Add our selection effects
pop_synth.set_flux_selection(flux_selector)
pop_synth.add_spatial_selector(gp)

# Add our auxiliary param - spectral index
pop_synth.add_observed_quantity(index)
```

Lets run it! The last parameter to set is adding some uncertainty to our observed flux values.

```python
population = pop_synth.draw_survey(flux_sigma=0.1)
```

We can now have a look at the properties of this simulated population, such as the detected and undetected fluxes and distances.

```python
population.display_fluxes(true_color=purple, obs_color=yellow, with_arrows=False, s=5);
```

```python
fig, ax = plt.subplots()
ax.hist(population.distances, color=purple, histtype="step", lw=3, label="All")
ax.hist(population.distances[population.selection], histtype="step", lw=3, 
        color=yellow, label="Detected")
ax.set_xlabel("$z$")
ax.legend()
```

We can also check out the spectral index distribution.

```python
fig, ax = plt.subplots()
ax.hist(population.index, color=purple, histtype="step", lw=3, label="All")
ax.hist(population.index[population.selection], histtype="step", lw=3, 
        color=yellow, label="Detected")
ax.set_xlabel("Spectral index")
ax.legend()
```

Let's see the distribution of objects on the sky in Galactic coordinates:

```python tags=["nbsphinx-thumbnail"]
c_all = SkyCoord(population.ra, population.dec, unit="deg", frame="icrs")
c_sel = SkyCoord(population.ra[population.selection], 
                 population.dec[population.selection], unit="deg", frame="icrs",)

fig, ax = plt.subplots(subplot_kw={"projection": "hammer"})
ax.scatter(c_all.galactic.l.rad-np.pi, c_all.galactic.b.rad, alpha=0.1, 
           color=purple, label="All")
ax.scatter(c_sel.galactic.l.rad-np.pi, c_sel.galactic.b.rad, alpha=0.8, 
           color=yellow, label="Detected")
ax.axhline(0, color="k")
ax.legend()
```

We can now imagine that by changing the input parameters, we can fit our model to the observations in order to have an optimal representation of the true BL Lac blazar population with this parameterizations.

```python

```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Auxiliary Samplers

Along with sampling the spatial and luminosity distributions, auxiliary properties and be sampled that both depend on and/or influence the luminosity as well each other. This allows you to build up arbitrailiy complex dependencies between parameters which can lead to diverse populations.


```python

import networkx as nx
import numpy as np
import matplotlib.pyplot as plt

%matplotlib inline
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)

purple = "#B833FF"
yellow = "#F6EF5B"

import warnings

warnings.simplefilter("ignore")


import popsynth

popsynth.update_logging_level("INFO")
```


## Built in auxiliary samplers

There are several built in auxiliary samplers that allow you to quickly add on auxiliary parameters.


```python
popsynth.list_available_auxiliary_samplers()
```

We can add these on to the populations, but let's have a look at how to use them.

```python
x = popsynth.NormalAuxSampler(name="aux_param", observed=True)

x.mu = 0
x.sigma = 1

# draws the observed values from normal distribution with std equal to tau
x.tau = 1
```


If value of x is observed (generates data), then we can set the width of the normal distribtuion from which the observed values are sampled from the latent values. Otherwise, only the latent values are stored. This applies to any of the built in auxiliary samplers. However, this can all be customized by adding our own:


## Creating a custom auxiliary sampler
Let's create two auxiliary samplers that sample values from normal distributions with some dependency on each other.

First, we specify the main population. This time, we will chose a SFR-like redshift distribution and a Schecter luminosity function


```python
pop_gen = popsynth.populations.SchechterSFRPopulation(
    r0=100,a=0.0157, rise=1.0, decay=1.0, peak=1.0, Lmin=1e50, alpha=2.0
)
```

<!-- #region -->
Suppose we have a property "demo" that we want to sample as well. For this property, we do not observe it directly. We will get to that. This means that our property latent and could influence other parameters but we can not measure it directly. If you are familiar with Bayesian hierarchical models, this concept may be more familiar to you. As an example, this could be the temperature of a star, which influences its spectrum. The spectrum creates an observable, but the tempreature is imply a random latent variable sampled from a distribution. 


We create an ```AuxiliarySampler``` child class, and define the *true_sampler* for the latent values:
<!-- #endregion -->

```python
class DemoSampler(popsynth.AuxiliarySampler):
    _auxiliary_sampler_name = "DemoSampler"

    mu = popsynth.auxiliary_sampler.AuxiliaryParameter(default=2)
    tau = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)

    def __init__(self):

        # pass up to the super class
        super(DemoSampler, self).__init__("demo", observed=False)

    def true_sampler(self, size):

        # sample the latent values for this property

        self._true_values = np.random.normal(self.mu, self.tau, size=size)
```

Now we instantiate it and then assign it our pop_gen object. Then we draw out survey

```python
demo1 = DemoSampler()

pop_gen.add_observed_quantity(demo1)


flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-9

pop_gen.set_flux_selection(flux_selector)

population = pop_gen.draw_survey()



## plot it
options = {"node_color": purple, "node_size": 3000, "width": 0.5}
pos = nx.drawing.nx_agraph.graphviz_layout(population.graph, prog="dot")
nx.draw(population.graph, with_labels=True, pos=pos, **options)
```


```python
fig = population.display_fluxes(obs_color=purple, true_color=yellow, s=15)
```

We can see that the population has stored out demo auxiliary property.

```python
all_demo = population.demo

obs_demo = population.demo_obs

selected_demo = population.demo_selected
```

We can also see that our demo sampler is now known which is important when creating populations from YAML files. This registering happens when we add the property ```_auxiliary_sampler_name = "DemoSampler" ``` which must be name of the class!

```python
popsynth.list_available_auxiliary_samplers()
```

## Observed auxiliary properties and dependent parameters

Suppose now we want to simulate a property that is observed by an instrument but depends on latent parameters.

We will create a second demo sampler and tell it what the observational error is as well as how to read from a secondary sampler:

```python
class DemoSampler2(popsynth.AuxiliarySampler):
    _auxiliary_sampler_name = "DemoSampler2"
    mu = popsynth.auxiliary_sampler.AuxiliaryParameter(default=2)
    tau = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)
    sigma = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)

    def __init__(
        self,
    ):

        # this time set observed=True
        super(DemoSampler2, self).__init__("demo2", observed=True, uses_distance=True)

    def true_sampler(self, size):

        # we access the secondary sampler dictionary. In this
        # case "demo". This itself is a sampler with
        # <>.true_values as a parameter
        secondary = self._secondary_samplers['demo']

        # now we sample the demo2 latent values and add on the dependence of "demo"

        tmp = np.random.normal(self.mu, self.tau, size=size)

        # for fun, we can substract the log of the distance as all
        # auxiliary samples know about their distances

        self._true_values = tmp + secondary.true_values - np.log10(1 + self._distance)

    def observation_sampler(self, size):

        # here we define the "observed" values, i.e., the latened values
        # with observational error

        self._obs_values = self._true_values + np.random.normal(
            0, self.sigma, size=size
        )
```

We recreate our base sampler:

```python
pop_gen = popsynth.populations.SchechterSFRPopulation(
    r0=100, a=0.0157, rise=1.0, decay=1.0, peak=1.0, Lmin=1e50, alpha=2.0
)
```


Now, make a new *demo1*, but this time we do not have to attach it to the base sampler. Instead, we will assign it as a secondary sampler to *demo2* and **popsynth** is smart enough to search for it when it draws a survey.

```python
demo1 = DemoSampler()


demo2 = DemoSampler2()

demo2.set_secondary_sampler(demo1)

# attach to the base sampler
pop_gen.add_observed_quantity(demo2)
```


```python
pos = nx.drawing.nx_agraph.graphviz_layout(pop_gen.graph, prog="dot")


fig, ax = plt.subplots()


nx.draw(pop_gen.graph, with_labels=True, pos=pos, ax=ax, **options)
```


```python

flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-8

pop_gen.set_flux_selection(flux_selector)
population = pop_gen.draw_survey(flux_sigma=0.1)
```

```python
fig, ax = plt.subplots()

ax.scatter(population.demo2_selected, population.demo_selected, c=purple, s=40)

ax.scatter(population.demo2, population.demo, c=yellow, s=20)
```


## Derived Luminosity sampler

Sometimes, the luminosity does not come directly from a distribution. Rather, it is computed from other quantities. In these cases, we want to use the **DerivedLumAuxSampler** class.

This allows you to sample auxiliary parameters and compute a luminosity from those.

```python
class DemoSampler3(popsynth.DerivedLumAuxSampler):
    _auxiliary_sampler_name = "DemoSampler3"
    mu = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1)
    tau = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)

    def __init__(self, mu=2, tau=1.0, sigma=1):

        # this time set observed=True
        super(DemoSampler3, self).__init__("demo3", uses_distance=False)

    def true_sampler(self, size):

        # draw a random number
        tmp = np.random.normal(self.mu, self.tau, size=size)

        self._true_values = tmp

    def compute_luminosity(self):

        # compute the luminosity
        secondary = self._secondary_samplers["demo"]

        return (10 ** (self._true_values + 54)) + secondary.true_values
```

```python
pop_gen = popsynth.populations.SchechterSFRPopulation(
    r0=100,a=0.0157, rise=1.0, decay=1.0, peak=1.0, Lmin=1e50, alpha=2.0
)
```


```python
demo1 = DemoSampler()


demo3 = DemoSampler3()

demo3.set_secondary_sampler(demo1)

# attach to the base sampler
pop_gen.add_observed_quantity(demo3)


pos = nx.drawing.nx_agraph.graphviz_layout(pop_gen.graph, prog="dot")

fig, ax = plt.subplots()

nx.draw(pop_gen.graph, with_labels=True, pos=pos, **options, ax=ax)
```

```python
flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-5
pop_gen.set_flux_selection(flux_selector)
population = pop_gen.draw_survey(flux_sigma=0.1)
```

```python
population.display_fluxes(obs_color=purple, true_color=yellow, s=15)
```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Selections

Selections on parameters including flux, distance and any auxiliary variables, can be performed in arbitrarily complex way.
We are familiar now with how to add selections onto fluxes and distances, now we will examine in more detail.




## built in selection functions

There are several available selection functions:

```python
import matplotlib.pyplot as plt
import numpy as np

%matplotlib inline

from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import popsynth

popsynth.loud_mode()
popsynth.list_available_selection_functions()
```

We can use these to set selections on parameters. Let's add a dummy parameter that is sampled from a normal distribution:

```python
aux_parameter = popsynth.NormalAuxSampler(name="dummy", observed=False)
aux_parameter.mu = 0
aux_parameter.sigma = 1
```


Now we will use the built in Box selection function. Here, we will assign it to an auxiliary sampler, so we need to tell it to select on the observed value:

```python
box_select = popsynth.BoxSelection(name="aux_selector", use_obs_value=True)
box_select.vmin = 0
box_select.vmax = 0.5
```

We can also add on a selection function for the flux

```python
flux_select = popsynth.HardFluxSelection()
flux_select.boundary = 1e-6
```

Now, we can put it all together and create a survey:

```python
ps = popsynth.SchechterZPowerCosmoPopulation(
    Lambda=50, delta=-2, Lmin=1e52, alpha=1.5, seed=1234
)

aux_parameter.set_selection_probability(box_select)

ps.set_flux_selection(flux_select)

ps.add_auxiliary_sampler(aux_parameter)

pop = ps.draw_survey()
```


```python
fig, ax = plt.subplots()

ax.scatter(
    np.log10(pop.fluxes_observed), pop.dummy, color="purple", alpha=0.7, label="total"
)
ax.scatter(
    np.log10(pop.selected_fluxes_observed),
    pop.dummy_selected,
    color="yellow",
    alpha=0.7,
    label="selected",
)

ax.set(xlabel="log10 fluxes", ylabel="dummy")
ax.legend()
```

<!-- #region -->
## custom selections

we can also create our own custom selection functions.


First, we will look at simply creating a selection. For simplicity, we will look at the Bernoulli selection class built in:
<!-- #endregion -->

```python
class BernoulliSelection(popsynth.SelectionProbability):
    
    # required to register class!
    _selection_name = "BernoulliSelection"

    # define the parameters to be used
    probability = popsynth.SelectionParameter(vmin=0, vmax=1, default=0.5)

    def __init__(self) -> None:

        super(BernoulliSelection, self).__init__(name="Bernoulli")

    def draw(self, size: int) -> None:
        """
        The draw function takes an integer for the size of the 
        samples and sets the private variable _selections which 
        should be an array of boolean values
        
        """
        
        self._selection = stats.bernoulli.rvs(
                self._probability, size=size).astype(bool)  # type: np.ndarray

```

The procedure can become arbitraliy complex. It is important to note that selections will know about several private variables:

```_observed_flux```
```_observed_value```
```_distance```
```_luminosity```


which enables you to use these values in your selection function.

Because of this, several of the build in selections can be used to select on these variables (though some of this is done in the background for you.)


```python
my_box_selection = popsynth.BoxSelection(name="box_flux_selection", use_flux=True)
my_box_selection.vmin = 1E-4
my_box_selection.vmax = 1E-2

```

Setting this as the flux selector will select only the fluxes above and below the limits
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---
# Distributions


The basic required object to create a population synth are a spatial
and (optional if a derived luminosity sampler is create) luminosity
distribution.


```python
%matplotlib inline


import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import networkx as nx
import numpy as np
import warnings

warnings.simplefilter("ignore")
```

`popsynth` comes with several built in distributions included

```python
import popsynth
popsynth.update_logging_level("INFO")


popsynth.list_available_distributions()
```

## Creating a simple population synth

First we create a spatial distribution, in the case, a Spherical distribution with a power law density.


```python
spatial_distribution = popsynth.ZPowerSphericalDistribution()

spatial_distribution.Lambda = 30
spatial_distribution.delta = -2
spatial_distribution.r_max = 10

```

And now we create a powerlaw luminosity distribution

```python
luminosity_distribution = popsynth.ParetoDistribution()

luminosity_distribution.alpha = 1.5
luminosity_distribution.Lmin = 1

```

Combining these together with a random seed, we have a population synthesis object

```python
pop_gen = popsynth.PopulationSynth(luminosity_distribution=luminosity_distribution, 
                                   spatial_distribution = spatial_distribution,
                                   seed=1234
                                  
                                  
                                  )
```

```python
pop_gen.display()
```

```python
population = pop_gen.draw_survey()
```

```python
fig=population.display_obs_fluxes_sphere(background_color="black",size=0.7);
```

## Cosmological Distributions

If we want to create cosmological spatial distributions, we can use
some of those that are built in.

```python
spatial_distribution = popsynth.ZPowerCosmoDistribution()
spatial_distribution.Lambda = 100
spatial_distribution.delta = -2
spatial_distribution.r_max = 10

```

These distributions know about the cosmological Universe and have
their fluxes computed using the luminosity distance rather than linear
distace.

```python
luminosity_distribution = popsynth.SchechterDistribution()

luminosity_distribution.alpha = 1.5
luminosity_distribution.Lmin = 1


```

```python
pop_gen = popsynth.PopulationSynth(luminosity_distribution=luminosity_distribution, 
                                   spatial_distribution = spatial_distribution,
                                   seed=1234
                                  
                                  
                                  )
```

```python
pop_gen.display()
```

```python
population = pop_gen.draw_survey()
```

```python
fig=population.display_obs_fluxes_sphere(cmap="viridis", background_color="black",size=0.7);
```

The cosmological parameters used when simulating are stored in the cosmology object:

```python
popsynth.cosmology.Om
```

```python
popsynth.cosmology.h0
```

```python
popsynth.cosmology.Ode
```


<div class="alert alert-info">

**Note:** The values of Om and h0 can be changed and will change the values of all cosmological calculations

</div>




```python
popsynth.cosmology.Om=0.7
```

```python
popsynth.cosmology.Ode
```

Let's re run the last simulation to see how this changes things

```python
pop_gen.clean()
```

```python
population = pop_gen.draw_survey()
```

```python
fig=population.display_obs_fluxes_sphere(background_color="black",size=0.7);
```

```python

```
---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Contributing 

Contributions to ```popsynth``` are always welcome. They can come in the form of:

## Bug reports

Please use the [Github issue tracking system for any
bugs](https://github.com/grburgess/popsynth/issues), for questions,
and or feature requests.

## Code and more distributions

While it is easy to create custom distributions in your local setup,
if you would like to add them to popsynth directly, go ahead. Please
include tests to ensure that your contributions are compatible with
the code and can be maintained in the long term.

## Documentation

Additions or examples, tutorials, or better explanations are always
welcome. To ensure that the documentation builds with the current
version of the software, I am using
[jupytext](https://jupytext.readthedocs.io/en/latest/) to write the
documentation in Markdown. These are automatically converted to and
executed as jupyter notebooks when changes are pushed to Github.


---
jupyter:
  jupytext:
    formats: ipynb,md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.11.1
  kernelspec:
    display_name: Python 3
    language: python
    name: python3
---

# Quick start

First, lets just run through some examples to see where we are going
by simulating a simple example population which we observe as a
survey. Let's say we are in a giant sphere surrounded by fire flies
that fill the volume homogeneously. Furthermore, the light they emit
follows a Pareto distribution (power law) in luminosity. Of course,
this population can be anything; active galactic nuclei (AGN),
gamma-ray bursts (GRBs), etc. The framework provided in popsynth is
intended to be generic.

```python
%matplotlib inline


import matplotlib.pyplot as plt
from jupyterthemes import jtplot

jtplot.style(context="notebook", fscale=1, grid=False)
purple = "#B833FF"
yellow = "#F6EF5B"

import popsynth

popsynth.update_logging_level("INFO")

import networkx as nx
import numpy as np
import warnings

warnings.simplefilter("ignore")
```

```python nbsphinx="hidden"
class DemoSampler(popsynth.AuxiliarySampler):
    _auxiliary_sampler_name = "DemoSampler"
    mu = popsynth.auxiliary_sampler.AuxiliaryParameter(default=2)
    tau = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)

    def __init__(self):

        super(DemoSampler, self).__init__("demo", observed=False)

    def true_sampler(self, size):

        self._true_values = np.random.normal(self.mu, self.tau, size=size)


class DemoSampler2(popsynth.DerivedLumAuxSampler):
    _auxiliary_sampler_name = "DemoSampler2"
    mu = popsynth.auxiliary_sampler.AuxiliaryParameter(default=2)
    tau = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)
    sigma = popsynth.auxiliary_sampler.AuxiliaryParameter(default=1, vmin=0)

    def __init__(self):

        super(DemoSampler2, self).__init__("demo2")

    def true_sampler(self, size):

        secondary = self._secondary_samplers["demo"]

        self._true_values = (
            (np.random.normal(self.mu, self.tau, size=size))
            + secondary.true_values
            - np.log10(1 + self._distance)
        )

    def observation_sampler(self, size):

        self._obs_values = self._true_values + np.random.normal(
            0, self.sigma, size=size
        )

    def compute_luminosity(self):

        secondary = self._secondary_samplers["demo"]

        return (10 ** (self._true_values + 54)) / secondary.true_values
```



## A spherically homogenous population of fire flies with a pareto luminosity function

**popsynth** comes with several types of populations included, though
you can easily [construct your
own](https://popsynth.readthedocs.io/en/latest/notebooks/custom.html). To
access the built in population synthesizers, one simply instantiates
the population from the **popsynth.populations** module. Here, we will
simulate a survey that has a homogenous spherical spatial distribution
and a pareto distributed luminosity.

```python
homogeneous_pareto_synth = popsynth.populations.ParetoHomogeneousSphericalPopulation(
    Lambda=5, Lmin=1, alpha=2.0  # the density normalization  # lower bound on the LF
)  # index of the LF

print(homogeneous_pareto_synth)
```

```python
homogeneous_pareto_synth.display()

```


If you have [networkx](https://networkx.org) and
[graviz](https://graphviz.readthedocs.io/en/stable/), you can plot a
graph of the connections.

```python
# we can also display a graph of the object


options = {"node_color": purple, "node_size": 3000, "width": 0.5}

pos = nx.drawing.nx_agraph.graphviz_layout(homogeneous_pareto_synth.graph, prog="dot")

nx.draw(homogeneous_pareto_synth.graph, with_labels=True, pos=pos, **options)
```


## Creating a survey

We can now sample from this population with the **draw_survey**
function, but first we need specify how the flux is selected by adding
a flux selection function. Here, we will use a hard selection function
in this example, but you [can make your
own](https://popsynth.readthedocs.io/en/latest/notebooks/selections.html#custom-selections). The
selection function will mark objects with **observed** fluxes below
the selection boundary as "hidden", but we will still have access to
them in our population. 

```python
flux_selector = popsynth.HardFluxSelection()
flux_selector.boundary = 1e-2

homogeneous_pareto_synth.set_flux_selection(flux_selector)
```
And by observed fluxes, we mean those where the latent flux is obscured by observational error, here we sample the observational error from a log normal distribution with $\sigma=1$. In the future, ```popsynth``` will have more options.

```python
population = homogeneous_pareto_synth.draw_survey(flux_sigma=0.1)
```

We now have created a population survey. How did we get here?

* Once the spatial and luminosity functions are specified, we can integrate out to a given distance and compute the number of expected objects.

* A Poisson draw with this mean is made to determine the number of total objects in the survey.

* Next all quantities are sampled (distance, luminosity)

* If needed, the luminosity is converted to a flux with a given observational error

* The selection function (in this case a hard cutoff) is applied

* A population object is created

We could have specified a soft cutoff (an inverse logit) with logarithmic with as well:

```python
homogeneous_pareto_synth.clean()
flux_selector = popsynth.SoftFluxSelection()
flux_selector.boundary = 1e-2
flux_selector.strength = 20


homogeneous_pareto_synth.set_flux_selection(flux_selector)

population = homogeneous_pareto_synth.draw_survey(flux_sigma=0.1)
```

More detail on the [process behind the
simulation](https://popsynth.readthedocs.io/en/latest/notebooks/distributions.html#Core-Concept)
can be found deeper in the documentation

## The Population Object

The population object stores all the information about the sampled
survey. This includes information on the latent parameters, measured
parameters, and distances for both the selected and non-selected
objects.


We can have a look at the flux-distance distribution from the
survey. Here, yellow dots are the *latent* flux value, i.e., without
observational noise, and purple dots are the *measured values for the
*selected* objects. Arrows point from the latent to measured values.

```python
fig = population.display_fluxes(obs_color=purple, true_color=yellow)
```

For fun, we can display the fluxes on in a simulated universe in 3D

```python
fig = population.display_obs_fluxes_sphere(background_color="black")
```

The population object stores a lot of information. For example, an array of selection booleans:

```python
population.selection
```

We can retrieve selected and non-selected distances:

```python
distances = population.selected_distances
```

```python
hidden_distances = population.hidden_distances
```

```python
fig, ax = plt.subplots()

bins = np.linspace(0, 6, 20)


ax.hist(hidden_distances, bins=bins, fc=yellow, ec="k",lw=1)
ax.hist(distances, bins=bins, fc=purple, ec="k",lw=1)
ax.set_xlabel("z")

```

## Saving the population
We can record the results of a population synth to an HDF5 file that
maintains all the information from the run. The true values of the
population parameters are always stored in the truth dictionary:


```python
population.truth
```

```python
population.writeto("saved_pop.h5")
```

```python
reloaded_population = popsynth.Population.from_file("saved_pop.h5")
```

```python
reloaded_population.truth
```

## Creating populations from YAML files

It is sometimes easier to quickly write down population in a YAML file
without having to create all the objects in python. Let's a take a
look at the format:

```yaml

# the seed
seed: 1234

# specifiy the luminosity distribution
# and it's parmeters
luminosity distribution:
    ParetoDistribution:
        Lmin: 1e51
        alpha: 2

# specifiy the flux selection function
# and it's parmeters
flux selection:
    HardFluxSelection:
        boundary: 1e-6

# specifiy the spatial distribution
# and it's parmeters

spatial distribution:
    ZPowerCosmoDistribution:
        Lambda: .5
        delta: -2
        r_max: 5

# specify the distance selection function
# and it's parmeters
distance selection:
    BernoulliSelection:
        probability: 0.5

# a spatial selection if needed
spatial selection:
    # None


# all the auxiliary functions
# these must be known to the
# registry at run time if
# the are custom!

auxiliary samplers:
    stellar_mass
        type: NormalAuxSampler
        observed: False
        mu: 0
        sigma: 1
        selection:
        secondary:
        init variables:

    demo:
        type: DemoSampler
        observed: False
        selection:
            UpperBound:
                boundary: 20

    demo2:
        type: DemoSampler2
        observed: True
        selection:
        secondary: [demo, stellar_mass] # other samplers that this sampler depends on


```

We can load this yaml file into a population synth. We use a saved file to demonstrate:

```python
my_file = popsynth.utils.package_data.get_path_of_data_file("pop.yml")

ps = popsynth.PopulationSynth.from_file(my_file)

print(ps)
```

```python
ps.display()
```

```python
options = {"node_color": purple, "node_size": 3000, "width": 0.5}

pos = nx.drawing.nx_agraph.graphviz_layout(ps.graph, prog="dot")

nx.draw(ps.graph, with_labels=True, pos=pos, **options)
```


<!-- #region -->
We can see that our population was created correctly for us.


Now, this means we can easily pass populations around to our collaborators for testing
<!-- #endregion -->

```python
pop = ps.draw_survey(flux_sigma=0.5)
```

Now, since we can read the population synth from a file, we can also write one we have created with classes to a file:

```python
ps.to_dict()
```

```python
ps.write_to("/tmp/my_pop_synth.yml")
```

but our population synth is also serialized to our population!

```python
pop.pop_synth
```

Therefore we always know exactly how we simulated our data.
.. popsynth documentation master file, created by
   sphinx-quickstart on Mon Aug 19 11:40:04 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to popsynth's documentation!
====================================
.. image:: ../external/logo.png

This framework provides an abstract way to generate survey populations from
arbitrary luminosity functions and redshift distributions. Additionally,
auxiliary quantities can be sampled and stored.

Populations can be saved and restored via an HDF5 files for later
use. Population synthesis routines can be created via classes or
structured YAML files.

Users can construct their own classes for spatial, luminosity,
etc. distributions which can all be connected to arbitrarily complex
selection functions.

.. note:: This is *not* Synth Pop. If you were expecting that… I suggest you check out Depeche Mode. Though, it is possible to combine coding and good music_.


.. image:: ../external/pop.gif

.. _music: https://open.spotify.com/playlist/601WLbJ3Vj91XIugGUJNUe?si=JuXYC9roSxm2PS51aqVaJw



.. toctree::
    :maxdepth: 2
    :hidden:

    notebooks/installation.ipynb
    notebooks/quickstart.ipynb
    notebooks/concept.ipynb
    notebooks/distributions.ipynb
    notebooks/custom.ipynb
    notebooks/selections.ipynb
    notebooks/aux.ipynb
    notebooks/contribute.ipynb
    api/API

.. nbgallery::
   :caption: Examples:

   notebooks/short_grbs.ipynb
   notebooks/bl_lacs.ipynb
   notebooks/stellar_mass.ipynb
   notebooks/milkyway.ipynb
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

.. auto{{ objtype }}:: {{ objname }}
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
popsynth
========

.. toctree::
   :maxdepth: 4

   popsynth
API
===

Here you can find the documentation of all classes and methods:

.. include:: modules.rst
