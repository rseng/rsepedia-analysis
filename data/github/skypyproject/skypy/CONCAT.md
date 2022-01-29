## Description

## Checklist
- [ ] Follow the [Contributor Guidelines](https://github.com/skypyproject/skypy/blob/main/CONTRIBUTING.rst)
- [ ] Write unit tests
- [ ] Write documentation strings
- [ ] Assign someone from your working team to review this pull request
- [ ] Assign someone from the infrastructure team to review this pull request
---
name: New model proposal
about: Description of new SkyPy model
title: ''
labels: enhancement
assignees: ''

---

## Description

## Inputs
 -

## Outputs
 -

## References
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.

ADR 3: Position sampling and patches of the sky
===============================================

Author: Nicolas Tessore  
Date: 5 February 2021  
Status: Accepted


Context
-------

This ADR addresses two related open problems for SkyPy simulations:

- How can we describe the regions in which to sample the positions of e.g.
  galaxies or supernovae?
- How can we break up the sampling over large fractions of the sky into smaller
  chunks?

The second point is particularly relevant when the simulation becomes so large
that it no longer fits into memory, as well as for parallelisation.


Extension to pipelines
----------------------

This ADR proposes to introduce a new top-level keyword `regions` [alt: `sky`,
`geometry`, `patches`] into pipelines that describes the geometry of the
simulation. For example, to include a "rectangular" and a "circular" region:

```yaml
regions:
- !rectangular [ ... ]
- !circular [ ... ]
```

The `regions` list can be traversed by the pipeline runner to create what are
effectively independent parallel simulations. The list items are objects with
the following interface.


Region interface
----------------

The regions need to support two sets of operations:

- Information queries: For example, there should be a `.area` [alt:
  `.solid_angle`] attribute that returns the solid angle of the region.
- Random sampling: There needs to be at least a `random_point()` [alt:
  `random_position()`, `uniform_in()`] function that can uniformly sample a
  random position from within the region.

When the pipeline runner traverses the list of regions, it can keep track of
the current region in a `$region` reference that can be used where necessary.
For example, to sample from a luminosity function with positions:

```yaml
tables:
  galaxies:
    z, M: !schechter_lf
      ...
      sky_area: $region.area
    ra, dec: !random_point [ $region ]
```


Support for HEALPix maps
------------------------

The above proposal is powerful enough to support advanced features such as
regions that are described by HEALPix maps. There may be a `healpix()` function
that generates a list of regions from HEALPix pixels:

```yaml
regions: !healpix
  nside: 8
```

The resulting list would contain `12 * nside**2 = 768` regions corresponding
to the HEALPix pixels of a map with `nside = 8`.

The function is easily extensible. For example, instead of using all HEALPix
pixels, there might be a footprint that describes a specific survey:

```yaml
regions: !healpix
  mask: survey-footprint-n512.fits
```

The `mask` keyword can be combined with the `nside` parameter to change the
resolution of the mask if requested.

If the HEALPix maps become finely resolved, it may be desirable to combine
several pixels into a single region. There may be a `batch` [alt: `combine`]
keyword for this purpose:

```yaml
regions: !healpix
  nside: 8
  batch: 16
```

The resulting list of regions will contain `768/16 = 48` regions. The `batch`
keyword may also take a quantity of solid angle and automatically choose the
number of pixels to combine accordingly.


Map making
----------

This ADR does not address the problem of how maps will be generated from the
list of regions. For example, a very real use case would be to generate
populations of galaxies and simply count the total number in each HEALPix pixel
to generate a density map. This will be addressed in a separate ADR.


Consequences
------------
The existing `Pipeline` class must be extended to support iterating regions. No
existing interfaces are affected.
# ADR 2: Mpc or Mpc/h
February 6, 2020

## Context
We need to decide on a unit convention as to whether units include the factor /h or not (for instance Mpc or Mpc/h as a unit of distance). For further discussion see e.g. 10.1017/pasa.2013.31

## Decision Drivers
- Flexibility: Mpc/h allows results to be easily propagated across the 'unknown' value of h (0.68 or 0.74 or something else).
- Consistency / least surprise: the default for astropy is Mpc

## Considered Options
- Mpc
- Mpc/h

## Decision Outcome
After [discussion](https://github.com/skypyproject/skypy/issues/23) and offline, Mpc has been chosen to ensure the closest integration and least surprise for astropy.# ADR 1: Considering options for the SkyPy `Model`
January 22, 2020

## Context
Within SkyPy all functions used to create a "simulation" will in practice be taking in some values (either parameters or columns from a table) and creating new column(s) in an output table *or* selecting specific rows from an input table.

The inputs and outputs of these functions are clearly defined so a directed acyclic graph (DAG) can be constructed to determine what order the functions should be run in.

To aid in the creation of the tables and the DAG a helper class or decorator should be used so the person writing the function does not have to worry about the implementation details. This class or decorator is what we are currently referring to as the `Model`.

For clarity in the options below we will assume the following example function:
```python
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

## Decision Drivers
- Ease of use: if there is too much boiler plate `Model`s will be annoying to write
- Clarity of implementation: the base `Model` should be easy to read, understand, and debug

## Considered Options

### A base `Model` class
In this implementation all functions must be written inside a class that inherits from the base `Model` class.  A different base class should be used depending on if the function adds a column to a table or selects rows.

The `__init__` method would define all the inputs and outputs and the inherited `__init__` can add this to the DAG.

The `compute` method will contain the custom function.

The `execute` method will call the `compute` method and add the results to the table/mask out rows.

- Ease of use: medium (lots of boiler plate)
- Clarity of implementation: high (Classes are well understood by most developers)

Example:
```python
import BaseModel
import numpy as np

class RedshiftGamma(BaseModel):
    def __init__(self):
        self.inputs = ["shape", "scale", "size"]
        self.outputs = ["redshift"]
        super(RedshiftGamma, self).__init__(self.inputs, self.outputs)
    
    def compute(shape, scale, size):
        """Gamma-distributed redshifts (Smail et al. 1994).

        Sample `size` redshifts from a gamma distribution with the
        given `shape` and `scale` parameters. See `numpy.random.gamma`.
        """

        # redshift distribution
        redshift = np.random.gamma(shape=shape, scale=scale, size=size)

        return redshift
```

### A `Model` decorator
In this implementation all functions must use an `@Model(inputs=[], outputs=[])` decorator. A different decorator should be written for adding columns and selecting rows. The decorator will:

1. Add the `inputs` and `outputs` to the DAG
2. Return a callable function that executes the wrapped function and add the results to the table/mask out rows.

- Ease of use: easy (one line added above a function)
- Clarity of implementation: medium (decorators are functions that return function that return function... This particular implementation will be at least 3 wrappers deep)

Example:
```python
import ModelWrapper
import numpy as np

@ModelWrapper(inputs=["shape", "scale", "size"], outputs=["redshift"])
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

### Use the DAG directly
Packages such as [pyungo](https://pypi.org/project/pyungo/) have APIs for most of the functionality we need here with decorators that define `inputs` and `outputs`.  When the compute graph is called all `inputs` and `outputs` are stored in the returned `results` data structure.  Once computed we can write a function that turns this into the final data table.

Also the actual function wrapping only needs to happen for functions contained in the configuration file preventing any un-needed nodes being added to the graph.

We kind of get masking for free here as the DAG does not care if/when the number of rows changes, we just have to be careful when constructing the final table out of the `results`.

- Ease of use: easy (one line added above a function)
- Clarity of implementation: high (we off load this to an existing package that we don't have to maintain)

Example:
```python
def redshift_gamma(shape, scale, size):
    """Gamma-distributed redshifts (Smail et al. 1994).

    Sample `size` redshifts from a gamma distribution with the
    given `shape` and `scale` parameters. See `numpy.random.gamma`.
    """

    # redshift distribution
    redshift = np.random.gamma(shape=shape, scale=scale, size=size)

    return redshift
```

After reading in the config we can wrap all the functions that we need:
```python
from pyungo import Graph

graph = Graph()
graph.register()(redshift_gamma)
res = graph.calculate(data={'shape': 1, 'scale': 1, 'size': 5})
# res is a dict with `shape`, `scale`, `size`, and `redshift`
# this dict can be turned into a table
```

## Decision Outcome
After [discussion](https://github.com/skypyproject/skypy/pull/38) option 3 has been picked.  This will be easiest for developers to write new functions and write clean unit tests.  Within the example given above `pyungo` was just used as an example, other DAG frameworks exist and picking one should be the topic of a different ADR.