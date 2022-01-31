# OMUSE #

OMUSE stands for Oceanographic MUltipurpose Software Environment. It is a 
package to conduct numerical experiments in oceanography and other Earth 
sciences.

### What is this repository for? ###

This repository contains the source tree for OMUSE. 

### How do I get set up? ###

Easiest way is to use a pip developer install:

- setup a python environment, e.g. using virtualenv, and activate it.
- optionally (see instructions [below](https://github.com/omuse-geoscience/omuse/blob/master/README.md#amuse-developer-install)) first do a develop install of [AMUSE](http://www.amusecode.org), 
- clone the OMUSE repository: `git clone https://github.com/omuse-geoscience/omuse`,
- go into the source directory `cd omuse` and set the environment variable `DOWNLOAD_CODES`, e.g. `export DOWNLOAD_CODES=latest`.
- now, do `pip install -e .` from the root of the package
- type `python setup.py build_codes --inplace` to build the codes. 
- the file `build.log` will report any errors in the build process.

The community codes of OMUSE can be build manually by going into each directory:

 + src/omuse/community/adcirc
 + src/omuse/community/swan
 + etc

and typing: first `make download` (for some) and then `make`

OMUSE has been tested on OSX and linux machines, with ifort and gfortran 
compilers, on desktop machines and on the Carthesius supercomputer.

In addition to the AMUSE dependencies, OMUSE needs/ can use the following 
packages:

 + matplotlib basemap
 + netCDF and netCDF for fortran and the python bindings
 + GRIB_API

#### AMUSE developer install ####

OMUSE depends on the "amuse-framework" package, and the instructions above will install this automatically from pypi. 
If you want to also have a developer install (from the repository source) for AMUSE you should take care that you install the amuse-framework package:

- clone the AMUSE [repository](https://github.com/amusecode/amuse): `git clone https://github.com/amusecode/amuse`
- go into the amuse-framework package directory: `cd amuse/packages/amuse-framework`
- do the developer install from here: `pip install -e .` 

### Documentation ###

Documentation can be found [here](https://omuse.readthedocs.io). In addition the base  [AMUSE documentation](https://amuse.readthedocs.io) can be consulted.

### Reporting issues ###

Issues can be reported at the OMUSE issue tracker; for framework issues, 
report them at the AMUSE [repository](https://github.com/amusecode/amuse).

### Contribution guidelines ###

Contributions are welcome. Note that most framework development happens at 
the AMUSE [repository](https://github.com/amusecode/amuse) A primer for 
writing code interfaces and other documentation can be found on the amuse 
website (www.amusecode.org).
### DALES python wrapper test suite

* To run the test suite type `nosetests -v` within this directory (using your python environment where OMUSE is installed)
* To test the dales low-level fortran interface to OMUSE, run the `run_dales` binary for an appropriate test case (see `dales-repo/cases`)## examples ##

* bubble-notebook.ipynb: Jupyter notebook containing a DALES convective bubble experiment

The Dales code directory contains a number of examples for its use (src/omuse/community/dales/example).

More examples can be found at: https://github.com/omuse-geoscience/omuse-examples
This package installs the i-emic community code for OMUSE.
This package installs the ERA5 interface for OMUSE.

In order to use the ERA5 interface, UID and API key from the CDS portal is 
needed (see the documentation for [CDSAPI](https://pypi.org/project/cdsapi/) 
for more details).
This package installs the SWAN community code for OMUSE.

SWAN is a third-generation wave model, developed at Delft University of Technology, that computes random, short-crested wind-generated waves in coastal regions and inland waters.

The latest version of SWAN can be found at [here](http://swanmodel.sourceforge.net)
This package installs the base OMUSE framework.
This package installs the qgmodel community code for OMUSE.
This package installs all available OMUSE packages.
DALES interface
===============

This is the OMUSE interface to DALES, the Dutch Atmospheric Large-Eddy Simulation.

Example::
  
  from omuse.community.dales.interface import Dales
  from omuse.units import units

  d = Dales(workdir='dales-workdir', channel_type='sockets', number_of_workers=1)
  
  # set parameters
  d.parameters_DYNAMICS.iadv_mom = 6  # 6th order advection for momentum
  d.parameters_DYNAMICS.iadv_thl = 5  # 5th order advection for scalars, less overshoots than 6th order
  d.parameters_DYNAMICS.iadv_qt = 5
  d.parameters_DYNAMICS.iadv_tke = 5

  # access model variables
  d.fields[:, :, :].U = 0 | units.m / units.s

  # evolve the model
  d.evolve_model(120 | units.s, exactEnd=True)
  

OMUSE models use variables with units, see :ref:`units`.

The parameters are doumented in the DALES documentation: options_. Setting
parameters is possible before the model is time stepped or the
model variables have been accessed, but not after.

Variable grids
--------------

various variables of the DALES model can be accessed from the OMUSE interface. The variables
are grouped in grids, according to their dimensionality.
For example, the U-velocity component in the model component in the model can
be accessed as ``d.fields[i,j,k].U``.

The following grids are defined:

================ ===================================================================
Grid             Description
================ ===================================================================
fields           3D grid
profiles         horizontal averages of the 3D grid (vertical profiles)
nudging_profiles vertical profiles of U, V, THL, QT to nudge the model towards.
forcing_profiles external forcings of U, V, THL, QT - vertical profiles
scalars          surface fluxes and roughness lengths, scalars.
surface_fields   2D surface fields, e.g. liquid water path.
================ ===================================================================

The grids contain the following variables:

================ ===================================================================
Grid             Variables
================ ===================================================================
fields           U, V, W, THL, QT, QL, QL_ice, QR, E12, T, pi, rswd, rswdir, rswdif,
\                rswu, rlwd, rlwu, rswdcs, rswucs, rlwdcs, rlwucs
profiles         U, V, W, THL, QT, QL, QL_ice, QR, E12, T, P, rho
nudging_profiles U, V, THL, QT
forcing_profiles U, V, THL, QT
scalars          wt, wq, z0m, z0h, Ps, QR
surface_fields   LWP, RWP, TWP, ustar, z0m, z0h, tskin, qskin, LE, H, obl
================ ===================================================================


Description of the variables:

============== ====== ===========================================================================
Variable       Unit   Description
============== ====== ===========================================================================
U,V,W          m/s    velocity components in x, y, and z directions
THL            K      liquid water potential temperature
QT             kg/kg  specific total humidity (i.e. vapor and cloud condensate but excluding precipitation)
QL             kg/kg  specific cloud condensate
QL_ice         kg/kg  specific cloud condensate in the form of ice
QR             kg/kg  precipitation
E12            m/s    sqrt(subgrid-scale turbulent kinetic energy)
T              K      temperature
P              Pa     pressure
Ps             Pa     surface pressure
rho            kg/m^3 density
pi             Pa     modified air pressure
rswu,rswd      W/m^2  up-/down-welling short-wave radiative flux
rlwu,rlwd      W/m^2  up-/down-welling long-wave radiative flux
r{s,l}w{u,d}cs W/m^2  up/down-welling long/short-wave radiative flux for clear sky
rswdir         W/m^2  downwelling shortwave direct radiative flux
rswdif         W/m^2  downwelling shortwave diffusive radiative flux
wt             K m/s  surface flux of heat
wq             m/s    surface flux of humidity
z0m, z0h       m      surface roughness length
ustar          m/s    friction velocity
LWP, RWP, TWP  kg/m^2 liquid-, rain-, total water path
tskin          K      skin temperature
qskin          kg/kg  skin humidity
H, LE          W/m^2  sensible and latent heat fluxes   
obl            m      Obukhov length
============== ====== ===========================================================================


The forcing and nudging profiles, and the 3D grids for prognostic variables can be read and written.
Diagnosed quantities, horizontal averages, and water paths can only be read.

Note: the indexing brackets can also be placed on the variable name instead of
on the grid name, e.g. ``d.fields.U[:, :, :]`` vs ``d.fields[:, :, :].U``.
Avoid using the first form, since assigning new values to the grid this way
does not work.



Example scripts
---------------

Some examples of using Dales with OMUSE are bundled included in
omuse/src/omuse/community/dales/example/ .

bubble.py
^^^^^^^^^

A Dales experiment with a warm air bubble. Shows model initialization, setting model parameters, setting initial conditions, evolving the model, and retreiving the state.

async.py
^^^^^^^^

Shows how to use OMUSE's asynchronous function calls, to let several model instances run concurrently.
See also :ref:`asynchronous`.


.. autoclass:: omuse.community.dales.interface.Dales
   :members:

   
      
Links and references
--------------------

* The official DALES git repository_
* DALES manual_
* DALES namelist options_
* Dales model description paper: Formulation of the Dutch Atmospheric Large-Eddy Simulation (DALES) and overview of its applications, T. Heus et al, `Geosci. Model Dev., 3, 415-444, 2010 <https://doi.org/10.5194/gmd-3-415-2010>`_

.. _repository: https://github.com/dalesteam/dales
.. _manual: https://github.com/dalesteam/dales/blob/master/utils/doc/input/dales-manual.pdf
.. _options: https://github.com/dalesteam/dales/blob/master/utils/doc/input/Namoptions.pdf
Installing OMUSE
================

Install required programs and libraries:

* make
* cmake
* python 3
* gcc
* gfortran
* mpi
* netCDF including Fortran bindings
* git
* mercurial
  
Set up and activate a virtual Python environment::
  
    python3 -m venv omuse_env
    source omuse_env/bin/activate

Get the OMUSE source code, install its Python dependencies and set up a development build of OMUSE, where the codes will be built in place::
  
    git clone https://github.com/omuse-geoscience/omuse/
    cd omuse/
    pip install -e .
    export DOWNLOAD_CODES=1

Build codes, select the ones needed::
  
    python setup.py build_code --code-name dales   --inplace
    python setup.py build_code --code-name cdo     --inplace
    python setup.py build_code --code-name qgcm    --inplace
    python setup.py build_code --code-name qgmodel --inplace
    python setup.py build_code --code-name swan    --inplace
    python setup.py build_code --code-name pop     --inplace
   
or try to build all of them::
  
    python setup.py develop_build

Install Jupyter in the virtual environment, and make the virtual environment's Python available as a kernel in Jupyter 

    python -m pip install ipykernel matplotlib
    python -m ipykernel install --user --name=omuse-env


Alternatively, see :ref:`Singularity-section` for instructions for setting up and using a Singularity container with
OMUSE and Jupyter.


Code versions
-------------

For DALES, there are additional options controlling which version of the code is used:
setting `DOWNLOAD_CODES=1` performs a shallow checkout of a single tag, while `DOWNLOAD_CODES="all"`
clones the whole DALES git repository, which is useful for development.
The environment variable `DALES_GIT_TAG` can be used to control which
branch or version tag to check out.
By default the variable points to a version tag in the DALES repository, which is tested to
work with the current OMUSE. 
.. _asynchronous:

Asynchronous calls
==================

Asynchronous function calls and requests is a mechanism in OMUSE to
let multiple models or code instances do work concurrently.

Any method call on an OMUSE object can be made asynchronous, by
appending ``.asynchronous`` to the method's name.
When an asynchronous call is made to a code, the Python script
continues to run while the call is being made.  The Python script can
then perform computations or calls to other codes in OMUSE.
Subsequent calls to the same code instance can be made, but they
will be performed in sequence.


An asynchronous OMUSE function call returns  immediately,
i.e. it does not wait for the worker code to finish the function.
Instead of the normal return value, the asynchronous function returns a request
object, which can be used to access the result of the function call later.
The request can also be added to a request pool, which manages
several concurrent requests.
::

   request1 = d1.evolve_model.asynchronous(target_time, exactEnd=True)
   # ... do something else ...
   print(request1.result()) 

Example
-------

Running two models simultaneously::

    from omuse.community.dales.interface import Dales
    from omuse.units import units
    from amuse.rfi.async_request import AsyncRequestsPool

    # create Dales objects
    d1 = Dales(workdir='dales1', channel_type='sockets', number_of_workers=1, case='bomex')
    d2 = Dales(workdir='dales2', channel_type='sockets', number_of_workers=1, case='bomex')

    # create a pool for managing asynchronous requests
    pool = AsyncRequestsPool()

    # add requests to the two codes to the pool
    request1 = d1.evolve_model.asynchronous(target_time, exactEnd=True)
    pool.add_request(request1)

    request2 = d2.evolve_model.asynchronous(target_time, exactEnd=True)
    pool.add_request(request2)

    # wait for the requests to finish
    pool.waitall()


Asynchronous variable access::
  
    # setting grid data
    # normal synchronous call
    d1.fields[:,:,3:6].U = 1 | units.m / units.s

    # asynchronous call
    d2.fields[:,:,3:6].request.U = 1 | units.m / units.s

    # getting grid data
    # synchronous call, returns an array
    uprofile = d1.fields[5,5,1:10].U

    # asynchronous call, returns a request.
    request = d2.fields[5,5,1:10].request.U
    uprofile = request.result() # retrieve result. Implicit wait for the request to finish

 

Caveats
-------

* Mixing asynchronous and synchronous calls produces correct results,
  but has consequences for the sequencing of the work: Performing a
  normal, synchronous call to a code after one or more asynchronous
  calls, causes an implicit wait for the asynchronous calls to complete.

* Accessing ``result()`` on a request causes a wait for the request to
  complete.

* When making several calls to the same code instance, the first call
  begins immediately. The second call begins only when the code is waited on,
  not automatically when the first call completes.
.. _Singularity-section:

Singularity image
=================

A Singularity recipe is included with OMUSE. When building it, the result is a Singularity image
which is a portable way to test OMUSE and the included models.

The image includes Jupyter for interactive Python notebooks and the following OMUSE models:

 * cdo
 * dales
 * qgcm (hogg2006,ocean_only,hogg2006_20km)
 * qgmodel
 * swan

   

Building the image::

    sudo singularity build omuse.img Singularity 

Create a directory which the container can use at runtime for storing notebooks::

    mkdir run

Launch the container to start Jupyter server inside::

    singularity run --contain -B examples:/opt/notebooks,run:/run/user omuse.img 
    
Then visit the reported localhost:8888 with a browser.

In the singularity command above, the --contain option disables mounting your home directory whereas the -B option specifies paths that will be mounted inside the container, in this case the examples distributed with OMUSE. If you have a folder with notebooks using OMUSE, you can mount this instead of the examples directory to execute them with the singularity container. The container above can also be used to run a regular python script that uses OMUSE functionality by piping the contents to the container python command::

    cat myscript.py | singularity exec --contain omuse.img python3 

Finally, it is also possible to launch a shell inside the container::

    singularity shell --contain -B run:/run/user omuse.img

to execute your python code with all OMUSE dependencies findable. We advise to use the --contain option whenever you
have OMUSE installed on your host system in $HOME/.local .
.. _units:

Units in OMUSE
==============

OMUSE code interfaces use quantities with units.
This has the advantages that the units of any
quantitiy is explicitly specified, and that
automatic conversion can be done between different units.

Unit example::

  from omuse.units import units

  # units are attached to a number with the | operator:
  velocity = 3 | units.m / units.s
  mass = 0.002 | units.kg

  print('velocity', velocity)

  # a quantity can be separated into a number and a unit
  print('number:', velocity.number, 'unit:', velocity.unit)
  print()

  print('mass', mass)

  # get value in a different unit
  print(mass.value_in(units.g), 'g')
  print()

  # arithmetic on quantities:
  print('momentum:', mass * velocity)

output::
  
  velocity 3 m / s
  number: 3 unit: m / s

  mass 0.002 kg
  2.0 g

  momentum: 0.006 m * kg * s**-1



.. OMUSE documentation master file, created by
   sphinx-quickstart on Tue Jul 30 14:59:01 2019.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to OMUSE's documentation!
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installing
   dales/dales
   asynchronous
   units
   singularity


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
