# HOOMD-blue code architecture

This document provides a high level overview of the codebase for developers. Details that are reported here relate to the internal API design. For detail on the user-facing public API, see the
[user documentation][hoomd_documentation].

[hoomd_documentation]: https://hoomd-blue.readthedocs.io/en/latest/

## Hardware and software support

Each minor and major release of HOOMD-blue at a minimum supports:

* x86_64 CPUs released in the four prior years.
* NVIDIA GPUs released in the four prior years.
* Compilers and software dependencies available on the oldest maintained Ubuntu LTS release.
* The two most recent major **CUDA** toolkit versions

## Testing

### Continuous integration

[Github Actions] performs continuous integration testing on HOOMD-blue. GitHub
Actions compiles HOOMD-blue, runs the unit tests, and and reports the
status to GitHub pull requests. A number of parallel builds test a variety of
compiler and build configurations as defined above.

Visit the [workflows] page to find recent builds. The pipeline configuration
files are in [.github/workflows/] and are built from Jinja templates in
[.github/workflows/templates/] using `make_workflows.py` which is automatically
run by `pre-commit`. To make changes to the workflows, edit the templates.

[GitHub Actions]: https://docs.github.com/en/actions
[workflows]: https://github.com/glotzerlab/hoomd-blue/actions
[.github/workflows/]: .github/workflows/
[.github/workflows/templates/]: .github/workflows/templates/

## Build system

HOOMD-blue consists of C++ code, Python code, and the CMake configuration scripts necessary to
compile and assemble a functioning Python module. The CMake configuration copies the `.py` files
into the build directory so that developers can iteratively build and test without needing to make
the `install` target and modify files outside the build directory. For more details on using the
build system, see the [user documentation][hoomd_documentation].

HOOMD-blue's CMake configuration follows the most modern CMake standards possible given the software
support constraint given above. For example, it uses `find_package(... CONFIG)` to find package
config files. It manages the linked libraries and additional include directories with the
appropriate visibility in `target_link_libraries` to pass these dependencies on to external
components. HOOMD-blue itself produces a CMake config file to use with `find_package`.

HOOMD has many optional dependencies (e.g. LLVM) and developers can build with or without components
or features (e.g. HPMC). These are set in CMake `ENABLE_*` and `BUILD_*` variables and passed
into the C++ code as preprocessor definitions. New code must observe these definitions so that
the code compiles correctly (or is excluded as needed) when a given option is set or not set.

## C++

The majority of HOOMD-blue's simulation engine is implemented in C++ with a design that strikes a
balance between performance, readability, and maintenance burden. In general, most classes in HOOMD
operate on the entire system of particles so that they can implement loops over the entire system
efficiently. To the extent possible, each class is responsible for a single isolated task and is
composable with other classes. Where needed classes provide a signal/slot mechanism to escape the
isolation and provide notification to client classes when relevant data changes. For example,
`ParticleData` emits a signal when the system box size changes.

This document provides a high level overview of the design, describing how the elements
interoperate. For full details on these classes, see the documentation in the source code comments.
With few exceptions, the C++ code for a class `ClassName` is in `ClassName.h`, `ClassName.cc`,
`ClassName.cuh`, and/or `ClassName.cu`. These files are in the directory corresponding to the Python
package where they reside.

## Execution model

HOOMD-blue implements all operations on the CPU and GPU. The CPU implementation is vanilla C++, and
the GPU implementation uses HIP to support both AMD and NVIDIA GPUs. The `ExecutionConfiguration`
class selects the device (CPU, GPU, or multiple GPUs) and configures global execution options. Each
operation class that needs to know the device configuration is given a shared pointer to the
`ExecutionConfiguration` - most classes store this in the member variable `m_exec_conf`.

To minimize code duplication and to provide a common interface for both CPU and GPU code paths,
HOOMD defines the CPU implementation of an operation in `ClassName` and the GPU implementation in a
subclass `ClassNameGPU`. The base class defines the data structures, parameters, getter/setter
methods, initialization, and other common tasks. The GPU class overrides key methods to perform the
expensive part of the computation on the GPU. The GPU subclass may use alternate data structures for
performance if needed, but this increases the code maintenance burden.

HOOMD-blue uses MPI for domain decomposition simulations on multiple CPUs or GPUs. The
`MPIConfiguration` class (held by `ExecutionConfiguration`) defines the MPI partition and ranks.
Many classes, such as `ParticleData` provide separate methods to access *local* properties on the
current rank and *global* properties (e.g. `ParticleData::getBox` and `ParticleData::getGlobalBox`).
The `Communicator` class is responsible for communicating and migrating particles and bonds between
neighboring ranks.

## Data model

If you think HOOMD-blue's data model is unnecessarily complex, you are correct. Understand that it
is a product of continual development since the year 2007 and has grown along with improving CUDA
functionality. No developer has the time or inclination to completely refactor the entire codebase
to be consistent with the current features. New code should follow the guidelines documented here.

Base data types:

* `Scalar` - Base floating point data type for particle properties. Configurable to either `double`
  or `float` at compile time.
* `Scalar2`, `Scalar3`, `Scalar4`, `int2`, `int3`, ... - 2,3, and 4-vectors of values. These map to
  the CUDA vector types which are aligned properly to enable efficient vector load instructions on
  the GPU. Use these types to store arrays of vector data. Prefer the `2` and `4` size vectors as
  they require fewer memory transactions to read/write than 3-vectors.
* `vec2<Real>` `vec3<Real>`, `quat<Real>` - Templated vector and quaternion types defined in
  `VectorMath.h`. Use these types and the corresponding methods (e.g. `dot`, `operator+`) to perform
  vector and quaternion math with a clean and readable syntax. Convert from and to the `ScalarN`
  vector types when reading inputs and writing outputs to arrays.

Array data:

* `GPUArray<T>` - Template array data type that stores two copies of the data, one on the CPU and
  one on the GPU. Use `ArrayHandle` to request a pointer to the data, which will copy the most
  recently written data to the requested device when needed. New code should use `ArrayHandle` to
  access existing data structures that use `GPUArray`. New code **should not** define new `GPUArray`
  arrays, use `GlobalArray` or `std::vector` with a managed allocator.
* `GlobalArray<T>` - Template array data type that stores one copy of the data in CUDA's unified
  memory in single process, multi-GPU execution. When using a single GPU per process, falls back on
  `GPUArray`. Use `ArrayHandle` to access data in `GlobalArray`.
* `std::vector<T, hoomd::detail::managed_allocator<T>>` - Store array data in a `std::vector` in
  CUDA's unifed memory. This data type is useful for parameter arrays that are exposed to Python.

When using `GlobalArray` or `std::vector<T, hoomd::detail::managed_allocator<T>>`, call
`cudaMemadvise` to set the appropriate memory hints for the array. Small parameter arrays should be
set to `cudaMemAdviseSetReadMostly`. Larger arrays accessed in portions in single-process multi-GPU
execution should be set to `cudaMemAdviseSetPreferredLocation` appropriately for the different
portions of the array.

System data:

* `ParticleData` - Stores the particle positions, velocities, masses, and other per-particle
  properties.
* `ParticleGroup` - Stores the indices of a subset of particles in the system.
* `BondedGroupData` - Stores bonds, angles, and dihedrals.
* `SystemDefinition` - Combines particle data and all bond data.
* `SnapshotSystemData` - Stores a copy of the global system state. Used for file I/O and user
  initialization/analysis.

## Class overview

Users configure HOOMD simulations by defining lists of **Operations** that act on the system state
and schedule when they occur with `Trigger`. The `System` class manages these lists and executes the
simulation in `System::run`. There are numerous types of operations, each with their own base class:

* `Compute` - Compute properties of the system state without modifying it. May provide results
  directly to the user and/or another operation class (e.g. `PotentialPair` uses `NeighborList`).
  A single compute instance may be used by multiple operations, so it must avoid recomputing results
  when `compute` is called multiple times during a single timestep.
* `Updater` - Change the system state.
* `Integrator` - Move the system state forward in time. There is only one `Integrator` in a
  `System`.
* `Analyzer` (named `Writer` in Python) - Computes properties of the system state without modifying
  it and writes them to an output stream or file.
* `Tuner` - Modify parameters of other operations without changing the system state or the
  correctness of the simulation. For example, `SFCPackTuner` reorders the particles in memory to
  improve performance by reducing cache misses.

### HPMC

The integrator `HPMCIntegrator` defines and stores the core parameters of the simulation, such as
the particle shape. All HPMC specific operations (such as `UpdaterClusters` and `UpdaterBoxMC`) take
a shared pointer to the HPMC integrator to access this information.

### MD

There are two MD integrators: `IntegratorTwoStep` implements normal MD simulations and
`FIREENergyMinimizer` implements energy minimization. In MD, the integrator maintains the list of
user-supplied forces (`ForceCompute`) to apply to particles. Both integrators also maintain a list
of user-supplied integration methos (`IntegrationMethodTwoStep`). Each method instance operates
on a single particle group (`ParticleGroup`) and is solely responsible for integrating the equations
of motion of all particles in that group.

## Templates evaluator

Many operations in HOOMD-blue provide similar functionality with different functional forms. For
example pair potentials with many different V(r) and HPMC integration with many different particle
shape classes. To reduce code duplication while maintaining high performance, HOOMD-blue uses
template evaluator classes combined with a single implementation of the general method. This allows
the method (e.g. pair potential evaluation) to be implemented only twice (once on the GPU and once
on the CPU) while each specific evaluator (e.g. V(r)) is also implemented only once. With the
functional form defined in a template class, the compiler is free to inline the evaluation of that
function into the inner loop generated in each template instantiation.

To add a new functional form to the code, a developer must:

1. Implement the evaluator class.
2. Instantiate the GPU kernel driver with the evaluator.
3. Instantiate the method class with the evaluator and export them to Python.
4. Add the Python Operation class to wrap the C++ implementation.

See existing examples in the codebase (e.g. grep for `EvaluatorPairLJ>`) for details.

## GPU kernel driver functions

Early versions of CUDA could compile only a minimal subset of C++ code. While the modern CUDA
compilers are much improved, there are still occasional cases where including complex C++ code like
`pybind11.h` (or even using some standard library features) causes compile errors. This can occur
even when the use of that code is used only in host code. To work around these cases, all GPU
kernels in HOOMD-blue are called via minimal driver functions. These driver functions are not member
functions of their respective class, and therefore must take a long C-style argument list consisting
of bare pointers to data arrays, array sizes, etc... The ABI for these calls is not strictly C, as
driver functions may be templated on functor classes and/or accept lightwight C++ objects as
parameters (such as `BoxDim`).

## Python

The Python code in HOOMD-blue is mostly written to wrap core functionality
written in C++.  Priority is given to ease of use for users in Python even
at the cost of code complexity (within reason).

**Note**: Most internal functions, classes, and methods are documented for
developers (where they are not feel free to add documentation or request it).
Thus, this section will not go over individual functions or methods in detail
except where necessary.

### Base Data Model

#### Definitions
- added: An object is associated with a `Simulation`.
- attached: An object has its corresponding C++ class instantiated and is
  connected to a `Simulation`.
- unattached: An object's data resides in pure Python and may be or not be
  added.
- detach: The transition from attached to unattached.
- removed: The removal of an unattached object from a `Simulation`.
- sync: The process of making a C++ container match a corresponding Python
  container. See the section on `SyncedList` for a concrete example.
- type parameter: An attribute that must be specified for all types or groups of
  types of a set length.

In Python, many base classes exist to facilitate code reuse. This leads to very
lean user facing classes, but requires more fundamental understanding of the model
to address cases where customization is necessary, or changes are required. The
first aspect of the data model discussed is that of most HOOMD operations or
related classes.

- `_HOOMDGetSetAttrBase`
- `_DependencyRelation`
    - `_HOOMDBaseObject`
        - `Operation`

#### `_DependencyRelation`

`_DependencyRelation` helps define a dependent dependency relationship between
two objects. The class is inherited by `_HOOMDBaseObject` to handle dependent
relationships between operations in Python. The class defines *dependents* and
*dependencies* of an object whose removal from a simulation (detaching) can be
handled by overwriting specific methods defined in `_DependencyRelation` in
`hoomd/operation.py`. See the interface of neighbor lists to pair potentials as
an example of this in `hoomd/md/nlist.py` and `hoomd/md/pair/pair.py`.

#### `_HOOMDGetSetAttrBase`

`_HOOMDGetSetAttrBase` provides hooks for exposing object attributes through two
internal (underscored) attributes `_param_dict` and `_typeparam_dict`.
`_param_dict` is an instance of type `hoomd.data.parameterdicts.ParameterDict`,
and `_typeparam_dict` is a dictionary of attribute names to
`hoomd.data.typeparam.TypeParameter`. This serves as the fundamental way of
specifying object attributes in Python. The class provides a multitude of hooks
to enable custom attribute querying and setting when necessary. See
`hoomd/operation.py` for the source code and the sections on `TypeParameter`'s
and `ParameterDict`'s for more information.

**Note**: This class allows the use of `ParameterDict` and `TypeParameter`
(described below) instances without C++ syncing or attaching. Internal custom
actions use this see `hoomd/custom/custom_action.py` for more information.

#### `_HOOMDBaseObject`

`_HOOMDBaseObject` combines `_HOOMDGetSetAttrBase` and `_DependencyRelation` to
provide dependency handling, validated and processed attribute setting,
pickling, and pybind11 C++ class syncing in an automated and structured way.
Most methods unique to or customized from base classes revolve around allowing
`_param_dict` and `_typeparam_dict` to sync to and from C++ when attached and
not when unattached. See `hoomd/operation.py` for source code. The `_add`,
`_attach`, `_remove`, and `_detach` methods handle the process of adding,
attaching, detaching, and removal.

### Attribute Validation and Defaults

In Python, four fundamental classes exist for value validation, processing, and
syncing when attached: `hoomd.data.parameterdicts.ParameterDict`,
`hoomd.data.parameterdicts.TypeParameterDict`,
`hoomd.data.typeparam.TypeParameter`, and
`hoomd.data.syncedlist.SyncedList`.
These can be/are used by `_HOOMDBaseObject` subclasses as well as others.
In addition, classes provided by `hoomd.data.collection` allow for nested
editing of Python objects while maintaining a correspondence to C++.

#### `ParameterDict`

`ParameterDict` provides a mapping (`dict`) interface from attribute names to
validated and processed values. Each instance of `ParameterDict` has its own
specification defining the validation logic for its keys. `ParameterDict`
contains logic when given a pybind11 produced Python object can sync between C++
and Python. See `ParameterDict.__setitem__` for this logic. Attribute specific
logic can be created using the `_getters` and `_setters` attributes. The logic
requires (outside custom getters and setters that all `ParameterDict` keys be
available as properties of the C++ object using the pybind11 property
implementation
(https://pybind11.readthedocs.io/en/stable/classes.html#instance-and-static-fields).
Properties can be read-only which means they will never be set through
`ParameterDict`, but can be through the C++ class constructor. Attempting to set
such a property after attaching will result in an `MutabiliyError` being thrown.

This class should be used to define all attributes shared with C++ member
variables that are on a per-object basis (i.e. not per type). Examples of
this can be seen in many HOOMD classes.  One good example is
`hoomd.md.update.ReversePerturbationFlow`.

#### `TypeParameterDict`

These classes work together to define validated mappings from types or groups of
types to values. The type validation and processing logic is identical to that
used of `ParameterDict`, but on a per key basis. In addition, these classes
support advanced indexing features compared to standard Python `dict` instances.
The class also supports smart defaults. These features can come with some
complexity, so looking at the code with `hoomd.data.typeconverter`,
`hoomd.data.smart_default`, and `hoomd.data.parameterdicts`, should help.

The class is used with `TypeParameter` to specify quantities such as pair
potential `params` (see `hoomd/md/pair/pair.py`) and HPMC shape specs (see
`hoomd/hpmc/integrate.py`).

#### `TypeParameter`

This class is a wrapper for `TypeParameterDict` to work with `_HOOMDBaseObject`
subclass instances. It provides very little independent logic and can be found
in `hoomd/data/typeparam.py`. This class primarily serves as the source of user
documentation for type parameters. `_HOOMDBaseObject` expects the values of
`_typeparam_dict` keys to be `TypeParameter` instances. See the methods
`hoomd.operation._HOOMDBaseObject._append_typeparm` and
`hoomd.operation._HOOMDBaseObject._extend_typeparm` for more information.

This class automatically handles attaching and providing the C++ class the
necessary information through a setter interface (retrieving data is similar).
The name of the setter/getter is the camel cased version of the name given to
the `TypeParameter` (e.g. if the type parameter is named `shape_def` then the
methods are `getShapeDef` and `setShapeDef`). These setters and getters need to
be exposed by the internal C++ class `obj._cpp_obj` instance.

#### `SyncedList`

`SyncedList` implements an arbitrary length list that is value validated and
synced with C++. List objects do not need to have a C++ direct counterpart, but
`SyncedList` must be provided a transform function from the Python object to the
expected C++ one. `SyncedList` can also handle the added and attached status of
its items automatically. An example of this class in use is in the MD integrator
for forces and methods (see `hoomd/md/integrate.py`).

#### Value Validation

The API for specifying value validation and processing in most cases is fairly
simple. The spec `{"foo": float, "bar": str}` does what would be expected,
`"foo"` must be a float and `"bar"` a string. In addition, both value do not
have a default. The basis for the API is that the container type `{}` for `dict`
and `set` (currently not supported), `[]` for `list`, and `()` for `tuple`
defines the object type (tuples validators are considered fixed size) while the
value(s) interior to them define the type to expect or a callable defining
validation logic. For instance,
`{"foo": [(float, float)], "bar": {"baz": lambda x: x % 5}` is perfectly valid
and validates or processes as expected. For more information see
`hoomd/data/typeconverter.py`.

#### Type Parameter Defaults

The effects of these defaults is found in `hoomd.data.TypeParameter`'s API
documentation; however the implementation can be found in
`hoomd/data/smart_default.py`. The defaults are similar to the value validation
in ability to be nested but has less customization due to terminating in
concrete values.

#### HOOMD Collection Types

In `hoomd.data.collection` classes exist which mimic Python `dict`s, `list`s,
and `tuples` (`set`s could easily be added). These classes keep a reference to
their owning object (e.g. a `ParameterDict` or `TypeParameterDict` instance).
Through this reference and read-modify-write approach, these classes facilitate
nested editing of Python attributes while maintaining a synced status in C++. In
general, developers should not need to worry about this as the use of these
classes is automated through previously described classes.

### Internal Python Actions

HOOMD-blue version 3 allows for much more interaction with Python objects within
its C++ core. One feature is custom actions (see the tutorial
https://hoomd-blue.readthedocs.io/en/latest/tutorial/03-Custom-Actions-In-Python/00-index.html
or API documentation for more introductory information). When using custom
actions internally for HOOMD, the classes `hoomd.custom._InternalAction` and one
of the `hoomd.custom._InternalOperation` subclasses are to be used. They allow
for the same interface of other HOOMD operations while having their logic
written in Python. See the examples in `hoomd/write/table.py` and
`hoomd/hpmc/tune/move_size.py` for more information.

### Pickling

By default all Python subclasses of `hoomd.operation._HOOMDBaseObject` support
pickling. This is to facilitate restartability and reproducibility of
simulations. For understanding what *pickling* and Python's supported magic
methods regarding it are see https://docs.python.org/3/library/pickle.html. In
general we prefer using `__getstate__` and `__setstate__` if possible to make
class's picklable.  For the implementation of the default pickling support for
`hoomd.operation._HOOMDBaseObject` see the class's `__getstate__` method.
*Notice* that we do not implement a generic `__setstate__`. We rely on Python's
default generally which is somewhat equivalent to `self.__dict__ =
self.__getstate__()`. Adding a custom `__setstate__` method is fine if necessary
(see [hoomd/write/table.py](hoomd/write/table.py)).  However, using `__reduce__`
is an appropriate alternative if is significantly reduces code complexity or has
demonstrable advantages; see [hoomd/filter/set\_.py](hoomd/filter/set_.py) for
an example of this approach.  _Note_ that `__reduce__` requires that a function
be able to fully recreate the current state of the object (this means that often
times the constructor will not work). Also, note `_HOOMDBaseObject`'s support a
class attribute `_remove_for_pickling` that allows attributes to be removed
before pickling (such as `_cpp_obj`).

**Testing**

To test the pickling of objects see the helper methods in
[hoomd/confest.py](hoomd/conftest.py), `pickling_check` and
`operation_pickling_check`. All objects that are expected to be picklable and
this is most objects in HOOMD-blue should have a pickling test.

A full test suite for collection-like objects can be found in
[hoomd/confest.py](hoomd/conftest.py). This suite is used by all HOOMD
collection like classes.

**Pybind11 Pickling**

For some simple objects like variants or triggers which have very thin Python
wrappers, supporting pickling using pybind11 (see
https://pybind11.readthedocs.io/en/stable/advanced/classes.html#pickling-support)
is acceptable as well. Care just needs to be made that users are not exposed to
C++ classes that "slip" out of their Python subclasses which can happen if no
reference in Python remains to a unpickled object. See
[hoomd/Trigger.cc](hoomd/Trigger.cc) for examples of using pybind11.

**Supporting Class Changes**

Supporting pickling with semantic versioning leads to the need to add support
for objects pickled in version 3.x to work with 3.y, y > x. If new parameters
are added in version "y", then a `__setstate__` method needs to be added if
`__getstate__` and `__setstate__` is the pickling method used for the object.
This `__setstate__` needs to add a default attribute value if one is not
provided in the `dict` given to `__setstate__`. If `__reduce__` was used for
pickling, if a function other than the constructor is used to reinstantiate the
object then any necessary changes should be made (if the constructor is used
then any new arguments to constructor must have defaults via semantic
versioning and no changes should be needed to support pickling). The removal of
internal attributes should not cause problems as well.

## Directory structure

The top level directories are:

* `CMake` - CMake scripts.
* `example_plugin` - External developers to copy this to start developing an external component.
* `hoomd` - Source code for the `hoomd` Python package. Subdirectories under `hoomd` follow the same
  layout as the final Python package.
* `sphinx-doc` - Sphinx configuration and input files for the user-facing documentation.

## Documentation

## User

The user facing documentation is compiled into a human readable document by Sphinx. The
documentation consists of `.rst` files in the `sphinx-doc` directory and the docstrings of
user-facing Python classes in the implementation (imported by the Sphinx autodoc extension).
HOOMD-blue's Sphinx configuration defines mocked imports so that the documentation may be built from
the source directory without needing to compile the C++ source code. This is greatly beneficial when
building the documentation on readthedocs.

The tutorial portion of the documentation is written in Jupyter notebooks housed in the
[hoomd-examples][hoomd_examples] repository. HOOMD-blue includes these in the generated Sphinx
documentation using [nbsphinx][nbsphinx].

[hoomd_examples]: https://github.com/glotzerlab/hoomd-examples
[nbsphinx]: https://nbsphinx.readthedocs.io/

## Detailed developer documentation

Like the user facing classes, internal Python classes document themselves with docstrings.
Similarly, C++ classes provide developer documentation in Javadoc comments. Browse the developer
documentation by viewing the source directly as HOOMD-blue provides no configuration for C++
documentation generation tools.
# HOOMD-blue

[![Citing HOOMD](https://img.shields.io/badge/cite-hoomd-blue.svg)](https://hoomd-blue.readthedocs.io/en/latest/citing.html)
[![conda-forge](https://img.shields.io/conda/vn/conda-forge/hoomd.svg?style=flat)](https://anaconda.org/conda-forge/hoomd)
[![conda-forge Downloads](https://img.shields.io/conda/dn/conda-forge/hoomd.svg?style=flat)](https://anaconda.org/conda-forge/hoomd)
[![GitHub Actions](https://github.com/glotzerlab/hoomd-blue/actions/workflows/test.yml/badge.svg?branch=trunk-patch)](https://github.com/glotzerlab/hoomd-blue/actions/workflows/test.yml)
[![Read the Docs](https://img.shields.io/readthedocs/hoomd-blue/latest.svg)](https://hoomd-blue.readthedocs.io/en/latest/?badge=latest)
[![Contributors](https://img.shields.io/github/contributors-anon/glotzerlab/hoomd-blue.svg?style=flat)](https://hoomd-blue.readthedocs.io/en/latest/credits.html)
[![License](https://img.shields.io/badge/license-BSD--3--Clause-green.svg)](LICENSE)

**HOOMD-blue** is a Python package that runs simulations of particle systems on CPUs and GPUs. It
performs hard particle Monte Carlo simulations of a variety of shape classes and molecular dynamics
simulations of particles with a range of pair, bond, angle, and other potentials. Many features are
targeted at the soft matter research community, though the code is general and capable of many
types of particle simulations.

## Resources

- [Documentation](https://hoomd-blue.readthedocs.io/):
  Tutorial, full Python API description, and usage information.
- [Citing HOOMD-blue](https://hoomd-blue.readthedocs.io/en/latest/citing.html)
  How to cite the code.
- [Installation guide](INSTALLING.rst):
  Instructions for installing **HOOMD-blue** binaries.
- [Compilation guide](BUILDING.rst):
  Instructions for compiling **HOOMD-blue**.
- [hoomd-users mailing list](https://groups.google.com/d/forum/hoomd-users):
  Send messages to the **HOOMD-blue** user community.
- [HOOMD-blue website](https://glotzerlab.engin.umich.edu/hoomd-blue/):
  Additional information and publications.
- [HOOMD-blue benchmark scripts](https://github.com/glotzerlab/hoomd-benchmarks):
  Scripts to evaluate the performance of HOOMD-blue simulations.

## Related tools

- [freud](https://freud.readthedocs.io/):
  Analyze HOOMD-blue simulation results with the **freud** Python library.
- [signac](https://signac.io/):
  Manage your workflow with **signac**.

## Example scripts

These examples demonstrate some of the Python API.

Hard particle Monte Carlo:
```python
import hoomd

mc = hoomd.hpmc.integrate.ConvexPolyhedron()
mc.shape['octahedron'] = dict(vertices=[
    (-0.5, 0, 0),
    (0.5, 0, 0),
    (0, -0.5, 0),
    (0, 0.5, 0),
    (0, 0, -0.5),
    (0, 0, 0.5),
])

cpu = hoomd.device.CPU()
sim = hoomd.Simulation(device=cpu, seed=20)
sim.operations.integrator = mc
# See HOOMD tutorial for how to construct an initial configuration 'init.gsd'
sim.create_state_from_gsd(filename='init.gsd')

sim.run(1e5)
```

Molecular dynamics:
```python
import hoomd

cell = hoomd.md.nlist.Cell()
lj = hoomd.md.pair.LJ(nlist=cell)
lj.params[('A', 'A')] = dict(epsilon=1, sigma=1)
lj.r_cut[('A', 'A')] = 2.5

integrator = hoomd.md.Integrator(dt=0.005)
integrator.forces.append(lj)
nvt = hoomd.md.methods.NVT(kT=1.5, filter=hoomd.filter.All(), tau=1.0)
integrator.methods.append(nvt)

gpu = hoomd.device.GPU()
sim = hoomd.Simulation(device=gpu)
sim.operations.integrator = integrator
# See HOOMD tutorial for how to construct an initial configuration 'init.gsd'
sim.create_state_from_gsd(filename='init.gsd')
sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.5)

sim.run(1e5)
```

## Change log

[CHANGELOG.rst](CHANGELOG.rst) contains the full change log.

## Contributing to HOOMD-blue

Contributions are welcomed via [pull requests](https://github.com/glotzerlab/hoomd-blue/pulls).
Please report bugs and suggest feature enhancements via the [issue
tracker](https://github.com/glotzerlab/hoomd-blue/issues). See [CONTRIBUTING.rst](CONTRIBUTING.rst)
and [ARCHITECTURE.md](ARCHITECTURE.md) for more information.

## License

**HOOMD-blue** is available under the [3-clause BSD license](LICENSE).
# HOOMD-blue Contributor Agreement

These terms apply to your contribution to the HOOMD-blue Open Source Project ("Project") owned or managed by the Regents of the University of Michigan ("Michigan"), and set out the intellectual property rights you grant to Michigan in the contributed materials. If this contribution is on behalf of a company, the term "you" will also mean the company you identify below. If you agree to be bound by these terms, fill in the information requested below and provide your signature.

1. The term "contribution" means any source code, object code, patch, tool, sample, graphic, specification, manual, documentation, or any other material posted or submitted by you to a project.
2. With respect to any worldwide copyrights, or copyright applications and registrations, in your contribution:
    * you hereby assign to Michigan joint ownership, and to the extent that such assignment is or becomes invalid, ineffective or unenforceable, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free, unrestricted license to exercise all rights under those copyrights. This includes, at Michigan's option, the right to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements;
    * you agree that both Michigan and you can do all things in relation to your contribution as if each of us were the sole owners, and if one of us makes a derivative work of your contribution, the one who makes the derivative work (or has it made) will be the sole owner of that derivative work;
    * you agree that you will not assert any moral rights in your contribution against us, our licensees or transferees;
    * you agree that we may register a copyright in your contribution and exercise all ownership rights associated with it; and
    * you agree that neither of us has any duty to consult with, obtain the consent of, pay or render an accounting to the other for any use or distribution of your contribution.
3. With respect to any patents you own, or that you can license without payment to any third party, you hereby grant to Michigan a perpetual, irrevocable, non-exclusive, worldwide, no-charge, royalty-free license to:
    * make, have made, use, sell, offer to sell, import, and otherwise transfer your contribution in whole or in part, alone or in combination with or included in any product, work or materials arising out of the project to which your contribution was submitted; and
    * at Michigan's option, to sublicense these same rights to third parties through multiple levels of sublicensees or other licensing arrangements.
4. Except as set out above, you keep all right, title, and interest in your contribution. The rights that you grant to Michigan under these terms are effective on the date you first submitted a contribution to Michigan, even if your submission took place before the date you sign these terms. Any contribution Michigan makes available under any license will also be made available under a suitable Free Software Foundation or Open Source Initiative approved license.
5. With respect to your contribution, you represent that:
    * it is an original work and that you can legally grant the rights set out in these terms;
    * it does not to the best of your knowledge violate any third party's copyrights, trademarks, patents, or other intellectual property rights; and
you are authorized to sign this contract on behalf of your company (if identified below).
6. The terms will be governed by the laws of the State of Michigan and applicable U.S. Federal Law. Any choice of law rules will not apply.

**By making contribution, you electronically sign and agree to the terms of the HOOMD-blue Contributor Agreement.**

![by-sa.png](https://licensebuttons.net/l/by-sa/3.0/88x31.png)

Based on the Sun Contributor Agreement - version 1.5.
This document is licensed under a Creative Commons Attribution-Share Alike 3.0 Unported License
http://creativecommons.org/licenses/by-sa/3.0/
<!-- Please confirm that your work is based on the correct branch. -->
<!-- Bug fixes should be based on *trunk-patch*. -->
<!-- Backwards compatible new features should be based on *trunk-minor*. -->
<!-- Incompatible API changes should be based on *trunk-major*. -->

## Description

<!-- Describe your changes in detail. -->

## Motivation and context

<!--- Why is this change required? What problem does it solve? -->

<!-- Replace ??? with the issue number that this pull request resolves. -->
Resolves #???

## How has this been tested?

<!--- Please describe in detail how you tested your changes. -->

<!--- Please build the sphinx documentation and check that any changes to
      documentation display properly. -->

## Change log

<!-- Propose a change log entry. -->
```

```

## Checklist:

- [ ] I have reviewed the [**Contributor Guidelines**](https://github.com/glotzerlab/hoomd-blue/blob/trunk-minor/CONTRIBUTING.md).
- [ ] I agree with the terms of the [**HOOMD-blue Contributor Agreement**](https://github.com/glotzerlab/hoomd-blue/blob/trunk-minor/ContributorAgreement.md).
- [ ] My name is on the list of contributors (`sphinx-doc/credits.rst`) in the pull request source branch.
---
name: Task
about: A task that needs to be completed
title: ''
labels: 'task'
assignees: ''

---

## Description

<!-- What needs to be done? -->

## Motivation and context

<!-- What additional information is helpful to understand this task? -->
---
name: Support
about: Request support
title: ''
labels: ''
assignees: ''

---

## Please request support of the mailing list

Please ask user support questions on the hoomd-users mailing list:
https://groups.google.com/d/forum/hoomd-users
---
name: Release checklist
about: '[for maintainer use]'
title: 'Release v3.x.y'
labels: ''
assignees: 'joaander'

---

Release checklist:

- [ ] Update actions versions.
  - See current actions usage with: `rg --no-filename --hidden uses: | awk '{$1=$1;print}' | sort | uniq`
  - Use global search and replace to update them to the latest tags
- [ ] Check for new or duplicate contributors since the last release:
  `comm -13 <(git log LAST_TAG --format="%aN <%aE>" | sort | uniq) <(git log --format="%aN <%aE>" | sort | uniq)`.
  Add entries to `.mailmap` to remove duplicates.
- [ ] Run *bumpversion*.
- [ ] Update change log.
  - ``git log --format=oneline --first-parent `git log -n 1 --pretty=format:%H -- CHANGELOG.rst`...``
  - [milestone](https://github.com/glotzerlab/hoomd-blue/milestones)
- [ ] Check readthedocs build, especially change log formatting.
  - [Build status](https://readthedocs.org/projects/hoomd-blue/builds/)
  - [Output](https://hoomd-blue.readthedocs.io/en/latest/)
- [ ] Tag and push.
- [ ] Update conda-forge recipe.
- [ ] Update glotzerlab-software.
- [ ] Announce release on mailing list.
---
name: Feature request
about: Suggest a new HOOMD-blue feature
title: ''
labels: 'enhancement'
assignees: ''

---

## Description

<!-- What new capability would you like in HOOMD-blue? -->

## Proposed solution

<!-- How should this capability be implemented? -->
<!-- What might the user API look like? -->

## Additional context

<!-- What additional information is helpful to understand this request? -->

## Developer

<!-- Who should implement the new functionality? We would welcome your contribution! -->
---
name: Bug report
about: Report a problem with HOOMD-blue
title: ''
labels: 'bug'
assignees: ''

---

## Description

<!-- Describe the problem. -->

## Script

```python
# Include a minimal script that reproduces the problem

```

<!-- Attach any input files needed to execute the script. -->

## Output

<!-- What output did you get? -->

## Expected output

<!-- What output did you expect? -->

## Configuration

<!-- What is your system configuration? -->

<!-- Remove items that do not apply. -->

Platform:
- Mac
- Linux
- GPU
- CPU

Installation method:
- Conda package
- Compiled from source
- glotzerlab-software container

<!-- What software versions do you have? -->

Versions:
- Python version: ?.?.?
- **HOOMD-blue** version: ?.?.?

## Developer

<!-- Who should implement the fix? We would welcome your contribution! -->
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Building from source
====================

To build the **HOOMD-blue** Python package from source:

1. `Install prerequisites`_::

   $ <package-manager> install cmake eigen git python numpy pybind11

2. `Obtain the source`_::

   $ git clone --recursive https://github.com/glotzerlab/hoomd-blue

3. `Configure`_::

   $ cmake -B build/hoomd -S hoomd-blue

4. `Build the package`_::

   $ cmake --build build/hoomd

5. `Install the package`_ (optional)::

   $ cmake --install build/hoomd

To build the documentation from source (optional):

1. `Install prerequisites`_::

   $ <package-manager> install sphinx sphinx_rtd_theme nbsphinx ipython

.. note::

   ``nbsphinx`` requires ``pandoc>=1.12.1``, which you may need to install separately.

2. `Build the documentation`_::

   $ sphinx-build -b html hoomd-blue/sphinx-doc build/hoomd-documentation

The sections below provide details on each of these steps.

.. _Install prerequisites:

Install prerequisites
---------------------

**HOOMD-blue** requires a number of tools and libraries to build. The options ``ENABLE_MPI``,
``ENABLE_GPU``, ``ENABLE_TBB``, and ``ENABLE_LLVM`` each require additional libraries when enabled.

.. note::

    This documentation is generic. Replace ``<package-manager>`` with your package or module
    manager. You may need to adjust package names and/or install additional packages, such as
    ``-dev`` packages that provide headers needed to build hoomd.

.. tip::

    Create a `virtual environment`_, one place where you can install dependencies and
    **HOOMD-blue**::

        $ python3 -m venv hoomd-venv

    You will need to activate your environment before configuring **HOOMD-blue**::

        $ source hoomd-venv/bin/activate

.. note::

    Some package managers (such as *pip*) and many clusters are missing some or all of **pybind11**,
    **eigen**, and **cereal**. ``install-prereq-headers.py`` will install these packages into your
    virtual environment::

    $ python3 hoomd-blue/install-prereq-headers.py

**General requirements:**

- C++17 capable compiler (tested with ``gcc`` 7, 8, 9, 10, 11 / ``clang`` 6, 7, 8, 9, 10, 11, 12, 13)
- Python >= 3.6
- NumPy >= 1.7
- pybind11 >= 2.2
- Eigen >= 3.2
- CMake >= 3.9

**For MPI parallel execution** (required when ``ENABLE_MPI=on``):

- MPI (tested with OpenMPI, MVAPICH)
- cereal >= 1.1

**For GPU execution** (required when ``ENABLE_GPU=on``):

- NVIDIA CUDA Toolkit >= 9.0

  *OR*

- AMD ROCm >= 3.5.0 with additional dependencies:

  - HIP [with ``hipcc`` and ``hcc`` as backend]
  - rocFFT
  - rocPRIM
  - rocThrust
  - hipCUB, included for NVIDIA GPU targets, but required as an
    external dependency when building for AMD GPUs
  - roctracer-dev
  - Linux kernel >= 3.5.0

  For **HOOMD-blue** on AMD GPUs, the following limitations currently apply.

   1. Certain kernels trigger an `unknown HSA error <https://github.com/ROCm-Developer-Tools/HIP/issues/1662>`_.
   2. The ``mpcd`` component is disabled on AMD GPUs.
   3. Multi-GPU execution via unified memory is not available.

**For threaded parallelism on the CPU** (required when ``ENABLE_TBB=on``):

- Intel Threading Building Blocks >= 4.3

**For runtime code generation** (required when ``ENABLE_LLVM=on``):

- LLVM >= 10.0, < 13

**To build the documentation:**

- sphinx
- sphinx_rtd_theme
- nbsphinx
- ipython

.. _virtual environment: https://docs.python.org/3/library/venv.html

.. _Obtain the source:

Obtain the source
-----------------

Clone using Git_::

   $ git clone --recursive https://github.com/glotzerlab/hoomd-blue

Release tarballs are also available on the `downloads page`_.

.. seealso::

    See the `git book`_ to learn how to work with Git repositories.

.. warning::

    **HOOMD-blue** uses Git submodules. Clone with the ``--recursive`` to clone the submodules.

    Execute ``git submodule update --init`` to fetch the submodules each time you switch branches
    and the submodules show as modified.

.. _downloads page: https://glotzerlab.engin.umich.edu/Downloads/hoomd
.. _git book: https://git-scm.com/book
.. _Git: https://git-scm.com/

.. _Configure:

Configure
---------

Use CMake_ to configure a **HOOMD-blue** build in the given directory. Pass
``-D<option-name>=<value>`` to ``cmake`` to set options on the command line. When modifying code,
you only need to repeat the build step to update your build - it will automatically reconfigure
as needed.

.. tip::

    Use Ninja_ to perform incremental builds in less time::

        $ cmake -B build/hoomd -S hoomd-blue -GNinja

.. tip::

    Place your build directory in ``/tmp`` or ``/scratch`` for faster builds. CMake_ performs
    out-of-source builds, so the build directory can be anywhere on the filesystem.

.. tip::

    Pass the following options to ``cmake`` to optimize the build for your processor:
    ``-DCMAKE_CXX_FLAGS=-march=native -DCMAKE_C_FLAGS=-march=native``.

.. important::

    When using a virtual environment, activate the environment and set the cmake prefix path
    before running CMake_: ``$ export CMAKE_PREFIX_PATH=<path-to-environment>``.

**HOOMD-blue**'s cmake configuration accepts a number of options.

Options that find libraries and executables only take effect on a clean invocation of CMake. To set
these options, first remove ``CMakeCache.txt`` from the build directory and then run ``cmake`` with
these options on the command line.

- ``PYTHON_EXECUTABLE`` - Specify which ``python`` to build against. Example: ``/usr/bin/python3``.

  - Default: ``python3.X`` detected on ``$PATH``.

- ``CMAKE_CUDA_COMPILER`` - Specify which ``nvcc`` or ``hipcc`` to build with.

  - Default: location of ``nvcc`` detected on ``$PATH``.

- ``MPI_HOME`` (env var) - Specify the location where MPI is installed.

  - Default: location of ``mpicc`` detected on the ``$PATH``.

- ``<package-name>_ROOT`` - Specify the location of a package.

  - Default: Found on the `CMake`_ search path.

Other option changes take effect at any time:

- ``BUILD_HPMC`` - When enabled, build the ``hoomd.hpmc`` module (default: ``on``).
- ``BUILD_MD`` - When enabled, build the ``hoomd.md`` module (default: ``on``).
- ``BUILD_METAL`` - When enabled, build the ``hoomd.metal`` module (default: ``on``).
- ``BUILD_TESTING`` - When enabled, build unit tests (default: ``on``).
- ``CMAKE_BUILD_TYPE`` - Sets the build type (case sensitive) Options:

  - ``Debug`` - Compiles debug information into the library and executables. Enables asserts to
    check for programming mistakes. **HOOMD-blue** will run slow when compiled in ``Debug`` mode,
    but problems are easier to identify.
  - ``RelWithDebInfo`` - Compiles with optimizations and debug symbols.
  - ``Release`` - (default) All compiler optimizations are enabled and asserts are removed.
    Recommended for production builds.

- ``CMAKE_INSTALL_PREFIX`` - Directory to install **HOOMD-blue**. Defaults to the root path of the
  found Python executable.
- ``ENABLE_GPU`` - When enabled, compiled GPU accelerated computations (default: ``off``).
- ``SINGLE_PRECISION`` - Controls precision (default: ``off``).

  - When set to ``on``, all calculations are performed in single precision.
  - When set to ``off``, all calculations are performed in double precision.

- ``ENABLE_HPMC_MIXED_PRECISION`` - Controls mixed precision in the ``hpmc`` component. When on,
  single precision is forced in expensive shape overlap checks.
- ``ENABLE_MPI`` - Enable multi-processor/GPU simulations using MPI.

  - When set to ``on``, multi-processor/multi-GPU simulations are supported.
  - When set to ``off`` (the default), always run in single-processor/single-GPU mode.

- ``ENABLE_MPI_CUDA`` - Enable CUDA-aware MPI library support.

  - Requires a MPI library with CUDA support to be installed.
  - When set to ``on``, **HOOMD-blue** will make use of the capability of the MPI library to
    accelerate CUDA-buffer transfers.
  - When set to ``off``, standard MPI calls will be used.

- ``ENABLE_TBB`` - Enable support for Intel's Threading Building Blocks (TBB).

  - When set to ``on``, **HOOMD-blue** will use TBB to speed up calculations in some classes on
    multiple CPU cores.
- ``PYTHON_SITE_INSTALL_DIR`` - Directory to install ``hoomd`` to relative to
  ``CMAKE_INSTALL_PREFIX``. Defaults to the ``site-packages`` directory used by the found Python
  executable.

These options control CUDA compilation via ``nvcc``:

- ``CUDA_ARCH_LIST`` - A semicolon-separated list of GPU architectures to compile.

.. _CMake: https://cmake.org/
.. _Ninja: https://ninja-build.org/

.. _Build the package:

Build the package
-----------------

The command ``cmake --build build/hoomd`` will build the **HOOMD-blue** Python package in the given
build directory. After the build completes, the build directory will contain a functioning Python
package.

.. _Install the package:

Install the package
-------------------

The command ``cmake --install build/hoomd`` installs the given **HOOMD-blue** build to
``${CMAKE_INSTALL_PREFIX}/${PYTHON_SITE_INSTALL_DIR}``. CMake autodetects these paths, but you can
set them manually in CMake.

.. _Build the documentation:

Build the documentation
-----------------------

Run `Sphinx`_ to build the documentation with the command
``sphinx-build -b html hoomd-blue/sphinx-doc build/hoomd-documentation``. Open the file
:file:`build/hoomd-documentation/index.html` in your web browser to view the documentation.

.. tip::

    When iteratively modifying the documentation, the sphinx options ``-a -n -W -T --keep-going``
    are helpful to produce docs with consistent links in the side panel and to see more useful error
    messages::

        $ sphinx-build -a -n -W -T --keep-going -b html \
            hoomd-blue/sphinx-doc build/hoomd-documentation

.. _Sphinx: https://www.sphinx-doc.org/
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Change Log
==========

v3.x
----

v3.0.0 (2022-03-22)
^^^^^^^^^^^^^^^^^^^

*Overview*

HOOMD-blue v3.0.0 is the first production release with the new API that has been developed and
implemented over more than 2 years. Those still using v2.x will need to make changes to their
scripts to use v3. See the `migrating` page for an overview and individual class and method
documentation for more information. To summarize, the new API is object oriented, allows HOOMD-blue
to work effectively as a Python package, and provides more hooks for Python code to directly
interface with the simulation.

*New features in v3 since v2.9.7:*

* Zero-copy data access through numpy and cupy.
* Triggers determine what timesteps operations execute on.
* User-defined operations, triggers, particle filters, variants, and forces.
* Logging subsystem supports array quantities.
* Implicit depletants for 2D shapes in HPMC.
* Harmonically mapped averaging for MD thermodynamic quantities of crystals.
* TWF and OPP pair potentials.
* Tether bond potential.
* Manifold constraints for MD integration methods (using RATTLE) and active forces.
* Document code architecture in ``ARCHITECTURE.md``.
* Overdamped viscous MD integration method.
* User-defined pair potentials work with HPMC on the GPU.
* Long range tail correction for Lennard-Jones potential.
* Anisotropic Lennard-Jones-like pair potential for polyhedra and ellipsoids.
* Newtownian event chain Monte Carlo for spheres and convex polyhedra.

See the full change log below for all v3 beta releases.

Changes from v3.0.0-beta.14:

*Added*

* ``hoomd.hpmc.tune.BoxMCMoveSize`` - Tune ``BoxMC`` move sizes to meet target acceptance ratios.
* ``hoomd.hpmc.nec.integrate.Sphere`` - Newtonian event chain Monte Carlo for hard spheres.
* ``hoomd.hpmc.nec.integrate.ConvexPolyhedron`` - Newtonian event chain Monte Carlo for hard convex
  polyhedra.
* ``hoomd.hpmc.nec.tune.ChainTime`` - Tune chain times in newtonian event chain Monte Carlo method.

*Changed*

* Improve documentation.
* [breaking] Renamed the ``hoomd.md.bond.Table`` energy parameter from ``V`` to ``U``.
* [breaking] Renamed the ``hoomd.md.pair.Table`` energy parameter from ``V`` to ``U``.
* [breaking] Renamed the ``hoomd.md.angle.Table`` energy parameter from ``V`` to ``U``.
* [breaking] Renamed the ``hoomd.md.dihedral.Table`` energy parameter from ``V`` to ``U``.
* [breaking] Renamed ``hoomd.md.nlist.Nlist`` to ``hoomd.md.nlist.NeighborList``.
* [developer] ``Updater`` and ``Analyzer`` in C++ have a ``m_trigger`` member now.
* [developer] ``_TriggeredOperation`` has been moved to ``TriggeredOperation`` and custom trigger
  setting and getting logic removed.

*Fixed*

* ``FIRE.converged`` may be queried before calling ``Simulation.run``.
* Bug where using ``__iadd__`` to certain attributes would fail with an exception.
* Bug where ``hoomd.md.pair.LJ.additional_energy`` is ``NaN`` when ``tail_correction`` is enabled
  and some pairs have ``r_cut=0``.
* Compile error with CUDA 11.7.
* Compile errors on native ubuntu 20.04 systems.
* Compile errors with ``ENABLE_GPU=on`` and ``clang`` as a host compiler.

*Removed*

* [developers] Removed ``IntegratorData`` class. It is replaced by structs that are defined in the
  integrator classes.
* ``get_ordered_vertices`` from ``hoomd.md.pair.aniso.ALJ``.
* Removed optional coxeter dependency.
* The ``limit`` parameter from ``hoomd.md.methods.NVE``.
* The ``limit`` parameter from ``hoomd.md.methods.rattle.NVE``.
* The ``diameter_shift`` parameter from ``hoomd.md.nlist.NeighborList``.
* The ``max_diameter`` parameter from ``hoomd.md.nlist.NeighborList``.

v3.0.0-beta.14 (2022-02-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

* ``hoomd.hpmc.external.field.Harmonic`` - harmonic potential of particles to specific sites in
  the simulation box and orientations.
* Support ``cereal`` 1.3.1
* Guide on how to model molecular systems.
* ``version.floating_point_precision`` - Floating point width in bits for the particle
  properties and local calculations.
* ``hoomd.md.pair.LJ.tail_correction`` - Option to enable the isotropic integrated long range tail
  correction.
* ``hoomd.md.Integrator.linear_momentum`` - Compute the total system linear momentum. Loggable.
* ``hoomd.md.bond.Table`` - Tabulated bond potential.
* ``hoomd.md.angle.Table`` - Tabulated angle potential.
* ``hoomd.md.dihedral.Table`` - Tabulated dihedral potential.
* ``hoomd.md.improper.Harmonic`` - Compute the harmonic improper potential and forces.
* Tutorial on Organizing and executing simulations.
* C++ and build system overview in ``ARCHITECTURE.md``.
* ``hoomd.hpmc.external.wall`` - Overlap checks between particles and wall surfaces.
* ``hoomd.md.pair.ansio.ALJ`` - an anisotropic Lennard-Jones-like pair potential for polyhedra and
  ellipsoids.
* New optional dependency: ``coxeter``, needed for some ``ALJ`` methods.

*Changed*

* Support variant translational and rotational spring constants in
  ``hoomd.hpmc.external.field.Harmonic``.
* [breaking] Renamed ``hoomd.md.angle.Cosinesq`` to ``hoomd.md.angle.CosineSquared``.
* [breaking] ``hoomd.Box`` no longer has a ``matrix`` property use ``to_matrix`` and
  ``from_matrix``.

*Fixed*

* Compilation errors on FreeBSD.
* ``TypeError`` when instantiating special pair forces.
* Inconsistent state when using the ``walls`` setter of a ``hoomd.md.external.wall.WallPotential``.

*Removed*

* [breaking] Removed ``hoomd.md.pair.SLJ`` potential and wall. Use ``hoomd.md.pair.ExpandedLJ``.
* [breaking] ``hoomd.Box.lattice_vectors`` property no longer exists.

v3.0.0-beta.13 (2022-01-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

* ``md.pair.ExpandedLJ`` - A Lennard-Jones potential where ``r`` is replaced with ``r-delta``.
* Support nested modification of operation parameters.
* ``wall`` - Define wall surfaces in the simulation box.
* ``md.external.wall`` - Pair interactions between particles and wall surfaces.
* ``Communicator.walltime`` - the wall clock time since creating the ``Communicator``.
* ``md.force.Custom`` - user defined forces in Python.

*Changed*

* Call ``update_group_dof`` implicitly in ``set_snapshot``, when changing integrators or integration
  methods, and on steps where ``FilterUpdater`` acts on the system.
* [breaking] ``update_group_dof`` defers counting the degrees of freedom until the next timestep or
  the next call to ``Simulation.run``.
* [breaking] Renamed ``md.bond.FENE`` to ``md.bond.FENEWCA``.
* ``md.bond.FENEWCA`` takes a user provided ``delta`` parameter and ignores the particle diameters.
* [breaking] ``md.pair.DLVO`` takes user provided ``a1`` and ``a2`` parameters and ignores the
  particle diameters.
* Removed invalid linker options when using gcc on Apple systems.
* Removed the ``r_on`` attribute and ``default_r_on`` constructor argument from pair potentials that
  do not use it.
* Building from source requires a C++17 compatible compiler.

*Fixed*

* Compile error with ``Apple clang clang-1300.0.29.30``.
* Incorrect OPLS dihedral forces when compiled with ``Apple clang clang-1300.0.29.30``.

*Deprecated*

* ``md.pair.SLJ`` - Replaced with ``md.pair.ExpandedLJ``.

*Removed*

* Leftover ``state`` logging category.

v3.0.0-beta.12 (2021-12-14)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

* Support simulations with arbitrarily large or small scales (within the limits of the floating
  point representation).

*Changed*

* Report full error details in the exception message.
* Improved documentation.
* [breaking]: ``buffer`` is now a required argument when constructing a neighbor list.
* [breaking]: ``force_tol``, ``angmom_tol``, and ``energy_tol`` are now required arguments to
  ``md.minimize.FIRE``

*Fixed*

* Allow neighbor lists to store more than ``2**32-1`` total neighbors.
* Return expected parameter values instead of ``NaN`` when potential parameters are set to 0.

v3.0.0-beta.11 (2021-11-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- Support Python 3.10.
- Support clang 13.

*Changed*

- [developers] Place all all HOOMD C++ classes in the ``hoomd`` and nested namespaces.
- [developers] Use official pre-commit clang-format repository.

v3.0.0-beta.10 (2021-10-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``md.minimize.FIRE`` - MD integrator that minimizes the system's potential energy.
- Include example AKMA and MD unit conversion factors in the documentation.
- ``BUILD_LLVM`` CMake option  (defaults off) to enable features that require LLVM.
- ``hpmc.pair.user.CPPPotential`` - user-defined pair potentials between particles in HPMC.
- ``hpmc.pair.user.CPPPotentialUnion`` - user-defined site-site pair potentials between shapes
  in HPMC.
- ``hpmc.external.user.CPPExternalPotential`` - user-defined external potentials in HPMC.
- Support user-defined pair potentials in HPMC on the GPU.

*Changed*

- Improved documentation.
- Improved error messages when setting operation parameters.
- Noted some dependencies of dependencies for building documentation.
- [developers] Removed ``m_comm`` from most classes. Use ``m_sysdef->isDomainDecomposed()`` instead.
- Add support for LLVM 12
- ``ENABLE_LLVM=on`` requires the clang development libraries.
- [breaking] Renamed the Integrator attribute ``aniso`` to ``integrate_rotational_dof`` and removed
  the ``'auto'`` option. Users must now explicitly choose ``integrate_rotational_dof=True`` to
  integrate the rotational degrees of freedom in the system.

*Fixed*

- Calling ``Operations.__len__`` no longer raises a ``RecursionError``.
- RATTLE integration methods execute on the GPU.
- Include ``EvaluatorPairDLVO.h`` in the installation for plugins.
- Bug in setting zero sized ``ManagedArrays``.
- Kernel launch errors when one process uses different GPU devices.
- Race condition that lead to incorrect simulations with ``md.pair.Table``.
- Bug where some particle filers would have 0 rotational degrees of freedom.

*Removed*

- The ``BUILD_JIT`` CMake option.
- Support for LLVM <= 9.

v3.0.0-beta.9 (2021-09-08)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``Communicator.num_partitions`` - the number of partitions in the communicator.
- ``domain_decomposition`` argument to ``State`` factory methods - set the parameters of the MPI
  domain decomposition
- ``State.domain_decomposition`` - number of domains in the x, y, and z directions in the domain
  decomposition.
- ``State.domain_decomposition_split_fractions`` - the fractional positions of the split planes in
  the domain decomposition.
- ``hoomd.update.FilterUpdater`` - an updater that evaluates the particles associated with a
  `hoomd.filter.ParticleFilter` instance.
- ``hoomd.update.RemoveDrift`` - Remove the average drift from a system restrained on a lattice.
- Developer documentation for HOOMD-blue's Python object data model in ``ARCHITECTURE.md``.
- Autocomplete support for interactive notebooks.
- ``hoomd.md.methods.OverdampedViscous`` - Overdamped integrator with a drag force but no random
  force .
- ``MutabilityError`` exception when setting read-only operation parameters.

*Changed*

- Improved documentation.
- [breaking] Moved ``manifold_constrant`` to separate integration method classes in
  ``hoomd.md.methods.rattle``.
- [breaking] Moved ``trigger`` to first argument position in `hoomd.update.BoxResize`,
  `hoomd.write.DCD`, and `hoomd.write.GSD`.
- [breaking] ``hoomd.data.LocalSnapshot`` particle data API now matches ``Snapshot``. Changes to
  angular momentum, moment of intertia, and rigid body id attributes.
- ``hoomd.write.CustomWriter`` now exposes action through the ``writer`` attribute.
- [breaking] Active force rotational diffusion is managed by
  ``hoomd.md.update.ActiveRotationalDiffusion``.

*Fixed*

- ``TypeParameter`` can set multiple parameters after calling ``hoomd.Simulation.run``.
- ``tune.LoadBalancer`` can be used in a simulation.
- ``hoomd.md.pair.Pair`` ``r_cut`` type parameter can be set to 0.
- MD integration methods can be removed from the integrator's method list.
- Neighborlist exclusions update when the number of bonds change.
- Errors related to equality checks between HOOMD operations.
- The integrator can be removed from a simulation after running.
- ``hoomd.md.constrain.Rigid.create_bodies`` method correctly assigns the body attribute.
- Setting rigid attribute of a MD integrator to ``None`` is allowed.

*Deprecated*

*Removed*

- ``Snapshot.exists`` - use ``snapshot.communicator.rank == 0``
- ``State.snapshot`` - use ``get_snapshot`` and ``set_snapshot``
-   The ``State.box`` property setter - use ``State.set_box``

v3.0.0-beta.8 (2021-08-03)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- Consistent documentation of parameter dimensions and units reference documentation.
- ``md.update.ReversePerturbationFlow`` - implementation of ``mueller_plathe_flow`` from v2.
- ``md.pair.ExpandedMie`` - Mie potential where ``r`` is replaced with ``r - delta``.
- ``md.pair.Table`` - Pair potential evaluated using the given tabulated values.
- ``md.constrain.Distance`` - fix distances between pairs of particles.
- ``hpmc.compute.SDF`` - compute the pressure of convex hard particle systems.
- ``Snapshot.wrap()`` - wrap snapshot particles back into the box.
- Support gcc11.
- ``md.bond.Tether`` - A bond with minimum and maximum lengths.
- ``State.get_snapshot`` and ``State.set_snapshot`` - methods to access the global snapshot.
- ``State.set_box`` set a new simulation box without modifying particle properties.
- ``md.long_range.pppm.make_pppm_coulomb_forces`` - Long range electrostatics evaluated by PPPM.
- ``md.long_range.pppm.Coulomb`` - The reciprocal part of PPPM electrostatics.
- ``md.force.ActiveOnManifold`` - Active forces constrained to manifolds.

*Changed*

- Improved documentation.
- [breaking] Constructor arguments that set a default value per type or pair of types now have
  default in their name (e.g. ``r_cut`` to ``default_r_cut`` for pair potentials and ``a`` to
  ``default_a`` for HPMC integrators).
- [developer] Support git worktree checkouts.
- [breaking] Rename ``nrank`` to ``ranks_per_partition`` in ``Communicator``.
- rowan is now an optional dependency when running unit tests.
- ``Snapshot`` and ``Box`` methods that make in-place modifications return the object.

*Fixed*

- Bug where ``ThermdynamicQuantities.volume`` returned 0 in 2D simulations.
- Update neighbor list exclusions after the number of particles changes.
- Test failures with the CMake option ``BUILD_MD=off``.
- ``write.Table`` can now display MD pressures.

*Deprecated*

- ``State.snapshot`` - use ``get_snapshot`` and ``set_snapshot``.
- The ability to set boxes with the property ``State.box`` - use ``set_box``.

*Removed*

- [breaking] ``Simulation.write_debug_data``.
- [breaking] ``shared_msg_file`` option to ``Device``. ``msg_file`` now has the same behavior as
  ``shared_msg_file``.
- [developers] C++ and Python implementations of ``constraint_ellipsoid``, from ``hoomd.md.update``
  and ``sphere`` and ``oneD`` from ``hoomd.md.constrain``.
- [developers] Doxygen configuration files.


v3.0.0-beta.7 (2021-06-16)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``md.constrain.Rigid`` - Rigid body constraints.
- ``dem_built``, ``hpmc_built``, ``md_built``, and ``mpcd_built`` to ``hoomd.version`` - flags that
  indicate when optional submodules have been built.
- ``GPU.compute_capability`` property.
- [developers] pre-commit enforced style guidelines for the codebase.
- [developers] Validation tests for MD Lennard-Jones simulations.
- [developers] Unit tests for bond, angle, and dihedral potentials.

*Changed*

- Improved documentation on compiling HOOMD.
- Operations raise a ``DataAccessError`` when accessing properties that are not available because
  ``Simulation.run`` has not been called.
- ``TypeConversionError`` is now in the ``hoomd.error`` package.
- ``from_gsd_snapshot`` only accesses the GSD snapshot on MPI rank 0.

*Fixed*

- Some broken references in the documentation.
- Missing documentation for ``md.pair.TWF``.
- Inconsistent documentation in ``md.pair``.
- Correctly identify GPUs by ID in ``GPU.devices``.
- Don't initialize contexts on extra GPUs on MPI ranks.
- Support 2D inputs in ``from_gsd_snapshot``.

*Deprecated*

- ``Snapshot.exists`` - use ``Snapshot.communicator.rank == 0`` instead.

*Removed*

- [developers] C++ implementations of ``rescale_temp`` and ``enforce2d``.
- [developers] Unused methods of ``Integrator``.

v3.0.0-beta.6 (2021-05-17)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``md.pair.LJ0804`` - 8,4 Lennard-Jones pair potential.
- ``md.nlist.Stencil`` - Stencil algorithm to generate neighbor lists.
- ``md.nlist.Tree`` - BVH algorithm to generate neighbor lists.
- ``hoomd.md.Force``, ``hoomd.md.Operation``, and ``hoomd.md.Operations`` objects are now picklable.
- Manifold constraints using RATTLE with ``md.methods.NVE``, ``md.methods.Langevin`` and
  ``md.methods.Brownian``
  - Supporting sphere, ellipsoid, plane, cylinder, gyroid, diamond, and primitive manifolds.
- ``md.compute.HarmonicAveragedThermodynamicQuantities`` - More accurate thermodynamic quantities
  for crystals

*Changed*

- Raise an exception when initializing systems with invalid particle type ids.

*Fixed*

- Setting the operations attribute in ``Simulation`` objects in specific circumstances.
- Misc documentation updates.
- ``'sim' is not defined`` error when using ``md.dihedral`` potentials.

*Removed*

- C++ implemtation of v2 logging infrastructure.

v3.0.0-beta.5 (2021-03-23)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``filter`` parameter to ``update.BoxResize`` - A ``ParticleFilter`` that identifies the particles
  to scale with the box.
- ``Simulation.seed`` - one place to set random number seeds for all operations.
- ``net_force``, ``net_torque``, and ``net_energy`` per-particle arrays in local snapshots.
- Support ``hpmc.update.Clusters`` on the GPU.
- ``hpmc.update.MuVT`` - Gibbs ensemble simulations with HPMC.
- ``md.update.ZeroMomentum`` - Remove linear momentum from the system.
- ``hpmc.compute.FreeVolume`` - Compute free volume available to test particles.
- Custom action tutorials.

*Changed*

- [breaking]  Removed the parameter ``scale_particles`` in ``update.BoxResize``
- [internal] Modified signature of ``data.typeconverter.OnlyTypes``
- Remove use of deprecated numpy APIs.
- Added more details to the migration guide.
- Support timestep values in the range [0,2**64-1].
- [breaking] Removed *seed* argument from ``State.thermalize_particle_momenta``
- [breaking] Removed *seed* argument from ``md.methods.NVT.thermalize_thermostat_dof``
- [breaking] Removed *seed* argument from ``md.methods.NPT.thermalize_thermostat_and_barostat_dof``
- [breaking] Removed *seed* argument from ``md.methods.NPH.thermalize_barostat_dof``
- [breaking] Removed *seed* argument from ``md.methods.Langevin``
- [breaking] Removed *seed* argument from ``md.methods.Brownian``
- [breaking] Removed *seed* argument from ``md.force.Active``
- [breaking] Removed *seed* argument from ``md.pair.DPD``
- [breaking] Removed *seed* argument from ``md.pair.DPDLJ``
- [breaking] Removed *seed* argument from all HPMC integrators.
- [breaking] Removed *seed* argument from ``hpmc.update.Clusters``
- [breaking] Removed *seed* argument from ``hpmc.update.BoxMC``
- [breaking] Removed *seed* argument from ``hpmc.update.QuickCompress``
- Use latest version of getar library.
- Improve documentation.
- Improve performance of ``md.pair.Mie``.
- [breaking] ``hpmc.update.Clusters`` re-implemented with a rejection free, but not ergodic,
  algorithm for anisotropic particles. The new algorithm does not run in parallel over MPI ranks.
- [breaking] HPMC depletion algorithm rewritten.
- [breaking, temporary] HPMC depletant fugacity is now set for type pairs. This change will be
  reverted in a future release.
- Tutorials require fresnel 0.13.
- Support TBB 2021.

*Fixed*

- Install ``ParticleFilter`` header files for external plugins.
- ``md.force.Active`` keeps floating point values set for ``active_force`` and ``active_torque``.
- ``create_state_from_snapshot`` accepts ``gsd.hoomd.Snapshot`` objects without error.
- HOOMD compiles on Apple silicon macOS systems.
- Memory leak in PPPM force compute.
- Segmentation fault that occurred when dumping GSD shapes for spheropolygons and spheropolyhedra
  with 0 vertices.
- Incorrect MD neighbor lists in MPI simulations with more than 1 rank.
- ``md.bond.FENE`` accepts parameters.

*Removed*

- Testing with CUDA 9, GCC 4.8, GCC 5.x, GCC 6.x, clang 5

v3.0.0-beta.4 (2021-02-16)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``hoomd.write.DCD`` - DCD trajectory writer.
- ``hoomd.md.many_body`` - RevCross, SquareDensity, and Tersoff triplet
  potentials.
- ``hoomd.md.methods.Berendsen`` - Berendsen integration method.
- ``hoomd.md.methods.NPH`` - Constant pressure constant enthalpy integration
  method.
- ``hoomd.md.pair.TWF`` - Potential for modeling globular proteins by Pieter
  Rein ten Wolde and Daan Frenkel.
- Custom particle filters in Python via ``hoomd.filter.CustomFilter``.

*Changed*

- Documentation improvements.

*Fixed*

- Correctly determine the maximum ``r_cut`` in simulations with more than one
  pair potential and more than one type.

v3.0.0-beta.3 (2021-01-11)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- ``hoomd.variant.Variant`` objects are picklable.
- ``hoomd.filter.ParticleFilter`` objects are picklable.
- ``hoomd.trigger.Trigger`` objects are picklable.
- ``hoomd.Snapshot.from_gsd_snapshot`` - Convert GSD snapshots to HOOMD.
- ``hoomd.md.pair.aniso.GayBerne`` - Uniaxial ellipsoid pair potential.
- ``hoomd.md.pair.aniso.Dipole`` - Dipole pair potential.
- ``hoomd.md.pair.OPP`` - Oscillating pair potential.

*Changed*

- Improved compilation docs.
- Box equality checking now returns ``NotImplemented`` for non-``hoomd.Box``
  objects.
- ``Simulation.create_state_from_snapshot`` now accepts ``gsd.hoomd.Snapshot``
  objects.
- Attempting to run in a local snapshot context manager will now raise a
  ``RuntimeError``.
- Attempting to set the state to a new snapshot in a local snapshot context
  manager will now raise a ``RuntimeError``.

*Fixed*

- ``hoomd.variant.Power`` objects now have a ``t_ramp`` attribute as documented.
- Enable memory buffers larger than 2-4 GiB.
- Correctly write large image flags to GSD files.
- Support more than 26 default type names.
- Correctly represent fractional degrees of freedom.
- Compute the minimum image in double precision.

v3.0.0-beta.2 (2020-12-15)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Added*

- Support pybind11 2.6.0
- Exclusive creation file mode for ``write.GSD``.
- ``hpmc.update.BoxMC``.
- ``walltime`` and ``final_timestep`` loggable properties in ``Simulation``.
- ``Null`` particle filter.
- Logging tutorial.

*Changed*

- [breaking] Replace ``write.GSD`` argument ``overwrite`` with ``mode``.
- [breaking] Rename ``flags`` to ``categories`` in ``Logger``
- ``hoomd.snapshot.ConfigurationData.dimensions`` is not settable and is
  determined by the snapshot box. If ``box.Lz == 0``, the dimensions are 2
  otherwise 3.
- Building from source requires a C++14 compatible compiler.
- Improved documentation.
- ``hpmc.integrate.FacetedEllipsoid``'s shape specification now has a default
  origin of (0, 0, 0).
- Document loggable quantities in property docstrings.
- Skip GPU tests when no GPU is present.
- ``write.Table`` writes integers with integer formatting.

*Fixed*

- ``Simulation.run`` now ends with a ``KeyboardInterrupt`` exception when
  Jupyter interrupts the kernel.
- Logging the state of specific objects with nested attributes.
- Broken relative RPATHs.
- Add missing documentation for ``version.version``
- Error when removing specific operations from a simulation's operations
  attribute.
- Find CUDA libraries on additional Linux distributions.
- ``hpmc.update.Clusters`` now works with all HPMC integrators.
- ``Simulation.timestep`` reports the correct value when analyzers are called.
- ``Logger`` names quantities with the documented namespace name.

v3.0.0-beta.1 (2020-10-15)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Overview*

v3 has a completely new Python API. See the tutorials, migration guide and new
API documentation learn about it. The API documentation serves as the complete
list of all features currently implemented in v3.0.0-beta.1. Not all features in
v2 have been ported in v3.0.0-beta.1. Future beta releases will add additional
functionality.

*Added*

- Zero-copy data access through numpy (CPU) and cupy (GPU).
- User-defined operations in Python.
- User-defined triggers determine what time steps operations execute on.
- New logging subsystem supports array quantities and binary log files.
- Implicit depletants are now supported by any **hpmc** integrator through
  ``mc.set_fugacity('type', fugacity)``.
- Enable implicit depletants for two-dimensional shapes in **hpmc**.
- ``jit.patch.user()`` and ``jit.patch.user_union()`` now support GPUs via
  NVRTC.
- Add harmonically mapped averaging.
- Add Visual Studio Code workspace

*Changed*

- The ``run`` method has minimal overhead
- All loggable quantities are directly accessible as object properties.
- Operation parameters are always synchronized.
- Operations can be instantiated without a device or MPI communicator.
- Writers write output for ``step+1`` at the bottom of the ``run`` loop.
- HOOMD writes minimal output to stdout/stderr by default.
- *CMake* >=3.9, *cereal*, *eigen*, and *pybind11* are required to compile
  HOOMD.
- Plugins must be updated to build against v3.
- By default, HOOMD installs to the ``site-packages`` directory associated with
  the ``python`` executable given, which may be inside a virtual environment.
- Refactored CMake code.
- ``git submodule update`` no longer runs when during CMake configuration.
- Use ``random123`` library for implicit depletants in **hpmc**.
- HOOMD requires a GPU that supports concurrent managed memory access (Pascal
  or newer).

*Bug fixes*

- Improved accuracy of DLVO potential on the GPU.
- Improved performance of HPMC simulations on the CPU in non-cubic boxes.

*Removed*

- HOOMD-blue no longer parses command line options.
- Type swap moves in ``hpmc.update.muvt()`` are no longer supported
  (``transfer_ratio`` option to ``muvt.set_params()``)
- The option ``implicit=True`` to ``hpmc.integrate.*`` is no longer available
  (use ``set_fugacity``).
- ``static`` parameter in ``dump.gsd``
- ``util.quiet_status`` and ``util.unquiet_status``.
- ``deprecated.analyze.msd``.
- ``deprecated.dump.xml``.
- ``deprecated.dump.pos``.
- ``deprecated.init.read_xml``.
- ``deprecated.init.create_random``.
- ``deprecated.init.create_random_polymers``.
- **hpmc** ``ignore_overlaps`` parameter.
- **hpmc** ``sphere_union::max_members`` parameter.
- **hpmc** ``convex_polyhedron_union``.
- **hpmc** ``setup_pos_writer`` method.
- **hpmc** ``depletant_mode='circumsphere'``.
- **hpmc** ``max_verts`` parameter.
- **hpmc** ``depletant_mode`` parameter.
- **hpmc** ``ntrial`` parameter.
- **hpmc** ``implicit`` boolean parameter.
- ``group`` parameter to ``md.integrate.mode_minimize_fire``
- ``cgcmm.angle.cgcmm``
- ``cgcmm.pair.cgcmm``
- ``COPY_HEADERS`` *CMake* option.
- Many other python modules have been removed or re-implemented with new names.
  See the migration guide and new API documentation for a complete list.
- Support for NVIDIA GPUS with compute capability < 6.0.

v2.x
----

v2.9.7 (2021-08-03)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Support CUDA 11.5. A bug in CUDA 11.4 may result in the error
  ``__global__ function call is not configured`` when running HOOMD.

v2.9.6 (2021-03-16)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Support TBB 2021.

v2.9.5 (2021-03-15)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Support macos-arm64.
* Support TBB 2021.
* Fix memory leak in PPPM.

v2.9.4 (2021-02-05)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Support thrust 1.10
* Support LLVM11
* Fix Python syntax warnings
* Fix compile errors with gcc 10

v2.9.3 (2020-08-05)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Fix a compile error with CUDA 11

v2.9.2 (2020-06-26)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Fix a bug where repeatedly using objects with ``period=None`` would use
  significant amounts of memory.
* Support CUDA 11.
* Reccomend citing the 2020 Computational Materials Science paper
  10.1016/j.commatsci.2019.109363.

v2.9.1 (2020-05-28)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Fixed a minor bug where the variable period timestep would be off by one when
  the timestep got sufficiently large.
* Updated collections API to hide ``DeprecationWarning``.
* Fix scaling of cutoff in Gay-Berne potential to scale the current maximum
  distance based on the orientations of the particles, ensuring ellipsoidal
  energy isocontours.
* Misc documentation fixes.


v2.9.0 (2020-02-03)
^^^^^^^^^^^^^^^^^^^

*New features*

* General

  * Read and write GSD 2.0 files.

    * HOOMD >=2.9 can read and write GSD files created by HOOMD <= 2.8 or GSD
      1.x. HOOMD <= 2.8 cannot read GSD files created by HOOMD >=2.9 or GSD >=
      2.0.
    * OVITO >=3.0.0-dev652 reads GSD 2.0 files.
    * A future release of the ``gsd-vmd`` plugin will read GSD 2.0 files.

* HPMC

  * User-settable parameters in ``jit.patch``.
  * 2D system support in muVT updater.
  * Fix bug in HPMC where overlaps were not checked after adding new particle
    types.

* MD

  * The performance of ``nlist.tree`` has been drastically improved for a
    variety of systems.

v2.8.2 (2019-12-20)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Fix randomization of barostat and thermostat velocities with
  ``randomize_velocities()`` for non-unit temperatures.
* Improve MPCD documentation.
* Fix uninitialized memory in some locations which could have led to
  unreproducible results with HPMC in MPI, in particular with
  ``ALWAYS_USE_MANAGED_MEMORY=ON``.
* Fix calculation of cell widths in HPMC (GPU) and ``nlist.cell()`` with MPI.
* Fix potential memory-management issue in MPI for migrating MPCD particles and
  cell energy.
* Fix bug where exclusions were sometimes ignored when ``charge.pppm()`` is
  the only potential using the neighbor list.
* Fix bug where exclusions were not accounted for properly in the
  ``pppm_energy`` log quantity.
* Fix a bug where MD simulations with MPI start off without a ghost layer,
  leading to crashes or dangerous builds shortly after ``run()``.
* ``hpmc.update.remove_drift`` now communicates particle positions after
  updating them.

v2.8.1 (2019-11-26)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

* Fix a rare divide-by-zero in the ``collide.srd`` thermostat.
* Improve performance of first frame written by ``dump.gsd``.
* Support Python 3.8.
* Fix an error triggering migration of embedded particles for MPCD with MPI +
  GPU configurations.

v2.8.0 (2019-10-30)
^^^^^^^^^^^^^^^^^^^

*New Features*

- MD:

  - ``hoomd.md.dihedral.harmonic`` now accepts phase offsets, ``phi_0``, for CHARMM-style periodic dihedrals.
  - Enable per-type shape information for anisotropic pair potentials that complements the existing pair parameters struct.

- HPMC:

  - Enable the use of an array with adjustable parameters within the user defined pair potential.
  - Add muVT updater for 2D systems.


*Bug fixes*

- Fix missing header in external plugin builds.
- Enable ``couple='none'`` option to ``md.integrate.npt()`` when randomly initializing velocities.
- Documentation improvements.
- Skip gsd shape unit test when required modules are not compiled.
- Fix default particle properties when new particles are added to the system (e.g., via the muVT updater).
- Fix ``charge.pppm()`` execution on multiple GPUs.
- Enable ``with SimulationContext() as c``.
- Fix a bug for ``mpcd.collide.at`` with embedded particles, which may have given incorrect results or simulation crashes.

v2.7.0 (2019-10-01)
^^^^^^^^^^^^^^^^^^^

*New features*

- General:

  - Allow components to use ``Logger`` at the C++ level.
  - Drop support for python 2.7.
  - User-defined log quantities in ``dump.gsd``.
  - Add ``hoomd.dump.gsd.dump_shape`` to save particle shape information in GSD files.

- HPMC:

  - Add ``get_type_shapes`` to ``ellipsoid``.

- MPCD:

  - ``mpcd.stream.slit_pore`` allows for simulations through parallel-plate (lamellar) pores.
  - ``mpcd.integrate`` supports integration of MD (solute) particles with bounce-back rules in MPCD streaming geometries.

*Bug fixes*

- ``hoomd.hdf5.log.query`` works with matrix quantities.
- ``test_group_rigid.py`` is run out of the ``md`` module.
- Fix a bug in ``md.integrate.langevin()`` and ``md.integrate.bd()`` where on the GPU the value of ``gamma`` would be ignored.
- Fix documentation about interoperability between ``md.mode_minimize_fire()`` and MPI.
- Clarify ``dump.gsd`` documentation.
- Improve documentation of ``lattice_field`` and ``frenkel_ladd_energy`` classes.
- Clarify singularity image download documentation.
- Correctly document the functional form of the Buckingham pair potential.
- Correct typos in HPMC example snippets.
- Support compilation in WSL.

v2.6.0 (2019-05-28)
^^^^^^^^^^^^^^^^^^^

*New features*

- General:

  - Enable ``HPMC`` plugins.
  - Fix plug-in builds when ``ENABLE_TBB`` or ``ALWAYS_USE_MANAGED_MEMORY`` CMake parameters are set.
  - Remove support for compute 3.0 GPUs.
  - Report detailed CUDA errors on initialization.
  - Document upcoming feature removals and API changes.

- MD:

  - Exclude neighbors that belong to the same floppy molecule.
  - Add fourier potential.

- HPMC:

  - New shape class: ``hpmc.integrate.faceted_ellipsoid_union()``.
  - Store the *orientable* shape state.

- MPCD:

  - ``mpcd.stream.slit`` allows for simulations in parallel-plate channels. Users can implement other geometries as a plugin.
  - MPCD supports virtual particle filling in bounded geometries through the ``set_filler`` method of ``mpcd.stream`` classes.
  - ``mpcd.stream`` includes an external ``mpcd.force`` acting on the MPCD particles. A block force, a constant force, and a sine force are implemented.

*Bug fixes*

- Fix compile errors with LLVM 8 and ``-DBUILD_JIT=on``.
- Allow simulations with 0 bonds to specify bond potentials.
- Fix a problem where HOOMD could not be imported in ``mpi4py`` jobs.
- Validate snapshot input in ``restore_snapshot``.
- Fix a bug where rigid body energy and pressure deviated on the first time step after ``run()``.
- Fix a bug which could lead to invalid MPI simulations with ``nlist.cell()`` and ``nlist.stencil()``.

*C++ API changes*

- Refactor handling of ``MPI_Comm`` inside library
- Use ``random123`` for random number generation
- CMake version 2.8.10.1 is now a minimum requirement for compiling from source

v2.5.2 (2019-04-30)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

- Support LLVM 9 in ``jit``
- Fix error when importing ``jit`` before ``hpmc``
- HPMC integrators raise errors when ``restore_state=True`` and state information is missing
- Send messages to replaced ``sys.stdout`` and ``sys.stderr`` streams
- Add ``hpmc.update.clusters`` to documentation index
- Fix a bug in the MPCD Gaussian random number generator that could lead to NaN values
- Fix issue where an initially cubic box can become non-cubic with ``integrate.npt()`` and ``randomize_velocities()``
- Fix illegal memory access in NeighborListGPU with ``-DALWAYS_USE_MANAGED_MEMORY=ON`` on single GPUs
- Improve ``pair.table`` performance with multi-GPU execution
- Improve ``charge.pppm`` performance with multi-GPU execution
- Improve rigid body performance with multi-GPU execution
- Display correct cell list statistics with the ``-DALWAYS_USE_MANAGED_MEMORY=ON`` compile option
- Fix a sporadic data corruption / bus error issue when data structures are dynamically resized in simulations that use unified memory (multi-GPU, or with -DALWAYS_USE_MANAGED_MEMORY=ON compile time option)
- Improve ``integrate.nve`` and ``integrate.npt`` performance with multi-GPU execution
- Improve some angular degrees of freedom integrators with multi-GPU execution
- Improve rigid body pressure calculation performance with multi-GPU execution

v2.5.1 (2019-03-14)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

- fix out-of-range memory access in ``hpmc.integrate.convex_polyheron``
- Remove support for clang3.8 and 4.0
- Documentation improvements
- Fix a segfault when using ``SLURM_LOCALID``

v2.5.0 (2019-02-05)
^^^^^^^^^^^^^^^^^^^

*New features*

-  General:

   -  Fix BondedGroupData and CommunicatorGPU compile errors in certain
      build configurations

-  MD:

   -  Generalize ``md.integrate.brownian`` and ``md.integrate.langevin``
      to support anisotropic friction coefficients for rotational
      Brownian motion.
   -  Improve NVLINK performance with rigid bodies
   -  ``randomize_velocities`` now chooses random values for the
      internal integrator thermostat and barostat variables.
   -  ``get_net_force`` returns the net force on a group of particles
      due to a specific force compute

-  HPMC:

   -  Fix a bug where external fields were ignored with the HPMC
      implicit integrator unless a patch potential was also in use.

-  JIT:

   -  Add ``jit.external.user`` to specify user-defined external fields
      in HPMC.
   -  Use ``-DHOOMD_LLVMJIT_BUILD`` now instead of ``-DHOOMD_NOPYTHON``

v2.4.2 (2018-12-20)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Miscellaneous documentation updates
-  Fix compile error with ``with -DALWAYS_USE_MANAGED_MEMORY=ON``
-  Fix MuellerPlatheFlow, cast input parameter to int to avoid C++
   constructor type mismatch
-  Improve startup time with multi-GPU simulations
-  Correctly assign GPUs to MPI processes on Summit when launching with
   more than one GPU per resource set
-  Optimize multi-GPU performance with NVLINK
-  Do not use mapped memory with MPI/GPU anymore
-  Fix some cases where a multi-GPU simulation fails with an alignment
   error
-  Eliminate remaining instance of unsafe ``__shfl``
-  Hide CMake warnings regarding missing CPU math libraries
-  Hide CMake warning regarding missing MPI<->CUDA interoperability
-  Refactor memory management to fix linker errors with some compilers

*C++ API Changes*

-  May break some plug-ins which rely on ``GPUArray`` data type being
   returned from ``ParticleData`` and other classes (replace by
   ``GlobalArray``)

v2.4.1 (2018-11-27)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Install ``WarpTools.cuh`` for use by plugins
-  Fix potential violation of detailed balance with anisotropic
   particles with ``hpmc.update.clusters`` in periodic boundary
   conditions
-  Support llvm 7.0

v2.4.0 (2018-11-07)
^^^^^^^^^^^^^^^^^^^

*New features*

-  General:

   -  Misc documentation updates
   -  Accept ``mpi4py`` communicators in ``context.initialize``.
   -  CUDA 10 support and testing
   -  Sphinx 1.8 support
   -  Flush message output so that ``python -u`` is no longer required
      to obtain output on some batch job systems
   -  Support multi-GPU execution on dense nodes using CUDA managed
      memory. Execute with ``--gpu=0,1,..,n-1`` command line option to
      run on the first n GPUs (Pascal and above).

      -  Node-local acceleration is implemented for a subset of kernels.
         Performance improvements may vary.
      -  Improvements are only expected with NVLINK hardware. Use MPI
         when NVLINK is not available.
      -  Combine the ``--gpu=..`` command line option with mpirun to
         execute on many dense nodes

   -  Bundle ``libgetar`` v0.7.0 and remove ``sqlite3`` dependency
   -  When building with ENABLE_CUDA=on, CUDA 8.0 is now a minimum
      requirement

-  MD:

   -  *no changes*.

-  HPMC:

   -  Add ``convex_spheropolyhedron_union`` shape class.
   -  Correctly count acceptance rate when maximum particle move is is
      zero in ``hpmc.integrate.*``.
   -  Correctly count acceptance rate when maximum box move size is zero
      in ``hpmc.update.boxmc``.
   -  Fix a bug that may have led to overlaps between polygon soups with
      ``hpmc.integrate.polyhedron``.
   -  Improve performance in sphere trees used in
      ``hpmc.integrate.sphere_union``.
   -  Add ``test_overlap`` method to python API

-  API:

   -  Allow external callers of HOOMD to set the MPI communicator
   -  Removed all custom warp reduction and scan operations. These are
      now performed by CUB.
   -  Separate compilation of pair potentials into multiple files.
   -  Removed compute 2.0 workaround implementations. Compute 3.0 is now
      a hard minimum requirement to run HOOMD.
   -  Support and enable compilation for sm70 with CUDA 9 and newer.

-  Deprecated:

   -  HPMC: The implicit depletant mode ``circumsphere`` with
      ``ntrial > 0`` does not support compute 7.0 (Volta) and newer GPUs
      and is now disabled by default. To enable this functionality,
      configure HOOMD with option the ``-DENABLE_HPMC_REINSERT=ON``,
      which will not function properly on compute 7.0 (Volta) and newer
      GPUs.
   -  HPMC: ``convex_polyhedron_union`` is replaced by
      ``convex_spheropolyhedron_union`` (when sweep_radii are 0 for all
      particles)

v2.3.5 (2018-10-07)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Document ``--single-mpi`` command line option.
-  HPMC: Fix a bug where ``hpmc.field.lattice_field`` did not resize 2D
   systems properly in combination with ``update.box_resize``.

v2.3.4 (2018-07-30)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  ``init.read_gsd`` no longer applies the *time_step* override when
   reading the *restart* file
-  HPMC: Add ``hpmc_patch_energy`` and ``hpmc_patch_rcut`` loggable
   quantities to the documentation

v2.3.3 (2018-07-03)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix ``libquickhull.so`` not found regression on Mac OS X

v2.3.2 (2018-06-29)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix a bug where gsd_snapshot would segfault when called without an
   execution context.
-  Compile warning free with gcc8.
-  Fix compile error when TBB include files are in non-system directory.
-  Fix ``libquickhull.so`` not found error on additional platforms.
-  HOOMD-blue is now available on **conda-forge** and the **docker
   hub**.
-  MPCD: Default value for ``kT`` parameter is removed for
   ``mpcd.collide.at``. Scripts that are correctly running are not
   affected by this change.
-  MPCD: ``mpcd`` notifies the user of the appropriate citation.
-  MD: Correct force calculation between dipoles and point charge in
   ``pair.dipole``

*Deprecated*

-  The **anaconda** channel **glotzer** will no longer be updated. Use
   **conda-forge** to upgrade to v2.3.2 and newer versions.

v2.3.1 (2018-05-25)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix doxygen documentation syntax errors
-  Fix libquickhull.so not found error on some platforms
-  HPMC: Fix bug that allowed particles to pas through walls
-  HPMC: Check spheropolyhedra with 0 vertices against walls correctly
-  HPMC: Fix plane wall/spheropolyhedra overlap test
-  HPMC: Restore detailed balance in implicit depletant integrator
-  HPMC: Correctly choose between volume and lnV moves in
   ``hpmc.update.boxmc``
-  HPMC: Fix name of log quantity ``hpmc_clusters_pivot_acceptance``
-  MD: Fix image list for tree neighbor lists in 2d

v2.3.0 (2018-04-25)
^^^^^^^^^^^^^^^^^^^

*New features*

-  General:

   -  Store ``BUILD_*`` CMake variables in the hoomd cmake cache for use
      in external plugins.
   -  ``init.read_gsd`` and ``data.gsd_snapshot`` now accept negative
      frame indices to index from the end of the trajectory.
   -  Faster reinitialization from snapshots when done frequently.
   -  New command line option ``--single-mpi`` allows non-mpi builds of
      hoomd to launch within mpirun (i.e. for use with mpi4py managed
      pools of jobs)
   -  For users of the University of Michigan Flux system: A ``--mode``
      option is no longer required to run hoomd.

-  MD:

   -  Improve performance with ``md.constrain.rigid`` in multi-GPU
      simulations.
   -  New command ``integrator.randomize_velocities()`` sets a particle
      group’s linear and angular velocities to random values consistent
      with a given kinetic temperature.
   -  ``md.force.constant()`` now supports setting the force per
      particle and inside a callback

-  HPMC:

   -  Enabled simulations involving spherical walls and convex
      spheropolyhedral particle shapes.
   -  Support patchy energetic interactions between particles (CPU only)
   -  New command ``hpmc.update.clusters()`` supports geometric cluster
      moves with anisotropic particles and/or depletants and/or patch
      potentials. Supported move types: pivot and line reflection
      (geometric), and AB type swap.

-  JIT:

   -  Add new experimental ``jit`` module that uses LLVM to compile and
      execute user provided C++ code at runtime. (CPU only)
   -  Add ``jit.patch.user``: Compute arbitrary patch energy between
      particles in HPMC (CPU only)
   -  Add ``jit.patch.user_union``: Compute arbitrary patch energy
      between rigid unions of points in HPMC (CPU only)
   -  Patch energies operate with implicit depletant and normal HPMC
      integration modes.
   -  ``jit.patch.user_union`` operates efficiently with additive
      contributions to the cutoff.

-  MPCD:

   -  The ``mpcd`` component adds support for simulating hydrodynamics
      using the multiparticle collision dynamics method.

*Beta feature*

-  Node local parallelism (optional, build with ``ENABLE_TBB=on``):

   -  The Intel TBB library is required to enable this feature.
   -  The command line option ``--nthreads`` limits the number of
      threads HOOMD will use. The default is all CPU cores in the
      system.
   -  Only the following methods in HOOMD will take advantage of
      multiple threads:

      -  ``hpmc.update.clusters()``
      -  HPMC integrators with implicit depletants enabled
      -  ``jit.patch.user_union``

Node local parallelism is still under development. It is not enabled in
builds by default and only a few methods utilize multiple threads. In
future versions, additional methods in HOOMD may support multiple
threads.

To ensure future workflow compatibility as future versions enable
threading in more components, explicitly set –nthreads=1.

*Bug fixes*

-  Fixed a problem with periodic boundary conditions and implicit
   depletants when ``depletant_mode=circumsphere``
-  Fixed a rare segmentation fault with ``hpmc.integrate.*_union()`` and
   ``hpmc.integrate.polyhedron``
-  ``md.force.active`` and ``md.force.dipole`` now record metadata
   properly.
-  Fixed a bug where HPMC restore state did not set ignore flags
   properly.
-  ``hpmc_boxmc_ln_volume_acceptance`` is now available for logging.

*Other changes*

-  Eigen is now provided as a submodule. Plugins that use Eigen headers
   need to update include paths.
-  HOOMD now builds with pybind 2.2. Minor changes to source and cmake
   scripts in plugins may be necessary. See the updated example plugin.
-  HOOMD now builds without compiler warnings on modern compilers (gcc6,
   gcc7, clang5, clang6).
-  HOOMD now uses pybind11 for numpy arrays instead of ``num_util``.
-  HOOMD versions v2.3.x will be the last available on the anaconda
   channel ``glotzer``.

v2.2.5 (2018-04-20)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Pin cuda compatible version in conda package to resolve ``libcu*.so``
   not found errors in conda installations.

v2.2.4 (2018-03-05)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix a rare error in ``md.nlist.tree`` when particles are very close
   to each other.
-  Fix deadlock when ```init.read_getar``` is given different file names
   on different ranks.
-  Sample from the correct uniform distribution of depletants in a
   sphere cap with ``depletant_mode='overlap_regions'`` on the CPU
-  Fix a bug where ternary (or higher order) mixtures of small and large
   particles were not correctly handled with
   ``depletant_mode='overlap_regions'`` on the CPU
-  Improve acceptance rate in depletant simulations with
   ``depletant_mode='overlap_regions'``

v2.2.3 (2018-01-25)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Write default values to gsd frames when non-default values are
   present in frame 0.
-  ``md.wall.force_shifted_lj`` now works.
-  Fix a bug in HPMC where ``run()`` would not start after
   ``restore_state`` unless shape parameters were also set from python.
-  Fix a bug in HPMC Box MC updater where moves were attempted with zero
   weight.
-  ``dump.gsd()`` now writes ``hpmc`` shape state correctly when there
   are multiple particle types.
-  ``hpmc.integrate.polyhedron()`` now produces correct results on the
   GPU.
-  Fix binary compatibility across python minor versions.

v2.2.2 (2017-12-04)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  ``md.dihedral.table.set_from_file`` now works.
-  Fix a critical bug where forces in MPI simulations with rigid bodies
   or anisotropic particles were incorrectly calculated
-  Ensure that ghost particles are updated after load balancing.
-  ``meta.dump_metadata`` no longer reports an error when used with
   ``md.constrain.rigid``
-  Miscellaneous documentation fixes
-  ``dump.gsd`` can now write GSD files with 0 particles in a frame
-  Explicitly report MPI synchronization delays due to load imbalance
   with ``profile=True``
-  Correctly compute net torque of rigid bodies with anisotropic
   constituent particles in MPI execution on multiple ranks
-  Fix ``PotentialPairDPDThermoGPU.h`` for use in external plugins
-  Use correct ghost region with ``constrain.rigid`` in MPI execution on
   multiple ranks
-  ``hpmc.update.muvt()`` now works with
   ``depletant_mode='overlap_regions'``
-  Fix the sampling of configurations with in ``hpmc.update.muvt`` with
   depletants
-  Fix simulation crash after modifying a snapshot and re-initializing
   from it
-  The pressure in simulations with rigid bodies
   (``md.constrain.rigid()``) and MPI on multiple ranks is now computed
   correctly

v2.2.1 (2017-10-04)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Add special pair headers to install target
-  Fix a bug where ``hpmc.integrate.convex_polyhedron``,
   ``hpmc.integrate.convex_spheropolyhedron``,
   ``hpmc.integrate.polyedron``, ``hpmc.integrate.faceted_sphere``,
   ``hpmc.integrate.sphere_union`` and
   ``hpmc.integrate.convex_polyhedron_union`` produced spurious overlaps
   on the GPU

v2.2.0 (2017-09-08)
^^^^^^^^^^^^^^^^^^^

*New features*

-  General:

   -  Add ``hoomd.hdf5.log`` to log quantities in hdf5 format. Matrix
      quantities can be logged.
   -  ``dump.gsd`` can now save internal state to gsd files. Call
      ``dump_state(object)`` to save the state for a particular object.
      The following objects are supported:

      -  HPMC integrators save shape and trial move size state.

   -  Add *dynamic* argument to ``hoomd.dump.gsd`` to specify which
      quantity categories should be written every frame.
   -  HOOMD now inter-operates with other python libraries that set the
      active CUDA device.
   -  Add generic capability for bidirectional ghost communication,
      enabling multi body potentials in MPI simulation.

-  MD:

   -  Added support for a 3 body potential that is harmonic in the local
      density.
   -  ``force.constant`` and ``force.active`` can now apply torques.
   -  ``quiet`` option to ``nlist.tune`` to quiet the output of the
      embedded ``run()`` commands.
   -  Add special pairs as exclusions from neighbor lists.
   -  Add cosine squared angle potential ``md.angle.cosinesq``.
   -  Add ``md.pair.DLVO()`` for evaluation of colloidal dispersion and
      electrostatic forces.
   -  Add Lennard-Jones 12-8 pair potential.
   -  Add Buckingham (exp-6) pair potential.
   -  Add Coulomb 1-4 special_pair potential.
   -  Check that composite body dimensions are consistent with minimum
      image convention and generate an error if they are not.
   -  ``md.integrate.mode.minimize_fire()`` now supports anisotropic
      particles (i.e. composite bodies)
   -  ``md.integrate.mode.minimize_fire()`` now supports flexible
      specification of integration methods
   -  ``md.integrate.npt()/md.integrate.nph()`` now accept a friction
      parameter (gamma) for damping out box fluctuations during
      minimization runs
   -  Add new command ``integrate.mode_standard.reset_methods()`` to
      clear NVT and NPT integrator variables

-  HPMC:

   -  ``hpmc.integrate.sphere_union()`` takes new capacity parameter to
      optimize performance for different shape sizes
   -  ``hpmc.integrate.polyhedron()`` takes new capacity parameter to
      optimize performance for different shape sizes
   -  ``hpmc.integrate.convex_polyhedron`` and
      ``convex_spheropolyhedron`` now support arbitrary numbers of
      vertices, subject only to memory limitations (``max_verts`` is now
      ignored).
   -  HPMC integrators restore state from a gsd file read by
      ``init.read_gsd`` when the option ``restore_state`` is ``True``.
   -  Deterministic HPMC integration on the GPU (optional):
      ``mc.set_params(deterministic=True)``.
   -  New ``hpmc.update.boxmc.ln_volume()`` move allows logarithmic
      volume moves for fast equilibration.
   -  New shape: ``hpmc.integrate.convex_polyhedron_union`` performs
      simulations of unions of convex polyhedra.
   -  ``hpmc.field.callback()`` now enables MC energy evaluation in a
      python function
   -  The option ``depletant_mode='overlap_regions'`` for
      ``hpmc.integrate.*`` allows the selection of a new depletion
      algorithm that restores the diffusivity of dilute colloids in
      dense depletant baths

*Deprecated*

-  HPMC: ``hpmc.integrate.sphere_union()`` no longer needs the
   ``max_members`` parameter.
-  HPMC: ``hpmc.integrate.convex_polyhedron`` and
   ``convex_spheropolyhedron`` no longer needs the ``max_verts``
   parameter.
-  The *static* argument to ``hoomd.dump.gsd`` should no longer be used.
   Use *dynamic* instead.

*Bug fixes*

-  HPMC:

   -  ``hpmc.integrate.sphere_union()`` and
      ``hpmc.integrate.polyhedron()`` missed overlaps.
   -  Fix alignment error when running implicit depletants on GPU with
      ntrial > 0.
   -  HPMC integrators now behave correctly when the user provides
      different RNG seeds on different ranks.
   -  Fix a bug where overlapping configurations were produced with
      ``hpmc.integrate.faceted_sphere()``

-  MD:

   -  ``charge.pppm()`` with ``order=7`` now gives correct results
   -  The PPPM energy for particles excluded as part of rigid bodies now
      correctly takes into account the periodic boundary conditions

-  EAM:

   -  ``metal.pair.eam`` now produces correct results.

*Other changes*

-  Optimized performance of HPMC sphere union overlap check and
   polyhedron shape
-  Improved performance of rigid bodies in MPI simulations
-  Support triclinic boxes with rigid bodies
-  Raise an error when an updater is given a period of 0
-  Revised compilation instructions
-  Misc documentation improvements
-  Fully document ``constrain.rigid``
-  ``-march=native`` is no longer set by default (this is now a
   suggestion in the documentation)
-  Compiler flags now default to CMake defaults
-  ``ENABLE_CUDA`` and ``ENABLE_MPI`` CMake options default OFF. User
   must explicitly choose to enable optional dependencies.
-  HOOMD now builds on powerpc+CUDA platforms (tested on summitdev)
-  Improve performance of GPU PPPM force calculation
-  Use sphere tree to further improve performance of
   ``hpmc.integrate.sphere_union()``

v2.1.9 (2017-08-22)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix a bug where the log quantity ``momentum`` was incorrectly
   reported in MPI simulations.
-  Raise an error when the user provides inconsistent ``charge`` or
   ``diameter`` lists to ``md.constrain.rigid``.
-  Fix a bug where ``pair.compute_energy()`` did not report correct
   results in MPI parallel simulations.
-  Fix a bug where make rigid bodies with anisotropic constituent
   particles did not work on the GPU.
-  Fix hoomd compilation after the rebase in the cub repository.
-  ``deprecated.dump.xml()`` now writes correct results when particles
   have been added or deleted from the simulation.
-  Fix a critical bug where ``charge.pppm()`` calculated invalid forces
   on the GPU

v2.1.8 (2017-07-19)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  ```init.read_getar``` now correctly restores static quantities when
   given a particular frame.
-  Fix bug where many short calls to ``run()`` caused incorrect results
   when using ``md.integrate.langevin``.
-  Fix a bug in the Saru pseudo-random number generator that caused some
   double-precision values to be drawn outside the valid range [0,1) by
   a small amount. Both floats and doubles are now drawn on [0,1).
-  Fix a bug where coefficients for multi-character unicode type names
   failed to process in Python 2.

*Other changes*

-  The Saru generator has been moved into ``hoomd/Saru.h``, and plugins
   depending on Saru or SaruGPU will need to update their includes. The
   ``SaruGPU`` class has been removed. Use ``hoomd::detail::Saru``
   instead for both CPU and GPU plugins.

v2.1.7 (2017-05-11)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix PPM exclusion handling on the CPU
-  Handle ``r_cut`` for special pairs correctly
-  Fix tauP reference in NPH documentation
-  Fixed ``constrain.rigid`` on compute 5.x.
-  Fixed random seg faults when using sqlite getar archives with LZ4
   compression
-  Fixed XZ coupling with ``hoomd.md.integrate.npt`` integration
-  Fixed aspect ratio with non-cubic boxes in
   ``hoomd.hpmc.update.boxmc``

v2.1.6 (2017-04-12)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Document ``hpmc.util.tune_npt``
-  Fix dump.getar.writeJSON usage with MPI execution
-  Fix a bug where integrate.langevin and integrate.brownian correlated
   RNGs between ranks in multiple CPU execution
-  Bump CUB to version 1.6.4 for improved performance on Pascal
   architectures. CUB is now embedded using a git submodule. Users
   upgrading existing git repositories should reinitialize their git
   submodules with ``git submodule update --init``
-  CMake no longer complains when it finds a partial MKL installation.

v2.1.5 (2017-03-09)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fixed a compile error on Mac

v2.1.4 (2017-03-09)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fixed a bug re-enabling disabled integration methods
-  Fixed a bug where adding particle types to the system failed for
   anisotropic pair potentials
-  scipy is no longer required to execute DEM component unit tests
-  Issue a warning when a subsequent call to context.initialize is given
   different arguments
-  DPD now uses the seed from rank 0 to avoid incorrect simulations when
   users provide different seeds on different ranks
-  Miscellaneous documentation updates
-  Defer initialization message until context.initialize
-  Fixed a problem where a momentary dip in TPS would cause walltime
   limited jobs to exit prematurely
-  HPMC and DEM components now correctly print citation notices

v2.1.3 (2017-02-07)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fixed a bug where the WalltimeLimitReached was ignored

v2.1.2 (2017-01-11)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  (HPMC) Implicit depletants with spheres and faceted spheres now
   produces correct ensembles
-  (HPMC) Implicit depletants with ntrial > 0 now produces correct
   ensembles
-  (HPMC) NPT ensemble in HPMC (``hpmc.update.boxmc``) now produces
   correct ensembles
-  Fix a bug where multiple nvt/npt integrators caused warnings from
   analyze.log.
-  update.balance() is properly ignored when only one rank is available
-  Add missing headers to plugin install build
-  Fix a bug where charge.pppm calculated an incorrect pressure

-  Other changes \*

-  Drop support for compute 2.0 GPU devices
-  Support cusolver with CUDA 8.0

v2.1.1 (2016-10-23)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix ``force.active`` memory allocation bug
-  Quiet Python.h warnigns when building (python 2.7)
-  Allow multi-character particle types in HPMC (python 2.7)
-  Enable ``dump.getar.writeJSON`` in MPI
-  Allow the flow to change directions in
   ``md.update.mueller_plathe_flow``
-  Fix critical bug in MPI communication when using HPMC integrators

v2.1.0 (2016-10-04)
^^^^^^^^^^^^^^^^^^^

*New features*

-  enable/disable overlap checks between pairs of constituent particles
   for ``hpmc.integrate.sphere_union()``
-  Support for non-additive mixtures in HPMC, overlap checks can now be
   enabled/disabled per type-pair
-  Add ``md.constrain.oned`` to constrain particles to move in one
   dimension
-  ``hpmc.integrate.sphere_union()`` now takes max_members as an
   optional argument, allowing to use GPU memory more efficiently
-  Add ``md.special_pair.lj()`` to support scaled 1-4 (or other)
   exclusions in all-atom force fields
-  ``md.update.mueller_plathe_flow()``: Method to create shear flows in
   MD simulations
-  ``use_charge`` option for ``md.pair.reaction_field``
-  ``md.charge.pppm()`` takes a Debye screening length as an optional
   parameter
-  ``md.charge.pppm()`` now computes the rigid body correction to the
   PPPM energy

*Deprecated*

-  HPMC: the ``ignore_overlaps`` flag is replaced by
   ``hpmc.integrate.interaction_matrix``

*Other changes*

-  Optimized MPI simulations of mixed systems with rigid and non-rigid
   bodies
-  Removed dependency on all boost libraries. Boost is no longer needed
   to build hoomd
-  Intel compiler builds are no longer supported due to c++11 bugs
-  Shorter compile time for HPMC GPU kernels
-  Include symlinked external components in the build process
-  Add template for external components
-  Optimized dense depletant simulations with HPMC on CPU

*Bug fixes*

-  fix invalid mesh energy in non-neutral systems with
   ``md.charge.pppm()``
-  Fix invalid forces in simulations with many bond types (on GPU)
-  fix rare cases where analyze.log() would report a wrong pressure
-  fix possible illegal memory access when using
   ``md.constrain.rigid()`` in GPU MPI simulations
-  fix a bug where the potential energy is misreported on the first step
   with ``md.constrain.rigid()``
-  Fix a bug where the potential energy is misreported in MPI
   simulations with ``md.constrain.rigid()``
-  Fix a bug where the potential energy is misreported on the first step
   with ``md.constrain.rigid()``
-  ``md.charge.pppm()`` computed invalid forces
-  Fix a bug where PPPM interactions on CPU where not computed correctly
-  Match logged quantitites between MPI and non-MPI runs on first time
   step
-  Fix ``md.pair.dpd`` and ``md.pair.dpdlj`` ``set_params``
-  Fix diameter handling in DEM shifted WCA potential
-  Correctly handle particle type names in lattice.unitcell
-  Validate ``md.group.tag_list`` is consistent across MPI ranks

v2.0.3 (2016-08-30)
^^^^^^^^^^^^^^^^^^^

-  hpmc.util.tune now works with particle types as documented
-  Fix pressure computation with pair.dpd() on the GPU
-  Fix a bug where dump.dcd corrupted files on job restart
-  Fix a bug where HPMC walls did not work correctly with MPI
-  Fix a bug where stdout/stderr did not appear in MPI execution
-  HOOMD will now report an human readable error when users forget
   context.initialize()
-  Fix syntax errors in frenkel ladd field

v2.0.2 (2016-08-09)
^^^^^^^^^^^^^^^^^^^

-  Support CUDA Toolkit 8.0
-  group.rigid()/nonrigid() did not work in MPI simulations
-  Fix builds with ENABLE_DOXYGEN=on
-  Always add -std=c++11 to the compiler command line arguments
-  Fix rare infinite loops when using hpmc.integrate.faceted_sphere
-  Fix hpmc.util.tune to work with more than one tunable
-  Fix a bug where dump.gsd() would write invalid data in simulations
   with changing number of particles
-  replicate() sometimes did not work when restarting a simulation

v2.0.1 (2016-07-15)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix acceptance criterion in mu-V-T simulations with implicit
   depletants (HPMC).
-  References to disabled analyzers, computes, updaters, etc. are
   properly freed from the simulation context.
-  Fix a bug where ``init.read_gsd`` ignored the ``restart`` argument.
-  Report an error when HPMC kernels run out of memory.
-  Fix ghost layer when using rigid constraints in MPI runs.
-  Clarify definition of the dihedral angle.

v2.0.0 (2016-06-22)
^^^^^^^^^^^^^^^^^^^

HOOMD-blue v2.0 is released under a clean BSD 3-clause license.

*New packages*

-  ``dem`` - simulate faceted shapes with dynamics
-  ``hpmc`` - hard particle Monte Carlo of a variety of shape classes.

*Bug fixes*

-  Angles, dihedrals, and impropers no longer initialize with one
   default type.
-  Fixed a bug where integrate.brownian gave the same x,y, and z
   velocity components.
-  Data proxies verify input types and vector lengths.
-  dump.dcd no longer generates excessive metadata traffic on lustre
   file systems

*New features*

-  Distance constraints ``constrain.distance`` - constrain pairs of
   particles to a fixed separation distance
-  Rigid body constraints ``constrain.rigid`` - rigid bodies now have
   central particles, and support MPI and replication
-  Multi-GPU electrostatics ``charge.pppm`` - the long range
   electrostatic forces are now supported in MPI runs
-  ``context.initialize()`` can now be called multiple times - useful in
   jupyter notebooks
-  Manage multiple simulations in a single job script with
   ``SimulationContext`` as a python context manager.
-  ``util.quiet_status() / util.unquiet_status()`` allow users to
   control if line status messages are output.
-  Support executing hoomd in Jupyter (ipython) notebooks. Notice,
   warning, and error messages now show up in the notebook output
   blocks.
-  ``analyze.log`` can now register python callback functions as sources
   for logged quantities.
-  The GSD file format (http://gsd.readthedocs.io) is fully implemented
   in hoomd

   -  ``dump.gsd`` writes GSD trajectories and restart files (use
      ``truncate=true`` for restarts).
   -  ``init.read_gsd`` reads GSD file and initializes the system, and
      can start the simulation from any frame in the GSD file.
   -  ``data.gsd_snapshot`` reads a GSD file into a snapshot which can
      be modified before system initialization with
      ``init.read_snapshot``.
   -  The GSD file format is capable of storing all particle and
      topology data fields in hoomd, either static at frame 0, or
      varying over the course of the trajectory. The number of
      particles, types, bonds, etc. can also vary over the trajectory.

-  ``force.active`` applies an active force (optionally with rotational
   diffusion) to a group of particles
-  ``update.constrain_ellipsoid`` constrains particles to an ellipsoid
-  ``integrate.langevin`` and ``integrate.brownian`` now apply
   rotational noise and damping to anisotropic particles
-  Support dynamically updating groups. ``group.force_update()`` forces
   the group to rebuild according to the original selection criteria.
   For example, this can be used to periodically update a cuboid group
   to include particles only in the specified region.
-  ``pair.reaction_field`` implements a pair force for a screened
   electrostatic interaction of a charge pair in a dielectric medium.
-  ``force.get_energy`` allows querying the potential energy of a
   particle group for a specific force
-  ``init.create_lattice`` initializes particles on a lattice.

   -  ``lattice.unitcell`` provides a generic unit cell definition for
      ``create_lattice``
   -  Convenience functions for common lattices: sq, hex, sc, bcc, fcc.

-  Dump and initialize commands for the GTAR file format
   (http://libgetar.readthedocs.io).

   -  GTAR can store trajectory data in zip, tar, sqlite, or bare
      directories
   -  The current version stores system properties, later versions will
      be able to capture log, metadata, and other output to reduce the
      number of files that a job script produces.

-  ``integrate.npt`` can now apply a constant stress tensor to the
   simulation box.
-  Faceted shapes can now be simulated through the ``dem`` component.

*Changes that require job script modifications*

-  ``context.initialize()`` is now required before any other hoomd
   script command.
-  ``init.reset()`` no longer exists. Use ``context.initialize()`` or
   activate a ``SimulationContext``.
-  Any scripts that relied on undocumented members of the ``globals``
   module will fail. These variables have been moved to the ``context``
   module and members of the currently active ``SimulationContext``.
-  bonds, angles, dihedrals, and impropers no longer use the
   ``set_coeff`` syntax. Use ``bond_coeff.set``, ``angle_coeff.set``,
   ``dihedral_coeff.set``, and ``improper_coeff.set`` instead.
-  ``hoomd_script`` no longer exists, python commands are now spread
   across ``hoomd``, ``hoomd.md``, and other sub packages.
-  ``integrate.\*_rigid()`` no longer exists. Use a standard integrator
   on ``group.rigid_center()``, and define rigid bodies using
   ``constrain.rigid()``
-  All neighbor lists must be explicitly created using ``nlist.\*``, and
   each pair potential must be attached explicitly to a neighbor list. A
   default global neighbor list is no longer created.
-  Moved cgcmm into its own package.
-  Moved eam into the metal package.
-  Integrators now take ``kT`` arguments for temperature instead of
   ``T`` to avoid confusion on the units of temperature.
-  phase defaults to 0 for updaters and analyzers so that restartable
   jobs are more easily enabled by default.
-  ``dump.xml`` (deprecated) requires a particle group, and can dump
   subsets of particles.

*Other changes*

-  CMake minimum version is now 2.8
-  Convert particle type names to ``str`` to allow unicode type name
   input
-  ``__version__`` is now available in the top level package
-  ``boost::iostreams`` is no longer a build dependency
-  ``boost::filesystem`` is no longer a build dependency
-  New concepts page explaining the different styles of neighbor lists
-  Default neighbor list buffer radius is more clearly shown to be
   r_buff = 0.4
-  Memory usage of ``nlist.stencil`` is significantly reduced
-  A C++11 compliant compiler is now required to build HOOMD-blue

*Removed*

-  Removed ``integrate.bdnvt``: use ``integrate.langevin``
-  Removed ``mtk=False`` option from ``integrate.nvt`` - The MTK NVT
   integrator is now the only implementation.
-  Removed ``integrate.\*_rigid()``: rigid body functionality is now
   contained in the standard integration methods
-  Removed the global neighbor list, and thin wrappers to the neighbor
   list in ``nlist``.
-  Removed PDB and MOL2 dump writers.
-  Removed init.create_empty

*Deprecated*

-  Deprecated analyze.msd.
-  Deprecated dump.xml.
-  Deprecated dump.pos.
-  Deprecated init.read_xml.
-  Deprecated init.create_random.
-  Deprecated init.create_random_polymers.

v1.x
----

v1.3.3 (2016-03-06)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix problem incluing ``hoomd.h`` in plugins
-  Fix random memory errors when using walls

v1.3.2 (2016-02-08)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix wrong access to system.box
-  Fix kinetic energy logging in MPI
-  Fix particle out of box error if particles are initialized on the
   boundary in MPI
-  Add integrate.brownian to the documentation index
-  Fix misc doc typos
-  Fix runtime errors with boost 1.60.0
-  Fix corrupt metadata dumps in MPI runs

v1.3.1 (2016-1-14)
^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix invalid MPI communicator error with Intel MPI
-  Fix python 3.5.1 seg fault

v1.3.0 (2015-12-8)
^^^^^^^^^^^^^^^^^^

*New features*

-  Automatically load balanced domain decomposition simulations.
-  Anisotropic particle integrators.
-  Gay-Berne pair potential.
-  Dipole pair potential.
-  Brownian dynamics ``integrate.brownian``
-  Langevin dynamics ``integrate.langevin`` (formerly ``bdnvt``)
-  ``nlist.stencil`` to compute neighbor lists using stencilled cell
   lists.
-  Add single value scale, ``min_image``, and ``make_fraction`` to
   ``data.boxdim``
-  ``analyze.log`` can optionally not write a file and now supports
   querying current quantity values.
-  Rewritten wall potentials.

   -  Walls are now sums of planar, cylindrical, and spherical
      half-spaces.
   -  Walls are defined and can be modified in job scripts.
   -  Walls execute on the GPU.
   -  Walls support per type interaction parameters.
   -  Implemented for: lj, gauss, slj, yukawa, morse, force_shifted_lj,
      and mie potentials.

-  External electric field potential: ``external.e_field``

*Bug fixes*

-  Fixed a bug where NVT integration hung when there were 0 particles in
   some domains.
-  Check SLURM environment variables for local MPI rank identification
-  Fixed a typo in the box math documentation
-  Fixed a bug where exceptions weren’t properly passed up to the user
   script
-  Fixed a bug in the velocity initialization example
-  Fixed an openmpi fork() warning on some systems
-  Fixed segfaults in PPPM
-  Fixed a bug where compute.thermo failed after reinitializing a system
-  Support list and dict-like objects in init.create_random_polymers.
-  Fall back to global rank to assign GPUs if local rank is not
   available

*Deprecated commands*

-  ``integrate.bdnvt`` is deprecated. Use ``integrate.langevin``
   instead.
-  ``dump.bin`` and ``init.bin`` are now removed. Use XML files for
   restartable jobs.

*Changes that may break existing scripts*

-  ``boxdim.wrap`` now returns the position and image in a tuple, where
   it used to return just the position.
-  ``wall.lj`` has a new API
-  ``dump.bin`` and ``init.bin`` have been removed.

v1.2.1 (2015-10-22)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix a crash when adding or removing particles and reinitializing
-  Fix a bug where simulations hung on sm 5.x GPUs with CUDA 7.5
-  Fix compile error with long tests enabled
-  Issue a warning instead of an error for memory allocations greater
   than 4 GiB.
-  Fix invalid RPATH when building inside ``zsh``.
-  Fix incorrect simulations with ``integrate.npt_rigid``
-  Label mie potential correctly in user documentation

v1.2.0 (2015-09-30)
^^^^^^^^^^^^^^^^^^^

*New features*

-  Performance improvements for systems with large particle size
   disparity
-  Bounding volume hierarchy (tree) neighbor list computation
-  Neighbor lists have separate ``r_cut`` values for each pair of types
-  addInfo callback for dump.pos allows user specified information in
   pos files

*Bug fixes*

-  Fix ``test_pair_set_energy`` unit test, which failed on numpy < 1.9.0
-  Analyze.log now accepts unicode strings.
-  Fixed a bug where calling ``restore_snapshot()`` during a run zeroed
   potential parameters.
-  Fix segfault on exit with python 3.4
-  Add ``cite.save()`` to documentation
-  Fix a problem were bond forces are computed incorrectly in some MPI
   configurations
-  Fix bug in pair.zbl
-  Add pair.zbl to the documentation
-  Use ``HOOMD_PYTHON_LIBRARY`` to avoid problems with modified CMake
   builds that preset ``PYTHON_LIBRARY``

v1.1.1 (2015-07-21)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  ``dump.xml(restart=True)`` now works with MPI execution
-  Added missing documentation for ``meta.dump_metadata``
-  Build all unit tests by default
-  Run all script unit tests through ``mpirun -n 1``

v1.1.0 (2015-07-14)
^^^^^^^^^^^^^^^^^^^

*New features*

-  Allow builds with ninja.
-  Allow K=0 FENE bonds.
-  Allow number of particles types to change after initialization.

   .. code::

       system.particles.types.add('newtype')

-  Allow number of particles to change after initialization.

   .. code::

       system.particles.add(‘A’)
       del system.particles[0]

-  OPLS dihedral
-  Add ``phase`` keyword to analyzers and dumps to make restartable jobs easier.
-  ``HOOMD_WALLTIME_STOP`` environment variable to stop simulation runs before they hit a wall clock limit.
-  ``init.read_xml()`` Now accepts an initialization and restart file.
-  ``dump.xml()`` can now write restart files.
-   Added documentation concepts page on writing restartable jobs.
-   New citation management infrastructure. ``cite.save()`` writes ``.bib`` files with a list of references to
    features actively used in the current job script.
-   Snapshots expose data as numpy arrays for high performance access to particle properties.
-  ``data.make_snapshot()`` makes a new empty snapshot.
-  ``analyze.callback()`` allows multiple python callbacks to operate at different periods.
-  ``comm.barrier()``and`` comm.barrier_all()``allow users to insert barriers into their scripts.
-   Mie pair potential.
-  ``meta.dump_metadata()`` writes job metadata information out to a json file.
-  ``context.initialize()`` initializes the execution context.
-  Restart option for ``dump.xml()``

*Bug fixes*

-  Fix slow performance when initializing ``pair.slj()``\ in MPI runs.
-  Properly update particle image when setting position from python.
-  PYTHON_SITEDIR hoomd shell launcher now calls the python interpreter
   used at build time.
-  Fix compile error on older gcc versions.
-  Fix a bug where rigid bodies had 0 velocity when restarting jobs.
-  Enable ``-march=native`` builds in OS X clang builds.
-  Fix ``group.rigid()`` and ``group.nonrigid()``.
-  Fix image access from the python data access proxies.
-  Gracefully exit when launching MPI jobs with mixed execution
   configurations.

*Changes that may require updated job scripts*

-  ``context.initialize()`` **must** be called before any ``comm``
   method that queries the MPI rank. Call it as early as possible in
   your job script (right after importing ``hoomd_script``) to avoid
   problems.

*Deprecated*

-  ``init.create_empty()`` is deprecated and will be removed in a future
   version. Use ``data.make_snapshot()`` and ``init.read_snapshot()``
   instead.
-  Job scripts that do not call ``context.initialize()`` will result in
   a warning message. A future version of HOOMD will require that you
   call ``context.initialize()``.

*Removed*

-  Several ``option`` commands for controlling the execution
   configuration. Replaced with ``context.initialize``.

v1.0.5 (2015-05-19)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix segfault when changing integrators
-  Fix system.box to indicate the correct number of dimensions
-  Fix syntax error in comm.get_rank with –nrank
-  Enable CUDA enabled builds with the intel compiler
-  Use CMake builtin FindCUDA on recent versions of CMake
-  GCC_ARCH env var sets the -march command line option to gcc at
   configure time
-  Auto-assign GPU-ids on non-compute exclusive systems even with
   –mode=gpu
-  Support python 3.5 alpha
-  Fix a bug where particle types were doubled with boost 1.58.0
-  Fix a bug where angle_z=true dcd output was inaccurate near 0 angles
-  Properly handle lj.wall potentials with epsilon=0.0 and particles on
   top of the walls

v1.0.4 (2015-04-07)
^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fix invalid virials computed in rigid body simulations when
   multi-particle bodies crossed box boundaries
-  Fix invalid forces/torques for rigid body simulations caused by race
   conditions
-  Fix compile errors on Mac OS X 10.10
-  Fix invalid pair force computations caused by race conditions
-  Fix invalid neighbour list computations caused by race conditions on
   Fermi generation GPUs

*Other*

-  Extremely long running unit tests are now off by default. Enable with
   -DHOOMD_SKIP_LONG_TESTS=OFF
-  Add additional tests to detect race conditions and memory errors in
   kernels

v1.0.3 (2015-03-18)
^^^^^^^^^^^^^^^^^^^

**Bug fixes**

-  Enable builds with intel MPI
-  Silence warnings coming from boost and python headers

v1.0.2 (2015-01-21)
^^^^^^^^^^^^^^^^^^^

**Bug fixes**

-  Fixed a bug where ``linear_interp`` would not take a floating point
   value for *zero*
-  Provide more useful error messages when cuda drivers are not present
-  Assume device count is 0 when ``cudaGetDeviceCount()`` returns an
   error
-  Link to python statically when ``ENABLE_STATIC=on``
-  Misc documentation updates

v1.0.1 (2014-09-09)
^^^^^^^^^^^^^^^^^^^

**Bug fixes**

1.  Fixed bug where error messages were truncated and HOOMD exited with
    a segmentation fault instead (e.g. on Blue Waters)
2.  Fixed bug where plug-ins did not load on Blue Waters
3.  Fixed compile error with gcc4.4 and cuda5.0
4.  Fixed syntax error in ``read_snapshot()``
5.  Fixed a bug where ``init.read_xml throwing`` an error (or any other
    command outside of ``run()``) would hang in MPI runs
6.  Search the install path for hoomd_script - enable the hoomd
    executable to be outside of the install tree (useful with cray
    aprun)
7.  Fixed CMake 3.0 warnings
8.  Removed dependancy on tr1/random
9.  Fixed a bug where ``analyze.msd`` ignored images in the r0_file
10. Fixed typos in ``pair.gauss`` documentation
11. Fixed compile errors on Ubuntu 12.10
12. Fix failure of ``integrate.nvt`` to reach target temperature in
    analyze.log. The fix is a new symplectic MTK integrate.nvt
    integrator. Simulation results in hoomd v1.0.0 are correct, just the
    temperature and velocity outputs are off slightly.
13. Remove MPI from Mac OS X dmg build.
14. Enable ``import hoomd_script as ...``

*Other changes*

1. Added default compile flag -march=native
2. Support CUDA 6.5
3. Binary builds for CentOS/RHEL 6, Fedora 20, Ubuntu 14.04 LTS, and
   Ubuntu 12.04 LTS.

Version 1.0.0 (2014-05-25)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*New features*

-  Support for python 3
-  New NPT integrator capable of flexible coupling schemes
-  Triclinic unit cell support
-  MPI domain decomposition
-  Snapshot save/restore
-  Autotune block sizes at run time
-  Improve performance in small simulation boxes
-  Improve performance with smaller numbers of particles per GPU
-  Full double precision computations on the GPU (compile time option
   must be enabled, binary builds provided on the download page are
   single precision)
-  Tabulated bond potential ``bond.table``
-  Tabulated angle potential ``angle.table``
-  Tabulated dihedral potental ``dihedral.table``
-  ``update.box_resize`` now accepts ``period=None`` to trigger an
   immediate update of the box without creating a periodic updater
-  ``update.box_resize`` now replaces *None* arguments with the current
   box parameters
-  ``init.create_random`` and ``init.create_random_polymers`` can now
   create random configurations in triclinc and 2D boxes
-  ``init.create_empty`` can now create triclinic boxes
-  particle, bond, angle, dihedral, and impropers types can now be named
   in ``init.create_empty``
-  ``system.replicate`` command replicates the simulation box

*Bug fixes*

-  Fixed a bug where init.create_random_polymers failed when lx,ly,lz
   were not equal.
-  Fixed a bug in init.create_random_polymers and init.create_random
   where the separation radius was not accounted for correctly
-  Fixed a bug in bond.\* where random crashes would occur when more
   than one bond type was defined
-  Fixed a bug where dump.dcd did not write the period to the file

*Changes that may require updated job scripts*

-  ``integrate.nph``: A time scale ``tau_p`` for the relaxation of the
   barostat is now required instead of the barostat mass *W* of the
   previous release. The time scale is the relaxation time the barostat
   would have at an average temperature ``T_0 = 1``, and it is related
   to the internally used (Andersen) Barostat mass *W* via
   ``W = d N T_0 tau_p^2``, where *d* is the dimensionsality and *N* the
   number of particles.
-  ``sorter`` and ``nlist`` are now modules, not variables in the
   ``__main__`` namespace.
-  Data proxies function correctly in MPI simulations, but are extremely
   slow. If you use ``init.create_empty``, consider separating the
   generation step out to a single rank short execution that writes an
   XML file for the main run.
-  ``update.box_resize(Lx=...)`` no longer makes cubic box updates,
   instead it will keep the current **Ly** and **Lz**. Use the ``L=...``
   shorthand for cubic box updates.
-  All ``init.*`` commands now take ``data.boxdim`` objects, instead of
   ``hoomd.boxdim`` (or *3-tuples*). We strongly encourage the use of
   explicit argument names for ``data.boxdim()``. In particular, if
   ``hoomd.boxdim(123)`` was previously used to create a cubic box, it
   is now required to use ``data.boxdim(L=123)`` (CORRECT) instead of
   ``data.boxdim(123)`` (INCORRECT), otherwise a box with unit
   dimensions along the y and z axes will be created.
-  ``system.dimensions`` can no longer be set after initialization.
   System dimensions are now set during initialization via the
   ``data.boxdim`` interface. The dimensionality of the system can now
   be queried through ``system.box``.
-  ``system.box`` no longer accepts 3-tuples. It takes ``data.boxdim``
   objects.
-  ``system.dimensions`` no longer exists. Query the dimensionality of
   the system from ``system.box``. Set the dimensionality of the system
   by passing an appropriate ``data.boxdim`` to an ``init`` method.
-  ``init.create_empty`` no longer accepts ``n_*_types``. Instead, it
   now takes a list of strings to name the types.

*Deprecated*

-  Support for G80, G200 GPUs.
-  ``dump.bin`` and ``read.bin``. These will be removed in v1.1 and
   replaced with a new binary format.

*Removed*

-  OpenMP mult-core execution (replaced with MPI domain decomposition)
-  ``tune.find_optimal_block_size`` (replaced by Autotuner)

v0.x
----

Version 0.11.3 (2013-05-10)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*Bug fixes*

-  Fixed a bug where charge.pppm could not be used after init.reset()
-  Data proxies can now set body angular momentum before the first run()
-  Fixed a bug where PPPM forces were incorrect on the GPU

Version 0.11.2 (2012-12-19)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*New features*

-  Block sizes tuned for K20

*Bug fixes*

-  Warn user that PPPM ignores rigid body exclusions
-  Document that proxy iterators need to be deleted before init.reset()
-  Fixed a bug where body angular momentum could not be set
-  Fixed a bug where analyze.log would report nan for the pressure
   tensor in nve and nvt simulations

Version 0.11.1 (2012-11-2)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*New features*

-  Support for CUDA 5.0
-  Binary builds for Fedora 16 and OpenSUSE 12.1
-  Automatically specify /usr/bin/gcc to nvcc when the configured gcc is
   not supported

*Bug fixes*

-  Fixed a compile error with gcc 4.7
-  Fixed a bug where PPPM forces were incorrect with neighborlist
   exclusions
-  Fixed an issue where boost 1.50 and newer were not detected properly
   when BOOST_ROOT is set
-  Fixed a bug where accessing force data in python prevented
   init.reset() from working
-  Fixed a bug that prevented pair.external from logging energy
-  Fixed a unit test that failed randomly

Version 0.11.0 (2012-07-27)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*New features*

1.  Support for Kepler GPUs (GTX 680)
2.  NPH integration (*integrate.nph*)
3.  Compute full pressure tensor
4.  Example plugin for new bond potentials
5.  New syntax for bond coefficients: *bond.bond_coeff.set(‘type’,
    params)*
6.  New external potential: *external.periodic* applies a periodic
    potential along one direction (uses include inducing lamellar phases
    in copolymer systems)
7.  Significant performance increases when running *analyze.log*,
    *analyze.msd*, *update.box_resize*, *update.rescale_temp*, or
    *update.zero_momentum* with a small period
8.  Command line options may now be overwritten by scripts, ex:
    *options.set_gpu(2)*
9.  Added *–user* command line option to allow user defined options to
    be passed into job scripts, ex: *–user=“-N=5 -phi=0.56”*
10. Added *table.set_from_file* method to enable reading table based
    pair potentials from a file
11. Added *–notice-level* command line option to control how much extra
    information is printed during a run. Set to 0 to disable, or any
    value up to 10. At 10, verbose debugging information is printed.
12. Added *–msg-file* command line option which redirects the message
    output to a file
13. New pair potential *pair.force_shifted_lj* : Implements
    http://dx.doi.org/10.1063/1.3558787

*Bug fixes*

1. Fixed a bug where FENE bonds were sometimes computed incorrectly
2. Fixed a bug where pressure was computed incorrectly when using
   pair.dpd or pair.dpdlj
3. Fixed a bug where using OpenMP and CUDA at the same time caused
   invalid memory accesses
4. Fixed a bug where RPM packages did not work on systems where the CUDA
   toolkit was not installed
5. Fixed a bug where rigid body velocities were not set from python
6. Disabled OpenMP builds on Mac OS X. HOOMD-blue w/ openmp enabled
   crashes due to bugs in Apple’s OpenMP implementation.
7. Fixed a bug that allowed users to provide invalid rigid body data and
   cause a seg fault.
8. Fixed a bug where using PPPM resulted in error messages on program
   exit.

*API changes*

1.  Bond potentials rewritten with template evaluators
2.  External potentials use template evaluators
3.  Complete rewrite of ParticleData - may break existing plugins
4.  Bond/Angle/Dihedral data structures rewritten

    -  The GPU specific data structures are now generated on the GPU

5.  DPDThermo and DPDLJThermo are now processed by the same template
    class
6.  Headers that cannot be included by nvcc now throw an error when they
    are
7.  CUDA 4.0 is the new minimum requirement
8.  Rewrote BoxDim to internally handle minimum image conventions
9.  HOOMD now only compiles ptx code for the newest architecture, this
    halves the executable file size
10. New Messenger class for global control of messages printed to the
    screen / directed to a file.

*Testing changes*

1. Automated test suite now performs tests on OpenMPI + CUDA builds
2. Valgrind tests added back into automated test suite
3. Added CPU test in bd_ridid_updater_tests
4. ctest -S scripts can now set parallel makes (with cmake > 2.8.2)

Version 0.10.1 (2012-02-10)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Add missing entries to credits page
2. Add ``dist_check`` option to neighbor list. Can be used to force
   neighbor list builds at a specified frequency (useful in profiling
   runs with nvvp).
3. Fix typos in ubuntu compile documentation
4. Add missing header files to hoomd.h
5. Add torque to the python particle data access API
6. Support boost::filesystem API v3
7. Expose name of executing gpu, n_cpu, hoomd version, git sha1, cuda
   version, and compiler version to python
8. Fix a bug where multiple ``nvt_rigid`` or ``npt_rigid`` integrators
   didn’t work correctly
9. Fix missing pages in developer documentation

Version 0.10.0 (2011-12-14)
^^^^^^^^^^^^^^^^^^^^^^^^^^^

*New features*

1.  Added *pair.dpdlj* which uses the DPD thermostat and the
    Lennard-Jones potential. In previous versions, this could be
    accomplished by using two pair commands but at the cost of reduced
    performance.
2.  Additional example scripts are now present in the documentation. The
    example scripts are cross-linked to the commands that are used in
    them.
3.  Most dump commands now accept the form:
    *dump.ext(filename=“filename.ext”)* which immediately writes out
    filename.ext.
4.  Added *vis* parameter to dump.xml which enables output options
    commonly used in files written for the purposes of visulization.
    dump.xml also now accepts parameters on the instantiation line.
    Combined with the previous feature, *dump.xml(filename=“file.xml”,
    vis=True)* is now a convenient short hand for what was previously

    .. code::

       xml = dump.xml()
       xml.set_params(position = True, mass = True, diameter = True,
                             type = True, bond = True, angle = True,
                             dihedral = True, improper = True, charge = True)
       xml.write(filename="file.xml")

5.  Specify rigid bodies in XML input files
6.  Simulations that contain rigid body constraints applied to groups of
    particles in BDNVT, NVE, NVT, and NPT ensembles.

    -  *integrate.bdnvt_rigid*
    -  *integrate.nve_rigid*
    -  *integrate.nvt_rigid*
    -  *integrate.npt_rigid*

7.  Energy minimization of rigid bodies
    (*integrate.mode_minimize_rigid_fire*)
8.  Existing commands are now rigid-body aware

    -  update.rescale_temp
    -  update.box_resize
    -  update.enforce2d
    -  update.zero_momentum

9.  NVT integration using the Berendsen thermostat
    (*integrate.berendsen*)
10. Bonds, angles, dihedrals, and impropers can now be created and
    deleted with the python data access API.
11. Attribution clauses added to the HOOMD-blue license.

*Changes that may break existing job scripts*

1. The *wrap* option to *dump.dcd* has been changed to *unwrap_full* and
   its meaning inverted. *dump.dcd* now offers two options for
   unwrapping particles, *unwrap_full* fully unwraps particles into
   their box image and *unwrap_rigid* unwraps particles in rigid bodies
   so that bodies are not broken up across a box boundary.

*Bug/fixes small enhancements*

1.  Fixed a bug where launching hoomd on mac os X 10.5 always resulted
    in a bus error.
2.  Fixed a bug where DCD output restricted to a group saved incorrect
    data.
3.  force.constant may now be applied to a group of particles, not just
    all particles
4.  Added C++ plugin example that demonstrates how to add a pair
    potential in a plugin
5.  Fixed a bug where box.resize would always transfer particle data
    even in a flat portion of the variant
6.  OpenMP builds re-enabled on Mac OS X
7.  Initial state of integrate.nvt and integrate.npt changed to decrease
    oscillations at startup.
8.  Fixed a bug where the polymer generator would fail to initialize
    very long polymers
9.  Fixed a bug where images were passed to python as unsigned ints.
10. Fixed a bug where dump.pdb wrote coordinates in the wrong order.
11. Fixed a rare problem where a file written by dump.xml would not be
    read by init.read_xml due to round-off errors.
12. Increased the number of significant digits written out to dump.xml
    to make them more useful for ad-hoc restart files.
13. Potential energy and pressure computations that slow performance are
    now only performed on those steps where the values are actually
    needed.
14. Fixed a typo in the example C++ plugin
15. Mac build instructions updated to work with the latest version of
    macports
16. Fixed a bug where set_period on any dump was ineffective.
17. print_status_line now handles multiple lines
18. Fixed a bug where using bdnvt tally with per type gammas resulted in
    a race condition.
19. Fix an issue where ENABLE_CUDA=off builds gave nonsense errors when
    –mode=gpu was requested.
20. Fixed a bug where dumpl.xml could produce files that init.xml would
    not read
21. Fixed a typo in the example plugin
22. Fix example that uses hoomd as a library so that it compiles.
23. Update maintainer lines
24. Added message to nlist exclusions that notifies if diameter or body
    exclusions are set.
25. HOOMD-blue is now hosted in a git repository
26. Added bibtex bibliography to the user documentation
27. Converted user documentation examples to use doxygen auto
    cross-referencing ``\example`` commands
28. Fix a bug where particle data is not released in dump.binary
29. ENABLE_OPENMP can now be set in the ctest builds
30. Tuned block sizes for CUDA 4.0
31. Removed unsupported GPUS from CUDA_ARCH_LIST

Version 0.9.2 (2011-04-04)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1. *New exclusion option:* Particles can now be excluded from the
   neighbor list based on diameter consistent with pair.slj.
2. *New pair coeff syntax:* Coefficients for multiple type pairs can be
   specified conveniently on a single line.

   .. code::

      coeff.set(['A', 'B', 'C', 'D'], ['A', 'B', 'C', 'D'], epsilon=1.0)

3. *New documentation:* HOOMD-blue’s system of units is now fully
   documented, and every coefficient in the documentation is labeled
   with the appropriate unit.
4. *Performance improvements:* Performance has been significantly
   boosted for simulations of medium sized systems (5,000-20,000
   particles). Smaller performance boosts were made to larger runs.
5. *CUDA 3.2 support:* HOOMD-blue is now fully tested and performance
   tuned for use with CUDA 3.2.
6. *CUDA 4.0 support:* HOOMD-blue compiles with CUDA 4.0 and passes
   initial tests.
7. *New command:* tune.r_buff performs detailed auto-tuning of the
   r_buff neighborlist parameter.
8. *New installation method:* RPM, DEB, and app bundle packages are now
   built for easier installation
9. *New command:* charge.pppm computes the full long range electrostatic
   interaction using the PPPM method

*Bug/fixes small enhancements*

1.  Fixed a bug where the python library was linked statically.
2.  Added the PYTHON_SITEDIR setting to allow hoomd builds to install
    into the native python site directory.
3.  FIRE energy minimization convergence criteria changed to require
    both energy *and* force to converge
4.  Clarified that groups are static in the documentation
5.  Updated doc comments for compatibility with Doxygen#7.3
6.  system.particles.types now lists the particle types in the
    simulation
7.  Creating a group of a non-existant type is no longer an error
8.  Mention XML file format for walls in wall.lj documentation
9.  Analyzers now profile themselves
10. Use ``\n`` for newlines in dump.xml - improves
    performance when writing many XML files on a NFS file system
11. Fixed a bug where the neighbor list build could take an
    exceptionally long time (several seconds) to complete the first
    build.
12. Fixed a bug where certain logged quantities always reported as 0 on
    the first step of the simulation.
13. system.box can now be used to read and set the simulation box size
    from python
14. Numerous internal API updates
15. Fixed a bug the resulted in incorrect behavior when using
    integrate.npt on the GPU.
16. Removed hoomd launcher shell script. In non-sitedir installs,
    ${HOOMD_ROOT}/bin/hoomd is now the executable itself
17. Creating unions of groups of non-existent types no longer produces a
    seg fault
18. hoomd now builds on all cuda architectures. Modify CUDA_ARCH_LIST in
    cmake to add or remove architectures from the build
19. hoomd now builds with boost#46.0
20. Updated hoomd icons to maize/blue color scheme
21. hoomd xml file format bumped to#3, adds support for charge.
22. FENE and harmonic bonds now handle 0 interaction parameters and 0
    length bonds more gracefully
23. The packaged plugin template now actually builds and installs into a
    recent build of hoomd

Version 0.9.1 (2010-10-08)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1. *New constraint*: constrain.sphere constrains a group of particles to
   the surface of a sphere
2. *New pair potential/thermostat*: pair.dpd implements the standard DPD
   conservative, random, and dissipative forces
3. *New pair potential*: pair.dpd_conservative applies just the
   conservative DPD potential
4. *New pair potential*: pair.eam implements the Embedded Atom Method
   (EAM) and supports both *alloy* and *FS* type computations.
5. *Faster performance*: Cell list and neighbor list code has been
   rewritten for performance.

   -  In our benchmarks, *performance increases* ranged from *10-50%*
      over HOOMD-blue 0.9.0. Simulations with shorter cutoffs tend to
      attain a higher performance boost than those with longer cutoffs.
   -  We recommended that you *re-tune r_buff* values for optimal
      performance with 0.9.1.
   -  Due to the nature of the changes, *identical runs* may produce
      *different trajectories*.

6. *Removed limitation*: The limit on the number of neighbor list
   exclusions per particle has been removed. Any number of exclusions
   can now be added per particle. Expect reduced performance when adding
   excessive numbers of exclusions.

*Bug/fixes small enhancements*

1.  Pressure computation is now correct when constraints are applied.
2.  Removed missing files from hoomd.h
3.  pair.yukawa is no longer referred to by “gaussian” in the
    documentation
4.  Fermi GPUs are now prioritized over per-Fermi GPUs in systems where
    both are present
5.  HOOMD now compiles against CUDA 3.1
6.  Momentum conservation significantly improved on compute#x hardware
7.  hoomd plugins can now be installed into user specified directories
8.  Setting r_buff=0 no longer triggers exclusion list updates on every
    step
9.  CUDA 2.2 and older are no longer supported
10. Workaround for compiler bug in 3.1 that produces extremely high
    register usage
11. Disabled OpenMP compile checks on Mac OS X
12. Support for compute 2.1 devices (such as the GTX 460)

Version 0.9.0 (2010-05-18)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1.  *New pair potential*: Shifted LJ potential for particles of varying
    diameters (pair.slj)
2.  *New pair potential*: Tabulated pair potential (pair.table)
3.  *New pair potential*: Yukawa potential (pair.yukawa)
4.  *Update to pair potentials*: Most pair potentials can now accept
    different values of r_cut for different type pairs. The r_cut
    specified in the initial pair.**\* command is now treated as the
    default r_cut, so no changes to scripts are necessary.
5.  *Update to pair potentials*: Default pair coeff values are now
    supported. The parameter alpha for lj now defaults to#0, so there is
    no longer a need to specify it for a majority of simulations.
6.  *Update to pair potentials*: The maximum r_cut needed for the
    neighbor list is now determined at the start of each run(). In
    simulations where r_cut may decrease over time, increased
    performance will result.
7.  *Update to pair potentials*: Pair potentials are now specified via
    template evaluator classes. Adding a new pair potential to hoomd now
    only requires a small amount of additional code.
8.  *Plugin API* : Advanced users/developers can now write, install, and
    use plugins for hoomd without needing to modify core hoomd source
    code
9.  *Particle data access*: User-level hoomd scripts can now directly
    access the particle data. For example, one can change all particles
    in the top half of the box to be type B:

    .. code::

       top = group.cuboid(name="top", zmin=0)
       for p in top:
           p.type = 'B'

    . *All* particle data including position, velocity, type, ‘’et
    cetera’’, can be read and written in this manner. Computed forces
    and energies can also be accessed in a similar way.
10. *New script command*: init.create_empty() can be used in conjunction
    with the particle data access above to completely initialize a
    system within the hoomd script.
11. *New script command*: dump.bin() writes full binary restart files
    with the entire system state, including the internal state of
    integrators.

    -  File output can be gzip compressed (if zlib is available) to save
       space
    -  Output can alternate between two different output files for safe
       crash recovery

12. *New script command*: init.read_bin() reads restart files written by
    dump.bin()
13. *New option*: run() now accepts a quiet option. When True, it
    eliminates the status information printouts that go to stdout.
14. *New example script*: Example 6 demonstrates the use of the particle
    data access routines to initialize a system. It also demonstrates
    how to initialize velocities from a gaussian distribution
15. *New example script*: Example 7 plots the pair.lj potential energy
    and force as evaluated by hoomd. It can trivially be modified to
    plot any potential in hoomd.
16. *New feature*: Two dimensional simulations can now be run in hoomd:
    #259
17. *New pair potential*: Morse potential for particles of varying
    diameters (pair.morse)
18. *New command*: run_upto will run a simulation up to a given time
    step number (handy for breaking long simulations up into many
    independent jobs)
19. *New feature*: HOOMD on the CPU is now accelerated with OpenMP.
20. *New feature*: integrate.mode_minimize_fire performs energy
    minimization using the FIRE algorithm
21. *New feature*: analyze.msd can now accept an xml file specifying the
    initial particle positions (for restarting jobs)
22. *Improved feature*: analyze.imd now supports all IMD commands that
    VMD sends (pause, kill, change trate, etc.)
23. *New feature*: Pair potentials can now be given names, allowing
    multiple potentials of the same type to be logged separately.
    Additionally, potentials that are disabled and not applied to the
    system dynamics can be optionally logged.
24. *Performance improvements*: Simulation performance has been
    increased across the board, but especially when running systems with
    very low particle number densities.
25. *New hardware support*: 0.9.0 and newer support Fermi GPUs
26. *Deprecated hardware support*: 0.9.x might continue run on compute#1
    GPUs but that hardware is no longer officially supported
27. *New script command*: group.tag_list() takes a python list of
    particle tags and creates a group
28. *New script command*: compute.thermo() computes thermodynamic
    properties of a group of particles for logging
29. *New feature*: dump.dcd can now optionally write out only those
    particles that belong to a specified group

*Changes that will break jobs scripts written for 0.8.x*

1. Integration routines have changed significantly to enable new use
   cases. Where scripts previously had commands like:

   .. code::

      integrate.nve(dt=0.005)

   they now need

   .. code::

      all = group.all()
      integrate.mode_standard(dt=0.005)
      integrate.nve(group=all)

   . Integrating only specific groups of particles enables simulations
   to fix certain particles in place or integrate different parts of the
   system at different temperatures, among many other possibilities.
2. sorter.set_params no longer takes the ‘’bin_width’’ argument. It is
   replaced by a new ‘’grid’’ argument, see the documentation for
   details.
3. conserved_quantity is no longer a quantity available for logging.
   Instead log the nvt reservoir energy and compute the total conserved
   quantity in post processing.

*Bug/fixes small enhancements*

1.  Fixed a bug where boost#38 is not found on some machines
2.  dump.xml now has an option to write particle accelerations
3.  Fixed a bug where periods like 1e6 were not accepted by updaters
4.  Fixed a bug where bond.fene forces were calculated incorrectly
    between particles of differing diameters
5.  Fixed a bug where bond.fene energies were computed incorrectly when
    running on the GPU
6.  Fixed a bug where comments in hoomd xml files were not ignored as
    they aught to be: #331
7.  It is now possible to prevent bond exclusions from ever being added
    to the neighbor list: #338
8.  init.create_random_polymers can now generate extremely dense systems
    and will warn the user about large memory usage
9.  variant.linear_interp now accepts a user-defined zero (handy for
    breaking long simulations up into many independent jobs)
10. Improved installation and compilation documentation
11. Integration methods now silently ignore when they are given an empty
    group
12. Fixed a bug where disabling all forces resulted in some forces still
    being applied
13. Integrators now behave in a reasonable way when given empty groups
14. Analyzers now accept a floating point period
15. run() now aborts immediately if limit_hours=0 is specified.
16. Pair potentials that diverge at r=0 will no longer result in invalid
    simulations when the leading coefficients are set to zero.
17. integrate.bdnvt can now tally the energy transferred into/out of the
    “reservoir”, allowing energy conservation to be monitored during bd
    simulation runs.
18. Most potentials now prevent NaN results when computed for
    overlapping particles
19. Stopping a simulation from a callback or time limit no longer
    produces invalid simulations when continued
20. run() commands limited with limit_hours can now be set to only stop
    on given timestep multiples
21. Worked around a compiler bug where pair.morse would crash on Fermi
    GPUs
22. ULF stability improvements for G200 GPUs.

Version 0.8.2 (2009-09-10)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1.  Quantities that vary over time can now be specified easily in
    scripts with the variant.linear_interp command.
2.  Box resizing updater (update.box_resize) command that uses the time
    varying quantity command to grow or shrink the simulation box.
3.  Individual run() commands can be limited by wall-clock time
4.  Angle forces can now be specified
5.  Dihedral forces can now be specified
6.  Improper forces can now be specified
7.  1-3 and 1-4 exclusions from the cutoff pair force can now be chosen
8.  New command line option: –minimize-cpu-usage cuts the CPU usage of
    HOOMD down to 10% of one CPU core while only decreasing overall
    performance by 10%
9.  Major changes have been made in the way HOOMD chooses the device on
    which to run (all require CUDA 2.2 or newer)

    -  there are now checks that an appropriate NVIDIA drivers is
       installed
    -  running without any command line options will now correctly
       revert to running on the CPU if no capable GPUs are installed
    -  when no gpu is explicitly specified, the default choice is now
       prioritized to choose the fastest GPU and one that is not
       attached to a display first
    -  new command line option: –ignore-display-gpu will prevent HOOMD
       from executing on any GPU attached to a display
    -  HOOMD now prints out a short description of the GPU(s) it is
       running on
    -  on linux, devices can be set to compute-exclusive mode and HOOMD
       will then automatically choose the first free GPU (see the
       documentation for details)

10. nlist.reset_exclusions command to control the particles that are
    excluded from the neighbor list

*Bug/fixes small enhancements*

1.  Default block size change to improve stability on compute#3 devices
2.  ULF workaround on GTX 280 now works with CUDA 2.2
3.  Standalone benchmark executables have been removed and replaced by
    in script benchmarking commands
4.  Block size tuning runs can now be performed automatically using the
    python API and results can be saved on the local machine
5.  Fixed a bug where GTX 280 bug workarounds were not properly applied
    in CUDA 2.2
6.  The time step read in from the XML file can now be optionally
    overwritten with a user-chosen one
7.  Added support for CUDA 2.2
8.  Fixed a bug where the WCA forces included in bond.fene had an
    improper cutoff
9.  Added support for a python callback to be executed periodically
    during a run()
10. Removed demos from the hoomd downloads. These will be offered
    separately on the webpage now to keep the required download size
    small.
11. documentation improvements
12. Significantly increased performance of dual-GPU runs when build with
    CUDA 2.2 or newer
13. Numerous stability and performance improvements
14. Temperatures are now calculated based on 3N-3 degrees of freedom.
    See #283 for a more flexible system that is coming in the future.
15. Emulation mode builds now work on systems without an NVIDIA card
    (CUDA 2.2 or newer)
16. HOOMD now compiles with CUDA 2.3
17. Fixed a bug where uninitialized memory was written to dcd files
18. Fixed a bug that prevented the neighbor list on the CPU from working
    properly with non-cubic boxes
19. There is now a compile time hack to allow for more than 4 exclusions
    per particle
20. Documentation added to aid users in migrating from LAMMPS
21. hoomd_script now has an internal version number useful for third
    party scripts interfacing with it
22. VMD#8.7 is now found by the live demo scripts
23. live demos now run in vista 64-bit
24. init.create_random_polymers can now create polymers with more than
    one type of bond

Version 0.8.1 (2009-03-24)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1.  Significant performance enhancements
2.  New build option for compiling on UMich CAC clusters:
    ENABLE_CAC_GPU_ID compiles HOOMD to read in the *$CAC_GPU_ID*
    environment variable and use it to determine which GPUs to execute
    on. No –gpu command line required in job scripts any more.
3.  Particles can now be assigned a *non-unit mass*
4.  *init.reset()* command added to allow for the creation of a looped
    series of simulations all in python
5.  *dump.pdb()* command for writing PDB files
6.  pair.lj now comes with an option to *shift* the potential energy to
    0 at the cutoff
7.  pair.lj now comes with an opiton to *smoothly switch* both the
    *potential* and *force* to 0 at the cutoff with the XPLOR smoothing
    function
8.  *Gaussian pair potential* computation added (pair.gauss)
9.  update and analyze commands can now be given a function to determine
    a non-linear rate to run at
10. analyze.log, and dump.dcd can now append to existing files

*Changes that will break scripts from 0.8.0*

1. *dump.mol2()* has been changed to be more consistent with other dump
   commands. In order to get the same result as the previous behavior,
   replace

   .. code::

       dump.mol2(filename="file.mol2")

   with

   .. code::

       mol2 = dump.mol2()
       mol2.write(filename="file.mol2")

2. Grouping commands have been moved to their own package for
   organizational purposes. *group_all()* must now be called as
   *group.all()* and similarly for tags and type.

*Bug/fixes small enhancements*

1.  Documentation updates
2.  DCD file writing no longer crashes HOOMD in windows
3.  !FindBoost.cmake is patched upstream. Use CMake 2.6.3 if you need
    BOOST_ROOT to work correctly
4.  Validation tests now run with –gpu_error_checking
5.  ULF bug workarounds are now enabled only on hardware where they are
    needed. This boosts performance on C1060 and newer GPUs.
6.  !FindPythonLibs now always finds the shared python libraries, if
    they exist
7.  “make package” now works fine on mac os x
8.  Fixed erroneously reported dangerous neighbor list builds when using
    –mode=cpu
9.  Small tweaks to the XML file format.
10. Numerous performance enhancements
11. Workaround for ULF on compute#1 devices in place
12. dump.xml can now be given the option “all=true” to write all fields
13. total momentum can now be logged by analyze.log
14. HOOMD now compiles with boost#38 (and hopefully future versions)
15. Updaters can now be given floating point periods such as 1e5
16. Additional warnings are now printed when HOOMD is about to allocate
    a large amount of memory due to the specification of an extremely
    large box size
17. run() now shows up in the documentation index
18. Default sorter period is now 100 on CPUs to improve performance on
    chips with small caches

Version 0.8.0 (2008-12-22)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

*New features*

1. Addition of FENE bond potential
2. Addition of update.zero_momentum command to zero a system’s linear
   momentum
3. Brownian dynamics integration implemented
4. Multi-GPU simulations
5. Particle image flags are now tracked. analyze.msd command added to
   calculate the mean squared displacement.

*Changes that will break scripts from 0.7.x*

1. analyze.log quantity names have changed

*Bug/fixes small enhancements*

1.  Performance of the neighbor list has been increased significantly on
    the GPU (overall performance improvements are approximately 10%)
2.  Profile option added to the run() command
3.  Warnings are now correctly printed when negative coefficients are
    given to bond forces
4.  Simulations no longer fail on G200 cards
5.  Mac OS X binaries will be provided for download: new documentation
    for installing on Mac OS x has been written
6.  Two new demos showcasing large systems
7.  Particles leaving the simulation box due to bad initial conditions
    now generate an error
8.  win64 installers will no longer attempt to install on win32 and
    vice-versa
9.  neighborlist check_period now defaults to 1
10. The elapsed time counter in run() now continues counting time over
    multiple runs.
11. init.create_random_polymers now throws an error if the bond length
    is too small given the specified separation radii
12. Fixed a bug where a floating point value for the count field in
    init.create_random_polymers produced an error
13. Additional error checking to test if particles go NaN
14. Much improved status line printing for identifying hoomd_script
    commands
15. Numerous documentation updates
16. The VS redistributable package no longer needs to be installed to
    run HOOMD on windows (these files are distributed with HOOMD)
17. Now using new features in doxygen#5.7 to build pdf user
    documentation for download.
18. Performance enhancements of the Lennard-Jones pair force
    computation, thanks to David Tarjan
19. A header prefix can be added to log files to make them more gnuplot
    friendly
20. Log quantities completely revamped. Common quantities (i.e. kinetic
    energy, potential energy can now be logged in any simulation)
21. Particle groups can now be created. Currently only analyze.msd makes
    use of them.
22. The CUDA toolkit no longer needs to be installed to run a packaged
    HOOMD binary in windows.
23. User documentation can now be downloaded as a pdf.
24. Analyzers and updaters now count time 0 as being the time they were
    created, instead of time step 0.
25. Added job test scripts to aid in validating HOOMD
26. HOOMD will now build with default settings on a linux/unix-like OS
    where the boost static libraries are not installed, but the dynamic
    ones are.

Version 0.7.1 (2008-09-12)
^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Fixed bug where extremely large box dimensions resulted in an
   argument error - ticket:118
2. Fixed bug where simulations ran incorrectly with extremely small box
   dimensions - ticket:138

Version 0.7.0 (2008-08-12)
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Note:* only major changes are listed here.

1.  Stability and performance improvements.
2.  Cleaned up the hoomd_xml file format.
3.  Improved detection of errors in hoomd_xml files significantly.
4.  Users no longer need to manually specify HOOMD_ROOT, unless their
    installation is non-standard
5.  Particle charge can now be read in from a hoomd_xml file
6.  Consistency changes in the hoomd_xml file format: HOOMD 0.6.0 XML
    files are not compatible. No more compatibility breaking changes are
    planned after 0.7.0
7.  Enabled parallel builds in MSVC for faster compilation times on
    multicore systems
8.  Numerous small bug fixes
9.  New force compute for implementing walls
10. Documentation updates
11. Support for CUDA 2.0
12. Bug fixed allowing simulations with no integrator
13. Support for boost#35.0
14. Cleaned up GPU code interface
15. NVT integrator now uses tau (period) instead of Q (the mass of the
    extra degree of freedom).
16. Added option to NVE integration to limit the distance a particle
    moves in a single time step
17. Added code to dump system snapshots in the DCD file format
18. Particle types can be named by strings
19. A snapshot of the initial configuration can now be written in the
    .mol2 file format
20. The default build settings now enable most of the optional features
21. Separated the user and developer documentation
22. Mixed polymer systems can now be generated inside HOOMD
23. Support for CMake 2.6.0
24. Wrote the user documentation
25. GPU selection from the command line
26. Implementation of the job scripting system
27. GPU can now handle neighbor lists that overflow
28. Energies are now calculated
29. Added a logger for logging energies during a simulation run
30. Code now actually compiles on Mac OS X
31. Benchmark and demo scripts now use the new scripting system
32. Consistent error message format that is more visible.
33. Multiple types of bonds each with the own coefficients are now
    supported
34. Added python scripts to convert from HOOMD’s XML file format to
    LAMMPS input and dump files
35. Fixed a bug where empty xml nodes in input files resulted in an
    error message
36. Fixed a bug where HOOMD seg faulted when a particle left the
    simulation , vis=True)\* is now a convenient short hand for what was
    previously box now works fine on mac os x
37. Fixed erroneously reported dangerous neighbor list builds when using
    –mode=cpu
38. Small tweaks to the XML file format.
39. Numerous performance enhancements
40. Workaround for ULF on compute#1 devices in place
41. dump.xml can now be given the option
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Contributing
============

Contributions are welcomed via `pull requests on GitHub
<https://github.com/glotzerlab/hoomd-blue/pulls>`__. Contact the **HOOMD-blue** developers before
starting work to ensure it meshes well with the planned development direction and standards set for
the project.

Features
--------

Implement functionality in a general and flexible fashion
_________________________________________________________

New features should be applicable to a variety of use-cases. The **HOOMD-blue** developers can
assist you in designing flexible interfaces.

Maintain performance of existing code paths
___________________________________________

Expensive code paths should only execute when requested.

Optimize for the current GPU generation
_______________________________________

Write, test, and optimize your GPU kernels on the latest generation of GPUs.

Version control
---------------

Base your work off the correct branch
_____________________________________

- Base backwards compatible bug fixes on ``trunk-patch``.
- Base additional functionality on ``trunk-minor``.
- Base API incompatible changes on ``trunk-major``.

Propose a minimal set of related changes
________________________________________

All changes in a pull request should be closely related. Multiple change sets that are loosely
coupled should be proposed in separate pull requests.

Agree to the Contributor Agreement
__________________________________

All contributors must agree to the Contributor Agreement before their pull request can be merged.

Source code
-----------

Use a consistent style
______________________

The **Code style** section of the documentation sets the style guidelines for **HOOMD-blue** code.

Document code with comments
___________________________

Use doxygen header comments for classes, functions, etc. Also comment complex sections of code so
that other developers can understand them.

Compile without warnings
________________________

Your changes should compile without warnings.

Tests
-----

Write unit tests
________________

Add unit tests for all new functionality.

Validity tests
______________

The developer should run research-scale simulations using the new functionality and ensure that it
behaves as intended.

User documentation
------------------

Write user documentation
________________________

Document public-facing API with Python docstrings in Google style.

Document version status
_______________________

Add `versionadded, versionchanged, and deprecated Sphinx directives
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-versionadded>`__
to each user-facing Python class, method, etc., so that users will be aware of how functionality
changes from version to version. Remove this when breaking APIs in major releases.

Add developer to the credits
____________________________

Update the credits documentation to list the name and affiliation of each individual that has
contributed to the code.

Propose a change log entry
__________________________

Propose a short concise entry describing the change in the pull request description.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Installing binaries
===================

**HOOMD-blue** binaries are available in the glotzerlab-software_ Docker_/Singularity_ images and in
packages on conda-forge_

.. _glotzerlab-software: https://glotzerlab-software.readthedocs.io
.. _Docker: https://hub.docker.com/
.. _Singularity: https://www.sylabs.io/
.. _conda-forge: https://conda-forge.org/docs/user/introduction.html

Singularity / Docker images
---------------------------

See the glotzerlab-software_ documentation for instructions to install and use the containers on
supported HPC clusters.

Conda package
-------------

**HOOMD-blue** is available on conda-forge_ on the *linux-64*, *osx-64*, and *osx-arm64* platforms.
Install the ``hoomd`` package from the conda-forge_ channel into a conda environment::

    $ conda install -c conda-forge hoomd

Recent versions of ``conda`` auto-detect whether your system has a GPU and installs the appropriate
package. Override this and force GPU package installation with::

    $ conda install -c conda-forge hoomd=*=*gpu*

.. tip::

    Use miniforge_ or miniconda_ instead of the full Anaconda distribution to avoid package
    conflicts with conda-forge_ packages.

.. _miniforge: https://github.com/conda-forge/miniforge
.. _miniconda: http://conda.pydata.org/miniconda.html
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Citing HOOMD-blue
===================

**Please cite this publication in any work that uses HOOMD-blue:**
    J. A. Anderson, J. Glaser, and S. C. Glotzer. HOOMD-blue: A Python package
    for high-performance molecular dynamics and hard particle Monte Carlo
    simulations Computational Materials Science 173: 109363, Feb 2020.
    `10.1016/j.commatsci.2019.109363 <https://dx.doi.org/10.1016/j.commatsci.2019.109363>`_

----------

*The following publications document significant contributions to features in
HOOMD-blue. We encourage you to cite these, if possible, when you make use of
these specific functionalities.*

HPMC:
    J. A. Anderson, M. E. Irrgang, and S. C. Glotzer. Scalable Metropolis Monte
    Carlo for simulation of hard shapes Computer Physics Communications 204:
    21-30, July 2016. `10.1016/j.cpc.2016.02.024 <https://dx.doi.org/10.1016/j.cpc.2016.02.024>`_

Implicit depletants in HPMC:
    J. Glaser, A. S. Karas, and S. C. Glotzer. A parallel algorithm for implicit
    depletant simulations The Journal of Chemical Physics 143: 184110, 2015.
    `10.1063/1.4935175 <https://dx.doi.org/10.1063/1.4935175>`_

MPI scaling:
    J. Glaser, T. D. Nguyen, J. A. Anderson, P. Lui, F. Spiga, J. A. Millan,
    D. C. Morse, S. C. Glotzer. Strong scaling of general-purpose molecular
    dynamics simulations on GPUs Computer Physics Communications 192: 97-107,
    July 2015. `10.1016/j.cpc.2015.02.028 <https://dx.doi.org/10.1016/j.cpc.2015.02.028>`_

Intra-node scaling on multiple GPUs:
    J. Glaser, P. S. Schwendeman, J. A. Anderson, S. C. Glotzer. Unified memory
    in HOOMD-blue improves node-level strong scaling Computational Materials
    Science 173: 109359, Feb 2020.
    `10.1016/j.commatsci.2019.109359 <https://dx.doi.org/10.1016/j.commatsci.2019.109359>`_

When including historical development of HOOMD-blue, or noting that HOOMD-blue was first implemented on GPUs, please also cite:
    J. A. Anderson, C. D. Lorenz, and A. Travesset. General purpose molecular
    dynamics simulations fully implemented on graphics processing units Journal
    of Computational Physics 227(10): 5342-5359, May 2008.
    `10.1016/j.jcp.2008.01.047 <https://dx.doi.org/10.1016/j.jcp.2008.01.047>`_

DEM:
    M. Spellings, R. L. Marson, J. A. Anderson, and S. C. Glotzer. GPU
    accelerated Discrete Element Method (DEM) molecular dynamics for
    conservative, faceted particle simulations Journal of Computational Physics
    334: 460-467, Apr 2017. `10.1016/j.jcp.2017.01.014 <https://dx.doi.org/10.1016/j.jcp.2017.01.014>`_

The tree or stencil MD neighbor list:
    M. P. Howard, J. A. Anderson, A. Nikoubashman, S. C. Glotzer, and A. Z.
    Panagiotopoulos. Efficient neighbor list calculation for molecular
    simulation of colloidal systems using graphics processing units Computer
    Physics Communications 203: 45-52, Mar 2016.
    `10.1016/j.cpc.2016.02.003 <https://dx.doi.org/10.1016/j.cpc.2016.02.003>`_

    M. P. Howard, A. Statt, F. Madutsa, T. M. Truskett, and A. Z.
    Panagiotopoulos. Quantized bounding volume hierarchies for neighbor search
    in molecular simulations on graphics processing units Computational
    Materials Science 164(15): 139-146, June 2019.
    `10.1016/j.commatsci.2019.04.004 <https://dx.doi.org/10.1016/j.commatsci.2019.04.004>`_

MPCD:
    M. P. Howard, A. Z. Panagiotopoulos, and A. Nikoubashman. Efficient
    mesoscale hydrodynamics: Multiparticle collision dynamics with massively
    parallel GPU acceleration Computer Physics Communications 230: 10-20, Sep.
    2018.
    `10.1016/j.cpc.2018.04.009 <https://dx.doi.org/10.1016/j.cpc.2018.04.009>`_

Rigid bodies in MD:
    T. D. Nguyen, C. L. Phillips, J. A. Anderson, and S. C. Glotzer. Rigid body
    constraints realized in massively-parallel molecular dynamics on graphics
    processing units Computer Physics Communications 182(11): 2313-2307,
    June 2011.
    `10.1016/j.cpc.2011.06.005 <https://dx.doi.org/10.1016/j.cpc.2011.06.005>`_

    J. Glaser, X. Zha, J. A. Anderson, S. C. Glotzer, A. Travesset. Pressure in
    rigid body molecular dynamics, Computational Materials Science Computational
    Materials Science 173: 109430, Feb 2020.
    `10.1016/j.commatsci.2019.109430 <https://dx.doi.org/10.1016/j.commatsci.2019.109430>`_

DPD:
    C. L. Phillips, J. A. Anderson, and S. C. Glotzer. Pseudo-random number
    generation for Brownian Dynamics and Dissipative Particle Dynamics
    simulations on GPU devices Journal of Computational Physics 230(19):
    7191-7201, Aug. 2011.
    `10.1016/j.jcp.2011.05.021 <https://dx.doi.org/10.1016/j.jcp.2011.05.021>`_

EAM:
    I.V. Morozov, A.M. Kazennova, R.G. Bystryia, G.E. Normana, V.V. Pisareva,
    and V.V. Stegailova. Molecular dynamics simulations of the relaxation
    processes in the condensed matter on GPUs Computer Physics Communications
    182(9): 1974-1978, 2011.
    `10.1016/j.cpc.2010.12.026 <https://dx.doi.org/10.1016/j.cpc.2010.12.026>`_

    L. Yang, F. Zhang, C. Wang, K. Ho, and A. Travesset. Implementation of
    metal-friendly EAM/FS-type semi-empirical potentials in HOOMD-blue: A
    GPU-accelerated molecular dynamics software Journal of Computational
    Physics 359(15): 352-360, 2018.
    `10.1016/j.jcp.2018.01.015 <https://dx.doi.org/10.1016/j.jcp.2018.01.015>`_

PPPM:
    D. N. LeBard, B. G. Levine, P. Mertmann, S. A. Barr, A. Jusufi, S. Sanders,
    M. L. Klein, and A. Z. Panagiotopoulos. Self-assembly of coarse-grained
    ionic surfactants accelerated by graphics processing units Soft Matter 8:
    2385-2397, 2012.
    `10.1039/c1sm06787g <https://dx.doi.org/10.1039/c1sm06787g>`_

CGCMM potential:
    B. G. Levine, D. N. LeBard, R. DeVane, W. Shinoda, A. Kohlmeyer, and M. L.
    Klein. Micellization studied by GPU-accelerated coarse-grained molecular
    dynamics Journal of Chemical Theory and Computation 7(12): 4135-4145, Oct.
    2011. `10.1021/ct2005193 <https://dx.doi.org/10.1021/ct2005193>`_
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

.. include:: ../CHANGELOG.rst
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc
==========

.. rubric:: Details

.. automodule:: hoomd.hpmc
    :synopsis: HPMC package.
    :members:

.. rubric:: Modules

.. toctree::
    :maxdepth: 3

    module-hpmc-compute
    module-hpmc-integrate
    module-hpmc-nec
    module-hpmc-tune
    module-hpmc-update
    module-hpmc-pair
    module-hpmc-external
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.pair
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.pair

.. autosummary::
    :nosignatures:

    Buckingham
    DLVO
    DPD
    DPDLJ
    DPDConservative
    Ewald
    ExpandedLJ
    ExpandedMie
    ForceShiftedLJ
    Fourier
    Gauss
    LJ
    LJ1208
    LJ0804
    Mie
    Morse
    Moliere
    OPP
    Pair
    ReactionField
    Table
    TWF
    Yukawa
    ZBL

.. rubric:: Details

.. automodule:: hoomd.md.pair
    :synopsis: Pair potentials.
    :members: Pair,
        Buckingham,
        DLVO,
        DPD,
        DPDLJ,
        DPDConservative,
        Ewald,
        ExpandedMie,
        ForceShiftedLJ,
        Fourier,
        Gauss,
        LJ,
        LJ1208,
        LJ0804,
        Mie,
        Morse,
        Moliere,
        OPP,
        ReactionField,
        ExpandedLJ,
        Table,
        TWF,
        Yukawa,
        ZBL
    :show-inheritance:

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-md-pair-aniso
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.md
========

.. rubric:: Details

.. automodule:: hoomd.md
    :synopsis: Molecular Dynamics.
    :members: Integrator

.. rubric:: Modules

.. toctree::
    :maxdepth: 3

    module-md-angle
    module-md-bond
    module-md-constrain
    module-md-compute
    module-md-data
    module-md-dihedral
    module-md-external
    module-md-force
    module-md-improper
    module-md-long_range
    module-md-manifold
    module-md-many_body
    module-md-methods
    module-md-minimize
    module-md-nlist
    module-md-pair
    module-md-special_pair
    module-md-update
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.update
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.update

.. autosummary::
    :nosignatures:

    ActiveRotationalDiffusion
    ReversePerturbationFlow
    ZeroMomentum


.. rubric:: Details

.. automodule:: hoomd.md.update
    :synopsis: Updaters.
    :members: ActiveRotationalDiffusion,
              ReversePerturbationFlow,
              ZeroMomentum
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.data
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.data

.. autosummary::
    :nosignatures:

    LocalSnapshot
    LocalSnapshotGPU
    AngleLocalAccessBase
    BondLocalAccessBase
    ConstraintLocalAccessBase
    DihedralLocalAccessBase
    ImproperLocalAccessBase
    PairLocalAccessBase
    ParticleLocalAccessBase
    hoomd.data.typeparam.TypeParameter

.. rubric:: Details

.. automodule:: hoomd.data
    :synopsis: Provide access in Python to data buffers on CPU or GPU.
    :members: AngleLocalAccessBase,
              BondLocalAccessBase,
              ConstraintLocalAccessBase,
              DihedralLocalAccessBase,
              ImproperLocalAccessBase,
              PairLocalAccessBase,
              ParticleLocalAccessBase

    .. autoclass:: LocalSnapshot
        :inherited-members:

    .. autoclass:: LocalSnapshotGPU
        :inherited-members:

    .. autoclass:: hoomd.data.typeparam.TypeParameter
        :members: default

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-hoomd-array
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.tune
---------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.tune

.. autosummary::
    :nosignatures:

    BoxMCMoveSize
    MoveSize

.. rubric:: Details

.. automodule:: hoomd.hpmc.tune
    :synopsis: Tuners for HPMC.
    :members:

    .. autoclass:: BoxMCMoveSize(trigger, moves, target, solver, max_move_size=None)
        :members: secant_solver, scale_solver

        .. method:: tuned()
            :property:

            Whether or not the moves sizes have converged to the desired acceptance rate.

            :type: bool


    .. autoclass:: MoveSize(trigger, moves, target, solver, types=None, max_move_size=None)
        :members: secant_solver, scale_solver

        .. method:: tuned()
            :property:

            Whether or not the moves sizes have converged to the desired acceptance rate.

            :type: bool
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.filter
------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.filter

.. autosummary::
    :nosignatures:

    ParticleFilter
    All
    CustomFilter
    Intersection
    Null
    SetDifference
    Tags
    Type
    Union

.. rubric:: Details

.. automodule:: hoomd.filter
    :synopsis: Particle selection filters.
    :no-members:

    .. autoclass:: ParticleFilter()
        :special-members: __call__, __hash__, __eq__, __str__
    .. autoclass:: All()
    .. autoclass:: CustomFilter()
        :special-members: __call__
    .. autoclass:: Intersection(f, g)
    .. autoclass:: Null()
    .. autoclass:: SetDifference(f, g)
    .. autoclass:: Tags(tags)
        :members: tags
    .. autoclass:: Type(types)
        :members: types
    .. autoclass:: Union(f, g)
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.external.field
-------------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.external.field

.. autosummary::
    :nosignatures:

    ExternalField
    Harmonic

.. rubric:: Details

.. automodule:: hoomd.hpmc.external.field
    :synopsis: External fields.
    :members: ExternalField, Harmonic
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.constrain
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.constrain

.. autosummary::
    :nosignatures:

    Constraint
    Distance
    Rigid


.. rubric:: Details

.. automodule:: hoomd.md.constrain
    :synopsis: Constraints.
    :undoc-members:
    :members: Constraint, Distance, Rigid
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.external
--------------

.. py:currentmodule:: hoomd.md.external

.. automodule:: hoomd.md.external
    :synopsis: External potentials for MD.

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-md-external-field
   module-md-external-wall
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.bond
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.bond

.. autosummary::
    :nosignatures:

    Bond
    FENEWCA
    Harmonic
    Table
    Tether

.. rubric:: Details

.. automodule:: hoomd.md.bond
    :synopsis: Bond potentials.
    :members: Bond,
              FENEWCA,
              Harmonic,
              Table,
              Tether
    :no-inherited-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.trigger
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.trigger

.. autosummary::
    :nosignatures:

    After
    And
    Before
    Not
    On
    Or
    Periodic
    Trigger

.. rubric:: Details

.. automodule:: hoomd.trigger
    :synopsis: Trigger events at specific time steps.
    :no-members:

    .. autoclass:: After(timestep)
        :show-inheritance:
    .. autoclass:: And(triggers)
        :show-inheritance:
    .. autoclass:: Before(timestep)
        :show-inheritance:
    .. autoclass:: Not(trigger)
        :show-inheritance:
    .. autoclass:: On(timestep)
        :show-inheritance:
    .. autoclass:: Or(triggers)
        :show-inheritance:
    .. autoclass:: Periodic(period, phase)
        :show-inheritance:
    .. autoclass:: Trigger()
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.many_body
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.many_body

.. autosummary::
    :nosignatures:

    Triplet
    RevCross
    SquareDensity
    Tersoff

.. rubric:: Details

.. automodule:: hoomd.md.many_body
    :synopsis: Many-body potentials.
    :members: Triplet,
        RevCross,
        SquareDensity,
        Tersoff
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.update
-----------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.update

.. autosummary::
    :nosignatures:

    BoxMC
    Clusters
    MuVT
    QuickCompress

.. rubric:: Details

.. automodule:: hoomd.hpmc.update
    :synopsis: HPMC updaters.
    :members: BoxMC, Clusters, MuVT, QuickCompress
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Units
+++++

HOOMD-blue does not adopt a particular system of units, nor does it offer a variety of systems
to choose from. Instead, it follows a self-consistent system of units where all derived units
(e.g. force) are defined in terms of base units (e.g. energy / length). To adopt a system of units
for your simulations, choose a set of base units (e.g. meters versus centimeters for length), and
then determine what the derived units are.

Base Units
==========

The base units are:

- :math:`[\mathrm{energy}]`
- :math:`[\mathrm{length}]`
- :math:`[\mathrm{mass}]`

Unit Conversion
===============

Example unit conversions between derived units and base units:

.. list-table::
   :header-rows: 1

   * - Derived units
     - Relation to base units
   * - :math:`[\mathrm{area}]`
     - :math:`[\mathrm{length}]^2`
   * - :math:`[\mathrm{volume}]`
     - :math:`[\mathrm{length}]^3`
   * - :math:`[\mathrm{time}]`
     - :math:`[\mathrm{energy}]^{-1/2} \cdot [\mathrm{length}] \cdot [\mathrm{mass}]^{1/2}`
   * - :math:`[\mathrm{velocity}]`
     - :math:`[\mathrm{energy}]^{1/2} \cdot [\mathrm{mass}]^{-1/2}`
   * - :math:`[\mathrm{force}]`
     - :math:`[\mathrm{energy}] \cdot [\mathrm{length}]^{-1}`
   * - :math:`[\mathrm{pressure}]`
     - :math:`[\mathrm{energy}] \cdot [\mathrm{length}]^{-3}`
   * - :math:`[\mathrm{charge}]`
     - :math:`\left(4 \pi \epsilon_{0} \cdot [\mathrm{energy}] \cdot [\mathrm{length}] \right)^{1/2}`
       - where :math:`\epsilon_{0}` is permittivity of free space

.. note::

    Most of the units on this page apply to MD simulations.

    In HPMC, the primary unit is that of length. Mass is factored out of the partition function and
    does not enter into the simulation. In addition, the energy scale is irrelevant in athermal
    HPMC systems where overlapping energies are infinite and valid configurations have
    zero potential energy. However, energy does appear implicitly in derived units like
    :math:`[\mathrm{pressure}] = [\mathrm{energy}] \cdot [\mathrm{length}]^{-3}`.  In
    HPMC, :math:`kT` is set to 1 :math:`\mathrm{energy}`.

Common unit systems
===================

Example base and derived units for common MD unit systems.

.. note::

    All conversion factors given here are computed with Wolfram Alpha using the provided links.

.. list-table::
   :header-rows: 1

   * - Unit
     - AKMA
     - MD
   * - :math:`[\mathrm{energy}]`
     - kcal/mol
     - kJ/mol
   * - :math:`[\mathrm{length}]`
     - Å
     - nm
   * - :math:`[\mathrm{mass}]`
     - atomic mass unit
     - atomic mass unit
   * - :math:`[\mathrm{area}]`
     - :math:`\mathrm{Å}^2`
     - :math:`\mathrm{nm}^2`
   * - :math:`[\mathrm{volume}]`
     - :math:`\mathrm{Å}^3`
     - :math:`\mathrm{nm}^3`
   * - :math:`[\mathrm{time}]`
     - `48.8882129 fs <https://www.wolframalpha.com/input/?i=angstrom+*+amu%5E%281%2F2%29+*+%28kcal%2FAvogadro+number%29%5E%28%E2%88%921%2F2%29>`__
     - `1 ps <https://www.wolframalpha.com/input/?i=nanometer+*+amu%5E%281%2F2%29+*+%28kilojoule%2FAvogadro+number%29%5E%28%E2%88%921%2F2%29>`__
   * - :math:`[\mathrm{velocity}]`
     - `0.02045482828 Å/fs <https://www.wolframalpha.com/input/?i=%28kcal%2FAvogadro+number%29%5E%281%2F2%29+*+amu%5E%28-1%2F2%29+in+angstrom%2Ffs>`__
     - 1 nm/ps
   * - :math:`[\mathrm{force}]`
     - kcal/mol/Å
     - kJ/mol/nm
   * - :math:`[\mathrm{pressure}]`
     - `68568.4230 atm <https://www.wolframalpha.com/input/?i=%28kcal%2FAvogadro+number%29+*+angstrom%5E%28-3%29+in+atmospheres>`__
     - `16.3882464 atm <https://www.wolframalpha.com/input/?i=%28kilojoule%2FAvogadro+number%29+*+nanometer%5E%28-3%29+in+atmospheres>`__
   * - :math:`[\mathrm{charge}]`
     - `0.05487686461 e <https://www.wolframalpha.com/input/?i=sqrt%284+*+pi+*+permittivity+of+free+space+*+1+%28kcal%2FAvogadro%27s+number%29+*+1+angstrom%29+%2F+proton+charge>`__
     - `0.0848385920 e <https://www.wolframalpha.com/input/?i=sqrt%284+*+pi+*+permittivity+of+free+space+*+1+%28kilojoule%2FAvogadro%27s+number%29+*+1+nanometer%29+%2F+proton+charge>`__
   * - :math:`k` (Boltzmann's constant)
     - `0.00198720426 kcal/mol/K <https://www.wolframalpha.com/input/?i=boltzmann%27s+constant+in+kcal%2FAvogadro+number%2FK>`__
     - `0.00831446262 kJ/mol/K <https://www.wolframalpha.com/input/?i=boltzmann%27s+constant+in+kilojoues%2FAvogadro+number%2FK>`__
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.operation
---------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.operation

.. autosummary::
    :nosignatures:

    Compute
    Integrator
    Operation
    Tuner
    TriggeredOperation
    Updater
    Writer

.. rubric:: Details

.. automodule:: hoomd.operation
    :synopsis: Classes define the interfaces and types for HOOMD-blue operations.
    :members: Compute, Integrator, Tuner, TriggeredOperation, Updater, Writer
    :show-inheritance:

    .. autoclass:: Operation
        :inherited-members:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.external.wall
-----------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.external.wall

.. autosummary::
    :nosignatures:

    ForceShiftedLJ
    Gauss
    LJ
    Mie
    Morse
    Yukawa
    WallPotential

.. rubric:: Details

.. automodule:: hoomd.md.external.wall
    :synopsis: MD wall potentials.
    :members: ForceShiftedLJ,
        Gauss,
        LJ,
        Mie,
        Morse,
        Yukawa,
        WallPotential
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Components
==========

Extend **HOOMD-blue** with a **component** implemented in C++ for performance-critical tasks, such
as pair potential evaluation. A component provides a number of related functionalities. For example,
the :py:mod:`hoomd.hpmc` component enables hard particle Monte Carlo methods with **HOOMD-blue**.

Compile and install components **built-in** or as **external** components. The **HOOMD-blue** build
process compiles all core and built-in components together, requiring one only one set of configure,
make, and install commands. External components compile and link against a **HOOMD-blue**
installation from a separate build directory with their own set of configure, make, and install
commands. You may compile a component either way. When the end user is compiling **HOOMD-blue** and
components from source, built-in components compile and install everything at once which minimizes
chances for errors (e.g. building **HOOMD-blue** against python 3.6, but the component against
python 3.7). External components provide more flexibility for packaging purposes.

The **HOOMD-Blue** source provides an example component template in the ``example_plugin``
subdirectory. ``example_plugin`` demonstrates how to add a new ``update`` command with both CPU and
GPU implementations. Use this as a template when developing your component.

Built-in components
-------------------

You can fork **HOOMD-blue** and add your component directly, or you can create a separate source
repository for your component. Create a symbolic link to the component in the ``hoomd`` source
directory to compile it as a built-in component::

  $ ln -s <path-to-component>/<component> hoomd-blue/hoomd/<component>

.. note::

    Built-in components may be used directly from the build directory or installed.

External components
-------------------

To compile an external component, you must first install **HOOMD-blue**. Then, configure your component
with CMake and install it into the hoomd python library. Point ``CMAKE_PREFIX_PATH`` at your virtual
environment (if needed) so that cmake can find **HOOMD-blue**::

  $ cmake -B build/<component> -S <path-to-component>
  $ cmake --build build/<component>
  $ cmake --install build/<component>

The component build environment, including the compiler, CUDA, MPI, python, and other libraries,
must match exactly with those used to build **HOOMD-blue**.

.. note::

    External components must be installed before use.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.manifold
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.manifold

.. autosummary::
    :nosignatures:

    Manifold
    Cylinder
    Diamond
    Ellipsoid
    Gyroid
    Plane
    Primitive
    Sphere

.. rubric:: Details

.. automodule:: hoomd.md.manifold
    :synopsis: Manifold constraints.
    :members: Manifold,
              Cylinder,
              Diamond,
              Ellipsoid,
              Gyroid,
              Plane,
              Primitive,
              Sphere
    :no-inherited-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.external
-------------------

.. automodule:: hoomd.hpmc.external
    :synopsis: External potentials for HPMC.

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-hpmc-external-user
   module-hpmc-external-field
   module-hpmc-external-wall
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.integrate
--------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.integrate

.. autosummary::
    :nosignatures:

    ConvexPolygon
    ConvexPolyhedron
    ConvexSpheropolygon
    ConvexSpheropolyhedron
    ConvexSpheropolyhedronUnion
    Ellipsoid
    FacetedEllipsoid
    FacetedEllipsoidUnion
    HPMCIntegrator
    Polyhedron
    SimplePolygon
    Sphere
    SphereUnion
    Sphinx

.. rubric:: Details

.. automodule:: hoomd.hpmc.integrate
    :synopsis: HPMC integrators.

    .. autoclass:: ConvexPolygon
        :show-inheritance:
        :members:
    .. autoclass:: ConvexPolyhedron
        :show-inheritance:
        :members:
    .. autoclass:: ConvexSpheropolygon
        :show-inheritance:
        :members:
    .. autoclass:: ConvexSpheropolyhedron
        :show-inheritance:
        :members:
    .. autoclass:: ConvexSpheropolyhedronUnion
        :show-inheritance:
        :members:
    .. autoclass:: Ellipsoid
        :show-inheritance:
        :members:
    .. autoclass:: FacetedEllipsoid
        :show-inheritance:
        :members:
    .. autoclass:: FacetedEllipsoidUnion
        :show-inheritance:
        :members:
    .. autoclass:: HPMCIntegrator
        :inherited-members:
        :show-inheritance:
    .. autoclass:: Polyhedron
        :show-inheritance:
        :members:
    .. autoclass:: SimplePolygon
        :show-inheritance:
        :members:
    .. autoclass:: Sphere
        :show-inheritance:
        :members:
    .. autoclass:: SphereUnion
        :show-inheritance:
        :members:
    .. autoclass:: Sphinx
        :show-inheritance:
        :members:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.methods.rattle
-----------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.methods.rattle

.. autosummary::
    :nosignatures:

    MethodRATTLE
    Brownian
    Langevin
    NVE
    OverdampedViscous

.. rubric:: Details

.. automodule:: hoomd.md.methods.rattle
    :synopsis: RATTLE integration methods.
    :members: MethodRATTLE,
        Brownian,
        Langevin,
        NVE,
        OverdampedViscous
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Deprecated
==========

Features deprecated in v3.x may be removed in a future v4.0.0 release.

v3.x
----

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
     - Deprecated in
   * - Particle diameters
     - Potentials such as `md.pair.ExpandedLJ`.
     - v3.0.0
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

.. include:: ../BUILDING.rst
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Migrating to HOOMD v3
=====================

HOOMD v3 introduces many breaking changes for both users and developers
in order to provide a cleaner Python interface, enable new functionalities, and
move away from unsupported tools. This guide highlights those changes.

Overview of API changes
-----------------------

HOOMD v3 introduces a completely new API. All classes have been renamed to match
PEP8 naming guidelines and have new or renamed parameters, methods, and
properties. See the tutorials and the Python module documentation for full
class-level details.

Here is a module level overview of features that have been moved or removed:

.. list-table::
   :header-rows: 1

   * - v2 module, class, or method
     - Replaced with
   * - ``hoomd.analyze.log``
     - `hoomd.logging`
   * - ``hoomd.benchmark``
     - *Removed.* Use Python standard libraries for timing.
   * - ``hoomd.cite``
     - *Removed.* See `citing`.
   * - ``hoomd.dump``
     - `hoomd.write`
   * - ``hoomd.compute.thermo``
     - `hoomd.md.compute.ThermodynamicQuantities`
   * - ``hoomd.context.initialize``
     - `hoomd.device.CPU` and `hoomd.device.GPU`
   * - ``hoomd.data``
     - `hoomd.State`
   * - ``hoomd.group``
     - `hoomd.filter`
   * - ``hoomd.init``
     - `hoomd.Simulation` ``create_state_from_`` factory methods
   * - ``hoomd.lattice``
     - *Removed.* Use an external tool.
   * - ``hoomd.meta``
     - `hoomd.logging.Logger`.
   * - ``hoomd.option``
     - *Removed.* Use Python standard libraries for option parsing.
   * - ``hoomd.update``
     - Some classes have been moved to `hoomd.tune`.
   * - ``hoomd.util``
     -  Enable GPU profiling with `hoomd.device.GPU.enable_profiling`.
   * - ``hoomd.hpmc.analyze.sdf``
     - `hoomd.hpmc.compute.SDF`
   * - ``hoomd.hpmc.data``
     - `hoomd.hpmc.integrate.HPMCIntegrator` properties.
   * - ``hoomd.hpmc.util``
     - `hoomd.hpmc.tune`
   * - ``hoomd.md.integrate.mode_standard``
     - `hoomd.md.Integrator`
   * - ``hoomd.md.update.rescale_temp``
     - `hoomd.State.thermalize_particle_momenta`
   * - ``hoomd.md.update.enforce2d``
     - *Removed.* This is not needed.
   * - ``hoomd.md.constrain.sphere``
     - `hoomd.md.manifold.Sphere`
   * - ``hoomd.md.constrain.oneD``
     - *Removed.*
   * - ``hoomd.md.update.constraint_ellipsoid``
     - `hoomd.md.manifold.Ellipsoid`
   * - ``hoomd.jit.patch``
     - `hoomd.hpmc.pair.user`
   * - ``hoomd.jit.external``
     - `hoomd.hpmc.external.user`

Removed functionality
---------------------

HOOMD v3 removes old APIs, unused functionality, and features better served by other codes:

:py:mod:`hoomd`:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - Python 2.7
     - Python >= 3.6
   * - Compute < 6.0 GPUs
     - Compute >= 6.0 GPUs
   * - ``static`` parameter in ``hoomd.dump.gsd``
     - ``dynamic`` parameter
   * - ``set_params`` and other ``set_*`` methods
     - Parameters and type parameters accessed by properties.
   * - ``context.initialize``
     - `device.CPU` / `device.GPU`
   * - ``util.quiet_status`` and ``util.unquiet_status``
     - No longer needed.

``hoomd.deprecated``:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - ``deprecated.analyze.msd``
     - Offline analysis: e.g. `Freud's msd module <https://freud.readthedocs.io>`_.
   * - ``deprecated.dump.xml``
     - `hoomd.write.GSD`
   * - ``deprecated.dump.pos``
     - `hoomd.write.GSD` with on-demand conversion to ``.pos``.
   * - ``deprecated.init.read_xml``
     - `Simulation.create_state_from_gsd`
   * - ``deprecated.init.create_random``
     - `mBuild <https://mosdef-hub.github.io/mbuild/>`_, `packmol <https://www.ime.unicamp.br/~martinez/packmol/userguide.shtml>`_, or user script.
   * - ``deprecated.init.create_random_polymers``
     - `mBuild <https://mosdef-hub.github.io/mbuild/>`_, `packmol <https://www.ime.unicamp.br/~martinez/packmol/userguide.shtml>`_, or user script.

:py:mod:`hoomd.hpmc`:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - ``sphere_union::max_members`` parameter
     - no longer needed
   * - ``convex_polyhedron_union``
     - :py:class:`ConvexSpheropolyhedronUnion <hoomd.hpmc.integrate.ConvexSpheropolyhedronUnion>`, ``sweep_radius=0``
   * - ``setup_pos_writer`` member
     - n/a
   * - ``depletant_mode='circumsphere'``
     - no longer needed
   * - ``max_verts`` parameter
     - no longer needed
   * - ``depletant_mode`` parameter
     - no longer needed
   * - ``ntrial`` parameter
     - no longer needed
   * - ``implicit`` boolean parameter
     - set ``fugacity`` non-zero

:py:mod:`hoomd.md`:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - ``group`` parameter to ``integrate.mode_minimize_fire``
     - Pass group to integration method.
   * - ``alpha`` parameter to ``pair.lj`` and related classes
     - n/a
   * - ``f_list`` and ``t_list`` parameters to ``md.force.active``
     - Per-type ``active_force`` and ``active_torque``
   * - ``md.pair.SLJ``
     - `md.pair.ExpandedLJ`

``hoomd.cgcmm``:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - ``cgcmm.angle.cgcmm``
     - no longer needed
   * - ``cgcmm.pair.cgcmm``
     - no longer needed

``hoomd.dem``:

.. list-table::
   :header-rows: 1

   * - Feature
     - Replace with
   * - DEM pair potentials
     - ALJ pair potential in `hoomd.md.pair.aniso`.

Not yet ported
--------------

The following v2 functionalities have not yet been ported to the v3 API. They may be added in a
future 3.x release:

- HPMC box volume move size tuner.

These contributed functionalities rely on the community for support. Please
contact the developers if you have an interest in porting these in a future release:

- ``hoomd.hdf5``
- ``hoomd.metal``
- ``hoomd.mpcd``


Compiling
---------

* CMake 3.8 or newer is required to build HOOMD v3.0.
* To compile with GPU support, use the option ``ENABLE_GPU=ON``.
* ``UPDATE_SUBMODULES`` no longer exists. Users and developers should use
  ``git clone --recursive``, ``git submodule update`` and ``git submodule sync``
  as appropriate.
* ``COPY_HEADERS`` no longer exists. HOOMD will pull headers from the source directory when needed.
* ``CMAKE_INSTALL_PREFIX`` is set to the Python ``site-packages`` directory (if
  not explicitly set by the user).
* **cereal**, **eigen**, and **pybind11** headers must be provided to build
  HOOMD. See :doc:`installation` for details.
* ``BUILD_JIT`` is replaced with ``ENABLE_LLVM``.

Components
----------

* HOOMD now uses native CUDA support in CMake. Use ``CMAKE_CUDA_COMPILER`` to
  specify a specific ``nvcc`` or ``hipcc``. Plugins will require updates to
  ``CMakeLists.txt`` to compile ``.cu`` files.

  - Remove ``CUDA_COMPILE``.
  - Pass ``.cu`` sources directly to ``pybind11_add_module``.
  - Add ``NVCC`` as a compile definition to ``.cu`` sources.

* External components require additional updates to work with v3. See
  ``example_plugin`` for details:

  - Remove ``FindHOOMD.cmake``.
  - Replace ``include(FindHOOMD.cmake)`` with
    ``find_package(HOOMD 3.Y REQUIRED)`` (where 3.Y is the minor version this
    plugin is compatible with).
  - Always force set ``CMAKE_INSTALL_PREFIX`` to ``${HOOMD_INSTALL_PREFIX}``.
  - Replace ``PYTHON_MODULE_BASE_DIR`` with ``PYTHON_SITE_INSTALL_DIR``.
  - Replace all ``target_link_libraries`` and ``set_target_properties`` with
    ``target_link_libraries(_${COMPONENT_NAME} PUBLIC HOOMD::_hoomd)`` (can link
    ``HOOMD::_md``, ``HOOMD::_hpmc``, etc. if necessary).

* Numerous C++ class APIs have changed, been removed, or renamed. Review the
  header files to see new class signatures. These changes may require you to
  update your component accordingly. Some of the more notable changes include:

  - ``Variant`` has been completely rewritten.
  - ``Trigger`` replaces periodic and variable period scheduling.
  - ``NeighborList`` has a ``addRCutMatrix`` method clients must use to specify
    the maximum cutoff radii per type pair.
  - ``timestep`` is now of type ``uint64_t``.
  - ``Saru`` has been removed. Use ``RandomGenerator``.
  - ``RandomGenerator`` is now constructed with a ``Seed`` and ``Counter``
    object that support 64-bit timesteps.
  - ``m_seed`` is no longer present in individual operation objects. Use the
    global seed provided by ``SystemDefinition``.
  - The HPMC integrators have been heavily refactored.
  - HPMC GPU kernels are now instantiated by template .cu files that are generated by CMake at
    configure time.
  - ``ParticleGroup`` instances are now constructed from immutable, reusable,
    and user-customizable ``ParticleFilter`` instances.
  - All GPU code is now written with HIP to support NVIDIA and AMD GPUs.
  - ``ActiveForceCompute`` always uses particle orientation in combination with
    per-type active forces and torques.
  - ``getProvidedLogQuantities`` and ``getLogQuantities`` have been removed. Provide loggable
    properties instead.
  - Removed the Sphere, Ellipsoid, and oneD constraints. Replaced with the more general RATTLE
    integration methods and Manifold classes.
  - Removed the Enforce2D and TempRescale Updaters. Enforce2D is not needed for 2D simulations,
    and TempRescale has been replaced by ``thermalize_`` methods.
  - Removed Doxygen configuration scripts. View the document for classes in the source files.
  - Particle types may no longer be added after a Simulation is initialized. Classes no longer
    need to subscribe to the types added signal and reallocate data structures when the number of
    types changes.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.dihedral
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.dihedral

.. autosummary::
    :nosignatures:

    Dihedral
    Harmonic
    OPLS
    Table

.. rubric:: Details

.. automodule:: hoomd.md.dihedral
    :synopsis: Dihedral potentials.
    :show-inheritance:
    :members: Dihedral,
              Harmonic,
              OPLS,
              Table
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.nec.tune
-------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.nec.tune

.. autosummary::
    :nosignatures:

    ChainTime

.. rubric:: Details

.. automodule:: hoomd.hpmc.nec.tune
    :synopsis: Tune Newtonian event chain parameters.
    :members: ChainTime
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.tune
-------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.tune

.. autosummary::
    :nosignatures:

    CustomTuner
    LoadBalancer
    ManualTuneDefinition
    ParticleSorter
    ScaleSolver
    SecantSolver
    SolverStep

.. rubric:: Details

.. automodule:: hoomd.tune
    :synopsis: Tuner simulation hyperparameters.
    :members: CustomTuner,
              LoadBalancer,
              ParticleSorter,
              ScaleSolver,
              SecantSolver,
              SolverStep
    :show-inheritance:

    .. autoclass:: ManualTuneDefinition
        :inherited-members:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.device
------------

.. py:currentmodule:: hoomd.device

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    CPU
    Device
    GPU
    auto_select

.. rubric:: Details

.. automodule:: hoomd.device
    :synopsis: Devices used for simulation runs
    :members: Device, CPU, GPU, auto_select
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

==========
HOOMD-blue
==========

.. only:: html

    |Citing-HOOMD|
    |conda-forge|
    |conda-forge-Downloads|
    |GitHub Actions|
    |Contributors|
    |License|


    .. |Citing-HOOMD| image:: https://img.shields.io/badge/cite-hoomd-blue.svg
        :target: https://hoomd-blue.readthedocs.io/en/latest/citing.html
    .. |conda-forge| image:: https://img.shields.io/conda/vn/conda-forge/hoomd.svg?style=flat
        :target: https://anaconda.org/conda-forge/hoomd
    .. |conda-forge-Downloads| image:: https://img.shields.io/conda/dn/conda-forge/hoomd.svg?style=flat
        :target: https://anaconda.org/conda-forge/hoomd
    .. |GitHub Actions| image:: https://github.com/glotzerlab/hoomd-blue/actions/workflows/test.yml/badge.svg?branch=trunk-patch
        :target: https://github.com/glotzerlab/hoomd-blue/actions/workflows/test.yml
    .. |Contributors| image:: https://img.shields.io/github/contributors-anon/glotzerlab/hoomd-blue.svg?style=flat
        :target: https://hoomd-blue.readthedocs.io/en/latest/credits.html
    .. |License| image:: https://img.shields.io/badge/license-BSD--3--Clause-green.svg
        :target: https://hoomd-blue.readthedocs.io/en/latest/license.html

**HOOMD-blue** is a Python package that runs simulations of particle systems on CPUs and GPUs. It
performs hard particle Monte Carlo simulations of a variety of shape classes and molecular dynamics
simulations of particles with a range of pair, bond, angle, and other potentials. Many features are
targeted at the soft matter research community, though the code is general and capable of many
types of particle simulations.

Resources
=========

- `GitHub Repository <https://github.com/glotzerlab/hoomd-blue>`_:
  Source code and issue tracker.
- :doc:`Citing HOOMD-blue </citing>`:
  How to cite the code.
- :doc:`Installation guide </installation>`:
  Instructions for installing **HOOMD-blue** binaries.
- :doc:`Compilation guide </building>`:
  Instructions for compiling **HOOMD-blue**.
- `hoomd-users mailing list <https://groups.google.com/d/forum/hoomd-users>`_:
  Send messages to the **HOOMD-blue** user community.
- `HOOMD-blue website <https://glotzerlab.engin.umich.edu/hoomd-blue/>`_:
  Additional information and publications.
- `HOOMD-blue benchmark scripts <https://github.com/glotzerlab/hoomd-benchmarks>`_:
  Scripts to evaluate the performance of HOOMD-blue simulations.

Related tools
=============

- `freud <https://freud.readthedocs.io/>`_:
  Analyze HOOMD-blue simulation results with the **freud** Python library.
- `signac <https://signac.io/>`_:
  Manage your workflow with **signac**.
- `Molecular Simulation Design Framework (MoSDeF)`_ tools:

  - `mbuild`_: Assemble reusable components into complex molecular systems.
  - `foyer`_: perform atom-typing and define classical molecular modeling force fields.

.. _Molecular Simulation Design Framework (MoSDeF): https://mosdef.org/
.. _mbuild: https://mbuild.mosdef.org/
.. _foyer: https://foyer.mosdef.org/

Example scripts
===============

These examples demonstrate some of the Python API.

Hard particle Monte Carlo:

.. code:: python

    import hoomd

    mc = hoomd.hpmc.integrate.ConvexPolyhedron()
    mc.shape['octahedron'] = dict(vertices=[
        (-0.5, 0, 0),
        (0.5, 0, 0),
        (0, -0.5, 0),
        (0, 0.5, 0),
        (0, 0, -0.5),
        (0, 0, 0.5),
    ])

    cpu = hoomd.device.CPU()
    sim = hoomd.Simulation(device=cpu, seed=20)
    sim.operations.integrator = mc
    # See HOOMD tutorial for how to construct an initial configuration 'init.gsd'
    sim.create_state_from_gsd(filename='init.gsd')

    sim.run(1e5)

Molecular dynamics:

.. code:: python

    import hoomd

    cell = hoomd.md.nlist.Cell()
    lj = hoomd.md.pair.LJ(nlist=cell)
    lj.params[('A', 'A')] = dict(epsilon=1, sigma=1)
    lj.r_cut[('A', 'A')] = 2.5

    integrator = hoomd.md.Integrator(dt=0.005)
    integrator.forces.append(lj)
    nvt = hoomd.md.methods.NVT(kT=1.5, filter=hoomd.filter.All(), tau=1.0)
    integrator.methods.append(nvt)

    gpu = hoomd.device.GPU()
    sim = hoomd.Simulation(device=gpu)
    sim.operations.integrator = integrator
    # See HOOMD tutorial for how to construct an initial configuration 'init.gsd'
    sim.create_state_from_gsd(filename='init.gsd')
    sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=1.5)

    sim.run(1e5)

.. toctree::
    :maxdepth: 1
    :caption: Getting started

    features
    installation
    building
    migrating
    changelog
    citing

.. toctree::
    :maxdepth: 1
    :caption: Tutorials

    tutorial/00-Introducing-HOOMD-blue/00-index
    tutorial/01-Introducing-Molecular-Dynamics/00-index
    tutorial/02-Logging/00-index
    tutorial/03-Parallel-Simulations-With-MPI/00-index
    tutorial/04-Custom-Actions-In-Python/00-index
    tutorial/05-Organizing-and-Executing-Simulations/00-index

.. toctree::
    :maxdepth: 1
    :caption: How to guides

    howto/molecular

.. toctree::
   :maxdepth: 3
   :caption: Python API

   package-hoomd
   package-hpmc
   package-md

.. toctree::
    :maxdepth: 1
    :caption: Developer guide

    contributing
    style
    testing
    components

.. toctree::
   :maxdepth: 1
   :caption: Reference

   notation
   units
   deprecated
   license
   credits
   indices
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.pair.user
--------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.pair.user

.. autosummary::
    :nosignatures:

    CPPPotentialBase
    CPPPotential
    CPPPotentialUnion

.. rubric:: Details

.. automodule:: hoomd.hpmc.pair.user
    :synopsis: User defined pair potentials for Monte Carlo.

    .. autoclass:: CPPPotentialBase
    .. autoclass:: CPPPotential
        :show-inheritance:
    .. autoclass:: CPPPotentialUnion
        :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Notation
==========

The HOOMD-blue documentation uses the following mathematical notation.

General notation:

.. list-table::

    * - :math:`x`
      - Scalar.
    * - :math:`\vec{a}`
      - Vector.
    * - :math:`a`
      - Magnitude of the vector :math:`\vec{a}`.
    * - :math:`\hat{a}`
      - Unit vector in the direction of :math:`\vec{a}`.
    * - :math:`\mathbf{A}`
      - Matrix.
    * - :math:`\mathbf{b}`
      - Quaternion.
    * - :math:`\vert \mathbf{b} \vert`
      - Magnitude of a quaternion.
    * - :math:`i`, :math:`j`, :math:`k`
      - Indices.
    * - :math:`i`
      - :math:`\sqrt{-1}` in some contexts.
    * - :math:`\mathrm{function}`
      - A mathematical function.
    * - a :math:`[\mathrm{length}]^2`
      - The quantity "a" has dimensions of :math:`\mathrm{length}^2`.

Symbol definitions:

.. list-table::

    * - :math:`\vec{v}_{ij}`
      - :math:`\vec{v}_j - \vec{v}_j`
    * - :math:`\vec{r}`
      - Position.
    * - :math:`\vec{v}`
      - Velocity.
    * - :math:`\mathbf{q}`
      - Orientation.
    * - :math:`m`
      - Mass.
    * - :math:`\vec{p}`
      - Momentum.
    * - :math:`I`
      - Moment of inertia.
    * - :math:`q`
      - Charge.
    * - :math:`N`
      - Number.
    * - :math:`L`
      - Length.
    * - :math:`V`
      - Volume.
    * - :math:`P`
      - Pressure.
    * - :math:`\rho`
      - Number density.
    * - :math:`\vec{a}_1`, :math:`\vec{a}_2`, :math:`\vec{a}_3`
      - Box unit cell vectors.
    * - :math:`\vec{F}`
      - Force.
    * - :math:`\vec{\tau}`
      - Torque.
    * - :math:`U`
      - Potential energy.
    * - :math:`K`
      - Kinetic energy.
    * - :math:`E`
      - Total energy.
    * - :math:`T`
      - Temperature.
    * - :math:`k`
      - Boltzmann's constant in some contexts. Typically appears multiplying temperature:
        :math:`kT`.
    * - :math:`k`
      - Spring constant in some contexts.
    * - :math:`\beta`
      - :math:`\frac{1}{kT}`
    * - :math:`t`
      - Time.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Index
=====

* :ref:`genindex`
* :ref:`modindex`
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.external.user
------------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.external.user

.. autosummary::
    :nosignatures:

    CPPExternalPotential

.. rubric:: Details

.. automodule:: hoomd.hpmc.external.user
    :synopsis: User defined fields for Monte Carlo.
    :members: CPPExternalPotential
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Code style
==========

All code in HOOMD-blue follows a consistent style to ensure readability. We
provide configuration files for linters (specified below) so that developers can
automatically validate and format files.

These tools are configured for use with `pre-commit`_ in
``.pre-commit-config.yaml``. You can install pre-commit hooks to validate your
code. Checks will run on pull requests. Run checks manually with::

    pre-commit run --all-files

.. _pre-commit: https://pre-commit.com/

Python
------

Python code in HOOMD-blue should follow `PEP8`_ with the formatting performed by
`yapf`_ (configuration in ``setup.cfg``). Code should pass all **flake8** tests
and formatted by **yapf**.

.. _PEP8: https://www.python.org/dev/peps/pep-0008
.. _yapf: https://github.com/google/yapf

Tools
^^^^^

* Linter: `flake8 <http://flake8.pycqa.org/en/latest/>`_

  * With these plugins:

    * `pep8-naming <https://github.com/PyCQA/pep8-naming>`_
    * `flake8-docstrings <https://gitlab.com/pycqa/flake8-docstrings>`_
    * `flake8-rst-docstrings <https://github.com/peterjc/flake8-rst-docstrings>`_

  * Configure flake8 in your editor to see violations on save.

* Autoformatter: `yapf <https://github.com/google/yapf>`_

  * Run: ``pre-commit run --all-files`` to apply style changes to the whole
    repository.

Documentation
^^^^^^^^^^^^^

Python code should be documented with docstrings and added to the Sphinx
documentation index in ``doc/``. Docstrings should follow `Google style`_
formatting for use in `Napoleon`_.

.. _Google Style: https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html#example-google
.. _Napoleon: https://www.sphinx-doc.org/en/master/usage/extensions/napoleon.html

Simulation operations should unambiguously document what calculations they perform using formal
mathematical notation and use a consistent set of symbols and across the whole codebase.
HOOMD-blue documentation should follow standard physics and statistical mechanics notation with
consistent use of symbols detailed in `notation`.

When referencing classes, methods, and properties in documentation, use ``name`` to refer to names
in the local scope (class method or property, or classes in the same module). For classes outside
the module, use the fully qualified name (e.g. ``numpy.ndarray`` or
``hoomd.md.compute.ThermodynamicQuantities``).

C++/CUDA
--------

* Style is set by **clang-format**

  * Whitesmith's indentation style.
  * 100 character line width.
  * Indent only with spaces.
  * 4 spaces per indent level.
  * See :file:`.clang-format` for the full **clang-format** configuration.

* Naming conventions:

  * Namespaces: All lowercase ``somenamespace``
  * Class names: ``UpperCamelCase``
  * Methods: ``lowerCamelCase``
  * Member variables: ``m_`` prefix followed by lowercase with words
    separated by underscores ``m_member_variable``
  * Constants: all upper-case with words separated by underscores
    ``SOME_CONSTANT``
  * Functions: ``lowerCamelCase``

Tools
^^^^^

* Autoformatter: `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_.

Documentation
^^^^^^^^^^^^^

Documentation comments should be in Javadoc format and precede the item they document for
compatibility with many source code editors. Multi-line documentation comment blocks start with
``/**`` and single line ones start with
``///``.

.. code:: c++

    /** Describe a class
     *
     *  Note the second * above makes this a documentation comment. Some
     *  editors like to add the additional *'s on each line. These may be
     * omitted
    */
    class SomeClass
        {
        public:
            /// Single line doc comments should have three /'s
            Trigger() { }

            /** This is a brief description of a method

                This is a longer description.

                @param arg This is an argument.
                @returns This describes the return value
            */
            virtual bool method(int arg)
                {
                return false;
                }
        private:

            /// This is a member variable
            int m_var;
        };

See ``Trigger.h`` for a good example.

Other file types
----------------

Use your best judgment and follow existing patterns when styling CMake,
restructured text, markdown, and other files. The following general guidelines
apply:

* 100 character line width.
* 4 spaces per indent level.
* 4 space indent.

Editor configuration
--------------------

`Visual Studio Code <https://code.visualstudio.com/>`_ users: Open the provided
workspace file (``hoomd.code-workspace``) which provides configuration
settings for these style guidelines.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.external.field
-----------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.external.field

.. autosummary::
    :nosignatures:

    Field
    Electric
    Periodic

.. rubric:: Details

.. automodule:: hoomd.md.external.field
    :synopsis: External field potentials.
    :members: Field,
        Electric,
        Periodic
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.communicator
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.communicator

.. autosummary::
    :nosignatures:

    Communicator

.. rubric:: Details

.. automodule:: hoomd.communicator
    :synopsis: MPI run interface.
    :members: Communicator
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.custom
------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.custom

.. autosummary::
    :nosignatures:

    Action
    CustomOperation

.. rubric:: Details

.. automodule:: hoomd.custom
    :synopsis: Classes for custom Python actions that allow injecting code into the run loop.
    :members: Action, CustomOperation
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.pair
---------------

.. automodule:: hoomd.hpmc.pair
    :synopsis: Pair potentials for HPMC.

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-hpmc-pair-user
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.long_range
-------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.long_range

.. rubric:: Details

.. automodule:: hoomd.md.long_range
    :synopsis: Long-range potentials for molecular dynamics.

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-md-long_range-pppm
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.compute
----------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.compute

.. autosummary::
    :nosignatures:

    HarmonicAveragedThermodynamicQuantities
    ThermodynamicQuantities

.. rubric:: Details

.. automodule:: hoomd.md.compute
    :synopsis: Compute system properties.
    :members: HarmonicAveragedThermodynamicQuantities, ThermodynamicQuantities
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd
=====

.. rubric:: Overview

.. py:currentmodule:: hoomd

.. autosummary::
    :nosignatures:

    Box
    Operations
    Simulation
    Snapshot
    State

.. rubric:: Details

.. automodule:: hoomd
    :synopsis: HOOMD-blue main package.
    :undoc-members:
    :imported-members:
    :members: Simulation,
              State,
              Snapshot,
              Operations,
              Box

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-hoomd-communicator
   module-hoomd-custom
   module-hoomd-data
   module-hoomd-device
   module-hoomd-error
   module-hoomd-filter
   module-hoomd-logging
   module-hoomd-operation
   module-hoomd-triggers
   module-hoomd-tune
   module-hoomd-update
   module-hoomd-variant
   module-hoomd-version
   module-hoomd-wall
   module-hoomd-write
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.version
-------------

.. py:currentmodule:: hoomd.version

.. automodule:: hoomd.version
    :synopsis: Version and build information.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.external.wall
------------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.external.wall

.. autosummary::
    :nosignatures:

    WallPotential

.. rubric:: Details

.. automodule:: hoomd.hpmc.external.wall
    :synopsis: HPMC walls.
    :members: WallPotential
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.data.array
----------------

.. py:currentmodule:: hoomd.data.array

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    HOOMDArray
    HOOMDGPUArray

.. rubric:: Details

.. automodule:: hoomd.data.array
    :synopsis: Array classes that expose internal data buffers through the local snapshot `hoomd.State` API.
    :members: HOOMDGPUArray

    .. autoclass:: HOOMDArray
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.nlist
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.nlist

.. autosummary::
    :nosignatures:

    NeighborList
    Cell
    Stencil
    Tree

.. rubric:: Details

.. automodule:: hoomd.md.nlist
    :synopsis: Neighbor list acceleration structures.
    :members: Cell, Stencil, Tree
    :no-inherited-members:
    :show-inheritance:

    .. autoclass:: NeighborList
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.pair.aniso
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.pair.aniso

.. autosummary::
    :nosignatures:

    ALJ
    AnisotropicPair
    Dipole
    GayBerne

.. rubric:: Details

.. automodule:: hoomd.md.pair.aniso
    :synopsis: Anisotropic pair potentials.
    :members:
        ALJ,
        AnisotropicPair,
        Dipole,
        GayBerne,
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.compute
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.compute

.. autosummary::
    :nosignatures:

    FreeVolume
    SDF

.. rubric:: Details

.. automodule:: hoomd.hpmc.compute
    :synopsis: Compute system properties.
    :members: FreeVolume, SDF
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.nec.integrate
------------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.hpmc.nec.integrate

.. autosummary::
    :nosignatures:

    ConvexPolyhedron
    HPMCNECIntegrator
    Sphere

.. rubric:: Details

.. automodule:: hoomd.hpmc.nec.integrate
    :synopsis: Integrators for Newtonian event chain Monte Carlo for hard particles.
    :members: ConvexPolyhedron,
              HPMCNECIntegrator,
              Sphere
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Testing
=======

All code in HOOMD must be tested to ensure that it operates correctly.

**Unit tests** check that basic functionality works, one class at a time. Unit
tests assume internal knowledge about how classes work and may use unpublished
APIs to stress test all possible input and outputs of a given class in order to
exercise all code paths. For example, test that the box class properly wraps
vectors back into the minimum image. Unit tests should complete in a fraction of
a second.

**System integration tests** check that many classes work together to produce
correct output. These tests are black box tests and should only use user-facing
APIs to provide inputs and check for correct outputs. For example, test that the
hard sphere HPMC simulation executes for several steps. System integration tests
may take several seconds.

**Validation tests** rigorously check that HOOMD simulations sample the correct
statistical ensembles. For example, validate the a Lennard-Jones simulation at a
given density matches the pressure in the NIST reference. Validation tests
should run long enough to ensure reasonable sampling, but not too long. These
test run in a CI environment on every pull request. Individual validation tests
should execute in less than 10 minutes.

Requirements
------------

The following Python packages are required to execute tests. Some tests will be skipped when
optional requirements are missing.

- gsd (optional)
- mpi4py (optional)
- pytest
- rowan (optional)
- CuPy (optional)

Running tests
-------------

Change to the build directory and execute the following commands to run the tests:

* ``ctest`` - Executes C++ tests
* ``python3 -m pytest hoomd``

pytest_ may be run outside the build directory by:

* Passing a full path to the build: ``python3 -m pytest <build-directory>/hoomd``
* After installing to an environment: ``python3 -m pytest --pyargs hoomd``

.. note::

    ``python3 -m pytest --pyargs hoomd`` tests the hoomd installation it finds by ``import hoomd``,
    which may not be the one you just built. You must also change to a directory outside the
    source, otherwise ``import hoomd`` attempts to import the uncompiled source.

.. seealso::

    See the pytest_ documentation for information on how to control output, select specific tests,
    and more.

.. _CTest: https://cmake.org/cmake/help/latest/manual/ctest.1.html
.. _pytest: https://docs.pytest.org/

Running tests with MPI
----------------------

When ``ENABLE_MPI=ON``, CTest_ will execute some tests with ``mpirun -n 1``, some with ``-n 2``
and some with ``-n 8``. Make sure your test environment (e.g. interactive cluster job) is correctly
configured before running ``ctest``.

pytest_ tests may also be executed with MPI with 2 ranks. pytest_ does not natively support
MPI. Execute it with the provided wrapper script in the build directory::

    mpirun -n 2 build/hoomd/hoomd/pytest/pytest-openmpi.sh -v -x build/hoomd

The wrapper script displays the outout of rank 0 and redirects rank 1's output to a file. Inspect
this file when a test fails on rank 1. This will result in an ``MPI_ABORT`` on rank 0 (assuming the
``-x`` argument is passed)::

    cat pytest.out.1

.. warning::

    Pass the ``-x`` option to prevent deadlocks when tests fail on only 1 rank.

.. note::

    The provided wrapper script supports OpenMPI_.

.. _OpenMPI: https://www.open-mpi.org/

Running validation tests
------------------------

Longer running validation tests do not execute by default. Run these with the ``--validate`` command
line option to pytest::

    $ python3 -m pytest build/hoomd --validate -m validate
    $ mpirun -n 2 hoomd/pytest/pytest-openmpi.sh build/hoomd -v -x -ra --validate -m validate

.. note::

    The ``-m validate`` option selects *only* the validation tests.

.. note::

    To run validation tests on an installed ``hoomd`` package, you need to specify additional
    options::

        python3 -m pytest --pyargs hoomd -p hoomd.pytest_plugin_validate -m validate --validate

Implementing tests
------------------

Most tests should be implemented in pytest_. HOOMD's test rig provides a ``device`` fixture that
most tests should use to cache the execution device across multiple tests and reduce test execution
time.

.. important::

    Add any new ``test_*.py`` files to the list in the corresponding ``CMakeLists.txt`` file.

Only add C++ tests for classes that have no Python interface or otherwise require low level testing.
If you are unsure, please check with the lead developers prior to adding new C++ tests.
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

.. include:: ../CONTRIBUTING.rst
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.error
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.error

.. autosummary::
    :nosignatures:

    DataAccessError
    TypeConversionError
    MutabilityError

.. rubric:: Details

.. automodule:: hoomd.error
    :synopsis: HOOMD errors
    :members: DataAccessError,
              TypeConversionError,
              MutabilityError
    :no-inherited-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.improper
-----------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.improper

.. autosummary::
    :nosignatures:

    Improper
    Harmonic

.. rubric:: Details

.. automodule:: hoomd.md.improper
    :synopsis: Improper potentials.
    :show-inheritance:
    :members: Improper,
              Harmonic,
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

.. include:: ../INSTALLING.rst
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.data
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.data

.. autosummary::
    :nosignatures:

    ForceLocalAccess
    ForceLocalAccessGPU

.. rubric:: Details

.. automodule:: hoomd.md.data
    :synopsis: Provide access in Python to force data buffers on CPU or GPU.
    :members:

    .. autoclass:: ForceLocalAccess
        :inherited-members:

    .. autoclass:: ForceLocalAccessGPU
        :inherited-members:

.. rubric:: Modules
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.long_range.pppm
------------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.long_range.pppm

.. autosummary::
    :nosignatures:

    Coulomb
    make_pppm_coulomb_forces

.. rubric:: Details

.. automodule:: hoomd.md.long_range.pppm
    :synopsis: Long-range potentials evaluated using the PPPM method.
    :members: Coulomb, make_pppm_coulomb_forces
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.angle
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.angle

.. autosummary::
    :nosignatures:

    Angle
    Harmonic
    CosineSquared
    Table

.. rubric:: Details

.. automodule:: hoomd.md.angle
    :synopsis: Angle potentials.
    :members: Angle,
              Harmonic,
              CosineSquared,
              Table
    :no-inherited-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Credits
=======

The following people have contributed to the to HOOMD-blue:

* Aaron Keys, University of Michigan
* Alain Kadar, University of Michigan
* Alex Travesset, Iowa State University and Ames Laboratory
* Alex Yang, Vanderbilt University
* Alexander Hudson
* Allen LaCour, University of Michigan
* Alyssa Travitz, University of Michigan
* Andrew Schultz, University at Buffalo
* Andrew Mark, Max Planck Institute
* Andrey Kazennov, Joint Institute for High Temperatures of RAS
* Antonio Osorio, University of Michigan
* Avisek Das, University of Michigan
* Axel Kohlmeyer, Temple University
* Ben Levine, Temple University
* Ben Swerdlow, University of Michgan
* Benjamin Schultz, University of Michgan
* Bjørnar Jensen, University of Bergen
* Bradley Dice, University of Michigan
* Brandon Butler, University of Michigan
* Brandon Denis Smith, University of Michigan
* Bryan VanSaders, University of Michgan
* Carl Simon Adorf, University of Michigan
* Carolyn Phillips, University of Michigan
* Charlie Slominski, Caltech
* Chengyu Dai, University of Michigan
* Chris Jones, Boise State University
* Christoph Junghans
* Christoph Klein, Vanderbilt University
* Chrisy Du, University of Michigan
* Cong Qiao, Brandeis University
* Corwin Kerr, University of Michigan
* Charlotte Zhao, University of Michigan
* Dan Evans, University of Michigan
* David LeBard, Temple University
* Elizabeth R Chen, University of Michigan
* Eric Harper, University of Michigan
* Eric Irrgang, University of Michigan
* Eric Jankowski, Boise State University
* Erin Teich, University of Michigan
* Fengyi Gao, University of Michigan
* Gabrielle Jones, University of Michigan
* Geert Kapteijns, University of Amsterdam
* Greg van Anders, University of Michigan
* Grey Garrett, University of Michigan
* Igor Morozov, Joint Institute for High Temperatures of RAS
* Isaac Bruss, University of Michigan
* Jakin B. Delony, University of South Florida
* James Antonaglia, University of Michigan
* James Proctor, University of Michigan
* James W. Swan, Massachusetts Institute of Technology
* Jenny Fothergill, Boise State University
* Jens Glaser, Oak Ridge National Laboratory
* Joseph Berleant, University of Michigan
* Joshua A. Anderson, University of Michigan
* Kelly Wang, University of Michigan
* Kevin Daly, Princeton University
* Kevin Kohlstedt, University of Michigan
* Kevin Silmore, Princeton University
* Khalid Ahmed, University of Michigan
* Kristi Pepa, University of Michigan
* Kwanghwi Je, University of Michigan
* Lin Yang, Iowa State University
* Ludwig Schneider, Georg-August Univeristy Goettingen
* Luis Y. Rivera-Rivera, University of Michigan
* Malcolm Ramsay
* Marco Klement, Friedrich-Alexander-Universität Erlangen-Nürnberg (FAU)
* Matthew Spellings, University of Michgan
* Michael Howard, University of Texas
* Mike Henry, Boise State University
* Nathan Horst
* Nipuli Gunaratne, University of Michigan
* Patrick Lawton, University of Michigan
* Paul Dodd, University of Michigan
* Pavani Medapuram Lakshmi Narasimha, University of Minnesota
* Pengji Zhou, University of Michgan
* Peter Palm, Leipzig University
* Peter Schwendeman, University of Michigan
* Philipp Mertmann, Ruhr University Bochum
* Philipp Schönhöfer, University of Michigan
* Rastko Sknepnek, Northwestern
* Raymond Asare, University of Michigan
* Richmond Newman, University of Michigan
* Roman Bystryi, Joint Institute for High Temperatures of RAS
* Ross Smith, University of Michigan
* Ryan Marson, University of Michigan
* Sam Nola, University of Michigan
* Simone Ciarella, Eindhoven University of Technology
* Shannon Moran, University of Michigan
* Sophie YouJung Lee, University of Michigan
* Stephen Thomas, Boise State University
* Steve Barr, Princeton University
* Sumedh R. Risbud, Massachusetts Institute of Technology
* Thi Vo, University of Michigan
* Tim Moore, University of Michigan
* Tobias Dwyer, University of Michigan
* Tommy Waltmann, University of Michigan
* Trung Dac Nguyen, University of Michigan
* Vyas Ramasubramani, University of Michigan
* Wenbo Shen, University of Michigan
* William Zygmunt, University of Michigan
* Wouter Ellenbroek, Eindhoven University of Technology
* Yuan Zhou, University of Michigan
* Åsmund Ervik, SINTEF
* Nathan Barrett, Pritzker School of Molecular Engineering
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.minimize
-----------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.minimize

.. autosummary::
    :nosignatures:

    FIRE


.. rubric:: Details

.. automodule:: hoomd.md.minimize
    :synopsis: FIRE energy minimizer.
    :members: FIRE
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.methods
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.methods

.. autosummary::
    :nosignatures:

    Method
    Brownian
    Langevin
    NPH
    NPT
    NVE
    NVT
    OverdampedViscous


.. rubric:: Details

.. automodule:: hoomd.md.methods
    :synopsis: Integration methods.
    :members: Method,
              Brownian,
              Langevin,
              NPH,
              NPT,
              NVE,
              NVT,
              OverdampedViscous
    :show-inheritance:

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-md-methods-rattle
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.hpmc.nec
---------------

.. automodule:: hoomd.hpmc.nec
    :synopsis: Newtownian event chain Monte Carlo.

.. rubric:: Modules

.. toctree::
   :maxdepth: 3

   module-hpmc-nec-integrate
   module-hpmc-nec-tune
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.force
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.force

.. autosummary::
    :nosignatures:

    Force
    Custom
    Active
    ActiveOnManifold

.. rubric:: Details

.. automodule:: hoomd.md.force
    :synopsis: Apply forces to particles.

    .. autoclass:: Force
        :members:

    .. autoclass:: Custom
        :members:

    .. autoclass:: Active
        :show-inheritance:
        :no-inherited-members:
        :members: create_diffusion_updater

    .. autoclass:: ActiveOnManifold
        :show-inheritance:
        :members: create_diffusion_updater
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.wall
----------

.. py:currentmodule:: hoomd.wall

.. rubric:: Overview

.. autosummary::
    :nosignatures:

    Cylinder
    Plane
    Sphere
    WallGeometry

.. rubric:: Details

.. automodule:: hoomd.wall
    :synopsis: Wall data classes that specify wall geometries for use in other hoomd submodules.
    :members: Cylinder, Plane, Sphere, WallGeometry
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

Features
========

Hard particle Monte Carlo
-------------------------

HOOMD-blue can perform simulations of hard particles using the Monte Carlo method (`hpmc`). Hard
particles are defined by their shape, and the `HPMC integrator <hpmc.integrate>` supports
polygons, spheropolygons, polyhedra, spheropolyhedra, ellipsoids, faceted ellipsoids, spheres,
indented spheres, and unions of shapes. HPMC can make both constant volume and constant pressure
box moves (`hpmc.update.BoxMC`), perform cluster moves (`hpmc.update.Clusters`)
and can compute the pressure during constant volume simulations (`hpmc.compute.SDF`).

HPMC can also apply external and pair potentials to the particles. Use
`hpmc.external.field.Harmonic` to restrain particles to a lattice (e.g. for Frenkel-Ladd
calculations) or `hpmc.external.user.CPPExternalPotential` to implement arbitrary external fields
(e.g. gravity). Use `hpmc.pair.user` to define arbitrary pairwise interactions between particles.
At runtime, `hoomd.version.hpmc_built` indicates whether the build supports HPMC simulations.

.. seealso::

    Tutorial: :doc:`tutorial/00-Introducing-HOOMD-blue/00-index`

Molecular dynamics
------------------

HOOMD-blue can perform molecular dynamics simulations (`md`) with NVE, NVT, NPT, NPH, Langevin,
Brownian, overdamped viscous integration methods (`md.methods`), and energy minimization
(`md.minimize`). Unless otherwise stated in the documentation, all integration methods integrate
both translational and rotational degrees of freedom. Some integration methods support manifold
constraints (`md.methods.rattle`). HOOMD-blue provides a number of pair potentials (`md.pair`)
including pair potentials that depend on particle orientation (`md.pair.aniso`) and many body
potentials (`md.many_body`). HOOMD-blue also provides bond potentials and distance constraints
commonly used in atomistic and coarse-grained force fields (`md.angle`, `md.bond`,
`md.constrain.Distance`, `md.dihedral`, `md.improper`, `md.special_pair`) and can model rigid bodies
(`md.constrain.Rigid`). External fields `md.external.field` apply potentials based only on the
particle's position and orientation, including walls (`md.external.wall`) to confine particles in a
specific region of space. `md.long_range` provides long ranged interactions, including the PPPM
method for electrostatics. HOOMD-blue enables active matter simulations with `md.force.Active` and
`md.update.ActiveRotationalDiffusion`. At runtime, `hoomd.version.md_built` indicates whether the
build supports MD simulations.

.. seealso::

    Tutorial: :doc:`tutorial/01-Introducing-Molecular-Dynamics/00-index`

Python package
--------------

HOOMD-blue is a Python package and is designed to interoperate with other packages in the scientific
Python ecosystem and to be extendable in user scripts. To enable interoperability, all operations
provide access to useful computed quantities as properties in native Python types or numpy arrays
where appropriate. Additionally, `State <hoomd.State>` and `md.force.Force` provide direct access to
particle properties and forces using Python array protocols. Users can customize their simulation or
extend HOOMD-blue with functionality implemented in Python code by subclassing `trigger.Trigger`,
`variant.Variant`, `hoomd.update.CustomUpdater`, `hoomd.write.CustomWriter`,
`hoomd.tune.CustomTuner`, or by using the HOOMD-blue API in combination with other Python packages
to implement methods that couple simulation, analysis, and multiple simulations (such as umbrella
sampling).

.. seealso::

    Tutorial: :doc:`tutorial/04-Custom-Actions-In-Python/00-index`

CPU and GPU devices
-------------------

HOOMD-blue can execute simulations on CPUs or GPUs. Typical simulations run more efficiently on
GPUs for system sizes larger than a few thousand particles, although this strongly depends on the
details of the simulation. The provided binaries support NVIDIA GPUs. Build from source to enable
preliminary support for AMD GPUs. CPU support is always enabled. GPU support must be enabled at
compile time with the ``ENABLE_GPU`` CMake option (see :doc:`building`). Select the device to use at
run time with the `device <hoomd.device>` module. Unless otherwise stated in the documentation,
**all** operations and methods support GPU execution. At runtime, `hoomd.version.gpu_enabled` indicates
whether the build supports GPU devices.

MPI
---

HOOMD-blue can use the message passing interface (MPI) to execute simulations in less time using
more than one CPU core or GPU. Unless otherwise stated in the documentation, **all** operations and
methods support MPI parallel execution. MPI support is optional, requires a compatible MPI library,
and must be enabled at compile time with the ``ENABLE_MPI`` CMake option (see :doc:`building`).
At runtime, `hoomd.version.mpi_enabled` indicates whether the build supports MPI.

.. seealso::

    Tutorial: :doc:`tutorial/03-Parallel-Simulations-With-MPI/00-index`

Threading
---------

Some operations in HOOMD-blue can use multiple CPU threads in a single process. Control this with
the `device.Device.num_cpu_threads` property. In this release, threading support in HOOMD-blue is
very limited and only applies to implicit depletants in `hpmc.integrate.HPMCIntegrator`, and
`hpmc.pair.user.CPPPotentialUnion`. Threading must must be enabled at compile time with the
``ENABLE_TBB`` CMake option (see :doc:`building`). At runtime, `hoomd.version.tbb_enabled` indicates
whether the build supports threaded execution.

Run time compilation
--------------------

Some operations allow the user to provide arbitrary C++ code that HOOMD-blue compiles at run time
and executes during the simulation. `hpmc.pair.user` and `hpmc.external.user` enable users to apply
arbitrary pair and external potentials to particles in HPMC simulations. `hpmc.pair.user`
supports both CPUs and NVIDIA GPUs while `hpmc.external.user` only supports CPUs. Run time
compilation must be enabled at compile time with the ``ENABLE_LLVM`` CMake option (see
:doc:`building`). At runtime, `hoomd.version.llvm_enabled` indicates whether the build supports run
time compilation.

Mixed precision
---------------

HOOMD-blue performs computations with mixed floating point precision. There is a **high precision**
type and a **reduced precision** type. All particle properties are stored in the high precision
type, and most operations also perform all computations with high precision. Operations that do not
mention "Mixed precision" in their documentation perform all calculations in high percision. Some
operations use reduced precision when possible to improve performance, as detailed in the
documentation for each operation. In this release, only `hpmc` implements mixed precision.

The precision is set at compile time with the ``SINGLE_PRECISION`` and
``ENABLE_HPMC_MIXED_PRECISION`` CMake options (see :doc:`building`). By default, the high precision
width is 64 bits and the reduced precision width is 32 bits. At runtime,
`hoomd.version.floating_point_precision` indicates the width of the floating point types.

Plugins
-------

Plugin code that provides additional functionality to HOOMD-blue may be implemented in pure Python
or as a package with C++ compiled libraries.

.. seealso::

    :doc:`components`
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.write
------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.write

.. autosummary::
    :nosignatures:

    DCD
    CustomWriter
    GSD
    Table

.. rubric:: Details

.. automodule:: hoomd.write
    :synopsis: Write data out.
    :members: DCD, CustomWriter, GSD
    :show-inheritance:

    .. autoclass:: Table(trigger, logger, output=stdout, header_sep='.', delimiter=' ', pretty=True, max_precision=10, max_header_len=None)
        :members:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

md.special_pair
---------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.md.special_pair

.. autosummary::
    :nosignatures:

    SpecialPair
    LJ
    Coulomb

.. rubric:: Details

.. automodule:: hoomd.md.special_pair
    :synopsis: Pair potentials between special pairs of particles
    :members: SpecialPair,
              LJ,
              Coulomb
    :no-inherited-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.logging
--------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.logging

.. autosummary::
    :nosignatures:

    log
    Logger
    LoggerCategories

.. rubric:: Details

.. automodule:: hoomd.logging
    :synopsis: Classes for logging data.
    :members: log, Logger
    :undoc-members:

    .. autoclass:: LoggerCategories
        :members: any
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.update
------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.update

.. autosummary::
    :nosignatures:

    BoxResize
    CustomUpdater
    FilterUpdater
    RemoveDrift

.. rubric:: Details

.. automodule:: hoomd.update
    :synopsis: Modify the system state periodically.
    :members: BoxResize, CustomUpdater, FilterUpdater, RemoveDrift
    :imported-members:
    :show-inheritance:
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

hoomd.variant
-------------

.. rubric:: Overview

.. py:currentmodule:: hoomd.variant

.. autosummary::
    :nosignatures:

    Constant
    Cycle
    Power
    Ramp
    Variant

.. rubric:: Details

.. automodule:: hoomd.variant
    :synopsis: Values that vary as a function of time step.
    :no-members:

    .. autoclass:: Constant(value)
        :members: __eq__
        :show-inheritance:
    .. autoclass:: Cycle(A, B, t_start, t_A, t_AB, t_B, t_BA)
        :members: __eq__
        :show-inheritance:
    .. autoclass:: Power(A, B, power, t_start, t_ramp)
        :members: __eq__
        :show-inheritance:
    .. autoclass:: Ramp(A, B, t_start, t_ramp)
        :members: __eq__
        :show-inheritance:
    .. autoclass:: Variant()
        :members: min, max, __getstate__, __setstate__
.. Copyright (c) 2009-2022 The Regents of the University of Michigan.
.. Part of HOOMD-blue, released under the BSD 3-Clause License.

How to model molecular systems
==============================

To model molecular systems using molecular dynamics in HOOMD-blue:

1. Define the bonds, angles, dihedrals, impropers, and constraints as needed to the
   system `State <hoomd.State>`.
2. Add the needed potentials and/or constraints to the integrator (see `md.bond`, `md.angle`,
   `md.dihedral`, `md.improper`, and `md.constrain`).

This code demonstrates bonds by modelling a Gaussian chain:

.. literalinclude:: molecular.py
    :language: python

Consider using a tool to build systems, such as the `Molecular Simulation Design Framework
(MoSDeF)`_. For example, `mBuild`_ can assemble reusable components into complex molecular systems,
and `foyer`_ can perform atom-typing and define classical molecular modeling force fields. The
`mosdef-workflows`_ `demonstrates how to use these tools with HOOMD-blue
<https://github.com/mosdef-hub/mosdef-workflows/blob/master/hoomd_lj/multistate_hoomd_lj.ipynb>`_

.. _Molecular Simulation Design Framework (MoSDeF): https://mosdef.org/
.. _mbuild: https://mbuild.mosdef.org/
.. _foyer: https://foyer.mosdef.org/
.. _mosdef-workflows: https://github.com/mosdef-hub/mosdef-workflows
