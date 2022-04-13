spglib
======
[![Version Badge](https://anaconda.org/atztogo/spglib/badges/version.svg)](https://anaconda.org/atztogo/spglib)
[![Downloads Badge](https://anaconda.org/atztogo/spglib/badges/downloads.svg)](https://anaconda.org/atztogo/spglib)
[![PyPI](https://img.shields.io/pypi/dm/spglib.svg?maxAge=2592000)](https://pypi.python.org/pypi/spglib)
[![Build Status](https://travis-ci.org/atztogo/spglib.svg?branch=master)](https://travis-ci.org/atztogo/spglib)

C library for finding and handling crystal symmetries

How to compile
---------------

    % aclocal
    % autoheader
    % libtoolize # or glibtoolize with macport etc
    % touch INSTALL NEWS README AUTHORS
    % automake -acf
    % autoconf
    % ./configure
* [Lsize]: Linear size---
Title: "Code interna"
Weight: 2
Description: "How to contribute to the code."
---

In this part of the web, some more details on the implementation are discussed.

{{% children depth=1 %}}---
Title: "Object orientation in Octopus"
Weight: 1
---


### Motivation 

Since version 10, {{< octopus >}} is gradually being converted into a fully object oriented code, using the OOP features of Fortran 2003.

The main features of OOP used in {{< octopus >}} are:

* Encapsulation
* Inheritance

The benefits of OOP design are:

* less code duplication
* better readable code
* code can be closer to the physics
* low level aspects can be hidden


### What is an object?

An object refers to a data-structure, which is bundled with the functions (methods) acting on that data. Usually, it is good practice to declare the data as private, 
and allow access only through so-called access functions (getters and setters). This is called encapsulation, and allows later to change details of the implementation
without affecting the rest of the code.

Objects are instances of classes, which define the data structures and the methods for its objects.


### Inheritance and class hierarchy 

Object oriented programming allows to build more complex classes by inheriting the structure and behaviour of other classes.
This reduces duplicated code, and ensures consistency between different classes, which share some functionality.

{{% expand "Example: linked list" %}}
Many classes require the functionality to iterate over a list of objects. These objects, however depend on the nature of the specific class.
In the diagram below, you see the inheritance tree of the {{<code linked_list_t >}} class, which is the parent for many more specialized classes.
{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}
{{% /expand %}}



### Templating

Languages like C++ allow to use templates to write functions or algorithms independent of the data type, they are to be applied to.
Unfortunately, this is not directly possible in Fortran. {{< octopus >}} uses a poor-mans approach to templates, by using macros and including 'templates' of functions, where the explicit types are replaced by macros. These routines are usually defined in files which have '_inc' in their name.

To create a 'templated' function, one can either provide an explicit interface
```Fortran
interface my_function
  module procedure dmy_function, zmy_function
end interface
```
of the function needs to be called using the `X(...)` macro. 

The body of the function is usually defined in a separate file (here {{< code "my_funtion_inc.F90" >}}) which needs to be included as shown here:

```Fortran

#include "undef.F90"
#include "real.F90"
#include "my_function_inc.F90"

#include "undef.F90"
#include "complex.F90"
#include "my_function_inc.F90"
...
```

In {{< code my_function_inc.F90 >}} we would have something like:
```Fortran
function X(my_function)(arg1, arg2) result(res))

  R_TYPE, intent(in)  :: arg1
  R_TYPE, intent(in)  :: arg2
  R_TYPE, intent(out) :: res

  ...

end function X(my_function)
```

The used macros are defined in these include files:

{{% expand "include/real.F90" %}}
```Fortran
#include_file src/include/real.F90
``` 
{{% /expand %}}

{{% expand "include/complex.F90" %}}
```Fortran
#include_file src/include/complex.F90
``` 
{{% /expand %}}

{{% expand "include/integer.F90" %}}
```Fortran
#include_file src/include/integer.F90
``` 
{{% /expand %}}

{{% expand "include/undef.F90" %}}
```Fortran
#include_file src/include/undef.F90
``` 
{{% /expand %}}
---
Title: "Abstract Hamiltonian class"
Weight: 1
Description: "How to contribute to the code."
---

This abstract class does not contains any real data or computational terms, but rather defines the abstract interface, which is inherited by the specific Hamiltonian classes, and some information, which is common to all systems, such as whether the Hamiltonian is Hermitian, and some variables describing the spectral range.

```Fortran
#include_type_def hamiltonian_abst_t
```

This abstract class is then specialized for the various systems: 
{{% graphviz-file "static/graph_data/hamiltonian_abst_t.viz" %}}
---
Title: "Hamiltonian"
Weight: 20
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The Grid in Octopus
====================

Introduction
------------


The fundamental element of Octopus is the real space grid, on which
all quantities, such as wafe functions, densities and potentials are
defined.

In general, a grid consists of

* the mesh
* the simulation box

The mesh is based on an regular mesh in _D_ dimensions, which can be distorted 
using so-called curvi-linear coordinates.

The underlying internal structure is a mesh of integer coordinates and the important
mapping arrays which 
* map the index to the _D_ integer coordinates
* map the integer coordinates to the index.

The grid is constructed in several steps:

1) Generation of a regular mesh filling a rectangular volume enclosing the simulation box.
2) Selecting inner points, enlargment points and ghost points. 
3) re-ordering the points and regeneration of the mapping arrays.



---
Title: "Electron Hamiltonian class"
Weight: 2
Description: "How to contribute to the code."
---

The electronic Hamiltonian is derived from the abstract Hamiltonian. 

{{% expand "Definition of \"hamiltonian_elec_t\"" %}}
```Fortran
#include_type_def hamiltonian_elec_t
```
{{% /expand %}}

It contains the 'physical' quantities, such as 
- information about the dimensionality,
- the contributions to the potential, 
- a pointer to the ions,
- electronic mass,
- etc.

and an instance of {{< code "hamiltonian_elec_base_t" >}}, which bundles some lower level variables.
The separation of quantities into these two classes is mostly historic, and currently has no deeper systematics.

{{% expand "Definition of \"hamiltonian_elec_base_t\"" %}}
```Fortran
#include_type_def hamiltonian_elec_base_t
```
{{% /expand %}}

The method for application of the Hamiltonian, is implemented as a wrapper routine, which checks the compatibility of the supplied wave functions,
and then calls the batched version of the routine.
{{% expand "Applying the Hamiltonian (wrapper)" %}}
```Fortran
#include_subroutine X(hamiltonian_elec_apply)
```
{{% /expand %}}

This batched routine now takes care of the actual application of the various terms of the Hamiltonian. The different terms to be applied
can be selected using the optional {{< code "terms" >}} argument.
{{% expand "Applying the Hamiltonian (batched)" %}}
```Fortran
#include_subroutine X(hamiltonian_elec_apply_batch)
```
{{% /expand %}}

{{% notice "note" %}}
The Hamiltonian is only allowed to modify the wave functions. The reason why also the initial state {{< code "psib" >}} has {{< code "intent(INOUT)" >}} is that the Hamiltonian can pack the wave functions, if requested. 
{{% /notice %}}---
Title: Startup
Weight: 1
---


General Startup procedure of the Octopus Code
=============================================

The flow chart for most calculations has many steps in common, most of them related to setting up the basic data structures of the code.


The 'main' program only performs tasks which are independent of the
actual calculation and the physical system:

[`main/main.F90:`](https://gitlab.com/octopus-code/octopus/blob/develop/src/main/main.F90)

In a strongly simplified way, we have:
```Fortran
program main

  [...] 

  ! "constructors":           ! start code components

  call global_init()               ! initialize the mpi, clocks, etc.
  call parser_init()               ! initialize the input parser
  call messages_init()             ! initialize the message system
  call walltimer_init()            ! initialize the timer module
  call io_init()                   ! initialize the I/O subsystem
  call calc_mode_par_init()        ! initialize parallelization strategy
  call profiling_init()            ! initialize and start the profiling system

  call run(inp_calc_mode)          ! pass control to the 'actual code' running the calculation

  ! "destructors":            ! stop code components in reverse order

  call profiling_end()
  call calc_mode_par_end()
  call io_end()
  call walltimer_end() 
  call messages_end()
  call parser_end() 
  call global_end()

end programme
```


The actual calculation is started from the routine `run()`, which initialized the data structures representing the actual system and calculation mode, before starting the corresponding calculation:

[`main/run.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/main/run.F90) 
defines the module `run_oct_m`:

This module does not contain own data (apart from constants).

{{% expand "Definition of run()" %}}
```Fortran
#include_subroutine run
```
{{% /expand %}}


`scf/ground_state.F90`:

```Fortran
#include_subroutine ground_state_run
```

```Fortran
#include_subroutine ground_state_run_legacy
```

```Fortran
#include_subroutine electrons_ground_state_run
```


---
Title: Main routines
Weight: 60
------
Title: "States"
Weight: 12
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}



{{% children depth=1 %}}---
Title: Wave functions
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

Wavefunctions in Octopus
========================

The wave functions in Octopus are referred to as the states.


They are handled by the module `states_abst_oct_m`, which is defined in [states/states_abst.F90](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states_abst.F90).

The exact way how states are stored in memory is flexible and depends on optimization  (and accelerator) settings in the input file. In particular, the order of indices depends on the `PACKED` setting.

States are stored in a hierarchy of 'containers'. Concepts of this hierarchy 
include _groups_, _batches_ and _blocks_.

#### The abstract class

The top level data structure, describing states is:

```Fortran
#include_type_def states_abst_t
```

This structure contains mainly metadata about states, describing _how_ states are represented, and defines the interface to the class.
As it is an abstract class, it cannot contain any information about the actual system. 

#### The states class for electrons

The class {{< code "states_elec_t" >}}, specialized the abstract class and contains more data specific to the electron system,
as well as pointers to other quantities, which are common to all states, such as the density, the current, etc.

{{% expand "Definition of \"states_elec_t\"" %}}
```Fortran
#include_type_def states_elec_t
```
{{% /expand %}}



The dimensions object contains a number of variables. The most relevant for this discussion is the `dim` variable, which denotes the dimension of one state, being `1` for spin-less states and `2` for spinors.

```Fortran
#include_type_def states_elec_dim_t
```


The wave functions themselves are stored in
```Fortran
    type(states_elec_group_t)     :: group
```
which, in turn, is defined in the module `states_elec_group_oct_m` in [`src/states/states_elec_group.F90`](https://gitlab.com/octopus-code/octopus/blob/develop/src/states/states_elec_group.F90):


The `group` contains all wave functions, grouped together in blocks or batches.
They are organised in an array of `batch_t` structures.

```Fortran
#include_type_def states_elec_group_t
```

```Fortran
    type(wfs_elec_t), pointer   :: psib(:, :)            !< A set of wave-functions
```
The indexing is as follows: `psib(ib,iqb)` where `ib` is the block index, and `iqn` the **k**-point. See below for the routine `states_init_block(st, mesh, verbose)` which creates the `group` object. On a given node, only wave functions of local blocks are available.
The `group` object does contain all information on how the batches are distributed over nodes.


```Fortran
#include_type_def wfs_elec_t
```




Creating the wave functions
---------------------------
A number of steps in initializing the `states_t` object are called from the `system_init()` routine:

`states_init()`:  
parses states-related input variables, and allocates memory for some book keeping variables. It does not allocate any memory for the states themselves.

`states_distribute_nodes()`:  
...


`states_density_init()`:  
allocates memory for the density (`rho`) and the core density (`rho_core`).

`states_exec_init()`:  
1. Fills in the block size (`st\%d\%block_size`);
2. Finds out whether or not to pack the states (`st\%d\%pack_states`);
3. Finds out the orthogonalization method (`st\%d\%orth_method`).
  


Memory for the actual wave functions is allocated in `states_elec_allocate_wfns()` which is called from the corresponding `*_run()` routines, such as `scf_run()` or `td_run()`, etc.


```Fortran
#include_subroutine states_elec_allocate_wfns
```

The routine `states_init_block` initializes the data components in `st` that describe how the states are distributed in blocks:

`st%nblocks`: this is the number of blocks in which the states are divided. 
Note that this number is the total number of blocks,
regardless of how many are actually stored in each node.  
`block_start`: in each node, the index of the first block.  
`block_end`: in each node, the index of the last block.
If the states are not parallelized, then `block_start` is 1 and `block_end` is `st%nblocks`.  
`st%iblock(1:st%nst, 1:st%d%nik)`: it points, for each state, to the block that contains it.  
`st%block_is_local()`: `st%block_is_local(ib)` is `.true.` if block `ib` is stored in the running node.  
`st%block_range(1:st%nblocks, 1:2)`: Block ib contains states fromn st\%block_range(ib, 1) to st\%block_range(ib, 2)  
`st%block_size(1:st%nblocks)`: Block ib contains a number st\%block_size(ib) of states.  
`st%block_initialized`: it should be .false. on entry, and .true. after exiting this routine.  

The set of batches `st%psib(1:st%nblocks)` contains the `block`s themselves.
  
```Fortran
#include_subroutine states_elec_init_block
```



The allocation of memory for the actual wave functions is performed in `batch_init_empty()` and `X(batch_allocate)()`. This routine, and the related `X(batch_add_state)()` show most clearly how the different memory blocks are related.

[`batch_init_empty()`](https://gitlab.com/octopus-code/octopus/blob/develop/src/grid/batch.F90) allocates the memory for `batch_state_t` `states` and `batch_states_l_t` `states_linear` and _nullifies_ the pointers within this types. Note that no memory for the actual wave functions has been allocated yet.

```Fortran
#include_subroutine batch_init_empty
```




Questions:
----------

How are the different objects pointing to states related?

* The usual storage for states is in `states_linear` which can be shadowed by `pack` in case 
  packed states are used.

What is the difference between `batch_add_state` and `batch_add_state_linear`?
---
Title: Batches
Weight: 2
---


In many situations, we need to perform the same operations over many mesh functions, such as the electronic wave functions. 
It is therefore advantageous to group those functions into one object. This can ensure that different mesh functions are contiguous in memory.

Due to the nature of stencil operations, which constitute a large part of the low level operations on mesh functions, it is often more efficient to perform the same stencil operation over different mesh functions (i.e. using the state index as fast index), than looping first over the mesh index, which would, in general, require a different stencil for each mesh point. This is, in particular, the case for calculations utilizing GPUs.

Therefore, we store mesh functions in linear or in so-called packed form. The former refers to the 'natural' ordering where the mesh index is the fastest moving, while the latter is transposed. 

The abstract class {{< code "batch_t" >}} is the parent class for batches, such as electronic {{< versioned-link "Developers/Code_Documentation/States/States" "wave functions" >}}. 

{{% expand "Definition of \"batch_t\"" %}}
```Fortran
#include_type_def batch_t
```
{{% /expand %}}

This class includes information about the dimensions of the functions (number of states, spatial dimension and number of mesh points), but also internal book-keeping variables, 
keeping track of the status of the batch. Furthermore, the {{< code "batch_t" >}} data type contains pointers to the actual data arrays, and defines the methods for interacting with a batch.

Empty batches can be initialized with:
{{% expand "Initializing empty batches" %}}
```Fortran
#include_subroutine X(batch_init)
```
{{% /expand %}}



{{% expand "Initializing batches with memory" %}}
```Fortran
#include_subroutine X(batch_init_with_memory_1)
```

```Fortran
#include_subroutine X(batch_init_with_memory_2)
```

```Fortran
#include_subroutine X(batch_init_with_memory_3)
```
{{% /expand %}}


{{% graphviz-file "static/graph_data/batch_t.viz" %}}


{{% expand "Definition of \"wfs_elec_t\"" %}}
```Fortran
#include_type_def wfs_elec_t
```
{{% /expand %}}
---
Title: Mesh operations
Weight: 5
---



Accessing functions
===================

Function, e.g. the density are given by an array rho(_i_,_spin_), where _i_ ranges over the mesh points associated with the given domain (see domain decomposition), representing

rho( **r**<sub>_i_</sub>,  _spin_ )


What is `mesh_x_global(mesh, ip)` ?

`mesh_x_global(mesh, ip)` returns a `FLOAT` vector, corresponding to the coordinates of point `ip` of the global mesh.

```Fortran
#include_function mesh_x_global
```


Integration
-----------

Integration is implemented as a straightforward summation:

```Fortran
 do ip = 1, mesh%np
      dd = dd + ff(ip)
 end do 
```

Differentiation
---------------

Taking derivatives is done by finite differences. The points involved
in a derivative are defined by the stencil (see below).

Derivatives are discussed in a separate document [Derivatives.md](../Derivatives).


Note on packed states:
----------------------



Note on curvilinear meshes:
---------------------------

The `mesh::x(:,:)` array always contains a regular mesh, which gets 'distorded' to a curvilinear mesh by additional function calls.


Domain decomposition
--------------------

See [Domain_Decomposition.md](Domain_Decomposition.md).



---
title: Coordinate systems
weight: 40
---

```Fortran
#include_type_def coordinate_system_t
```---
Title: Mesh Setup
Weight: 3
---


Setting up a grid
==================




In particular, we have:


* the "first stage":

  This sets up derivatives, double-grid and calls stages 1 and 2 of mesh_init.
At this stage we are still unaware or parallelism and domain decomposition.

{{% expand %}}
```Fortran
#include_subroutine grid_init_stage_1
```
{{% /expand %}}

* the second stage:

{{% expand %}}
```Fortran
#include_subroutine grid_init_stage_2
```
{{% /expand %}}


Mesh initialization:
--------------------

* First stage: set up the mesh_t::idx index structure, defining the lower and upper bounds.

{{% expand %}}
```Fortran
#include_subroutine mesh_init_stage_1
```
{{% /expand %}}


* Second stage: 

{{% expand %}}
```Fortran
#include_subroutine mesh_init_stage_2
```
{{% /expand %}}



* Third stage:

  We can finally generate the list of mesh points mesh_t::x.
At this stage, domain decomposition will be considered and x(:,:) contains
only the domain local points.

{{% expand %}}
```Fortran
! ---------------------------------------------------------
!> When running parallel in domains, stencil and np_stencil
!! are needed to compute the ghost points.
!! mpi_grp is the communicator group that will be used for
!! this mesh.
! ---------------------------------------------------------

#include_subroutine mesh_init_stage_3
```
{{% /expand %}}

which 'serializes' the mesh by a choice of space-filling curve methods.


---
Title: Mesh implementation
Weight: 2
---



The data structures:
--------------------

{{% expand "Definition of basis_set_abst_t" %}}
```Fortran
#include_type_def basis_set_abst_t
```
{{% /expand %}}

{{% graphviz-file "static/graph_data/basis_set_abst_t.viz" %}}


The mesh descriptor is defined in 
`grid/mesh.F90`:

```Fortran
#include_type_def mesh_t
```
which uses the structure `index_t`, containing the range and the mapping arrays:
```Fortran
#include_type_def index_t
```

About lxyz and lxyz_inv maps:
-----------------------------

The direct map:

> lxyz(_d_, _i_) = R<sup>_(d)_</sup><sub>_i_</sub>

where _d_ denotes the dimension (i.e. x, y or z) and _i_ is the index of the point.

The inverse map:

> lxyz_inv( R<sup>x</sup><sub>_i_</sub> , R<sup>y</sup><sub>_i_</sub> , R<sup>z</sup><sub>_i_</sub> ) 
= _i_

The points defined by lxyz define a rectangular box of _d_ dimensions. 
The real mesh vectors are related to R<sub>_i_</sub> by multiplication with spacing
and possibly distortion for curvilinear coordinates.

Note that this index array is the same in all parallel domains, and relates the integer coordinates to the global index of a given point.




```Fortran
#include_type_def mesh_cube_map_t
```







Note on packed states:
----------------------



Note on curvilinear meshes:
---------------------------

The `mesh::x(:,:)` array always contains a regular mesh, which gets 'distorded' to a curvilinear mesh by additional function calls.


Domain decomposition
--------------------



---
Title: "Grid"
Weight: 2
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The Grid in Octopus
====================

Introduction
------------


The fundamental element of Octopus is the real space grid, on which
all quantities, such as wave functions, densities and potentials are
defined.

In general, a grid consists of

* the mesh
* the simulation box

The mesh is based on an regular mesh in _D_ dimensions, which can be distorted 
using so-called curvilinear coordinates.

The underlying internal structure is a mesh of integer coordinates and the important
mapping arrays which 
* map the index to the _D_ integer coordinates
* map the integer coordinates to the index.

The grid is constructed in several steps:

1) Generation of a regular mesh filling a rectangular volume enclosing the simulation box.
2) Selecting inner points, enlargement points and ghost points. 
3) re-ordering the points and regeneration of the mapping arrays.

The top level data structure, describing a grid is defined in
`grid/grid.F90`:

```Fortran
#include_type_def grid_t
```



{{% children depth=5 %}}---
Title: Simulation Box
Weight: 1
---


The simulation box is defined by the class:
```Fortran
#include_type_def simul_box_t
```
which is derived from {{< developers "Code_documentation/ions/box" "box_t" >}}.

The main function of the simulation box is {{< code "simul_box_contains_points" >}}, which determines for a set of points whether they are inside or outside the simulation box.

{{% expand "Definition of simul_box_contains_points()" %}}
```Fortran
#include_function simul_box_contains_points
```
{{% /expand %}}
---
Title: Derivatives
Weight: 10
---

{{< octopus >}} use finite differences to evaluate derivatives, such as the Laplacian (kinetic energy).
The derivative at a point of the mesh is a weighted sum over neighboring points. For instance, 
the general form for the Laplacian is:
$$
    \nabla^2f(n_xh,\,n_yh) = \sum_i^n\sum_j^n\frac{c_{ij}}h\,f(n_xh + ih,\,n_yh+jh)
$$
The coefficients $\(c_{ij}\)$ depend on the mesh and number of points used and define the {{< emph stencil >}}.

The derivatives object contains operators for the gradient and the Laplacian, as well as information about the order and the stencil type, and information about boundaries:
```Fortran
#include_type_def derivatives_t
```

The derivative operators themselves are represented as non-local operators:
```Fortran
#include_type_def nl_operator_t
```

```Fortran
#include_type_def stencil_t
```

```Fortran
#include_type_def stargeneral_arms_t
```
---
Title: "The Self-consistency cycle"
Weight: 1
---

{{% expand "Definition of scf_t" %}}
```Fortran
#include_type_def scf_t
```
{{% /expand %}}

{{% expand "Definition of scf_run()" %}}
```Fortran
#include_subroutine scf_run
```
{{% /expand %}}---
Title: "Electrons"
Weight: 10
Description: "How to contribute to the code."
---

---
Title: "Systems in Octopus"
Weight: 11
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


The {{<developers "Code_Documentation:Multisystem_framework" "multisystem framework">}} allows to define and couple different physical systems.

The following diagram represents the family tree of the system classes. Rounded boxes denote abstract classes, while rectangular boxes are normal classes, which can be instantiated.

{{% graphviz-file "static/graph_data/interaction_partner_t.viz" %}}

### Abstract classes


#### {{< code interaction_partner_t >}}

An interaction_partner in {{< octopus >}} is anything, which can have an interaction with any other system.
For instance electrons and ions are interaction partners, but also the photonic system described by the Maxwell system, or external potentials.
Therefore, it is the base class of all possible systems and multisystems.

{{% expand "Definition of interaction_partner_t" %}}
```Fortran
#include_type_def interaction_partner_t
```
{{% /expand %}}

Each {{< name interaction_partner >}} is associated with a {{< name namespace >}}, owns a {{< name clock >}}, as well as a list of {{< developers "code_documentation/multisystems/interactions" "interactions" >}} in which it can be a partner, and a list of physical quantities, which can be exposed to other systems (through the interactions). See {{< developers "Code_Documentation:Multisystems:quantity" "here" >}} for the list of possible quantities.

It also provides the basic functions to update exposed quantities, and copy them to the interaction. More details about this mechanism are described in the section on 
{{< developers "code_documentation/multisystems/interactions" "interactions" >}}.

#### {{< code system_t >}}

The {{< code system_t >}} type is the abstract type for all systems. 
As all possible systems are potential partners of some interaction, the {{< code system_t >}} type itself extends the abstract {{< code interaction_partner_t >}} type.

{{% expand "Definition of system_t" %}}
```Fortran
#include_type_def system_t
```
{{% /expand %}}

The {{< code system_t >}} class adds information about the physical space in which the system exists, the propagator, and the list of interactions, which are owned by the system. (Check the section {{< developers "Code_Documentation:Multisystems/Interactions" "interactions" >}} for details on who owns an interaction.)


#### {{< code multisystem_t >}}

The {{< code multisystem_t >}}, finally adds a list of systems to the type definietion.

{{% expand "Definition of multisystem_t" %}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}

{{< notice note >}}
{{< code multisystem_t >}} is an abstract class and cannot be used as such to describe a set of systems in the code.
The {{< code type >}} to be used for combined systems is {{< code multisystem_basic_t >}}, which extends {{< code multisystem_t >}}.
{{< /notice >}}




### Specific classes

{{% children depth=2 %}}


Currently, the following system types are defined:
```Fortran
#include_code_doc system_types
```
{{% notice note %}}
When using these system types, **always** use the parameters, and not their numerical values, as they might change over time.
{{% /notice %}}




---
Title: "Multisystem: Containers"
Weight: 2
---

### Introduction

In {{<octopus>}} there is the possibility to group several systems into containers, which are implemented in the {{<code "multisystem_basic_t">}} class.

Containers can have different purposes and applications. In the simplest case, containers are simply a collection of other systems, and do not have their own interactions with anything else. In this case, containers do not introduce any different physics (or approximations), but simply help in the book-keeping of the problem.

{{% notice note %}}
Containers do not correspond to a given region in space, but only to a selection of systems. In many cases, these systems *might* be confined to a certain region in space, but this is not a property of the container. In many other cases, e.g. combining matter and maxwell fields, both systems occupy the same space, or have a substantial overlap.
{{% /notice %}}

Another use case might be to group systems into a container, and then only interact with the whole container, instead of the individual systems. This, however, is an approximation, and furthermore (at least, at the moment) has some limitations due to the implementation.

{{% notice warning %}}
Note, that containers themselves do not move. This has consequences to the definition of the energy contributions. In particular, a container does not have it's own kinetic energy, and all kinetic energy contributions of the constituents are accounted for in the internal energy. For more information, see 
{{<versioned-link "Developers/Code_Documentation/Multisystem_framework/energies" "Calculating Energies">}}
{{% /notice %}}


### Implementation

Most of the functionality is implemented at the level of the abstract class {{<code "multisystem_t">}}:
{{% expand "Definition of multisystem_t"%}}
```Fortran
#include_type_def multisystem_t
```
{{% /expand %}}

The specific {{<code "multisystem_basic_t">}} class only adds the finalizer:
{{% expand "Definition of multisystem_basic_t"%}}
```Fortran
#include_type_def multisystem_basic_t
```
{{% /expand %}}
---
Title: "System: Electrons"
Weight: 2
---

{{% notice warning %}}
The {{<code "electron_t">}} is currently in the process of being refactored. In the end, this class shall describe the electrons alone, which are interacting
with ions, external fields, etc. through interactions. Also the electron-electron interaction, described according to the various theory levels (see {{<variable "TheoryLevel">}}) are implemented through a so-called intra-interaction.
{{% /notice %}}

{{% expand "Definition of \"electrons_t\"" %}}
```Fortran
#include_type_def electrons_t
```
{{% /expand %}}

{{% expand "Definition of the constructor" %}}
```Fortran
#include_function electrons_constructor
```
{{% /expand %}}
---
Title: Abstract box class
Weight: 1
---

```Fortran
#include_type_def box_t
```

{{% graphviz-file "static/graph_data/box_t.viz" %}}
---
Title: "Ions"
Weight: 15
Description: "How to contribute to the code."
---

{{% children depth=5 %}}---
Title: Propagator class
Weight: 1
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


### The propagator class

{{% expand "Definition of propagator_t" %}}
```Fortran
#include_type_def propagator_t
```
{{% /expand %}}

The type {{< code propagator_t >}} is an extension of {{< code algorithm_t >}}. Therefore, it contains the list of operations, which define the propagator.
The elements of the propagator algorithm are defined as {{< developers "Code_Documentation:Multisystems:Algorithms#algorithmic-operations" "algorithmic operations" >}}.
The complete propagation algorithm is then defined by adding each step of the algorithm to the propagator (see the examples).

{{< notice note >}}
Note, that at this level, progators are independent of the actual implementation of each step. These have to be implemented within the system, for which the propagator will be applied.
{{< /notice >}}

Here, we define the following operations:
```Fortran
#include_code_doc general_propagation_operations
```
These operations are general and not bound to a specific propagator, or a specific system. Therefore, they are implemented in the {{< code "system_t" >}} class. For a discussion, see the section on {{< developers "Code_Documentation:Multisystems:Time_propagation" "time propagation" >}}.


The class procedures of {{< code "propagator_t" >}} are those. handling the internal state of the propagator.

Specific propagators are defined as classes extending {{< code "propagator_t" >}}. The necessary specific algorithmic steps are to be defined in the scope of the module file, containing the extending class.

Examples are:
* {{< developers "Code_Documentation:Propagators:Verlet" "Verlet algorithm" >}}
* {{< developers "Code_Documentation:Propagators:Beeman" "Beeman algorithm" >}}
* {{< developers "Code_Documentation:Propagators:Exp_mid" "exponential midpoint algorithm" >}}

{{% expand "Placement in the class hierarchy" %}}
{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}
{{% /expand %}}---
Title: "Propagators"
Weight: 12
Description: "How to contribute to the code."
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Propagators are implemented as a [state machine](https://en.wikipedia.org/wiki/Finite-state_machine). In every iteration step of the main td-loop ({{< developers "Code_Documentation:Multisystems:Time_propagation" "see Time Propagation" >}}), each the propagator is progressed by one algorithmic step.

* {{< developers "Code_Documentation:Propagators:Propagator_class" "Description of the propagator class " >}} 
* {{< developers "Code_Documentation:Propagators:Algorithms" "Description of the algorithm class " >}} 
* Examples are:
  * {{< developers "Code_Documentation:Propagators:Verlet" "Verlet algorithm" >}} (details on how to implement propagators)
  * {{< developers "Code_Documentation:Propagators:Beeman" "Beeman algorithm" >}} (example of predictor-corrector propagators)
  * {{< developers "Code_Documentation:Propagators:Exp_mid" "exponential midpoint algorithm" >}}---
Title: "Algorithms"
Weight: 3
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In {{< octopus >}}, many operations, such as time propagations, geometry optimization, etc. are implemented in terms of algorithms.
An algorithm, in general, contains a set of instructions, which are performed in a well-defined order. 


### Algorithm container

The {{< code algorithm_t >}} class itself, is only an extention of the {{< developers "Miscellanea:Linked_list" "linked list" >}}.

{{% expand "Definition of algorithm_t" %}}
```Fortran
#include_type_def algorithm_t
```
{{% /expand %}}

[Algorithmic operations](#algorithmic-operations) can be added with the function {{< code "add_operation()" >}}. Examples are discussed in the section 
{{< developers "Code_Documentation:Propagators" "Propagators" >}}.

{{% expand "Placement in the class hierarchy" %}}
{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}
{{% /expand %}}

### Algorithm iterator

{{% expand "Definition of algorithm_iterator_t" %}}
```Fortran
#include_type_def algorithm_iterator_t
```
{{% /expand %}}

### Algorithmic operations 

{{% expand "Definition of algorithmic_operation_t" %}}
```Fortran
#include_type_def algorithmic_operation_t
```
{{% /expand %}}

Some global algorithmic steps are defined in {{< source "multisystem/propagator.F90" >}}:
```Fortran
#include_code_doc general_propagation_operations
```

Derived propagators can then add their own steps, such as e.g. in {{< source "multisystem/propagator_verlet.F90" >}}:
```Fortran
#include_code_doc verlet_propagation_operations
```

{{% notice note %}}
It is important to stress here, that one algorithmic step, in general, does not advance any clock. Clocks are only advanced in steps which update the corresponding entity, which could be a system, an exposed quantity or an interaction. It is therefore quite common, that at a specific state of the algorithm, clocks or different entities have different values.
{{% /expand %}}---
Title: "Exponential Midpoint"
Weight: 12
---

For the exponential midpoint propagator, we need to define the following operations:
```Fortran
#include_code_doc exp_mid_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_exp_mid_constructor
```
### The timeline explained



{{< d3-sequence file="graph_data/propagation-exp-mid-3body-equal.json" viewContainers="yes" viewGhosts="yes" >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.
---
Title: "Static Propagator"
Weight: 9
---

The "static propagator" is just a dummy propagator, which does not propagate the system. The only implemented operation is to 
update the interactions, as this is necessary in a multi-system calculation, where other systems are propagated by a "real" propagator.


```Fortran
#include_type_def propagator_static_t

```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_static_constructor
```


---
Title: "Detailled example: Verlet"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### The Verlet algorithm

According to [Wikipedia](https://en.wikipedia.org/wiki/Verlet_integration), the Verlet algorithm is defined as:

- Calculate $\vec{x}(t + \Delta t) = \vec{x}(t) + \vec{v}(t) \Delta t + \tfrac12 \vec{a}(t) \Delta t^2$.
- Derive $\vec{a}(t + \Delta t)$ from the interaction potential using $\vec{x}(t + \Delta t)$.
- Calculate $\vec{v}(t + \Delta t) = \vec{v}(t) + \tfrac12 \big(\vec{a}(t) + \vec{a}(t + \Delta t)\big)\Delta t$.



The Verlet operator is represented by the type {{< code "propagator_verlet_t" >}} which extends {{< code "propagator_t" >}}:

{{% expand "Definition of propagator_verlet_t" %}}
```Fortran
#include_type_def propagator_verlet_t
```
{{% /expand %}}

{{< notice note >}}
Note, that the class definition does not add anything to the {{< code "propagator_t" >}} class. The only differences are the definition of the operations, and the overloaded constructor. 
{{< /notice >}}

For the Verlet propagator, we need to define the following operations:
```Fortran
#include_code_doc verlet_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_verlet_constructor
```
The algorithm also uses steps, which are not specific to the Verlet algorithm, and are defined in the {{< code "propagator_oct_m" >}} module.

{{< notice note >}}
Note, the difference between {{< code "OP_VERLET_FINISH" >}} and {{< code "OP_FINISHED" >}}. The former denotes a specific step, to be taken at the end of one time step, 
while the latter generally denotes the end of a time step.
{{< /notice >}}

As can be seen in this example, the definition of the propagator in terms of the algorithmic operations is a direct translation of the algorithm, shown above.

### Implementation of the steps

So far, the propagator is only defined in an abstract way. It is down to the individual systems (or better, their programmers) to implement the individual steps of the algorithm, in terms of the dynamic variables for that system. An easy examample to demonstrate this is the classical particle, as implemented in {{< code "classical_particles_t" >}}

{{% expand "definition of classical_particles_t" %}}
```Fortran
#include_type_def classical_particles_t
```
{{% /expand %}}
This describes an extremely simple system, consisting of a set of classical, neutral particles. The dynamic variables (i.e. the state) of the system are the positions, velocities and accelerations.

One ''tick'' of the propagator is defined in the function {{< code "classical_particle_do_td" >}}. As can be seen in the code below, this function implements all possible algorithmic steps for all propagators, allowed for that system.


{{% expand "Example implementation for classical particles" %}}
```Fortran
#include_subroutine classical_particles_do_td
```
{{% /expand %}}


### The timeline explained

The Verlet algorithm is a good (because simple) example to illustrate how the systems and the interaction are updated as the code progresses through the main loop.


{{< d3-sequence file="graph_data/propagation-3body-verlet-equal-step.json"  viewContainers="yes" viewGhosts="yes" >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.


</br>

{{% expand "Example with different time steps" %}}
This example shows the propagation of the solar system, where different time steps are used for the three systems.

{{< d3-sequence file="graph_data/propagation-3body-verlet-different-step.json" viewContainers="yes" viewGhosts="yes" >}}

In contrast to the previous example we can see here that the slowest system (i.e. ''sun'') is waiting in ''Updating interaction'' for many computational steps, as it needs to wait for the other systems to complete their time steps. These waiting steps are computationally very cheap, as no grid operations are performed.

{{% /expand %}}
---
Title: "Custom Diagram"
Weight: 20
---

This page allows custom propagation diagrams, for arbitrary time propagation runs.

To generate the required input files, currently a special branch of Octopus needs to be compiled and run. Later, the functionality will be included in main branches.
You need to checkout the ''multisystem_debug'' branch from the {{< octopus-git >}} repository.

In order to enable the output, you need to set
{{< code-block >}}
{{< variable "Debug">}} = propagation_graph
{{< /code-block >}}

{{< notice note >}}
Please, keep in mind that a lot of output is produced, and lower the number of time steps to the minimum to complete a few steps for all subsystems
{{< /notice >}}

Use the button below to navigate to the {{< file "multisystem_propagation.log" >}} of your calculation.

{{< d3-sequence file="load"  viewContainers="yes" viewGhosts="yes" >}}
---
Title: "Beeman"
Weight: 11
---

[Beeman's algorithm](https://en.wikipedia.org/wiki/Beeman%27s_algorithm) is one of the simplest which allows a predictor-corrector implementation, and demonstrates an algorithm, including a self-consistent loop within one propagator time-step. The self-contistent predictor-corrector feature is implemented as optional.


For the Beeman propagator, we need to define the following operations:
```Fortran
#include_code_doc beeman_propagation_operations
```
These are used to define the algorithm, which is done in the constructor of the propagator:
```Fortran
#include_function propagator_beeman_constructor
```

### The timeline explained



{{< d3-sequence file="graph_data/propagation-beeman-3body-equal.json" viewContainers="yes" viewGhosts="yes"  >}}

This graph illustrates how the state machine is stepping through the algorithm. Each system is picking the next algorithmic step from the propagator. For the containers (i.e. ''root'' and ''earth''), the only steps are ''Updating interactions'' and ''Finished''. The {{< emph real >}} systems, on the other hand, are progressing orderly through the operations, defined in the propagator.

{{% expand "Example with different time steps" %}}
This example shows the propagation of the solar system, where different time steps are used for the three systems.

{{< d3-sequence file="graph_data/propagation-beeman-3body-different.json" viewContainers="yes" viewGhosts="yes"  >}}
 {{% /expand %}}
---
Title: "Multisystem framework"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Support for multisystem in {{< octopus >}} is implemented through an object-oriented framework.

Currently, two major modes are implemented: 

1. The legacy mode, which only allows for one ''matter'' system, consisting of electrons, ions and external fields.
2. The new multisystem framework, which allows for several coupled systems, e.g. maxwell, charged particles, etc. 

At the time of writing (Feb. 2021), electrons and ions are not yet available as separate systems.

### Legacy mode

The ''legacy'' mode of the code is used whenever the input file does _not_ have a {{< variable "Systems">}} block.
In this case, the top level system is initialized to be of (the old) {{< code "electrons_t" >}} type. This type describes the combined electron-ion system.
It is planned for te future, that this will be split into the new {{< code electrons_t >}} and {{< code ions_t >}}, which will descibe the electrons and ions as separate systems.
The current {{< code electrons_t >}} is to be replaced by the {{< code matter_t >}}, which is a multisystem, containing electrons and ions.
### Multisystem mode

If the input file contains the {{< variable "Systems" >}} block, the code uses the new multisystems mode.
In this multisystem mode, from the user perspective, the highest level system is a ''multisystem''. Multisystems are containers which can host other system types, including other multisystems. From the code perspective, the {{< code multisystem_t >}} type is a special case of the {{< code system_t >}} type (i.e. it {{< emph extends >}}  {{< code system_t >}}).

The following chapters will discuss in more detail:

* [the multisystem classes](system_classes/)
* [interactions](interactions/)
* [exposed quantities](quantities/)
* [clocks](clocks/)
* [initialization](initialization/)
* [time propagation](time_propagation/)


---
Title: Initialization
Weight: 5
---

### Constructing the systems

In the multisystem framework, it is important to detail how the various components are initialized.
The top level system is created in {{< source "main/run.F90" >}} by calling the {{< code multisystem_basic_t >}} constructor
{{% expand "Definition of multisystem_basic_constructor()" %}}
```Fortran
#include_function multisystem_basic_constructor
```
{{% /expand %}}

which itself calls {{< code "multisystem_basic_init()" >}} and {{< code "multisystem_init()" >}}.

{{% expand "Definition of multisystem_basic_init()" %}}
```Fortran
#include_subroutine multisystem_basic_init
```
```Fortran
#include_subroutine multisystem_init
```
{{% /expand %}}

Finally, {{< code "multisystem_create_system()" >}} loops over the subsystems and calls the system factory to create each of them.
{{% expand "Definition of multisystem_create_system()" %}}
```Fortran
#include_subroutine multisystem_create_system
```
{{% /expand %}}

### Initializing the propagators

The propagators for the system(s) are initialized in {{< code "time_dependent_run_multisystem()" >}} by a call to {{< code "systems%init_propagators()" >}}.
As {{< code "time_dependent_run_multisystem()" >}} is called with {{< code "systems" >}} being of type {{< code "multisystem_basic_t" >}}, this translates to a call to
{{< code "multisystem_init_propagator()" >}}.
{{% expand "Definition of multisystem_init_propagator()" %}}
```Fortran
#include_subroutine multisystem_init_propagator
```
{{% /expand %}}

This routine parses the input for the propagator type, and then creates a propagator of the given type using the respective constructor.

Then, {{< code "multisystem_init_propagator()" >}} loops over the subsystems, and calls {{< code "system%propagator_init()" >}} for each.

In general (i.e. if the function is not overloaded by a derived class), systems use the function {{< code "system_propagator_init()" >}}: 
{{% expand "Definition of system_init_propagator()" %}}
```Fortran
#include_subroutine system_init_propagator
```
{{% /expand %}}

{{< code "system_propagator_init()" >}} parses the input file, looking for the variable {{< variable "TDSystemPropagator" >}}. Here it is important to remember
the order in which the parser treats namespaces.

If the system is specified in the input file as:
```text
%Systems
"System_A" | <type_A>
"System_B" | <type B>
%
```
the top-level multisystem will internally be called {{< code "." >}} and we have in total three namespaces:
```text
. 
./System_A
./System_B
```

The {{< code "system%propagator_init()" >}} routines for {{< code "System_A" >}} and {{< code "System_B" >}} will first look whether {{< variable "TDSystemPropagator" >}} is defined in their respective namespace. If this is _not_ the case, they will look for the variable in the parent namespace, which here is the global namespace {{< code "." >}}.

{{% notice warning %}}
The code does allow to specify different propagators for the multisystem, and the subsystems. While this _might_ work if the subsystems do _not_ interact, it will most likely fail for interacting systems. Therefore, it is highly recommended **not** to specify the propagators for the subsystems separately, unless one knows exactly what one is doing.
{{% /notice %}}


---
Title: "Time propagation"
Weight: 10
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


In the new multisystem framework, the time propagation is handled by the routine {{< code "time_dependent_run_multisystem()" >}}

### The main loop

{{% expand "Implementation of time_dependent_run_multisystem()" %}}
```Fortran
#include_subroutine time_dependent_run_multisystem
```
{{% /expand %}}

The code is quite self-explanatory. 
The central part of the time propagation is the {{< code "do ... while" >}} loop (see {{< code "\"! The full TD loop\"" >}}).
Here the function {{< code "system%dt_operation()" >}} is called successively, triggering one algorithmic step of the multisystem propagator.
{{< notice note >}}
The propagators in general will need several calls of this to finish one time step of the propagator!
{{< /notice >}} 
{{< notice info >}}
At the top level, there is only one system (in general called {{< code "." >}}), which is of type {{< code "multisystem_basic_t" >}}, which then contains other systems as subsystems. These subsystems can be multisystems themselves.
{{< /notice >}}
### Updating the system


In the multisystem framework, the highest level system is of type {{< code "multisystem_basic_t" >}}, which inherits the {{< code "dt_operation()" >}} routine from 
the abstract {{< code "multisystem_t" >}} type.  

{{% expand "Implementation of multisystem_dt_operation" %}}
```Fortran
#include_subroutine multisystem_dt_operation
```
{{% /expand %}}

This routine first calls the general {{< code "system_dt_operation()" >}} subroutine from the parent class {{< code "system_t" >}} and then loops over all
subsystems, which are part of the multisystem, and calls their specific {{< code "system%dt_operations()" >}} routine.

The general {{< code "system_dt_operation()" >}} subroutine is handling all general algorithmic operations, such as starting and ending the SCF loop, or updating interaction. 

{{% expand "Implementation of system_dt_operation()" %}}
```Fortran
#include_subroutine system_dt_operation
```
{{% /expand %}}

This routine fetches the next operation from the propagator and executes it.
The following operations are implemented here:

* {{< code "UPDATE_INTRERACTIONS" >}}
  - increment the **propagator**-clock by one tick (if all interactions can be updated)
  - try to update all interactions (this might not always be possible, if some interaction partners are still lagging behind)
  - if all interactions are updated succesfully, progress the propagator to the next step.
* {{< code "FINISHED" >}}
  - increment the **system**-clock by one tick.
  - rewind the propagator
* {{< code "START_SCF_LOOP" >}} and {{< code "END_SCF_LOOP" >}}
  - handles self-consistency and convergence for predictor-corrector algorithms.
  - no clocks are updated.

All remaining operations are passed down to the system specific routine {{< code "system%do_dt_operation()" >}}, which is the place where specific systems need to implement all algorithmic steps, required by the propagators, allowed for that system.

{{% notice note %}}
In a multisystem, the multisystem itself as well as all its subsystems are propagating according to their _own_ propagators.
In principle, they could be different, which could lead to unexpected behaviour. If propagators of the subsystems are not specified explicitely in the 
input file, they inherit automatically the propagator of the parent multisystem. Therefore, they should *not* be specified individually.
The multisystem container itself always uses an 'empty' propagator, which only implements the operations {{< code OP_UPDATE_INTERACTIONS >}} and {{< code OP_FINISHED >}}.
For {{< code multisystem_t >}} objects, the {{< code "do_td_operations()" >}} routine only implements the {{< code SKIP >}} operation. 
{{% /notice %}}

This routine is triggering the update of the interactions, via the call to {{< code "system%update_interations()" >}},
which loops over the interactions associated with the system, and all required exposed quantities.
### Updating the interactions

Each of the systems (i.e. the multisystem and its subsystems) are attempting to update their interactions with {{< code "system%update_interactions()" >}}.

{{% expand "Implementation of system%update_interations()" %}}
```Fortran
#include_function system_update_interactions
```
{{% /expand %}}

The first part makes sure that the {{< code "update_interactions_start()" >}} routine is only called when no interaction has been updated yet in this time step, iu.e. if their clocks are behind the system clock.

{{< notice note >}}
It is assumed that the interaction can never be ahead of time, compared to the propagator. It is therefore sufficient to ask whether the interaction time equals the propagator time to determine whether an interaction is up-to-date.
{{< /notice >}}

In the second part, we actually attempt to update the interactions, if needed. If an interaction is *not* up-to-date, it first needs to be ensured that the quantities, on which the interaction depends, are up-to-date. Once all quantities are updated, the {{< code "update()" >}} routine of the interaction can be called with the current propagator time as argument. 

Finally, if all interactions are updated successfully, there is a step {{< code "update_interactions_finish()" >}}, which can be necessary for some systems.

{{% expand "Implementation of interaction_with_partner_update()" %}}
```Fortran
#include_function interaction_with_partner_update
```
{{% /expand %}}

This function is to be called with for a specific {{< code requested_time >}}, which is the time at which the interaction is needed.

The function first tries to update the exposed quantities of the partner. If this is successfull, the function calls the {{< code "calculate()" >}} routine of the interaction, sets the interaction clock to the requested time and returns with the result {{< code ".true." >}}, indicating a successful update. In case the function was blocked by one of the exposed quantities, the interaction is not calculated and the function returns {{< code ".false." >}}.


As can be seen, it is, possible that an interaction cannot be updated at a given time. This can be the case when the interaction also depends on quantities of the other interaction partner, and that partner is not yet at the requested time.
 

{{% expand "Definition of system_update_exposed_quantities()" %}}
```Fortran
#include_function system_update_exposed_quantities
```
{{% /expand %}}

This routine demands some more comments:

Firstly, the routine is only allowed if the interaction is of type {{< code "interaction_with_partner_t" >}}, as otherwise there is no partner to update exposed quantities.
It is important to remember who called the function, and which quantities are supposed to be updated:

{{< code "system_update_exposed_quantities()" >}} is called from {{< code "Implementation of interaction_with_partner_update()" >}} for each interaction of a system for the interaction partner of that system with respect to a given interaction.

In the following situations, the routine is **not** allowd to update the exposed quantities:
* we are in an internal SCF loop of the propagation step
* one of the quantities cannot be updated becuase it is either
  - too far behind in time (where two far is more than one clock tick)
  - it is one clock tick behind, but the user requested {{< code "OPTION__INTERACTIONTIMING__TIMING_EXACT" >}}
In these cases, the routine will return a {{< code ".false." >}}. The propagator will continue, but the system in question will be held at this step.
Other systems, however, can continue to progress, and might be able to get to the time, where then this quantity can be updated.

If the timing conditions are fulfilled (and we are not in an internal SCF loop) the system will call the specific routines to update the exposed quantities, and copies their values to the interaction.

### Deadlocks and how to avoid them


---
Title: "Clocks"
Weight: 4
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Clocks are essential to keep things synchronized. This is also the case in {{<octopus>}}, where instances of the {{< code clock_t >}} type are used to ensure that different systems stay in sync during time propagation. Clocks are used to determine whether quantities or interactions need to be updated in a propagation step.

The smallest unit of time in {{< octopus >}} is one {{< code tick >}}.


{{% expand "Definition of clock_t" %}}
```Fortran
#include_type_def clock_t
```
{{% /expand %}}
---
Title: Exposed quantities
Weight: 3
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


Systems can expose quantities that can be used to calculate interactions
with other systems.

Some quantities are dynamical variables of the system. Such quantities are
usually updated by the propagation algorithm and cannot be calculated
on-demand. Such quantities must be marked as "protected".

The module {{< source "multisystem/quantity.F90" "quantity.F90" >}} defines the parameters, which act as index to an exposed quantity within a system.
```Fortran
#include_code_doc quantity
```

Any system, through its base class {{< code "interaction_partner_t" >}} owns an array of type {{< code "quantity_t" >}}

```Fortran
#include_type_def quantity_t
```

This determines whether a quanity is required for a given system, and also associates a specific clock with each quantity.

{{< notice info >}}
''Protected'' quantities deserve some extra mention: Protected quantities are updated in the propagation step of the system, and do not need to be (cannot be) updated explicitely by {{< code "update_exposed_quantity()" >}}. 
{{< /notice >}}---
Title: Interactions
section: Developers
Weight: 2
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}

### Introduction

In physics, interactions are in general symmetric. However, the magnitude with which it affects the interaction partners can be very different.
Also, the way the interactions have to be implemented can differ depending on the ''direction'' of the interaction.

Therefore, in {{< octopus >}} interactions are broken down into their forward and their backward parts, which are treated seperately.


{{% graphviz %}}
digraph {
  "partner A" -> "partner B" [label=" forward   "];
  "partner B" -> "partner A" [label=" backward   "];
 }
{{% /graphviz %}}

As we will see, in some approximations, only one of them (usually the forward interaction) are used.

The ''partners'' are instances of {{< code interaction_partner_t >}} derived types. The interactions (from now on uni-directional) are implemented as derived types of the abstract {{< code interaction_t >}} class.

{{% graphviz %}}
digraph {
  "partner A" -> "partner B" [label=" interaction "];
 }
{{% /graphviz %}}

An interaction ''belongs'' to partner A and acts on partner B.

### Abstract classes

#### {{< code interaction_t >}}

{{% expand "Definition of interaction_t" %}}
```Fortran
#include_type_def interaction_t
```
{{% /expand %}}

This is the minimal data structure for interactions, which only contains information relating to partner A, who owns the interacion.
In particular, the interaction has a list of quantities, which partner A needs to expose, so that the interaction can be calculated.
The IDs for the exposed quantities are defined in the section {{< developers "Code_Documentation:Multisystems:quantity" "Exposed Quantities" >}}.

Furthermore, this abstract type already contains the clock for the interaction, and defines the interfaces for the deferred {{< code "update()" >}} 
and {{< code "calculate()" >}} routines.

{{% graphviz-file "static/graph_data/interaction_t.viz" %}}


#### {{< code interaction_with_partner_t >}}

Curretly, all interactions are derived from {{< code interaction_with_partner_t >}}, which extends {{< code interaction_t >}} and adds the information about partner B, from here on only called ''partner''.

{{% expand "Definition of interaction_with_partner_t" %}}
```Fortran
#include_type_def interaction_with_partner_t
```
{{% /expand %}}

This type adds a pointer to the interaction partner and a list of IDs of the quantities, the partner exposes to the interaction.
Furthermore, at this level, we can implement the {{< code "update()" >}} procedure.


### Abstract children

#### {{< code force_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a force on the partner. Here we only add the actual force vector, acting on the partner system.
{{% expand "Definition of force_interaction_t" %}}
```Fortran
#include_type_def force_interaction_t
```
{{% /expand %}}

#### {{< code potential_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a potential for the partner. 

{{% expand "Definition of potential_interaction_t" %}}
```Fortran
#include_type_def potential_interaction_t
```
{{% /expand %}}


#### {{< code density_interaction_t >}}

The next level of specialization of interaction, are all interactions which create a potential for the partner. 

{{% expand "Definition of density_interaction_t" %}}
```Fortran
#include_type_def density_interaction_t
```
{{% /expand %}}

</br>

### Specific classes:

Specific interaction classes extend the abstract ones. The most important element they add to the abstract classes is the information about the quantities, required to calculate the interaction. In case of the system. owning the interaction (system A), it is sufficient to keep pointers to the data, stored in the system itself. Thr reason is that the interaction is always updated by the propagator of the system A. For the partner system (system B), however, the interaction keeps a copy of the exposed quantities. This allows the partner system to continue the propagation beyond the time for which the quantities are requested, which might happen if the two systems are using different time steps.




#### Ghost interaction

{{% expand "Definition of ghost_interaction_t" %}}
```Fortran
#include_type_def ghost_interaction_t
```
{{% /expand %}}

#### Gravity

{{% expand "Definition of gravity_t" %}}
```Fortran
#include_type_def gravity_t
```
{{% /expand %}}

#### Coulomb force

{{% expand "Definition of coulomb_force_t" %}}
```Fortran
#include_type_def coulomb_force_t
```
{{% /expand %}}

#### Lorentz force

{{% expand "Definition of lorentz_force_t" %}}
```Fortran
#include_type_def lorentz_force_t
```
{{% /expand %}}


### Interaction factory

Instances of {{< code "interaction_t" >}} or derived types are, like systems, generated using a {{< versioned-link "developers/code_documentation/multisystems/factories/#interaction-factories" "factory" >}}.


Currently, the following interaction types are defined:
```Fortran
#include_code_doc interaction_types
```
{{% notice note %}}
When using these system types, **always** use the parameters, and not their numerical values, as they might change over time.
{{% /notice %}}


{{% graphviz-file "/static/graph_data/interaction_graph.dot" %}}---
Title: "Calculating energies"
Weight: 15
---

In the multisysten framework, the total energy ({{<code "total_energy">}}) of a calculation consists of various contributions, where we follow the standard definitions of
thermodynamics (see [here](https://en.wikipedia.org/wiki/Internal_energy))

* The kinetic energy ({{<code "kinetic_energy">}})
* The internal energy ({{<code "internal_energy">}})
* The potential energy ({{<code "potential_energy">}})

### Kinetic energy

We use the term *kinetic energy* in the usual sense. The only slight exception is for container systems:
{{% notice note %}}
Ideally, for containers, the kinetic energy should be the kinetic energy with respect to the centre of mass of the container.
The kinetic energies of the constituents, relative to the centre of mass frame, should be accounted for as part of the internal energy of the container.
{{% /notice %}}

{{% expand "multisystem_update_kinetic_energy()" %}}
```Fortran
#include_subroutine multisystem_update_kinetic_energy
```
{{% /expand %}}

Specific systems need their specific routines to calculate the kinetic energy.
{{% notice note %}}
For the Maxwell system, this is the so-called electromagnetic energy.
{{% /notice %}}

### Interaction energies

Everything else is treated as an interaction energy. Here we have to distinguish between the interactions between two different systems (which can be of the same type), and the interaction of a system with itself (intra-interaction). 

#### inter-interactions

The interactions between two different systems do not require any further thought, as by definition no physical self-interaction can occur. 
As the method to calculate the interactions are part of the interaction class, it is independent of the actual system and can be implemented 
at the level of the class{{<code "system_t">}}.

{{% expand "system_update_potential_energy()" %}}
```Fortran
#include_subroutine system_update_potential_energy
```
{{% /expand %}}

The exception are containers. Here we need to loop over the constituents. In order to distinguish inter- from intra-interactions, we need to
query each interaction for its interaction partner, and skip the interaction, if the partner is part of the container.

{{% expand "multisystem_update_potential_energy()" %}}
```Fortran
#include_subroutine multisystem_update_potential_energy
```
{{% /expand %}}

#### intra-interactions

Systems may contain more than one physical particle (e.g. the electrons, a set of ions or container systems). In order to account for the interaction of these particles with other particles of the same system, we decided to treat this case as a system interacting with itself, which we call intra-interaction. 

In some cases, such as a single particle, this intra interaction has to be zero, while in other cases with many particles, the interactions have to be calculated, where -- of course -- the interaction of one particle with itself has to be removed (at least approximatively).


Another important aspect of the implementation is that {{<octopus>}} deals with _one-sided_ interactions. This has the implication that there is no double counting when calculating
the interaction energies. Both contributions have to be counted: 
**for each system, we add the the energy of the system in the field of the partner.**

{{% expand "system_update_internal_energy()" %}}
```Fortran
#include_subroutine system_update_internal_energy
```
{{% /expand %}}


{{% notice note %}}
For containers, this means that the interaction energy contains the complete interaction energies of all subsystems in the container, _plus_ the energies of the subsystems in the field of all other systems (not part of that container).
{{% /notice %}}

{{% notice warning %}}
Note that the internal_energy of a container is _not_ the sum of the self_energies of the constituents, and also the external energy of a container is _not_
the sum of the external energies of the constituents!
{{% /notice %}}

{{% expand "multisystem_update_internal_energy()" %}}
```Fortran
#include_subroutine multisystem_update_internal_energy
```
{{% /expand %}}





### Total energy

The total energy of the whole simulation (top level container) is clearly defined. It contains the sum of all internal and interaction energies.

For a specific system (e.g. electrons), the total energy also is the sum of the internal energy and the interaction energies (including the intra interaction energy).

{{% expand "system_update_total_energy()" %}}
```Fortran
#include_subroutine system_update_total_energy
```
{{% /expand %}}

 

{{% notice warning %}}
Total energies and interaction energies depend on the state of other systems, and therefore should not be used for convergence criteria for a given system. Only the internal energy is safe to use here.
{{% /notice %}}

---
Title: Factories
Weight: 6
---

Octopus uses the so-called [factory pattern](https://en.wikipedia.org/wiki/Factory_method_pattern) to create instances of the systems and interaction classes.
The abstract factory classes are introduced to avoid the problem of circular dependencies.

The function of the factories is to create an object of a (dynamically) given type and return a pointer to it. This is done by calling the respective 
constructors of the classes, of which an instance is to be created.

### System factories

```Fortran
#include_type_def system_factory_abst_t
```

{{% graphviz-file "static/graph_data/system_factory_abst_t.viz" %}}

```Fortran
#include_type_def system_factory_t
```

{{% expand "Definition of system_factory_create()" %}}
```Fortran
#include_function system_factory_create
```
{{% /expand %}}





### Interaction factories


```Fortran
#include_type_def interactions_factory_abst_t
```

{{% graphviz-file "static/graph_data/interactions_factory_abst_t.viz" %}}

```Fortran
#include_type_def interactions_factory_t
```

{{% expand "Definition of interactions_factory_abst_create_interactions()" %}}
```Fortran
#include_subroutine interactions_factory_abst_create_interactions
```
{{% /expand %}}


{{% expand "Definition of interactions_factory_create()" %}}
```Fortran
#include_function interactions_factory_create
```
{{% /expand %}}

---
Title: "Miscellanea"
Weight: 5
Description: "How to contribute to the code."
---

{{% children depth=5 %}}---
Title: "Linked list"
Weight: 1
---


The {{< code "linked_list_t" >}} type implements a linked list of unlimited polymorphic
values. This allows the storage of any type of data. 

{{% expand "Definition of linked_list_t" %}}
```Fortran
#include_type_def linked_list_t
```
{{% /expand %}}

where the data is stored in {{< code "list_node_t" >}} type objects.

{{% expand "Definition of list_node_t" %}}
```Fortran
#include_type_def list_node_t
```

{{% notice note %}}
The {{<code "class(*), pointer :: value => null()">}} allows storage of a pointer to {{<emph any>}} data type. 
{{% /notice %}}
{{% /expand %}}


Iterating over the list is done using the associated iterator. 
{{% expand "Definition of linked_list_iterator_t" %}}
```Fortran
#include_type_def linked_list_iterator_t
```
{{% /expand %}}


These classes are not meant to used as is, but rather to be extended and by providing an add method to the
list and a get_next method to the iterator.

{{% graphviz-file "static/graph_data/linked_list_t.viz" %}}

Also the {{< code "linked_list_iterator_t" >}} is the parent of a number of specialized iterators for the various systems:

{{% graphviz-file "static/graph_data/linked_list_iterator_t.viz" %}}
---
Title: "Eigensolver"
Weight: 30
---

{{< notice warning >}}
Work in progress!
{{< /notice >}}


---
Title: "Eigensolver"
Weight: 30
---

```Fortran
#include_type_def eigensolver_t
```


```Fortran
#include_subroutine eigensolver_run
```

---
title: "Graphical User Interface"
#series: "Manual"
Weight: 101
hidden: true
draft: true
---



---------------------------------------------
---
title: "About Octopus"
#series: "Manual"
weight: 2
---


{{< octopus >}} is a software package for {{< name "density-functional theory" >}} (DFT), and {{< name "time-dependent density functional theory" >}} (TDDFT).

## Developers

The main development team of this program is composed of:

* Joseba Alberdi
* Xavien Andrade
* Heiko Appel
* Florian Buchholz
* Alberto Castro
* Tilman Dannert
* Umberto De Giovannini
* Alain Delgado
* Nicole Helbig
* Hannes Huebener
* Joaquim Jornet-Somoza
* Ask Larsen
* Irina Lebedeva
* Martin L&uuml;ders
* Miguel A. L. Marques
* Fernando Nogueira
* Micael Oliveira
* Carlo Andrea Rozzi
* Angel Rubio
* Ravindra Shinde
* Jose R. F. Sousa
* David Strubbe
* Iris Theophilou
* Alejandro Varas
* Matthieu Verstraete
* Philipp Wopperer

Former developers:

* Fulvio Berardi
* Johanna Fuks
* David Kammerlander
* Kevin Krieger
* Florian Lorenzen
* Danilo Nitsche
* Roberto Olivares-Amaya
* Arto Sakko
* Axel Thimm
* Jessica Walkenhorst
* Jan Werschnik

Other contributors are:

* Sebastien Hamel: parallel version of oct-excite
* Eugene S. Kadantsev: linear response code

## Introduction

{{< octopus >}} is a pseudopotential real-space package aimed at the simulation of the electron-ion dynamics of one-, two-, and three-dimensional ﬁnite systems subject to time-dependent
electromagnetic ﬁelds. The program is based on time-dependent density-functional theory (TDDFT) in the Kohn-Sham scheme. All quantities are expanded in a regular mesh
in real space, and the simulations are performed in real time. The program has been
successfully used to calculate linear and non-linear absorption spectra, harmonic spectra,
laser induced fragmentation, etc. of a variety of systems. The fundamentals of DFT and
TDDFT can be found, e.g., in the books 
[^footnote-1]
and[^footnote-2]
. 
All information about the octopus
package can be found in its homepage, {{< octopus-web >}}, and in the articles 
[^footnote-3] and [^footnote-4].

The main advantage of real-space methods is the simplicity and
intuitiveness of the whole procedure. First of all, quantities like the density or
the wave-functions are very simple to visualize in real space. Furthermore,
the method is fairly simple to implement numerically for 1-, 2-, or 3-dimensional
systems, and for a variety of different boundary conditions. For example, one
can study a finite system, a molecule, or a cluster without the need of a super-cell,
simply by imposing that the wave-functions are zero at a surface far enough from the system.
In the same way, an infinite system, a polymer, a surface, or bulk material can be
studied by imposing the appropriate cyclic boundary conditions. Note also that
in the real-space method there is only one convergence parameter, namely the grid-spacing, and that decreasing the grid spacing always improves the result.

Unfortunately, real-space methods suffer from a few drawbacks. For example,
most of the real-space implementations are not variational, i.e., we may
find a total energy lower than the true energy, and if we reduce the grid-spacing
the energy can actually increase. Moreover, the grid breaks translational
symmetry, and can also break other symmetries that the system may possess. This can
lead to the artificial lifting of some degeneracies, to the appearance of
spurious peaks in spectra, etc. Of course all these problems can be minimized
by reducing the grid-spacing.

## History

{{< octopus >}} is based on a fixed-nucleus code written by George F. Bertsch and K. Yabana to perform real-time dynamics in clusters 
[^footnote-5]

and on a condensed matter real-space plane-wave based code written by A. Rubio, X. Blase and S.G. Louie 
[^footnote-6]
. 
The code was afterwards extended to handle periodic systems by G.F. Bertsch, J.I. Iwata, A. Rubio, and K. Yabana 
[^footnote-7]
. 
Contemporaneously there was a major rewrite of the original cluster code to handle a vast majority of finite systems. At this point the cluster code was named {{< name "tddft" >}}.

This version was consequently enhanced and beautified by A. Castro (at the time Ph.D. student of A. Rubio), originating a fairly verbose 15,000 lines of {{< name "Fortran 90/77" >}}. In the year 2000, M. Marques (aka Hyllios, aka António de Faria, corsário português), joined the A. Rubio group in Valladolid as a postdoc. Having to use {{< name "tddft" >}} for his work, and being petulant enough to think he could structure the code better than his predecessors, he started a major rewrite of the code together with A. Castro, finishing version 0.2 of {{< name "tddft" >}}. But things were still not perfect: due to their limited experience in {{< name "Fortran 90" >}}, and due to the inadequacy of this language for anything beyond a {{< name "HELLO WORLD" >}} program, several parts of the code were still clumsy. Also the idea of GPLing the almost 20,000 lines arose during an alcoholic evening. So after several weeks of frantic coding and after getting rid of the {{< name "Numerical Recipes" >}} code that still lingered around, {{< octopus >}} was born.

The present released version has been completely rewritten and keeps very little relation to the old version (even input and output files) and has been enhanced with major new flags to perform various excited-state dynamics in finite and extended systems. The code will be updated frequently and new versions can be found here.

If you find the code useful for you research we would appreciate if you give reference to this work and previous ones.

## Contributing to Octopus

If you have some free time, and if you feel like taking a joy ride with {{< name "Fortran 90" >}}, just drop us an email. You can also send us patches, comments, ideas, wishes, etc. They will be included in new releases of octopus.

If you found a have a bug, please report it to our Bug Tracking System: {{< octopus-bugs >}}

## The {{< octopus >}} Copying Conditions

This program is “free”; this means that everyone is free to use it and free to redistribute it on a free basis. What is not allowed is to try to prevent others from further sharing any version of this program that they might get from you.

Specifically, we want to make sure that you have the right to give away copies of the program, that you receive source code or else can get it if you want it, that you can change this program or use pieces of them in new free programs, and that you know you can do these things.

To make sure that everyone has such rights, we have to forbid you to deprive anyone else of these rights. For example, if you distribute copies of the program, you must give the recipients all the rights that you have. You must make sure that they, too, receive or can get the source code. And you must tell them their rights.

Also, for our own protection, we must make certain that everyone finds out that there is no warranty for this program. If these programs are modified by someone else and passed on, we want their recipients to know that what they have is not what we distributed, so that any problems introduced by others will not reflect on our reputation.

The precise conditions of the license are found in the General Public Licenses that accompany it.

Please note that {{< octopus >}} distribution normally comes with some external libraries that are not covered by the GPL license, please see the 
{{< manual "Appendix:Copying" "Copying" >}} Appendix for the copying conditions or these packages.


[^footnote-1]: {{< book title="A Primer in Density Functional Theory" author="C. Fiolhais, F. Nogueira, and M.A.L. Marques (editors)" publisher="Springer Berlin Heidelberg New York" isbn="3-540-03082-2" series="Lecture Notes in Physics" year="2006" issn="0075-8450" >}}

[^footnote-2]: {{< book title="Time-dependent Density Functional Theory" author="M. A. L. Marques and C. A. Ullrich and F. Nogueira and A. Rubio and K. Burke and E. K. U. Gross (editors)" publisher="Springer Berlin Heidelberg New York" isbn="3-540-35422-0" series="Lecture Notes in Physics" year="2006" issn="0075-8450" >}}

[^footnote-3]: {{< article title="octopus: a first principles tool for excited states electron-ion dynamics" authors="M.A.L. Marques, A. Castro, G. F. Bertsch, and A. Rubio" journal="Comp. Phys. Comm." volume="151" pages="60 " year="2003" >}}

[^footnote-4]: {{< article title="octopus: a tool for the application of time-dependent density functional theory" authors="A. Castro, H. Appel, M. Oliveira, C. A. Rozzi, X. Andrade, F. Lorenzen, M.A.L. Marques, E. K. U. Gross, and A. Rubio" journal="Phys. Stat. Sol. (b)" volume="243" pages="2465" year="2006" >}}

[^footnote-5]: {{< article authors="G.F. Bertsch and K. Yabana" title="Time-dependent local-density approximation in real time " journal="Phys. Rev. B" volume="54" pages="4484" year="1996" >}}

[^footnote-6]: {{< article authors="A. Rubio, X. Blase, and S.G. Louie" title="Ab Initio Photoabsorption Spectra and Structures of Small Semiconductor and Metal Clusters" journal="Phys. Rev. Lett." volume="77" pages="247" year="1996" >}}

[^footnote-7]: {{< article authors="G.F. Bertsch, J.I. Iwata, A. Rubio, and K. Yabana" title="Real-space, real-time method for the dielectric function" journal="Phys. Rev. B" volume="62" pages="7998" year="2000" >}}


{{< manual-foot prev="" next="Installation" >}}
---
title: "Installation"
#series: "Manual"
weight: 4
---


Maybe somebody else installed {{< octopus >}} for you. In that case, the files
should be under some directory that we can call {{< inst_file >}}, the executables
in {{< inst_file "bin/" >}} (e.g. if {{< inst_file >}}={{< file "/usr/local/" >}}, the main {{< octopus >}}
executable is then in {{< file "/usr/local/bin/octopus" >}}); the pseudopotential files that {{< octopus >}} will
need in {{< inst_file "share/octopus/PP/" >}}, etc.

However, you may be unlucky and that is not the case. In the following we will try
to help you with the still rather unfriendly task of compiling and installing the {{< octopus >}}.

### Instructions for Specific Architectures

See {{< manual "installation/Specific architectures" "step-by-step_instructions" >}} for some specific supercomputers and generic configurations (including Ubuntu and Mac OSX). If the system you are using is in the list, this will be the easiest way to install. Add entries for other supercomputers or generic configurations on which you have successfully built the code.

## Downloading

Download the latest {{< octopus >}} version here: {{< octopus-download >}}.

<!-- Currently, we do not have a binary distribution -->
<!--
### Binaries

If you want to install {{< octopus >}} in a {{< name "Debian" >}}-based {{< name "Linux" >}} box ({{< name "Debian" >}} or {{< name "Ubuntu" >}}, you might not need to compile it; we release binary packages for some platforms. Keep in mind that these packages are intended to run on different systems and are therefore only moderately optimized. If this is an issue you have to compile the package yourself with the appropriate compiler flags and libraries (see below).

Download the appropriate {{< file ".deb" >}} file from the {{< octopus-download >}}. Install it (using root access) with the command below, using the appropriate filename you downloaded:

{{< command-line "dpkg -i octopus_package.deb" "<nowiki>-</nowiki>" >}} 
-->

### Source code

If you have a different system from those mentioned above or you want to compile {{< octopus >}} you need to get the source code file ({{< file ".tar.gz" >}} file) and follow the compilation instructions below.

## Building
### Quick instructions

For the impatient, here is the quick-start: 
<!-- 
{{< command-line "tar xzf octopus-{{< octopus_version >" >}}.tar.gz}} 
{{< command-line "cd octopus-{{< octopus_version >" >}}}}
{{< command-line "./configure" >}}
{{< command-line "make" >}}
{{< command-line "make install" >}}
-->

```bash
tar xzf octopus-{{< octopus-version >}}.tar.gz
cd octopus-{{< octopus-version >}}
./configure
make
make install
```

This will probably {{< emph "not" >}} work, so before giving up, just read
the following paragraphs.

### Slow instructions

There is an {{< manual "installation/Building from scratch" "appendix" >}} with detailed instructions on how to compile {{< octopus >}} and the required libraries from scratch -- you only need to do this if you are unable to install from a package manager or use per-built libraries.

### Long instructions

The code is written in standard {{< name "Fortran 2003" >}}, with some routines written in {{< name "C" >}} (and in {{< name "bison" >}}, if we count the input parser). To build it you will need both a {{< name "C" >}} compiler ({{< name "gcc" >}} works just fine and it is available for almost every piece of silicon), and a
{{< name "Fortran 2003" >}} compiler. You can check in the {{< manual "Installation/Porting_Octopus_and_Platform_Specific_Instructions" "Compilers Appendix" >}} which compilers {{< octopus >}} has been tested with. This appendix also contains hints on potential problems with certain platform/compiler combinations and how to fix them.

#### Requirements

Besides the compiler, you will also need:

* {{< name "make" >}}: most computers have it installed, otherwise just grab and install the {{< name "GNU make" >}}.

* {{< name "cpp" >}}: The {{< name "C" >}} preprocessor is heavily used in {{< octopus >}} to preprocess {{< name "Fortran" >}} code. It is used for both C (from the CPP variable) and Fortran (FCCPP). {{< name "GNU cpp" >}} is the most convenient but others may work too. For more info, see {{< versioned-link "Developers/Preprocessors" "Preprocessors" >}}.

* {{< name "Libxc" >}}: The library of exchange and correlation functionals. It used to be a part of {{< octopus >}}, but since version 4.0.0 it is a standalone library and needs to be installed independently. For more information, see the [libxc page](https://www.tddft.org/programs/Libxc). {{< octopus >}} 4.0.0 and 4.0.1 require version 1.1.0 (not 1.2.0 or 1.0.0). {{< octopus >}} 4.1.2 requires version 2.0.x or 2.1.x, and won't compile with 2.2.x. (Due to bugfixes from libxc version 2.0 to 2.1, there will be small discrepancies in the testsuite for {{< code "functionals/03-xc.gga_x_pbea.inp" >}} and {{< code "periodic_systems/07-tb09.test" >}}). {{< octopus >}} 5.0.0 supports libxc versions 2.0.x, 2.1.x and 2.2.x. Please note: The Libxc testsuite prior to 2.1 will report some errors in most cases. This is not something to worry about.

* {{< name "FFTW" >}}: We have relied on this great library to perform Fast Fourier Transforms ({{< name "FFTs" >}}). You may grab it from the [{{< name "FFTW" >}} site](https://www.fftw.org/). You require {{< name "FFTW" >}} version 3.

* {{< name "LAPACK/BLAS" >}}: Our policy is to rely on these two libraries as much as possible on these libraries for linear-algebra operations. If you are running {{< name "Linux" >}}, there is a fair chance they are already installed in your system. The same goes to the more heavy-weight machines ({{< name "alphas" >}}, {{< name "IBMs" >}}, {{< name "SGIs" >}}, etc.). Otherwise, just grab the source from [{{< name "netlib" >}} site](https://www.netlib.org).

* {{< name "GSL" >}}: Finally someone had the nice idea of making a public scientific library! {{< name "GSL" >}} still needs to grow, but it is already quite useful and impressive. {{< octopus >}} uses splines, complex numbers, special functions, etc. from {{< name "GSL" >}}, so it is a must! If you don't have it already installed in your system, you can obtain {{< name "GSL" >}} from the [{{< name "GSL" >}} site](https://www.gnu.org/software/gsl/). You will need version 1.9 or higher. Version 4.0 of {{< octopus >}} (and earlier) can only use GSL 1.14 (and earlier). A few tests will fail if you use GSL 1.15 or later. Version 5.0.0 of {{< octopus >}} (and earlier) can only use GSL 1.16 or earlier, due to a bug in our configure script.

* {{< name "Perl" >}}: During the build process {{< octopus >}} runs several scripts in this language. It's normally available in every modern {{< name "Unix" >}} system.


#### Optional libraries


There are also some optional packages; without them some parts of {{< octopus >}} won't work:

* {{< name "MPI" >}}: If you want to run {{< octopus >}} in multi-tentacle (parallel) mode, you will need an implementation of {{< name "MPI" >}}. [MPICH](https://www-unix.mcs.anl.gov/mpi/mpich/) or [Open MPI](https://www.open-mpi.org/) work just fine in our {{< name "Linux" >}} boxes.

* {{< name "PFFT" >}}: We rely on this great library for highly scalable parallel Poisson solver, based on Fast Fourier Transforms ({{< name "FFTs" >}}). You may grab it from the [{{< name "M. Pippig's" >}} site](https://www-user.tu-chemnitz.de/~potts/workgroup/pippig/software.php.en). You also require {{< name "FFTW" >}} version 3.3 compiled with MPI and with a small patch by M. Pippig (also available there).

* {{< name "NetCDF" >}}: The [Network Common Dataform](https://www.unidata.ucar.edu/software/netcdf/) library is needed for writing the binary files in a machine-independent, well-defined format, which can also be read by {{< manual "Visualization" "visualization programs" >}} such as [OpenDX](https://www.opendx.org/)

* {{< name "GDLib" >}}: A library to read graphic files. See {{< tutorial "Particle in an octopus" "Particle in an octopus" >}}. (The simulation box in 2D can be specified via {{< variable "BoxShapeImage" >}}.) Available from [GDLib](https://www.libgd.org).

* {{< name "SPARSKIT" >}}: [Library for sparse matrix calculations](https://www-users.cs.umn.edu/~saad/software/SPARSKIT/). Used for one propagator technique.

* {{< name "ETSF I/O" >}}: An input/output library implementing the ETSF standardized formats, requiring NetCDF, available at [libraries_and_tools](https://www.etsf.eu/resources/software/). Versions 1.0.2, 1.0.3, and 1.0.4 are compatible with {{< octopus >}} (though 1.0.2 will produce a small discrepancy in a filesize in the testsuite). It must have been compiled with the same compiler you are using with {{< octopus >}}. To use ETSF_IO, include this in the <tt>configure</tt> line, where <tt>$DIR</tt> is the path where the library was installed:

{{< code-line "--with-etsf-io-prefix=\"$DIR\"">}}
 

* {{< name "LibISF" >}}: (version 5.0.0 and later) To perform highly scalable parallel Poisson solver, based on BigDFT 1.7.6, with a cheap memory footprint. You may grab it from the [{{< name "BigDFT" >}} site](https://bigdft.org/). You require {{< name "BigDFT" >}} version 1.7.6 compiled with MPI, following these instructions: [installation instructions](https://bigdft.org/Wiki/index.php?title=Installation-Building_the_Poisson_Solver_library_only). Probably, you have to manually copy the files "libwrappers.a" and "libflib.a" to the installation "/lib" directory. To configure {{< octopus >}}, you have to add this configure line:

{{< code-line "--with-isf-prefix=\"$DIR\"" >}}

{{< min-version "10.0" >}}
From version 10.0 on, the ISF library is replaced by its successor [PSOLVER](https://gitlab.com/l_sim/psolver).
The corresponding configure line is:

{{< code-block >}}
  --with-psolver-prefix=DIR
                          Directory where PSolver was installed.
  --with-psolver-include=DIR
                          Directory where PSolver Fortran headers were
                          installed.
{{< /code-block >}}

{{< /min-version >}}
#### Unpacking the sources

Uncompress and untar it ({{< command "gzip -cd octopus-{{< octopus_version >" >}}.tar.gz | tar -xvf -}}). In the following, {{< file "OCTOPUS-SRC/" >}} denotes the source directory of {{< octopus >}}, created by the {{< command "tar" >}} command.


The {{< file "OCTOPUS-SRC/" >}} contains the following subdirectories of interest to users:

* {{< file "doc/" >}}: The documentation of {{< octopus >}}, mainly in HTML format.

* {{< file "liboct_parser/" >}}: The C library that handles the input parsing.

* {{< file "share/PP/" >}}: Pseudopotentials. In practice now it contains the Troullier-Martins (PSF and UPF formats) and Hartwigsen-Goedecker-Hutter pseudopotential files.

* {{< file "share/util/" >}}: Currently, the {{< emph "utilities" >}} include a couple of IBM OpenDX networks ({{< file "mf.net" >}}), to visualize wavefunctions, densities, etc.

* {{< file "testsuite/" >}}: Used to check your build. You may also use the files in here as samples of how to do various types of calculations.

* {{< file "src/" >}}: Fortran90 and C source files. Note that the Fortran90 files have to be preprocessed before being fed to the Fortran compiler, so do not be scared by all the - directives.

#### Development version

You can get the development version of {{< octopus >}} by downloading it from the {{< octopus >}} project on {{< octopus-git >}}.

You can also get the current version with the following command (you need the {{< name "git" >}} package):

{{< command-line "git clone git@gitlab.com:octopus-code/octopus.git" >}}

Before running the configure script, you will need to run the GNU autotools. This may be done by executing:

{{< command-line "autoreconf -i" >}}

Note that you need to have working recent versions of the {{< name "automake" >}} and {{< name "autoconf" >}}. In particular, the configure script may fail in the part <tt>checking for Fortran libraries of mpif90</tt> for <tt>autoconf</tt> version 2.59 or earlier. The solution is to update <tt>autoconf</tt> to 2.60 or later, or manually set <tt>FCLIBS</tt> in the <tt>configure</tt> command line to remove a spurious apostrophe.

If autoreconf is failing with "aclocal: warning: couldn't open directory 'm4': No such file or directory", create an empty folder named m4 inside external_libs/spglib-1.9.9/.

Please be aware that the development version may contain untested changes that can affect the execution and the results of {{< octopus >}}, especially if you are using new and previously unreleased features. So if you want to use the development version for production runs, you should at least contact {{< octopus >}} developers.

#### Configuring

Before configuring you can (should) set up a couple of options. Although
the {{< name "configure" >}} script tries to guess your system settings for you, we recommend 
that you set explicitly the default Fortran compiler and the compiler options.
Note that {{< name "configure" >}} is a standard tool for Unix-style programs and you can find a lot of generic documentation on how it works elsewhere.

For example, in {{< name "bash" >}} you would typically do:

{{< command-line "export FC=ifort" >}}
{{< command-line "export FCFLAGS=\"-O2 -xHost\"" >}}

if you are using the Intel Fortran compiler on a linux machine.

Also, if you have some of the required libraries in some unusual directories,
these directories may be placed in the variable {{< code "LDFLAGS" >}} (e.g.,
{{< command "export LDFLAGS=$LDFLAGS:/opt/lib/" >}}).

The configuration script will try to find out which compiler you are using.
Unfortunately, and due to the nature of the primitive language that {{< octopus >}}
is programmed in, the automatic test fails very often. Often it is better to set
the variable {{< code "FCFLAGS" >}} by hand, check the {{< manual "Installation/Porting_Octopus_and_Platform_Specific_Instructions" "Compilers Appendix" >}} page for which flags have been reported to work with different Fortran compilers.

You can now run the configure script 
{{< command-line "./configure" >}}

You can use a fair amount of options to spice {{< octopus >}} to your own taste. To obtain a full list just type {{< command "./configure --help" >}}. Some 
commonly used options include:

* {{< flag "--prefix=" >}}{{< inst_file >}}: Change the base installation dir of {{< octopus >}} to {{< inst_file >}}. {{< inst_file >}} defaults to the home directory of the user who runs the {{< name "configure" >}} script.

* {{< flag "--with-fft-lib=" >}}{{< file "<lib>" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "FFTW" >}} library exactly in the way that it is specified in the {{< file "<lib>" >}} argument. You can also use the {{< code "FFT_LIBS" >}} environment variable.

* {{< flag "--with-pfft-prefix=" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "PFFT" >}} library.

* {{< flag "--with-pfft-lib=" >}}{{< file "<lib>" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "PFFT" >}} library exactly in the way that it is specified in the {{< file "<lib>" >}} argument. You can also use the {{< code "PFFT_LIBS" >}} environment variable.

* {{< flag "--with-blas=" >}}{{< file "<lib>" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "BLAS" >}} library in the way that it is specified in the {{< file "<lib>" >}} argument.

* {{< flag "--with-lapack=" >}}{{< file "<lib>" >}}: Instruct the {{< name "configure" >}} script to look for the {{< name "LAPACK" >}} library in the way that it is specified in the {{< file "<lib>" >}} argument.

* {{< flag "--with-gsl-prefix=" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "GSL" >}} library. The libraries are expected to be in {{< file "DIR/lib/" >}} and the include files in {{< file "DIR/include/" >}}. The value of {{< file "DIR/" >}} is usually found by issuing the command {{< command "gsl-config --prefix" >}}.

* {{< flag "--with-libxc-prefix=" >}}{{< file "DIR/" >}}: Installation directory of the {{< name "Libxc" >}} library.

If you have problems when the {{< name "configure" >}} script runs, you can find more details of what happened in the file {{< name "config.log" >}} in the same directory.

#### Compiling and installing

Run {{< command "make" >}} and then {{< command "make install" >}}. The compilation may take some time, so you might want to speed it up by running {{< command "make" >}} in parallel ({{< command "make -j" >}}). If everything went fine, you should now be able to taste {{< octopus >}}. 

Depending on the value given to the {{< flag "--prefix" >}}={{< inst_file >}} given, the executables will reside in {{< inst_file "bin/" >}}, and the auxiliary files will be copied to {{< inst_file "share/octopus" >}}.

#### Testing your build

After you have successfully built {{< octopus >}}, to check that your build works as expected there is a battery of tests that you can run. They will check that {{< octopus >}} executes correctly and gives the expected results (at least for these test cases). If the parallel version was built, the tests will use up to 6 MPI processes, though it should be fine to run on only 4 cores. (MPI implementations generally permit using more tasks than actual cores, and running tests this way makes it likely for developers to find race conditions.)

To run the tests, in the sources directory of {{< octopus >}} use the command

{{< command-line "make check" >}}

or if you are impatient,

{{< command-line "make check-short" >}}

which will start running the tests, informing you whether the tests are passed or not. For examples of job scripts to run on a machine with a scheduler, please see {{< manual "Installation/Specific_architectures" "Specific_architectures" >}}.

If all tests fail, maybe there is a problem with your executable (like a missing shared library). 

If only some of the tests fail, it might be a problem when calling some external libraries (typically blas/lapack). Normally it is necessary to compile all Fortran libraries with the same compiler. If you have trouble, try to look for help in the [{{< octopus >}} mailing list](https://www.tddft.org/mailman/listinfo/octopus-users).

#### Fast recompilation

NOTE: This feature is currently only available in the development version and the plan is to include it in the {{< octopus >}} 9 release.

If you have already compiled the code and if you are changing only one file, you can run

{{< command-line "make NODEP=1" >}}

to ignore the dependencies due to Fortran module files. This will only compile the files that have changed and link the executables; therefore, it is much faster. If you change, e.g., interfaces of modules or functions, you need to to run {{< command "make" >}} without {{< command "NODEP=1" >}} to ensure a correct handling of the dependencies.


{{< manual-foot prev="Manual:About Octopus" next="Manual:Input file" >}}
---------------------------------------------
---
Title: "Manual"
weight: 20
description: "the user manual"
menu: "top"
---

Welcome to the Octopus online manual.





{{% children depth=1 sort="Weight" %}}---
title: "About this manual"
#series: "Manual"
weight: 1
---


## Copying

'''This manual is for {{< octopus >}} {{< octopus-version >}}, a first principles, electronic structure, excited states, time-dependent density functional theory program.'''

''Copyright © 2006 Miguel A. L. Marques, Xavier Andrade, Alberto Castro and Angel Rubio''

Permission is granted to copy, distribute and/or modify this document under the terms of the GNU Free Documentation License, Version 1.1 or any later version published by the Free Software Foundation; with no Invariant Sections, no Front-Cover Texts, and no Back-Cover Texts.

You can find a copy of the Free Documentation License [here](https://www.gnu.org/copyleft/fdl.html)

## Reading this manual

If you are new to {{< octopus >}} it is recommended that you read this manual sequentially. To read this manual you need some basic knowledge of {{< name "Density Functional Theory" >}} ({{< name "DFT" >}}) and {{< name "Time Dependent Density Functional Theory" >}} ({{< name "TDDFT" >}}), also a bit of numerical methods may be useful.

The way to tell what to do to {{< octopus >}} is through input variables, each one of these variables has its own proper documentation describing what it does and what are the possible values that it can take. All the documentation of the variables is contained in the {{< versioned-link "Variables" "Variables Reference" >}}. That information is complementary to this manual, here we will only give a brief description of the functionality of each variable and we will leave the details for the Reference Manual. If you are reading an electronic version of this manual, all references to input variables will be marked as hyperlinks. Opening that link will take you to a page with the reference documentation for that variable. 

There are some styles used in this documents: "{{< name "proper names" >}}" (like {{< name "Fortran" >}}), "{{< file "files and directories" >}}" (like {{< file "inp" >}} or {{< file "/tmp/" >}}), "{{< command "commands" >}}" (like {{< command "ls" >}}) and "{{< code "input file text" >}}" (like {{< code "a=0.5" >}}).

<!--
## Getting help

If you need help using {{< octopus >}}, you have any doubts about the content of this manual or you want to contact us for any reason, please do it so through the [octopus users mailing list](https://tddft.org/mailman/listinfo/octopus-users).

## Contributing

This manual is developed using a Wiki, this means that you can directly cooperate with it. You can always find the original version in https://www.tddft.org/programs/octopus/wiki/index.php/Manual.

-->


{{< manual-foot prev="" next="About Octopus" >}}
---
title: "Pseudopotentials"
#series: "Manual"
weight: 100
description: " "
---


{{< octopus >}} is capable of using different formats for pseudopotentials (see {{< variable Species >}} for more information). We include in the distribution some LDA Troullier-Martins pseudopotentials for some common atoms, but it is likely that you will need other pseudopotentials.

#### Repositories on the web 

These are places on the web that include pseudopotential files in a variety of formats. Note that not all formats are read by {{octopus}} but many are.

* [NNIN Virtual Vault for Pseudopotentials](https://www.nnin.org/nnin_comp_psp_vault.html) Contains an extensive list of pseudo-potential databases and a search engine. It also has links to pseudo-potential related tools.
* [ABINIT](https://www.abinit.org) has a quite complete [set of pseudopotentials](https://www.abinit.org/atomic-data-files). At this moment only <tt>.fhi</tt> files are supported.
* [QuantumEspresso](https://www.quantum-espresso.org/)'s format is also accepted. So you can use the <i>norm-conserving</i> pseudopotentials from their [repository](https://www.quantum-espresso.org/pseudopotentials).

#### Pseudopotential generators 

In order to generate your own pseudopotentials, we advise the use of the following programs:
* [APE](https://www.tddft.org/programs/APE) 
* [opium](https://opium.sourceforge.net/) 
* [FHI98PP](https://www.fhi-berlin.mpg.de/th/fhi98md/fhi98PP/) 

---
title: "Curvilinear Coordinates"
#series: "Manual"
---

---
title: "Deprecated Utilities"
#series: "Manual"
Weight: 100
---


Some utilities present in older version of octopus have superseded by other utilities or integrated into the main code, in this section you can find how to replace them.

### octopus_cmplx

Not exactly a utility, now there is no complex version of the code, the main executable {{< command "octopus" >}} can calculate real or complex wavefunctions depending of the requirements of input file.

### oct-sf

This utility was replaced by [[Manual:External utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]].

### oct-cross-section

This utility was replaced by [[Manual:External utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]].

### oct-hs-mult and oct-hs-acc

These utilities were replaced by [[Manual:External utilities:oct-harmonic-spectrum | {{< command "oct-harmonic-spectrum" >}}]].

### oct-excite

This utility was used to calculate the excitation spectrum within linear response. It no longer exists because the linear-response formalism calculations are now done by the main code, by making use of the "casida" run mode. See Manual:Calculation Modes:Casida

### oct-make-st

This utility was used to replace some of the Kohn-Sham states in the restart directory by Gaussians wave packets. This sort of operation can now be done using the {{< variable "UserDefinedStates" >}} variable.

### oct-rotatory_strength and oct-rsf

This utility has been removed. It is replaced by a new {{< variable "PropagationSpectrumType" >}} in {{< file "oct-propagation_spectrum" >}}.

{{< manual-foot prev="Manual:External utilities:oct-xyz-anim" next="Manual:Examples:Hello world" >}}
---------------------------------------------
---
title: "Advanced ways of running Octopus"
#series: "Manual"
---


## Parallel {{< octopus >}}

{{< octopus >}} can be run in multi-tentacle (parallel) mode if you have configured it that way when you compiled and installed it. There are four possible parallelization strategies: domains, states, spin/''k''-points, and electron-hole pairs (only for Casida). You do not need to specify in the input file what kind of parallelization strategy you want to use as {{< octopus >}} can find reasonable defaults. However, the experienced user may want to tune the corresponding variables.

### Parallel in States 
The first strategy is parallel in states. This means that each processor takes some states and propagates them independently of the others. 

To use this include 

{{< code-block >}}
 {{< variable "ParStates" >}} = auto
{{< /code-block >}}

in your input file.

This is clearly the most efficient (and simplest) way of parallelizing the code if the number of states is much
larger than the number of processors available (which is usually the case for us). You can expect a more or less linear speed-up with number of processors until you reach
around 5 states/processor. Beyond this it starts to saturate.

### Parallel in Domains 

When running parallel in "domains", {{< octopus >}} divides the simulation region (the box) into separate regions (domains) and assigns each of of these to a different processor. This allows it not only to speed up the calculation, but also to divide the memory among the different processors. The first step of this process, the splitting of the box, is in general a very complicated process. Note that we are talking about a simulation box of an almost arbitrary shape, and of an arbitrary number of processors. Furthermore, as the communication between processors grows proportionally to the surface of the domains, one should use an algorithm that divides the box such that each domain has the same number of points, and at the same time minimizes the total area of the domains. {{< octopus >}} uses the [METIS](https://www-users.cs.umn.edu/~karypis/metis/) library to make this domain decomposition of the space, METIS is included in the source code of {{< octopus >}} and will be used by default.

So for example, if you are planning to run the ground state of Na there's only one choice:

{{< code-block >}}
 {{< variable "ParDomains" >}} = auto
{{< /code-block >}}

{{< octopus >}} will then partition the real-space mesh and compute the resulting
pieces on different nodes. Later when you plan to do time propagations
you have more choices at hand. You could run only parallel in states

{{< code-block >}}
 {{< variable "ParStates" >}} = auto
 {{< variable "ParDomains" >}} = no
{{< /code-block >}}

or in domains (as above), or you can decide to employ both parallelization
strategies

{{< code-block >}}
 {{< variable "ParStates" >}} = auto
 {{< variable "ParDomains" >}} = auto
{{< /code-block >}}

### Problems with parallelization 

The parallelization is only effective if a sufficiently large amount of
load is given to each node. Otherwise communication dominates. In that
case one should do the calculation in serial. This happens, for example,
if you attempt to do parallelization in real-space domains, but your grid
does not have too many points. If the code detects that the number of
points per node is not large enough, it gives the warning:

```text
 ** Warning:
 ** From node =    0
 ** I have less elements in a parallel group than recommended.
 ** Maybe you should reduce the number of nodes 
```

### Compiling and running in parallel 
The steps to compile and use {{< octopus >}} in parallel are:

1) use the {{< code "--enable-mpi" >}} in the configuration script 

2) give to {{< octopus >}} the libraries, etc. necessary to compile a parallel
code. When you install mpich in your system, it creates a series of
programs called mpicc, mpif77, mpif90, etc. These are wrappers that call
the specific compiler with the appropriate flags, libraries, etc. This
means that to compile {{< octopus >}} in parallel you'll need mpif90. Please check
that this is present in your system (perhaps in {{< file "/usr/local/bin" >}} if you
compiled mpich from the sources). {{< octopus >}} will try to find this compiler,
but it may not be able to find it for several reasons:

* {{< file mpif90 >}} is _not_ in the path.

* there is another f90 compiler specified by the user/system. This means that if you have the variable FC  defined when you run the configure script, {{< octopus >}} will use the fortran compiler specified by that variable, and will not try to find mpif90. Then, the configuration fails.

```bash 
 export FC=<whatever path>/mpif90
 export FCFLAGS=<whatever flags>
 ./configure --enable-mpi
```

It is wise to also set the {{< code FCFLAGS >}}, as {{< octopus >}} cannot tell which fortran
compiler you are using if you compile with mpif90. 

3) Now you just have to run {{< octopus >}} in parallel (this step depends on your actual system, you may have to use mpirun or mpiexec to accomplish it). 

In some run modes (e.g., {{< code td >}}), you can use the multi-level parallelization, i.e., to run in parallel in more than one way at the same time. In the {{< code td >}} case, you can run parallel in states and in domains at the same time. In order to fine-tune this behavior, please take a look at the variables {{< variable "ParStates" >}}, {{< variable "ParDomains" >}}, {{< variable "ParKPoints" >}}, and {{< variable "ParOther" >}}. In order to check if everything is OK, take a look at the output of {{< octopus >}} in the "Parallelization" section. This is an example:

```text
 ************************** Parallelization ***************************
 Octopus will run in *parallel*
 Info: Number of nodes in par_states  group:     8 (      62)
 Info: Octopus will waste at least  9.68% of computer time
 **********************************************************************
```

In this case, {{< octopus >}} runs in parallel only in states, using 8 processors (for 62 states). Furthermore, some of the processors will be idle for 9.68% of the time (this is not that great, so maybe a different number of processors would be better in this case).

## Passing arguments from environment variables

{{< octopus >}} can also read input variables from the environment variables of the shell. This is especially useful if you want to call {{< octopus >}} from scripts or if you want to pass one-time arguments when running {{< octopus >}}.

To activate this feature you must first set {{< code "OCT_PARSE_ENV" >}} to some value, for example in sh/bash

{{< command-line "export OCT_PARSE_ENV=1" >}}

After that, to pass an input variable to {{< octopus >}}, you must prepend the name of the variable by {{< code "'OCT_'" >}}, for example:

{{< command-line "export OCT_Radius=10" >}}

After that you can run {{< octopus >}} as usual.

You can also pass the variables in the same command line:

{{< command-line "OCT_PARSE_ENV=1 OCT_Radius=10 octopus" >}}

If you set a variable both from the input file and as a environment variable, the environment variable takes precedence. Be careful with this behaviour. Blocks are not supported (suggestions as how to elegantly put a block inside of environment variables are welcome).

There is an additional environment variable {{< code "OCTOPUS_SHARE" >}} that can be set to specify the location of the {{< file variables >}} file used by the parser. If this environment variable is not set, the default is as set by the configure script, namely {{< file "prefix/share/octopus" >}}. You can set this to something else if necessary, for example to handle moving an executable from one machine to another, where the original path does not exist on the new machine.

### Examples

For example this is a simple script that helps to determine the optimal Radius of the simulation box (an {{< file "inp" >}} file with the rest of the parameters has to be provided):

```bash
export OCT_PARSE_ENV=1
for OCT_Radius in `seq 5 16`
do
  export OCT_Radius
  octopus >& out-$OCT_Radius
  energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
  echo $OCT_Radius $energy
done
```


{{< manual-foot prev="Basics:Visualization" next="External utilities:oct-analyze_projections" >}}
---
title: "Visualization"
#series: "Manual"
Weight: 102
---


Every given number of time iterations, or after ground-state calculations, some of the functions that characterise the system may be written to disk so that they may be analized. Files are written within {{< file "static/" >}} output directory after the self-consistent field, or within {{< file "td.x/" >}} directories, during evolution, where “x” stands for the iteration number at which each write is done. 

The function that you want to plot is selected by the {{< variable "Output" >}} variable and the output format is chosen by the 
{{< variable "OutputFormat" >}}.

### dx

This is an OpenDX network, aimed at the visualization of wavefunctions. To be able to use it, you need to have properly installed the [OpenDX](https://opendx.org) program, as well as the Chemistry extensions developed at the Cornell Theory Center. Unfortunately, since this software is old and unmaintained, you are likely to have trouble finding and installing it.

Once these are working, you may follow the {{< versioned-link "tutorial/basics/visualization/#benzene-molecule" "tutorial for the benzene molecule" >}}
### XCrySDen

Atomic coordinates (finite or periodic), forces, and functions on a grid can be plotted with the free program [https://www.xcrysden.org/ XCrySDen]. Its XSF format also can be read by [V_sim](https://inac.cea.fr/sp2m/L_Sim/V_Sim/index.en.html) and [Vesta](https://www.geocities.jp/kmo_mma/crystal/en/vesta.html). Beware, these all probably assume that your output is in Angstrom units (according to the [specification](https://www.xcrysden.org/doc/XSF.html)), so use UnitsOutput = eV_Angstrom, or your data will be misinterpreted by the visualization software.

### CUBE

The Gaussian cube format (see https://gaussian.com/cubegen/, https://www.ks.uiuc.edu/Research/vmd/plugins/molfile/cubeplugin.html and https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html for a more detailed description of the format) can be output, and can be read by VMD, XCrysDen, Avogadro, and other software. Note that CUBE files are always in atomic units, so the UnitsOutput input option will be ignored.

### PDB

Everything is supposed to be in Angstroms: https://deposit.rcsb.org/adit/docs/pdb_atom_format.html

### XYZ

Generally considered to be in Angstroms: https://openbabel.org/wiki/XYZ_%28format%29, https://en.wikipedia.org/wiki/XYZ_file_format, https://www.molpro.net/info/2012.1/doc/manual/node100.html, https://departments.icmab.es/leem/siesta/Documentation/Manuals/siesta-3.1-manual/node32.html

{{< manual-foot prev="Manual:Geometry Optimization" next="Manual:Advanced ways of running Octopus" >}}
---------------------------------------------
---
title: "Running Octopus"
section: "Manual"
weight: 2
description: " "
---


### Input

In order to run, {{< octopus >}} requires an input file that ''must'' be called {{< file "inp" >}}. Depending on your input file there are other files that you might need, like pseudopotentials or coordinate files (we will discuss this later in this manual).

The rest of the files that are required are part of the {{< octopus >}} installation; if {{< octopus >}} is correctly installed they should be available and {{< octopus >}} should know where to find them. With {{< octopus >}} you can't just copy the binary between systems and expect it to work.

### Executing

To run {{< octopus >}} you just have to give the {{< command "octopus" >}} command in the directory where you have your input file. While running, {{< octopus >}} will display information on the screen. If you want to capture this information you can send the output to a file, let's say {{< file "out.log" >}}, by executing it like this:

```bash
octopus > out.log
```

This captures only the normal output. If there is a warning or an error, it will be still printed on the screen. To capture everything to a file, run

```bash
octopus >& out.log
```

If you want to run {{< octopus >}} in the background, append {{< command "&" >}} to the last command.

### Output

While running, {{< octopus >}} will create several output files, all of them inside subdirectories in the same directory where it was run. The files that contain the physical information depend on what {{< octopus >}} is doing and they will be discussed in the next chapter. 

One directory that is always created is {{< file "exec/" >}}, this file contains information about the {{< octopus >}} run. Here you will find the {{< file "parser.log" >}} file, a text file that contains {{< emph "all" >}} the input variables that were read by {{< octopus >}}, both the variables that are in the input file and the variables that took default values; in the second case they are marked by a comment as {{< code "-default" >}}. This file is very useful if you want to check that {{< octopus >}} is correctly parsing a variable or what are the default values that it is taking.

### Clean Stop

{{< octopus >}} allows to stop the code gracefully after finishing the current iteration. These may not work for gcm, invert_ks, casida run modes. You can use this to prevent possible corruption of restart information when your job is killed by a scheduler, by preemptively asking {{< octopus >}} to quit automatically in a job script like this:

This can be done by creating a file called {{< file "stop" >}} in the running directory:
```text
 -PBS -l walltime=4:10:00
 mpirun $HOME/octopus/bin/octopus &> output &
 sleep 4h
 touch stop
 wait
```

or more sophisticatedly like this:

```text
 MIN=`qstat -a $PBS_JOBID | awk '{wall=$9} END {print $9}' | awk -F: '{min=($1*60)+$2; print min}'`
 sh ~/sleep_stop.sh stop $((MIN-10))m > sleepy &
 mpirun $HOME/octopus/bin/octopus &> output
```

with auxiliary script {{< code "sleep_stop.sh" >}}:

```text
 -!/bin/bash
 if [ $- -ne 2 ]; then
    echo "Usage: sleep_stop.sh FILENAME TIME"
 else
    echo "Time until $1 created: $2"
    rm -f $1
    sleep $2
    touch $1
 fi
```
{{< min-version 10 >}}
Since version 10, a second way is to use the input {{< variable Walltime >}} variable, which sets the time (in minutes) before which the code is gracefully stopped.
{{< /min-version >}}
### Restarting

Another directory that will created is {{< file "restart/" >}}; in this directory {{< octopus >}} saves the information from the calculation that it is doing. This information can be used in the following cases:

* If {{< octopus >}} is stopped without finishing by some reason, it can restart from where it was without having to do all work again.
* If after the calculation is done (or even if it was stopped), the user wants to do the same simulation with some different parameters, {{< octopus >}} can save some work by starting from the restart information.
* There are some calculations that require the results of other type of calculation as an input; in this case it uses the files written in {{< file "restart/" >}} by the previous calculation (we will discuss this case later, when we talk about the different calculation modes).

Sometimes it's not desired to restart a calculation, but to start it from the very beginning. {{< octopus >}} can be instructed to do so by setting the input variable {{< variable "fromScratch" >}} to {{< code "yes" >}}.

{{< manual-foot prev="Basics:Input file" next="Basics:Units" >}}
---
title: "Hamiltonian"
section: "Manual"
weight: 5
description: " "
manuals: "Basics"
---


{{< octopus >}} is based upon Density Functional Theory in the Kohn-Sham formulation. The Kohn-Sham Hamiltonian is the main part of this formulation; in this section we describe how the Hamiltonian is treated in {{< octopus >}} and what are the variables that control it.

## Introduction

Although the Hohenberg-Kohn theorem states that DFT is exact, the Kohn-Sham
method of reducing an interacting many-particle problem to a non-interacting single-particle problem introduces an approximation: the exchange-correlation term.

### Ground-State DFT

The Kohn-Sham method of DFT assumes that, for each interacting ground-state
density $n(r)$, there exists a non-interacting electron system with the same ground-state density. The interacting ground state is obtainable through the solution of the
Kohn-Sham equations

$$
  \\left\[-\\frac{\\nabla^2}{2}+v\_{\\rm KS}\[n\](r)\\right\]\\varphi\_i(r) = \\varepsilon\_i \\varphi\_i(r)
$$

The notation $v_{\rm KS}[n]$ means that the Kohn-Sham potential, 
$v_{\rm KS}$, has a functional dependence on $n$,
the electronic density, which is defined in terms of the Kohn-Sham wave-functions by

$$
  n(r) = \\sum\_i^{\\rm occ} | \\varphi\_i(r) |^2\\;.
$$

The potential $v_{\rm KS}$ is defined as the sum of the external potential (normally the 
potential generated by the nuclei), the Hartree term, and the exchange-correlation (xc) potential

$$
  v\_{\\rm KS}\[n\](r) = v\_{\\rm ext}(r) + v\_{\\rm Hartree}\[n\](r) + v\_{\\rm xc}\[n\](r)\\;.
$$

Due to the functional dependence on the density, these equations form a set of nonlinear 
coupled equations. The standard procedure to solve it is iterating until 
self-consistency is achieved.

The total energy of the electronic system is given by

$$
E\_{\\rm elec}\[n\] = T\_s\[n\] + \\int\\!\\!d^3r\\;n(r)v\_{\\rm ext}(r) + E\_{\\rm Hartree}\[n\] + E\_{\\rm x}\[n\] + E\_{\\rm c}\[n\]
$$

where $T_s[n]$ is the non-interacting kinetic energy, $v_{\rm ext}$ is the external potential, and $E_{\rm x,c}[n]$ are the exchange (x) and correlation (c) energies. The second term is called the "external energy" in the code.
In practice, from the solution of the Kohn-Sham equations, we can evaluate the energy via the eigenvalues as

$$
E\_{\\rm elec}\[n\] = \\sum \\epsilon\_i - \\int\\!\\!d^3r\\;n(r)v\_{\\rm xc}(r) - E\_{\\rm Hartree}\[n\] + E\_{\\rm x}\[n\] + E\_{\\rm c}\[n\]
$$

where $v_{\rm xc}[n]$ is the exchange-correlation potential. To find the total energy of the entire system, we additionally include ion-ion interaction and ionic kinetic energy.

$E_{\rm xc} = E_{\rm x} + E_{\rm c}$ is an unknown object and includes all the non-trivial many-body eﬀects required to make KS theory exact. Several approximations to $E_{\rm xc}$ have been proposed. The most used is the local density approximation (LDA). In this approximation $E_{\rm xc}[n(r)]$ is taken to be the exchange and correlation energy of a homogeneous electron gas with density $n = n(r)$.
Although there exists an exact expression for the exchange energy in this model, the exact value of the correlation energy is known only in the limit of very high densities. Ceperley and Alder did a Monte Carlo simulation of the homogeneous electron gas at several densities. Several parameterizations of the correlation energy for any density were then obtained interpolating the Monte Carlo results. One particularly simple parameterization was proposed by Perdew and Zunger, and this option may be used by {{< octopus >}}. You can, of course, choose other xc functionals ({{< variable "XCFunctional" >}}), via the extensive library provided by {{< libxc >}}.

### Time-dependent DFT

Time-dependent density-functional theory (TDDFT) extends the basic ideas
of ground-state density-functional theory (DFT) to the treatment of
excitations and of more general
time-dependent phenomena. TDDFT can be viewed as an alternative formulation
of time-dependent quantum mechanics but, in contrast to the normal approach
that relies on wave-functions and on the many-body Schrödinger equation, its
basic variable is the one-body electron density, $n(r,t)$. The advantages are
clear: The many-body wave-function, a function in a $3N$-dimensional
space (where $N$ is the number of electrons in the system), is a very complex 
mathematical object, while the density is a simple function that depends solely
on the 3-dimensional vector $r$. The standard way to obtain $n(r,t)$ 
is with the help of a fictitious system of non-interacting electrons,
the Kohn-Sham system.
The final equations are simple to tackle numerically, and are routinely
solved for systems with a large number of atoms. These electrons
feel an effective potential, the time-dependent Kohn-Sham potential.
The exact form of this potential is unknown, and has therefore to be
approximated.

The time-dependent Kohn-Sham equations are

$$
  i \\frac{\\partial}{\\partial t}\\varphi\_i(r,t) = \\left\[
    -\\frac{\\nabla^2}{2} + v\_{\\rm KS}(r,t)
  \\right\]\\varphi\_i(r,t)
$$

The density of the interacting system can be obtained from the 
time-dependent Kohn-Sham orbitals

$$
  n(r, t) = \\sum\_i^{\\rm occ} \\left|\\varphi\_i(r,t)\\right|^2
  \\;.
$$

The time-dependent Kohn-Sham equations, having the form of a one-particle equation,
is fairly easy to solve numerically. We stress, however, that
the Kohn-Sham equation ''is not'' a mean-field approximation: If we
knew the exact Kohn-Sham potential, $v_{\rm KS}$, we would obtain the exact Kohn-Sham orbitals, and from these the correct density of the system.
The Kohn-Sham potential is conventionally separated in the following way

$$
  v\_{\\rm KS}(r,t) = v\_{\\rm ext}(r,t) + v\_{\\rm Hartree}(r,t) + v\_{\\rm xc}(r,t)
  \\,.
$$

The first term is again the external potential. The Hartree potential
accounts for the classical electrostatic interaction between the electrons

$$
  v\_{\\rm Hartree}(r,t) = \\int\\!\\!d^3r' \\frac{n(r,t)}{\\left|r-r'\\right|}
  \\;.
$$

The time-dependence of the exchange and correlation potential introduces the need
for an approximation beyond the one made in the time-independent case. The simplest method of obtaining a time-dependent xc potential consists in assuming that
the potential is the time-independent xc potential evaluated at the time-dependent
density, i.e.,

$$
  v^{\\rm adiabatic}\_{\\rm xc}(r, t) = 
  \\left.\\tilde v\_{\\rm xc}\[n\](r)\\right|\_{n=n(t)}
  \\;,
$$

This is called the adiabatic approximation. If the time-independent xc potential chosen is the LDA, then we obtain the so-called adiabatic local density approximation
(ALDA). This approximation gives remarkably good excitation energies but suﬀers
from the same problems as the LDA, most notably the exponential fall-oﬀ of the
xc potential. If a strong laser pushes the electrons to regions far from the nucleus,
ALDA should not be expected to give an accurate description of the system.
Other options for the time-dependent xc potential are orbital-dependent
potentials like the exact exchange functional (EXX) (usually in the Krieger-Li-Iafrate (KLI) approximation).

## Occupations

You can set the occupations of the orbitals by hand, using the variable {{< variable "Occupations" >}}.

## External Potential

You can add an external uniform and constant electric or magnetic field, set by the {{< variable "StaticElectricField" >}} or {{< variable "StaticMagneticField" >}}. If you want to add a more general potential, you can do it using a user-defined {{< variable "Species" >}}.

If you got the coordinates from a PDB file, you can add the potential generated by the point charges defined there by setting the {{< variable "ClassicalPotential" >}} variable to yes.

## Electron-Electron Interaction

You can neglect this term by setting {{< code-inline >}}{{< variable "TheoryLevel" >}} = independent_particles{{< /code-inline >}} variable. This implies that both the Hartree and exchange-correlation terms will not be calculated.

### Exchange and correlation potential
{{< octopus >}} does not distinguish, for the moment, between ground-state DFT xc functionals, and time-dependent DFT functionals. In other words, in all cases the *adiabatic* approximation is assumed. In mathematical terms, we may formalize this in the following way: let $n$ be the time-dependent density, a function living in the four-dimensional space-time world. We will call $n_t(\vec{r}) = n(\vec{r}, t)\,$, the electronic density at time $t\,$, a function in the three-dimensional space. An exchange and/or correlation *energy* functional $E^{\alpha}\,$ in ground state DFT may then be used to build an *adiabatic* exchange and/or correlation *action* functional in the following way:

$$
  A^{\\alpha}\[n\] = \\int\_{t\_0}^{t\_f}{\\rm d}\\tau E^{\\alpha}\[n\_\\tau\]\\,,
$$

The time-dependent potential functional is then:

$$
  v^{\\alpha}\[n\](\\vec{r},t) \\equiv { \\delta A^{\\alpha} \\over \\delta
  n(\\vec{r}, t)} = { \\delta E^{\\alpha} \\over \\delta n\_{t}(\\vec{r}) } =
  v^{\\alpha \\rm (GS)}\[n\_t\](\\vec{r})\\,.
$$

We use the distinct notation $v^{\\alpha}\[n\](\\vec{r},t)\\,$ and $v^{\alpha \rm (GS)}\[n_t]\(\vec{r})\,$ to stress that the exchange and correlation potential in TDDFT -- the former -- and the exchange and correlation potential in GS-DFT -- the latter -- are in principle two conceptually different objects, which coincide only thanks to the adiabatic approximation. This is the reason why
we may actually only refer to the functionals in the GS-DFT context.

We may classify the xc functionals contained in the {{< octopus >}} code following John Perdew's Jacob's Ladder scheme:

* LDA rung: A functional $\alpha\,$ belonging to the LDA rung depends only on the electronic density (on the spin density in spin-polarized or spinors cases). Moreover, it has a <i>local</i> dependency on the density, i.e.:

$$
  E^{\\alpha}\_{\\rm LDA} = E^{\\alpha}\_{\\rm LDA}\[n\] = \\int {\\rm d}^3r e^{\\alpha}\_{\\rm LDA}(n(\\vec{r}))\\,.
$$

The potential may be then derived by functional derivation:

$$
  v^{\\alpha}\_{LDA}\[n\](\\vec{r}) = 
  {\\delta {E^{\\alpha}\_{\\rm LDA}} \\over {\\delta n(\\vec{r})}} 
  = {{{\\rm d}e^{\\alpha}\_{\\rm LDA} } \\over {{\\rm d}n }}(n(\\vec{r}))\\,.
$$

* GGA rung: A functional $\alpha\,$ belonging to the GGA rung depends on the electronic density, and also on its gradient. Moreover, it also has a *local* dependency (actually, the GGA is very often called a ''semi-local'' functional due to this).

$$
  E^{\\rm \\alpha (GGA)} = E^{\\rm \\alpha (GGA)}\[n, \\vec{\\nabla}n\] =
  \\int {\\rm d}^3r e^{\\alpha}\_{\\rm GGA}(n(\\vec{r}), \\vec{\\nabla}n(\\vec{r}))\\,.
$$

* meta-GGA rung:

* OEP rung: This is the family of functionals which are defined in terms of the occupied Kohn-Sham orbitals, $\lbrace \varphi_i \rbrace_{i=1}^{N}\,$. These are in fact the only *non-local* functionals. The name of the rung, OEP, stands for "optimized-effective-potential", the reason being that in general this is the method used to derive the potential from the energy functional (direct functional derivation is in this case not possible). A more suitable name would be orbital-dependent functionals.


{{< octopus >}} comes with several Exchange and Correlation potentials, including several flavours of LDA, GGA and OEP. You can choose your favorites by setting the variable {{< variable "XCFunctional" >}}. (Note that until {{< octopus >}} <= 2.1 the exchange-correlation functional was chosen with the variables <tt>XFunctional</tt> and <tt>CFunctional</tt>.)

When using OEP Exchange, the variable {{< variable "OEPLevel" >}} controls the level of approximation required.

You can also include Self-Interaction Correction (SIC), controlled by the variable {{< variable "SICCorrection" >}}.

## Relativistic Corrections

The variable {{< variable "RelativisticCorrection" >}} allows one to choose the relativistic correction to be used. Up to now only spin-orbit coupling ({{< variable "RelativisticCorrection" >}} = {{< code "spin_orbit" >}}) is implemented.


### Spin-orbit Coupling 
The spin-orbit coupling as it is implemented in {{< octopus >}} is included in the pseudo-potentials. These can either be HGH pseudo-potentials or Troullier-Martins-like pseudo-potentials. In the latter case the pseudo-potentials need to be generated from fully relativistic calculations and their fully separable form is given by:

$$
  \\hat{v}\_{KB}=v\_{local}+\\sum\_{l,j,m\_j} \\frac{|\\varphi^{PP}\_{ljm\_j}\\delta v\_{lj}^{PP}\>\<\\varphi^{PP}\_{ljm\_j}\\delta v\_{lj}^{PP}|}{\<\\varphi^{PP}\_{ljm\_j}| \\delta v\_{lj}^{PP}|\\varphi^{PP}\_{ljm\_j}\>}
$$

Since the angular part of the pseudo wave-functions $\varphi^{PP}_{ljm_j} $ are spherical spinors the wave-functions should be complex spinors and so the {{< variable "SpinComponents" >}} needs to be set to {{< code "non_collinear" >}}. This is also true for HGH pseudo-potentials.

Note that currently {{< octopus >}} is only able to read j-dependent Troullier-Martins-like pseudo-potentials that are provided in the UPF file format.


{{< manual-foot prev="Basics:Physical System" next="Basics:Discretization" >}}
---
Title: "Basics"
Weight: 10
Description: " "
section: "Manual"
---

{{% children depth=3 %}}---
title: "Physical System"
section: "Manual"
weight: 4
description: " "
---


The first thing that {{< octopus >}} has to know is the physical system you want to treat. To do this you have specify a group of species and their positions. 

## Dimensions

{{< octopus >}} can work in a space with 1, 2 or 3 dimensions. You can select the dimension of your system with the {{< variable "Dimensions" >}} variable.

## Species

An {{< octopus >}} species is very generic and can be a nucleus (represented by pseudopotentials or by the full Coulomb potential), a jellium sphere or even a user-defined potential. The information regarding the atomic species goes into the {{< variable "Species" >}} block.

### Pseudopotentials
The many-electron Schroedinger equation can be greatly simpliﬁed if electrons are
divided in two groups: valence electrons and inner core electrons. The electrons
in the inner shells are strongly bound and do not play a signiﬁcant role in the
chemical binding of atoms, thus forming with the nucleus an inert core. Binding
properties are almost completely due to the valence electrons, especially in metals
and semiconductors. This separation implies that inner electrons can be ignored,
reducing the atom to an inert ionic core that interacts with the valence electrons.
This suggests the use of an eﬀective interaction, a pseudopotential, that gives an
approximation to the potential felt by the valence electrons due to the nucleus and
the core electrons. This can signiﬁcantly reduce the number of electrons that have
to be dealt with. Moreover, the pseudo wave functions of these valence electrons
are much smoother in the core region than the true valence wave functions, thus
reducing the computational burden of the calculations.

Modern pseudopotentials are obtained by inverting the free atom Schroedinger equation 
for a given reference electronic conﬁguration, and forcing the pseudo wave
functions to coincide with the true valence wave functions beyond a certain cutoﬀ
distance. The pseudo wave functions are also forced to have the same norm as the
true valence wave functions, and the energy pseudo eigenvalues are matched to the
true valence eigenvalues. Diﬀerent methods of obtaining a pseudo eigenfunction
that satisﬁes all these requirements lead to diﬀerent non-local, angular momentum
dependent pseudopotentials. Some widely used pseudopotentials are the Troullier
and Martins potentials, the Hamann potentials, the Vanderbilt potentials and the 
Hartwigsen-Goedecker-Hutter potentials. The default potentials used by {{< octopus >}} 
are of the Troullier and Martins type, although you can also opt for the HGH potentials.

{{< octopus >}} comes with a package of pseudopotentials and the parameters needed to use them. If you want to have a look you can find them under {{< inst_file "share/octopus/pseudopotentials" >}}. These pseudopotentials serve to define many species that you might wish to use in your coordinates block, e.g. a helium atom, "He". If you are happy to use these predefined pseudopotentials, you do not need to write a species block.

However it is also possible to define new species to use in your coordinates block by adding a {{< variable "Species" >}} block. You can check the documentation of that variable for the specific syntax. With this block you may specify the format of your file that contains a pseudopotential and parameters such as the atomic number, or you may define an algebraic expression for the potential with the user-defined potential. A user defined potential should be finite everywhere in the region where your calculation runs.

If you want to search other repositories on the web or create your own pseudopotentials, check the {{< manual Pseudopotentials pseudopotentials>}} page. Save the pseudopotential file in the same directory as the inp file and specify its format with the species block.

### All-Electron Nucleus

The potential of this species is the full Coulomb potential 

$$
V\\left( \\vec{r}\\right)=\\frac1{\\left|\\vec{R}\_i-\\vec{r}\\right|}\\ .
$$

The main problem to represent this potential is the discontinuity over $\vec{R}_i$. To overcome this problem we do the following: 

* First we assume that atoms are located over the closest grid point.
* Then we calculate the charge density associated with the nucleus: a delta distribution with the value of the charge at this point and zero elsewhere.
* Now we solve the Poisson equation for this density.

In this way we get a potential that is the best representation of the Coulomb potential for our grid (we will discuss about grids later) and is continuous in $\vec{R}_i$ (the value is the average of the potential over the volume associated with the grid point).

The main problem is that the requirement of having atoms over grid points is quite strong: it is only possible for a few systems with simple geometries and you can't move the atoms. Also the Coulomb potential is very hard, which means you will need a very small spacing, and as you have to consider both core and valence electrons, this species is only suitable for atoms or very small molecules.

### User Defined

It is also possible to define an external, user-defined potential in the input file. All {{< versioned-link "manual/basics/input_file#mathematical-expressions" "functions" >}} accepted by the parser can be used. Besides that, one can use the symbols $x$, $y$, $z$, and $r$. In this way it is trivial to calculate model systems, like harmonic oscillators, quantum dots, etc.

## Coordinates

For each instance of a species (even for user-defined potentials), you have to specify its position inside the simulation box. To do this you can use the {{< variable "Coordinates" >}} block which describes the positions inside of the input file or one of the {{< variable "XYZCoordinates" >}} or {{< variable "PDBCoordinates" >}} variables, that specify an external file, in {{< name "xyz" >}} or {{< name "PDB" >}} format respectively, from where the coordinates will be read.

Before using a geometry with {{< octopus >}} we recommend that you center it. For this you can use the {{< manual "External utilities:oct-center-geom" "oct-center-geom" >}} utility.

## Velocities

If you are going to do ion dynamics you may want to have an initial velocity for the particles. You have several choices for doing this:

* Don't put anything in the input file; particles will have zero initial velocity.
* Give them a random velocity according to a temperature (in degrees Kelvin) given by the {{< variable "RandomVelocityTemp" >}} variable.
* Explicitly give the initial velocity for each particle, either through the {{< variable "Velocities" >}} block or from a pseudo-xyz file detailed by the variable {{< variable "XYZVelocities" >}}.

## Number of Electrons

Each species adds enough electrons to make the system neutral. If you want to add or remove electrons you can specify the total charge of your system with the {{< variable "ExcessCharge" >}} variable (a negative charge implies to add electrons).

{{< manual-foot prev="Basics:Units" next="Basics:Hamiltonian" >}}

---
title: "Troubleshooting"
section: "Manual"
weight: 9
description: " "
---


If {{< octopus >}} works properly on your system (i.e. you can recreate the results in the {{< versioned-link "/tutorials" "tutorials" >}} but you have troubles using it for your own work, here are some things to try. Please feel free to add your own ideas here.

#### Read parser.log 
If you look at the file {{< file "exec/parser.log" >}}, it will tell you the value of the variables that you set with the {{< file "inp" >}} file, as well as all the variables which are taking their default value. This can sometimes be helpful in understanding the behavior of the program.

#### Use OutputDuringSCF 
If you add {{< variable "OutputDuringSCF" >}}{{< code " = yes" >}} to your input file, you can examine the results of each iteration in the Self Consistent Field calculation. So if you also have the variable {{< variable "Output" >}} set to {{< code "Output = density + potential" >}}, both the electron density and the Kohn-Sham, bare, exchange-correlation and Hartree potentials will be written to a folder (called, e.g., {{< file "scf.0001" >}}) after each SCF iteration.

#### Set Debug 
Set the variable {{< variable "Debug" >}} to {{< code info >}} for some extra diagnostic info and a stack trace with any fatal error, {{< code trace >}} to add a full stack strace, and {{< code trace_file >}} to get a stack trace from each MPI task when running in parallel.

{{< manual-foot prev="Basics:Symmetry" next="Calculations:Ground State" >}}
---
title: "Output"
section: "Manual"
weight: 7
description: " "
---


At first you may be quite happy that you have mastered the input file, and {{< octopus >}} runs without errors. However, eventually you (or your thesis advisor) will want to learn something about the system you have managed to describe to {{< octopus >}}.

#### Ground-State DFT 
{{< octopus >}} sends some relevant information to the standard output (which you may have redirected to a file). Here you will see energies and occupations of the eigenstates of your system. These values and other information can also be found in the file {{< file "static/info" >}}.

However {{< octopus >}} also calculates the wavefunctions of these states and the positions of the nuclei in your system. Thus it can tell you the density of the dipole moment, the charge density, or the matrix elements of the dipole moment operator between different states. Look at the values that the {{< variable "Output" >}} variable can take to see the possibilities.

For example, if you include 
{{< code-line "Output = wfs_sqmod + potential" >}} 

in your {{< file "inp" >}} file, {{< octopus >}} will create separate text files in the directory {{< file "static" >}} with the values of the square modulus of the wave function and the local, classical, Hartree, and exchange/correlation parts of the Kohn-Sham potential at the points in your mesh.

You can specify the formatting details for these input files with the {{< variable "OutputFormat" >}} variable and the other variables in the {{< variable "Output" >}} section of the Reference Manual. For example, you can specify that the file will only contain values along the x, y, or z axis, or in the plane x=0, y=0, or z=0. You can also set the format to be readable by the graphics programs [OpenDX](https://en.wikipedia.org/wiki/IBM_OpenDX)[^opendx], [gnuplot](https://www.gnuplot.info/) or [MatLab](https://www.mathworks.com/products/matlab/). OpenDX can make plots of iso-surfaces if you have data in three-dimensions. However gnuplot can only make a 3-d plot of a function of two variables, i.e. if you have the values of a wavefunction in a plane, and 2-d plots of a function of one variable, i.e. the value of the wavefunction along an axis.

[^opendx]: Unfortunately, OpenDX is no longer actively maintained. Some Linux distributions, e.g. Debian or Ubuntu still have opendx packages, but we recommend to move to different visualization programmes, such as [visit](https://wci.llnl.gov/simulation/computer-codes/visit) or [paraview](https://www.paraview.org/).

#### Time-Dependent DFT 

##### Optical Properties 

A primary reason for using a time-dependent DFT program is to obtain the optical properties of your system.  You have two choices for this, linear-response theory <i>a la</i> [Jamorski, Casida & Salahub](https://dx.doi.org/10.1063/1.471140)[^footnote-1], 
or explicit time-propagation of the system after a perturbation, a la [Yabana & Bertsch](https://prola.aps.org/abstract/PRB/v54/i7/p4484_1)[^footnote-2]. 
You may wish to read more about these methods in the paper by [Castro et al.](https://prola.aps.org/abstract/PRL/v76/i8/p1212_1)[^footnote-3].

##### Linear-Response Theory 
Linear-response theory is based on the idea that a small (time-dependent) perturbation in an externally applied electric potential 
$\delta v (r, \omega )$ will result in a (time-dependent) perturbation of the electronic density 
$ \delta \rho (r, \omega )$ which is linearly related to the size of the perturbation:
$\delta \rho (r, \omega ) = \int d^{3} r' \chi (r, r'; \omega) \delta v (r, \omega )$. Here, obviously, the time-dependence is Fourier-transformed into a frequency-dependence, $ \omega $. The susceptibility,
$\chi(r, r'; \omega) $, is a density-density response function, because it is the response of the charge density to a potential that couples to the charge density of the system. Because of this, it has poles at the excitation energies of the many-body system, meaning that the induced density also has these poles. One can use this analytical property to find a related operator whose eigenvalues are these many-body excitation energies. The matrix elements of the operator contain among other things: 1) occupied and unoccupied Kohn-Sham states and energies (from a ground state DFT calculation) and 2) an exchange-correlation kernel, 
$ f_{xc}(r, r', \omega) = {{\delta v_{xc}[n(r,\omega)]}\over{\delta n(r',\omega)}}\mid_{\delta v_{ext} = 0}$.

Casida's equations are a full solution to this problem (for real wavefunctions). The Tamm-Dancoff approximation uses only occupied-unoccupied transitions. The Petersilka approximation uses only the diagonal elements of the Tamm-Dancoff matrix, <i>i.e.</i> there is no mixing of transitions.
[^footnote-4] It takes only a little more time to calculate the whole matrix, so Petersilka is provided mostly for comparison.

These methods are clearly much faster (an order of magnitude) than propagating
in time, but it turns out that they are very sensitive to the quality
of the unoccupied states. This means that it is very hard to converge the
excitation energy, because one requires a very large simulation box (much
larger than when propagating in real time).

##### Electronic Excitations by Means of Time-Propagation 

See [Time-Dependent: Delta-kick: Calculating an absorption spectrum](../../calculations/time-dependent#delta-kick-calculating-an-absorption-spectrum)

{{< manual-foot prev="Basics:Discretization" next="Basics:Troubleshooting" >}}

#### References 

[^footnote-1]: {{< article title="Dynamic polarizabilities and excitation spectra from a molecular implementation of time-dependent density-functional response theory: N2 as a case study" authors="Christine Jamorski, Mark E. Casida, and Dennis R. Salahub " journal="J. Chem. Phys." year="1996" vol="104" issue="13" pages="5134-5147" url="https://scitation.aip.org/getabs/servlet/GetabsServlet?prog=normal&id=JCPSA6000104000013005134000001&idtype=cvips&gifs=yes" doi="10.1063/1.471140 " >}}

[^footnote-2]: {{< article title="Time-dependent local-density approximation in real time " authors="K. Yabana, G. F. Bertsch" journal="Phys. Rev. B" year="1996" vol="54" issue="7" pages="4484 - 4487" url="https://link.aps.org/abstract/PRB/v54/p4484" doi="10.1103/PhysRevB.54.4484" >}}

[^footnote-3]: {{< article title="octopus: a tool for the application of time-dependent density functional theory" authors="Alberto Castro, Heiko Appel, Micael Oliveira, Carlo A. Rozzi, Xavier Andrade, Florian Lorenzen, M. A. L. Marques, E. K. U. Gross, Angel Rubio " journal="physica status solidi (b)" year="2006" volume="243" issue="11" pages="2465-2488" url="https://www3.interscience.wiley.com/cgi-bin/abstract/112650991/ABSTRACT?CRETRY=1&SRETRY=0" doi="10.1002/pssb.200642067" >}}

[^footnote-4]: {{< article title="Excitation Energies from Time-Dependent Density-Functional Theory" authors="Petersilka, M. and Gossmann, U. J. and Gross, E. K. U." journal="Phys. Rev. Lett." volume="76" number="8" pages="1212--1215" numpages="3" year="1996" month="Feb" doi="10.1103/PhysRevLett.76.1212" publisher="American Physical Society" url="https://prola.aps.org/abstract/PRL/v76/i8/p1212_1" >}}



---
title: "Units"
section: "Manual"
weight: 3
description: " "
---


Before entering into the physics in {{< octopus >}} we have to address a very important issue: {{< emph "units" >}}. There are different unit systems that can be used at the atomic scale: the most used ones are atomic units and what we call "convenient" units. Here we present both unit systems and explain how to use them in {{< octopus >}}.

### Atomic Units

Atomic units are a Gaussian system of units (by "Gaussian" it means that the
vacuum dielectric constant has no dimensions and is set to be
$\epsilon_0 = {1 \over {4\pi}}$),
in which the numerical values of the Bohr radius, the electronic
charge, the electronic mass, and the reduced Planck's constant are set to one:

$$
(1)\\qquad a\_0 = 1; e^2 = 1; m\_e = 1; \\hbar = 1.
$$

This simplifies formulae (although some may feel it presents a serious hazard for dimensional analysis,
interpretation and understanding of formulae, and physics in general.
But this is just a personal taste).
This sets directly two fundamental units, the atomic units of length and of mass:

$$
(2)\\qquad {\\rm au}\_{\\rm length} = a\_0 = 5.2917721\\times 10^{-11}~{\\rm m};\\quad 
{\\rm au}\_{\\rm mass} = m\_e = 9.1093819\\times 10^{-31}~{\\rm kg}.
$$

Since the squared charge must have units of energy times length, we can thus
set the atomic unit of energy

$$
(3)\\qquad {\\rm au}\_{\\rm energy} = {e^2 \\over a\_0} = 4.3597438\\times 10^{-18}~{\\rm J},
$$

which is called Hartree, {{< hartree >}}. And, since the energy has units of mass times
length squared per time squared, this helps us get the atomic unit of time:

$$
(4)\\qquad {\\rm Ha} = m\_e { a\_0^2 \\over {\\rm} {\\rm au}\_{\\rm time}^2} \\to 
{\\rm au}\_{\\rm time} = a\_0 \\sqrt{m\_e \\over {\\rm Ha}} = {a\_0 \\over e} \\sqrt{m\_e a\_0}
= 2.4188843\\times 10^{-17}~{\\rm s}.
$$

Now the catch is: what about Planck's constant? Its dimensions are of energy
times time, and thus we should be able to derive its value by now. But at the
beginning we set it to one! The point is that the four physics constants
used ($a_0, m_e, e^2, \hbar$) are not independent, since:

$$
(5)\\qquad a\_0 = { \\hbar^2 \\over {m\_e \\; {e^2 \\over {4 \\pi \\epsilon\_0} } } }.
$$

In this way, we could actually have derived the atomic unit of time in an
easier way, using Planck's constant:

$$
(6)\\qquad \\hbar = 1\\; {\\rm Ha}\\,{\\rm au}\_{\\rm time} \\Rightarrow {\\rm au}\_{\\rm time} = { \\hbar \\over {\\rm Ha}} = 
{ {\\hbar a\_0} \\over e^2}\\,.
$$

And combining (6) and (5) we retrieve (4).

### Convenient Units

Much of the literature in this field is written using Ångströms and electronvolts
as the units of length and of energy, respectively. So it may be "convenient"
to define a system of units, derived from the atomic system of units, in which
we make that substitution. And so we will call it "convenient".

The unit mass remains the same, and thus the unit of time must change, being
now $\hbar /{\rm eV}\,$,
with $\hbar = 6.582\,1220(20)\times 10^{-16}~\rm eV\,s$.

### Units in {{< octopus >}}

Except where otherwise noted, {{< octopus >}} expects all values in the input file to be in atomic units. If you prefer to use other units in the input file, the code provides some handy conversion factors. For example, to write some length value in Ångströms, you can simply multiply the value by {{< code "angstrom" >}}:

{{< code-line "Spacing = 0.5*angstrom" >}}
 
A complete list of units {{< octopus >}} knows about can be found in the {{< variable "Units" >}} variable description. 

By default {{< octopus >}} writes all values atomic units. You can switch to convenient units by setting the variable {{< variable "UnitsOutput" >}} to {{< value "ev_angstrom" >}}.

### Mass Units

An exception for units in {{< octopus >}} is mass units. When dealing with the mass of ions, [atomic mass units (amu)](https://en.wikipedia.org/wiki/Amu) are always used. This unit is defined as $1/12$ of the mass of the <sup>12</sup>{{< name "C" >}} atom. In keeping with standard conventions in solid-state physics, effective masses of electrons are always reported in units of the electron mass (''i.e.'' the atomic unit of mass), even in the eV-Å system.

### Charge Units

In both unit systems, the charge unit is the electron charge ''e'' (''i.e.'' the atomic unit of charge).

### Unit Conversions

Converting units can be a very time-consuming and error-prone task when done by hand, especially when there are implicit constants set to one, as in the case of atomic units. That is why it's better to use as specialized software like [GNU Units](https://www.gnu.org/software/units/units.html).

In some fields, a very common unit to express the absorption spectrum is Mb. To convert a strength function from 1/eV to Mb, multiply by $\pi h c r_e\,$, with $r_e=e^2/(m_e c^2)\,$. The numerical factor is 109.7609735.


{{< manual-foot prev="Basics:Running Octopus" next="Basics:Physical System" >}}

---
title: "Symmetry"
section: "Manual"
weight: 8
description: " "
---


There is not much use of symmetry in Octopus.

In finite systems, you will get an analysis just for your information, but which will not be used anywhere in the code. It is of this form (''e.g.'' for silane):

```text
 ***************************** Symmetries *****************************
 Symmetry elements : 4*(C3) 3*(C2) 3*(S4) 6*(sigma)
 Symmetry group    : Td
 **********************************************************************
```

Many symmetries will in fact be broken by the real-space mesh. Since it is always orthogonal, it will break three-fold rotational symmetry of benzene, for example.

In periodic systems, you will also get an analysis, _e.g._ like this for bulk silicon in its 8-atom convention cell:

```text
 ***************************** Symmetries *****************************
 Space group No.227
  International: Fd -3 m
  International(long): Fd -3 m _1
  Schoenflies: Oh^7
  Multiplicity: 192
 Point group
  International: m -3 m
  Schoenflies: Oh
 Identity has a fractional translation     0.500000    0.500000    0.000000
 Identity has a fractional translation     0.500000    0.000000    0.500000
 Identity has a fractional translation     0.000000    0.500000    0.500000
 Disabling fractional translations. System appears to be a supercell.
 Info: The system has    24 symmetries that can be used.
 **********************************************************************
```

The analysis is done by the library [spglib](https://spglib.github.io/spglib/). The comments on fractional translations and supercell are due to the use of the conventional cell, since the 2-atom primitive cell does not have orthogonal lattice vectors and therefore cannot be used in Octopus currently.

For large systems, the symmetry analysis might be very time-consuming; in rare cases, the symmetry analysis might crash and stop the calculation. In either situation, you can use the variable {{< variable "SymmetriesCompute" >}} to turn off the symmetry analysis. You can manually identify a direction in which symmetry is broken with {{< variable "SymmetryBreakDir" >}}, which is appropriate for the case that an external field is applied in a particular direction.

Symmetries are used in two ways for periodic systems: first, to reduce the set of k-points needed. This behavior is controlled by {{< variable "KPointsUseSymmetries" >}}. For high-symmetry systems, this can make a dramatic reduction in the time required for a calculation. For example, in our silicon example, if we set

```text
 %KPointsGrid
  4   | 4   | 4
  0.5 | 0.5 | 0.5
 %
```

we will obtain not 4×4×4 = 64 k-points but only 4:

```text
    4 k-points generated from parameters :
  ---------------------------------------------------
     n =    4    4    4      s =  0.50  0.50  0.50
 
  index |    weight    |             coordinates              |
      1 |     0.125000 |    0.125000    0.125000    0.125000  |
      2 |     0.375000 |    0.125000    0.125000    0.375000  |
      3 |     0.375000 |    0.375000    0.375000    0.125000  |
      4 |     0.125000 |    0.375000    0.375000    0.375000  |
```

The density and other scalar quantities are straightforwardly computed from the density due to each k-point and the weight; the same has been done for some vector quantities such as the gradients used in GGA functionals and the forces.. Note, however, that the proper calculation of more complicated tensorial properties such as the response functions in the Sternheimer calculation modes have not been implemented with symmetry. You can also use the symmetry operations to symmetrize the density, via {{< variable "SymmetrizeDensity" >}}. In general, you should beware of use of symmetries for partially periodic systems (_e.g._ wire or sheet geometry) for which some problems have been found.

The use of symmetry for reducing the number of calculations needed for obtaining absorption spectra by time-propagation is discussed here:  
{{< tutorial "response/use_of_symmetries_in_optical_spectra_from_time-propagation" "Use of symmetries in optical spectra from time-propagation" >}}


{{< manual-foot prev="Basics:Output" next="Basics:Troubleshooting" >}}
---
title: "Input file"
section: "Manual"
weight: 1
description: " "
---


{{< octopus >}} uses a single input file from which to read user instructions to know what to calculate and how. This page explains how to generate that file and what is the general format. The Octopus parser is a library found in the {{< file "liboct_parser" >}} directory of the source, based on bison and C. You can find two (old) separate release versions of it at the bottom of the {{< octopus-releases >}} page.

### Input file

Input options should be in a file called {{< file "inp" >}}, in the directory {{< octopus >}} is run from. This is a plain {{< name "ASCII" >}} text file, to create or edit it you can use any text editor like {{< name "emacs" >}}, {{< name "vi" >}}, {{< name "jed" >}}, {{< name "pico" >}}, {{< name "gedit" >}}, etc. For a fairly comprehensive example, just look at the tutorial page {{< tutorial "Basics/Basic_input_options" "Basic input options" >}}.

At the beginning of the program, the parser reads the input file, parses it, and generates a list of variables that will be read by {{< octopus >}} (note that the input is case-independent). There are two kind of variables: scalar values (strings or numbers), and blocks (that you may view as matrices).

### Scalar Variables

A scalar variable {{< code "var" >}} can be defined by:

{{< code-line " var = exp" >}}

{{< code "var" >}} can contain any alphanumeric character plus _, and {{< code "exp" >}} can be a quote-delimited string, a number (integer, real, or complex), a variable name, or a mathematical expression. Complex numbers are defined as <tt>{real, imag}</tt>. Real numbers can use scientific notation with <tt>e</tt> or <tt>E</tt> (no <tt>d</tt> or <tt>D</tt>, Fortran people), such as <tt>6.02e23</tt>. Variable names are not case-sensitive, and you must not redefine a previously defined symbol -- especially not the reserved variables <tt>x, y, z, r, w, t</tt> which are used in space- or time-dependent expressions, where <tt>w</tt> is the 4th space coordinate when operating in 4D. A list of predefined variables follows [below](#predefined-variables)

### Mathematical expressions

The parser can interpret expressions in the input file, either to assign the result to a variable, or for defining functions such as a potential in the {{< variable "Species" >}} block or a time-dependent function in the {{< variable "TDFunctions" >}} block. The arguments can be numbers or other variables.

#### Arithmetic operators
* {{< code "a+b" >}}: addition
* {{< code "a-b" >}}: subtraction
* {{< code "-a" >}}: unary minus
* {{< code "a*b" >}}: multiplication
* {{< code "a/b" >}}: division
* {{< code "a^b" >}}: exponentiation

#### Logical operators

Logical operation will return 0 for false or 1 for true. You can exploit this to define a piecewise expression, e.g. <tt>"2 * (x <= 0) - 3 * (x > 0)"</tt> (although the <tt>step</tt> function may also be used). The comparison operators (except <tt>==</tt>) use only the real part of complex numbers.

* {{< code "a < b" >}}: less than
* {{< code "a <= b" >}}: less than or equal to ($\le$)
* {{< code "a > b" >}}: greater than
* {{< code "a >= b" >}}: greater than or equal to ($\ge$)
* {{< code "a == b" >}}: equal to
* {{< code "a && b" >}}: logical and
* {{< code "a || b" >}}: logical or
* {{< code "!a" >}}: logical not

#### Functions

* {{< code "sqrt(x)" >}}: The square root of {{< code "x" >}}.
* {{< code "exp(x)" >}}: The exponential of {{< code "x" >}}.
* {{< code "log(x)" >}} or {{< code "ln(x)" >}}: The natural logarithm of {{< code "x" >}}.
* {{< code "log10(x)" >}}: Base 10 logarithm of {{< code "x" >}}.
* {{< code "logb(x, b)" >}}: Base {{< code "b" >}} logarithm of {{< code "x" >}}.
* {{< code "{x, y&-125;" >}}: The complex number $x + iy$.
* {{< code "arg(z)" >}}: Argument of the complex number {{< code "z" >}}, $\arg(z)$, where $-\pi < \arg(z) <= \pi$.
* {{< code "abs(z)" >}}: Magnitude of the complex number {{< code "z" >}}, $|z|$.
* {{< code "abs2(z)" >}}: Magnitude squared of the complex number {{< code "z" >}}, $|z|^2$.
* {{< code "logabs(z)" >}}: Natural logarithm of the magnitude of the complex number {{< code "z" >}}, $\log|z|$. It allows an accurate evaluation of $\log|z|$ when $|z|$ is close to one. The direct evaluation of {{< code "log(abs(z))" >}} would lead to a loss of precision in this case.
* {{< code "conjg(z)" >}}: Complex conjugate of the complex number {{< code "z" >}}, $z^* = x - i y$.
* {{< code "inv(z)" >}}: Inverse, or reciprocal, of the complex number {{< code "z" >}}, $\frac{1}{z} = \frac{x - i y}{x^2 + y^2}$.
* {{< code "sin(x)" >}}, {{< code "cos(x)" >}}, {{< code "tan(x)" >}}, {{< code "cot(x)" >}}, {{< code "sec(x)" >}}, {{< code "csc(x)" >}}: The sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
* {{< code "asin(x)" >}}, {{< code "acos(x)" >}}, {{< code "atan(x)" >}}, {{< code "acot(x)" >}}, {{< code "asec(x)" >}}, {{< code "acsc(x)" >}}: The inverse (arc-) sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
* {{< code "atan2(x,y)" >}}: = $\mathrm{atan}(y/x)$.
* {{< code "sinh(x)" >}}, {{< code "cosh(x)" >}}, {{< code "tanh(x)" >}}, {{< code "coth(x)" >}}, {{< code "sech(x)" >}}, {{< code "csch(x)" >}}: The hyperbolic sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.
* {{< code "asinh(x)" >}}, {{< code "acosh(x)" >}}, {{< code "atanh(x)" >}}, {{< code "acoth(x)" >}}, {{< code "asech(x)" >}}, {{< code "acsch(x)" >}}: The inverse hyperbolic sine, cosine, tangent, cotangent, secant and cosecant of {{< code "x" >}}.

* {{< code "min(x, y)" >}}: The minimum of {{< code "x" >}} and {{< code "y" >}}.
* {{< code "max(x, y)" >}}: The maximum of {{< code "x" >}} and {{< code "y" >}}.
* {{< code "step(x)" >}}: The Heaviside step function in {{< code "x" >}}. This can be used for piecewise-defined functions.
* {{< code "erf(x)" >}}: The error function $\mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x dt e^{-t^2}$.

* {{< code "realpart(z)" >}}: The real part of the complex number {{< code "z" >}}.
* {{< code "imagpart(z)" >}}: The imaginary part of the complex number {{< code "z" >}}.
* {{< code "floor(x)" >}}: The largest integer less than the real number {{< code "x" >}}.
* {{< code "ceiling(x)" >}}: The smallest integer greater than the real number {{< code "x" >}}.

These mathematical operations are all based on the GSL library and are defined in {{< code "symbols.c" >}} and {{< code "grammar.y" >}}.

#### References
* https://www.gnu.org/software/gsl/manual/html_node/Properties-of-complex-numbers.html
* https://www.gnu.org/software/gsl/manual/html_node/Complex-arithmetic-operators.html
* https://www.gnu.org/software/gsl/manual/html_node/Error-Function.html
* https://www.gnu.org/software/libc/manual/html_node/Rounding-Functions.html
* https://www.gnu.org/software/gsl/manual/html_node/Representation-of-complex-numbers.html

### Predefined variables

There are some predefined constants for your convenience:

* {{< code "pi" >}}: {{< code "3.141592653589793" >}}.
* {{< code "e" >}}: The base of the natural logarithms.
* {{< code "false" >}} or {{< code "no" >}}: False.
* {{< code "true" >}} or {{< code "yes" >}}: True.
* {{< code "i" >}}: The imaginary unit $i$, ''i.e.'' {{< code "{0, 1" >}}}

Since version 7.0 there are also some predefined units that can be found by searching for the decimal point {{< code "." >}} in {{< code "share/variables" >}}:

* {{< code "angstrom" >}}: {{< code "1.8897261328856432" >}}.
* {{< code "pm" >}} or {{< code "picometer" >}}: {{< code "0.018897261328856432" >}}.
* {{< code "nm" >}} or {{< code "nanometer" >}}: {{< code "18.897261328856432" >}}.
* {{< code "ry" >}} or {{< code "rydberg" >}}: {{< code "0.5" >}}.
* {{< code "eV" >}} or {{< code "electronvolt" >}}: {{< code "0.03674932539796232" >}}.
* {{< code "invcm" >}}: {{< code "4.5563353e-06" >}}.
* {{< code "kelvin" >}}: {{< code " 3.1668105e-06" >}}.
* {{< code "kjoule_mol" >}}: {{< code " 0.00038087988" >}}.
* {{< code "kcal_mol" >}}: {{< code " 0.0015936014" >}}.
* {{< code "as" >}} or {{< code "attosecond" >}}: {{< code "0.0413413737896" >}}.
* {{< code "fs" >}} or {{< code "femtosecond" >}}: {{< code "41.3413737896" >}}.
* {{< code "ps" >}} or {{< code "picosecond" >}}: {{< code "41341.3737896" >}}.
* {{< code "c" >}}: {{< code "137.035999139" >}}.

### Blocks

Blocks are defined as a collection of values, organised in row and column format.
The syntax is the following:
{{< code-block >}}
%var
 exp | exp | exp | ... 
 exp | exp | exp | ... 
 ...
%
{{< /code-block >}}

Rows in a block are separated by a newline, while columns are separated by the character | or by a tab. There may be any number of lines and any number of columns in a block. Note also that each line can have a different number of columns. Values in a block don't have to be of the same type.

### Comments

Everything following the character {{< code "-" >}} until the end of the line is considered a comment and is simply cast into oblivion.

### Includes

With <code>include FILENAME</code> it is possible to include external files into the input file. To illustrate the usage of this command we can split the input file from the {{< versioned-link "tutorial/basics/total_energy_convergence#methane-molecule" "Methane tutorial" >}} in two.

In the input file {{< file "inp" >}} we have: 
```text
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom

 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.22*angstrom
 
 
 include geometry.oct
```
with the geometry being defined in {{< file "geometry.oct" >}}
```text
 CH = 1.2*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
```

### Environment variables


You can also {{< versioned-link "manual/advanced_ways_of_running_octopus/#passing-arguments-from-environment-variables" "set variables using the environment" >}}, which can be helpful in scripting.

### Default values

If {{< octopus >}} tries to read a variable that is not defined in the input file, it automatically assigns to it a default value (there are some cases where {{< octopus >}} cannot find a sensible default value and it will stop with an error). All variables read (present or not in the input file) are output to the file {{< file "exec/parser.log" >}}. The variable that are not defined in the input file will have a <tt>-default</tt> comment to it. If you are not sure of what the program is reading, just take a look at it.

We recommend you to keep the variables in the input file to a minimum: ''do not write a variable that will be assigned its default value''. The default can change in newer versions of {{< octopus >}} and old values might cause problems. Besides that, your input files become difficult to read and understand.

### Documentation

Each input variable has (or should have) its own documentation explaining what it does and the valid values it may take. This documentation can be obtained {{< versioned-link "/Variables" "online" >}} or it can also be accessed by the {{< manual "External_utilities:oct-help" "oct-help" >}} command.

### Experimental features

Even in the stable releases of {{< octopus >}} there are many features that are being developed and are not suitable for production runs. To protect users from inadvertly using these parts they are declared as ''Experimental''.

When you try to use one of these experimental functionalities {{< octopus >}} will stop with an error. If you want to use it you need to set the variable {{< variable "ExperimentalFeatures" >}} to {{< code "yes" >}}. Now {{< octopus >}} will only emit a warning.

By setting {{< variable "ExperimentalFeatures" >}} to {{< code "yes" >}} you will be allowed to use parts of the code that are not complete or not well tested and most likely produce wrong results. If you want to use them for production runs you should contact the Octopus developers first.

### Good practices

In order to ensure compatibility with newer versions of {{< octopus >}} and avoid problems, keep in mind the following rules of good practice when writing input files:

* Although input variables that take an option as an input can also take a number, the number representation makes the input file less readable and it is likely to change in the future. So '''avoid using numbers instead of values'''. For example
```text
 UnitsOutput = ev_angstrom
```
''must always'' be used instead of
```text
 UnitsOutput = 3
```
* '''Do not include variables that are not required in the input file''', especially declarations of values that are just a copy of the default value. This makes the input file longer, less readable and, since defaults are likely to change, it makes more probable that your input file will have problems with newer versions of {{< octopus >}}. Instead rely on default values.

* '''Avoid duplicating information in the input file'''. Use your own variables and the mathematical-interpretation capabilities for that. For example, you should use:
```text
 m = 0.1
 c = 137.036
 E = m*c^2
```
instead of
```text
 m = 0.1
 c = 137.036
 E = 1877.8865
```
In the second case, you might change the value of {{< code "m" >}} (or {{< code "c" >}} if you are a cosmologist) while forgetting to update {{< code "E" >}}, ending up with an inconsistent file.

{{< manual-foot prev="Basics:Installation" next="Basics:Running Octopus" >}}
---
title: "Discretization"
#series: "Manual"
weight: 6
description: " "
manuals: "Basics"
---


Besides all the approximations we have to do, no computer can solve an infinite continuous problem. We have to discretize our equations somehow. {{< octopus >}} uses a grid in real space to solve the Kohn-Sham equations. That is, functions are represented by their value over a set of points in real space. Normally the grid is equally spaced, but also non-uniform grids can be used. The shape of the simulation region may also be tuned to suit the geometric configuration of the system.

## Grid

In {{< octopus >}} functions are represented by a set of values that correspond to the value of the function over a set of points in real space. By default these points are distributed in a uniform grid, which means that the distance between points is a constant for each direction. It is possible to have grids that are not uniform, but as this is a bit more complex we will discuss it later. 

In this scheme, the separation between points, or spacing, is a critical value. When it becomes large the representation of functions gets worse and when it becomes small the number of points increases, increasing memory use and calculation time. This value is equivalent to the {{< name "energy cutoff" >}} used by plane-wave representations.

In {{< octopus >}} you can choose the Spacing of your simulation by the {{< variable "Spacing" >}} variable. If you set this as a single value it will give the same spacing for each directions (for example {{< code "Spacing=0.3" >}}). If you want to have different spacings for each coordinate you can specify {{< variable "Spacing" >}} as a block of three real values.

If you are working with the default pseudopotential species of {{< octopus >}} they come with a recommended {{< variable "Spacing" >}} value and you don't need to give it in the input file. Normally this default values are around 0.4 {{< bohr >}} (~0.2 {{< angstrom >}}), but you may need smaller spacing in some cases. Do not rely on these default values for production runs.

### Double grid (experimental)

The double-grid technique is a method to increase the precision of the representation of the pseudopotentials in the grid that has been recently integrated in {{< octopus >}} (not available in {{< octopus >}} 2.1 and previous versions). To activate this technique, set the variable {{< variable "DoubleGrid" >}} to {{< value "yes" >}}. The use of a double grid increases the cost of the transfer of the potential to the grid, but as in most cases this is done only a few times per run, the overhead to the total computation time is negligible. The only exception is when the atoms are displaced while doing a time-dependent simulation, where the double grid can severely increase the computation time.

## Box

We also have to select a finite domain of the real space to run our simulation, which is known as the simulation box. {{< octopus >}} can use several kinds of shapes of box. This is controlled by the variable {{< variable "BoxShape" >}}. Besides standard shapes {{< octopus >}} can take shapes given by a user-defined function or even by an image.

The way to give the size of the simulation box changes for each shape, but for most of them it is given by the variable {{< variable "Radius" >}}.

### Zero boundary conditions

By default {{< octopus >}} assumes zero boundary conditions, that is, wavefunctions and density are zero over the boundary of the domain. This is the natural boundary condition when working with finite systems.

In this case choosing an adequate box size is very important: if the box is too small the wavefunctions will be forced to go to zero unnaturally, but if the box is too large, a larger number of points is needed, increasing calculation time and memory requirements.

{{< manual-foot prev="Basics:Hamiltonian" next="Basics:Output" >}}
---
title: "Phonons"
#series: "Manual"
weight: 7
hidden: true
---

to be continued...

---------------------------------------------
---
title: "Optimal Control"
#series: "Manual"
weight: 5
---


{{< tutorial "Unsorted/Basic_QOCT" "Tutorial on basic Optimal Control Theory" >}}

{{< manual-foot prev="Calculations:Linear Response" next="Calculations:Geometry Optimization" >}}
---
title: "Calculations"
Weight: 20
---

{{% children depth=3 %}}---
title: "Ground State"
#series: "Manual"
weight: 1
description: "Converge the ground state of a system"
---


The ground-state electronic density in a Kohn-Sham (KS)-based electronic-structure code such as {{< octopus >}} is obtained after a self-consistent process that attempts to solve the KS equations. 

## Kohn-Sham Ground State

In essence, the problem is the following: at a given iteration step, one departs from an approximate solution – some KS eigenfunctions $\psi^{inp}_j$, eigenvalues $\epsilon^{inp}_j$ and density $\rho^{inp}$, which determines a KS “input” Hamiltonian. By diagonalizing this Hamiltonian, one obtains the corresponding “output” eigenfunctions, eigenvalues, and density. This density (or, alternatively, the corresponding Kohn-Sham potential) is then used to build a new input Hamiltonian, that will be diagonalized in the next iteration step. This cycle is considered to be closed, and the solution achieved, when the input and output are similar enough that some convergence criterion is fulfilled. In our case, we have allowed for four different criteria, to be defined below. The self-consistent procedure will stop either when the first of the convergence criterions is fulfilled, or when a maximum number of iterations has been performed.

### Mixing

The output density (or potential) of a given iteration is not used directly to construct the Kohn-Sham potential for the following iteration. Instead, it is "mixed" with some densities (or potentials) of previous iteration steps. The manner in which this mixing is produced is determined by the variables {{< variable "MixingScheme" >}}, {{< variable "Mixing" >}}, {{< variable "MixField" >}} and {{< variable "MixNumberSteps" >}}.

### Convergence
After each iteration {{< octopus >}} checks whether some convergence criterion is met. One criterion is that the error in the electron density be smaller than some threshold. Of course, the true electron density is not known, so this "error" is really the change in the density since the last iteration:
$\epsilon = \int {\rm d}^3r |\rho^{out}(\mathbf r) -\rho^{inp}(\mathbf r)|.$ We call this criterion {{< variable "ConvAbsDens" >}}.

However, since the density is proportional to the number of electrons $N$, this absolute criterion is not very transferable between different system sizes. Therefore the default criterion to use is {{< variable "ConvRelDens" >}}, in which the relative density, $ {\rho \over N}$, is used rather than the absolute density. $\epsilon = {1\over N} \int {\rm d}^3r |\rho^{out}(\mathbf r) -\rho^{inp}(\mathbf r)|.$ By default we use a value of <tt>1e-5</tt>.

The other convergence variables {{< variable "ConvAbsDens" >}}, {{< variable "ConvAbsEv" >}}, and {{< variable "ConvRelEv" >}} are set to <tt>0</tt>, indicating that they will not be used as criteria for convergence.

These other available criteria use errors defined as follows:

{{< variable "ConvAbsEv" >}}: The change in each eigenvalue is found and the sum of these changes must be smaller than the threshold.
$ \epsilon = \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert. $

{{< variable "ConvRelEv" >}}: The eigenvalues are scaled by total eigenvalue sum, otherwise as above.
<!--$\epsilon = {1 \over N} \sum_{j=1}^{N_{occ}} \vert \epsilon_j^{out}-\epsilon_j^{inp}\vert.$-->
$\epsilon = \frac{ \left| \sum_{j=1}^{N_{occ}} ( \epsilon_j^{out} - \epsilon_j^{inp} ) \right|}
    {\left| \sum_{j=1}^{N_{occ}} \epsilon_j^{out} \right|}$


To use them, set the relevant variable equal to a number indicating your desired error. Set the other three convergence variables to zero.

### Eigensolver

In each iteration of the self-consistency problem described above, {{< octopus >}} must diagonalize the Kohn-Sham input Hamiltonian to obtain the output eigenfunctions and eigenvalues. This diagonalization is done by the eigensolver, also an iterative procedure. There are several options for the iterative scheme used to diagonalize the Hamiltonian, which are specified in the documentation for the variable 
{{< variable "Eigensolver" >}}. You may specify the threshhold for considering this iterative diagonalization finished with the variable {{< variable "EigensolverTolerance" >}}. The variable {{< variable "EigensolverMaxIter" >}} sets a maximum number of steps in the diagonalization, so that if it is reached, the diagonalization is considered finished even if some of the eigenvectors are not fully converged.

During each self-consistent field cycle iteration {{< octopus >}} reports the eigenvalues it has obtained by diagonalizing the Hamiltonian, and how many of those eigenvectors are fully converged.

```text
                                                                              
*********************** SCF CYCLE ITER -    3 ************************
 etot =  4.52092631E+00 abs_ev   =  7.91E+02 rel_ev   =  4.28E+01
                        abs_dens =  5.16E-01 rel_dens =  1.43E-02
Matrix vector products:   3669
Converged eigenvectors:      6
Eigenvalues [H]
 -st  Spin   Eigenvalue     Occupation       Error
   1   --    -1.288198       2.000000      (7.2E-07)
   2   --    -0.830676       2.000000      (1.0E-06)
   3   --    -0.826885       2.000000      (8.8E-07)
   4   --    -0.808297       2.000000      (6.2E-07)
...
```

It is not too important whether the eigenvectors are not converged in the SCF steps, only whether they are converged at the end.

### LCAO

Since the solution of the ground-state problem is done iteratively, we need an initial guess; a set of initial Kohn-Sham orbitals. If we are doing the calculation with pseudopotentials (as opposed to model potentials defined by the user), we can use the pseudo-orbitals that are used to generate the pseudopotential. By default, the program will fill the initial guess states with pseudo-orbitals. Whether or not this is done is determined by the variable {{< variable "LCAOStart" >}}. The guess density is the sum of the atomic densities.

Note, however, that those pseudo-orbitals are not passed directly to the iterative cycle. Instead, the code performs an initial diagonalization with the Hamiltonian generated by the guess density. Therefore, the SCF cycle is started with the linear combination of those atomic orbitals, with the coefficients that result of that diagonalization (LCAO stands for linear combination of atomic orbitals). In other words, the first step of the SCF cycle is performed inside the LCAO subspace, whereas the following steps are performed in the full space.

This diagonalization will typically be done with a number of pseudo-orbitals that is larger than the number that will be used later in the KS SCF cycle. Once we diagonalize that LCAO matrix, we take the lowest lying eigenstates to proceed with the calculation. There is some default number of pseudo-orbitals that will be used, but one can change it making use of variable {{< variable "LCAODimension" >}}.

If you set {{< variable "SCFinLCAO" >}}, the LCAO calculation will be performed self-consistently. Or, in other words, the whole SCF cycle will be done inside the LCAO subspace.


## Unoccupied states

This is {{< variable "CalculationMode" >}} {{< code "= unocc" >}}. The purpose of this run mode is to calculate higher lying Kohn-Sham orbitals. For that purpose, it reads the restart information from a converged previous ground-state calculation, and builds the corresponding Hamiltonian. Then, it calculates the unoccupied eigenvalues and eigenfunctions. The number of unoccupied orbitals calculated is given by the {{< variable "ExtraStates" >}}.

{{< manual-foot prev="Calculations:Troubleshooting" next="Calculations:Time-Dependent" >}}
---
title: "Linear Response"
#series: "Manual"
weight: 4
description: "Polarizabilities"
---


{{< octopus >}} can calculate dynamic polarizabilities and first-order hyperpolarizabilites in a linear-response scheme using the Sternheimer equation. It is also possible to calculate optical spectra with this technique, but it is slower than time-evolution.

##### Ground state

The first thing we will need for linear response is a {{< manual "Calculations/Ground State" "Ground State" >}} calculation. Unlike the Casida approach, when using the Sterheimer equation you needn't do a unoccupied-states calculation. To improve the convergence of the linear-response calculation, it is better to use tightly converged wavefunctions. For example, you can add these parameters to your gs calculation:

```text

EigenSolverFinalTolerance = 1e-10
ConvRelDens = 1e-9
```

##### Input

The {{< variable "CalculationMode" >}} for polarizability calculations is {{< code "em_resp" >}}. The main parameter you have to specify is the frequency of the perturbation, given by the {{< variable "EMFreqs" >}} block. You can also add an imaginary part to the frequency by setting the variable {{< variable "EMEta" >}}. Adding a small imaginary part is required if you want to get the imaginary part of the polarizability or to calculate polarizabilities near resonance; a reasonable value is {{< code "0.1 eV" >}}.

To get the hyperpolarizabilties, you also have to specify the variable {{< variable "EMHyperpol" >}} with the three coefficients with respect to the base frequency; the three values must sum to zero.

##### Output

After running, for each frequency in the input file, {{< octopus >}} will generate a subdirectory under {{< file "em_resp/" >}}. In each subdirectory there is a file called {{< file "alpha" >}} that contains the real part of the polarizability tensor $\alpha_{ij}$ and the average polarizability

$$
\\bar{\\alpha}=\\frac13\\sum\_{i=1}^3\\alpha\_{ii}\\,
$$

The imaginary part $\eta$ is written to file {{< file "eta" >}}. If $\eta > 0$, there is also a file called {{< file "cross_section_tensor" >}} that contains the photo-absorption cross section tensor for that frequency, related to the imaginary part of the polarizability ($\sigma = \frac{4 \pi \omega}{c} \mathrm{Im} \alpha $).

The hyperpolarizability will be in a file called {{< file "beta" >}} at the base frequency, containing all the 27 components and some reduced quantities:

$$
\\beta\_{||\\,i} = \\frac15 \\sum\_{j=1}^3(\\beta\_{ijj}+\\beta\_{jij}+\\beta\_{jji})\\ .
$$

Optionally, Born charges can also be calculated.

#### Finite differences

In this mode only static polarizability can be obtained. The calculation is done by taking the numerical derivative of the energy with respect to an external static and uniform electric field. To use this, run with {{< variable "ResponseMethod" >}} {{< code "=finite_differences" >}}. {{< octopus >}} will run several ground-state energy calculations and then calculate the polarizability using a finite-differences formula for the derivative. The results will be in the {{< file "em_resp_fd" >}} directory. Hyperpolarizability and Born charges can also be calculated.

{{< manual-foot prev="Calculations:Casida" next="Calculations:Optimal Control" >}}
---
title: "Casida"
#series: "Manual"
weight: 3
description: "Alternative to linear response"
---


Mark Casida's formulation of Linear-Response TDDFT allows calculations of the excitation energies of a finite system. For small molecules, this is normally the fastest way to calculate them.

To perform a Casida calculation you first need a {{< manual "Calculations/Ground State" "ground-state" >}} calculation and a calculation of unoccupied states; then you run the code with {{< variable "CalculationMode" >}}{{< code "=casida" >}}.

Here are the input files for a linear-response calculation of the nitrogen dimer, found at {{< file "share/testsuite/linear_response/01-casida.*" >}}.

#### Ground state calculation:

The first step is to converge the ground-state:

{{< expand "Ground State input file" >}}
{{< code-block >}}
#include_input testsuite/linear_response/01-casida.01-gs.inp
{{< /code-block >}}
{{< /expand >}}

#### Unoccupied states

Once the ground state calculation is converged, it is necessary to add and converge more unoccupied states.

{{< expand "Unoccupied States input file" >}}
{{< code-block>}}
#include_input testsuite/linear_response/01-casida.0?-unocc.inp
{{< /code-block >}}
{{< /expand >}}

#### Casida calculation for excited states

{{< expand "Casida input file" >}}
{{< code-block>}}
#include_input testsuite/linear_response/01-casida.0?-casida.inp
{{< /code-block >}}
{{< /expand >}}

See also the {{< tutorial "Response/Optical_spectra_from_Casida" "tutorial on Casida calculations">}}.

{{< manual-foot prev="Calculations:Time-Dependent" next="Calculations:Linear_Response" >}}
---
title: "Geometry Optimization"
#series: "Manual"
weight: 6
---


To perform a geometry optimization with {{< octopus >}} one should set the {{< variable "CalculationMode" >}} to {{< code "go" >}}. The method to be used should be set using the variable {{< variable "GOMethod" >}}. Note that most of the methods use the minimization routines from [GSL](https://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html).

The stopping criteria can be set with the variables {{< variable "GOTolerance" >}} and {{< variable "GOMinimumMove" >}}. Then minimization will be stopped when all forces on ions are smaller than {{< variable "GOTolerance" >}} or when all species coordinates change less than {{< variable "GOMinimumMove" >}} during one minimization step. If none of the previous criteria is matched after a number of minimization steps equal to {{< variable "GOMaxIter" >}}, then the minimization will be stopped with an error message.

After each minimization step taken by the [GSL](https://www.gnu.org/software/gsl/) algorithms {{< octopus >}} will write the current geometry to a file named {{< file "go.XXXX.xyz" >}} (XXXX is the iteration number) in the directory {{< file "geom" >}}. As an extra information, the title of the xyz file (second line) contains the total energy.

Currently, the default method used for geometry optimizations is {{< variable "GOMethod" >}} = {{< code "FIRE" >}}. Otherwise, many times, if the stopping criteria are too small, the GSL routines will return an error message stating that they were unable to improve the solution. 

At the end of the run, if the minimization was successful, the minimized geometry will be written to a file named {{< file "min.xyz" >}} in the working directory.

See also the tutorial on {{< tutorial "Unsorted/Geometry_optimization" "geometry optimization">}}

{{< manual-foot prev="Calculations:Optimal Control" next="Visualization" >}}
---
title: "Time-Dependent"
#series: "Manual"
weight: 2
description: "Time evolution"
---


### Time evolution

When {{< variable "CalculationMode" >}} = {{< code "td" >}}, the code performs the time-propagation of the electronic orbitals and &ndash; if required &ndash; the ionic positions. The latter task does not pose major algorithmical problems (the usual Verlet algorithms deal with that task); however the best way to propagate a Schr&ouml;dinger-like equation is still unclear. Due to this fact, we provide with a rather excessive selection of possibilities for that purpose. Before describing the set of variables necessary to specify the way in which the time evolution is to be performed, it is worth making a brief introduction to the problem.

We are concerned with a set of Schrödinger-like equations for the electronic orbitals $\psi_j(t)$:

$$
i{\\partial \\over \\partial t}\\psi\_j(t) = H(t)\\psi\_j(t) \\,\\!
$$

$$
\\psi\_j(t=0) = \\psi\_j^{(0)}\\,\\!
$$

Because this equation is linear (the time derivative and the Hamiltonian are both linear operators), one may formally define a linear “evolution” operator, $U(T, t)$, which transforms the initial vector into the solution at time $T$:

$$
\\psi_j(T) = U(T, 0)\\psi\_j^{(0)}
$$

Moreover, there is the formally exact expression for the evolution operator

$$
\\psi\_j(T) = \\mathcal{T} \\! \\! \\exp \\left\\{ -i \\int\_0^{T} d \\tau H(\\tau) \\right\\} \\psi\_j^{(0)}
$$

where $\mathcal{T} \\! \\! \exp$ is the [time-ordered exponential](https://en.wikipedia.org/wiki/Ordered_exponential), which is a short-hand for:

$$
\\psi\_j(T) = \\left\\{\\sum\_{n=0}^{\\infty} \\frac{\\left(-i\\right)^n}{n!} \\int\_0^t d\\tau\_1 \\cdots \\int\_0^t d\\tau\_n \\mathcal{T} H(\\tau\_1) \\cdots H(\\tau\_n)\\right\\}\\psi\_j^{(0)}\\,\\!
$$

If the Hamiltonian commutes with itself at different times, we can drop the time-ordering product, and leave a simple exponential. If the Hamiltonian is time-independent, which makes it trivially self commuting, the solution is simply written as:

$$
\\psi\_j(T) = \\exp\\left\\{ -iTH\\right\\} \\psi\_j^{(0)}\\,\\!.
$$

Unfortunately, this is not the case for TDDFT when the system is exposed to external time-dependent perturbations like electric and magnetic fields or pulsed lasers. But even without an external time-dependency, there remains the intrinsic time-dependency of the Kohn-Sham Hamiltonian, which is built &ldquo;self-consistently&rdquo; from the varying electronic density.

The first step to tackle this problem is to split the propagation of the long interval $[0, T]\,\\!$ into $N\\,\\!$ smaller steps by utilizing the group property

$$
U(T, t) = U(T, t')U(t', t)\\,\\!
$$

of the time-evolution operator. This yields the following time discretization:

$$
U(T, 0) = \\prod\_{i=0}^{N-1}U(t\_i+\\Delta t, t\_i)\\,\\!,
$$

where $t_0=0\\,\\!$, $t_N=T\\,\\!$, $\Delta t = T/N\,\\!$. So at each time step we are dealing with the problem of performing the short-time propagation:

$$
\\psi\_j(t+\\Delta t) = U(t+\\Delta t, t)\\psi\_j(t) = \\mathcal{T}\\!\\!\\exp\\left\\{ -i\\int\_{t}^{t+\\Delta t}d\\tau H(\\tau)\\right\\} \\psi\_j(t)\\,\\!.
$$

In this way, one can monitor the evolution in the interior of $[0, T]\\,\\!$. In fact, the possibility of monitoring the evolution is generally a requirement.

This requirement imposes a natural restriction on the maximum size of $\Delta t\\,\\!$: if $\omega_{\mathrm{max}}\\,\\!$ is the maximum frequency that we want to discern, $\Delta t\\,\\!$ should be no larger than $\approx 1/\omega_{\mathrm{max}}\\,\\!$. Below this $\Delta t_{\mathrm{max}}\\,\\!$, we are free to choose $\Delta t\\,\\!$ considering performance reasons: technically, the reason for the discretization is two-fold: the time-dependence of $H\\,\\!$ is alleviated, and the norm of the exponential argument is reduced (the norm increases linearly with $\Delta t\\,\\!$).

Since we cannot drop the time-ordering product, the desired algorithm cannot be reduced, in principle, to the calculation of the action of the exponential of an operator over the initial vector. Some algorithms tailored to approximate the evolution operator, in fact, do not even need to peform such operator exponentials. Most of them, however, do rely on the calculation of one or more exponentials, such as the ones used by {{< octopus >}}. This is why in principle we need to specify two different issues: the &ldquo;evolution method&rdquo;, and the &ldquo;exponential method&rdquo;. In other words: we need an algorithm to approximate the evolution operator $U(t+\Delta{t},t)\\,\\!$ &ndash; which will be specified by variable {{< variable "TDPropagator" >}} &ndash; and, if this algorithm requires it, we will also need an algorithm to approximate the exponential of a matrix operator $\\exp\\left\\{ A\right\\}\\,\\!$ &ndash; which will be specified by variable {{< variable "TDExponentialMethod" >}}.

#### Propagators for the time-dependent Kohn-Sham equations

In the following, we describe the propagators available in {{< octopus >}} for the time-dependent Kohn-Sham equations. These propagators solve the problem of approximating the orbitals $\psi_j(t+\Delta t)\\,\\!$ from the knowledge of $\psi_j(\tau)\\,\\!$ and $H(\tau)\\,\\!$for $0\le\tau\le t\\,\\!$. Some methods require the knowledge of the Hamiltonian at some points $\tau\\,\\!$ in time between $t\\,\\!$ and $t+\Delta t\\,\\!$.

This quantity can be approximated by the following procedure:

- Approximate $H(\tau)\\,\\!$ through extrapolation of a polynomial fit to a certain number of previous time steps.
- Propagate $\psi_j(t)\\,\\!$ to get $\psi_j(t+\Delta t)\\,\\!$.
- Calculate $H(t+\Delta t)\\,\\!$ from the orbitals $\psi_j(t+\Delta t)\\,\\!$.
- Interpolate the required $H(\tau)\\,\\!$ from $H(t)\\,\\!$ and $H(t+\Delta t)\\,\\!$.
- Repeat the steps 2 to 4 until self-consistency is reached.

In {{< octopus >}}, however, the above scheme is dropped for performance reasons and only step 1 is implemented, via a second-order extrapolation, except for the first two steps where the extrapolation obviously cannot be trusted. Instead, we rely on a sufficiently small $\Delta t\\,\\!$.

##### Midpoint rules

The implicit midpoint rule or [Crank-Nicolson method](https://en.wikipedia.org/wiki/Crank%E2%80%93Nicolson_method) ({{< variable "TDPropagator" >}}={{< code "crank_nicolson" >}}) calculates the exponential of the Hamiltonian by a first-order Pad&eacute; approximation. The Hamiltonian is evaluated at $t + \Delta t/2\\,\\!$ and the integral is dropped:

$$
U\_{\\mathrm{CN}}(t + \\Delta t, t) = \\frac{1- i \\frac{\\Delta t}{2}H(t+\\Delta t/2)}{1+ \\frac{\\Delta t}{2}H(t+\\Delta t/2)}\\,\\!
$$

The calculation of the matrix fraction is transformed into the solution of the linear system

$$
L \\psi\_j(t+\\Delta t) = b\\,\\!
$$

with the known quantities

$$
L = 1+i\\frac{\\Delta t}{2}H(t+\\Delta t/2)\\,\\!
$$

and

$$
b = \\left\\{1-i\\frac{\\Delta t}{2}H(t+\\Delta t/2)\\right\\}\\psi\_j(t).\\,\\!
$$

The Crank-Nicolson scheme is unitary and preserves time-reversal symmetry.

A simpler but similar scheme is the exponential midpoint rule
({{< variable "TDPropagator" >}}={{< code "exp_mid" >}}), which is unitary and time-reversal-symmetry-preserving, but has the drawback that it requires fairly small time steps:

$$
U\_{\\mathrm{EM}}(t + \\Delta t, t) = \\exp\\left\\{-i\\Delta t H(t + \\Delta t/2)\\right\\}\\,\\!
$$

##### Magnus expansions

[Magnus expansion](https://en.wikipedia.org/wiki/Magnus_expansion) ({{< variable "TDPropagator" >}}={{< code "magnus" >}}) is the most sophisticated propagation implemented in {{< octopus >}}. According to Magnus, there exists an operator $\Omega(t+\Delta t, t)\\,\\!$ for which holds

$$
U(t+\\Delta t, t) = \\exp\\left\\{\\Omega(t+\\Delta t, t)\\right\\}\\,\\!
$$

given by the series

$$
\\Omega(t + \\Delta t, t) = \\sum\_{k=1}^\\infty \\Omega\_k(t + \\Delta t, t)\\,\\!
$$

that converges at least for some local environment of $t\\,\\!$.

The operators $\Omega_k(t+\Delta t, t)\\,\\!$ are generated with the recursive procedure

$$
\\Omega\_k(t+\\Delta t, t) = \\sum\_{l=0}^{k-1} \\frac{B\_l}{l!}\\int\_t^{t+\\Delta t} S^l\_k(\\tau)d\\tau,\\,\\!
$$

$$
S\_1^0(\\tau) = -iH(\\tau), \\quad S\_k^0=0\\ \\mathrm{for}\\ k\>1,\\,\\!
$$

$$
S\_k^j(\\tau) = \\sum\_{m=1}^{k-l}\\left\[\\Omega\_m(t+\\Delta t, t), S\_{k-m}^{l-1}(\\tau)\\right\]\\ \\mathrm{for}\\ 1\\le l\\le k-1,\\,\\!
$$

where $B_l\\,\\!$ are [Bernoulli numbers](https://en.wikipedia.org/wiki/Bernoulli_number).

In {{< octopus >}} we have implemented a fourth-order Magnus expansion, which means that the operator series is truncated to second order and the integrals are calculated by a second-order quadrature formula. The resulting operator is

$$
\\Omega\_{M(4)}(t + \\Delta t, t) = -i\\frac{\\Delta t}{2}\\left\[H(t\_1)+ H(t\_2)\\right\] - \\frac{\\sqrt{3}\\Delta t^2}{12}\\left\[H(t\_2), H(t\_1)\\right\]\\,\\!
$$

with the quadrature sampling points $t_{1,2} = t + \left[\frac{1}{2}\mp \frac{\sqrt{3}}{6}\right]\Delta t\\,\\!$.

With a modified Hamiltonian

$$
H\_{M(4)}(t, \\Delta t) = \\overline{H}(t, \\Delta t) + i \\left\[T + v\_{\\mathrm{ext}}^{\\mathrm{nonlocal}}(t), \\overline{\\Delta V}(t, \\Delta t)\\right\],\\,\\!
$$

$$
\\overline{H}(t, \\Delta t) = T + \\frac{1}{2}\\left\[V\_{\\mathrm{KS}}(t\_1)+V\_{\\mathrm{KS}}(t\_2)\\right\],\\,\\!
$$

$$
\\overline{\\Delta V}(t, \\Delta t) = \\frac{\\sqrt{3}}{12}\\Delta t\\left\[V\_{\\mathrm{KS}}(t\_1)-V\_{\\mathrm{KS}}(t\_1)\\right\]\\,\\!
$$

the propagator takes the form

$$
U\_{M(4)}(t+\\Delta t, t) = \\exp\\left\\{-i \\Delta t H\_{M(4)}(t, \\Delta t)\\right\\}.\\,\\!
$$

Note that only the nonlocal components of the Kohn-Sham Hamiltonian contribute to the commutator and, furthermore, that we make the assumption that the nonlocal potential $v_{\mathrm{ext}}^{\mathrm{nonlocal}}(t)\\,\\!$ does not vary significantly in the interval $[t, t+\Delta t]\\,\\!$. This nonlocal component stems from the ionic pseudopotentials. Consequently, its variation is caused by the ionic movement, which is negligible on the electronic time scale determining $\Delta t\\,\\!$.

##### Time-reversal-symmetry based propagation

Due to time-reversal symmetry, propagating backwards by $\Delta t/2\\,\\!$ starting from $\psi_j(t+\Delta t)\\,\\!$ should lead to the same result as propagating forwards by $\Delta t/2\\,\\!$ starting from $\psi_j(t)\\,\\!$. Using the simplest approximation to the evolution operator, we end up with the condition

$$
\\exp\\left\\{+i\\frac{\\Delta t}{2}H(t +\\Delta t)\\right\\} \\psi\_j(t+\\Delta t) = \\exp\\left\\{-i\\frac{\\Delta t}{2}H(t)\\right\\} \\psi\_j(t)\\,\\!
$$

Rearranging the terms gives an approximation for the propagator:

$$
U\_{\\mathrm{ETRS}}(t + \\Delta t, t) = \\exp\\left\\{-i \\frac{\\Delta t}{2} H(t + \\Delta t)\\right\\}
\\exp\\left\\{-i \\frac{\\Delta t}{2}H(t)\\right\\}\\,\\!
$$

This ''enforced time-reversal symmetry'' method is available in {{< octopus >}} by setting ({{< variable "TDPropagator" >}}={{< code "etrs" >}}).

It is worthwhile to give a sidenote on the approximation of $H(t+\Delta t)\\,\\!$: for the {{< code "etrs" >}} method the Hamiltonian at time $t+\Delta t\\,\\!$ is not extrapolated, but rather built from the density $n' = \sum_j|\psi'_j(t+\Delta t)|^2\\,\\!$, which, in turn, is calculated from the orbital estimate

$$
\\psi'\_j(t+\\Delta t) = \\exp\\left\\{-i\\Delta t H(t)\\right\\}\\psi\_j(t).\\,\\!
$$

There exists, however, also a variant {{< variable "TDPropagator" >}}={{< code "aetrs" >}} ''(approximated enforced time-reversal symmetry)'' that performs a second-order polynomial extrapolation of $H(t-k\Delta t)\\,\\!$, $k=0,1,2\\,\\!$, to get $H(t+\Delta t)\\,\\!$, which is about 40&nbsp;% faster than {{< code "etrs" >}}.

#### Approximations to the exponential of an operator

Most of the evolution methods described in the previous section require the calculation of the exponential of a matrix. Before going into the details of the {{< variable "TDExponentialMethod" >}} that {{< octopus >}} provides, one general difficulty deserves mentioning.

In principle, we would like to calculate $\\exp\\left\\{A\\right\\}v\\,\\!$ by first calculating $\\exp\\left\\{A\\right\\}\\,\\!$ and then applying this result to any vector $v\\,\\!$. Unfortunately, the size of our Hamiltonian matrix is of the order $\approx 10^5\\,\\!$ which forbids its full storage in matrix form. Apart from that, methods to calculate the exponential of a matrix explicitly are limited to a few thousand matrix elements. For this reason, we have to turn to iterative solutions that calculate $\\exp\\left\\{A\\right\\}v\\,\\!,$ for a particular vector $v$.

##### Polynomial expansions

The simplest possibility is to use the definition of the exponential of a matrix

$$
\\exp\\left\\{A\\right\\} = \\sum\_{k=0}^\\infty \\frac{1}{k!}A^k\\,\\!
$$

to get the approximation of order $k\\,\\!$ ({{< variable "TDExponentialMethod" >}}={{< code "taylor" >}}):

$$
\\mathrm{taylor}\_N\\left\\{A, v\\right\\} = \\sum\_{k=0}^N\\frac{1}{k!}A^kv\\,\\!
$$

$\\mathrm{taylor}\_N\\left\\{A, v\\right\\}\\,\\!$ amounts to the expansion of the exponential in the standard polynomial base $\\left\\{1, x, x^2, \\ldots\\right\\}\\,\\!$.

Experience shows that the choice $k=4\\,\\!$ gives particularly good results for our TDDFT implementation.

The standard polynomial basis is not the only possibility; {{< octopus >}} also implements the expansion of the exponential in the Chebyshev basis ({{< variable "TDExponentialMethod" >}}={{< code "chebyshev" >}}):

$$
\\mathrm{cheb}\_N\\left\\{A, v\\right\\} = \\sum\_{k=0}^Nc\_kT\_k(A)kv\\,\\!
$$

with $T_k\\,\\!$ being the [Chebyshev polynomial](https://en.wikipedia.org/wiki/Chebyshev_polynomials) of order $k\\,\\!$.

The truncation $N\\,\\!$ for both these expansions can be set by the input variable {{< variable "TDExpOrder" >}}.

##### Krylov-subspace projection

The [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace) of order $N\\,\\!$ for a given operator $A\\,\\!$ and vector $v\\,\\!$ is defined as

$$
\\mathcal{K}\_N\\left\\{A, v\\right\\} = \\mathrm{span}\\left\\{v, Av, A^2v, \\ldots, A^{N-1}v\\right\\}.\\,\\!
$$

The idea is to find the element of $\\mathcal{K}\_N \\left\\{A, v \\right\\}$, that optimally approximates $\\exp\\left\\{A\\right\\}\\,\\!$ in a least-squares sense. It can be proven that this is given by

$$
\\beta V\_N(V\_N^T\\exp\\left\\{A\\right\\}V\_N)e\_1\\,\\!
$$

with $V_N = [v_1, \ldots, v_N]\\,\\!$ being an orthonormal basis of $\\mathcal{K}\_N\\left\\{A, v\\right\\}\\,\\!$ calculated with the Lanczos procedure. To get rid of the exponential, the approximation

$$
V\_N^T\\exp\\left\\{A\\right\\}V\_N \\approx \\exp\\left\\{V\_N^TAV\_N\\right\\}\\,\\!
$$

is introduced. The object $H_N=V_N^TAV_N\\,\\!$ is also a result of the Lanczos procedure and of the small dimension $k\times k\\,\\!$ which allows for direct calculation of the remaining exponential. So, we have

$$
\\mathrm{lanczos}\_N\\left\\{A, v\\right\\} = V\_N\\exp\\left\\{H\_N\\right\\}e\_1,\\,\\!
$$

which can be selected by setting {{< variable "TDExponentialMethod" >}}={{< code "lanczos" >}}. When actually calculating the exponential, {{< octopus >}} recursively generates the quantities $H\_N\\,\\!$ and $V\_N\\,\\!$ until some convergence criterion &ndash; controlled by {{< variable "TDLanczosTol" >}} &ndash; is met or the order exceeds the maximum allowed order controlled by {{< variable "TDExpOrder" >}}. In the latter case, a warning is written. Care should be taken in choosing the convergence parameter small enough to avoid faulty evolutions.

###  Performing a time-propagation  

To actually perform a time-propagation, a necessary prerequisite is to have the ground-state density (see {{< manual "Calculations:Ground_State" "Ground State" >}}) of the system under consideration.

Then, the {{< variable "CalculationMode" >}} has to be set to {{< code "td" >}}.

Besides the already mentioned input variables concerning the propagation method, two more important parameters have to be set: {{< variable "TDTimeStep" >}} and {{< variable "TDMaxSteps" >}}. The first one sets $\\Delta t\\,\\!$, the second sets the number of iterations. It is convenient to define these two variables like this:


 T  = 0.1     - Length of propagation.
 dt = 0.002   - Length of time step.
 
 {{< variable "TDMaxSteps" >}} = T/dt
 {{< variable "TDTimeStep" >}} = dt


In the absence of a perturbation, the system's total energy should be a constant. This check to determine a reasonable time step has to be done before any other calculation.

The relevant parts of {{< octopus >}}' output might look like

```text
  Iter           Time        Energy     Elapsed Time
      1       0.002000   -541.542073         4.800
      2       0.004000   -541.542073         4.560
      3       0.006000   -541.542073         3.880
      .
      .
      .
     48       0.096000   -541.542073         7.420
     49       0.098000   -541.542073         7.620
     50       0.100000   -541.542073         7.490
```

The energy is printed in the third column and should remain constant.
In general, larger time steps are desirable to shorten computational propagation time but keep in mind that the size of the time step limits to the maximum frequency you will be able to observe, and, consequently, also any external perturbation to which you might wish to expose the system.

### External fields 



##### Delta kick: Calculating an absorption spectrum

To obtain the linear optical absorption spectrum of the system, we follow the scheme proposed by Yabana and Bertsch, and
excite all frequencies of the system by giving some small
momentum (${\\mathcal K}$) to the electrons. This is achieved by transforming the ground-state wavefunctions according to:

$$
  \\varphi\_i(r, \\delta t) = e^{i {\\mathcal K} z} \\varphi\_i(r, 0)\\,,
$$

and then propagating these wavefunctions for some (finite) time.
The spectrum can then be obtained from the expression for the dipole
strength function $S(\omega)$:

$$
  S(\\omega) = \\frac{2 \\omega}{\\pi} {\\mathcal Im} \\alpha(\\omega)\\,,
$$

where the dynamical polarizability, $\alpha(\omega)$, is
essentially the Fourier transform of the dipole moment of the system
$d(t)$:

$$
  \\alpha(\\omega) = \\frac{1}{\\mathcal{K}}\\int\\!\\! dt \\; 
  e^{i\\omega t} \\left\[d(t)-d(0)\\right\]\\,.
$$

With this definition, the Thomas-Reiche-Kuhn ''f''-sum rule 
for the number of electrons, $N$, is given by the
integral:

$$
 N = \\int\\!\\! d\\omega\\;S(\\omega)\\,.
$$

This sum rule can be used to check the quality of the calculations. Another check
is energy conservation, which TDDFT respects when no external field is applied.

To obtain a spectrum with such a recipe in {{< octopus >}} one has to follow the steps:

- Choose a system of linearly independent (can be non-orthogonal) axes (defined by the variable {{< variable "TDPolarization" >}}). The defaults are the 3 Cartesian axes.
- Add an electric kick using the variable {{< variable "TDDeltaStrength" >}}. Its value should be small so that we remain in the linear regime, but large enough to avoid numerical errors.
- Run a time-dependent simulation for the 3 independent directions. This involves running {{< octopus >}} 3 times, changing the value of {{< variable "TDPolarizationDirection" >}} each time. After each run, move the file {{< code "td.general/multipoles" >}} to {{< code "td.general/multipoles.X" >}}, where {{< code X >}} is 1, 2, or 3.
- Run the utility {{< manual "External_utilities:oct-propagation_spectrum" "oct-propagation_spectrum" >}}.

##### Lasers

To calculate non-linear optical properties, we follow
the evolution of the system under the influence of a laser field that is treated 
in the dipole approximation (although this constraint can be removed). 
The harmonic emission spectrum can then be calculated from
the acceleration of the dipole moment:

$$
  H(\\omega) \\propto \\left| \\int\\!\\! dt \\; e^{i\\omega t} \\frac{d^2}{dt^2} d(t) \\right|^2\\,.
$$

During the propagation, charge density is absorbed at the boundaries
of the simulation region, either by an imaginary absorbing potential
or a mask function. In the first case, we add to the Kohn-Sham potential a term

$$
  V\_{\\rm eff}(r, t) = V\_{\\rm KS}(r,t) - i V\_{\\rm abs}(r)
  \\,,
$$

where $V_{\rm abs}$ is zero in the inner region of the simulation box, and
rises smoothly till the edges. By adjusting both the height and 
the shape of the potential, we can select which momenta are absorbed
and prevent the unwanted reflections at the boundary.
When using a mask, the wavefunction is multiplied in each time-step
by a function which is 1 in the inner simulation region and gradually 
goes to 0 at the borders:

$$
  \\varphi(r, t) \\rightarrow M(r) \\varphi(r, t)
  \\,.
$$

The absorbed charge can be interpreted as an ionization probability
and can be used to estimate the photo-electron spectra.
The box size has to be big enough so that the physical system is not
perturbed by the absorbing boundaries. Note that the wavefunctions are no 
longer normalized as the system slowly gets charged.


### Symmetries

The dynamic polarizability (trivially related to optical absorption) is, in its most general form, a 3x3 tensor. The reason is that we can shine light on the system polarized in any of the three Cartesian axes, and for each of these three cases measure how the dipole of the molecule oscillates along the three Cartesian axes. This usually means that to obtain the full dynamic polarizability of the molecule we usually need to apply 3 different perturbations along $x, y, z\,$, by setting {{< variable "TDPolarizationDirection" >}} to 1, 2, or 3.

However, if the molecule has some symmetries, it is in general possible to reduce the total number of calculations from 3 to 2, or even 1. This is explained in detail [here](https://nano-bio.ehu.es/node/576). To use this formalism in {{< octopus >}}, you can use the variables {{< variable "TDPolarization" >}}, {{< variable "TDPolarizationDirection" >}}, {{< variable "TDPolarizationEquivAxes" >}}, and {{< variable "TDPolarizationWprime" >}}. The block {{< variable "TDPolarization" >}} defines a basis (not necessarily orthogonal) chosen to maximize the number of equivalent axes, and {{< variable "TDPolarizationDirection" >}}, and {{< variable "TDPolarizationWprime" >}} are vectors specified in that basis.
Let us give a couple of examples.

The methane molecule has $T_d$ symmetry, which means that the response is identical for all directions. This means that we only need one propagation to obtain the whole tensor. This propagation can be performed in any direction we wish. So we could use the input

```text
 %{{< variable "TDPolarization" >}}
  1 | 0 | 0
  0 | 1 | 0
  0 | 0 | 1
 %
 {{< variable "TDPolarizationDirection" >}} = 1
 {{< variable "TDPolarizationEquivAxes" >}} = 3
 %{{< variable "TDPolarizationWprime" >}}
  0 | 0 | 1
 %
```

Note that we could have omitted the blocks {{< variable "TDPolarization" >}} and {{< variable "TDPolarizationWprime" >}} in the previous input file, as these are their default values. 

Now let us look at a linear molecule. In this case, you might think that we need two calculations to obtain the whole tensor, one for the direction along the axis of the molecule, and another for the axis perpendicular to the molecule. The fact is that we need only one, in a specially chosen direction, so that our field has components both along the axis of the molecule and perpendicular to it. Let us assume that the axis of the molecule is oriented along the $x$-axis. Then we can use

```text
 %{{< variable "TDPolarization" >}}
  1/sqrt(2) | -1/sqrt(2) | 0
  1/sqrt(2) |  1/sqrt(2) | 0
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
 {{< variable "TDPolarizationDirection" >}} = 1
 {{< variable "TDPolarizationEquivAxes" >}} = 3
 %{{< variable "TDPolarizationWprime" >}}
  0 | 0 | 1
 %
```

You should try to convince yourself that the three axes are indeed equivalent and linearly independent.

Finally, let us look at a general planar molecule (in the ''xy'' plane). In principle we need only two calculations (reduced to one if more symmetries are present, as for benzene). In this case we chose one of the polarization axes on the plane, and the other two rotated 45 degrees out of plane:

```text
 %{{< variable "TDPolarizationDirection" >}}
  1/sqrt(2) | 0 | 1/sqrt(2)
  1/sqrt(2) | 0 |-1/sqrt(2)
  0         | 1 | 0
 %
 
 {{< variable "TDPolarizationEquivAxes" >}} = 2
```

In this case, we need two runs, one for {{< variable "TDPolarizationDirection" >}} equal to 1, and another for it equal to 3. Note that if there are fewer than 3 equivalent axes, {{< variable "TDPolarizationWprime" >}} is irrelevant.


{{< manual-foot prev="Calculations:Ground_State" next="Calculations:Casida" >}}
---
title: "Benzene"
#series: "Manual"
---


Well, the sodium atom is a bit too trivial. Let's try something harder: benzene. you will just need the geometry for benzene to be able to play. Here it is (in Å):

```text
      C  0.000  1.396  0.000
      C  1.209  0.698  0.000
      C  1.209 -0.698  0.000
      C  0.000 -1.396  0.000
      C -1.209 -0.698  0.000
      C -1.209  0.698  0.000
      H  0.000  2.479  0.000
      H  2.147  1.240  0.000
      H  2.147 -1.240  0.000
      H  0.000 -2.479  0.000
      H -2.147 -1.240  0.000
      H -2.147  1.240  0.000
```

Follow now the steps of the previous example. Carbon and Hydrogen have a much harder pseudo-potential than Sodium, so you will probably have to use a tighter mesh. It also takes much more time... 

{{< manual-foot series="" prev="Examples:Hello world" next="Updating to a new version" >}}
---------------------------------------------
---
title: Examples
weight: 30

---

{{% children depth=3 %}}

---
title: "Hello world"
#series: "Manual"
---


As a first example, we will take a sodium atom. With your favourite text editor, create the file {{< file "inp" >}}.

```text
 {{< variable "CalculationMode" >}} = gs
 %{{< variable "Coordinates" >}}
     'Na' | 0.0 | 0.0 | 0.0 
 %
```

This input file should be essentially self-explanatory. 

Note that when a species is not specified in the {{< variable "Species" >}} block, octopus reads the information of pseudopotentials from the {{< file "defaults" >}} file (located under {{< inst_file "share/octopus/PP/" >}}. This file also contains default values for {{< variable "Radius" >}} and {{< variable "Spacing" >}}.

Then run octopus – for example, do 

{{< command-line "octopus > out" >}}

so that the output is stored in {{< file "out" >}} file. If everything goes OK, {{< file "out" >}} should look like:
```text

    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                ___
                             .-'   `'.
                            /         \
                            |         ;
                            |         |           ___.--,
                   _.._     |0) ~ (0) |    _.---'`__.-( (_.
            __.--'`_.. '.__.\    '--. \_.-' ,.--'`     `""`
           ( ,.--'`   ',__ /./;   ;, '.__.'`    __
           _`) )  .---.__.' / |   |\   \__..--""  """--.,_
          `---' .'.''-._.-'`_./  /\ '.  \ _.-~~~````~~~-._`-.__.'
                | |  .' _.-' |  |  \  \  '.               `~---`
                 \ \/ .'     \  \   '. '-._)
                  \/ /        \  \    `=.__`~-.
             jgs  / /\         `) )    / / `"".`\
            , _.-'.'\ \        / /    ( (     / /
             `--~`   ) )    .-'.'      '.'.  | (
                    (/`    ( (`          ) )  '-;
                     `      '-;         (-'

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA

    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                           Running octopus

Version                : 5.0.1
Revision               : 15042
Build time             : Tue Jan 12 11:04:57 EST 2016
Configuration options  : max-dim=3 mpi sse2
Optional libraries     : arpack berkeleygw etsf_io gdlib metis mpi2 netcdf newuoa parmetis parpack pfft scalapack sparskit
Architecture           : x86_64
C compiler             : /opt/local/bin/mpicc-mpich-mp (/usr/bin/clang)
C compiler flags       :  -pipe -O3 -arch x86_64
Fortran compiler       : /opt/local/bin/mpif90-mpich-mp (/opt/local/bin/gfortran-mp-5)
Fortran compiler flags : -pipe -O3

  The octopus is swimming in dhcp-18-189-27-45.dyn.MIT.EDU (Darwin)


            Calculation started on 2016/01/13 at 20:24:58


************************** Calculation Mode **************************
Input: [CalculationMode = gs]
**********************************************************************

Reading Coordinates from Coordinates block

****************************** Species *******************************
Reading pseudopotential from file:
      '/opt/local/share/octopus/PP/PSF/Na.psf'
      Calculating atomic pseudo-eigenfunctions for species Na....
Info: l =  0 component used as local potential.
Info: l =  0 is maximum angular momentum considered.
Number of orbitals: total =     16, bound =      4
**********************************************************************


***************************** Symmetries *****************************
Symmetry elements : (i) (Cinf) (sigma)
Symmetry group    : Kh
**********************************************************************

Input: [SpinComponents = unpolarized]
Input: [SmearingFunction = semiconducting]
Input: [SymmetrizeDensity = no]

******************************* States *******************************
Total electronic charge  =        1.000
Number of states         =        1
States block-size        =        1
**********************************************************************

Info: Using default spacing(1) [b] =  0.567
Info: Using default spacing(2) [b] =  0.567
Info: Using default spacing(3) [b] =  0.567
Input: [CurvMethod = curv_uniform]
Input: [DerivativesStencil = stencil_star]

************************** Parallelization ***************************
Octopus will run in *serial*
**********************************************************************

Info: Generating weights for finite-difference discretization of x-gradient
Info: Generating weights for finite-difference discretization of y-gradient
Info: Generating weights for finite-difference discretization of z-gradient
Info: Generating weights for finite-difference discretization of Laplacian

******************************** Grid ********************************
Simulation Box:
  Type = minimum
  Species =    Na     Radius =  13.228 b
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [b] = ( 0.567, 0.567, 0.567)    volume/point [b^3] =      0.18221
  - inner mesh =      52971
  - total mesh =      79699
  Grid Cutoff [H] =    15.354165    Grid Cutoff [Ry] =    30.708329
**********************************************************************

Info: states-block size = 0.6 MiB
Input: [StatesOrthogonalization = gram_schmidt]

****************************** Hartree *******************************
The chosen Poisson solver is 'interpolating scaling functions'
**********************************************************************


**************************** Theory Level ****************************
Input: [TheoryLevel = dft]

Exchange-correlation:
  Exchange
    Slater exchange (LDA)
    [1] PAM Dirac, Proceedings of the Cambridge Philosophical Society 26, 376 (1930)
    [2] F Bloch, Zeitschrift fuer Physik 57, 545 (1929)
  Correlation
    Perdew & Zunger (Modified) (LDA)
    [1] Perdew and Zunger, Phys. Rev. B 23, 5048 (1981)
    [2] Modified to improve the matching between the low- and high-rs parts

Input: [SICCorrection = sic_none]
**********************************************************************

Input: [FilterPotentials = filter_none]
Info: Pseudopotential for Na
  Radii for localized parts:
    local part     =  3.5 b
    non-local part =  0.0 b
    orbitals       = 19.9 b

Input: [RelativisticCorrection = non_relativistic]
Input: [AbsorbingBoundaries = not_absorbing]

****************** Approximate memory requirements *******************
Mesh
  global  :       1.5 MiB
  local   :       1.8 MiB
  total   :       3.4 MiB

States
  real    :       0.6 MiB (par_kpoints + par_states + par_domains)
  complex :       1.2 MiB (par_kpoints + par_states + par_domains)

**********************************************************************

Info: Generating external potential
      done.
Info: Octopus initialization completed.
Info: Starting calculation mode.
Info: Allocating ground state wave-functions
Info: Blocks of states
      Block       1 contains       1 states:       1 -       1
Info: Ground-state allocation done.

** Warning:
**   Could not find 'restart/gs' directory for restart.
**   No restart information will be read.


** Warning:
**   Unable to read wavefunctions.
**   Starting from scratch!

Input: [MixField = density] (what to mix during SCF cycles)
Input: [TypeOfMixing = broyden]

**************************** Eigensolver *****************************
Input: [Eigensolver = cg]
Input: [Preconditioner = pre_filter]
Input: [SubspaceDiagonalization = standard]
**********************************************************************

Input: [LCAOStart = lcao_full]
Input: [LCAOScaleFactor = 1.000]
Input: [LCAOMaximumOrbitalRadius = 20.00 b]
Info: Single-precision storage for     1 extra orbitals will be allocated.
Info: Unnormalized total charge =      0.999665
Info: Renormalized total charge =      1.000000
Info: Setting up Hamiltonian.
Info: Performing initial LCAO calculation with      2 orbitals.
Info: Getting Hamiltonian matrix elements.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation
   1   --    -0.103066       1.000000
Info: Ground-state restart information will be written to 'restart/gs'.
Info: SCF using real wavefunctions.
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    1 ************************
 etot = -1.84238245E-01 abs_ev   =  3.88E-04 rel_ev   =  3.75E-03
                        abs_dens =  6.54E-03 rel_dens =  6.54E-03
Matrix vector products:     27
Converged eigenvectors:      0
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103455       1.000000      (5.5E-05)

Elapsed time for SCF step     1:          0.08
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    2 ************************
 etot = -1.84238838E-01 abs_ev   =  1.02E-04 rel_ev   =  9.88E-04
                        abs_dens =  4.55E-03 rel_dens =  4.55E-03
Matrix vector products:     27
Converged eigenvectors:      0
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103353       1.000000      (1.3E-06)

Elapsed time for SCF step     2:          0.10
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    3 ************************
 etot = -1.84239600E-01 abs_ev   =  2.23E-04 rel_ev   =  2.16E-03
                        abs_dens =  4.35E-04 rel_dens =  4.35E-04
Matrix vector products:     19
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103130       1.000000      (9.3E-07)

Elapsed time for SCF step     3:          0.09
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    4 ************************
 etot = -1.84239592E-01 abs_ev   =  2.81E-07 rel_ev   =  2.73E-06
                        abs_dens =  6.62E-04 rel_dens =  6.62E-04
Matrix vector products:     12
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103130       1.000000      (7.2E-07)

Elapsed time for SCF step     4:          0.07
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER -    5 ************************
 etot = -1.84239609E-01 abs_ev   =  6.42E-06 rel_ev   =  6.23E-05
                        abs_dens =  9.12E-06 rel_dens =  9.12E-06
Matrix vector products:     15
Converged eigenvectors:      1
Eigenvalues [H]
 -st  Spin   Eigenvalue      Occupation     Error
   1   --    -0.103123       1.000000      (8.7E-07)

Elapsed time for SCF step     5:          0.07
**********************************************************************


             Info: Writing states. 2016/01/13 at 20:24:59


        Info: Finished writing states. 2016/01/13 at 20:24:59

Info: SCF converged in    5 iterations

Info: Finished writing information to 'restart/gs'.

             Calculation ended on 2016/01/13 at 20:24:59

                          Walltime:  01.  6s

Octopus emitted 2 warnings.
```

Take now a look at the working directory. Besides the initial file ({{< file "inp" >}}) and the {{< file "out" >}} file, three new directories appear. In {{< file "static/" >}}, you will find the file {{< file "info" >}}, with information about the static calculation (it should be hopefully self-explanatory, otherwise please complain to the authors...). In {{< file "restart/" >}}, you will find the {{< file "gs" >}} directory that contains restart information about the ground-state, which is used if, for example, you want to start a time-dependent calculation afterwards. Finally, the {{< file "exec" >}} directory has information about the run of octopus; inside the {{< file "parser.log" >}} contains all the input variables parsed by octopus.

### Exercises

* Study how the total energy and eigenvalue of the sodium atom improve with the mesh spacing.
* Calculate the static polarizability of the sodium atom ({{< variable "CalculationMode" >}} = em_resp). A {{< file "em_resp/freq_0.0000/alpha" >}} will be created containing the static polarizability tensor. 
* Calculate a few unoccupied states ({{< variable "CalculationMode" >}} = unocc). The eigenspectrum will be in the file {{< file "static/eigenvalues" >}}. Why don't we find a Rydberg series in the eigenspectrum?
* Repeat the previous calculation with different exchange and correlation functionals like PBE, LB94, and exact exchange (see {{< variable "XCFunctional" >}}).
* Perform a time-dependent evolution ({{< variable "CalculationMode" >}} = td), to calculate the optical spectrum of the Na atom. Use a {{< variable "TDDeltaStrength" >}} = 0.05, polarised in the ''x''-direction. The multipole moments of the density are output to the file {{< file "td.general/multipoles" >}}. You can process this file with the utility [[Manual:External_utilities:oct-propagation_spectrum | {{< command "oct-propagation_spectrum" >}}]] to obtain the optical spectrum. If you have computer time to waste, re-run the time-dependent simulation for some other xc choices. 

{{< manual-foot series="Examples" prev="Deprecated Utilities" next="Examples:Benzene" >}}
---
title: "oct-xyz-anim"
#series: "Manual"
---


### NAME 
oct-xyz-anim - Constructs an "animated" xyz file from a time-propagation file where the atoms were allowed to move

### SYNOPSIS 
{{< command "oct-xyz-anim" >}}

[oct-xyz-anim does not read the standard input: all standard input will
be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-xyz-anim accepts no command line arguments,
since there is not a standard way to do this with Fortran 90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities. It reads out the {{< file "td.general/coordinates" >}} file, and makes a movie in XYZ format, called 'movie.xyz'.

{{< manual-foot prev="Manual:External utilities:oct-vibrational_spectrum" next="Manual:Deprecated_Utilities" >}}
---------------------------------------------
---
title: "oct-run periodic table"
#series: "Manual"
---


### NAME 
oct-run_periodic_table - run octopus for atoms

### SYNOPSIS 
{{< command  " oct-run_periodic_table [ option ] ..." >}}

### DESCRIPTION 
This script is one of the octopus utilities.

Simple script to run octopus for all atoms in the periodic table for
which pseudopotentials are available.

{{< flag "-h" >}}
Show this message

{{< flag "-n" >}}
Run in dry run mode (show what would be executed)

{{< flag "-s species" >}}
Run octopus for given species.

{{< flag "-a" >}}
Run octopus for all atoms

{{< flag "-r" >}}
Do a restart, i.e. do not use {{< variable "fromScratch" >}}=yes.

{{< flag "-t temperature" >}}
Run with electronic temperature.

{{< flag "-e no_states" >}}
Use extra states

{{< flag "-x octopus_executable" >}}

{{< flag "-d defaults_file_for_pseudopotentials" >}}

### EXAMPLES 
The following command would run the octopus executable to
calculate the He atom.

{{< command-line " oct-run_periodic_table -s He -x octopus" >}}

{{< manual-foot prev="Manual:External utilities:oct-propagation_spectrum" next="Manual:External utilities:oct-run_regression_test.pl" >}}
---------------------------------------------
---
title: "oct-infrared spectrum"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-photoelectron spectrum"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-run regression test.pl"
#series: "Manual"
---


<!DOCTYPE html>
<html><head><title>Forbidden</title></head>
<body><h1>Forbidden</h1><p>Invalid file extension found in the path info or query string.</p></body></html>

---------------------------------------------
---
title: "oct-center-geom"
#series: "Manual"
---


### NAME 
oct-center-geom - Centers a molecule's geometry

### SYNOPSIS 

```bash
oct-center-geom
```

[oct-center-geom does not read the standard input: all standard input
will be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-center-geom accepts no command-line
arguments, since there is not a standard way to do this with Fortran
90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

It reads the coordinates defined in the {{< file "inp" >}} file, and constructs an output xyz file, that will be called
{{< file "adjusted.xyz" >}} file, that describes the same system but in which the
atomic coordinates are centered, and (optionally) has the axes aligned.

To control the orientation of the centered molecule there are two parameters: {{< variable "MainAxis" >}} and {{< variable "AxisType" >}}.

Do not forget then to change your input file to use this file instead of your old geometry (by changing {{< variable "XYZCoordinates" >}}={{< value "'adjusted.xyz'" >}}). 

Be careful with units, this utility honours the {{< variable "Units" >}}, {{< variable "UnitsInput" >}} and {{< variable "UnitsOutput" >}} variables, so the {{< file "adjusted.xyz" >}} file will be in the specified output units.

{{< manual-foot prev="Manual:External utilities:oct-casida_spectrum" next="Manual:External utilities:oct-check_deallocs" >}}
---------------------------------------------
---
title: "oct-vdW c6"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-local multipoles"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-oscillator-strength"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-check deallocs"
#series: "Manual"
---



---------------------------------------------
---
title: External Utilities
weight: 40
---

{{% children depth=3 %}}

---
title: "oct-vibrational spectrum"
#series: "Manual"
---


This utility calculates the vibrational spectrum from a molecular-dynamics run. 

What this utility does is to read the velocity from the {{< file "td.general/coordinates" >}} and calculate the Velocity Autocorrelation Function:

$$
C\_{v}(t)=\\sum\_{i=1}^{N\_{atoms}}\\vec{v}\_i(t)\\cdot\\vec{v}\_i(t\_0)\\ ,
$$

afterward a cosinusoidal envelope is added, to make the periodic extension of the function continuous and then the spectrum is calculated by taking the Fourier transform of the function. On exit, two files are generated {{< file "td.general/velocity_autocorrelation" >}} and {{< file "td.general/vibrational" >}}.

This utility honours the variables {{< variable "PropagationSpectrumStartTime" >}} and {{< variable "PropagationSpectrumEndTime" >}} to control the time of sampling. Note that the velocity in the initial time must be different from zero, or $C_{v}$ will be identically zero.

As a discrete Fourier tranform is used, this utility can take several minutes to process a large run. If {{< octopus >}} was compiled with OpenMP support, this utility can run in several threads.

{{< manual-foot prev="Manual:External utilities:oct-vdW_c6" next="Manual:External utilities:oct-xyz-anim" >}}
---------------------------------------------
---
title: "oct-propagation spectrum"
#series: "Manual"
---


### NAME 
oct-propagate_spectrum - Calculates the absorption cross section tensor from the results of a time-propagation run.

### SYNOPSIS 
{{< command "oct-propagate_spectrum" >}}

[oct-propagate_spectrum does not read the standard input: all standard input
will be simply ignored. An input file named {{< file "inp" >}} must be present in the
running directory. Also, oct-propagate_spectrum accepts no command line
arguments, since there is not a standard way to do this with Fortran
90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

This utility generates the dipole strength function of the given system. Its main input is the td.general/multipoles file. Output is written to a file called spectrum. This file is made of two columns: energy (in eV or a.u., depending on the units specified in the input file), and dipole strength function (in 1/eV, or 1/a.u., idem).

In the input file, the user may set the {{< variable "PropagationSpectrumTransform" >}} (this should be set to “sine” for proper use), the {{< variable "PropagationSpectrumDampMode" >}} (recommended value is “polynomial”, which ensures fulfilling of the N-sum rule), the {{< variable "PropagationSpectrumStartTime" >}}, the {{< variable "PropagationSpectrumEndTime" >}}, the {{< variable "PropagationSpectrumEnergyStep" >}}, and the {{< variable "PropagationSpectrumMaxEnergy" >}}.

{{< manual-foot prev="Manual:External utilities:oct-photoelectron_spectrum" next="Manual:External utilities:oct-run_periodic_table" >}}
---------------------------------------------
---
title: "oct-conductivity"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-convert"
#series: "Manual"
---


### Name 
oct-convert - Octopus utility to read obj files and write in many different file formats

### Description 

This executable gives the ability to read files written during the ground-state or time-dependent execution

### Example 

You can run ground-state and time-dependent execution of the [benzene example](../Benzene_molecule).

Then, we have to add this to the inp file, if we want to have the ground state density in DX format:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertFolder" >}} = 'restart/gs'
 {{< variable "ConvertFilename" >}} = 'density'
 {{< variable "ConvertIterateFolder" >}} = no
```

To convert the restart wave-functions (from 1 to 10) of a td run:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertIterateFolder" >}} = no
 {{< variable "ConvertFilename" >}} = ' '
 {{< variable "ConvertStart" >}} = 1
 {{< variable "ConvertEnd" >}}   = 10
 {{< variable "ConvertFolder" >}} = 'restart/td/'
```

If we want to convert the densities of the time-dependent executions, from files td.0000001 to td.0000010:

```bash
 {{< variable "Output" >}} = density
 {{< variable "OutputFormat" >}} = dx
 {{< variable "ConvertIterateFolder" >}} = yes
 {{< variable "ConvertStart" >}} = 1
 {{< variable "ConvertEnd" >}}   = 10
```


{{< manual-foot prev="Manual:External utilities:oct-conductivity" next="Manual:External utilities:oct-dielectric-function" >}}
---------------------------------------------
---
title: "oct-help"
#series: "Manual"
---


Oct-help is a utility to helps you to write input files. It has been available since version 3.0. Oct-help has two modes:

####  search  

```bash
oct-help search string
```

In this mode oct-help will print all variable names that contain ''string'' in their name.

####  show  

```bash
oct-help show variable
```

In this mode oct-help will print the documentation for the indicated variable.

{{< manual-foot prev="Manual:External utilities:oct-harmonic-spectrum" next="Manual:External utilities:oct-infrared_spectrum" >}}
---------------------------------------------
---
title: "oct-display partitions"
#series: "Manual"
---


### NAME 
oct-display_partitions - Graphical representation of the mesh partitions.

### SYNOPSIS 

```bash
oct-display_partitions
```

oct-display_partitions does not read the standard input: all standard
input will be simply ignored. Also, oct-display_partitions accepts no
command line arguments.

### DESCRIPTION 
This program is one of the octopus utilities.

It is a script to use gnuplot to plot the partitioning of the mesh. The
files {{< emph TODO >}} {{< file "mesh_partition.???" >}} as produced in {{< file "debug/mesh_partition" >}} by the
{{< octopus >}} debug mode have to be present in the current working directory
when this script is invoked. The plot can be found in the file
{{< file "mesh_partitions.png" >}}. This script generates a gnuplot-script called
{{< file "mesh_partitions_index.gp" >}} which is stored in the current working
directory and can be loaded into gnuplot manually. With:

```bash
gnuplot> load mesh_partitions_index.gp
gnuplot> set term x11
gnuplot> replot
```

the plot can be reproduced and shown on the screen so that rotating and
zooming is possible.

{{< manual-foot prev="Manual:External utilities:oct-dielectric-function" next="Manual:External utilities:oct-harmonic-spectrum" >}}
---------------------------------------------
---
title: "oct-harmonic-spectrum"
#series: "Manual"
---


### NAME 
oct-harmonic-spectrum - Calculates the harmonic spectrum.

### SYNOPSIS 

```bash
oct-harmonic-spectrum
```

oct-harmonic-spectrum does not read the standard input: all standard
input will be simply ignored. An input file named {{< file "inp" >}} must be present in
the running directory. Also, oct-harmonic-spectrum accepts no command
line arguments, since there is not a standard way to do this with
Fortran 90.

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

{{< notice note >}} 
TODO: MISSING SHORT DESCRIPTION OF THE PROGRAM
{{< /notice >}}

{{< manual-foot prev="Manual:External utilities:oct-display_partitions" next="Manual:External utilities:oct-help" >}}

---------------------------------------------
---
title: "oct-analyze projections"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-run testsuite.sh"
#series: "Manual"
---

### NAME

{{< code oct-run_testsuite >}} - Script to run the octopus testsuite.
### SYNOPSIS

```bash
oct-run_testsuite [ option ] ...
```

### DESCRIPTION

This script is one of the octopus utilities.

This script runs the octopus testsuite.


{{< flag "-h" >}}

Prints out a short help message.

{{< flag "-n" >}}

Run in dry mode (show what would be executed).

{{< flag "-f" >}}

Run full testsuite.

{{< flag "-l" >}}

Local run.

{{< flag "-p prefix" >}}

Installation prefix [default: /usr]

{{< flag "-m address" >}}

Mail report to address


### EXAMPLES


---------------------------------------------
---
title: "oct-atomic occupations"
#series: "Manual"
---


This script prints out to standard output the electronic configuration of a given atom, and the way in which this electronic configuration should be reflected in the {{< variable "Occupations" >}} block of a corresponding {{< octopus >}} input file.

### Options 

{{< flag "-s species" >}}
species should be the atomic symbol (e.g. Na, Au, etc).

{{< flag "-h" >}}
Show a brief summary of command line options.

### Examples 
```bash
oct-atomic_occupations -s Na
```

```bash
oct-atomic_occupations -s Ti_sc
```

```bash
$ for x in \$(cat /usr/share/octopus/PP/defaults | awk '{print \$1}')
> do oct-atomic_occupations -s $x
> done
```

{{< manual-foot prev="Manual:External utilities:oct-analyze_projections" next="Manual:External utilities:oct-casida_spectrum" >}}
---------------------------------------------
---
title: "oct-dielectric-function"
#series: "Manual"
---



---------------------------------------------
---
title: "oct-casida spectrum"
#series: "Manual"
---


### NAME 
oct-casida_spectrum - Broadens the linear-response spectra from the Casida run mode.

### SYNOPSIS 

```bash
oct-casida_spectrum
```

[oct-casida_spectrum does not read standard input: all standard input will be
simply ignored. An input file named {{< file "inp" >}} must be present in the running
directory. Also, oct-casida_spectrum accepts no command line arguments, since
there is not a standard way to do this with Fortran 90.]

### DESCRIPTION 
This program is one of the {{< octopus >}} utilities.

It post-processes the files {{< file "casida/casida" >}}, {{< file "casida/petersilka" >}} and
{{< file "casida/eps_diff" >}}, that contain the excitation energies and oscillator strengths.

The parameters of the spectrum can be set using the variables {{< variable "CasidaSpectrumBroadening" >}}, {{< variable "CasidaSpectrumMinEnergy" >}}, {{< variable "CasidaSpectrumMaxEnergy" >}}, and {{< variable "CasidaSpectrumEnergyStep" >}}.

{{< manual-foot prev="Manual:External utilities:oct-atomic_occupations" next="Manual:External utilities:oct-center-geom" >}}
---------------------------------------------
---
title: Octopus Manual
section: "Manual"
---

{{% children depth=3 %}}

---
title: "Copying"
section: "Manual"
---


###  INTRODUCTION  

The real-space TDDFT code octopus ("octopus") is provided under the 
GNU General Public License ("GPL"), Version 2, with exceptions for 
external libraries that are contained solely for convenience in this 
distribution. 

You can find a copy of the GPL license [https://www.gnu.org/copyleft/gpl.html  here. ]

A copy of the exceptions and licenses follow this introduction.

###  LICENSE EXCEPTIONS  

Octopus provides several external libraries, which are located in the 
"external_libs" subdirectory of the octopus source distribution. The 
GNU General Public License does not apply to these libraries. Separate 
copyright notices can be found for each library in the respective 
subdirectories of "external_libs". Copyright notices are also contained 
in this document.

Currently the following external libraries are provided:

###  Expokit  
*Roger B. Sidje
*Department of Mathematics, University of Queensland 
*Brisbane, QLD-4072, Australia
*Email: rbs@maths.uq.edu.au
*WWW: https://www.maths.uq.edu.au/expokit/ 
*Copyright: https://www.maths.uq.edu.au/expokit/copyright
```text
 
```

###  Metis 4.0  

*George Karypis 
*Department of Computer Science & Engineering
*Twin Cities Campus
*University of Minnesota, Minneapolis, MN, USA
*Email: karypis@cs.umn.edu
*WWW: https://www-users.cs.umn.edu/~karypis/metis/index.html

####  METIS COPYRIGHT NOTICE  

The ParMETIS/METIS package is copyrighted by the Regents of the
University of Minnesota. It can be freely used for educational and
research purposes by non-profit institutions and US government agencies
only.  Other organizations are allowed to use ParMETIS/METIS only for
evaluation purposes, and any further uses will require prior approval.
The software may not be sold or redistributed without prior approval.
One may make copies of the software for their use provided that the
copies, are not sold or distributed, are used under the same terms
and conditions.

As unestablished research software, this code is provided on an
``as is'' basis without warranty of any kind, either expressed or
implied. The downloading, or executing any part of this software
constitutes an implicit agreement to these terms. These terms and
conditions are subject to change at any time without prior notice.

###  qshep  

*Robert Renka
*University of North Texas
*(817) 565-2767


####  QSHEP COPYRIGHT NOTICE  

ACM Software Copyright Notice

[https://www.acm.org/pubs/toc/CRnotice.html]

Copyright   (c)  1998  Association   for  Computing   Machinery,  Inc.
Permission to  include in application  software or to make  digital or
hard copies  of part or all of  this work is subject  to the following
licensing agreement.

ACM Software License Agreement

All software, both binary and  source published by the Association for
Computing  Machinery  (hereafter,  Software)  is  copyrighted  by  the
Association  (hereafter, ACM) and  ownership of  all right,  title and
interest in and to the Software remains with ACM.  By using or copying
the Software, User agrees to abide by the terms of this Agreement.

Noncommercial Use

The ACM  grants to you (hereafter, User)  a royalty-free, nonexclusive
right  to execute,  copy, modify  and distribute  both the  binary and
source  code   solely  for   academic,  research  and   other  similar
noncommercial uses, subject to the following conditions:

1. User  acknowledges that the  Software is  still in  the development
stage  and that  it is  being supplied  "as is,"  without  any support
services   from  ACM.    Neither  ACM   nor  the   author   makes  any
representations or warranties,  express or implied, including, without
limitation, any  representations or warranties  of the merchantability
or fitness for any particular  purpose, or that the application of the
software, will not infringe on any patents or other proprietary rights
of others.

2. ACM shall  not be held  liable for direct, indirect,  incidental or
consequential  damages arising  from any  claim by  User or  any third
party with respect  to uses allowed under this  Agreement, or from any
use of the Software.

3. User agrees  to fully  indemnify and hold  harmless ACM  and/or the
author(s) of  the original work from  and against any  and all claims,
demands, suits, losses, damages, costs and expenses arising out of the
User's use of the Software, including, without limitation, arising out
of the User's modification of the Software.

4. User may modify  the Software and distribute that  modified work to
third  parties provided  that: (a)  if posted  separately,  it clearly
acknowledges  that it  contains  material copyrighted  by  ACM (b)  no
charge is associated  with such copies, (c) User  agrees to notify ACM
and the Author(s)  of the distribution, and (d)  User clearly notifies
secondary users that such modified work is not the original Software.

5. User agrees that  ACM, the authors of the  original work and others
may enjoy  a royalty-free, non-exclusive license to  use, copy, modify
and redistribute these modifications to  the Software made by the User
and  distributed to  third parties  as  a derivative  work under  this
agreement.

6. This agreement will terminate immediately upon User's breach of, or
non-compliance with,  any of its terms.   User may be  held liable for
any  copyright   infringement  or   the  infringement  of   any  other
proprietary rights  in the Software  that is caused or  facilitated by
the User's failure to abide by the terms of this agreement.

7. This agreement  will be construed  and enforced in  accordance with
the law  of the  state of New  York applicable to  contracts performed
entirely  within the State.   The parties  irrevocably consent  to the
exclusive jurisdiction of  the state or federal courts  located in the
City of New York for all disputes concerning this agreement.

Commercial Use

Any User wishing to make a commercial use of the Software must contact
ACM  at   permissions@acm.org  to  arrange   an  appropriate  license.
Commercial use  includes (1) integrating or incorporating  all or part
of the source code into a product for sale or license by, or on behalf
of, User to third parties, or (2) distribution of the binary or source
code  to third  parties  for use  with  a commercial  product sold  or
licensed by, or on behalf of, User.

Revised 6/98

{{< manual-foot prev="Manual:Appendix:Reference Manual" next="Manual" >}}
---------------------------------------------
---
title: Spack
draft: true
---
---
title: EasyBuild
draft: true
---
---
Title: "Direct installation using system compilers"
Weight: 5
---

---
title: "Specific architectures"
#series: "Manual"
---

{{% notice note %}}
This page needs updating!
{{% /notice %}}

Here the {{< manual "Install" "general instructions" >}} for building Octopus are supplemented by instructions for some specific large supercomputers. Included also is how to configure the optional ETSF_IO library, and how you can run the testsuite in parallel.

### Generic Ubuntu 10.10

Install the following packages: 
```bash
 * Required: gfortran, gcc, make, automake, m4, libtool, libgsl0-dev, libblas-dev, liblapack-dev, libfftw3-dev
 * Optional: mpi-default-dev, libgd2-xpm-dev, libsparskit-dev, libnetcdf-dev, libtrilinos-dev
```

[Unfortunately you should not use <tt>libblacs-mpi-dev, libscalapack-mpi-dev</tt> packages because they can give incorrect results. To use BLACS and ScaLAPACK (optional) you must build them yourself.]

First [download libxc](https://www.tddft.org/programs/libxc/download), extract with {{< code "tar xzf" >}}, enter the resulting directory, and execute:

```bash
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3
 make install
 make check
```

For octopus, download and extract, and execute the following (inserting the path where you put libxc):

```bash
 ./configure --prefix=`pwd` CC=gcc CXX=g++ FC=gfortran FCCPP="/lib/cpp -ansi -freestanding" FCFLAGS="-O3 -ffree-line-length-none" CFLAGS=-O3 --with-libxc-prefix=''LIBXC_PATH''
 make install
```

If you are using MPI, replace the compilers by <tt>CC=/usr/bin/mpicc CXX=/usr/bin/mpicxx FC=/usr/bin/mpif90</tt> and add <tt>--enable-mpi</tt>.

To build ETSF_IO (optional), <tt>libnetcdf-dev</tt> is required.

```bash
 ./configure --prefix=`pwd` CC=gcc FC=gfortran FCFLAGS="-O3 -ffree-line-length-none" CFLAGS="-O3" --with-netcdf-ldflags="-lnetcdff"
 make install
```

Add <tt>--with-etsf-io-prefix="$DIR/etsf_io-1.0.4"</tt> to the Octopus configure line, where <tt>$DIR/etsf_io-1.0.4</tt> is where you installed ETSF_IO.

(Aug 2018)

<!--
### Generic MacOS

Octopus is now available in MacPorts.
For more info, and how to build by hand with MacPorts libraries, see [[Compiling_octopus_in_OS_X]].
-->

### United States


{{% expand "Sequoia/Vulcan (IBM Blue Gene/Q)" %}}
Supercomputers at Lawrence Livermore National Laboratory based on the IBM Blue Gene/Q architecture. This was tested in Vulcan, but it should work on Sequoia as well.

https://computing.llnl.gov/tutorials/bgq/

```bash
 export LIBS_BLAS="-L/usr/local/tools/essl/5.1/lib/ -lesslsmpbg"
 export LIBS_LAPACK="/usr/local/tools/lapack/lib/liblapack.a"
 export LIBS_FFT="-L/usr/local/tools/fftw-3.3.3/lib/ -lfftw3_omp"
 export FC_INTEGER_SIZE=4
 export CC_FORTRAN_INT=int
 export CC=mpixlc_r
 export FC=mpixlf95_r
 export LDFLAGS="-qsmp=omp"
 export CFLAGS="-g -O3"
 export FCFLAGS=$CFLAGS" -qxlf90=autodealloc -qessl"
 ./configure --with-libxc-prefix=$HOME --prefix=$HOME --host=powerpc32-unknown-linux-gnu --build=powerpc64-unknown-linux-gnu --with-gsl-prefix=$HOME --enable-openmp --enable-mpi
```

(This build does not use Scalapack)
{{% /expand %}}


{{% expand "Stampede" %}}
10 PFLOPS (PF) Dell Linux Cluster based on 6,400+ Dell PowerEdge server nodes, each outfitted with 2 Intel Xeon E5 (Sandy Bridge) processors and an Intel Xeon Phi Coprocessor (MIC Architecture). Texas Advanced Computing Center, University of Texas, Austin, USA.

https://portal.tacc.utexas.edu/user-guides/stampede

```bash
 module swap mvapich2 impi
 module load fftw3 gsl hdf5 netcdf
 MKL_DIR=$MKLROOT
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3" --enable-mpi \
 --with-blas="-L$MKL_DIR/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread" --with-fftw-prefix="$TACC_FFTW3_DIR" \
 --with-netcdf-prefix="$TACC_NETCDF_DIR" --enable-newuoa --with-libxc-prefix=`pwd`/../libxc-2.2.2 --with-blacs="$MKL_DIR/lib/intel64/libmkl_blacs_intelmpi_lp64.a" \
 --with-scalapack="$MKL_DIR/lib/intel64/libmkl_scalapack_lp64.a"
 make -j 6 install
```

For libxc:
```bash
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 FCFLAGS="-O3" CFLAGS="-O3"
 make -j 6 install
 make check
```

testsuite script:

```bash
 #!/usr/bin/env bash
 
 #SBATCH -n 16
 #SBATCH -p development
 #SBATCH -t 01:00:00
 #SBATCH --export=ALL
 
 module load perl
 WORKDIR=$PWD
 export TEMPDIRPATH=$SCRATCH/tmp
 cd $HOME/octopus
 export OCT_TEST_NJOBS=8
 export MPIEXEC=`which ibrun`
 make check &> $WORKDIR/makecheck
```

(13 September 2016)
{{% /expand %}}

{{% expand "Edison" %}}
Cray XC30 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

https://www.nersc.gov/nusers/systems/

```bash
 module load cray-fftw gsl cray-netcdf-hdf5parallel parpack
 
 ./configure --prefix=`pwd` --with-libxc-prefix=/usr/common/usg/libxc/3.0.0/intel/ivybridge --enable-mpi \
 CC=cc CXX=CC FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo" \
 FCCPP="/lib/cpp -ffreestanding" --with-fftw-prefix=$FFTW_DIR/..
 --with-arpack="$ARPACK" --with-parpack=no \
 --with-berkeleygw-prefix=/usr/common/software/berkeleygw/2.0
 make -j 6 install
```

To run the testsuite in parallel:

```bash

#!/bin/bash -l
#SBATCH -J pulpo
#SBATCH -N 2
#SBATCH -p regular
#SBATCH -t 04:00:00
#SBATCH --export=ALL

cd $HOME/edison/octopus
export OMP_NUM_THREADS=1
export MPIEXEC=`which srun`
export TEMPDIRPATH=$SCRATCH/tmp
export OCT_TEST_NJOBS=15
make check &> $SLURM_SUBMIT_DIR/makecheck
```

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```bash
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo" \
 --with-netcdf-module-path="$NETCDF_DIR/include" --with-netcdf-prefix="$NETCDF_DIR"
 make -j 6 install
 make check
```

(version 8.1, July 2018)
{{% /expand %}}

{{% expand "Cori" %}}
Cray XC40 at National Energy Research Scientific Computing Center (NERSC), Lawrence Berkeley National Laboratory, USA

https://www.nersc.gov/users/computational-systems/cori/

```bash
 module load cray-fftw gsl cray-netcdf-hdf5parallel parpack
```

For libxc (you can use the installed version instead)
```bash
 ./configure --prefix=`pwd` CC=cc FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

For octopus (note: using parpack here causes scalapack problems):
```bash
 ./configure --prefix=`pwd` --with-libxc-prefix=/usr/common/software/libxc/4.2.3/intel/haswell --enable-mpi \
 CC=cc CXX=CC FC=ftn FCFLAGS="-fast -no-ipo" CFLAGS="-fast -no-ipo" \
 FCCPP="/lib/cpp -ffreestanding" --with-fftw-prefix=$FFTW_DIR/.. \
 --with-arpack="$ARPACK" --with-parpack=no --with-berkeleygw-prefix=/usr/common/software/berkeleygw/2.0
 make -j 6 install
```

To run the testsuite in parallel:

```bash

-!/bin/bash -l
-SBATCH -J pulpo
-SBATCH -n 64
-SBATCH -C haswell
-SBATCH -p regular
-SBATCH -t 04:00:00
-SBATCH --export=ALL

cd $HOME/cori/octopus
export OMP_NUM_THREADS=1
export MPIEXEC=`which srun`
export TEMPDIRPATH=$SCRATCH/tmp
export OCT_TEST_NJOBS=15
make check &> $SLURM_SUBMIT_DIR/makecheck
```

To build ''ETSF_IO'' 1.0.4 (optional): (some tests fail)

```bash
 module load cray-netcdf-hdf5parallel
 ./configure --prefix=`pwd` FC=ftn FCFLAGS="-fast -no-ipo" CC=cc CFLAGS="-fast -no-ipo"
 make -j 6 install
 make check
```

(30 July 2018)
{{% /expand %}}


{{% expand "Lonestar" %}}
Dell Linux cluster at Texas Advanced Computing Center (TACC), University of Texas, Austin, USA

Part of National Science Foundation's [https://teragrid.org/ TeraGrid]

https://services.tacc.utexas.edu/index.php/lonestar-user-guide

```bash
 module load gsl
 module load mkl
 module load netcdf
 module load fftw3
 cd libxc
 ./configure --prefix=`pwd`/.. --enable-mpi CC=mpicc FC=mpif90 F77=mpif90 FCFLAGS=-O3 CFLAGS=-O3
 make install
 cd ..
 export SLKPATH=/opt/apps/intel11_1/mvapich2_1_6/scalapack/1.8.0/lib
 ./configure --prefix=`pwd` --enable-mpi CC=mpicc FC=mpif90 \
 --with-blas="-Wl,--start-group $TACC_MKL_LIB/libmkl_intel_lp64.a $TACC_MKL_LIB/libmkl_sequential.a $TACC_MKL_LIB/libmkl_core.a -Wl,--end-group -lpthread" \
 --with-fft-lib="-L$TACC_FFTW3_LIB -lfftw3" --with-netcdf-prefix="$TACC_NETCDF_DIR" --disable-gdlib \
 --with-etsf-io-prefix="$HOME/etsf_io-1.0.3" --enable-newuoa --with-libxc-prefix=`pwd` \
 --with-blacs="$SLKPATH/libscalapack.a $SLKPATH/blacs_MPI-LINUX-0.a $SLKPATH/blacsF77init_MPI-LINUX-0.a $SLKPATH/blacs_MPI-LINUX-0.a"
 make install
```

To build ''ETSF_IO'':

```bash
 module load netcdf
 ./configure --prefix=`pwd` FC=mpif90 --with-netcdf-module-path=$TACC_NETCDF_INC --with-netcdf-ldflags=-L$TACC_NETCDF_LIB FCFLAGS=-O3
 make install
 make check
```

To run the testsuite in parallel:

```bash
 #!/bin/bash
 
 #$ -N test_pulpo
 #$ -cwd
 #$ -pe 6way 12
 #$ -q development
 #$ -l h_rt=00:45:00
 #$ -V
 
 cd $HOME/octopus
 export TEMPDIRPATH=$SCRATCH/tmp
 if [ ! -d $TEMPDIRPATH ]; then
     mkdir $TEMPDIRPATH
 fi
 
 export MPIEXEC=`which ibrun`
 export OCT_TEST_NPROCS=1
 make check &> makecheck
```

(2 May 2011)
{{% /expand %}}

{{% expand "Lawrencium" %}}
Dell Linux cluster at Lawrence Berkeley National Laboratory, USA

Cluster available only to Lawrence Berkeley National Laboratory researchers

https://lrc.lbl.gov/html/Lawrencium.html

```bash
 module load icc/10.1.018
 module load ifort/10.1.018
 module load mkl/2011.1.107
 module load openmpi/1.4.3-intel
 module load netcdf/4.0-intel
 module load fftw/3.1.2-intel
 module load autoconf/2.65
 module load gsl/1.13-intel
 module load gdlib/2.0.35
 module load libtool/2.2.6b
```

```bash
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3 -check bounds" --enable-mpi \
 --with-blas="-L$MKLROOT/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread" \
 --with-fft-lib="-L/global/software/centos-5.x86_64/modules/fftw/3.1.2-intel/lib -lfftw3" --with-netcdf-prefix="$NETCDF" --with-netcdf-include="$NETCDF/include" --enable-newuoa \
 --with-libxc-prefix=`pwd` --with-blacs="$MKL_DIR/lib/intel64/libmkl_blacs_openmpi_lp64.a" --with-scalapack="$MKL_DIR/lib/intel64/libmkl_scalapack_lp64.a"
```

To build ''ETSF_IO'':

```bash
 module load netcdf/4.2-intel-p
 ./configure --prefix=`pwd` CC=mpicc FC=mpif90 CFLAGS="-O3" FCFLAGS="-O3" --with-netcdf-module-path=$NETCDF_DIR/include/ --with-netcdf-ldflags="-L$NETCDF_DIR/lib/ -lnetcdf"
```

(21 Aug 2012)
{{% /expand %}}

### Europe

{{% expand "Curie" %}}

##### PFFT 1.0.5 and BigDFT 1.6.0

The Curie supercomputer, owned by GENCI and operated into the TGCC by CEA, is the first French Tier0 system open to scientists through the French participation into the PRACE research infrastructure.

Curie is offering 3 different fractions of x86-64 computing resources for addressing a wide range of scientific challenges and offering an aggregate peak performance of 2 PetaFlops.

General information: https://www-hpc.cea.fr/en/complexe/tgcc-curie.htm

Compilation with PFFT version 1.0.5 and LibISF taken from BigDFT 1.6.0 (previously compiled):

```bash
 #!/bin/bash
 
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 WPREFIX=$WORKDIR/software
 PREFIX=$HOME/software
 export ISF_HOME="$WORKDIR/software/libisf"
 
 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$PREFIX/fftw/3.3/include -I$ISF_HOME/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PREFIX/pfft/1.0.5/lib -L$PREFIX/fftw/3.3/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 export LDFLAGS=$LDFLAGS" -L$ISF_HOME/lib -lisfsolver -lrt"
  
 ../configure --prefix=$WPREFIX/octopus/superciliosus_pfft \
   --disable-openmp --with-libxc-prefix=$WPREFIX/libxc/2.0.x \
   --disable-gdlib --with-gsl-prefix=/usr/local/gsl-1.14 \
   --enable-mpi --enable-newuoa \
   --with-pfft-prefix=$PREFIX/pfft/1.0.5
```

BigDFT was configured like this:

```bash
 #!/bin/bash 
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
   --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
   --prefix=$WORKDIR/software/bigdft
```


##### PFFT 1.0.7 and BigDFT 1.7.0

```bash
 #!/bin/bash

 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 WPREFIX=$WORKDIR/poisson/software
 WLPREFIX=$WORKDIR/poisson/local
 export ISF_HOME="$WLPREFIX/libisf"
 export FFTW_HOME="$WLPREFIX/fftw-3.3.3"
 export PFFT_HOME="$WLPREFIX/pfft-1.0.7-alpha"

 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$FFTW_HOME/include -I$ISF_HOME/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PFFT_HOME/lib -L$FFTW_HOME/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 
 export LDFLAGS=" -L/usr/lib64 -L$ISF_HOME/lib -lPSolver-1 -lrt"
 export LD_LIBRARY_PATH=$ISF_HOME/lib:$LD_LIBRARY_PATH
 
 ../configure --prefix=$WPREFIX/octopus/superciliosus_pfft \
   --disable-openmp --with-libxc-prefix=$WPREFIX/libxc/11066 \
   --disable-gdlib --with-gsl-prefix=/usr/local/gsl-1.14 \
   --enable-mpi --enable-newuoa \
   --with-pfft-prefix=$PFFT_HOME
```

LibISF compilation:
```bash
 #!/bin/bash
 
 module load c/intel
 module load fortran/intel
 module load bullxmpi
 module load mkl
 module load gsl/1.14
 
 export FC=mpif90
 
 ../configure --with-ext-linalg="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lm" \
    --with-ext-linalg-path="-L/usr/local/Intel_compilers/c/composer_xe_2011_sp1.7.256/mkl/lib/intel64" \
    --prefix=$WORKDIR/poisson/local/libisf \
    --disable-libbigdft --disable-binaries --enable-dynamic-libraries
```

FFTW 3.3.3 compilation script is in https://www-user.tu-chemnitz.de/~mpip/software/install_fftw-3.3.3_gcc.sh and it is changed as follow:

```bash
 < myprefix=$HOME/local
 --- 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 23c29,30
 < wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 ---
 > - wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 > cp ../fftw-$FFTW_VERSION.tar.gz .
 876c883
 < ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared \
 ---
 > ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-shared --disable-shared \ 
```


PFFTW 1.0.7 compilation script is in https://www-user.tu-chemnitz.de/~mpip/software/install_pfft-1.0.7-alpha_gcc.sh and it is changed as follow:

```bash
 < myprefix=$HOME/local
 ---
 > 
 > module load c/intel
 > module load fortran/intel
 > module load bullxmpi
 > module load mkl
 > module load gsl/1.14
 > 
 > myprefix=$WORKDIR/poisson/local
 24c31,32
 < wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 ---
 > -wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 > cp ../pfft-$PFFT_VERSION.tar.gz .
```

{{% /expand %}}

{{% expand "Fermi/Juqueen" %}}
IBM Blue Gene/Q, Cineca, Italy

General information: https://www.hpc.cineca.it/content/ibm-fermi-user-guide

The exactly the same script has been tested in the Juqueen, Forschungszentrum Jülich, Jülich Supercomputing Centre (JSC), Jülich, Germany. General information: https://www.fz-juelich.de/ias/jsc/EN/Expertise/Supercomputers/JUQUEEN/JUQUEEN_node.html

```bash
 #/bin/bash
 
 module load gsl
 module load lapack
 module load blacs
 module load scalapack 
 
 export PFFT_HOME="$HOME/local/pfft-threads-1.0.5-alpha"
 export FFTW3_HOME="$HOME/local/fftw-threads-3.3.2"
 export LD_LIBRARY_PATH=$PFFT_HOME/lib:$FFTW3_HOME/lib:$LD_LIBRARY_PATH 
 
 export FC_INTEGER_SIZE=8
 export CC_FORTRAN_INT=int
 export CC=mpicc
 export CFLAGS="-g -Wl,-relax -O3 -I$PFFT_HOME/include -I$FFTW3_HOME/include"
 export FC=mpif90
 export FCFLAGS=$CFLAGS" -qxlf90=autodealloc -qessl -qsmp=omp "
 export LIBS_BLAS="-lesslsmpbg -L/opt/ibmmath/essl/4.4/lib -lesslbg"
 export LIBS_LAPACK="-L$LAPACK_LIB -llapack"
 
 export LIBS_FFT="-L$FFTW3_HOME/lib/ -lfftw3_mpi -lfftw3 -lm"
 echo "GSL DIR = $GSL_HOME "
 echo "PFFT DIR = $PFFT_HOME "
 export LDFLAGS=$LDFLAGS" -R/opt/ibmcmp/lib/bg/bglib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/powerpc-bgp-linux/lib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/gnu-linux/lib/gcc/powerpc-bgp-linux/4.1.2 \
    -L/opt/ibmcmp/xlf/bg/11.1/bglib \
    -L/opt/ibmcmp/xlmass/bg/4.4/bglib \
    -L/opt/ibmcmp/xlsmp/bg/1.7/bglib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/runtime/SPI \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/sys/lib \
    -L/bgsys/drivers/V1R4M2_200_2010-100508P/ppc/comm/default/lib \
    -L/usr/local/zlib/v1.2.3/lib \
    -L/bgsys/drivers/ppcfloor/runtime/SPI \
    -L$PFFT_HOME/lib -L$FFTW3_HOME/lib \
    -lpfft -lfftw3_mpi -lfftw3 -lm \
    -ldl -lxl -lxlopt -lrt -lpthread"
 
 SVN_VERS=$(svn info .. | grep Revision | awk '{print $2}')
 echo $SVN_VERS
 
 ../configure --prefix=$HOME/software/octopus_omp \
    --with-libxc-prefix=$HOME/software/libxc_new \
    --with-fft-lib="$FFTW3_HOME/lib/libfftw3.a  $FFTW3_HOME/lib/libfftw3_omp.a" \
    --host=powerpc32-unknown-linux-gnu \
    --build=powerpc64-unknown-linux-gnu \
    --disable-gdlib \
    --disable-f90-forall \
    --with-blas="-L/opt/ibmmath/lib64 -lesslbg" \
    --with-lapack="-L$LAPACK_LIB -llapack" \
    --with-gsl-prefix=$GSL_HOME \
    --enable-mpi --enable-openmp
```

To build FFTW3. Small changes are made to the script taken from M. Pippig web site https://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
It has a patch and is build with MPI and multithreaded (multithreaded may not work in Juqueen):

```bash
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw-threads-${FFTW_VERSION}"
 TMP="tmp-mine-fftw-${FFTW_VERSION}" 
 
 # bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)" 
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
  fi 
 
 mkdir $TMP && cd $TMP 
 
 wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 # fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
 	  sbo[pe] = (int) (pe * (b * p->tblock) * vn);
 	  rbs[pe] = (int) (db * bt * vn);
 	  rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -	  if (sbs[pe] != (b * p->tblock) * vn
 -	      || rbs[pe] != (p->block * bt) * vn)
 +	  if (dbt != p->tblock
 +	      || db != p->block)
 	       equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC=mpixlc_r F77=mpixlf90_r MPICC=mpixlc_r MPICXX=mpixlcxx_r MPIFC=mpixlf90_r MPILIBS=" " \
 CFLAGS='-O3 -g' \
 FCFLAGS='-O3 -g' \ 
 
 make -j 4
 make install
```

To build PFFT, also taken from https://www-user.tu-chemnitz.de/~mpip/software.php?lang=en
```bash
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft-threads-$PFFT_VERSION
 FFTWDIR=$myprefix/fftw-threads-$FFTW_VERSION
 TMP="tmp-mine-pfft-$PFFT_VERSION" 
 
 # bash check if directory exists
 if [ -d $TMP ]; then
         echo "Directory $TMP already exists. Delete it? (y/n)"
 	read answer 
 	if [ ${answer} = "y" ]; then 
 		rm -rf $TMP 
 	else
 		echo "Program aborted."
 		exit 1
 	fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --build=powerpc64-bgq-linux-gnu --host=powerpc64-bgq-linux \
 --prefix=$INSTDIR --with-fftw3=$FFTWDIR \
 MPICC=mpif90 MPICXX=mpicxx MPIFC=mpif90 CC=mpicc FC=mpif90 \
 CFLAGS='-O3 -g -qmaxmem=-1' \
 FCFLAGS='-O3 -g -qmaxmem=-1' \
 --disable-shared
 
 make -j 4
 make install
```
{{% /expand %}}

{{% expand "Jugene" %}}

FZJ-JSC IBM Blue Gene/P, JÜLICH SUPERCOMPUTING CENTRE (JSC), Germany

May also work on other Blue Gene/P computers

General information: https://www.fz-juelich.de/jsc/jugene

Usage information: https://www.fz-juelich.de/jsc/jugene/usage/
```bash
 module load gsl;
 module load lapack;
 autoreconf -i;
 export FC_INTEGER_SIZE=4;
 export CC_FORTRAN_INT=int;
 export CC=mpixlc_r;
 export CFLAGS='-g -O3 -qarch=450d';
 export FC=mpixlf90_r;
 export FCFLAGS='-g -O3 -qarch=450d -qxlf90=autodealloc -qessl -qsmp=omp';
 export LIBS_BLAS="-lesslsmpbg -L/opt/ibmmath/essl/4.4/lib -lesslbg";
 export LIBS_LAPACK="-L$LAPACK_LIB -llapack";
 export LIBS_FFT="-L/bgsys/local/fftw3/lib/ -lfftw3 -lm";
 ./configure --prefix=$outdir \
      --host=powerpc32-unknown-linux-gnu \
      --build=powerpc64-unknown-linux-gnu \
      --disable-gdlib \
      --with-gsl-prefix=$GSL_DIR \
      --enable-mpi --enable-openmp;
```

To run you have to change to $WORK directory and run LoadLeveler. Below is an example of the tutorial (https://www.tddft.org/programs/octopus/wiki/index.php/Tutorial:Benzene_molecule) executing in 32 nodes in the SMP mode (4 threads per chip). The execution mode is hybrid; OpenMP/MPI. 

Here is a example job input script:
```bash
 # @ job_name = Octopus_Sample_1
 # @ comment = "BGP Job by Size"
 # @ error = $(job_name).$(jobid).out
 # @ output = $(job_name).$(jobid).out
 # @ environment = COPY_ALL;
 # @ wall_clock_limit = 00:30:00
 # @ notification = error
 # @ notify_user = foo@gmail.com
 # @ job_type = bluegene
 # @ bg_size = 32
 # @ queue
 mpirun -exe $OCTOPUS_HOME/bin/octopus_mpi -mode SMP -verbose 1 
```
To execute the above script:
```bash
 llsubmit executing_script
```
{{% /expand %}}

{{% expand "MareNostrum II" %}}

MareNostrum is a supercomputer in Europe, the most powerful in Spain. This is one of the seven supercomputers of the Spanish Supercomputing Network. The supercomputer consists of 2,560 JS21 blade computing nodes, each with 2 dual-core IBM 64-bit PowerPC 970MP processors running at 2.3 GHz for 10,240 CPUs in total. The computing nodes of MareNostrum communicate primarily through Myrinet.

You may need to compile GSL and FFTW libraries before starting then installation.

Here is a script to build Octopus:
```bash
 #!/bin/bash
 
 export LD_LIBRARY_PATH=/usr/lib
 export CC=mpicc
 export FC=mpif90
 
 export CFLAGS="-q64 -O3 -qtune=ppc970 -qarch=ppc970 -qcache=auto -qnostrict -qignerrno -qinlglue"
 export FCFLAGS=$CFLAGS" -qessl -qxlf90=autodealloc" 
 
 export LIBS_BLAS="-L/usr/lib64 -lessl" 
 
 ./configure \
    --with-lapack=/gpfs/apps/LAPACK/lib64/liblapack.a \
    --disable-gsltest \
    --disable-gdlib \
    --with-gsl-prefix=$outdirGSL \
    --with-fft-lib=$outdirFFT/lib/libfftw3.a \
    --prefix=$outdir --enable-mpi --disable-f90-forall
```

As you can see in the script the --with-gsl-prefix has to point to your GSL and --with-fft-lib to your FFTW. You should compile those with the same CFLAGS and FCFLAGS. Here is an example of configuring GSL (we use GSL-1.11):
```bash
 ./configure --prefix=$outdirGSL CC=gcc -m64 --enable-static --disable-shared
```
We use FFTW-3.1.3. To configure FFTW3 just write:
```bash
 ./configure --prefix=$outdirFFT
```

Recently Marenostrum has changed to the SLURM manager. To execute a job you have to submit it with jobsubmit (more details in the [https://www.bsc.es/media/859.pdf User's Guide]).

{{% /expand %}}

{{% expand "MareNostrum III" %}}

MareNostrum III is a supercomputer based on Intel SandyBridge processors, iDataPlex Compute Racks, a Linux Operating System and an Infiniband interconnection. More information: https://www.bsc.es/marenostrum-support-services/mn3

```bash
 module load MKL
 module load GSL
 
 export PREFIX=$HOME/Software
 export FFTW3_HOME=$HOME/local/fftw-3.3.2
 export PFFT_HOME=$HOME/local/pfft-1.0.5-alpha
 
 export MKL_LIBS=$MKL_HOME/lib/intel64
 export MKLROOT=$MKL_HOME
 export MKL_LIBS="-L$MKLROOT/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm"
 export CC=mpicc
 export CFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$ISF_HOME/include -openmp -I$MKLROOT/include"
 export FC=mpif90
 export FCFLAGS="-xHost -O3 -m64 -ip -prec-div -static-intel -sox -I$FFTW3_HOME/include -openmp -I$ISF_HOME/include -I$MKLROOT/include"
 
 export LIBS_BLAS="${MKL_LIBS}"
 export LIBS_FFT="-Wl,--start-group -L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3 -lfftw3_mpi -Wl,--end-group"
 
 
 ../configure --prefix=$PREFIX/octopus/ \
  --enable-openmp --with-libxc-prefix=$PREFIX/libxc \
  --disable-gdlib --with-gsl-prefix=/apps/GSL/1.15 \
  --enable-mpi --enable-newuoa \
  --with-pfft-prefix=$PFFT_HOME
```
{{% /expand %}}

{{% expand "Corvo" %}}

Corvo is the computational cluster of the Nano-bio Spectroscopy Group. It consists of 1504 processor cores and 5284 GiB of RAM connected by an Infiniband network.

Script to build LIBFM library included in the Scafacos library (modified Scafacos svn version 1920, this is valid since version 9887 of Octopus):

```bash
 ./configure --enable-fcs-solvers=fmm,direct --enable-fcs-fmm-comm=armci \
   --enable-fcs-fmm-unrolled --enable-fcs-fmm-max-mpol=40 \
   CFLAGS=-O3 CXXFLAGS=-O3 FCFLAGS=-O3 --disable-doc CXX=mpic++ CC=mpicc \
   FC=mpif90 --prefix=$HOME/Software/scafacos \
   LDFLAGS=-L/opt/scalapack/1.8.0/lib -L/opt/blacs/1.1/lib \
   -L/opt/gotoblas2/1.08/lib -L/opt/netcdf/4.0.1/lib -L/opt/lapack/3.2/lib \
   -L/opt/etsf_io/1.0.2/lib --enable-fcs-int=int --enable-fcs-float=double \
   --enable-fcs-integer=integer --enable-fcs-real=real*8 
 
 make
 make install
```

The scripts to compile FFTW 3.3.2 and PFFT 1.0.5

```sh
 #!/bin/sh -e 
 
 myprefix=$HOME/local
 FFTW_VERSION="3.3.2"
 INSTALLDIR="$myprefix/fftw/${FFTW_VERSION}"
 TMP="tmp-fftw-${FFTW_VERSION}"
 
 # bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 
 wget https://www.fftw.org/fftw-$FFTW_VERSION.tar.gz
 gzip -dc fftw-$FFTW_VERSION.tar.gz | tar xvf -
 cd fftw-$FFTW_VERSION
 
 # fix bug in fftw alltoall routine that causes deadlocks
 patch -b -p0 <<\EOF
 --- mpi/transpose-alltoall.c
 +++ mpi/transpose-alltoall.c
 @@ -223,8 +223,8 @@
           sbo[pe] = (int) (pe * (b * p->tblock) * vn);
           rbs[pe] = (int) (db * bt * vn);
           rbo[pe] = (int) (pe * (p->block * bt) * vn);
 -         if (sbs[pe] != (b * p->tblock) * vn
 -             || rbs[pe] != (p->block * bt) * vn)
 +         if (dbt != p->tblock
 +             || db != p->block)
               equal_blocks = 0;
      }
      pln->send_block_sizes = sbs;
 EOF
 
 ./configure --prefix=$INSTALLDIR --enable-mpi --enable-openmp --enable-threads --disable-shared \
 CC="mpicc" F77="mpif90" MPICC="mpicc" \
 CFLAGS="-O3" FFLAGS="-O3" \
 MPILIBS=" "
 
 make -j 4
 make install
```


```bash
 #!/bin/sh -e
 
 myprefix=/home/local
 PFFT_VERSION=1.0.5-alpha
 FFTW_VERSION=3.3.2
 INSTDIR=$myprefix/pfft/$PFFT_VERSION
 FFTWDIR=$myprefix/fftw/$FFTW_VERSION
 TMP="tmp-pfft-$PFFT_VERSION"
 module unload fftw/3.2.2
 
 # bash check if directory exists
 if [ -d $TMP ]; then
        echo "Directory $TMP already exists. Delete it? (y/n)"
        read answer
        if [ ${answer} = "y" ]; then
                rm -rf $TMP
        else
                echo "Program aborted."
                exit 1
        fi
 fi
 
 mkdir $TMP && cd $TMP
 export CC=mpicc
 export FC=mpif90 
 
 wget https://www.tu-chemnitz.de/~mpip/software/pfft-$PFFT_VERSION.tar.gz
 gzip -dc pfft-$PFFT_VERSION.tar.gz | tar xvf -
 cd pfft-$PFFT_VERSION
 ./configure --prefix=$INSTDIR --with-fftw3=$FFTWDIR --disable-shared
 
 make install
```

This is the scrip to compile Octopus with PFFT and FMM support (through Scafacos library):

```bash
 module add gcc
 module add ifort
 module add gsl
 module add etsf_io
 module add lapack
 module add netcdf
 module add gotoblas2
 module add mpibull2 
 module add blacs
 module add scalapack
 
 export CC=mpicc
 export FC="mpif90"
 
 # It is important that "/opt/intel/Compiler//11.1/038/bin/intel64/" not to be in the PATH
 export PATH=/opt/mpi/mpibull2-1.3.9-14.s/bin:/opt/netcdf/4.0.1/bin:/opt/etsf_io/1.0.2/bin:/opt/gsl//1.13/bin:/opt/intel/composer_xe_2011_sp1.10.319//bin/intel64/:/opt/gcc//4.4/bin:/opt/slurm/bin:/usr/kerberos/bin:/opt/cuda//bin:/usr/local/bin:/bin:/usr/bin
 
 export LIBFM_HOME="$HOME/Software/scafacos"
 export PFFT_HOME="$HOME/local/pfft-1.0.5-alpha"
 export FFTW3_HOME="$HOME/local/fftw-3.3.2"
 export LIBS_PFFT="-L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3_mpi -lfftw3"
 export CFLAGS=" -I$PFFT_HOME/include -I$FFTW3_HOME/include -I$LIBFM_HOME/include"
 export LDFLAGS=$LDFLAGS" -L$LIBFM_HOME/lib -L$PFFT_HOME/lib -L$FFTW3_HOME/lib -lpfft -lfftw3_mpi -lfftw3 -lm "
 export FCFLAGS=$CFLAGS
 export CPPFLAGS=" -I$LIBFM_HOME/include "
 export LD_LIBRARY_PATH=$LIBFM_HOME/lib:$LD_LIBRARY_PATH

 ../configure --enable-mpi --with-libxc-prefix=$HOME/Software/libxc_9882 \
    --prefix=$HOME/Software/octopus \
    --with-blacs=/opt/blacs/1.1/lib/libblacs.a --with-scalapack \
    --with-libfm="-L$LIBFM_HOME/lib -lfcs4fortran -lfcs -lfcs_direct -lfcs_fmm -lfcs_near -lfcs_gridsort -lfcs_common" \
    --disable-openmp
```
{{% /expand %}}

{{% expand "LaPalma" %}}

LaPalma is a supercomputer of the IAC. It is a part of the old MareNostrum II computer. It was not possible to compile with XL compilers, so GNU 4.6.1 version is used.

Previously to Octopus, Libxc (2.2), FFTW (3.3.4), GSL (1.16) and OpenBLAS (0.2.13) were required to compile.

GSL
```bash
 export LD_LIBRARY_PATH=/usr/lib:/gpfs/apps/MPICH2/slurm-2.5.6/64/lib/:/gpfs/apps/GCC/4.9.2/lib64/
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.9.2/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.9.2/bin/gfortran"
 
 ./configure --prefix=$HOME/local/gsl/gnu --enable-static --disable-shared
```

OpenBLAS
```bash
 make PREFIX=$HOME/local/openblas/4.6.1 BINARY=64 CC=/gpfs/apps/GCC/4.6.1/bin/gcc FC=/gpfs/apps/GCC/4.6.1/bin/gfortran
```
Libxc
```bash
 export PATH=/gpfs/apps/MPICH2/mx/default/64/bin:$PATH
 export LD_LIBRARY_PATH=/gpfs/apps/MPICH2/mx/default/64/lib:/gpfs/apps/MPICH2/slurm/64/lib
 
 export LD_LIBRARY_PATH=/gpfs/apps/GCC/4.6.1/lib64/:$HOME/local/openblas/lib:$LD_LIBRARY_PATH
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.6.1/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.6.1/bin/gfortran"
 
 export CFLAGS="-m64 -O3"
 export FCFLAGS=$CFLAGS
 
 ./configure --prefix=$HOME/local/libxc/gnu/4.6.1/2.2
```

Octopus 

```bash
 export PATH=/gpfs/apps/MPICH2/mx/default/64/bin:$PATH
 export LD_LIBRARY_PATH=/gpfs/apps/MPICH2/mx/default/64/lib:/gpfs/apps/MPICH2/slurm/64/lib
 
 export LD_LIBRARY_PATH=/gpfs/apps/GCC/4.6.1/lib64/:$HOME/local/openblas/4.6.1/lib:$LD_LIBRARY_PATH
 export CC="/gpfs/apps/MPICH2/mx/default/64/bin/mpicc -cc=/gpfs/apps/GCC/4.6.1/bin/gcc"
 export FC="/gpfs/apps/MPICH2/mx/default/64/bin/mpif90 -fc=/gpfs/apps/GCC/4.6.1/bin/gfortran"
 
 export CFLAGS="-m64 -O2 -I $PWD/external_libs/isf/wrappers "
 export FCFLAGS=$CFLAGS
 
 ../configure \
    --with-blas=$HOME/local/openblas/4.6.1/lib/libopenblas.so.0 \
    --with-lapack=/gpfs/apps/LAPACK/lib64/liblapack.a \
    --disable-gdlib \
    --with-gsl-prefix=$HOME/local/gsl/gnu \
    --with-fft-lib=$HOME/local/fftw/3.3.4/lib/libfftw3.a \
    --with-libxc-prefix=$HOME/local/libxc/gnu/4.6.1/2.2 \
    --prefix=/gpfs/projects/ehu31/software/octopus/gnu/4.6.1 \
    --enable-mpi --disable-f90-forall
```
{{% /expand %}}

{{% expand "XBES" %}}

Located at the University of Hamburg in Hamburg. It has a login node and a 24 nodes with Intel Xeon. Here a basic configuration script:

```bash
 ./configure \
 PATH="$PATH:/home/common/lib/gsl-1.16-intel/bin" \
 CC="mpiicc -m64" CFLAGS=" -O3 -march=native -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include" \
 FC="mpiifort -m64" FCFLAGS="-O3 -I${MKLROOT}/include/intel64/lp64 -I${MKLROOT}/include" \
 GSL_CONFIG=/home/common/lib/gsl-1.16-intel/bin/gsl-config --with-gsl-prefix="/home/common/lib/gsl-1.16-intel/" \
 --with-libxc-prefix=$LIBXC_PREFIX \
 --with-blas="${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl" \
 --with-lapack="${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_sequential.a -Wl,--end-group -lpthread -lm -ldl" \
 --with-fft-lib=/home/common/lib/fftw-3.3.4-intel/lib/libfftw3.a
```
{{% /expand %}}

{{% expand "Magerit (CeSViMa)" %}}

Magerit is a cluster located in "Centro de Supercomputación y Visualización de Madrid". It consists on 245 nodes eServer BladeCenter PS702 with 2 Power7 processors each of 8 cores at 3'3 GHz (422,4 GFlops) and with 32 GB de RAM.
Here the "optimal" configuration script which uses GNU compilers with ESSL IBM library.

```bash
 module purge
 module load gcc/4.4
 module load openmpi/1.6.3
 
 export LIBSDIR='/sw/openmpi/'
 
 export SOFTWARE_DIR="$HOME/SOFTWARE/exec"
 
 export FC=mpif90
 export CC=mpicc
 export CXX=mpicxx
 
 export STDFLAGS="-O3 -mcpu=power7 -mtune=power7"
 
 export FCFLAGS=$STDFLAGS
 export CFLAGS=$STDFLAGS
 export CXXLAGS=$STDFLAGS
 
 export C_TYPE="gnu64-4.4-essl"
 
 export branch=`awk '{print $2}' ../.git/HEAD`
 export branch=`basename $branch`
 
 INSTALLDIR="$SOFTWARE_DIR/octopus/git-$C_TYPE/$branch"
 
 #GSL  1.16 --compiled with gcc 4.4.6
 export GSL_HOME="$LIBSDIR/GSL/1.16/"
 
 #LAPACK gcc 4.4
 export LAPACK_HOME="$LIBSDIR/LAPACK/3.4.2-gnu64-4.4"
 
 #PARMETIS gcc 4.7.2
 export PARMETIS_HOME="$LIBSDIR/METIS/PARMETIS-4.0.3/"
 
 #METIS gcc 4.7.2
 export METIS_HOME="$LIBSDIR/METIS/METIS-5.1.0/"
 
 #FFTW3.3.4 gcc 4.4
 export FFTW_HOME="$LIBSDIR/FFTW/3.3.4-gnu64-4.4"
 export LD_LIBRARY_PATH=$FFTW_HOME/lib:$LD_LIBRARY_PATH
 
 #SCALAPACK gcc 4.4.6
 export SCALAPACK_HOME="$LIBSDIR/SCALAPACK/2.0.2-gnu64-4.4/"
 export LD_LIBRARY_PATH=$SCALAPACK_HOME/lib:$LD_LIBRARY_PATH
 
 #LIBXC
 export LIBXC_HOME="$SOFTWARE_DIR/libxc/3.0.0-$C_TYPE/"
 
 #LIBVDWXC
 export LIBVDWXC_HOME="$SOFTWARE_DIR/libvdwxc/0.2.0-$C_TYPE/"
 export CPPFLAGS="-I$LIBVDWXC_HOME/include -DHAVE_LIBVDWXC"
 #export LIBS=-lvdwxcfort
 export LDFLAGS="-L$LIBVDWXC_HOME/lib -lvdwxcfort -Wl,-rpath=$LIBVDWXC_HOME/lib"
 export LD_LIBRARY_PATH=$LIBVDWXC_HOME/lib:$LD_LIBRARY_PATH
 
 export LD_LIBRARY_PATH=/opt/ibmcmp/lib64/:/opt/ibmcmp/xlsmp/2.1/lib64/:/opt/ibmcmp/xlf/13.1/lib64/:$LD_LIBRARY_PATH
 
 export LIBS_BLAS=" -L$SCALAPACK_HOME/lib -L$LAPACK_HOME/lib -lscalapack -llapack -lessl -lblacssmp -L$FFTW_HOME/lib -lfftw3 -lfftw3_threads -lfftw3_mpi / opt/gnu/gcc/4.4.6/lib64/libgfortran.a -L/opt/ibmcmp/xlf/13.1/lib64/ -lxlf90_r -lxlfmath -lxl -L/opt/ibmcmp/xlsmp/2.1/lib64/ -lxlomp_ser -Wl,--rpath /opt/ibmcmp/xlf/13.1/lib64/ /opt/ibmcmp/xlf/13.1/lib64/libxl.a -R/opt/ibmcmp/lib64/"
 
 make distclean
 make clean
 ../configure --prefix=$INSTALLDIR \
     --disable-gdlib --disable-gsltest \
     --with-libxc-prefix=$LIBXC_HOME \
     --with-fftw-prefix=$FFTW_HOME \
     --with-gsl-prefix=$GSL_HOME \
     --with-metis-prefix=$METIS_HOME\
     --with-parmetis-prefix=$PARMETIS_HOME \
     --enable-mpi \
     --enable-openmp
```
{{% /expand %}}


{{% expand "MPCDF systems (draco, cobra, raven)" %}}

Draco, Cobra and Raven are the HPC machines of the Max-Planck Compute and Data Facility in Garching. On these, Octopus can be built with the following script:

```bash
 #! /bin/bash -l
 # compilation script for octopus
 # needs to be executed in a subfolder (e.g. _build) of the root directory
 #
 # if you want to skip the configure step, call as ./build_octopus.sh noconfig
 # if you want to skip the configure step and the dependency checking,
 #   call as ./build_octopus.sh noconfig nodep
 
 if [[ ! -f ../configure.ac ]]; then
   echo "Error! Please execute this script in a subfolder of the root directory."
   exit 1
 fi
 
 compiler=intel
 #compiler=gnu
 cuda=yes
 
 echo Using $compiler compiler.
 [[ $cuda == yes ]] && echo Building with CUDA support.
 
 module purge
 if [[ $compiler == intel ]]; then
   # newest intel version
   module load intel/19.1.2 impi/2019.8 mkl/2020.2
   export CC=mpiicc
   export FC=mpiifort
   export CXX=mpiicpc
   export CFLAGS="-O3 -xCORE-AVX512 -qopt-zmm-usage=high -fma -ip"
   export FCFLAGS="$CFLAGS"
   export CXXFLAGS="$CFLAGS"
   export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -liomp5 -lpthread -lm -ldl"
 elif [[ $compiler == gnu ]]; then
   module load gcc/9 impi/2019.8 mkl/2020.2
   export CC=mpicc
   export FC=mpif90
   export CXX=mpicxx
   export CFLAGS="-O3 -march=skylake-avx512 -g"
   export FCFLAGS="$CFLAGS"
   export CXXFLAGS="$CFLAGS"
   export MKL="-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -lmkl_blacs_intelmpi_lp64 -lgomp -lpthread -lm -ldl"
 else
   echo "Compiler $compiler unknown."
   exit 1
 fi
 
 module load gsl hdf5-serial netcdf-serial libxc metis parmetis
 export LDFLAGS="-Xlinker -rpath=$MKL_HOME/lib/intel64:$GSL_HOME/lib:$NETCDF_HOME/lib:$ELPA_HOME/lib:$METIS_HOME/lib:$PARMETIS_HOME/lib"
 
 if [[ $cuda == yes ]]; then
   module load cuda/10.2
   CUDA_FLAGS="--enable-cuda --enable-nvtx --with-cuda-prefix=$CUDA_HOME"
   LDFLAGS="$LDFLAGS:$CUDA_HOME/lib64"
 else
   CUDA_FLAGS=""
 fi
 module list
 
 export INSTALLDIR=$PWD/installed
 
 if [[ "$1" != "noconfig" ]]; then
   pushd .. && autoreconf -i && popd
   ../configure $CUDA_FLAGS \
     FCFLAGS_FFTW="-I$MKLROOT/include/fftw" \
     FCCPP="cpp -ffreestanding" \
     --prefix=$INSTALLDIR \
     --enable-mpi --enable-openmp \
     --disable-gdlib \
     --with-gsl-prefix="$GSL_HOME" \
     --with-libxc-prefix="$LIBXC_HOME" \
     --with-blas="$MKL" \
     --with-lapack="$MKL" \
     --with-blacs="$MKL" \
     --with-scalapack="$MKL" \
     --with-netcdf-prefix="$NETCDF_HOME" \
     --with-metis-prefix="$METIS_HOME" \
     --with-parmetis-prefix="$PARMETIS_HOME" \
     || exit 1
 fi
 
 echo "\n\nBuilding octopus...\n"
 
 if [[ "$2" == "nodep" ]]; then
   make NODEP=1 -j20 && make NODEP=1 install || exit 1
 else
   make -j20 && make install || exit 1
 fi
 
 mkdir -p $INSTALLDIR/.build.doc/
 cp -f config.log $INSTALLDIR/.build.doc/
 cp -f $0 $INSTALLDIR/.build.doc/
 
 echo "... done"
```
{{% /expand %}}

{{< manual-foot prev="Tutorial:Running Octopus on Graphical Processing Units (GPUs)" next="Manual:Appendix:Porting Octopus and Platform Specific Instructions" >}}
---------------------------------------------
---
Title: "Installation"
Weight: 30
---

{{% children depth=3 %}}---
title: "Porting Octopus and Platform Specific Instructions"
#series: "Manual"
---

{{% notice note %}}
This page needs updating!
{{% /notice %}}

This page contains information about {{< octopus >}} portability, with specific information to compile and run octopus for many architectures. If you managed to compile octopus for a different system, please contribute.
Warning: this information is quite out of date and may no longer be valid.

### General information and tips about compilers

* {{< octopus >}} is developed in the {{< name "GNU" >}} environment and sometimes we depend on system features that are {{< name "GNU" >}} extensions without knowing it. These are bugs and we will try to fix them; please report any problem that you find.

* If you have problems with the C compiler, try to use {{< name "gcc" >}}. It normally works better than vendor compilers and it's available on most systems. However, be careful with locally installed versions of gcc: sometimes they don't work.

* The {{< name "Fortran" >}} {{< code "//" >}} concatenation operator is sometimes recognized as a {{< name "C++" >}}-style comment and the preprocessor gives erroneous results: sometimes it doesn't expand macros after it or simply eliminates what comes after. To solve this problem, use the preprocessor with the {{< name "-C" >}} (keep comments) and {{< name "-ansi" >}} or equivalent options (in {{< name "ANSI C" >}} {{< code "//" >}} is not a comment).

* If you are compiling in dual 32/64-bit architectures like {{< name "PowerPC" >}}, {{< name "UltraSparc" >}} or {{< name "AMD64" >}} systems here are some tips: 
** A 64-bit version of {{< octopus >}} is only needed if you are going to use more than 2-3 {{< name "Gb" >}} of physical {{< name "RAM" >}}.
** Some operating systems have 64 bits {{< name "kernels" >}} and 32 bits {{< name "userspace" >}} ({{< name "Solaris" >}}, {{< name "OS X" >}}); if you want a 64-bit {{< octopus >}} there, you have to compile all required libraries in 64-bit (normally a 64-bit {{< name "libc" >}} is available).
** Typically Linux distributions {{< emph "for" >}} {{< name "AMD64" >}} have a 64-bit {{< name "userspace" >}}, so you will get a 64-bit executable there.

### SSE2 support

* We have some SSE2 code written using compiler primitives that can give an important increase in perfomance. For this you need hardware support (AMD Opteron/Athlon 64/Sempron/Turion or Intel Pentium 4 or newer) and compiler support, supported compilers are GCC and pathcc. For gcc you need to put the correct {{< flag "-march" >}} flags (for example {{< code "<nowiki>-march=opteron</nowiki>" >}} or {{< code "<nowiki>-march=pentium4</nowiki>" >}}).

* Besides this, for {{< name "x86" >}} (32 bits) you have to link dynamically because we have to use a tweaked malloc function that doesn't work with static linking. For x86_64 (64 bits) this is not needed.

### Operating systems

#### Linux

The main development operating system for {{< octopus >}}.

#### Solaris

Octopus compiles correctly either with sun compilers or gcc/gfortran. By default {{< name "Solaris" >}} doesn't have {{< name "GNU coreutils" >}}, so some test won't run.

#### Tru 64

It works.

#### Mac OS X

It works. Don't try to compile static binaries, they are not supported by the OS.

#### Windows

Toy operating systems are not supported for the moment, sorry.

### Compilers

#### Intel Compiler for x86/x86_64

* status: ok
* Version 9 and version 7.1 Build 20040901Z are ok. Versions 8 and 8.1 can be problematic.
* Recommended flags: FCFLAGS="-u -zero -fpp1 -nbs -pc80 -pad -align -unroll -O3 -ip -tpp7 -xW"
* {{< name "Intel" >}} artificially blocks their compilers from using certain optimization in non-{{< name "Intel" >}} processors.
* With Octopus 3.2.0, use of the flags <tt>-check all -traceback</tt> with <tt>ifort</tt> 10.1.018 will cause an internal compiler segmentation fault while compiling <tt>src/grid/mesh_init.F90</tt>.

#### Intel Compiler for Itanium

* status: ok
* Version: 8.1.033 (older 8 releases and version 9 are reported to cause problems), version 10 works but it is much slower than 8.1.
* Recommended flags:
** FCFLAGS="-O3 -tpp2 -ip -IPF_fp_relaxed -ftz -align all -pad"
** CFLAGS="-O3 -tpp2 -ip -IPF_fp_relaxed -ftz"

#### [Open64](https://www.open64.net/)

This is an open source compiler based on the liberated code of SGI MIPSpro compiler. It is available for x86, x86_64 and Itanium architectures.

#### [Pathscale Fortran Compiler](https://www.pathscale.com/) 

* Versions tested: 2.2, 2.3 and 2.5
* Architecure: x86, AMD64
* Recommended flags:
```text
FCFLAGS="-Wall -O3 -march=auto -mcpu=auto -OPT:Ofast -fno-math-errno"
```
* Issues: 
** Everything works.
** It's necessary to compile blas/lapack with the same compiler.

#### NAG compiler

AMD64:
```text

FCFLAGS="-colour -kind=byte -mismatch_all -abi=64 -ieee=full -O4 -Ounroll=4"
```

#### [GNU C Compiler (gcc)](https://gcc.gnu.org/)

#### [GNU Fortran (gfortran)](https://gcc.gnu.org/fortran/)

* Status: ok.
* Version: gcc version 4.1.1 or newer. (4.1.0 {{< emph "does not" >}} work) For the parallel version you need at least gfortran 4.3.
* You may also need to compile blas, lapack and fftw3 using that specific gfortran version.
* Some recommended flags: <tt>-march=athlon64 -msse2 -mfpmath=sse -malign-double -funroll-loops -O3</tt>

#### [https://www.g95.org/ g95] 

* Status: works
* Tested architectures: x86/Linux, PowerPC/Darwin
* Version: version 4.0.3 (g95 0.91!) May 24 2007
* G95 doesn't recognize the linemarkers created by the preprocessor, so it's necessary to pass the -P flag to cpp.
* Flags:
```text

FC=g95
FCFLAGS="-O3 -funroll-loops -ffast-math"
FCCPP="cpp -ansi-P"
```

There may be problems with versions 0.92 or 0.93, depending on the underlying version of gcc. See [[G95]] for info on building version 0.94 with gcc 4.2.4.

#### Portland 6

Flags:

```text

FCFLAGS="-fast -mcmodel=medium -O4"
```


Known problems:

The following problem with the PGI compiler version 6.0 and MPICH version 1.2.6 on {{< name "x86_64" >}} has been reported:

The MPI detection during the {{< command "configure" >}} step does not work properly. This may lead to compilation failures on e. g. the file {{< file "par_vec.F90" >}}. This problem is considered a bug in either the PGI compiler or the MPICH implementation. Please apply the following change by hand after running {{< command "configure" >}}:

In the file {{< file "config.h" >}}, replace the line 

```text

/* -undef MPI_H */
```

by

```text

-define MPI_H 1
```

and remove the line

```text

-define MPI_MOD 1
```

#### Portland 7, 8, 9

Flags (tested on Cray XT4):

```text

FCFLAGS="-O2 -Munroll=c:1 -Mnoframe -Mlre -Mscalarsse -Mcache_align -Mflushz"
```

The configure script may fail in the part <tt>checking for Fortran libraries of mpif90</tt> for <tt>autoconf</tt> version 2.59 or earlier. The solution is to update <tt>autoconf</tt> to 2.60 or later, or manually set <tt>FCLIBS</tt> in the <tt>configure</tt> command line to remove a spurious apostrophe.

#### Portland 10

For Octopus 3.2.0, the file <tt>src/basic/lookup.F90</tt> is incorrectly optimized yielding many segmentation faults in the testsuite. With PGI 10.5 the optimization flag should be <tt>-O2</tt> or less; with PGI 10.8 the optimization flag should be <tt>-O1</tt> or less. Note that <tt>-fast</tt> and <tt>-fastsse</tt> are between <tt>-O2</tt> and <tt>-O3</tt>. For later versions of Octopus, a PGI pragma compels this file to be <tt>-O0</tt> regardless of what is specified in <tt>FCFLAGS</tt>, so you may safely set <tt>FCFLAGS</tt> to <tt>-fast</tt>.

####  Portland 11  

11.4 does not work and will crash with glibc memory corruption errors. 11.7 is fine.

#### Portland 12

12.5 and 12.6 cannot compile due to an internal compiler errors of this form:

```text
 PGF90-S-0000-Internal compiler error. sym_of_ast: unexpected ast    6034 (simul_box.F90: 1103)
```

12.4 and 12.9 are ok.

#### Absoft
Flags x86:
```text

FCFLAGS="-O3 -YEXT_NAMES=LCS -YEXT_SFX=_"
```

Flags amd64/em64t:
```text

FCFLAGS="-O3 -mcmodel=medium -m64 -cpu:host -YEXT_NAMES=LCS -YEXT_SFX=_"
```

#### Compaq compiler 

```text

FCFLAGS="-align dcommons -fast -tune host -arch host -noautomatic"
```

#### Xlf
* Status: works
-bmaxdata:0x80000000 -qmaxmem=-1 -qsuffix=f=f90 -Q -O5 -qstrict -qtune=auto -qarch=auto -qhot -qipa

* Because of the exotic mixture of {{< name "MAC OS" >}} and {{< name "BSD" >}}, this system is not very standard. Compiling {{< octopus >}} can be problematic.

* {{< name "OS X" >}} doesn't support static linking of binaries, so don't try.

#### SGI MIPS 

-O3 -INLINE -n32 -LANG:recursive=on

#### Sun Studio

You can download this compiler for free, it supports Linux and Solaris over x86, amd64 and sparc. A very fast compiler but quite buggy.

* Flags:
** CFLAGS="-fast -xprefetch -xvector=simd -D__SSE2__"
** FCFLAGS=$FLAGS

### MPI Implementations

#### OpenMPI

#### MPICH2

#### [SGI MPT](https://www.sgi.com/products/software/mpt/)

#### Intel MPI

#### [Sun HPC ClusterTools](https://www.sun.com/software/products/clustertools/) 

#### [MVAPICH](https://mvapich.cse.ohio-state.edu/)

### NetCDF

{{< octopus >}} uses the Fortran 90 interface of {{< name "netCDF" >}}, this means that it's likely that you will have to compile it using the same compiler you will use to compile {{< octopus >}}. You can get the sources and follow installation instructions from the [https://www.unidata.ucar.edu/software/netcdf/ {{< name "NetCDF" >}} site].

### BLAS and LAPACK

These are standard libraries that provide a series of common vector and matrix operations. {{< octopus >}} uses as much as possible this libraries. There are several version available depending on your hardware. Around 40% of {{< octopus >}} execution time is spend in BLAS level 1 routines, so getting a fast implementation for your hardware might be important. On the other hand, Lapack performance is not very important.

#### AMD ACML

This is the AMD Mathematical Library optimized to run in Athlon and Opteron processors. You can get a free copy from https://developer.amd.com/acml.jsp .

#### [ATLAS](https://math-atlas.sourceforge.net/) 

#### Compaq CXML

#### GOTO BLAS

Probably the fastest implementation of blas, source code is available and it can be compiled in many architectures.

#### Intel MKL

See https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor for MKL's advice on the proper way to link. Here is an example, in which {{< code "--with-lapack" >}} is left blank because it is included in {{< code "--with-blas" >}}.

```text
 MKL_DIR=/opt/intel/mkl/lib/lintel64
 --with-blas="-L$MKL_DIR -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread"
 --with-blacs="$MKL_DIR/libmkl_blacs_intelmpi_lp64.a" --with-scalapack="$MKL_DIR/libmkl_scalapack_lp64.a"
```

#### Netlib

The reference implementation of BLAS and Lapack. It is available in most linux distributions. You can get the source code from https://www.netlib.org/blas https://www.netlib.org/lapack .

{{< manual-foot prev="Manual:Specific architectures" next="Manual:Appendix:Reference Manual" >}}
---
title: "Building from scratch"
#series: "Manual"
---


Compiling {{< octopus >}} on some systems might be difficult, this is mainly because some libraries are required and they have to be all compiled with same Fortran compiler. This is a guide on how to compile {{< octopus >}} and the required libraries starting from scratch (a Fortran compiler is required) in a standard Unix system.

{{< notice note >}}
This page should be considered a last resort. On most systems you will be able to install many of these dependencies from a package manager or other pre-built binaries, which will save you a lot of effort.
{{< /notice >}}

###  Prerequisites  

You will need the following:

* A standard Unix/Linux system, where you want to install {{< octopus >}}.
* A basic knowledge on how to use the system (manage files, edit text files, extract tar archives, etc.).
* A directory where we will install the created binaries, this can be a specially created directory or even your {{< code "$HOME" >}}. From now, we will assume that the directory is {{< basedir >}}, whenever it appears you should replace it with the actual directory you are using. ''Note that the files that you will download and the directories with the sources, including the ones of {{< octopus >}}, don't need to be placed in this directory.''
* A working C compiler, {{< name "gcc" >}} is the ideal choice for {{< name "x86/x86_64" >}} systems. For other platforms you might want to choose a compiler that produces more optimized code. We will assume that the command to run the compiler is {{< cc >}}.
* A working Fortran 95 compiler, if you don't have one you can probably get [https://gcc.gnu.org/fortran/ {{< name "gfortran" >}}] or [https://www.g95.org/ {{< name "g95" >}}]. We will assume that the command to run the Fortran compiler is {{< f90 >}}.
* If you are compiling for a dual 32/64 bits architecture, it is a good idea to tell the compilers explicitly what type of binaries you want, many problems are caused by unadvertedly mixing 32 and 64 bits code. For {{< octopus >}} normally you will want a 64 bits binary ({{< code "-m64" >}} flag in most compilers). Since this flag should be passed always we recommend you to include it in the compiler command used in this guide (sometimes optimization flags are left out), so for example {{< cc >}} could be {{< command "gcc -m64" >}}.

###  Compiler flags  

Since probably you want {{< octopus >}} to run fast, you will probably would to set optimization flags for both the C and Fortran compilers on your system (at least -O3). In general {{< octopus >}} should run fine with aggressive optimization options. We will assume now that you have found a set of flags for {{< cc >}} and {{< f90 >}}, we will call them {{< cflags >}} and {{< fflags >}}. For example {{< cflags >}} could be {{< command "-O3 -march=native" >}}.

{{< notice note >}}
If {{< cc >}}, {{< f90 >}}, {{< cflags >}} or {{< fflags >}} contain spaces you should not enclose them in " ". We will put the " " explicitly when it is necessary.
{{< /notice >}}
###  Compilation of libraries  

First we will compile and install the libraries required by {{< octopus >}}.
####  BLAS  

The first library we will compile is {{< name "blas" >}}. We will use the generic reference implementation, this version is not highly optimized but it is free software and simple to install. So after you have a working version of {{< octopus >}} you may want to use a more optimized {{< name "blas" >}} implementation as [https://www.tacc.utexas.edu/resources/software/-blas {{< name "libgoto" >}}], [https://math-atlas.sourceforge.net {{< name "ATLAS" >}}], {{< name "MKL" >}}, {{< name "ACML" >}}, {{< name "ESSL" >}}, etc.

* Download the package from the {{< name "netlib" >}} site: https://www.netlib.org/blas/blas.tgz
* Extract the package and enter into the newly created {{< file "BLAS" >}} directory.
* Edit the {{< file "make.inc" >}} file and modify the definition of the fortran compiler that now should be:

```bash
 FORTRAN  = {{< f90 >}}
 OPTS     = {{< fcflags >}}
 DRVOPTS  = $(OPTS)
 NOOPT    =
 LOADER   = {{< f90 >}}
 LOADOPTS =
```
* Now type {{< command "make" >}} to compile, this will take a while.
* One compilation is finished, create the directory {{< basedir >}}{{< file "/lib" >}} 
```bash
 mkdir {{< basedir >}}/lib
```
and copy the file {{< file "blas_LINUX.a" >}} to {{< basedir >}}{{< file "/lib/libblas.a" >}}
```bash
 cp blas_LINUX.a {{< basedir >}}/lib/libblas.a
```
{{< notice note >}}
Note: the generated file librrary will be always called {{< file "blas_LINUX.a" >}} independently of the operating system you are using, this is just a name.
{{< /notice >}}
####  LAPACK  

We will use the open source reference implementation of Lapack.

* Get the Lapak package from netlib: https://www.netlib.org/lapack/lapack.tgz
* Extract the archive and enter the newly created lapack directory.
* Copy {{< file "make.inc.example" >}} to {{< file "make.inc" >}}:
```bash
 cp make.inc.example make.inc
```
* Edit {{< file "make.inc" >}} to indicate the compiler the will be used. The relevant part of the file should look like:
```bash
 FORTRAN  = {{< f90 >}}
 OPTS     = {{< fcflags >}}
 DRVOPTS  = $(OPTS)
 NOOPT    =
 LOADER   = {{< f90 >}}
 LOADOPTS =
```
* Now build the library with
```bash
 make lib
```
* Copy the newly created library {{< file "lapack_LINUX.a" >}} to {{< file "<basedir>/lib/liblapack.a" >}} .

####  GSL  

* Get the source of the latest version of GSL from ftp://ftp.gnu.org/gnu/gsl/ , currently it is GSL 1.16: ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz .
* Extract the archive and enter the newly created directory.
* Run the configure script with:
```bash
 ./configure CC="{{< cc >}}" --prefix={{< basedir >}} --disable-shared --enable-static
```
* Compile and install:
```bash
 make
 make install
```

####  FFTW 3  

* Get the sources for the latest version of [https://www.fftw.org/ FFTW], currently https://www.fftw.org/fftw-3.3.4.tar.gz .
* Extract the archive and enter the newly created directory.
* Configure:
```bash
 ./configure  --prefix={{< basedir >}} CC="{{< cc >}}" CFLAGS="{{< cflags >}}" F77="{{< f90 >}}" F77FLAGS="{{< fflags >}}"
```
* Compile and install:
```bash
 make
 make install
```

####  LibXC  

* Get {{< libxc >}}, currently https://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-{{< libxc-version >}}.tar.gz .
* Extract the archive and enter the newly created directory.
```bash
 tar -xvzf libxc-{{< libxc-version >}}.tar.gz
 cd libxc-{{< libxc-version >}}
```
* Configure:
```bash
 ./configure --prefix={{< basedir >}} CC="{{< cc >}}" CFLAGS="{{< cflags >}}" FC="{{< f90 >}}" FCFLAGS="{{< fcflags >}}"
```
* Compile and install:
```bash
 make
 make install
```

###  Compilation of Octopus  

After compiling the libraries, now we are ready to compile {{< octopus >}}. These are the steps you have to follow:

* Download the last version of {{< octopus >}}: https://www.tddft.org/programs/octopus/down.php?file={{< octopus-version >}}/octopus-{{< octopus-version >}}.tar.gz
* Extract the file and enter the newly created directory.
* Define the following environment variables to be used by the configure script (we assume that you are using {{< name "bash" >}}, if you are using another shell, the commands should be modified accordingly):
```bash
 export LIBS_BLAS={{< basedir >}}/lib/libblas.a
 export LIBS_LAPACK={{< basedir >}}/lib/liblapack.a
 export LIBS_FFT={{< basedir >}}/lib/libfftw3.a
```
* Now call the configure script:
```bash
 ./configure CC="{{< cc >}}" CFLAGS="{{< cflags >}}" FC="{{< f90 >}}" FCFLAGS="{{< fcflags >}}" --prefix={{< basedir >}} --with-gsl-prefix={{< basedir >}} --with-libxc-prefix={{< basedir >}}
```
* Compile and install
```bash
 make
 make install
```
* If everything went fine, you should have {{< octopus >}} installed in the {{< basedir >}} directory, this means that the executables are in {{< basedir >}}{{< file "/bin/" >}}. You may want to add this last directory to your path to run octopus commands directly, otherwise you will have to give the full path.
* To test the compilation is correct, you can run the testsuite of {{< octopus >}}, that compares the results obtained with the created version against reference values. To do it, run
```bash
 make check
```
All tests should be passed if the compilation is correct. If all of them fail, there is probably a problem with your executable, typically missing dynamic libraries. If
just some fail, there might be a more serious problem, we recommend you to look for help in the [[Mailing lists]].


{{< manual-foot prev="Manual:Updating to a new version" next="Tutorial:Running Octopus on Graphical Processing Units (GPUs)" >}}
---------------------------------------------
---
title: "Updating to a new version"
#series: "Manual"
---


This page lists the issues that may appear when updating to a new version of {{< octopus >}} and how to solve them.

#### Octopus 8

##### Variables

* The {{< variable "Species" >}} block can now take a {{< emph "set" >}} option to choose which pseudopotential set to use.
* The default for {{< variable "XCFunctional" >}} is now taken, when possible, from the pseudopotentials.

#### Octopus 7

##### Variables

* The UnitsInput and Units variables have been removed. Input file values are now always in atomic units, unless explicitly stated otherwise in the corresponding variable description. A few constants are available to help using eV/Angstrom units in the input file, such that one can write {{< code "1=Spacing = 0.25*angstrom" >}}.
* Now by default the code assumes XYZ files to be in Angstrom. This can be changed by using the {{< variable "UnitsXYZFiles" >}} variable.

#### Octopus 6

##### Installation

* A compiler supporting Fortran 2003 iso_c_binding is now required.
* The executable is always called {{< file "octopus" >}}, no longer {{< file "octopus_mpi" >}} when compiled in parallel.
* This version supports [[libxc 3.0.0]]. Although it is still possible to use [[libxc 2.0.0]] or any of the 2.x versions, updating to [[libxc 3.0.0]] is highly recommended.

##### Variables

* There is a new format for the {{< variable "Species" >}} block and the names of the species options have been revised.
* The variable OutputHow has been renamed to {{< variable "OutputFormat" >}}.
* The variable ParallelizationStrategy and ParallelizationGroupRanks have been replaced with {{< variable "ParDomains" >}}, {{< variable "ParStates" >}}, {{< variable "ParKPoints" >}}, and {{< variable "ParOther" >}}.

##### Utilities

The {{< file "oct-rotatory_strength" >}} utility has been removed. This utility is replaced by a new {{< variable "PropagationSpectrumType" >}} in {{< file "oct-propagation_spectrum" >}}. 

#### Octopus 5.0

##### Restart files

This version breaks the backwards compatibility of the restart files. Nevertheless, it is not too difficult to update the restart files generated with a previous version. Here is a summary of the changes and how to update the files:

* Casida restart data now goes to a directory {{< file "restart/casida" >}}, and a file {{< file "kernel" >}} or {{< file "kernel_triplet" >}}. One then needs to create the directory and rename the files.
* The format of the {{< file "occs" >}} file changed with the addition of an extra field for the imaginary part of the eigenvalues. To update the file one needs to add a column of zeros by hand after the column of eigenvalues.
* The line starting with {{< emph "fft_alpha" >}} was removed from the {{< file "mesh" >}} file.

Note that the mesh partition restart is also not compatible any more. In this case the partition needs to be recalculated, which is usually not a problem.

##### Variables

Some variables changed name or changed format. If an obsolete variable is found, {{< octopus >}} will stop and will tell you the variable that replaces it.

The datasets feature was removed. The input file of a calculation using datasets needs to be split into several independent input files.

#### Octopus 4.1

##### Installation

This version requires [[libxc 2.0.0]] or higher, so it might be necessary to update libxc before installing {{< octopus >}}. 

##### Restart files
Casida restart file is incompatible with that generated by version 4.0. The new format avoids problems with restarting with a different number of states.

#### Octopus 4.0

##### Installation

* Libxc is now an independent library. To compile {{< octopus >}} 4.0.0 you will have to compile [[libxc 1.1.0]] first. These are the short instructions to compile it (we assume that libxc will be installed in $DIR, this can be the same directory where {{< octopus >}} is going to be installed):
```text
 tar -xvzf libxc-1.1.0.tar.gz
 cd libxc-1.1.0 
 ./configure --prefix=$DIR
 make
 make install
```
Now, when configuring {{< octopus >}} pass the option --with-libxc-prefix=$DIR.

* The configure option for the location of netcdf and etsf_io are --with-netcdf-prefix and --with-etsf-io-prefix.

##### Variables

* TDEvolutionMethod is now called {{< variable "TDPropagator" >}}.

##### Utilities

* {{< file "oct-cross_section" >}} was renamed to {{< file "oct-propagation_spectrum" >}}.
* {{< file "oct-broad" >}} was renamed to {{< file "oct-casida_spectrum" >}}.
* The format for {{< file "oct-help" >}} has changed. Now {{< file "-s" >}} is used instead of {{< file "search" >}}, {{< file "-p" >}} instead of {{< file "show" >}}, and {{< file "-l" >}} instead of {{< file "list" >}}.

#### Octopus 3.0

##### Input variables

Some variables changed name or changed format. If an obsolete variable is found, {{< octopus >}} will stop and will tell you the variable that replaces it.

* {{< variable "Units" >}}, {{< variable "UnitsInput" >}} and {{< variable "UnitsOutput" >}} now take a named option as argument instead of a string. So

```text
 Units = "eVA"
```

should be replaced by
```text
 
 Units = eV_Angstrom
```

The code will stop if the old format is encountered.

* XFunctional and CFunctional were replaced by {{< variable "XCFunctional" >}}, for example

```text
 XFunctional = lda_x
 CFunctional = lda_c_vwn
```

must be replaced with

```text
 XCFunctional = lda_x + lda_c_vwn
```

* TDLasers was replaced by {{< variable "TDExternalFields" >}}.
* Some options for {{< variable "CalculationMode" >}} were renamed.

##### Output directories

Some directories in {{< octopus >}} output were renamed, restart files now are stored under {{< file "restart/" >}} instead of {{< file "tmp/" >}}. Files previously found under {{< file "status/" >}} are now located in {{< file "exec/" >}}.

##### Restart file format

{{< octopus >}} writes restart files in a binary format, this format has been updated and improved. As a result {{< octopus >}} is no longer able to restart automatically the calculation of a run performed with older versions. 

##### Recovering old restart files

If you really need to recover your old restart files, first you have to generate NetCDF restart files with your old version of {{< octopus >}}. Then you will have to rename the {{< file "tmp/" >}} directory to {{< file "restart/" >}} and remove the {{< file "restart_" >}} prefix from its subdirectories, for example {{< file "tmp/restart_gs/" >}} now should be called {{< file "restart/gs/" >}}.

Now you can run {{< octopus >}} 3.0 and it will find your restart information.


{{< manual-foot prev="Manual:Examples:Benzene" next="Manual:Building from scratch" >}}
---------------------------------------------
---
Title: "Tutorials"
weight: 30
description: "Entry page for the Octopus tutorials"
menu: "top"
---

Here you will find a collection of tutorials covering a wide range of topics, from the basics of performing calculations with {{<octopus>}} to more advanced features. Several series of linked tutorials are proposed to guide you on how to perform certain types of calculations. These tutorials do not by any means cover all the things that {{<octopus>}} can do for you, but hopefully they cover the most common options. You can find more information in the {{< manual "" "online Manual" >}}.

If you have never used {{<octopus>}} before, then you should start with the {{<octopus>}} basics series of tutorials.

A complete list of all available tutorials can be found here. The tutorials are also organized by categories that you can browse in case you are interested in some specific type of system, feature, or calculation mode.

## Main tutorial series

* {{< tutorial "Basics" "Octopus Basics" >}} - getting started with {{< octopus >}}.
* {{< tutorial "Response" "Optical Response" >}} - how to calculate several types of optical response with different methods.
* {{< tutorial "Model" "Model systems" >}} - working with model systems, like quantum dots or quantum wells.
* {{< tutorial "Periodic Systems" "Periodic systems" >}} - periodic boundary conditions, band structures, etc.

## Difficulty level

{{< tutorial-list "difficulties" >}}

## Theory level

{{< tutorial-list "theories" >}}

## Calculation mode

{{< tutorial-list "calculation_modes" >}}

## System type

{{< tutorial-list "system_types" >}}

## Species type

{{< tutorial-list "species_types" >}}

## Feature

{{< tutorial-list "features" >}}

## Utility

{{< tutorial-list "utilities" >}}
---
title: "Periodic systems"
#tags: ["Basic", "Ground State", "Unoccupied", "Bulk", "Pseudopotentials", "DFT", "Band Structure", "DOS"]
#series: "Tutorial"
tutorials: ["Octopus Basics", "Periodic Systems"]
difficulties: "basic"
theories: "DFT"
calculation_modes: ["Ground state", "Unoccupied"]
system_types: "bulk"
species_types: "Pseudopotentials"
features: ["Band structure", "DOS"]
description: "How to perform some basic calculation using bulk silicon as an example."
weight: 6
---


The extension of a ground-state calculation to a periodic system is quite straightforward in {{< octopus >}}. In this tutorial we will explain how to perform some basic calculation using bulk silicon as an example.

## Input
As always, we will start with a simple input file. In this case we will use a primitive cell of Si, composed of two atoms. 

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/periodic_systems/1.start/inp
{{< /code-block >}}


Lets see more in detail some of the input variables:

* {{< code-inline >}}{{< variable "PeriodicDimensions" >}} = 3{{< /code-inline >}}: this input variable must be set equal to the number of dimensions you want to consider as periodic. Since the system is 3D ({{< code-inline >}}{{< variable "Dimensions" >}} = 3{{< /code-inline >}} is the default), by setting this variable to 3 we impose periodic boundary conditions at all borders. This means that we have a fully periodic infinite crystal.

* {{< code-inline >}}{{< variable "LatticeVectors" >}}{{< /code-inline >}} and {{< code-inline >}}{{< variable "LatticeParameters" >}}{{< /code-inline >}}: these two blocks are used to define the primitive lattice vectors that determine the unit cell. {{< code-inline >}}{{< variable "LatticeVectors" >}}{{< /code-inline >}} defines the direction of the vectors, while {{< variable "LatticeParameters" >}} defines their length.

* {{< code-inline >}}{{< variable "ReducedCoordinates" >}}{{< /code-inline >}}: the position of the atoms inside the unit cell, in reduced coordinates.

* {{< code-inline >}}{{< variable "KPointsGrid" >}}{{< /code-inline >}}: this specifies the ''k''-point grid to be used in the calculation. Here we employ a 2x2x2 Monkhorst-Pack grid with four shifts. The first line of the block defines the number of ''k''-points along each axis in the Brillouin zone. Since we want the same number of points along each direction, we have defined the auxiliary variable {{< code "nk = 2" >}} This will be useful later on to study the convergence with respect to the number of ''k''-points. The other four lines define the shifts, one per line, expressed in reduced coordinates of the Brillouin zone. Alternatively, one can also define the reciprocal-space mesh by explicitly setting the position and weight of each ''k''-point using the {{< variable "KPoints" >}} or {{< variable "KPointsReduced" >}} variables.

* {{< code-inline >}}{{< variable "KPointsUseSymmetries" >}} = yes{{< /code-inline >}}: this variable controls if symmetries are used or not. When symmetries are used, the code shrinks the Brillouin zone to its irreducible portion and the effective number of ''k''-points is adjusted. 

* {{< code-inline >}}{{< variable "Output" >}} = dos{{< /code-inline >}}: we ask the code to output the density of states.

Here we have taken the value of the grid spacing to be 0.5 bohr. Although we will use this value throughout this tutorial, remember that in a real-life calculation the convergence with respect to the grid spacing must be performed for all quantities of interest.

Note that for periodic systems the default value for the {{< code-inline >}}{{< variable "BoxShape" >}}{{< /code-inline >}} variable is {{< code parallelepiped >}}, although in this case the name can be misleading, as the actual shape also depends on the lattice vectors. This is the only box shape currently available for periodic systems.

## Output
Now run {{< octopus >}} using the above input file. Here are some important things to note from the output.

{{< code-block >}}
******************************** Grid ********************************
Simulation Box:
  Type = parallelepiped
  Lengths [b] = (   3.627,   3.627,   3.627)
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 3 dimension(s).

  Lattice Vectors [b]
    0.000000    5.130000    5.130000
    5.130000    0.000000    5.130000
    5.130000    5.130000    0.000000
  Cell volume =           270.0114 [b^3]
  Reciprocal-Lattice Vectors [b^-1]
   -0.612396    0.612396    0.612396
    0.612396   -0.612396    0.612396
    0.612396    0.612396   -0.612396
Main mesh:
  Spacing [b] = ( 0.484, 0.484, 0.484)    volume/point [b^3] =      0.08000
  # inner mesh =       3375
  # total mesh =      10639
  Grid Cutoff [H] =    21.095389    Grid Cutoff [Ry] =    42.190778
**********************************************************************
{{< /code-block >}}

Here {{< octopus >}} outputs some information about the cell in real and reciprocal space.

{{< code-block >}}
***************************** Symmetries *****************************
Space group No.  227
International: Fd-3m
Schoenflies: Oh^7
  Index                Rotation matrix                      Fractional translations
    1 :     1   0   0     0   1   0     0   0   1      0.000000    0.000000    0.000000
    2 :     0  -1   1     0  -1   0     1  -1   0      0.000000    0.000000    0.000000
    3 :    -1   0   0    -1   0   1    -1   1   0      0.000000    0.000000    0.000000
    4 :     0   1  -1     1   0  -1     0   0  -1      0.000000    0.000000    0.000000
    5 :    -1   1   0    -1   0   0    -1   0   1      0.000000    0.000000    0.000000
    6 :     0   0   1     1   0   0     0   1   0      0.000000    0.000000    0.000000
    7 :     1  -1   0     0  -1   1     0  -1   0      0.000000    0.000000    0.000000
    8 :     0   0  -1     0   1  -1     1   0  -1      0.000000    0.000000    0.000000
    9 :     0  -1   0     1  -1   0     0  -1   1      0.000000    0.000000    0.000000
   10 :    -1   0   1    -1   1   0    -1   0   0      0.000000    0.000000    0.000000
   11 :     0   1   0     0   0   1     1   0   0      0.000000    0.000000    0.000000
   12 :     1   0  -1     0   0  -1     0   1  -1      0.000000    0.000000    0.000000
   13 :     0   0  -1     1   0  -1     0   1  -1      0.000000    0.000000    0.000000
   14 :    -1   1   0    -1   0   1    -1   0   0      0.000000    0.000000    0.000000
   15 :     1  -1   0     0  -1   0     0  -1   1      0.000000    0.000000    0.000000
   16 :     0   0   1     0   1   0     1   0   0      0.000000    0.000000    0.000000
   17 :     1   0  -1     0   1  -1     0   0  -1      0.000000    0.000000    0.000000
   18 :     0  -1   0     0  -1   1     1  -1   0      0.000000    0.000000    0.000000
   19 :     0   1   0     1   0   0     0   0   1      0.000000    0.000000    0.000000
   20 :    -1   0   1    -1   0   0    -1   1   0      0.000000    0.000000    0.000000
   21 :     0   1  -1     0   0  -1     1   0  -1      0.000000    0.000000    0.000000
   22 :     1   0   0     0   0   1     0   1   0      0.000000    0.000000    0.000000
   23 :    -1   0   0    -1   1   0    -1   0   1      0.000000    0.000000    0.000000
   24 :     0  -1   1     1  -1   0     0  -1   0      0.000000    0.000000    0.000000
Info: The system has    24 symmetries that can be used.
**********************************************************************
{{< /code-block >}}

This block tells us about the space-group and the symmetries found for the specified structure.

{{< code-block >}}
     2 k-points generated from parameters :
 ---------------------------------------------------
    n =    2    2    2
 
    s1  =  0.50  0.50  0.50
 
    s2  =  0.50  0.00  0.00
 
    s3  =  0.00  0.50  0.00
 
    s4  =  0.00  0.00  0.50
  
 index |    weight    |             coordinates              |
     1 |     0.250000 |    0.250000    0.000000    0.000000  |
     2 |     0.750000 |   -0.250000    0.250000    0.250000  |
{{< /code-block >}}

Next we get the list of the ''k''-points in reduced coordinates and their weights. Since symmetries are used, only two ''k''-points are generated. If we had not used symmetries, we would have 32 ''k''-points instead.

The rest of the output is much like its non-periodic counterpart. After a few iterations the code should converge:

{{< code-block >}}
*********************** SCF CYCLE ITER -    8 ************************
 etot  = -7.92843852E+00 abs_ev   =  4.39E-06 rel_ev   =  1.79E-05
 ediff =       -1.19E-06 abs_dens =  3.49E-05 rel_dens =  4.37E-06
Matrix vector products:     77
Converged eigenvectors:     10

#  State  KPoint  Eigenvalue [H]  Occupation    Error
      1       1       -0.257101    2.000000   (6.5E-07)
      1       2       -0.184777    2.000000   (9.3E-07)
      2       2       -0.079243    2.000000   (9.4E-07)
      2       1        0.010873    2.000000   (8.4E-07)
      3       2        0.023247    2.000000   (8.2E-07)
      4       2        0.073689    2.000000   (8.6E-07)
      3       1        0.128013    2.000000   (9.0E-07)
      4       1        0.128013    2.000000   (9.6E-07)
      5       2        0.209032    0.000000   (9.3E-07)
      5       1        0.232110    0.000000   (9.2E-07)

Density of states:

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
------------------------------------------------------%---------------
%---------%--------------%------------%-%------%------%-----------%--%
                                                      ^


Elapsed time for SCF step     8:          0.10
**********************************************************************


             Info: Writing states. 2018/08/06 at 17:15:24


        Info: Finished writing states. 2018/08/06 at 17:15:24

Info: SCF converged in    8 iterations
{{< /code-block >}}

As usual, the {{< file "static/info" >}} file contains the most relevant information concerning the calculation. Since we asked the code to output the density of states, we also have a few new files in the {{< file "static" >}} directory:

* {{< file "dos-XXXX.dat" >}}: the band-resolved density of states (DOS);

* {{< file "total-dos.dat" >}}: the total DOS (summed over all bands);

* {{< file "total-dos-efermi.dat" >}}: the Fermi Energy in a format compatible with {{< file "total-dos.dat" >}}.

Of course you can tune the output type and format in the same way you do in a finite-system calculation.

## Convergence in k-points

Similar to the convergence in spacing, a convergence must be performed for the sampling of the Brillouin zone. To do this one must try different numbers of ''k''-points along each axis in the Brillouin zone. This can easily be done by changing the value of the <tt>nk</tt> auxiliary variable in the previous input file. You can obviously do this by hand, but this is something that can also be done with a script. Here is such a script, which we will name {{< file "kpts.sh" >}}:

```bash
#include_file doc/tutorials/octopus_basics/periodic_systems/2.k-points/kpts.sh
```
 
After running the script ( {{< code "source kpts.sh" >}}), you should obtain a file called {{< file "kpts.log" >}} that should look like this:

{{< code-block >}}
#nk Total Energy
2 -7.92843852
4 -7.93452207
6 -7.93462084
8 -7.93462226
{{< /code-block >}}

As you can see, the total energy is converged to within 0.0001 hartree for {{< code "nk = 6" >}}.
For an extended range, e.g. from 2 to 12, the script should be

```bash
#!/bin/bash
echo "#nk Total Energy" > kpts.log
list="2 4 6 8 10 12"
for nkpt in $list
do
    sed -i  "s/nk = [0-9][0-9]*/nk = $nkpt/g" inp
    octopus >& out-$nkpt
    energy=`grep Total static/info  | head -2 | tail -1 | cut -d "=" -f 2`
    echo $nkpt $energy >> kpts.log
    rm -rf restart
done
```

Then, the {{< file "kpts.log" >}} file should look like this:

{{< code-block >}}
#nk Total Energy
2 -7.92843852
4 -7.93452207
6 -7.93462084
8 -7.93462226
10 -7.93462251
12 -7.93462217
{{< /code-block >}}

## Band-structure

We now proceed with the calculation of the band-structure of Si. In order to compute a band-structure, we must perform a non-self-consistent calculation, using the density of a previous ground-state calculation. So the first step is to obtain the initial ground-state. To do this, rerun the previous input file, but changing the number of k-points to {{< code "nk=6" >}}. Next, modify the input file such that it looks like this:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/periodic_systems/4.band-structure/inp
{{< /code-block >}}

Here are the things we changed:

* {{< code-inline >}}{{< variable "CalculationMode" >}} = unocc{{< /code-inline >}}: we are now performing a non-self-consistent calculation, so we use the {{< code unoccupied >}} calculation mode;

* {{< code-inline >}}{{< variable "ExtraStates" >}} = 10{{< /code-inline >}}: this is the number of unoccupied bands to calculate;

* {{< code-inline >}}{{< variable "ExtraStatesToConverge" >}} = 5{{< /code-inline >}}: the highest unoccupied states are very hard to converge, so we use this variable to specify how many unoccupied states are considered for the stopping criterion of the non-self-consistent run.

* {{< code-inline >}}{{< variable "KPointsPath" >}}{{< /code-inline >}}: this block is used to specify that we want to calculate the band structure along a certain path in the Brillouin zone. This replaces the {{< variable "KPointsGrid" >}} block. The first row describes how many ''k''-points will be used to sample each segment. The next rows are the coordinates of the ''k''-points from which each segment starts and stops. In this particular example, we chose the following path: L-Gamma, Gamma-X, X-Gamma using a sampling of 10-10-15 ''k''-points.

* {{< code-inline >}}{{< variable "KPointsUseSymmetries" >}} = no{{< /code-inline >}}: we have turned off the use of symmetries.

After running {{< octopus >}} with this input file, you should obtain a file named {{< file "bandstructure" >}} inside the {{< file "static" >}} directory. This is how the first few lines of the file should look like:

{{< code-block >}}
# coord. kx ky kz (red. coord.), bands:    14 [H]
     0.00000000    0.50000000    0.00000000    0.00000000    -0.19984750    -0.10463050     0.11059562     0.11059567     0.21320524     0.27708534     0.27708537     0.43517458     0.54466879     0.54466894     0.57089044     0.57496796     0.57496801     0.63339690
     0.02640129    0.45000000   -0.00000000   -0.00000000    -0.20506492    -0.09711271     0.11125581     0.11125586     0.21395854     0.27788925     0.27788927     0.43640709     0.51950978     0.51950989     0.55510464     0.60282729     0.60282738     0.64782337
     0.05280258    0.40000000   -0.00000000   -0.00000000    -0.21725636    -0.07804331     0.11323409     0.11323414     0.21620689     0.28007561     0.28007563     0.43967737     0.48696242     0.48696251     0.52612903     0.64322288     0.64322300     0.66832900
...
{{< /code-block >}}

#include_eps doc/tutorials/octopus_basics/periodic_systems/4.band-structure/Si_bandstructure.eps caption="Band structure of bulk silicon. The zero of energy has been shifted to the maximum of the occupied bands."

The first column is the coordinate of the ''k''-point along the path. The second, third, and fourth columns are the reduced coordinates of the ''k''-point. The following columns are the eigenvalues for the different bands. In this case there are 14 bands (4 occupied and 10 unoccupied).

On the right you can see the plot of the band structure. This plot shows the occupied bands (purple) and the first 5 unoccupied bands (green). If you are using gnuplot, you can obtain a similar plot with the following command:

```text
plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
```

Note that when using the {{< variable "KPointsPath" >}} input variable, {{< octopus >}} will run in a special mode, and the restart information of the previous ground-state calculation will not be altered in any way. The code informs us about this just before starting the unoccupied states iterations:

```text
Info: The code will run in band structure mode.
     No restart information will be printed.
```

{{< tutorial-footer >}}










---
title: "Recipe"
tags: ["Basic", "Recipe"]
#series: "Tutorial"
tutorials: ["Octopus Basics"]
description: "...finally, learn how to prepare a tasty octopus dish"
difficulties: "basic"
weight: 8
---


If the input file sets

{{< code-block >}}
 {{< variable "CalculationMode" >}} = recipe
{{< /code-block >}}

then the code will print you a tasty recipe, randomly selected from those available in the package. For example:

{{< code-block >}}
 ************************** Calculation Mode **************************
 Input: [CalculationMode = recipe]
 **********************************************************************
 
 Info: Octopus initialization completed.
 Info: Starting calculation mode.
 
 
 OCTOPUS a la GALLEGA:
 
 Ingredients: For 4 persons
   - 2 kg. of octopus
   - 2 kg. of potatoes
   - 100 grs. of paprika
   - 100 grs. thick salt
   - olive oil
 
 Preparation: Wash the octopus in cold water. Place a copper (Cu) pan with water
 on the fire, and when it starts boiling submerge the octopus, and remove
 it again. Repeat this procedure 3 times, always waiting for the water to
 start boiling. Then, boil the octopus for 20 minutes, remove the pan from
 the fire and let it rest for 5 minutes. Remove the octopus from the water,
 and cut it in thin slices with some scissors. It is served on wood plates,
 spiced in the following order: first the salt, then the paprika and the
 olive oil. It is served with potatoes.
 
 
 DISCLAIMER: The authors do not guarantee that the implementation of this
 recipe leads to an edible dish, for it is clearly "system-dependent".
{{< /code-block >}}

There are currently recipes available in English, Spanish, Italian, and Basque. Which one you get is controlled by the LANG environment variable, so if you run "LANG=es octopus", you may get:

{{< code-block >}}
 PULPO A FEIRA:
 
 Ingredientes: Para 4 personas
   - 2 kg. de pulpo
   - 2 kg. de papas / patatas
   - 100 grs. de ají / chile / pimentón picante
   - 100 grs. de sal gorda
   - aceite
 
 Preparación: Se lava el pulpo en agua fría, se pone una olla de cobre con
 agua al fuego y cuando rompa a hervir se coge el pulpo, se mete y se saca
 del agua tres veces dejando que en cada intervalo vuelva a hervir el agua.
 Se deja cocer el pulpo durante unos 20 minutos retirándolo del fuego y
 dejándolo reposar durante 5 minutos. A continuación, se quita del agua y
 se corta en trozos finos con unas tijeras. Para servirlo se pone en unos
 platos de madera condimentándolo por este orden: sal, pimentón, aceite y
 se añaden unos cachelos (Patatas).
{{< /code-block >}}

Contributions in other languages are welcome!

{{< tutorial-footer >}}


---
title: "Total energy convergence"
tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
#series: "Tutorial"
tutorials: ["Octopus Basics"]
difficulties: "basic"
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"
description: "Make sure that the results are converged"
weight: 3
---


In this tutorial we will show how to converge a quantity of interest (the total energy) with respect to the parameters that define the real-space grid. For this we will use two examples, the Nitrogen atom from the {{< tutorial "Basics:Basic_input_options" "Basic input options" >}} tutorial, and a methane molecule.

## Nitrogen atom: finding a good spacing

The key parameter of a real-space calculation is the spacing between the points of the mesh. The default option in {{< octopus >}} is to use a Cartesian grid. This is a regular grid where the spacing along each Cartesian direction is the same. The first step in any calculation should then be making sure that this spacing is good enough for our purposes. This should be done through a convergence study, very similar to the ones performed in plane-wave calculations.

The needed spacing essentially depends on the pseudopotentials that are being used. The idea is to repeat a series of ground-state calculations, with identical input files except for the grid spacing. There are many different ways of doing it, the simplest one being to change the input file by hand and run {{< octopus >}} each time. Another option is to use a little bash script like the one bellow. 

As a first example, we will use the Nitrogen atom from the previous tutorial, so make sure you have the corresponding {{< file "inp" >}} and {{< file "N.xyz" >}} files in your working directory. Next, create a file called {{< code spacing.sh >}} with the following content:


```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/spacing.sh
```

What this script does is quite simple: it runs {{< octopus >}} for a list of spacings (0.26, 0.24, ..., 0.14) and after each calculation it collects the values of the total energy and of the eigenvalues, and prints them to a file called {{< file "spacing.log" >}}. In order to change the value of the spacing between each calculation, it uses the feature that one can override variables in the {{< file "inp" >}} file by {{< manual "Advanced_ways_of_running_Octopus#passing-arguments-from-environment-variables" "defining them as environment variables" >}} in the shell. Note that we unset the variable OCT_Spacing at the end of the script, to avoid it to remaining defined in case you run the script directly in the shell. 

Now, to run the script type 
```bash
source spacing.sh
```
(**NOT** {{< command-line "sh spacing.sh" >}} or {{< command-line "bash spacing.sh" >}}!). 

Note that for this to work, the {{< file "octopus" >}} executable must be found in the shell path. If that is not the case, you can replace
```bash
  octopus >& out-$Spacing
```
with the actual path to the executable
```bash
  /path/to/octopus >& out-$Spacing
```

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/1.nitrogen_atom_spacing/NitrogenSpacing.eps caption="Convergence with spacing of N"

Once the script finishes running, the {{< file "spacing.log" >}} file should look something like this: 

```text
 #Sp    Energy        s_eigen   p_eigen
 0.26 -256.59213179 -19.851724 -6.757584
 0.24 -260.25614052 -18.814241 -7.081600
 0.22 -262.57033689 -18.192766 -7.311321
 0.20 -262.91105590 -18.098579 -7.357661
 0.18 -262.19593061 -18.286207 -7.290490
 0.16 -261.78224630 -18.392571 -7.247292
 0.14 -261.78536939 -18.389733 -7.248998
```

You can also plot the results from the file in your favorite plotting program, e.g. {{< command "gnuplot" >}}.

The results, for this particular example, are shown in the figure. In this figure we are actually plotting the error with respect to the most accurate results (smallest spacing). That means that we are plotting the '''difference''' between the values for a given spacing and the values for a spacing of 0.14 Å. So, in reading it, note that the most accurate results are at the left (smallest spacing). A rather good spacing for this nitrogen pseudopotential seems to be 0.18 Å. However, as we are usually not interested in total energies, but in energy differences, probably a larger one may also be used without compromising the results.

## Methane molecule

We will now move on to a slightly more complex system, the methane molecule CH<sub>4</sub>, and add a convergence study with respect to the box size.

### Input

As usual, the first thing to do is create an input file for this system. From our basic chemistry class we know that methane has a tetrahedral structure. The only other thing required to define the geometry is the bond length between the carbon and the hydrogen atoms. If we put the carbon atom at the origin, the hydrogen atoms have the coordinates given in the following input file:

{{< code-block>}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.22*angstrom
 
 CH = 1.2*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
{{< /code-block >}}

Here we define a variable CH that represents the bond length between the carbon and the hydrogen atoms, which simplifies the writing of the coordinates. We start with a bond length of 1.2 Å but that's not so important at the moment, since we can optimize it later (for more information see the [Geometry optimization](../Geometry optimization) tutorial). (You should not use 5 Å or so, but something slightly bigger than 1 Å is fine.)

Some notes concerning the input file:
* We do not tell {{< octopus >}} explicitly which {{< variable "BoxShape" >}} to use. The default is a union of spheres centered around each atom. This turns out to be the most economical choice in almost all cases.
* We also do not specify the %{{< variable "Species" >}} block. In this way, {{< octopus >}} will use default pseudopotentials for both Carbon and Hydrogen. This should be OK in many cases, but, as a rule of thumb, you should do careful testing before using any pseudopotential for serious calculations.

If you use the given input file you should find the following values in the resulting {{< file "static/info" >}} file.

{{< code-block >}}
Eigenvalues [eV]
 #st  Spin   Eigenvalue      Occupation
   1   --   -15.988934       2.000000
   2   --    -9.064039       2.000000
   3   --    -9.064039       2.000000
   4   --    -9.064039       2.000000

Energy [eV]:
      Total       =      -219.01537542
      Free        =      -219.01537542
      -----------
      Ion-ion     =       236.08498119
      Eigenvalues =       -86.36210292
      Hartree     =       393.44961913
      Int[n*v_xc] =      -105.38115372
      Exchange    =       -69.66918783
      Correlation =       -11.00060043
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         0.00000000
      -TS         =        -0.00000000
      Kinetic     =       163.96741296
      External    =      -931.84569058
      Non-local   =       -49.74424172
{{< /code-block >}}


Now the question is whether these values are converged or not. This will depend on two things: the {{< variable "Spacing" >}}, as seen above, and the {{< variable "Radius" >}}. Just like for the Nitrogen atom, the only way to answer this question is to try other values for these variables.

###  Convergence with the spacing  

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/Spacing_CH4.eps caption="Convergence with spacing of methane"

As before, we will keep all entries in the input file fixed except for the spacing that we will make smaller by 0.02 Å all the way down to 0.1 Å. So you have to run {{< octopus >}} several times. You can use a similar script as the one from the Nitrogen atom example:

```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/3.methane_spacing/spacing.sh
```

For this example we changed the list of spacings to try and, to keep things simple, we are only writing the total energy to the {{< file "spacing.log" >}} file.

Now run the script to get the following results:

{{< code-block >}}
#Sp    Energy
0.22 -219.01537351
0.20 -218.58291295
0.18 -218.20448388
0.16 -218.14599180
0.14 -218.14037189
0.12 -218.13666585
0.10 -218.13059765
{{< /code-block >}}


If you give these numbers to ''gnuplot'' (or other software to plot) you will get a curve like the one shown on the right.

As you can see from this picture, the total energy is converged to within 0.1 eV for a spacing of 0.18 Å. So we will use this spacing for the next calculations.

### Convergence with the radius

#include_eps doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/Radius_CH4.eps caption="Convergence with radius of methane"

Now we will see how the total energy changes with the {{< variable "Radius" >}} of the box. We will change the radius in steps of 0.5 Å. You can change the input file by hand and run Octopus each time, or again use a small script:

```bash
#include_file doc/tutorials/octopus_basics/total_energy_convergence/4.methane_radius/radius.sh
```


Before running the script, make sure that the spacing is set to 0.18 Å in the input file. You should then get a file {{< file "radius.log" >}} that looks like this:

```text
#Rad   Energy
2.5 -217.99419600
3.0 -218.16990473
3.5 -218.20448384
4.0 -218.21155936
4.5 -218.21258666
5.0 -218.21322714
```

On the right you can see how this looks like when plotted. Note that in this case, the most accurate results are on the right (larger radius).

If we again ask for a convergence up to 0.1 eV we should use a radius of 3.5 Å. How does this compare to the size of the molecule? Can you explain why the energy increases (becomes less negative) when one decreases the size of the box?

{{< tutorial-footer >}}

---
title: "Octopus Basics"
weight: 10
description: "Basic introduction to run Octopus"
---

This is a series of linked tutorials that will teach you the basics of using Octopus. This is where you should start if you have never used Octopus before. It assumes that you have a basic understanding of solid state physics and Density Functional Theory, that you have Octopus installed on your system, and that you know how to run commands on a shell. 


{{< tutorial-series "Octopus Basics" >}}---
title: "Time-dependent propagation"
tags: ["Basic", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Laser"]
#series: "Tutorial"
tutorials: ["Octopus Basics"]
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
difficulties: "basic"
features: "Laser"
weight: 7
description: "Time evolution of the density of the methane molecule."
---


OK, it is about time to do TDDFT. Let us perform a calculation of the time evolution of the density of the methane molecule. 

## Ground-state

The reason to do first a ground-state DFT calculation is that this will be the initial state for the real-time calculation.

Run the input file below to obtain a proper ground-state (stored in the {{< code restart >}} directory) as from the 
{{< tutorial "basics:Total Energy Convergence" "total energy convergence tutorial" >}}. (Now we are using a more accurate bond length than before.)

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/time-dependent_run/1.gs/inp
{{< /code-block >}}

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
{{< /code-block >}}

## Time-dependent

#### Input

Now that we have a starting point for our time-dependent calculation, we modify the input file so that it looks like this:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/time-dependent_run/2.td/inp
{{< /code-block >}}

{{< code-block >}}
 {{< variable "CalculationMode" >}} = td
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
  
 Tf  = 0.1/eV
 dt = 0.002/eV
 
 {{< variable "TDPropagator" >}} = aetrs
 {{< variable "TDMaxSteps" >}} = Tf/dt
 {{< variable "TDTimeStep" >}} = dt
{{< /code-block >}}

Here we have switched on a new run mode, {{< code-inline >}}{{< variable "CalculationMode" >}} = td{{< /code-inline >}} instead of {{< code gs >}}, to do a time-dependent run. We have also added three new input variables:

* {{< variable "TDPropagator" >}}: which algorithm will be used to approximate the evolution operator. {{< octopus >}} has a large number of possible propagators that you can use (see Ref.[^footnote-1] for an overview).


* {{< variable "TDMaxSteps" >}}: the number of time-propagation steps that will be performed.
* {{< variable "TDTimeStep" >}}: the length of each time step.

Note that, for convenience, we have previously defined a couple of variables, {{< code Tf >}} and {{< code dt >}}. We have made use of one of the possible propagators, {{< code aetrs>}}. The manual explains about the possible options; in practice this choice is usually good except for very long propagation where the {{< code etrs >}} propagator can be more stable.

#### Output

Now run {{< octopus >}}. This should take a few seconds, depending on the speed of your machine. The most relevant chunks of the standard output are

{{< code-block >}}
Input: [IonsConstantVelocity = no]
Input: [Thermostat = none]
Input: [MoveIons = no]
Input: [TDIonicTimeScale = 1.000]
Input: [TDTimeStep = 0.2000E-02 hbar/eV]
Input: [TDPropagationTime = 0.1000 hbar/eV]
Input: [TDMaxSteps = 50]
Input: [TDDynamics = ehrenfest]
Input: [TDScissor = 0.000]
Input: [TDPropagator = aetrs]
Input: [TDExponentialMethod = taylor]
{{< /code-block >}}
...
{{< code-block >}}
********************* Time-Dependent Simulation **********************
  Iter           Time        Energy   SC Steps    Elapsed Time
**********************************************************************

      1       0.002000   -218.785065         1         0.118
      2       0.004000   -218.785065         1         0.121
      3       0.006000   -218.785065         1         0.118
{{< /code-block >}}
...
{{< code-block >}}
     49       0.098000   -218.785065         1         0.123
     50       0.100000   -218.785065         1         0.130
{{< /code-block >}}

It is worthwhile to comment on a few things:

* We have just performed the time-evolution of the system, departing from the ground-state, under the influence of no external perturbation. As a consequence, the electronic system does not evolve. The total energy does not change (this you may already see in the output file, the third column of numbers), nor should any other observable change. However, this kind of run is useful to check that the parameters that define the time evolution are correct.
* As the evolution is performed, the code probes some observables and prints them out. These are placed in some files under the directory {{< file "td.general" >}}, which should show up in the working directory. In this case, only two files show up, the {{< file "td.general/energy" >}}, and the {{< file "td.general/multipoles" >}} files. The {{< file "td.general/multipoles" >}} file contains a large number of columns of data. Each line corresponds to one of the time steps of the evolution (except for the first three lines, that start with a {{< code "-" >}} symbol, which give the meaning of the numbers contained in each column, and their units). A brief overview of the information contained in this file follows:
** The first column is just the iteration number.
** The second column is the time.
** The third column is the dipole moment of the electronic system, along the $x$-direction: $\\langle \\Phi (t) \\vert \\hat{x} \\vert \\Phi(t)\\rangle = 
\\int\\!\\!{\\rm d}^3r\\; x\\,n(r)$. Next are the $y$- and $z$-components of the dipole moment.
* The file {{< file "td.general/energy" >}} contains the different components of the total energy.
* It is possible to restart a time-dependent run. Try that now. Just increase the value of {{< variable "TDMaxSteps" >}} and rerun {{< octopus >}}. If, however, you want to start the evolution from scratch, you should set the variable {{< variable "FromScratch" >}} to {{< code yes >}}.

## The time step

A key parameter is, of course, the time step, {{< variable "TDTimeStep" >}}. Before making long calculations, it is worthwhile spending some time choosing the largest time-step possible, to reduce the number of steps needed. This time-step depends crucially on the system under consideration, the spacing, on the applied perturbation, and on the algorithm chosen to approximate the evolution operator. 

In this example, try to change the time-step and to rerun the time-evolution. Make sure you are using {{< variable "TDMaxSteps" >}} of at least 100, as with shorter runs an instability might not appear yet. You will see that for time-steps larger than {{< code "0.0024" >}} the propagation gets unstable and the total energy of the system is no longer conserved. Very often it diverges rapidly or even becomes NaN.

Also, there is another input variable that we did not set explicitly, relying on its default value, {{< variable "TDExponentialMethod" >}}. Since most propagators rely on algorithms to calculate the action of the exponential of the Hamiltonian, one can specify which algorithm can be used for this purpose.

You may want to learn about the possible options that may be taken by {{< variable "TDExponentialMethod" >}}, and {{< variable "TDPropagator" >}} -- take a look at the manual. You can now try some exercises:

* Fixing the propagator algorithm (for example, to the default value), investigate how the several exponentiation methods work (Taylor, Chebyshev, and Lanczos). This means finding out what maximum time-step one can use without compromising the proper evolution.

* And fixing now the {{< variable "TDExponentialMethod" >}}, one can now play around with the various propagators.

## Laser fields

Now we will add a time-dependent external perturbation (a laser field) to the molecular Hamiltonian. For brevity, we will omit the beginning of the file, as this is left unchanged. The relevant part of the input file, with the modifications and additions is:

{{< code-block >}}
#include_input_snippet doc/tutorials/octopus_basics/time-dependent_run/3.laser/inp laser
{{< /code-block >}}

#include_eps doc/tutorials/octopus_basics/time-dependent_run/3.laser/Tutorial_TD_Laser.eps caption="Laser field used in the tutorial" 


The most important variables here are the {{< variable "TDExternalFields" >}} block and the associated {{< variable "TDFunctions" >}} block. You should carefully read the manual page dedicated to these variables: the particular laser pulse that we have employed is the one whose envelope function is a cosine.

Now you are ready to set up a run with a laser field. Be careful to set a total time of propagation able to accommodate the laser shot, or even more if you want to see what happens afterwards. You may also want to consult the meaning of the variable {{< variable "TDOutput" >}}.

A couple of important comments:
* You may supply several laser pulses: simply add more lines to the {{< variable "TDExternalFields" >}} block and, if needed, to the {{< variable "TDFunctions" >}} block.
* We have added the {{< code laser >}} option to {{< variable "TDOutput" >}} variable, so that the laser field is printed in the file {{< file "td.general/laser" >}}. This is done immediately at the beginning of the run, so you can check that the laser is correct without waiting. 
* You can have an idea of the response to the field by looking at the dipole moment in the {{< file "td.general/multipoles" >}}. What physical observable can be calculated from the response?
* When an external field is present, one may expect that unbound states may be excited, leading to ionization. In the calculations, this is reflected by density reaching the borders of the box. In these cases, the proper way to proceed is to setup absorbing boundaries (variables {{< variable "AbsorbingBoundaries" >}}, {{< variable "ABWidth" >}} and {{< variable "ABHeight" >}}).

[^footnote-1]: {{< article title="Propagators for the time-dependent Kohn–Sham equations" authors="A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Chem. Phys." volume="121" pages="3425-3433" year="2004" doi="10.1063/1.1774980" >}}


{{< tutorial-footer >}}


---
title: "Centering a geometry"
tags: ["Basic", "Molecule", "oct-center-geom"]
#series: "Tutorial"
tutorials: ["Octopus Basics"]
difficulties: "basic"
system_types: "Molecule"
utilities: "oct-center-geom"
description: "Translate the center of mass to the origin."
weight: 5
---


Before running an {{< octopus >}} calculation of a molecule, it is always a good idea to run the {{< manual "External utilities:oct-center-geom"  "oct-center-geom" >}} utility, which will translate the center of mass to the origin, and align the molecule, by default so its main axis is along the ''x''-axis. Doing this is often helpful for visualization purposes, and making clear the symmetry of the system, and also it will help to construct the simulation box efficiently in the code. The current implementation in {{< octopus >}} constructs a parallelepiped containing the simulation box and the origin, and it will be much larger than necessary if the system is not centered and aligned. For periodic systems, these considerations are not relevant (at least in the periodic directions).

## Input
For this example we need two files.

#### {{< file "inp" >}}
We need only a very simple input file, specifying the coordinates file.

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/centering_a_geometry/inp
{{< /code-block >}}

#### {{< file "tAB.xyz" >}}

We will use this coordinates file, for the molecule [https://en.wikipedia.org/wiki/Azobenzene ''trans''-azobenzene], which has the interesting property of being able to switch between ''trans'' and ''cis'' isomers by absorption of light.

{{< code-block >}}
#include_file doc/tutorials/octopus_basics/centering_a_geometry/tAB.xyz
{{< /code-block >}}

## Centering the geometry

{{< manual "External utilities:oct-center-geom" "oct-center-geom" >}} is a serial utility, and it will be found in the {{< code "bin" >}} directory after installation of {{< octopus >}}.

When you run the utility, you should obtain a new coordinates file {{< code "adjusted.xyz" >}} for use in your calculations with {{< octopus >}}. The output is not especially interesting, except for perhaps the symmetries and moment of inertia tensor, which come at the end of the calculation:
 
{{< code-block >}}
***************************** Symmetries *****************************
Symmetry elements : (sigma)
Symmetry group    : Cs
**********************************************************************

Using main axis        1.000000       0.000000       0.000000
Input: [AxisType = inertia]
Center of mass [b] =        5.410288       2.130154      -1.005363
Moment of inertia tensor [amu*b^2]
      6542962.320615     -4064048.884840     -3486905.049193
     -4064048.884840      8371776.892023     -2789525.719087
     -3486905.049193     -2789525.719087     10730451.796958
Isotropic average      8548397.003199
Eigenvalues:            1199576.606657          11623018.898148          12822595.504791
Found primary   axis        0.707107       0.565686       0.424264
Found secondary axis       -0.624695       0.780869      -0.000000
{{< /code-block >}}

You can now visualize the original and new coordinates files with {{< code "xcrysden" >}} or your favorite visualization program (''e.g.'' Jmol, Avogadro, VMD, Vesta, etc.), and see what transformations the utility performed. You can see more options about how to align the system at {{< variable "AxisType" >}} and {{< variable "MainAxis" >}}.

{{< tutorial-footer >}}

---
title: "Getting started"
#tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
tutorials: ["Octopus Basics"]
difficulties: "basic"
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"
weight: 1
description: "Learn how to run the code"
---


The objective of this tutorial is to give a basic idea of how {{< octopus >}} works.

### Generating the input file

With a text editor, create a text file called {{< file "inp" >}} containing the following text:

<!--
{{< code-block >}}
{{< variable "CalculationMode" >}} = gs

%{{< variable "Coordinates" >}}
 'H' | 0 | 0 | 0
%
{{< /code-block >}}
-->
{{< code-block >}}
#include_input testsuite/tutorials/01-octopus_basics-getting_started.01-H_atom.inp
{{< /code-block >}}


This is the simplest example of an {{< octopus >}} input file:

* {{< code-inline >}}{{< variable "CalculationMode" >}} = gs {{< /code-inline >}}: This variable defines the run mode -- please consult the manual for the full list of the possible run modes. In this case we set it to {{< code gs >}}, which instructs the code to start a ground-state calculation.

* {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}: The entry is not just the definition of a variable, but rather of a full set of them -- a "block" of variables. The beginning of a block is marked by the {{< code "%identifier" >}} line, and ended by a {{< code "%" >}} line. In this case the identifier is {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}, where we list the atoms or species in our calculation and its coordinates, one per line. In this case, we put a single hydrogen atom in the center of our simulation box. 

The reason this input file can be so simple is that {{< octopus >}} comes with default values for the simulation parameters, and a set of default pseudopotentials for several elements (for properly converged calculations you might need to adjust these parameters, though).

To get a general idea of the format of the {{< octopus >}} input file, go and read the page about the {{< manual "Basics/Input file" "Input file" >}} in the manual.

The documentation for each input variable can be found in the {{< versioned-link "Variables/" "variable reference" >}} online, and can also be accessed via the {{< manual "External utilities/oct-help" "oct-help" >}} utility.

### Running Octopus

Once you have written your input file, run the {{< command "octopus" >}} command (using {{< command "mpirun" >}} and perhaps a job script if you are using the parallel version). If everything goes correctly, you should see several lines of output in the terminal (if you don't, there must be a problem with your installation). As this is probably the first time you run {{< octopus >}}, we will examine the most important parts of the output.

{{% notice note %}}
Be aware that the precise values you find in the output might differ from the ones in the tutorial text. This can be due to updates in the code, or also changes in the compilation and run configuration.
{{% /notice %}}

* First there is an octopus drawn in ASCII art, the copyright notice and some information about the octopus version you are using and the system where you are running:

{{< code-block >}}
      <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
                                ___
                             .-'   `'.
                            /         \
                            |         ;
                            |         |           ___.--,
                   _.._     |0) ~ (0) |    _.---'`__.-( (_.
            __.--'`_.. '.__.\    '--. \_.-' ,.--'`     `""`
           ( ,.--'`   ',__ /./;   ;, '.__.'`    __
           _`) )  .---.__.' / |   |\   \__..--""  """--.,_
          `---' .'.''-._.-'`_./  /\ '.  \ _.-~~~````~~~-._`-.__.'
                | |  .' _.-' |  |  \  \  '.               `~---`
                 \ \/ .'     \  \   '. '-._)
                  \/ /        \  \    `=.__`~-.
             jgs  / /\         `) )    / / `"".`\
            , _.-'.'\ \        / /    ( (     / /
             `--~`   ) )    .-'.'      '.'.  | (
                    (/`    ( (`          ) )  '-;
                     `      '-;         (-'

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2, or (at your option)
    any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA

    <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

                           Running octopus

Version                : 8.1
Commit                 : 9f1e56678d2698471df70e38f7dd4e37fd3d9044
Build time             : Thu Jul 19 14:51:54 CEST 2018
Configuration options  : max-dim=4 sse2 avx
Optional libraries     : arpack berkeleygw etsf_io gdlib netcdf pspio sparskit nlopt
Architecture           : x86_64
C compiler             : gcc
C compiler flags       : -g -Wall -O2 -march=native
Fortran compiler       : gfortran (GCC version 6.4.0)
Fortran compiler flags : -g -Wall -ffree-line-length-none -O2 -march=native

             The octopus is swimming in fenugreek (Linux)


            Calculation started on 2018/07/19 at 15:04:27

{{< /code-block >}}

Note that it also gives you the revision number, the compiler, and the compiler flags used. You should always include this information when submitting a bug report!

* The type of calculation it was asked to perform:
{{< code-block >}}
************************** Calculation Mode **************************
Input: [CalculationMode = gs]
**********************************************************************
{{< /code-block >}}

* The species and pseudopotentials it is using:
{{< code-block >}}
****************************** Species *******************************
  Species 'H'
    type             : pseudopotential
    file             : '/opt/local/share/octopus/pseudopotentials/PSF/H.psf'
    file format      : PSF
    valence charge   : 1.0
    atomic number    :   1
    form on file     : semilocal
    orbital origin   : calculated
    lmax             : 0
    llocal           : 0
    projectors per l : 1
    total projectors : 0
    application form : local
    orbitals         : 16
    bound orbitals   :  1
 
**********************************************************************
{{< /code-block >}}


* After some other output, {{< octopus >}} prints information about the grid: as we didn't say anything in the input file, {{< octopus >}} used the parameters recommended for this pseupopotential:
{{< code-block >}}
******************************** Grid ********************************
Simulation Box:
  Type = minimum
  Species =     H     Radius =   7.559 b
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [b] = ( 0.435, 0.435, 0.435)    volume/point [b^3] =  0.08210
  # inner mesh =      22119
  # total mesh =      37759
 Grid Cutoff [H] =    26.123439    Grid Cutoff [Ry] =    52.246878
**********************************************************************
{{< /code-block >}}


* The level of theory and, in the case of (TD)DFT, the approximation to the exchange-correlation term:
{{< code-block >}}
**************************** Theory Level ****************************
Input: [TheoryLevel = dft]

Exchange-correlation:
  Exchange
    Slater exchange (LDA)
    [1] P. A. M. Dirac, Math. Proc. Cambridge Philos. Soc. 26, 376 (1930)
    [2] F. Bloch, Z. Phys. 57, 545 (1929)
  Correlation
    Perdew & Zunger (Modified) (LDA)
    [1] J. P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981), modified to improve the matching between the low- and high-rs

Input: [SICCorrection = sic_none]
**********************************************************************
{{< /code-block >}}


* At this point, {{< octopus >}} tries to read the wave-functions from a previous calculation. As there are none, it will give a warning.
{{< code-block >}}
** Warning:
**   Could not find 'restart/gs' directory for restart.
**   No restart information will be read.

** Warning:
**   Unable to read wavefunctions.
**   Starting from scratch!
{{< /code-block >}}

* Now {{< octopus >}} commences the calculation. To get a reasonable starting point for the DFT calculation, the initial wavefunctions are calculated as a {{< manual "Calculations:Ground_State#LCAO" "Linear Combination of Atomic Orbitals" >}} (LCAO).
{{< code-block >}}
Info: Performing initial LCAO calculation with      1 orbitals.
Info: Getting Hamiltonian matrix elements. 
ETA: .......1......2.......3......4......5.......6......7.......8......9......0
 
Eigenvalues [H]
 #st  Spin   Eigenvalue      Occupation
   1   --    -0.233314       1.000000
Info: Ground-state restart information will be written
{{< /code-block >}}

* After the LCAO, the real DFT calculation starts. For each self-consistency step some information is printed. When SCF {{< manual "Calculations:Ground_State#Convergence" "converges" >}}, the calculation is done.
{{< code-block >}}
*********************** SCF CYCLE ITER -    1 ************************
 etot  = -4.48042396E-01 abs_ev   =  1.09E-03 rel_ev   =  4.66E-03
 ediff =       -2.76E-03 abs_dens =  8.90E-03 rel_dens =  8.90E-03
Matrix vector products:     27
Converged eigenvectors:      0

#  State  Eigenvalue [H]  Occupation    Error
      1       -0.234405    1.000000   (6.0E-05)

Elapsed time for SCF step     1:          0.12
**********************************************************************
{{< /code-block >}}
...
{{< code-block >}}
*********************** SCF CYCLE ITER -    5 ************************
 etot  = -4.46377047E-01 abs_ev   =  3.50E-06 rel_ev   =  1.50E-05
 ediff =       -4.24E-06 abs_dens =  4.33E-06 rel_dens =  4.33E-06
Matrix vector products:      7
Converged eigenvectors:      1

#  State  Eigenvalue [H]  Occupation    Error
      1       -0.233013    1.000000   (9.6E-07)

Elapsed time for SCF step     5:          0.04
**********************************************************************


             Info: Writing states. 2018/07/19 at 15:13:58


        Info: Finished writing states. 2018/07/19 at 15:13:58

Info: SCF converged in    5 iterations

Info: Finished writing information to 'restart/gs'.

             Calculation ended on 2018/07/19 at 15:13:58

                          Walltime:  01. 10s

Octopus emitted 2 warnings.
{{< /code-block >}}

Just running the command {{< command "octopus" >}} will write the output directly to the terminal. To have a saved copy of the output, it is generally advisable to redirect the output into a file, and to capture the standard error stream as well, which can be done like this: {{< command "octopus &> log" >}}. That would create a file called {{< file "log" >}} containing all output including warnings and errors in their context.

###  Analyzing the results  

After finishing the calculation you will find a series of files in the directory you ran:

{{< code-block >}}
% ls
  exec inp restart static
{{< /code-block >}}

For the moment we will ignore the '''exec'''  and  '''restart''' directories and focus on the {{< file "static/info" >}} file, which contains the detailed results of the ground-state calculation. If you open that file, first you will see some parameters of the calculations (that we already got from the output) and then the calculated energies and eigenvalues in Hartrees:

{{< code-block >}}
 Eigenvalues [H]
  #st  Spin   Eigenvalue      Occupation
    1   --    -0.233013       1.000000
 
 Energy [H]:
       Total       =        -0.44637705
       Free        =        -0.44637705
       -----------
       Ion-ion     =         0.00000000
       Eigenvalues =        -0.23301327
       Hartree     =         0.28415332
       Int[n*v_xc] =        -0.30429841
       Exchange    =        -0.19375604
       Correlation =        -0.03975282
       vanderWaals =         0.00000000
       Delta XC    =         0.00000000
       Entropy     =         1.38629436
       -TS         =        -0.00000000
       Kinetic     =         0.41780616
       External    =        -0.91483022
       Non-local   =         0.00000000
{{< /code-block >}}


Since by default {{< octopus >}} does a spin-unpolarized density-functional-theory calculation with the local-density approximation, our results differ from the exact total energy of 0.5 H. Our exchange-correlation functional can be set by the variable {{< variable "XCFunctional" >}}, using the set provided by the {{% libxc %}} library.

### Extra

If you want to improve the LDA results, you can try to repeat the calculation with spin-polarization:

{{< code-block >}}
 {{< variable "SpinComponents" >}} = spin_polarized
{{< /code-block >}}

And if you want to obtain the exact Schödinger equation result (something possible only for very simple systems like this one) you have to remove the self-interaction error (a problem of the LDA). Since we only have one electron the simplest way to do it for this case is to use independent electrons:

{{< code-block >}}
 {{< variable "TheoryLevel" >}} = independent_particles
{{< /code-block >}}

A more general way would be to include self-interaction correction.

{{< tutorial-footer >}}

---
title: "Visualization"
tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Visualization"]
tutorials: ["Octopus Basics"]
difficulties: "basic"
theories: ["DFT"]
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Visualization"
weight: 4
description: "Example: Benzene"
---


In this tutorial we will explain how to visualize outputs from {{< octopus >}}. Several different file formats are available in {{< octopus >}} that are suitable for a variety of visualization software. We will not cover all of them here. See {{< manual "Visualization" "Visualization" >}} for more information.


## Benzene molecule

As an example, we will use the benzene molecule. At this point, the input should be quite familiar to you:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/visualization/inp
{{< /code-block >}}

Coordinates are in this case given in the {{< file "benzene.xyz" >}} file. This file should look like:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/visualization/benzene.xyz
{{< /code-block >}}

Here we have used two new input variables:

* {{< variable "Output" >}}: this tells {{< octopus >}} what quantities to output.
* {{< variable "OutputFormat" >}}: specifies the format of the output files.

In this case we ask {{< octopus >}} to output the wavefunctions ({{< code wfs >}}), the density, the electron localization function ({{< code elf >}}) and the Kohn-Sham potential. We ask {{< octopus >}} to generate this output in several formats, that are required for the different visualization tools that we mention below. If you want to save some disk space you can just keep the option that corresponds to the program you will use. Take a look at the documentation on variables {{< variable "Output" >}} and {{< variable "OutputFormat" >}} for the full list of possible quantities to output and formats to visualize them in.

If you now run {{< octopus >}} with this input file, you should get the following files in the {{< file "static" >}} directory:

```text
 %ls static
 convergence      forces      vh.z=0       wf-st0001.cube     wf-st0003.dx       wf-st0005.xsf      wf-st0007.y=0,z=0  wf-st0009.z=0      wf-st0012.cube     wf-st0014.dx
 density.cube     info        vks.cube     wf-st0001.dx       wf-st0003.xsf      wf-st0005.y=0,z=0  wf-st0007.z=0      wf-st0010.cube     wf-st0012.dx       wf-st0014.xsf
 density.dx       v0.cube     vks.dx       wf-st0001.xsf      wf-st0003.y=0,z=0  wf-st0005.z=0      wf-st0008.cube     wf-st0010.dx       wf-st0012.xsf      wf-st0014.y=0,z=0
 density.xsf      v0.dx       vks.xsf      wf-st0001.y=0,z=0  wf-st0003.z=0      wf-st0006.cube     wf-st0008.dx       wf-st0010.xsf      wf-st0012.y=0,z=0  wf-st0014.z=0
 density.y=0,z=0  v0.xsf      vks.y=0,z=0  wf-st0001.z=0      wf-st0004.cube     wf-st0006.dx       wf-st0008.xsf      wf-st0010.y=0,z=0  wf-st0012.z=0      wf-st0015.cube
 density.z=0      v0.y=0,z=0  vks.z=0      wf-st0002.cube     wf-st0004.dx       wf-st0006.xsf      wf-st0008.y=0,z=0  wf-st0010.z=0      wf-st0013.cube     wf-st0015.dx
 elf_rs.cube      v0.z=0      vxc.cube     wf-st0002.dx       wf-st0004.xsf      wf-st0006.y=0,z=0  wf-st0008.z=0      wf-st0011.cube     wf-st0013.dx       wf-st0015.xsf
 elf_rs.dx        vh.cube     vxc.dx       wf-st0002.xsf      wf-st0004.y=0,z=0  wf-st0006.z=0      wf-st0009.cube     wf-st0011.dx       wf-st0013.xsf      wf-st0015.y=0,z=0
 elf_rs.xsf       vh.dx       vxc.xsf      wf-st0002.y=0,z=0  wf-st0004.z=0      wf-st0007.cube     wf-st0009.dx       wf-st0011.xsf      wf-st0013.y=0,z=0  wf-st0015.z=0
 elf_rs.y=0,z=0   vh.xsf      vxc.y=0,z=0  wf-st0002.z=0      wf-st0005.cube     wf-st0007.dx       wf-st0009.xsf      wf-st0011.y=0,z=0  wf-st0013.z=0
 elf_rs.z=0       vh.y=0,z=0  vxc.z=0      wf-st0003.cube     wf-st0005.dx       wf-st0007.xsf      wf-st0009.y=0,z=0  wf-st0011.z=0      wf-st0014.cube
```

## Visualization

### Gnuplot


#include_eps doc/tutorials/octopus_basics/visualization/benzene_2D.eps caption="Electronic density of benzene along the z=0 plane. Generated in Gnuplot with the 'pm3d' option." width="700px"

Visualization of the fields in 3-dimensions is complicated. In many cases it is more useful to plot part of the data, for example the values of a scalar field along a line of plane that can be displayed as a 1D or 2D plot. 

The {{< code gnuplot >}} program is very useful for making these types of plots (and many more) and can directly read Octopus output. First type {{< code gnuplot >}} to enter the program's command line. Now to plot a the density along a line (the X axis in this case) type

```bash
  plot "static/density.y=0,z=0" w l
```

To plot a 2D field, for example the density in the plane z=0, you can use the command

```bash
  splot "static/density.z=0" w l
```

If you want to get a color scale for the values of the field use "{{< code "w pm3d" >}}" instead of "{{< code "w l" >}}".
{{% expand "Gnuplot script" %}}
The above plot has been produced with the script:
{{< code-block >}}
#include_file doc/tutorials/octopus_basics/visualization/plot.gp
{{< /code-block >}}
{{% /expand %}}

### VisIt

The program [VisIt](https://visit.llnl.gov/) is an open-source data visualizer developed by Lawrence Livermore National Laborartory. It supports multiple platforms including Linux, OS X and Windows, most likely you can download a precompiled version for your platform [here](https://wci.llnl.gov/simulation/computer-codes/visit/executables). 

For this tutorial you need to have Visit installed and running. VisIt can read several formats, the most convenient for Octopus data is the <tt>cube</tt> format: {{< variable "OutputFormat" >}} = <tt>cube</tt>. This is a common format that can be read by several visualization and analysis programs, it has the advantage that includes both the requested scalar field and the geometry in the same file. Visit can also read other formats that Octopus can generate like <tt>xyz</tt>, <tt>XSF</tt> (partial support), and <tt>VTK</tt>. 

The following video contains the instruction on how to open a <tt>cube</tt> file in VisIt. Note that the video has subtitles in case you don't have audio in your current location, you can enable them by clicking the "CC" button.

{{< youtube id="c9E4iObMQ-M" scale="70%">}}
<youtube width="800" height="600">https://www.youtube.com/watch?v=c9E4iObMQ-M</youtube>


### XCrySDen

Atomic coordinates (finite or periodic), forces, and functions on a grid can be plotted with the free program [XCrysDen](https://www.xcrysden.org/). Its XSF format also can be read by [V_Sim](https://inac.cea.fr/sp2m/L_Sim/V_Sim/index.en.html) and [VESTA](https://jp-minerals.org/vesta/en/download.html). Beware, these all probably assume that your output is in Angstrom units (according to the [specification](https://www.xcrysden.org/doc/XSF.html)), so use {{< variable "UnitsOutput" >}} = eV_Angstrom, or your data will be misinterpreted by the visualization software.



### UCSF Chimera

{{% figure src="/images/Tutorial_chimera_1.png" width="500px" caption="Ground-state density plotted with UCSF Chimera for the benzene molecule, colored according to the value of the electrostatic potential" %}}


Those {{< file ".dx" >}} files can be also be opened by [UCSF Chimera](http:s//www.cgl.ucsf.edu/chimera/) .

Once having UCSF Chimera installed and opened, open the benzene.xyz file from '''File > Open'''. 

To open the {{< file ".dx" >}} file, go to '''Tools > Volume Data > Volume Viewer'''. A new window opens. There go to '''File > Open''' and select {{< file "static/density.dx" >}}.

To "color" the density go to '''Tools > Surface Color'''. A new dialog opens and there choose '''Browse''' to open the vh.dx file. Then click '''Color''' and you should have an image similar to the below one.

{{< tutorial-footer >}}

---
title: "Basic input options"
#tags: ["Basic", "Ground State", "Molecule", "Pseudopotentials", "DFT", "Total Energy"]
weight: 2
tutorials: ["Octopus Basics"]
difficulties: "basic"
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Total energy"
description: "Obtain the ground state of the nitrogen atom."
---


Now we will move to a more complicated (and realistic) input file. We will obtain the ground state of the nitrogen atom. We will introduce several basic input variables and will give a more detailed description of the output for this example.

## The input files

This sample input file lets us obtain the ground state of the nitrogen atom, within the LDA approximation, in a closed-shell (unpolarized) configuration (as explained below, you need an auxiliary {{< file ".xyz" >}} input). Note that this is not the correct ground state of the nitrogen atom! However, it will permit us to describe some of the most important input variables:

{{< code-block >}}
#include_input doc/tutorials/octopus_basics/basic_input_options/inp
{{< /code-block >}}

We have introduced here several new variables:

* {{< code-inline >}}{{< variable "UnitsOutput" >}} = eV_Angstrom{{< /code-inline >}}: Two different unit systems may be used for output: the usual atomic units (which is the default, and the ones used internally in the code); and the system in which the Ångström is substituted for the atomic unit of length, and the electronvolt is substituted for the atomic unit of energy. You can find a more detailed description of units in {{< octopus >}} in the {{< manual "Basics/Units" "Units" >}} page of the manual.

* The following entry in the input file is not a variable that {{< octopus >}} will read directly, but rather illustrates the possibility of writing "user-defined" values and expressions to simplify the input file. In this case, we define the nitrogen mass ({{< code-inline >}}Nitrogen_mass = 14.0{{< /code-inline >}}) (note that in this case, as an exception, the value is expected to be in the so-called "atomic mass units", rather than in "atomic units"). This definition may be used elsewhere in the input file.

* The {{< variable "Species" >}} block should contain the list of species that are present in the system to be studied. In this case we have only one species: nitrogen. The first field is a string that defines the name of the species, "N" in this case. The second field defines the type of species, in this case {{< code-inline >}}species_pseudo{{< /code-inline >}}. Then a list of parameters follows. The parameters are specified by a first field with the parameter name and the field that follows with the value of the parameter. Some parameters are specific to a certain species while others are accepted by all species. In our example {{< code-inline >}}set{{< /code-inline >}} instructs {{< octopus >}} to use a pseudopotential for nitrogen from the {{< code-inline >}}standard{{< /code-inline >}} set. This happens to be a Troullier-Martins pseudopotential defined in the {{< file "N.psf" >}} file found in the directory {{< file "share/octopus/pseudopotentials/PSF" >}}. Then come maximum {{< code-inline >}}lmax{{< /code-inline >}} - component of the pseudopotential to consider in the calculation, and the {{< code-inline >}}lloc{{< /code-inline >}} - component to consider as local. Generally, you want to set the maximum ''l'' to the highest available in the pseudopotential and the local ''l'' equal to the maximum ''l''. Finally, the mass of the species can also be modified from the default values by setting {{< code-inline >}}mass{{< /code-inline >}} parameter.

* {{< code-inline >}}{{< variable "XYZCoordinates" >}} = 'N.xyz'{{< /code-inline >}}: The geometry of the molecule (in this case, a single atom in the grid origin) is described in this case in a file with the well known {{< code-inline >}}XYZ{{< /code-inline >}} format. The file for this outrageously simple case is given by:

{{< code-block >}}
 1
 This is a comment line
 N 0 0 0
 {{< /code-block >}}

* {{< code-inline >}}{{< variable "ExtraStates" >}} = 1{{< /code-inline >}}: By default, {{< octopus >}} performs spin-unpolarized calculations (restricted closed-shell, in Hartree-Fock terminology). It then places two electrons in each orbital. The number of orbitals, or Kohn-Sham states, is then calculated by counting the number of valence electrons present in the system, and dividing by two. In this case, since we have five valence electrons, the code would use three orbitals. However, we know beforehand that the HOMO orbital has a three-fold degeneracy, and as a consequence we need to put each one of the three _p_ electrons in a different orbital. We therefore need one more orbital, which we get with this line in the input file.

* {{< code-inline >}}%{{< variable "Occupations" >}}{{< /code-inline >}} block: Generally, the occupations of the Kohn-Sham orbitals are automatically decided by the code, filling the lowest-energy orbitals. However, if we have degeneracies in the LUMO as in this case, the user may want to accommodate the electrons in a certain predefined way. In this example, the obvious way to fill the orbitals of the nitrogen atom is to put two electrons in the first and deepest orbital (the _s_ orbital), and then one electron on each of the second, third and fourth orbitals (the _p_ orbitals, which should be degenerate).

* {{< code-inline >}}{{< variable "BoxShape" >}} = sphere{{< /code-inline >}}: This is the choice of the shape of the simulation box, which in this case is set to be a sphere (other possible choices are {{< code-inline >}}minimum{{< /code-inline >}}, {{< code-inline >}}cylinder{{< /code-inline >}}, or {{< code-inline >}}parallelepiped{{< /code-inline >}}).

* {{< code-inline >}}{{< variable "Radius" >}} = 5.0*angstrom{{< /code-inline >}}: The radius of the sphere that defines the simulation box.

* {{< code-inline >}}{{< variable "Spacing" >}} = 0.18*angstrom{{< /code-inline >}}: As you should know, {{< octopus >}} works in a real-space regular cubic mesh. This variable defines the spacing between points, a key numerical parameter, in some ways equivalent to the energy cutoff in plane-wave calculations.

## Output

Once you have constructed the input file and created the {{< file "N.xyz" >}} file, you may unleash {{< octopus >}} on it. Lets now go over some of the sections of the output.

#### Species

{{< code-block >}}
****************************** Species *******************************
  Species 'N'
    type             : pseudopotential
    file             : '/opt/share/octopus/pseudopotentials/PSF/N.psf'
    file format      : PSF
    valence charge   : 5.0
    atomic number    :   7
    form on file     : semilocal
    orbital origin   : calculated
    lmax             : 1
    llocal           : 0
    projectors per l : 1
    total projectors : 1
    application form : kleinman-bylander
    orbitals         : 16
    bound orbitals   : 16

**********************************************************************
{{< /code-block >}}

Here the code searches for the needed pseudopotential files, and informs the user about its success or failure. In this case, only the {{< file "N.psf" >}} file is required. Once that file has been processed, some information about it is written to the output. One of the most important pieces of information to be found here is the valence charge, which tells us how many electrons from this species will be considered in the calculation.

#### Grid

{{< code-block >}}
******************************** Grid ********************************
Simulation Box:
  Type = sphere
  Radius  [A] =   5.000
  Octopus will run in 3 dimension(s).
  Octopus will treat the system as periodic in 0 dimension(s).
Main mesh:
  Spacing [A] = ( 0.180, 0.180, 0.180)    volume/point [A^3] =      0.00583
  # inner mesh =      89727
  # total mesh =     127183
  Grid Cutoff [eV] =  1160.586810    Grid Cutoff [Ry] =    85.301565
**********************************************************************
{{< /code-block >}}

This step is about the construction of the mesh. As requested in the input file, a sphere of radius 5 Å is used, which contains a cubic regular real-space grid with spacing 0.18 Å. This implies 89727 points ({{< code-inline >}}inner mesh =  89727{{< /code-inline >}}). For the sake of comparison with plane-wave-based codes, this is more or less equivalent to a plane-wave calculation that imposes a density cutoff of 1160.595 eV = 42.6 Hartree (except that in this case there is no artificial periodic repetition of the system).

#### Mixing

{{< code-block >}}
Input: [MixField = potential] (what to mix during SCF cycles)
Input: [MixingScheme = broyden]
{{< /code-block >}}

During the self-consistent procedure one has to use a {{< manual "Calculations/Ground_State#Mixing" "mixing scheme" >}} to help convergence. One can mix either the density or the potential, and there are several mixing schemes available.

#### Eigensolver

{{< code-block >}}
**************************** Eigensolver *****************************
Input: [Eigensolver = cg]
Input: [Preconditioner = pre_filter]
Input: [PreconditionerFilterFactor = 0.5000]
Input: [SubspaceDiagonalization = standard]
**********************************************************************
{{< /code-block >}}

Here we see that the {{< manual "Calculations:Ground_State#Eigensolver" "eigensolver" >}} used will be simple conjugate gradients (cg), and a preconditioner is used to speed up its convergence.

#### LCAO
After some output you should see something like:

{{< code-block >}}
Info: Performing initial LCAO calculation with      4 orbitals.
Info: Getting Hamiltonian matrix elements.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

Eigenvalues [eV]
 #st  Spin   Eigenvalue      Occupation
   1   --   -17.398866       2.000000
   2   --    -6.414013       1.000000
   3   --    -6.414013       1.000000
   4   --    -6.414013       1.000000
{{< /code-block >}}

This is the first step of a ground-state calculation: obtaining a reasonably good starting density and Kohn-Sham orbitals to feed in the self-consistent (SCF) procedure. For this purpose, {{< octopus >}} performs an initial calculation restricted to the basis set of atomic orbitals ( {{< manual "Calculations:Ground_State#LCAO" "Linear Combination of Atomic Orbitals" >}}, LCAO). The resulting eigenvalues of this calculation are written to standard output.

#### Wavefunction kind

{{< code-block >}}
 Info: SCF using real wavefunctions.
{{< /code-block >}}

Very often one can work with real wave-functions. This is particularly helpful as calculations with real wave-functions are much faster than with complex ones. However, if a magnetic field is present, if the system is periodic, or if spin-orbit coupling is present, complex wave-functions are mandatory. But don't worry: the program is able to figure out by itself what to use.

#### SCF

{{< code-block >}}
*********************** SCF CYCLE ITER -    1 ************************
 etot  = -2.54485949E+02 abs_ev   =  4.52E-01 rel_ev   =  8.29E-03
 ediff =        2.69E-01 abs_dens =  2.82E-01 rel_dens =  5.65E-02
Matrix vector products:    108
Converged eigenvectors:      0

#  State  Eigenvalue [eV]  Occupation    Error
      1      -17.468243    2.000000   (5.8E-04)
      2       -6.518418    1.000000   (1.3E-04)
      3       -6.518418    1.000000   (1.3E-04)
      4       -6.518418    1.000000   (1.3E-04)

Density of states:

----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
----------------------------------------------------------------------
---------------------------------------------------------------------%
---------------------------------------------------------------------%
%--------------------------------------------------------------------%
                                                                     ^


Elapsed time for SCF step     1:          1.11
**********************************************************************
{{< /code-block >}}


Now the SCF cycle starts. For every step, {{< octopus >}} outputs several pieces of information:

* The values {{< code-inline >}}abs_dens{{< /code-inline >}} and {{< code-inline >}}rel_dens{{< /code-inline >}} are to monitor the absolute and relative convergence of the density, while {{< code-inline >}}rel_ev{{< /code-inline >}} and {{< code-inline >}}abs_ev{{< /code-inline >}} are two alternative measures of the convergence, based on measuring the difference between input and output eigenvalues. The SCF procedure, by default, is stopped when {{< code-inline >}}rel_dens{{< /code-inline >}} is smaller than $10^{-5}$. This may be altered with the appropriate input variables (see in the manual the variables {{< variable "ConvAbsDens" >}}, {{< variable "ConvRelDens" >}}, {{< variable "ConvAbsEv" >}} and {{< variable "ConvRelEv" >}}).

* The line {{< code-inline >}}Matrix vector products:    108{{< /code-inline >}} tells us that the Hamiltonian was applied 108 times. This gives us an idea of the computational cost.

* The line {{< code-inline >}}Converged eigenvectors:      0{{< /code-inline >}} tells us that upon completion of the diagonalization procedure, none of the orbitals met the required precision criterion for the wavefunctions. In a following example, we will modify this criterion in the input file.

* The list of eigenvalues is then printed, along with their errors: how much they deviate from "exact" eigenvalues of the current Hamiltonian. This number is the so-called "residue".

You can now take a look at the file {{< file "static/info" >}} that will hold a summary of the calculation.

## Restarting

Any ground-state calculation may be restarted later (to refine it if it did not converge properly, or with any other purpose), provided that the contents of the {{< code restart >}} directory are preserved. You can try this now, just by running {{< octopus >}} again. You will notice that {{< octopus >}} did not give any warning after the line

{{< code-block >}}
 Info: Loading restart information.
{{< /code-block >}}

This is useful if you change slightly the parameters of the simulation (for example the XC functional or the convergence criteria). If you change the grid parameters {{< octopus >}} will not be able to restart from the previous calculation. If you do not want {{< octopus >}} to try to restart a calculation, you can set the variable {{< variable "FromScratch" >}}.

In case you ware wondering what the restart information looks like, you can have a look at the contents of the {{< file "restart" >}} directory. This is where the files needed to restart a calculation are stored. It may contain several sub-directories depending on the calculations previously performed. In this case, it just contains one:

```bash
 % ls restart
 gs
 % ls restart/gs
 0000000001.obf  0000000003.obf  density      df_010101.obf  df_010103.obf  dv_010102.obf  f_old_0101.obf  lxyz.obf  mixing  states  vhxc.obf          wfns
 0000000002.obf  0000000004.obf  density.obf  df_010102.obf  dv_010101.obf  dv_010103.obf  grid            mesh      occs    vhxc    vin_old_0101.obf
```

{{< octopus >}} stores each individual state in a different binary (yet platform-independent) file. In this case, we only have four states (files {{< file "0000000001.obf" >}} to {{< file "0000000004.obf" >}}). Some other useful quantities, like the density, are also stored in binary form. The other files are text files that contain diverse control information. It is unlikely that you will ever have to work directly with these files, but you may take a look around if you are curious. 


{{< tutorial-footer >}}

---
Title: "Verlet Propagation"
tutorials: "Multisystem"
weight: 4
---

We are now ready to perform the time propagation of the system. 
The following input file 

{{% expand "input file" %}}
{{< code-block >}}
#include_input testsuite/multisystem/04-propagation_verlet.02-three_body.inp
{{< /code-block >}}
{{% /expand %}}

The input file specifies the propagator as:
{{< code-block >}}
{{< variable "TDSystemPropagator" >}} = verlet
{{< /code-block >}}

According to [Wikipedia](https://en.wikipedia.org/wiki/Verlet_integration), the Verlet algorithm is defined as:

- Calculate $\vec{x}(t + \Delta t) = \vec{x}(t) + \vec{v}(t) \Delta t + \tfrac12 \vec{a}(t) \Delta t^2$.
- Calculate interaction potential using $\vec{x}(t + \Delta t)$.
- Derive $\vec{a}(t + \Delta t)$ from the interaction potential using $\vec{x}(t + \Delta t)$.
- Calculate $\vec{v}(t + \Delta t) = \vec{v}(t) + \tfrac12 \big(\vec{a}(t) + \vec{a}(t + \Delta t)\big)\Delta t$.

These 4 steps have to be performed for each time step.

{{% expand "Detailled description of the algorithmic steps" %}}
The following graph depicts all algorithmic steps of the state machine (see the {{% developers "Code_Documentation/Propagators/" "description of the propagator implementation" %}}).
{{< d3-sequence file="/develop/graph_data/propagation-3body-verlet-equal-step.json" >}}

Each box in this diagram represents a function call of Octopus, which is related to the algorithmic steps of the propagation.
It can be seen how all three systems step through the above algorithmic steps of the Verlet algorithm. You can use the (+) and (-) buttons to manually step through the algorithm to see in which order the steps are performed, and in which steps the various clocks of the systems are advanced.
The interactive graph also allows to display the (unphysical) container systems and the so-called ''ghost'' interactions, which are always present and are an internal tool to ensure that systems stay synchronized, even if they do not have a physical interaction between them.

You can set the variable {{% code-inline %}}{{< variable "Debug" >}} = info{{% /code-inline %}} in the input file and (re-)run the example. Look for the ''Debug:'' lines in the output can compare to the graph above.
{{% /expand %}}


{{< tutorial-footer >}}
---
Title: "Interaction graph"
tutorials: "Multisystem"
weight: 3
---

The code can generate interaction graphs, which can help debugging, if results are not as expected. This can be enabled with the {{< variable "Debug" >}} variable set to {{< code "interaction_graph" >}}.

The following input file 

{{% expand "input file" %}}
{{< code-block >}}
#include_input testsuite/multisystem/02-interaction_graph.01-three_body.inp
{{< /code-block >}}
{{% /expand %}}

gives rise to the following graph:
{{% graphviz-file "/static/graph_data/interaction_graph.dot" %}}

This graph can be created using the {{< code "graphviz" >}} package.

{{% expand "How to install and use graphviz" %}}
On {{% name Debian %}}-like systems, {{< code graphviz >}} can simply be installed with
{{% code-block %}}
sudo apt-get install graphviz
{{% /code-block %}}
Then, the {{< code "interaction_graph.dot" >}} can be converted into a {{< code png >}} file by using:
{{% code-block %}}
dot -Tpng interaction_graph.dot > interaction_graph.png
{{% /code-block %}}

{{% /expand %}}

{{< tutorial-footer >}}
---
Title: "Multisystem"
weight: 40
---


{{< tutorial-series "Multisystem" >}}


---
Title: "Multisystem introduction"
tutorials: "Multisystem"
weight: 1
---

## Octopus with multisystem support

{{< octopus >}} now supports the calculation of compound systems. While the aim of this development is to allow coupling matter to Maxwell fields, or coupling a full DFT system to a system, treated by a tight-binding Hamiltonian, in this tutorial, we will start with a very different model, namely celestial bodies. Due to their simplicity, they are well suited to explain how to set up a system, containing various subsystems.



New variables and concepts in this tutorial: 

* {{< variable "Systems" >}}:
* {{< variable "Interactions" >}}
* Namespaces

### Systems

{{< variable "Systems" >}} describes the nature of a subsystem. Examples are

* classical particles
* charged particles
* a Maxwell system
* a tight binding system
* electrons (not yet implemnented in the multisystem framework)

{{< code-block >}}
%{{< variable "Systems">}}
"Name-1" | {{< emph "system-type-1" >}}
"Name-2" | {{< emph "system-type-2" >}}
...
%
{{< /code-block >}}
The possible system types can be found in the reference of {{< variable "Systems" >}}.
The names, given to a system, define namespaces, which will be used in the remainder of the input file, in order to define system specific values for variables.
This can be seen in the example below.

If a (sub-)system is a multisystem, the components need to be specified in another {{< variable "Systems" >}} block.


### Interactions

Different subsystems are pretty boring, if they cannot interact with each other. These mutual interactions are defined in the {{< variable "Interactions" >}} block.

{{< code-block >}}
%{{< variable "Interactions" >}}
{{< emph "interaction_type" >}} | {{< emph "interaction_mode" >}} [ | {{< emph "partners" >}} ]
...
%
{{< /code-block >}}

### Namespaces

Many possible systems share the same variables. In the following examples of celestial bodies, planets are all classical particles (albeit quite big ones), which are specified by their mass. Obviously, they need to be able to have different masses. For that reason, the system-specific variables are defined by their variable name, with the system name added as prefix (or namespace).

In our example, that would be:

{{< code-block >}}
%{{< variable "Systems">}}
"Earth" | classical_particle
...
%

Earth.{{< variable "ParticleMass" >}} = 5.97237e24
%Earth.{{< variable "ParticleInitialPosition" >}}
 -147364661998.16476 | -24608859261.610123 | 1665165.2801353487
%
%Earth.{{< variable "ParticleInitialVelocity" >}}
 4431.136612956525 | -29497.611635546345 | 0.343475566161544
%
{{< /code-block >}}


{{< tutorial-footer >}}
---
Title: "Example: Celestial dynamics"
tutorials: "Multisystem"
weight: 2
---


In the simplest example we define a system, which consists of similar subsystems. The input file below shows a model describing the motion of the earth, the moon and the sun:

Here, all three subsystems are treated equally. We have two levels of systems: 

* the solar system, which can be thought of as a container
* the celestial bodies (sub, earth and moon)

{{< expand "Expand for input file with two levels of nesting" >}}
{{< code-block >}}
#include_input testsuite/multisystem/01-nested_systems.01-two_levels.inp
{{< /code-block >}}
{{< /expand >}}

While this is a perfectly fine definition of the system, sometimes it is beneficial to group systems together.
For instance, we might define a "Earth" as the combined system of "Terra" and "Luna", which is orbiting the sun:

{{< expand "Expand for input file with three levels of nesting" >}}
{{< code-block >}}
#include_input testsuite/multisystem/01-nested_systems.02-three_levels.inp
{{< /code-block >}}
{{< /expand >}}


{{< tutorial-footer >}}
---
title: "Atomic Simulation Environment (ASE)"
tags: ["Tutorial", "Advanced"]
difficulties: "advanced"
#series: "Tutorial"
description: "Using Octopus as calculator in ASE"
---


## Overview

About ASE: https://wiki.fysik.dtu.dk/ase/about.html

Installation: https://wiki.fysik.dtu.dk/ase/install.html

The Octopus interface is distributed as part of the ASE and does not require any special configuration or recompilation of Octopus. It can be used to set up and run Octopus calculations from a Python script, and these can be combined with different algorithms for structure optimization, molecular dynamics, nudged-elastic-band method for saddle-point searches, vibration/phonon analysis, genetic algorithm, and other features of ASE.

Also, the ASE contains functions to generate structures like crystals, surfaces, or nanoparticles, as well as many standard molecule geometries.

The ASE interface works by writing input files to the disk and running Octopus as an external process. That means some steps like structure optimizations can be less efficient because they require disk I/O to restore quantities such as densities or wavefunctions. This disadvantage can be circumvented by using a ram-disk, if the amount of memory is not too large of course. One can also simply use ASE to script the generation of input files (documentation may be slightly inadequate).

## Exercises

See https://wiki.fysik.dtu.dk/ase/ase/calculators/octopus.html


Back to {{< tutorial "" "Tutorials" >}}



---------------------------------------------
---
title: "Basic QOCT"
tags: ["Tutorial", "Advanced", "Optimal Control", "Molecule", "User-defined Species", "Independent Particles", "Laser"]
difficulties: "advanced"
theories: "Independent particles"
system_types: "Molecule"
calculation_modes: ["Optimal control"]
species_types: "User-defined species"
features: "Laser"
#series: "Tutorial"
description: "Introduction into Quantum Optimal Control Theory."
---


If you set {{< code-inline >}}{{< variable "CalculationMode" >}} = opt_control{{< /code-inline >}}, you will be attempting to solve
the Quantum Optimal Control Theory (QOCT) equations. This tutorial introduces the basic theoretical concepts of
QOCT, as well as the procedure to perform basic QOCT calculations with {{< octopus >}}.

If you already know QOCT, you can skip the Theoretical Introduction.

## Theoretical Introduction

The QOCT can be formulated in the following way:

We consider a quantum mechanical system, governed in a temporal interval $[0,T]$ by the Hamiltonian
$\hat{H}[\varepsilon(t)]$. The Hamiltonian depends parametrically on a "control parameter" $\varepsilon(t)$.
In the cases we are concerned with, this is the temporal shape of an external electromagnetic field, typically
the temporal shape of a laser pulse. This Hamiltonian drives the system from one specified initial state
$\vert\Psi_0\rangle$ to a final state $\vert\Psi(T)\rangle$, according to Schr&ouml;dinger's equation:

$i\frac{\rm d}{{\rm d}t}\vert\Psi(t)\rangle = \hat{H}[\varepsilon(t)]\vert\Psi(t)\rangle\,.$

We now wish to find the ''optimal'' control parameter $\varepsilon$ such that it maximizes a given
objective, which in general should be expressed as a functional, $J_1[\Psi]$ of the full evolution
of the state. This search must be constrained for two reasons: One, the control parameter must be constrained
to ''physical'' values -- that is, the laser field cannot acquire unlimitedly large values. Two, the evolution
of the state must obviously follow Schr&ouml;dinger equation.

These considerations are translated into mathematical terms by prescribing the maximization of the following
functional:

$J[\Psi,\chi,\varepsilon] = J_1[\Psi] + J_2[\varepsilon] + J_3[\Psi,\chi,\varepsilon]\\,.$

* $J_1[\Psi]$ is the objective functional. In this basic tutorial, we will be concerned with the case in which it depends only on the value of the state at the end of the propagation, $\Psi(T)$. Moreover, the functional is the expectation value of some positive-definite operator $\hat{O}$: \
\
$J_1[\Psi] = \langle\Psi(T)\vert\hat{O}\vert\Psi(T)\rangle\\,.$

* $J_2[\varepsilon]$ should constrain the values of the control parameter. Typically, the idea is to maximize the objective $J_1$ while minimizing the ''fluence'': \
\
$J_2[\varepsilon] = -\alpha\int_0^T {\rm d}t \varepsilon^2(t)\\,.$ \
\
The constant $\alpha$ is usually called ''penalty''.

* The maximization of $J_3$ must ensure the verification of Schr&ouml;dinger's equation: \
\
$J_3[\Psi,\chi,\varepsilon] = -2 \Im \left[ \int_0^T {\rm d}t\langle\chi(t)\vert i\frac{\rm d}{{\rm d}t} - \hat{H}[\varepsilon(t)] \vert\Psi(t)\rangle \right]$ \
\
The auxiliary state $\chi$ plays the role of a Lagrangian multiplier.

In order to obtain the threesome $(\Psi,\chi,\varepsilon)$ that maximizes this functional we must find the solution to the corresponding Euler equation. The functional derivation leads to the QOCT equations:

$i\frac{\rm d}{{\rm d}t}\vert\Psi(t)\rangle = \hat{H}[\varepsilon(t)]\vert\Psi(t)\rangle\\,,$

$\vert\Psi(0)\rangle = \vert\Psi_0\rangle\\,.$

$i\frac{\rm d}{{\rm d}t}\vert\chi(t)\rangle = \hat{H}^\dagger[\varepsilon(t)]\vert\chi(t)\rangle\\,,$

$\vert\chi(T)\rangle = \hat{O}\vert\Psi(T)\rangle\\,,$

$\alpha\varepsilon(t) = \Im \left[ \langle\chi(t)\vert \frac{\partial \hat{H}}{\partial \varepsilon}\vert\Psi(t)\rangle\right]\\,.$


## Population target: scheme of Zhu, Botina and Rabitz

An important case in QOCT is when we want to achieve a large transition probability from a specific initial state into a final target state. That is, the operator $\hat{O}$ is the projection onto one given ''target'' state, $\Psi_t$:

$J_1[\Psi] = \vert \langle \Psi(T)\vert\Psi_t\rangle \vert^2 $

As observed by Zhu, Botina and Rabitz[^zhu-1998], for this case one
can use an alternative definition of $J_3$:

$ J_3[\Psi,\chi,\varepsilon] = -2\Im [ \langle\Psi(T)\vert\Psi_t\rangle \int_0^T {\rm d}t \langle\chi(t)\vert i\frac{\rm d}{{\rm d}t}-\hat{H}[\varepsilon(t)]\vert\Psi(t)\rangle ] $

This alternative definition has the virtue of decoupling the equations for $\chi$ and $\Psi$:

$i\frac{\rm d}{{\rm d}t}\vert\Psi(t)\rangle = \hat{H}[\varepsilon(t)]\vert\Psi(t)\rangle\\,,$

$\vert\Psi(0)\rangle = \vert\Psi_0\rangle\\,.$

$i\frac{\rm d}{{\rm d}t}\vert\chi(t)\rangle = \hat{H}^\dagger[\varepsilon(t)]\vert\chi(t)\rangle\\,,$

$\vert\chi(T)\rangle = \vert\Psi_t\rangle\\,,$

$\alpha\varepsilon(t) = \Im \left[ \langle\chi(t)\vert \frac{\partial \hat{H}}{\partial \varepsilon}\vert\Psi(t)\rangle\right]\\,.$

These are the equation that we are concerned with in the present example.

## The asymmetric double well: preparing the states

We are now going to solve the previous equations, for the case of attempting a transition from the ground state of a given system, to the first excited state. The system that we have chosen is a one-dimensional asymmetric double well, characterized by the potential:

$V(x) = \frac{\omega_0^4}{64B}x^4 - \frac{\omega_0^2}{4}x^2 + \beta x^3 $

Let us fix the parameter values to $\omega_0 = B = 1$ and $\beta=\frac{1}{256}$. The ground state of this system may be obtained with {{< octopus >}} running the following inp file:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs 
 {{< variable "Dimensions" >}} = 1 
 {{< variable "FromScratch" >}} = yes
 {{< variable "BoxShape" >}} = sphere
 {{< variable "Spacing" >}} = 0.1
 {{< variable "Radius" >}} = 8.0
 {{< variable "TheoryLevel" >}} = independent_particles
 %{{< variable "Species" >}}
 "AWD1D" | species_user_defined | potential_formula | "1/64*(x)^4-1/4*(x)^2+1/256*(x)^3" | valence | 1
 %
 %{{< variable "Coordinates" >}}
 "AWD1D" | 0 
 %
 {{< variable "ExtraStates" >}} = 1
 {{< variable "EigensolverTolerance" >}} = 1.0e-6
 {{< variable "EigensolverMaxIter" >}} = 1000
 {{< variable "Mixing" >}} = 1.0
 {{< variable "MixingScheme" >}} = linear
{{< /code-block >}}

{{< figure src="/images/Tutorial_basicqoct_potential.png" width="500px" caption="Fig. 1. Black curve: external potential defining the asymmetric double well used in this example. Red and green: ground state and excited state, respectively, corresponding to this system" >}}

Note that we have set {{< code-inline >}}{{< variable "ExtraStates" >}}=1{{< /code-inline >}} because we want to get the first excited state, in order to be able to do QOCT later (the objective will be the maximum population of this state). 

In this way we obtain the first two (single particle) states for this system. It is useful to make a plot with the results, together with the potential. If you do this, you should the curves depicted in Fig. 1. Observe how the ground state is localized in the left well (which is, not surprisingly, deeper), whereas the first excited state is localized in the right well. We are going to attempt to transfer the wavepacket from one side to the other. The energies of these two states are -0.6206 a.u. and -0.4638, respectively.

Note that, in order to plot the functions, you will have to add to the {{< file "inp" >}} file the variables {{< variable "Output" >}} and {{< variable "OutputFormat" >}}.

## Running {{< octopus >}} in QOCT mode

We now run the code with the following {{< file "inp" >}} file:

{{< code-block >}}
 {{< variable "ExperimentalFeatures" >}} = yes
 {{< variable "Dimensions" >}} = 1
 {{< variable "FromScratch" >}} = yes
 {{< variable "CalculationMode" >}} = opt_control
 ###################
 # Grid
 ###################
 {{< variable "BoxShape" >}} = sphere
 {{< variable "Spacing" >}} = 0.1
 {{< variable "Radius" >}} = 8.0
 ###################
 # System
 ###################
 %{{< variable "Species" >}}
 "ADW1D" | species_user_defined | potential_formula | "1/64*(x)^4-1/4*(x)^2+1/256*(x)^3" | valence | 1
 %
 %{{< variable "Coordinates" >}}
 "ADW1D" | 0
 %
 {{< variable "TheoryLevel" >}} = independent_particles
 ###################
 # TD RUN Parameters
 ###################
 stime = 100.0
 dt = 0.01
 {{< variable "TDPropagator" >}} = aetrs
 {{< variable "TDExponentialMethod" >}} = taylor
 {{< variable "TDExponentialOrder" >}} = 4
 {{< variable "TDLanczosTol" >}} = 5.0e-5
 {{< variable "TDMaxSteps" >}} = stime/dt
 {{< variable "TDTimeStep" >}} = dt
 ###################
 # OCT parameters
 ###################
 {{< variable "OCTPenalty" >}} =1.0
 {{< variable "OCTEps" >}} = 1.0e-6
 {{< variable "OCTMaxIter" >}} = 50
 {{< variable "OCTInitialState" >}} = oct_is_groundstate
 {{< variable "OCTTargetOperator" >}} = oct_tg_gstransformation
 %{{< variable "OCTTargetTransformStates" >}}
 0 | 1
 %
 {{< variable "OCTScheme" >}} = oct_zbr98
 {{< variable "OCTDoubleCheck" >}} = yes
 ###################
 # Laser field = Initial guess
 ###################
 ampl = 0.06
 freq = 0.157 
 %{{< variable "TDExternalFields" >}}
 electric_field | 1 | 0 | 0 | freq | "envelope_function"
 %
 %TDFunctions
 "envelope_function" | tdf_cw | ampl
 %
 ###################
 # Output
 ###################
 {{< variable "Output" >}} = wfs
 {{< variable "OutputFormat" >}} = axis_x
 {{< variable "TDOutput" >}} = laser + td_occup
{{< /code-block >}}

Let us go through some of the key elements of this {{< file "inp" >}} file:

* The first part of the file is just a repetition of the basic information already given in the {{< file "inp" >}} file for the ground state calculation: definition of the system through the {{< variable "Species" >}} block, specification of the mesh, etc.

* Each QOCT run consists of a series of backward-forward propagations for both $\Psi$ and $\chi$; the variables in the section "TD Run parameters" specify how those propagations are done. It is especially important to set the total time of the propagation, which fixes the interval $[0,T]$ mentioned above in the QOCT problem definition.

* The next section ("OCT parameters") define the particular way in which OCT will be done:

 1. {{< variable "OCTScheme" >}} is the particular OCT algorithm to be used. For this example, we will use the one suggested by Zhu, Botina and Rabitz[^zhu-1998]. This choice is "oct_zbr98".
 2. {{< variable "OCTMaxIter" >}} sets a maximum number of iterations for the the ZBR98 algorithm. Each iterations involves a forward and a backward propagation of the wavefunction.
 3. {{< variable "OCTEps" >}} sets a convergence threshold for the iterative scheme. This iterative scheme can be understood as a functional that processes an input control parameter $\varepsilon_i$ and produces an output control parameters $\varepsilon_o$, Upon convergence, input and output should converge and the iterations stopped. This will happen whenever<br>
$D[\varepsilon_i,\varepsilon_o] = \int_0^T \vert \varepsilon_o(t)-\varepsilon_i(t)\vert^2\,.$<br>
is smaller than {{< variable "OCTEps" >}}. Note that this may never happen; there is no convergence guarantee.
 4. {{< variable "OCTPenalty" >}} sets the $\alpha$ penalty parameter described previously.
 5. {{< variable "OCTInitialState" >}} should specify the initial states from which the propagation starts. In most cases this will just be the ground state of the system -- this is the case for this example. The value "oct_is_groundstate" means that.
 6. {{< variable "OCTTargetOperator" >}} sets which kind of target operator ($\hat{O}$, with the notation used above) will be used. As mentioned before, for this example we will be using a very common kind of target: the projection onto a target state: $\hat{O}=\vert\Psi_t\rangle\langle\Psi_t\vert$. This is set by doing {{< variable "OCTTargetOperator" >}} = oct_tg_gstransformation. This means that the target state is constructed by taking a linear combination of the spin-orbitals that form the ground state. The particular linear combination to be used is set by the block {{< variable "OCTTargetTransformStates" >}}.
 7. {{< variable "OCTTargetTransformStates" >}}. This is a block. Each line of the code refers to each of the spin orbitals that form the state the system; in this example we have a one-particle problem, and therefore one single spin-orbital. This spin orbital will be a linear combination of the spin-orbitals calculated in the previous ground-state run; each of the columns will contain the corresponding coefficient. In our case, we calculated two spin-orbitals in the ground state run, therefore we only need two columns. And our target state is simply the first excited state, so the first column (coefficient of the ground state orbital) is zero, and the second column is one.
 8. {{< variable "OCTDoubleCheck" >}} simply asks the code to perform a final propagation with the best field it could find during the QOCT iterative scheme, to check that it is really working as promised.

* The next section ("Laser field = Initial guess") simply specifies which is the initial laser field to be used as "initial guess". It is mandatory, since it not only provides with the initial guess for the "control parameter" (which is simply a real function $\varepsilon(t)$ that specifies the temporal shape of the time-dependent external perturbation), but also specifies what kind of external perturbation we are talking about: an electric field with x, y or z polarization, a magnetic field, etc.
* Finally, we would like to do cute plots with our results, so we add an "Output" section. This only applies if we have set {{< variable "OCTDoubleCheck" >}}=yes, since the code does not plot anything during the propagations performed in the QOCT iterations. The variables given will permit us to plot the shape of the wavefunction during the propagation, as well as the time-dependent electric field and the projection of the propagating wave function onto the ground state and first excited state.

Even if this is a simple one-dimensional system, the run will take a few minutes because we are doing a lot of propagations. You will see how, at each iteration, some information is printed to standard output, e.g.:

{{< code-block >}}
 ****************** Optimal control iteration -   12 ******************
        => J1       =      0.87422
        => J        =      0.54580
        => J2       =     -0.32842
        => Fluence  =      0.32842
        => Penalty  =      1.00000
        => D[e,e']  =     1.66E-04
 **********************************************************************
{{< /code-block >}}

The meaning of each number should be clear from the notation that we have introduced above. By construction, the algorithm is "monotonously convergent", meaning that the value of $J$ should increase. This does not guarantee that the final value of $J_1$ will actually approach one, which is our target (full overlap with the first excited state), but we certainly hope so.

At the end of the run it will do a "final propagation with the best field" (usually the one corresponding to the last iteration). For this example, you will get something like:

{{< code-block >}}
 *************** Final propagation with the best field ****************
        => J1       =      0.94926
        => J        =      0.61505
        => J2       =     -0.33420
        => Fluence  =      0.33420
 **********************************************************************
{{< /code-block >}}

So we ''almost'' achieved a full transfer: $J_1[\Psi]=\vert\langle\Psi(T)\vert\Psi_t\rangle\vert^2 = 0.94926$.

##  Analysis of the results  

{{< figure src="/images/Tutorial_basicqoct_convergence.png" width="500px" caption="Fig. 2. Convergence history for the $J$ (black) and $J_1$ (red) functionals" >}}
{{< figure src="/images/Tutorial_basicqoct_final.png" width="500px" caption="Fig. 3. Initial and target states (red and green). In blue, the final propagated state corresponding to the last iteration field." >}}

The run will create a directory called {{< file "opt-control" >}} (among others). In it there is a file called {{< file "convergence" >}} that contains the convergence history of the QOCT run, i.e. the values of the relevant functionals at each QOCT iteration. It is useful to make a plot with these values; you can see the result of this in Fig. 2. Observe the monotonous increase of $J$. Also, in this case, the functional $J_1$ (which is, in fact, the one we are really interested in) grows and approaches unity -- although the final limit will not be one, but in this case something close to 0.95.

You may want to get the laser pulses used at each iteration, you can do so by setting the variable {{< variable "OCTDumpIntermediate" >}}.

There is also a directory called {{< file "initial" >}}, useful to plot the initial state; a directory called {{< file "target" >}}, that contains the target state; and a directory called {{< file "final" >}} that contains the state at the end of the propagation corresponding to the last QOCT iteration. Fig. 3 is a plot constructed with the information contained in these directories. Note that while the initial and target states have no imaginary part, the final propagated wavefunction is complex, and therefore we have to plot both the real and the imaginary part.

In the directories {{< file "laser.xxxx" >}} you will find a file called {{< file "cp" >}} which contains the control parameter for the corresponding QOCT iteration. Another file, called {{< file "Fluence" >}} informs about the corresponding fluence. The format of file {{< file "cp" >}} is simple: first column is time, second and third column are the real and imaginary part of the control parameter. Since $\varepsilon$ is defined to be real in our formulation, the last column should be null. Finally, the directory {{< file "laser.bestJ1" >}} will obviously contain the control parameter corresponding to the iteration that produced the best $J_1$.



[^zhu-1998]: W. Zhu, J. Botina and H. Rabitz, J. Chem. Phys. '''108''', 1953 (1998)
---
title: "Sternheimer linear response"
#tags: ["Tutorial", "Advanced", "Electromagnetic Response", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "Sternheimer"]
difficulties: "advanced"
theories: "DFT"
calculation_modes: "Electromagnetic Response"
system_stypes: "Molecule"
species_types: "Pseudopotentials"
features: ["Optical Absorption", "Sternheimer"]
#series: "Tutorial"
description: "More details on the Sternheimer approach."
---


The Sternheimer approach to perturbation theory allows efficient calculations of linear and non-linear response properties.
[^footnote-1]

The basis of this method, just as in standard perturbation theory, is to calculate the variation of the wave-functions $\psi^{1}$ under a given perturbing potential. The advantage of the method is that the variations are obtained by solving the linear equation 

$$
(H^0-\\epsilon^0 + \\omega )|\\psi^{1}\>=-P\_{\\rm c} H^{1}|\\psi^{0}\>\\ ,
$$

that only depends on the occupied states instead of requiring an (infinite) sum over unoccupied states. In the case of (time-dependent) density functional theory the variation of the Hamiltonian includes a term that depends on the variation of the density, so this equation must be solved self-consistently.

To run a Sternheimer calculation with {{< octopus >}}, the only previous calculation you need is a ground-state calculation.
For this tutorial we will use a water molecule, with this basic input file for the ground state:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 
 %{{< variable "Coordinates" >}}
  'O'  |  0.000000  | -0.553586  |  0.000000
  'H'  |  1.429937  |  0.553586  |  0.000000
  'H'  | -1.429937  |  0.553586  |  0.000000
 %
 
 {{< variable "Radius" >}} = 10
 {{< variable "Spacing" >}} = 0.435
 {{< variable "ConvRelDens" >}} = 1e-6
{{< /code-block >}}

We use a tighter setting on SCF convergence (ConvRelDens) which will help the ability of the Sternheimer calculation to converge numerically, and we increase a bit the size of the box as response calculations tend to require more space around the molecule than ground-state calculations to be converged.
[^footnote-2]


After the ground-state calculation is finished, we change the run mode to {{< code "em_resp" >}}, to run a calculation of the electric-dipole response:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = em_resp
{{< /code-block >}}

Next, to specify the frequency of the response we use the {{< variable "EMFreqs" >}} block; in this case we will use three values 0.00, 0.15 and 0.30 {{< units "{{< Hartree >" >}}}}:

{{< code-block >}}
 %{{< variable "EMFreqs" >}} 
 3 | 0.0 | 0.3
 %
{{< /code-block >}}

and we will also specify a small imaginary part to the frequency of 0.1 {{< units "eV" >}}, which avoids divergence on resonance:

{{< code-block >}}
 {{< variable "EMEta" >}} = 0.1*eV
{{< /code-block >}}

and finally we add a specification of the linear solver, which will greatly speed things up compared to the default:

{{< code-block >}}
 {{< variable "LinearSolver" >}} = qmr_dotp
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}

In the run, you will see calculations for each frequency for the ''x'', ''y'', and ''z'' directions, showing SCF iterations, each having linear-solver iterations for the individual states' $\psi^{1}$, labelled by the k-point/spin (ik) and state (ist). The norm of $\psi^{1}$, the number of linear-solver iterations (iter), and the residual $\left|(H^0-\epsilon^0 + \omega )|\psi^{1}>+P_{\rm c} H^{1}|\psi^{0}>\right|$ are shown for each. First we see the static response:

{{< code-block >}}
****************** Linear-Response Polarizabilities ******************
Wavefunctions type: Complex
Calculating response for   3 frequencies.
**********************************************************************

Info: Calculating response for the x-direction and frequency 0.0000.
Info: EM Resp. restart information will be written to 'restart/em_resp'.
Info: EM Resp. restart information will be read from 'restart/em_resp'.

** Warning:
**   Unable to read response wavefunctions from 'wfs_x_f1+': Initializing to zero.

Info: Finished reading information from 'restart/em_resp'.
--------------------------------------------
LR SCF Iteration:   1
Frequency:             0.000000 Eta :             0.003675
   ik  ist                norm   iters            residual
    1    1            0.216316      21        0.110754E-03
    1    2            1.573042      21        0.286815E-02
    1    3            1.620482      21        0.973599E-02
    1    4            1.177410      21        0.540540E-03

             Info: Writing states. 2016/01/14 at 19:50:27


        Info: Finished writing states. 2016/01/14 at 19:50:27

SCF Residual:     0.122899E+01 (abs),     0.153623E+00 (rel)
{{< /code-block >}}

Later will come the dynamical response. The negative state indices listed indicate response for $-\omega$. For each frequency, the code will try to use a saved response density from the closest previously calculated frequency.

{{< code-block >}}
Info: Calculating response for the x-direction and frequency 0.1500.
Info: EM Resp. restart information will be written to 'restart/em_resp'.
Info: EM Resp. restart information will be read from 'restart/em_resp'.
Read response density 'rho_0.0000_1'.
Info: Finished reading information from 'restart/em_resp'.
--------------------------------------------
LR SCF Iteration:   1
Frequency:             0.150000 Eta :             0.003675
   ik  ist                norm   iters            residual
    1    1            0.181350       8        0.136183E-02
    1   -1            0.240803      19        0.794480E-04
    1    2            0.953447       8        0.508247E-02
    1   -2            1.708845      19        0.166997E-02
    1    3            0.864834       8        0.858907E-02
    1   -3            1.780595      19        0.984578E-02
    1    4            0.725681       8        0.529092E-02
    1   -4            1.421823      19        0.536713E-03
{{< /code-block >}}

At the end, you will have a directory called {{< file "em_resp" >}} containing a subdirectory for each frequency calculated, each in turn containing {{< file "eta" >}} (listing $\eta$ = 0.1 {{< units "eV" >}}), {{< file "alpha" >}} (containing the real part of the polarizability tensor), and {{< file "cross_section" >}} (containing the cross-section for absorption, based on the imaginary part of the polarizability).

For example, {{< file "em_resp/freq_0.0000/alpha" >}} says

{{< code-block >}}
 # Polarizability tensor [b^3]
            10.238694           -0.000000           -0.000000
             0.000000           10.771834           -0.000000
            -0.000000           -0.000000            9.677212
 Isotropic average           10.229247
{{< /code-block >}}

Exercise: compare results for polarizability or cross-section to a calculation from time-propagation or the Casida approach.










---------------------------------------------
[^footnote-1]: {{< article title="Time-dependent density functional theory scheme for efficient calculations of dynamic (hyper)polarizabilities" authors="Xavier Andrade, Silvana Botti, Miguel Marques and Angel Rubio" journal="J. Chem. Phys" volume="126" pages="184106" year="2007" doi="10.1063/1.2733666" >}}

[^footnote-2]: {{< article title="Basis set effects on the hyperpolarizability of CHCl<sub>3</sub>: Gaussian-type orbitals, numerical basis sets and real-space grids" authors="F. D. Vila, D. A. Strubbe, Y. Takimoto, X. Andrade, A. Rubio, S. G. Louie, and J. J. Rehr" journal="J. Chem. Phys." volume="133" pages="034111" year="2010" doi="10.1063/1.3457362" >}}

---
title: "Running Octopus on Graphical Processing Units (GPUs)"
#tags: ["Tutorial", "Expert"]
difficulties: "expert"
#series: "Tutorial"
description: ""
---


{{< notice warning >}}
This chapter needs to be updated, as the GPU implementation has changed since it was written!!!
{{< /notice >}}


Recent versions of Octopus support execution on graphical processing
units (GPUs). In this tutorial we explain how the GPU support works
and how it can be used to speed-up calculations.

{{< octopus >}} is in the process of being ported to machines with graphical processing units (GPUs). Not all operations are supported yet, but the number of supported features is growing with time. 

{{< notice note >}}
Note that the code might fall back to CPU operation for unsupported features. 
{{< /notice >}}

### Consideration before using a GPU

#### Calculations that will be accelerated by GPUs

A GPU is a massively parallel processor, to work efficiently it need
large amounts of data to process. This means that using GPUs to
simulate small systems, like molecules of less than 20 atoms or low
dimensionality systems, will probably will be slower than using the CPU.

Not all the calculations you can do with Octopus will be effectively
accelerated by GPUs. Essentially ground state calculations with the
rmmdiis eigensolver (and calculations based on ground state calculations, like geometry optimization) and time propagation with the etrs and aetrs
propagators are the simulations that will use the GPU more
effectively. For other types of calculations you might see some
improvements in perfomance. Moreover, some Octopus features do not
support GPUs at all. Currently these are:

* HGH or relativistic pseudopotentials.
* Non-collinear spin.
* Curvilinear coordinates.

* Periodic systems are partially supported, but there are limitations which might force the code to run on the CPU.

In these cases Octopus will be silently executed on the CPU.

#### Supported GPUs

Octopus GPU support is based on the OpenCL framework, so it is vendor
independent, currently it has been tested on GPUs from Nvidia and AMD
(ex-ATI). Recently, CUDA support is added using a CUDA wrapping layer.

{{< notice warning >}}
Actually, the non-CUDA version is currently untested, as our development systems are all NVidia CUDA based.
{{< /notice >}}

In order to run Octopus, the GPU must have double precision
support. For AMD/ATI cards this is only available in high-end models:
Radeon 58xx, 69xx, 78xx, 79xx, and Firepro/Firestream cards. All
recent Nvidia cards support double precision (4xx and newer), low- and
middle-end cards however have very limited processing power and are
probably slower running Octopus than a CPU.

The other important factor is GPU memory, you need a GPU with at least
1.5 GiB of RAM to run Octopus, but 3 GiB or more are recommended. 


#### Activating GPU support

Octopus has GPU support since version 4.0; to activate it you need a version compiled with OpenCL support enabled (see the {{< manual "Building_Octopus_with_GPU_support" "manual" >}} for details). If your version of Octopus was compiled with OpenCL support you should see that 'opencl' among the configuration options:

```text
                            Running octopus
 
 Version                : superciliosus
 Revision               :
 Build time             : Thu Feb 14 22:19:10 EST 2013
 Configuration options  : max-dim=3 openmp '''opencl''' sse2 avx
 Optional libraries     : netcdf '''clamdfft clamdblas'''
```

Ideally your Octopus version should also be compiled with the optional libraries 'clamdfft' and 'clamdblas', that are part of the [https://developer.amd.com/tools/heterogeneous-computing/amd-accelerated-parallel-processing-math-libraries/ AMD APPML] package (these libraries work on hardware from other vendors too, including Nvidia).

### Starting the tutorial: running without OpenCL

If you have a version of Octopus compiled with OpenCL support, by setting the {{< variable "DisableAccel" >}} to {{< code "yes" >}} you can tell Octopus not to use OpenCL.

### Selecting the OpenCL device

OpenCL is designed to work with multiple devices from different vendors. Most likely you will have a single OpenCL device, but in some cases you need to select the OpenCL device that Octopus should use.

In OpenCL a 'platform' identifies an OpenCL implementation, you can have multiple platforms installed in a system. Octopus will list the available platforms at the beginning of the execution, for example:

```text
 Info: Available CL platforms: 2
     * Platform 0 : NVIDIA CUDA (OpenCL 1.1 CUDA 4.2.1)
       Platform 1 : AMD Accelerated Parallel Processing (OpenCL 1.2 AMD-APP (1084.4))
```

By default Octopus will use the first platform. You can select a different platform with the 
{{< max-version 9 >}}{{< variable "OpenCLPlatform" >}}{{< /max-version >}}{{< min-version 10 >}}{{< variable "AccelPlatform" >}}{{< /min-version >}} 
variable. You can select the variable by name, for example:

{{< min-version "10" >}}
{{< code-block >}}
{{< variable "AccelPlatform" >}} = amd
{{< /code-block >}}
{{< /min-version >}}

{{< max-version "9" >}}
{{< code-block >}}
{{< variable "OpenCLPlatform" >}} = amd
{{< /code-block >}}
{{< /max-version >}}

or you can use numeric values that select the platform by its index in the list:

{{< min-version "10" >}}
{{< code-block >}}
{{< variable "AccelPlatform" >}} = 1
{{< /code-block >}}
{{< /min-version >}}

{{< max-version "9" >}}
{{< code-block >}}
{{< variable "OpenCLPlatform" >}} = 1
{{< /code-block >}}
{{< /max-version >}}

The first form should be preferred as it is independent of the order of the platforms.

Each OpenCL platform can contain several devices, the available devices on the selected platform are printed by Octopus, like this:

```text
 Info: Available CL devices: 2
       Device 0 : Tahiti
       Device 1 : Intel(R) Core(TM) i7-3820 CPU @ 3.60GHz
```

Each OpenCL device has a certain type that can be 'cpu', 'gpu', or 'accelerator'. By default Octopus tries to use the first device of type 'gpu', if this fails Octopus tries to use the default device for the system, as defined by OpenCL implementation. To select a different device you can use the 
{{< max-version 9 >}}{{< variable "OpenCLDevice" >}}{{< /max-version >}}{{< min-version 10 >}}{{< variable "AccelDevice" >}}{{< /min-version >}} variable, you can select the device by its index or by its type ('cpu', 'gpu', or 'accelerator').

{{< notice note >}}
the Octopus OpenCL implementation is currently designed to execute efficiently on GPUs. Octopus will run on CPUs using OpenCL, but it will be slower than using the CPU directly.
{{< /notice >}}

After Octopus has selected an OpenCL device it will print some information about it, for example:

```text
 Selected CL device:
      Device vendor          : Advanced Micro Devices, Inc.
      Device name            : Tahiti
      Driver version         : 1084.4 (VM)
      Compute units          : 32
      Clock frequency        : 925 GHz
      Device memory          : 2048 MiB
      Max alloc size         : 512 MiB
      Device cache           : 16 KiB
      Local memory           : 32 KiB
      Constant memory        : 64 KiB
      Max. workgroup size    : 256
      Extension cl_khr_fp64  : yes
      Extension cl_amd_fp64  : yes
```



### The block-size

To obtain good peformance on the GPU (and CPUs), Octopus uses blocks of states as the basic data object. The size of these blocks is controlled by the {{< variable "StatesBlockSize" >}} variable. For GPUs the value must be a power of 2 (the default is 32) and for optimal performance must be as large as possible, the catch is that for large values you might run out of GPU memory.

### StatesPack

### The Poisson solver

If you want to used a GPU-accelerated Poisson solver, you need to use the 'fft' solver and to have the 'clamdfft' library support compiled. To select the solver, just add

```text
 {{< variable "PoissonSolver" >}} = fft
 {{< variable "FFTLibrary" >}} = clfft
```

to your input file.


<span class=noprint><hr>
Back to [[Tutorials]]



---------------------------------------------
---
title: "Large systems: the Fullerene molecule"
tags: ["Tutorial", "Advanced", "Ground State", "Molecule", "Pseudopotentials", "DFT"]
difficulties: "advanced"
theories: "DFT"
calculation_modes: "Ground state"
system_types: "Molecule"
species_types: "Pseudopotentials"
#series: "Tutorial"
description: "How to set up calculations for larger systems."
---


In this tutorial, by solving the Fullerene molecule, we show you how to do simulations of systems with a large number of atoms with {{< octopus >}}. Most of the default parameters work fine for large systems, but there are a few things you might need to adjust.

### The eigensolver

The default eigensolver of {{< octopus >}} is simple and reliable, but for large systems (> 50 orbitals) it can be quite slow. For these systems it is better to use the RMMDIIS eigensolver. This is a very efficient diagonalizer, but it can have convergence problems. In {{< octopus >}} we follow the algorithm implementation of Kresse and Furthmüller [^footnote-1] that provides a very robust algorithm if the following considerations are taken:

* A good initial approximation of the eigenvalues and eigenvectors is required. In {{< octopus >}}, this means you have to do an LCAO before using the RMMDIIS eigensolver (this is the default). You might as well use an alternative eigensolver for one or two iterations.

* You have to provide a reasonable number of {{< variable "ExtraStates" >}}, it should be around 10% to 20% of the number of occupied states. You must include them by hand.

Let's see an example: this is an input file, {{< file "inp" >}}, to calculate the ground state of the fullerene molecule with the RMMDIIS eigensolver:

{{< code-block >}}
 {{< variable "XYZCoordinates" >}} = "C60.xyz"
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 {{< variable "CalculationMode" >}} = gs
 {{< variable "EigenSolver" >}} = rmmdiis
 {{< variable "ExtraStates" >}} = 20
{{< /code-block >}}

In this case, to use RMMDIIS the only thing we need to do is to change the {{< variable "Eigensolver" >}} variable and to include 20 {{< variable "ExtraStates" >}} (by default {{< octopus >}} does not include any unoccupied states).

This is the content of the coordinates file, {{< file "C60.xyz" >}}:

{{< code-block >}}
60

 C    -1.415441     3.015027    -1.167954
 C     2.583693     2.292084    -0.721147
 C     0.721155    -2.583700     2.292090
 C     3.015024    -1.167941     1.415429
 C    -3.015019    -1.167950     1.415434
 C    -0.721134    -2.583697     2.292085
 C    -2.583698     2.292079    -0.721148
 C    -3.461397    -0.000006     0.694328
 C     1.415423     3.015023    -1.167951
 C     3.461399     0.000011     0.694319
 C    -1.415430    -3.015023    -1.167948
 C    -2.292086    -0.721156    -2.583690
 C    -0.000003     0.694321    -3.461398
 C     2.292087    -0.721146    -2.583706
 C     1.415436    -3.015019    -1.167948
 C    -2.292087     0.721141    -2.583692
 C     1.167941     1.415428    -3.015025
 C     3.015020    -1.167945    -1.415437
 C     0.694330    -3.461397     0.000006
 C    -2.583701    -2.292092    -0.721145
 C     0.721144     2.583697     2.292082
 C     1.167948     1.415436     3.015018
 C    -0.000002     0.694327     3.461397
 C    -1.167951     1.415435     3.015031
 C    -0.721152     2.583693     2.292084
 C    -0.694329     3.461397    -0.000006
 C     2.583701     2.292093     0.721144
 C     2.292088    -0.721140     2.583692
 C    -1.167941    -1.415428     3.015025
 C    -3.015020     1.167945     1.415437
 C     1.167955    -1.415429     3.015019
 C    -2.292077    -0.721143     2.583701
 C    -2.583697     2.292080     0.721148
 C     0.694318     3.461398    -0.000005
 C     3.015030     1.167958     1.415430
 C    -0.694317    -3.461398     0.000005
 C    -3.015030    -1.167958    -1.415430
 C    -1.167955     1.415429    -3.015019
 C     2.292078     0.721143    -2.583701
 C     2.583697    -2.292079    -0.721148
 C     1.167951    -1.415434    -3.015031
 C     0.721152    -2.583693    -2.292084
 C    -0.721144    -2.583696    -2.292082
 C    -1.167948    -1.415436    -3.015019
 C     0.000002    -0.694327    -3.461397
 C     0.721134     2.583697    -2.292085
 C     3.461397     0.000007    -0.694329
 C     1.415440    -3.015027     1.167954
 C    -2.583694    -2.292084     0.721147
 C    -3.015025     1.167942    -1.415429
 C    -2.292087     0.721146     2.583706
 C    -1.415436     3.015019     1.167948
 C     1.415430     3.015023     1.167947
 C     2.292086     0.721156     2.583690
 C     0.000003    -0.694321     3.461398
 C    -1.415423    -3.015023     1.167952
 C    -3.461399    -0.000011    -0.694320
 C    -0.721155     2.583700    -2.292090
 C     3.015019     1.167950    -1.415434
 C     2.583698    -2.292079     0.721148
{{< /code-block >}}

Create both files and run {{< octopus >}}. It should take a few minutes (depending on the speed of your machine). There are some peculiarities of the RMMDIIS eigensolver that you might notice:

* The first self-consistency iterations take less time than the rest of the iterations. To ensure a proper initial approximation, the eigesolver does a different diagonalization during the first steps. 

* With respect to other eigensolvers, it takes many more self-consistency iterations to converge the calculation. In compensation each RMMDIIS step is faster. This is because at each step a small and fixed number of eigensolver steps are done, so in practice the eigensolver and self-consistency convergence are done together.

* At some iterations the eigenvalues might not be ordered -- this is normal. When convergence is reached, the eigenvalues will be properly ordered.

### Self-consistency convergence

Once your calculation is finished, you might notice that the final residue in the eigenvectors is rather large (of the order of $10^{-4}$). This is for two reasons: first, it is possible to converge the density without properly converging the eigenvectors sometimes. Try a tighter tolerance and rerun (with restart), adding

{{< code-block >}}
 {{< variable "ConvRelDens" >}} = 1e-6
{{< /code-block >}}

to your {{< file "inp" >}} file. Now your eigenvectors should be better converged.

Can you converge all the eigenvectors? How many states are actually properly converged? Suppose you wanted to calculate correct unoccupied states. How would you do it and be sure they were converged?

The second reason for poor convergence is that we are using a large (default) spacing -- check what it is in the output. The problem is just not as well-conditioned with a large spacing. So, you can try using a smaller spacing such as the one you found in the [methane tutorial](../Methane_molecule ).


[^footnote-1]: [Kresse and Furthmüller](https://dx.doi.org/10.1103/PhysRevB.54.11169) 





---------------------------------------------
---
title: "Vibrational modes"
#tags: ["Tutorial", "Advanced", "Vibrational Modes", "Geometry Optimization", "Molecule", "Pseudopotentials", "DFT", "Forces", "Vibrations"]
#series: "Tutorial"
theories: "DFT"
calculation_modes: ["Geometry Optimization"]
system_types: "Molecule"
species_types: "Pseudopotentials"
features: ["Vibrational Modes","Forces", "Vibrations"]
description: "Obtain vibrational modes by using Sternheimer linear response."
---


This tutorial will show you how to obtain vibrational modes in {{< octopus >}} by using Sternheimer linear response. We will use the water molecule as example. To start we will need the geometry. Copy the following in a file called {{< file "h2o.xyz" >}}:

{{< code-block >}}
 3
 Water molecule in Angstroms 
 O  0.000 0.000 0.0
 H  0.757 0.586 0.0
 H -0.757 0.586 0.0
{{< /code-block >}}

Since to obtain proper vibrational modes we need to be as close as possible to the minimum of energy, first we will optimize the geometry under the parameters we have defined for the simulation. For that we use the following input file:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = go
 
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 {{< variable "XYZCoordinates" >}} = "h2o.xyz"
 
 {{< variable "Spacing" >}} = 0.16*angstrom
 {{< variable "Radius" >}} = 4.5*angstrom
 {{< variable "BoxShape" >}} = minimum

 {{< variable "FilterPotentials" >}} = filter_ts
 
 {{< variable "GOMethod" >}} = fire
 
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}


As you see we will use a smaller than usual spacing and a larger spherical box. This is required since forces require more precision than other quantities. Of course for production runs these values should be properly converged. Also to increase the precision of the forces and the total energy we will use the variable {{< variable "FilterPotentials" >}}, which specifies that a filter should be applied to the ionic potentials that removes the Fourier components that are higher than what we can represent on our grid.

Now run {{< octopus >}}, once the calculation is finished you should get the optimized geometry in {{< file "min.xyz" >}}. Modify the value of {{< variable "XYZCoordinates" >}} to point to this file.

To perform a calculation of the vibrational modes we need a ground-state calculation with a very well converged forces. So change the {{< variable "CalculationMode" >}} to {{< code "gs" >}} and run {{< octopus >}} again, by adding {{< variable "ConvForce" >}} = 1.e-7 ( the criterion to stop the self-consistency process will be given by the change in the force). Once done, check that the forces in the {{< file "static/info" >}} file are reasonably small (<0.05 eV/A).

Now we are ready for the actual calculation of vibrational modes. To do it, change {{< variable "CalculationMode" >}} to {{< code "vib_modes" >}} and run {{< octopus >}}. This calculation will take a while, since $3N_{atoms}$ response calculations are required.

Finally, the results will be written into the {{< code "vib_modes" >}} folder.

- Check the {{< file "normal_frequencies_lr" >}} file. 

{{< code-block >}}
     1    3714.66865804
     2    3612.19180564
     3    1541.89763162
     4     282.95538165
     5     189.98486286
     6     155.62791393
     7    -197.74013582
     8    -246.65205687
     9    -258.63857751
{{< /code-block >}}

What kind of information does it give? Did we find the correct equilibrium geometry? 

Try now to visualize the vibrational eigenmodes. You can use XCrySDen to open the file {{< file "normal_modes_lr.axsf" >}}. In the "Display" menu set "Forces", and in "Modify" set the factor to 1. Try to classify the modes as translations, rotations, and true internal vibrations. Which ones have the "negative" (actually imaginary) frequencies?---
title: "Unsorted tutorials"
weight: 500
---

{{< children >}}---
title: "Polarizable Continuum Model (PCM)"
theory: "DFT"
calculation_modes: ["Ground state","Time-dependent"]
system_types: "Molecule"
species_types: "Pseudopotentials"
features: ["Polarizable Continuum Model"]
#series: "Tutorial"
description: "Calculate the solvation energy of the Hydrogen Fluoride molecule in water solution."
---



{{< notice note >}}
Carlo: How about splitting this tutorial in two parts, so it better fits the categorization or the other tutorials? I propose 
part 1: "Solvation energy of a molecule in solution within the Polarizable Continuum Model", and 
part 2: "Time propagation in a solvent within the Polarizable Continuum Model". 

I have already modified the intros, but I haven't split the page yet, so that Gabriel can finish it first
{{< /notice >}}


## Solvation energy of a molecule in solution within the Polarizable Continuum Model

In this tutorial we show how calculate the solvation energy of the Hydrogen Fluoride molecule in water solution by setting up a ground state calculation using the Polarizable Continnum Model (PCM).

At the time being PCM is only implemented for TheoryLevel = DFT.

### Input

The required input keywords to activate a gs PCM run are {{< variable "PCMCalculation" >}} and {{< variable "PCMStaticEpsilon" >}}.

{{< variable "PCMCalculation" >}} just enables the PCM part of the calculation. {{< variable "PCMStaticEpsilon" >}} contains the value of the relative static dielectric constant of the solvent medium (in this particular example water with $\epsilon = 78.39$).

The corresponding minimal input for the ground state calculation in solution is

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom

 {{< variable "PseudopotentialSet" >}} = hgh_lda
 
 FH = 0.917*angstrom
 %{{< variable "Coordinates" >}}
  "H" |   0 | 0 | 0
  "F" |  FH | 0 | 0
 %

 {{< variable "PCMCalculation" >}} = yes
 {{< variable "PCMStaticEpsilon" >}} = 78.39
{{< /code-block >}}

#### Warnings

* Please note that the latter input example is intended to make you familiar with the running of Octopus for PCM calculations, it does not provide any benchmark for HF molecule. 
* Actually, the geometry might not be optimal and the gs energy might not be converged with respect to the Radius and Spacing variables. A new convergence analysis w.r.t. Radius and Spacing (independent of the one in vacuo) is in principle due when PCM is activated. 
* Geometry optimization for solvated molecules within Octopus is still under development; it might work properly only when there are minor adjustments in the geometry when passing from vacuuum to the solvent. 
* The selection of the PseudopotentialSet and the implicit choice of LDA xc functional for this particular example is arbitrary.

### Output

The main output from such a gs PCM calculation is the solvation energy of the system, which can be extracted from the file

* {{< file "static/info" >}}

The actual info file produced from a calculation with the previous output contains the following energetic contributions:

{{< code-block >}}
 Energy [eV]:
       Total       =      -680.55034366
       Free        =      -680.55034366
       -----------
       Ion-ion     =       109.92094773
       Eigenvalues =      -120.75757382
       Hartree     =       708.30289945
       Int[n*v_xc] =      -176.84277019
       Exchange    =      -121.31897823
       Correlation =       -13.43829837
       vanderWaals =         0.00000000
       Delta XC    =         0.00000000
       Entropy     =         0.00000000
       -TS         =        -0.00000000
       E_e-solvent =         1.44141920
       E_n-solvent =        -2.34199070
       E_M-solvent =        -0.90057151
       Kinetic     =       556.30947697
       External    =     -1919.42310315
       Non-local   =        42.69395221
{{< /code-block >}}

E_M-solvent is the solvation energy, which is in turn split in a term due to electrons and another one due to the nuclei, E_e-solvent and E_n-solvent, respectively.

{{< notice note >}}
Carlo: add the first few lines of each output file
{{< /notice >}}

When we perform the calculation of HF in water using the previous input, a new folder called pcm is created inside the working directory where input file lies. This pcm folder contains the following files: 

* {{< file "cavity_mol.xyz" >}} - useful way to visualize the molecule inside the cavity tessellation. The cavity is plotted as a collection of fictitious Hydrogen atoms. Checking up this file is good to be sure that the molecule fits well inside the cavity and that the tessellation doesn't contain artifacts.

The first lines of the files are:

{{< code-block >}}
    62

   H        1.46477445     1.44399256    -0.00000000
   H        1.08627161     1.44399256     0.52096446
   H        0.47384116     1.44399256     0.32197374
{{< /code-block >}}


* {{< file "pcm_info.out" >}} - a summary of gs PCM calculation containing information on the tessellation and a table of PCM energetic contributions and total polarization charges due electrons and nuclei for all iterations of the SCF cycle. The PCM energetic contributions are split into the interaction energy of: 1) electrons and polarization charges due to electrons E<sub>ee</sub>, 2) electrons and polarization charges due to nuclei E<sub>en</sub>, 3) nuclei and polarization charges due to nuclei E<sub>nn</sub>, and finally, 4) nuclei and polarization charges due to electrons. 

Below, an example of the file:

{{< code-block >}}
# Configuration: Molecule + Solvent
################################### 
# Epsilon(Solvent) =       78.390
#
# Number of interlocking spheres =    1
#
# SPHERE   ELEMENT               CENTER  (X,Y,Z) (A)                    RADIUS (A)
# ------   -------    -------------------------------------------       ----------
#    1       F        0.91700000      0.00000000      0.00000000        1.54440000
#
# Number of tesserae / sphere = 60
#
# 
# Total number of tesserae =   60
# Cavity surface area (A^2) =    29.972947
# Cavity volume (A^3) =    15.430073
# 
#
#    iter              E_ee                     E_en                     E_nn                     E_ne                     E_M-solv                 Q_pol                    deltaQ^e                 Q_pol                    deltaQ^n
       1            -293.25112588             294.71995750            -296.89466500             294.71995750              -0.70587588              -7.87251916              -0.02542700               7.89664876              -0.00129741
       2            -293.13633571             294.54433718            -296.89466500             294.54433718              -0.94232634              -7.87358521              -0.02436096               7.89664876              -0.00129741
       3            -293.08877935             294.51502795            -296.89466500             294.51502795              -0.95338845              -7.87305159              -0.02489458               7.89664876              -0.00129741
{{< /code-block >}}

* {{< file "pcm_matrix.out" >}} - the PCM response matrix. It is printed for advanced debug and specially to check against GAMESS PCM matrix in the same conditions (same solvent dielectric constant, same geometry of the cavity).

* {{< file "ASC_e.dat" >}} - polarization charges due to the electrons.

* {{< file "ASC_n.dat" >}} - polarization charges due to the nuclei.

* {{< file "ASC_sol.dat" >}} - polarization charges due to the full molecule.

{{< notice note >}}
Carlo: a word of warning about the signs of these charges here
{{< /notice >}}

The format of the polarization charge file is:

{{< code-block >}}
 cavity position x | cavity position y | cavity position z | polarization charge | label 
{{< /code-block >}}

where | indicates column separation (actually, replace by spaces in the ASC_*.dat files).

### Exercises

{{< figure src="/images/PCM_gs_convergence.png" width="500px" caption="Convergence of the solvation energy of HF in water solution as a function of the number of SCF iterations." >}}

* Look first at the convergence of the total energy of the system. 
{{< notice note >}}
Carlo: clarify if you mean the energy of the molecule or the solvation energy, and make the figure accordingly
{{< /notice >}} 
It should look like in the Figure. What it is most important to notice is the fact that there is a stabilization of the system in solvent with respect to the ''in vacuo'' case. Find the final value and sign of the energy difference between the both cases obtained upon convergence.
{{< notice note >}}
Carlo: Please rename to something meaningful such as PCM_gs_convergence.png and make the figure larger. Mediawiki will thumbnail it automatically. Also I am not sure if the inset is useful. Maybe suggest plotting in log scale as an exercise
{{< /notice >}}


* Another important check is to ascertain whether Gauss' theorem is fulfilled for the total polarization charges, i.e., if $Q_{pol}=-(\epsilon-1)/\epsilon \times Q_M$ is valid or not, where $Q_M$ is the nominal charge of the molecule. If the latter fails, by setting {{< variable "PCMRenormCharges" >}} = yes and manipulating {{< variable "PCMQtotTol" >}} Gauss' theorem is recovered at each SCF iteration (or time-step).

### Advanced settings

* The cavity can be improved/refined in several ways.
  * You can increase the number of points in the cavity by changing the default value of the variable {{< variable "PCMTessSubdivider" >}} from 1 to 4. In practice, the tessellation of each sphere centered at each atomic position (but Hydrogen's) from 60 to 240 points. Check the convergence and the final value of the total energy with respect to different tessellations.
  * You can also relax the constrain on the minimum distance {{< variable "PCMTessMinDistance" >}} (=0.1 A by default) between tesserae so as to obtain a more smooth cavity. This might be useful when the geometry of the molecule is complicated. Check if and how the total energy changes.
  * You can construct the molecular cavity by putting spheres also in Hydrogen, i.e., by setting {{< variable "PCMSpheresOnH" >}} = yes. 
  * You can change the radius of the spheres used to build up the cavity by manipulating {{< variable "PCMRadiusScaling" >}}. 
  * Moreover, by using {{< variable "PCMCavity" >}} = 'full path to cavity file' keyword the cavity geometry can be read from a file (instead of being generated inside Octopus, which is the default).

## Time propagation in a solvent within the Polarizable Continuum Model

In this tutorial we show how perform a time propagation of the Methene molecule in a solvent by setting up a time-dependent calculation and using the Time-dependent Polarizable Continnum Model (TD-PCM)

As in the case of gas-phase calculations, you should perform always the gs calculation before the td one. At the time being PCM is only implemented for TheoryLevel = DFT.

{{< notice note >}}
Carlo: I put the different methods into different sections. Used the same naming as in the octopus paper. Pleas add and comment inputs and outputs.
{{< /notice >}}

There are several ways to set up a PCM td calculation in Octopus, corresponding to different kinds of approximation for the coupled evolution of the solute-solvent system.

## Equilibrium TD-PCM

The first and most elementary approximation would be to consider an evolution where the solvent is able to instantaneously equilibrate with the solute density fluctuations. Formally within PCM, this means that the dielectric function is constant and equal to its static (zero-frequency) value for all frequencies $\epsilon(\omega)=\epsilon_0$. In practice, this approximation requires minimal changes to the gs input file: just replacing gs by td in the CalculationMode leaving the other PCM related variables unchanged. <span style="color: red>Carlo: Gabriel, please check that this is corresponds with what you mean</span>

### Input

{{< code-block >}}
 {{< variable "CalculationMode" >}} = td
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom

 {{< variable "PseudopotentialSet" >}} = hgh_lda
 
 FH = 0.917*angstrom
 %{{< variable "Coordinates" >}}
  "H" |   0 | 0 | 0
  "F" |  FH | 0 | 0
 %

 {{< variable "PCMCalculation" >}} = yes
 {{< variable "PCMStaticEpsilon" >}} = 78.39
{{< /code-block >}}

### Output

## Inertial TD-PCM

The first nonequilibrium approximation is to consider that there are fast and slow degrees of freedom of the solvent, of which only the former are able to equilibrate in real time with the solute. The slow degrees of freedom of the solvent that are frozen in time, in equilibrium with the initial (gs) state of the propagation. Formally, it means to consider two dielectric functions for the zero-frequency and high-frequency regime (namely, static and dynamic dielectric constants). In practice, this approximation is finally activated by setting extra input keywords, i.e., {{< variable "PCMTDLevel" >}} = neq and {{< variable "PCMDynamicEpsilon" >}}, the latter to the dynamic dielectric constant value.

### Input

### Output

## Equation of motion PCM

Finally, in the last --and most accurate-- nonequilibrium approximation there is an equation of motion for the evolution of the solvent polarization charges, rendering it history-dependent. Formally, this requires to use a specific model for the frequency-dependence for the dielectric function.

### Input

### Output

---
title: "Berry"
tags:: ["Tutorial", "Advanced", "Bulk", "Berry phase"]
#series: "Tutorial"
difficulties: "advanced"
system_types: "bulk"
features: "Berry phase"
hidden: True
---


In case of a finite electric field applied to a periodic system, the Octopus code can compute the effect of a finite electric field of finite, using techniques from the Modern Theory of Polarization and the notion of Berry phases.
This tutorial gives more details on what is currently implemented and how to use it.

At the moment, the code can only be used for calculations performed with a single k-point, and is therefore valid in the limit of large supercells.

/!\ This tutorial is still under construction. /!\







---------------------------------------------
---
title: "Geometry optimization"
tags: ["Tutorial", "Beginner", "Geometry Optimization", "Molecule", "Pseudopotentials", "DFT", "Total Energy", "Forces"]
#series: "Tutorial"
difficulties: "beginner"
theories: "DFT"
calculation_modes: "Geometry optimization"
species_types: "Pseudopotentials" 
system_types: "Molecule"
features: ["Total energy", "Forces"]
description: "Determining the geometry of methane."
---


In this tutorial we will show how to perform geometry optimizations using {{< octopus >}}.

## Methane molecule 

###  Manual optimization  
We will start with the methane molecule. Before using any sophisticated algorithm for the geometry optimization, it is quite instructive to do it first "by hand". Of course, this is only practical for systems with very few degrees of freedom, but this is the case for methane, as we only need to optimize the CH bond-length. We will use the following files.

####  {{< file inp >}} file

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 3.5*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom
 
 CH = 1.2*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3) 
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3) 
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
{{< /code-block >}}

This input file is similar to ones used in other tutorials. We have seen in the {{< tutorial "Basics:Total energy convergence" "Total_energy_convergence" >}} tutorial that for these values of the spacing and radius the total energy was converged to within 0.1 eV.

####  {{< file ch.sh >}} script  
You can run the code by hand, changing the bond-length for each calculation, but you can also use the following script:

```bash
#!/bin/bash
echo "-CH distance   Total Energy" > ch.log
list="0.90 0.95 1.00 1.05 1.10 1.15 1.20 1.25 1.30"
for CH in $list
do
    sed -ie "/CH = /s|[0-9.]\+|$CH|" inp
    octopus >& out-$CH
    energy=`grep Total static/info  | head -1 | cut -d "=" -f 2`
    echo $CH $energy >> ch.log
    rm -rf restart
done
```

{{< figure src="/images/Methane_Energy_vs_CH.png" width="500px" caption="Energy versus CH bond length for methane" >}}

####  Results  

Once you run the calculations, you should obtain values very similar to the following ones:

{{< code-block >}}
#CH distance   Total Energy
0.90 -215.22444186
0.95 -217.00785398
1.00 -218.09421507
1.05 -218.64692827
1.10 -218.78607337
1.15 -218.61468911
1.20 -218.20448384
1.25 -217.61525239
1.30 -216.89630205
{{< /code-block >}}

You can see on the right how this curve looks like when plotted. By simple inspection we realize that the equilibrium CH distance that we obtain is around 1.1 Å. This should be compared to the experimental value of 1.094 Å. From this data can you calculate the frequency of the vibrational breathing mode of methane? 

### FIRE

Since changing the bond-lengths by hand is not exactly the simplest method to perform the geometry optimization, we will now use automated algorithms for doing precisely that. Several of these algorithms are available in {{< octopus >}}. The default and recommend one is the the FIRE method[^footnote-1]
. 
{{< octopus >}} can also use the [GSL optimization routines](https://www.gnu.org/software/gsl/manual/html_node/Multidimensional-Minimization.html), so all the methods provided by this library are available.

#### Input

First we will run a simple geometry optimization using the following input file:

{{< code-block >}}
 {{< variable "FromScratch" >}} = yes
 {{< variable "CalculationMode" >}} = go
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "BoxShape" >}} = sphere
 {{< variable "Radius" >}} = 5.0*angstrom
 {{< variable "Spacing" >}} = 0.18*angstrom
 
 CH = 1.2*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0 
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3) 
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3) 
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
{{< /code-block >}}

In this case we are using default values for all the input variables related to the geometry optimization. Compared to the input file above, we have changed the {{< variable "CalculationMode" >}} to {{< code go >}}. Because the default minimum box of spheres centered on the atoms has some caveats when performing a geometry optimization (see bellow for more details), we also changed the box shape to a single sphere with a 5 Å radius. Like always, convergence of the final results should be checked with respect to the mesh parameters. In this case you can trust us that a radius of 5 Å is sufficient.

The geometry in the input file will be used as a starting point. Note that although a single variable CH was used to specify the structure, the coordinates will be optimized independently.

#### Output

When running {{< octopus >}} in geometry optimization mode the output will be a bit different from the normal ground-state mode. Since each minimization step requires several SCF cycles, only one line of information for each SCF iteration will printed: 

{{< code-block >}}
...
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    1 : etot -2.28028257E+02 : abs_dens 2.03E+00 : force  1.51E+00 : etime     1.1s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    2 : etot -2.28357657E+02 : abs_dens 2.31E-01 : force  4.90E-01 : etime     1.0s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    3 : etot -2.20168027E+02 : abs_dens 2.44E-01 : force  1.94E-02 : etime     1.0s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    4 : etot -2.18151090E+02 : abs_dens 6.26E-02 : force  2.02E-01 : etime     0.9s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    5 : etot -2.18383137E+02 : abs_dens 5.06E-03 : force  2.01E-02 : etime     0.9s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    6 : etot -2.18212576E+02 : abs_dens 3.38E-03 : force  1.24E-02 : etime     0.8s
ETA: .......1......2.......3......4......5.......6......7.......8......9......0

iter    7 : etot -2.18212135E+02 : abs_dens 3.30E-05 : force  1.81E-07 : etime     0.5s

             Info: Writing states. 2018/08/06 at 11:01:41


        Info: Finished writing states. 2018/08/06 at 11:01:41

Info: SCF converged in    7 iterations



++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    1 ++++++++++++++++++++++
  Energy    = -218.2121352982 eV
  Max force =    2.5241058228 eV/A
  Max dr    =    0.0002435392 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

...
{{< /code-block >}}

At the end of each minimization step the code will print the total energy, the maximum force acting on the ions and the maximum displacement of the ions from the previous geometry to the next one. {{< octopus >}} will also print the current geometry to a file named {{< file "go.XXXX.xyz" >}} (XXXX is the minimization iteration number) in a directory named {{< file "geom" >}}. In the working directory there is also a file named {{< file "work-geom.xyz" >}} that accumulates the geometries of each SCF calculation. The results are also summarized in {{< file "geom/optimization.log" >}}. The most recent individual geometry is in {{< file "last.xyz" >}}.

After some minimization steps the code will exit. 

{{< code-block >}}
...

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:   32 ++++++++++++++++++++++
  Energy    = -218.7910677857 eV
  Max force =    0.0513974737 eV/A
  Max dr    =    0.0052755032 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



Writing final coordinates to min.xyz
Info: Finished writing information to 'restart/gs'.

...
{{< /code-block >}}

If the minimization was successful a file named {{< file "min.xyz" >}} will be written in the working directory containing the minimized geometry. In this case it should look like this:

{{< code-block >}}
   5
 units: A
      C                   -0.000000   -0.000000   -0.000000
      H                    0.633663    0.633663    0.633663
      H                   -0.633663   -0.633663    0.633663
      H                    0.633663   -0.633663   -0.633663
      H                   -0.633663    0.633663   -0.633663
{{< /code-block >}}

Has the symmetry been retained? What is the equilibrium bond length determined? Check the forces in {{< file "static/info" >}}. You can also visualize the optimization process as an animation by loading {{< file "work-geom.xyz" >}} as an XYZ file into XCrySDen.

## Sodium trimer

Let's now try something different: a small sodium cluster. This time we will use an alternate method related to conjugate gradients. 

#### Input

Here is the input file to use

{{< code-block >}}
 {{< variable "FromScratch" >}} = yes
 {{< variable "CalculationMode" >}} = go
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 9.0*angstrom
 {{< variable "Spacing" >}} = 0.3*angstrom
 
 {{< variable "XYZCoordinates" >}} = "inp.xyz"
 
 {{< variable "GOMethod" >}} = cg_bfgs
{{< /code-block >}}

and the corresponding {{< file "inp.xyz" >}} file

{{< code-block >}}
  3

  Na     -1.65104552837    -1.51912734276    -1.77061890357
  Na     -0.17818042320    -0.81880768425    -2.09945089196
  Na      1.36469690888     0.914070162416    0.40877139068
{{< /code-block >}}

In this case the starting geometry is a random geometry.

#### Output
After some time the code should stop successfully. The information about the last minimization iteration should look like this:

{{< code-block >}}
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    3 ++++++++++++++++++++++
  Energy    =  -16.9047778355 eV
  Max force =    0.0257973692 eV/A
  Max dr    =    0.1384047905 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
{{< /code-block >}}

You might have noticed that in this case the code might perform several SCF cycles for each minimization iteration. This is normal, as the minimization algorithm used here requires the evaluation of the energy and forces at several geometries per minimization iteration. The minimized geometry can be found in the {{< file "min.xyz" >}} file:

{{< code-block >}}
   3
 units: A
      Na                  -2.269359   -1.808026   -1.617806
      Na                   0.457604   -0.515049   -2.237456
      Na                   1.347047    0.899166    0.394155
{{< /code-block >}}

#### Caveats of using a minimum box

For this example we are using the minimum box ({{< code-inline >}}{{< variable "BoxShape" >}} = minimum{{< /code-inline >}}). This box shape is the one used by default, as it is usually the most efficient one, but it requires some care when using it for a geometry optimization. The reason is that the box is generated at the beginning of the calculation and is not updated when the atoms move. This means that at the end of the optimization, the centers of the spheres composing the box will not coincide anymore with the position of the atoms. This is a problem if the atoms move enough such that their initial and final positions differ significantly. If this is the case, the best is to restart the calculation using the previously minimized geometry as the starting geometry, so that the box is constructed again consistently with the atoms. In this example, this can be done by copying the contents of the {{< file "min.xyz" >}} file {{< file "inp.xyz" >}} and running {{< octopus >}} again without changing input file.

This time the calculation should converge after only one iteration:

{{< code-block >}}
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
+++++++++++++++++++++ MINIMIZATION ITER -:    1 ++++++++++++++++++++++
  Energy    =  -16.8896899780 eV
  Max force =    0.0447462812 eV/A
  Max dr    =    0.0450282495 A
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
{{< /code-block >}}

You can now compare the old {{< file "min.xyz" >}} file with the new one:

{{< code-block >}}
   3
 units: A
      Na                  -2.252443   -1.788415   -1.588626
      Na                   0.484193   -0.498523   -2.232081
      Na                   1.302019    0.862970    0.362209
{{< /code-block >}}
Although the differences are not too large, there are also not negligible.

#### Further improving the geometry
Imagine that you now want to improve the geometry obtained and reduce the value of the forces acting on the ions. The stopping criterion is controlled by the {{< variable "GOTolerance" >}} variable (and {{< variable "GOMinimumMove" >}}). Its default value is 0.051 eV/Å (0.001 H/b), so let's set it to 0.01 eV/Å:

{{< code-block >}}
 {{< variable "FromScratch" >}} = yes
 {{< variable "CalculationMode" >}} = go
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 9.0*angstrom
 {{< variable "Spacing" >}} = 0.3*angstrom
 
 {{< variable "XYZCoordinates" >}} = "inp.xyz"
 
 {{< variable "GOMethod" >}} = cg_bfgs
 {{< variable "GOTolerance" >}} = 0.01*eV/angstrom
{{< /code-block >}}

Copy the contents of file {{< file "min.xyz" >}} to file {{< file "inp.xyz" >}} and rerun the code. This time, after some minimization steps {{< octopus >}} should return the following error message:

{{< code-block >}}
 **************************** FATAL ERROR *****************************
 *** Fatal Error (description follows)
 *--------------------------------------------------------------------
 * Error occurred during the GSL minimization procedure:
 * iteration is not making progress towards solution
 **********************************************************************
{{< /code-block >}}

This means the minimization algorithm was unable to improve the solutions up to the desired accuracy. Unfortunately there are not many things that can be done in this case. The most simple one is just to restart the calculation starting from the last geometry (use the last {{< file "go.XXXX.xyz" >}} file) and see if the code is able to improve further the geometry.











---------------------------------------------
[^footnote-1]: {{< article title="Structural relaxation made simple" authors="Erik Bitzek, Pekka Koskinen, Franz Gähler, Michael Moseler, and Peter Gumbsch" journal="Phys. Rev. Lett." volume="97" pages="170201" year="2006" url="https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201" doi="10.1103/PhysRevLett.97.170201" link="APS" >}}

---
title: "Unfolding"
#tags: ["Tutorial", "Expert", "Ground State", "Bulk", "oct-unfold"]
difficulties: "expert"
system_types: "bulk"
calculation_modes: "Ground state"
utilities: "oct_unfold"
description: "How to perform a band-structure unfolding."
---



In this tutorial, we look at how to perform a band-structure unfolding using Octopus. 

This calculation is done in few steps, which are described below. First, we need to compute the ground state of a supercell.

## Input
Here we start with the input file for Octopus. This is a ground state calculation for bulk silicon, similar to the one of the tutorial {{< tutorial "Periodic_systems" "Periodic Systems" >}}.
In the present case however, we will use a supercell of Si, composed of 8 atoms. 

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 
 {{< variable "PeriodicDimensions" >}} = 3
 BoxShape = parallelepiped
 {{< variable "Spacing" >}} = 0.5
 
 a = 10.18
 %{{< variable "LatticeParameters" >}}
  a | a  | a
 90 | 90 |90
 %
 
 %{{< variable "ReducedCoordinates" >}}
  "Si" | 0.0         | 0.0       | 0.0 
  "Si" |   1/2       | 1/2       | 0.0
  "Si" |   1/2       | 0.0       | 1/2
  "Si" |   0.0       | 1/2       | 1/2
  "Si" |   1/4       | 1/4       | 1/4
  "Si" |   1/4 + 1/2 | 1/4 + 1/2 | 1/4
  "Si" |   1/4 + 1/2 | 1/4       | 1/4 + 1/2
  "Si" |   1/4       | 1/4 + 1/2 | 1/4 + 1/2
 %
 
 nk = 4
 %{{< variable "KPointsGrid" >}}
   nk |  nk |  nk
 %
 {{< variable "KPointsUseSymmetries" >}} = yes
 
 {{< variable "ConvRelDens" >}} = 1e-7
 {{< variable "EigensolverTolerance" >}} = 1e-8
{{< /code-block >}}


Most of these variables are already used in other tutorials, so the interested read could refer for instance to {{< tutorial "Periodic_systems" "Periodic Systems" >}} for a detailed description.

## Unfolding - step 1 

After running Octopus using the previous input file, we have the obtained the ground state of the supercell.
We now what to define the primitive cell on which we want to unfold or supercell, and the specific k-point path that we are interested in.

This is done by adding to the input file the following lines:

{{< code-block >}}
  {{< variable "UnfoldMode" >}} = unfold_setup
  %{{< variable "UnfoldLatticeParameters" >}}
    a | a | a
  %
  %{{< variable "UnfoldLatticeVectors" >}}
   0.  | 0.5 | 0.5
   0.5 | 0.  | 0.5
   0.5 | 0.5 | 0.0  
  %
  %{{< variable "UnfoldKPointsPath" >}}
   4 | 4 | 8
   0.5 | 0.0 | 0.0 - L point
   0.0 | 0.0 | 0.0 - Gamma point
   0.0 | 0.5 | 0.5 - X point
   1.0 | 1.0 | 1.0 - Another Gamma point
  %
{{< /code-block >}}

Lets see more in detail some of the input variables:

* <tt>{{< variable "UnfoldMode" >}} = unfold_setup</tt>: this variable instruct the utility <tt>oct-unfold</tt> in which mode we are running. As a first step, we are runing in the <tt>unfold_setup</tt> mode, which generates necessary files for getting the unfolded band structure.

* <tt>{{< variable "UnfoldLatticeParameters" >}}</tt> specifies the lattice parameters of the primitive cell. This variable is similar to {{< variable "LatticeParameters" >}}.

* <tt>{{< variable "UnfoldLatticeVectors" >}}</tt> specifies the lattice vectors of the primitive cell. This variable is similar to {{< variable "LatticeVectors" >}}.

* <tt>{{< variable "UnfoldKPointsPath" >}}</tt> specifies the k-point path for the unfolded. The coordinates are indicated as reduced coordinates for the primitive lattice. This variable is similar to {{< variable "KPointsPath" >}}.

By running <tt>oct-unfold</tt>, you will obtain two files. The first one (<tt>unfold_kpt.dat</tt>) contains the list of k-points of the specified k-point path, but expressed in reduced coordinates of the supercell. These are the k-points for which we need to evaluate the wavefunctions of the supercell.
The second file (<tt>unfold_gvec.dat</tt>) contains the list of reciprocal lattice vectors that relate the k-points of the primitive cell to the one in the supercell.

The unfolding still been an experimental feature, you also need to add to your input file the line
{{< code-block >}}
  ExperimentalFeatures = yes
{{< /code-block >}}

## Unfolding - step 2

In order to perform the unfolding, we must compute wavefunctions in the supercell at specific k-points. These points are described in the file <tt>unfold_kpt.dat</tt>.
To get them, we need to run an non self-consistent calculation using
{{< code-block >}}
  {{< variable "CalculationMode" >}} = unocc
{{< /code-block >}}

In order to use the list of k-points obtained in step 1, we remove from the input file the k-point related variables (%{{< variable "KPointsGrid" >}}, %{{< variable "KPointsPath" >}}, or %{{< variable "KPoints" >}}) and we replace them by the line
{{< code-block >}}
  include unfold_kpt.dat
{{< /code-block >}}

We then run {{< octopus >}} on this input file. 
Note that in order to get extra states, one should use the variables {{< variable "ExtraStates" >}} and {{< variable "ExtraStatesToConverge" >}}, as for a normal band structure calculation, 
see {{< tutorial "Periodic_systems" "Periodic Systems" >}}.

For this tutorial, we used

{{< code-block >}}
  {{< variable "ExtraStates" >}} = 6
  {{< variable "ExtraStatesToConverge" >}} = 4
  MaximumIter = 100
{{< /code-block >}}

Note that after you performed this step, you cannot run a TD calculation in this folder.
Indeed, the original GS states have been replaced by the one of this unocc calculation.

Here, we are requesting to converge 4 unoccupied states, which all correspond to the first valence band in the folded primitive cell, as the supercell is 4 times larger than the initial cell.

## Unfolding - step 3

Now that we have computed the states for the required k-points, we can finally compute the spectral function for each of the k-points of the specified k-point path.
This is done by changing unfolding mode to be
{{< code-block >}}
  {{< variable "UnfoldMode" >}} = unfold_run
{{< /code-block >}}

This will produce the spectral function for the full path (<tt>static/ake.dat</tt>) and for each individual points of the path (<tt>static/ake_XXX.dat</tt>). 
The unfolded bandstructure of silicon is shown below.

{{< figure src="/images/Si_unfolded.png" width="500px" caption="Unfolded band structure of bulk silicon." >}}

Note that in order to control the energy range and the energy resolution for the unfolded band structure, one can specify the variables {{< variable "UnfoldMinEnergy" >}}, {{< variable "UnfoldMaxEnergy" >}}, and {{< variable "UnfoldEnergyStep" >}}.

It is important to note here that the code needs to read all the unoccupied wavefunctions, and therefore will need to have enough memory.







---
title: "RDMFT"
tags: ["Tutorial", "Advanced", "RDMFT"]
difficulties: "advanced"
theories: ["independent-particles", "rdmft"]
calculation_modes: ["Ground state"]
system_types: "Molecule"
species_types: "Pseudopotentials"
description: "How to do a Reduced Density Matrix Functional Theory (RDMFT) calculation."

#series: "Tutorial"
---


In this tutorial, you will learn how to do a Reduced Density Matrix Functional Theory (RDMFT) calculation with {{< octopus >}}. 

In contrast to density functional theory or Hartree-Fock theory, here we do **not** try to find one optimal slater-determinant (single-reference) to describe the ground-state of a given many-body system, but instead approximate the one-body reduced density matrix (1RDM) of the system. Thus, the outcome of an RDMFT minimization is a set of eigen-vectors, so called ''natural orbitals'' (NOs,) and eigenvalues, so-called ''natural occupation numbers'' (NONs), of the 1RDM. One important aspect of this is that we need to use ''more'' orbitals than the number of electrons. This additional freedom allows to include static correlation in the description of the systems and hence, we can describe settings or processes that are very difficult for single-reference methods. One such process is the '''dissociation of a molecule''', which we utilize as example for this tutorial: You will learn how to do a properly converged dissociation study of H$_2$. For further information about RDMFT, we can recommend e.g. chapter ''Reduced Density Matrix Functional Theory (RDMFT) and Linear Response Time-Dependent RDMFT (TD-RDMFT)'' in the book of Ferre, N., Filatov, M. & Huix-Rotllant, M. ''Density-Functional Methods for Excited States'' (Springer, 2015). doi:10.1007/978-3-319-22081-9 

## Basic Steps of an RDMFT Calculation

The first thing to notice for RDMFT calculations is that, we will explicitly make use of a ''basis set'', which is in contrast to the minimization over the full real-space grid that is performed normally with {{< octopus >}}. Consequently, we always need to do a preliminary calculation to generate our basis set. We will exemplify this for the hydrogen molecule in equilibrium position. Create the following {{< file "inp" >}} file:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "TheoryLevel" >}} = independent_particles
 
 {{< variable "Dimensions" >}} = 3
 {{< variable "Radius" >}} = 8
 {{< variable "Spacing" >}} = 0.15
 
 # distance between H atoms
 d = 1.4172 - (equilibrium bond length H2)
 
 %{{< variable "Coordinates" >}}
  'H' | 0 | 0 | -d/2
  'H' | 0 | 0 | d/2
 %
 
 {{< variable "ExtraStates" >}} = 14
{{< /code-block >}}

You should be already familiar with all these variables, so there is no need to explain them again, but we want to stress the two important points here: 
- We chose {{< variable "TheoryLevel" >}} = independent_particles, which you should always use as a basis in an RDMFT calculation (see actual octopus paper?)
- We set the {{< variable "ExtraStates" >}} variable, which controls the size of the basis set that will be used in the RDMFT calculation. Thus in RDMFT, we have with the ''basis set'' a ''new numerical parameter'' besides the {{< variable "Radius" >}} and the {{< variable "Spacing" >}} that needs to be converged for a successful calculation. The size of the basis set M is just the number of the orbitals for the ground state (number of electrons divided by 2) plus the number of extra states.
- All the numerical parameters ''depend on each other'': If we want to include many {{< variable "ExtraStates" >}}  to have a sufficiently large basis, we will need also a larger {{< variable "Radius" >}} and especially a smaller {{< variable "Spacing" >}} at a certain point. The reason is that the additional states will have function values bigger than zero in a larger region than the bound states. Additionally, all states are orthogonal to each other, and thus the grid needs to resolve more nodes with an increasing number of {{< variable "ExtraStates" >}}.

With this basis, we can now do our first rdmft calculation. For that, you just need to change the {{< variable "TheoryLevel" >}} to 'rdmft' and ''add the part {{< variable "ExperimentalFeatures" >}} = yes''. In the standard out there will be some new information:
* Calculating Coulomb and exchange matrix elements in basis --this may take a while--
* Occupation numbers: they are not zero or one/two here but instead fractional!
* ''anything else important?''

## Basis Set Convergence

Now, we need to converge our calculation for the three numerical parameters. This needs to be done in some iterative way, which means we first set {{< variable "ExtraStates" >}} and then do convergence series for {{< variable "Radius" >}} and {{< variable "Spacing" >}} as we know it. Then, we utilize the converged values and perform a series for {{< variable "ExtraStates" >}}. Having converged {{< variable "ExtraStates" >}}, we start again with the series for {{< variable "Radius" >}} and {{< variable "Spacing" >}} and see if the energy (or better the electron density) still changes considerably. If yes, this meas that the higher lying {{< variable "ExtraStates" >}} where net well captured by the simulation box and we need to adjust the parameters again. We go on like this until convergence of all parameters. This sound like a lot of work, but don't worry, if you have done such convergence studies a few times, you will get a feeling for good values.

So let us do one iteration as example. We start with a guess for {{< variable "ExtraStates" >}}, which should not be too small, say 14 (thus $M=1+14=15$) and converge {{< variable "Radius" >}} and {{< variable "Spacing" >}}. To save time, we did this for you and found the values
{{< variable "Radius" >}} = 8
{{< variable "Spacing" >}} = 0.15.
If you feel insecure with such convergence studies, you should have a look in the corresponding tutorial: [[Tutorial:Total_energy_convergence]]

Now let us perform a basis-set series. For that, we first create reference basis by repeating the above calculation with say {{< variable "ExtraStates" >}} = 29. We copy the restart folder in a new folder that is called 'basis' and execute the following script (if you rename anything you need to change the respective part in the script):
```bash
 #!/bin/bash
 
 series=ES-series
 outfile="$series.log"
 echo "#ES    Energy    1. NON  2.NON" > $outfile
 list="4 9 14 19 24 29"
 
 export OCT_PARSE_ENV=1 
 for param in $list 
 do
   folder="$series-$param"
   mkdir $folder
   out_tmp=out-RDMFT-$param
   cd $folder
    # here we specify the basis folder
    cp -r ../basis/restart .
    # create inp file
    {
 cat <<-EOF
 calculationMode = gs
 TheoryLevel = rdmft
 ExperimentalFeatures = yes
 
 ExtraStates = $param
 
 Dimensions = 3
 
 Radius = 8
 Spacing = 0.15
 
 # distance between H atoms
 d = 1.4172 - (equilibrium) 
 
 %Coordinates
   'H' | 0 | 0 | -d/2
   'H' | 0 | 0 | d/2
 %
 
 Output = density
 OutputFormat = axis_x
 EOF
 } > inp
  
   # put the octopus dir here
   [octopus-dir]/bin/octopus > $out_tmp
 
    energy=`grep -a "Total energy" $out_tmp  | tail -1 | awk '{print $3}'`
    seigen=`grep -a " 1  " static/info |  awk '{print $2}'`
    peigen=`grep -a " 2  " static/info |  awk '{print $2}'`
    echo $param $energy $seigen $peigen >> ../$outfile
  cd ..
 done
```

If everything works out fine, you should find the following output in the 'ES-series.log' file:

{{< code-block >}}
 #ES    Energy    1. NON  2.NON
 4 -1.1476150018E+00 1.935750757008 0.032396191164
 9 -1.1498006205E+00 1.932619193929 0.032218954381
 14 -1.1609241676E+00 1.935215440985 0.032526426664
 19 -1.1610006378E+00 1.934832116929 0.032587169713
 24 -1.1622104536E+00 1.932699204653 0.032997371081
 29 -1.1630839347E+00 1.932088486112 0.032929702131
{{< /code-block >}}

So we see that even with M=30, the calculation is not entirely converged! Thus for a correctly converged result, one would need to further increase the basis and optimize the simulation box. Since the calculations become quite expensive, we would need to do them on a high-performance cluster, which goes beyond the scope of this tutorial. We will thus have to proceed with '''not entirely converged''' parameters in the next section.

## H2 Dissociation

In this third and last part of the tutorial, we come back to our original goal: the calculation of the H2 dissociation curve. Now, we want to change the ''d''-parameter in the input file and consequently, we should use an adapted simulation box. Most suited for a dissociation, is the {{< code-inline >}}{{< variable "BoxShape" >}}= cylinder{{< /code-inline >}}. However, since this makes the calculations again much more expensive, we just stick to the default choice of {{< code-inline >}}{{< variable "BoxShape" >}}=minimal {{< /code-inline >}}.

We choose {{< code-inline >}}{{< variable "ExtraStates" >}}=14{{< /code-inline >}} as a compromise between accuracy and numerical cost and do the d-series with the following script:

```bash
 #!/bin/bash
 
 series=dist-series
 outfile="$series.log"
 echo "#ES    Energy    1. NON  2.NON" > $outfile
 list="0.25 0.5 0.75 1.0 1.5 2.0 2.5 3.0 4.0 5.0 6.0 8.0"
 
 export OCT_PARSE_ENV=1 
 for param in $list 
 do
  folder="$series-$param"
  mkdir $folder
  out_tmp=out-$series-$param
  cd $folder
 
  # create inp file
 {
 cat <<-EOF
 calculationMode = gs
 TheoryLevel = rdmft
 ExperimentalFeatures = yes
 
 ExtraStates = 14
 
 Dimensions = 3
 
 Radius = 8
 Spacing = 0.15
 
 # distance between H atoms
 d = $param
 
 %Coordinates
   'H' | 0 | 0 | -d/2
   'H' | 0 | 0 | d/2
 %
 
 Output = density
 OutputFormat = axis_x
 EOF
 } > inp
  
  # put the octopus dir here
  oct_dir=~/octopus-developer/bin  
  
  # generate the basis
  #export OCT_TheoryLevel=independent_particles
  #${oct_dir}/octopus > ${out_tmp}_basis
  
  # do the actual rdmft calculation
  #export OCT_TheoryLevel=rdmft
  #${oct_dir}/octopus > ${out_tmp}
 
  energy=`grep -a "Total energy" $out_tmp  | tail -1 | awk '{print $3}'`
  non1=`grep -a " 1  " static/info |  awk '{print $2}'`
  non2=`grep -a " 2  " static/info |  awk '{print $2}'`
  echo $param $energy $non1 $non2 >> ../$outfile
  cd ..
 done
```

The execution of the script should take about 30-60 minutes on a normal machine, so now it's time for well-deserved coffee!

After the script has finished, you have find a new log file and if you plot the energy of the files, you should get something of the following shape:
[[File:H2 diss curve.png|thumb|Dissociation curve of the Hydrogen molecule, calculated with RDMFT using the Müller functional.]]

Comments:
* This dissociation limit is not correct. Since the basis set is not converged, we get a too large value of E_diss=-1.009 than the Müller functional obtains for a converged set (E_diss=-1.047 Hartree), which is paradoxically closer to the exact dissociation limit of E_diss=1.0 Hartree because of an error cancellation (too small bases always increase the energy).
* The principal shape is good. Especially, we see that the curve becomes constant for large d, which single-reference methods cannot reproduce. 





---------------------------------------------
---
title: "Wannier90"
#tags: ["Tutorial", "Expert", "Ground State", "Bulk", "Metals", "oct-wannier90"]
difficulties: "expert"
system_types: ["bulk", "metal"]
calculation_modes: "Ground state"
utilities: "oct-wannier90"
description: "Using the Wannier90 interface."
---


Starting from Octopus Selene (10.x), {{< octopus >}} provides an interface for Wannier90. In order to perform this tutorial, you need first to compile Wannier90. Detailed information can be found on the [Wannier90 website](https://www.wannier.org/) .

In this tutorial, we follow the second example from the Wannier90 tutorials, which consists in computing the Fermi surface of lead. This requires several steps, which are explained in detail bellow.


## Ground-state

### Input
We start with the input file for {{< octopus >}}. At this point there is nothing special to specify as {{< octopus >}} does not need to know that Wannier90 will be called later.

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "PeriodicDimensions" >}} = 3
 
 {{< variable "Spacing" >}} = 0.4
 
 %{{< variable "LatticeVectors" >}}
   0.0 | 0.5 | 0.5 
   0.5 | 0.0 | 0.5
   0.5 | 0.5 | 0.0
 %
 
 a = 9.3555
 %{{< variable "LatticeParameters" >}}
  a | a | a
 %
 
 {{< variable "PseudopotentialSet" >}} = hgh_lda
 %{{< variable "ReducedCoordinates" >}}
  "Pb" | 0.0 | 0.0 | 0.0 
 %
 
 %{{< variable "KPointsGrid" >}}
   4 |  4 |  4
 %
 
 {{< variable "Smearing" >}} = 0.1*eV
 {{< variable "SmearingFunction" >}} = fermi_dirac
 {{< variable "ExtraStates" >}} = 4
{{< /code-block >}}

The use of most of these variables is already explained in detail in other tutorials. See, for example, the [Periodic Systems](../Periodic_systems) tutorial. Nevertheless, some options deserve a more detailed explanation:
* <tt>{{< variable "SmearingFunction" >}} = fermi_dirac</tt>: as we are interested in lead, which is a metal, we need to use a smearing of the occupations in order to get a converged ground-state. Here we have chosen the smearing function to be a Fermi-Dirac distribution.
* <tt>{{< variable "Smearing" >}} = 0.1*eV</tt>: this defines the effective temperature of the smearing function.
* <tt>{{< variable "PseudopotentialSet" >}} = hgh_lda</tt>: since lead is not provided with the default pseudopotential set, we must select one containing Pb. We thus use the LDA HGH pseudopotential set here.

It is also important to note that the ''k''-point symmetries cannot be used with Wannier90. Also, only Monkhorst-Pack grids are allowed, so that only one ''k''-point shift is allowed.

### Output
Now run {{< octopus >}} using the above input file. The output should be similar to other ground-state calculations for periodic systems, so we are not going to look at it in detail. Nevertheless, to make sure there was no problem with the calculation and since we will need this information later, you should check the calculated Fermi energy. You can find this information in the {{< file "static/info" >}} file:
{{< code-block >}}

Fermi energy =     0.157425 H
{{< /code-block >}}


## Wannier90

The usage of Wannier90 as a post-processing tool requires several steps and several invocations of both Wannier90 and {{< octopus >}}. Note that all the interfacing between Wannier90 and {{< octopus >}} is done through a specialized utility called {{< file "oct-wannier90" >}}.

### Step one: create the Wannier90 input file

In order to run Wannier90, you must prepare an input file for it. To make this task easier, the {{< file "oct-wannier90" >}} utility can take the {{< octopus >}} input file and "translate" it for Wannier90. Here is how the input file for the utility should look like:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "PeriodicDimensions" >}} = 3
 
 {{< variable "Spacing" >}} = 0.4
 
 %{{< variable "LatticeVectors" >}}
   0.0 | 0.5 | 0.5 
   0.5 | 0.0 | 0.5
   0.5 | 0.5 | 0.0
 %
 
 a = 9.3555
 %{{< variable "LatticeParameters" >}}
  a | a | a
 %
 
 {{< variable "PseudopotentialSet" >}} = hgh_lda
 %{{< variable "ReducedCoordinates" >}}
  "Pb" | 0.0 | 0.0 | 0.0 
 %
 
 %{{< variable "KPointsGrid" >}}
   4 |  4 |  4
 %
 
 {{< variable "Smearing" >}} = 0.1*eV
 {{< variable "SmearingFunction" >}} = fermi_dirac
 {{< variable "ExtraStates" >}} = 4
 
 {{< variable "Wannier90Mode" >}} = w90_setup
{{< /code-block >}}

This is basically the same input file as above with the following option added:
* <tt>{{< variable "Wannier90Mode" >}} = w90_setup</tt>: this tells the {{< file "oct-wannier90" >}} utility that we want to setup an input file for Wannier90 based on the contents of our input file.

Now run the the {{< file "oct-wannier90" >}} utility. This will create a file named {{< file "w90.win" >}} that should look like this:
{{< code-block >}}

- this file has been created by the Octopus wannier90 utility
 
begin unit_cell_cart
Ang
   0.00000000   2.47535869   2.47535869
   2.47535869   0.00000000   2.47535869
   2.47535869   2.47535869   0.00000000
end unit_cell_cart
 
begin atoms_frac
Pb     0.00000000   0.00000000   0.00000000
end atoms_frac
 
use_bloch_phases = .true.
 
num_bands    6
num_wann    6
 
write_u_matrices = .true.
translate_home_cell = .true.
write_xyz = .true.
 
mp_grid    4   4   4
 
begin kpoints
 -0.00000000  -0.00000000  -0.00000000
 -0.25000000  -0.00000000  -0.00000000
...
{{< /code-block >}}

Before running Wannier90, we are going to slightly tweak this file. As we are interested in lead, we will replace the following lines:
{{< code-block >}}

use_bloch_phases = .true.
 
num_bands    6
num_wann    6
{{< /code-block >}}

by
{{< code-block >}}

num_bands    6
num_wann    4
begin projections
Pb:sp3
end projections
{{< /code-block >}}

By doing this, we are telling Wannier90 that we want to find 4 Wannier states out of the 6 bands contained in the DFT calculation. Moreover, we are using guess projections based on sp3 orbitals attached to the lead atom. More details about Wannier90, projections, and disentanglement can be found in the tutorials and manual of Wannier90.

### Step 2: generate the {{< file ".nnkp" >}} file

After editing the file {{< file "w90.win" >}}, run Wannier90 in post-processing mode:
{{< command-line "wannier90.x -pp w90.win" >}}

This will produce a file named {{< file "w90.nnkp" >}}. This file contains information telling {{< octopus >}} what further information does Wannier90 require for the following steps.

### Step 3: generate files for Wannier90

At this stage, we are going to tell {{< octopus >}} to generate some files for Wannier90 using the information from the {{< file "w90.nnkp" >}} file. To do so, we change the value of the {{< variable "Wannier90Mode" >}} variable in the above {{< octopus >}} input file to:
{{< code-block >}}
 {{< variable "Wannier90Mode" >}} = w90_output
{{< /code-block >}}
Then, run the {{< file "oct-wannier90" >}} utility again. This time the utility will produces the files {{< file "w90.eig" >}}, {{< file "w90.mmn" >}}, and {{< file "w90.amn" >}}. The variable that controls which files are produced is {{< variable "Wannier90Files" >}}.

### Step 4: run Wannier90

Finally, you call Wannier90 to perform the Wannierization:
{{< command-line "wannier90.x w90.win" >}}

## Computing the Fermi surface

{{< figure src="/images/FS_Lead_band2_band3.png" width="500px" caption="Fermi surface for bands 2 and 3 in lead." >}}

Once you performed all the previous steps, you know how to run Wannier90 on top of an {{< octopus >}} calculation.
Lets now compute the interpolated Fermi surface. To do that, we add the following lines to the {{< file "w90.win" >}} file:
{{< code-block >}}

restart = plot
fermi_energy = 4.283754645
fermi_surface_plot = true
{{< /code-block >}}

The value of the Fermi energy corresponds to the one obtained by {{< octopus >}}, converted to eV.

By calling again Wannier90, you will obtain a file {{< file "w90.bxsf" >}}, which you can plot using, for instance, xcrysden:
{{< command-line "xcrysden --bxsf w90.bxsf" >}}
From here you should be able to produce the plot of the Fermi surface for band 2 and band 3, as shown in the figure.

## Visualization of Wannier states using Octopus or Wannier90

Wannier90 allows one to compute the Wannier states. For this, you need to output the periodic part of the Bloch states on the real-space grid. This produces the {{< file ".unk" >}} files. For this, you need to specify the output of {{< file "oct-wannier90" >}} by adding the following line to the {{< file "inp" >}} file
{{< code-block >}}
 {{< variable "Wannier90Files" >}} = w90_unk
{{< /code-block >}}

Note that in order to get all the files you need in one run, you could instead use <tt>{{< variable "Wannier90Files" >}} = w90_amn + w90_mmn + w90_eig + w90_unk</tt>.

An alternative way of getting the states is to let {{< octopus >}} compute them. You need to get the <tt>U</tt> matrices from Wannier90 for this by adding the following line to the {{< file "w90.win" >}} file:
{{< code-block >}}
  write_u_matrices = .true.
{{< /code-block >}}

Then run the {{< file "oct-wannier90" >}} utility with the appropriate mode
{{< code-block >}}
  {{< variable "Wannier90Mode" >}} = w90_wannier
{{< /code-block >}}
This will create a folder named {{< file "wannier" >}} containing the calculated Wannier states. The main advantage of doing this is that the Wannier states can be obtained in any of the output formats supported by {{< octopus >}}. This is controlled, as usual, by the {{< variable "OutputFormat" >}} variable.








---
title: "Parallelization and performance"
tags: ["Tutorial", "Advanced"]
difficulties: "advanced"
#series: "Tutorial"
description: "Using Octopus on supercomputers."
---

{{< notice note >}}
This tutorial page was set up for the Benasque TDDFT school 2014. The specific references to the supercomputer used at that time will have to be adapted for others to use this tutorial.
{{< /notice >}}

In this tutorial we will see how to run a relatively big system in the [Hopper supercomputer](https://www.nersc.gov/users/computational-systems/hopper) (at NERSC in California), and how to measure its performance. There are a few key things you need to know about how to interact with the machine. To log in, run {{< code "ssh trainX@hopper.nersc.gov" >}} in your terminal, substituting the actual name of your training account. Be aware that since this machine is far away, you should not try running X-Windows programs! You submit jobs by the {{< code "qsub" >}} command, ''e.g.'' {{< code "qsub job.scr" >}}, which will put them in the queue for execution when there is free space. You can see what jobs you currently have in the queue by executing {{< code "qstat -u $USER" >}}, so you can see when your job finishes. A status code will be shown: Q = waiting in the queue, R = running, C = complete. You can cancel a job by {{< code "qdel" >}} + the job number, as written by {{< code "qstat" >}}. The job script (''e.g.'' {{< code "job.scr" >}}) specifies parameters to the PBS/Torque queuing system about how many cores to use, what commands to run, etc.

# Running the ground state

We will need the input file, job submission script, coordinates file, and a pseudopotential for Mg (for the other elements we will use the default ones that come with Octopus). The pseudopotential is available in the Octopus directory at {{< code "jube/input/Mg.fhi" >}}. You can copy the files on hopper directly from {{< code "/global/homes/d/dstrubbe/octopus_tutorial" >}} to your scratch directory as follows:

```bash
 cd $SCRATCH
 cp -r /global/homes/d/dstrubbe/octopus_tutorial .
```

The input file ({{< file inp >}}):
{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 
 #### System size and parameters
 
 {{< variable "Spacing" >}} = 0.20
 {{< variable "Radius" >}} = 4.0
 
 {{< variable "Units" >}} = ev_angstrom
 {{< variable "XYZCoordinates" >}} = "xyz"
 %{{< variable "Species" >}}
  "Mg" | 24.305 | spec_ps_fhi | 12 | 3 | 2
 %
 {{< variable "ExcessCharge" >}} = 0
 
 {{< variable "XCFunctional" >}} = gga_x_pbe + gga_c_pbe
 
 {{< variable "ExtraStates" >}} = 18
 {{< variable "Eigensolver" >}} = rmmdiis
 {{< variable "LCAOAlternative" >}} = yes
 {{< variable "SmearingFunction" >}} = fermi_dirac
 {{< variable "Smearing" >}} = 0.1
 {{< variable "Mixing" >}} = 0.15

 #### GS
 {{< variable "MaximumIter" >}} = 300
 {{< variable "EigensolverTolerance" >}} = 1e-8
 {{< variable "ConvRelDens" >}} = 5e-8
 
 #### Saving memory

 {{< variable "SymmetriesCompute" >}} = no
 {{< variable "PartitionPrint" >}} = no
 {{< variable "MeshPartitionPackage" >}} = metis
 
 # Additional options
 {{< variable "ExperimentalFeatures" >}} = yes
{{< /code-block >}}

Submission script {{< file job.scr >}} with 24 CPU processor cores: 
```bash
 #!/bin/bash
 #PBS -q regular
 #PBS -l mppwidth=24
 #PBS -l advres=benasque.348
 #PBS -l walltime=0:30:00
 #PBS -N testing_chl
 #PBS -V
 
 module load octopus/4.1.2
 cd $PBS_O_WORKDIR
 aprun -n 24 octopus_mpi &> output_gs_24
```

To run:
```bash
 qsub job.scr
```

Coordinates file {{< file xyz >}}. Take a look at it (on your local machine) with visualization software such as {{< code "xcrysden" >}} to see what kind of molecule we are dealing with.
{{% expand "Expand: xyz file" %}}
{{< code-block >}}
 148
 units: A
      C                  -10.053414   -2.206790   -1.658417
      C                   -9.302986   -0.895144   -1.860284
      O                   -8.207781   -0.752638   -1.308974
      H                  -10.995492   -2.174759   -2.206494
      H                  -10.255731   -2.349212   -0.596194
      H                   -9.446853   -3.033807   -2.027068
      N                   -9.793020    0.106019   -2.626896
      C                   -9.183715    1.417441   -2.564153
      H                  -10.613638    0.026372   -3.226328
      H                   -9.703953    2.093287   -3.243002
      H                   -8.135921    1.345699   -2.855849
      H                   -9.251150    1.801765   -1.546540
      Mg                  -7.589714   -0.406247    0.589093
      C                   -5.742517   -3.346473    0.598957
      C                   -4.800687    1.285199   -0.432991
      C                   -9.267821    2.551225    1.101757
      C                  -10.244120   -2.188907    2.138434
      N                   -5.505702   -0.986494    0.214199
      C                   -5.032092   -2.255827    0.120854
      C                   -3.658251   -2.282504   -0.557964
      C                   -3.334778   -0.768793   -0.749853
      C                   -4.616448   -0.057153   -0.306840
      C                   -2.157453   -0.329206    0.109120
      C                   -3.729999   -3.115235   -1.844687
      C                   -2.766872   -2.716411   -2.960238
      C                   -1.331215   -2.744537   -2.525131
      O                   -0.791040   -3.518163   -1.775651
      O                   -0.671785   -1.759611   -3.177730
      N                   -7.068068    1.624336    0.467327
      C                   -5.959997    2.072749   -0.116268
      C                   -6.119618    3.519474   -0.405906
      C                   -7.367643    3.886460    0.028358
      C                   -8.016810    2.676853    0.594467
      C                   -5.092331    4.342863   -1.061147
      C                   -7.969215    5.193689   -0.096865
      C                   -8.764346    5.768120    0.809053
      N                   -9.475374    0.106282    1.501210
      C                   -9.944387    1.380878    1.591879
      C                  -11.236629    1.354534    2.222881
      C                  -11.533164   -0.000127    2.475692
      C                  -10.409531   -0.760031    2.032458
      C                  -12.065705    2.482921    2.565847
      O                  -11.852095    3.614620    2.163993
      C                  -12.773981   -0.530453    3.073700
      C                  -13.971292   -0.223432    2.172826
      N                   -7.915114   -2.308227    1.405664
      C                   -9.126921   -2.901012    1.866373
      C                   -8.885865   -4.386110    1.948373
      C                   -7.612956   -4.594110    1.506344
      C                   -7.041245   -3.265509    1.152599
      C                   -9.873356   -5.349056    2.414656
      C                   -6.572822   -5.562522    1.244981
      O                   -6.595619   -6.753651    1.424254
      C                   -5.368183   -4.810959    0.605175
      C                   -4.115033   -5.056118    1.400777
      O                   -3.994084   -5.204225    2.587175
      O                   -3.071063   -5.083305    0.537626
      C                   -1.772257   -5.277910    1.082827
      C                    0.753800   -1.655738   -3.032461
      C                    1.062861   -0.433850   -2.240797
      C                    1.344919   -0.409138   -0.931698
      C                    1.342960   -1.602548   -0.036439
      C                    1.703393    0.883925   -0.254470
      C                    3.229775    1.030190   -0.125447
      C                    3.578204    2.129802    0.882753
      C                    5.078428    2.494313    0.882648
      C                    5.257557    3.825237    1.627432
      C                    5.895268    1.372849    1.548820
      C                    7.408164    1.625845    1.487162
      C                    8.167514    0.444435    2.106308
      C                    9.683652    0.673563    2.265625
      C                    9.978504    1.909940    3.121920
      C                   10.398360    0.775333    0.905797
      C                   10.942484   -0.577896    0.428481
      C                   12.283515   -0.895300    1.103652
      C                   12.802866   -2.293828    0.723383
      C                   14.303660   -2.390580    1.022546
      C                   12.056040   -3.380248    1.508281
      H                   -3.971743    1.871588   -0.839072
      H                   -9.894483    3.467231    1.143013
      H                  -11.139444   -2.723143    2.491943
      H                   -2.902664   -2.745806    0.121901
      H                   -3.127452   -0.534048   -1.818442
      H                   -1.989258    0.753032    0.043950
      H                   -1.225073   -0.813801   -0.218322
      H                   -2.303939   -0.577629    1.169379
      H                   -3.543942   -4.185902   -1.581612
      H                   -4.763611   -3.094413   -2.254413
      H                   -2.880990   -3.436020   -3.813951
      H                   -3.035086   -1.724063   -3.383246
      H                   -5.426480    5.379329   -1.240079
      H                   -4.814253    3.940870   -2.049085
      H                   -4.174428    4.410148   -0.457160
      H                   -7.707449    5.719082   -1.025637
      H                   -9.043716    5.319894    1.750247
      H                   -9.186822    6.753651    0.680885
      H                  -12.935051    2.315258    3.226144
      H                  -12.713878   -1.625588    3.247884
      H                  -12.955458   -0.101758    4.088044
      H                  -14.185704    0.848513    2.128895
      H                  -13.784329   -0.571622    1.147629
      H                  -14.870493   -0.715865    2.568269
      H                   -9.798172   -6.334093    1.921673
      H                  -10.917255   -5.017070    2.290652
      H                   -9.746714   -5.558010    3.501605
      H                   -5.239174   -5.189527   -0.447800
      H                   -1.114846   -4.997929    0.244201
      H                   -1.652845   -6.330646    1.359351
      H                   -1.603686   -4.627444    1.947904
      H                    1.089674   -1.567319   -4.088044
      H                    1.172188   -2.585399   -2.602590
      H                    1.049836    0.486560   -2.828784
      H                    0.921746   -2.501009   -0.512396
      H                    0.749437   -1.428138    0.873723
      H                    2.361521   -1.853834    0.298450
      H                    1.296005    1.752865   -0.810301
      H                    1.230975    0.930527    0.748394
      H                    3.680330    0.071964    0.195869
      H                    3.670214    1.253311   -1.114605
      H                    3.279266    1.810989    1.899143
      H                    2.984923    3.037266    0.660345
      H                    5.423686    2.621092   -0.169554
      H                    6.303818    4.153543    1.610281
      H                    4.955503    3.747690    2.675806
      H                    4.666802    4.624833    1.168771
      H                    5.581171    1.260899    2.604502
      H                    5.663558    0.409177    1.057908
      H                    7.657453    2.557380    2.028025
      H                    7.731488    1.788373    0.443677
      H                    8.007397   -0.459379    1.487656
      H                    7.741534    0.221822    3.106361
      H                   10.097389   -0.215807    2.810860
      H                    9.433721    1.873203    4.072707
      H                    9.711620    2.837351    2.607095
      H                   11.049875    1.968786    3.352603
      H                    9.720700    1.204097    0.146879
      H                   11.237355    1.495880    0.987191
      H                   10.208178   -1.378522    0.628804
      H                   11.075839   -0.562083   -0.668671
      H                   13.030624   -0.126520    0.820076
      H                   12.191260   -0.815694    2.203489
      H                   12.646233   -2.454757   -0.370577
      H                   14.870493   -1.631152    0.474361
      H                   14.499359   -2.243318    2.092463
      H                   14.704426   -3.370353    0.747273
      H                   10.974607   -3.312656    1.362529
      H                   12.253579   -3.281398    2.585541
      H                   12.380065   -4.380739    1.209889
{{< /code-block >}}
{{% /expand %}}

When your job finishes, take a look at the output to see what happened and make sure it completed successfully. Then we can do time-propagation.

# Running the time-dependent profiling

We change the input file accordingly. Change the {{< variable "CalculationMode" >}} from {{< code gs >}} to {{< code td >}}, and add the following lines:

{{< code-block >}}
 ##### TD

 T = 18
 dt = 0.003
 {{< variable "TDPropagator" >}} = aetrs
 {{< variable "TDTimeStep" >}} = dt
 
 # Profiling
 
 {{< variable "ProfilingMode" >}} = prof_memory
 {{< variable "TDMaxSteps" >}} = 30
 {{< variable "FromScratch" >}} = yes
{{< /code-block >}}

Now it is the time to do exactly the same TD run, changing the number of CPU processor cores. You have to change the XXX to powers of 2, 2^x. Start at 64 (which will be fastest) and divide by 2, in steps down to 4. (Running on 2 or 1 cores may not work.)

```bash
 #!/bin/bash
 #PBS -q regular
 #PBS -l mppwidth=XXX
 #PBS -l advres=benasque.348
 #PBS -l walltime=0:30:00
 #PBS -N testing_chl
 #PBS -V 
 
 module load octopus/4.1.2
 cd $PBS_O_WORKDIR
 aprun -n XXX octopus_mpi &> output_td_XXX
```

Different {{< code "profiling.000xxx" >}} folders will be created with each execution. We need to process them, mainly to be able to plot the information they contain. For that we can run the next script. It runs fine without any argument, but we can have more control in the files that it is going to process by using the following arguments: "analyze.sh 64 000004 2". The first argument is the biggest number of CPU processor cores that is going to be considered. The second optional argument is the number of the reference file that is going to be used. The third one is the starting number, i.e. the smallest number of CPU cores to consider.

```bash
 #!/bin/bash
 
 ## Copyright (C) 2012,2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## analyze.sh
 
 # Define the biggest number of processors.
 if [ -z "$1" ]; then
     last_proc=64
 else
     last_proc=$1
 fi
 
 # Define the reference file/folder
 if [ -z "$2" ]; then
     ref=000004
 else
     ref=$2
 fi
 
 # Define the starting value
 if [ -z "$3" ]; then
     start=1
 else
     start=$3
 fi
 
 #Initialise the output file
 echo "-  " > profile_$start
 for ((num=$start;num<=$last_proc;num*=2)); do
     echo $num >> profile_$start
 done
 
 rm -f tmp1
 
 # Analyze all profiling.XXXXXX/time.000000 to get the time per subroutine
 count=0
 for function_name in $(less profiling.$ref/time.000000  |  awk '{print $1}')
 do
     if [ $count -lt 4 ]; then
 	count=$((count + 1))
     else
 	echo $function_name >> tmp1
        -iterate over power of two profilings
 	for ((num=$start;num<=$last_proc;num*=2)); do
 	    folder=`printf 'profiling.%06d\n' $num `
 	    x=$(less $folder/time.000000 | grep "^$function_name " | awk '{print $3}' )
 	    zero=_"$x"_
 	    if [ "$zero" != "__" ]; then
 		echo $x >> tmp1
 	    else
 		echo "0" >> tmp1
 	    fi
 	done
 	paste profile_$start tmp1 > tmp2
 	rm tmp1
 	cp tmp2 profile_$start
     fi
 done
 
 echo "The result is in the \"profile_$start\" file"
```

At this point we should run "analyze.sh 64 000004 2". Thus, we will create files named "profile_2". You can take a look at the following columns in the profiling data:

* TIME_STEP; the iteration time. It has a good scaling.
* COMPLETE_DATASET; the whole time of the execution. In general it decreases, it is more obvious in a real execution, where the initialization time is the same and execution one is bigger.
* SYSTEM_INIT; initialization time. We were able to stop the increasing time, and now is almost constant independently of the number of processes.
* POISSON_SOLVER; execution time for the Poisson. It is somehow constant in this case, but now with the other solvers and domain parallelization.
* RESTART_WRITE; time for writing the restart files. It depends much in the system status, more than in the number of running processes. Could be heavily decreased if it is written to the local drive.

Now we can plot it using the following script:
```bash
 #!/bin/bash
 
 ## Copyright (C) 2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## plot_function.sh
 
 if [ $- -eq 0 ]; then
    function="TIME_STEP"
 else
    function=$1
 fi
 echo $function
 
 column_number=$( awk -v fun=$function '                                                                   
 { for(i=1;i<=NF;i++){                                                                                     
    if ($i == fun)                                                                                        
       {print i+1 }                                                                                       
    }                                                                                                     
 }' profile_2 )
 
 script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
 
 sed  "s/REF/$column_number/g" $script_dir/plot_ref > plot_base
 sed -i "s/FUNCTION/$function/g" plot_base
 gnuplot plot_base
```

We also need this auxiliary file:
```bash
 ## Copyright (C) 2014 J. Alberdi-Rodriguez
 ##
 ## This program is free software; you can redistribute it and/or modify
 ## it under the terms of the GNU General Public License as published by
 ## the Free Software Foundation; either version 2, or (at your option)
 ## any later version.
 ##
 ## This program is distributed in the hope that it will be useful,
 ## but WITHOUT ANY WARRANTY; without even the implied warranty of
 ## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 ## GNU General Public License for more details.
 ##
 ## You should have received a copy of the GNU General Public License
 ## along with this program; if not, write to the Free Software
 ## Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
 ## 02111-1307, USA.
 ##
 ## plot_ref
 
 set t postscript eps enhanced color solid
 set output "gnuplot.eps"
 
 set xlabel "MPI processes"
 set ylabel "t (s)"
 set logscale yx 2
 
 plot "profile_2" u 1:REF w linespoint t "FUNCTION 2^x"
```

Something else you can try is 12, 24, 48 and 96 cores, because each of the nodes has 24 cores. In this case, you would need "analyze.sh 96 000003 3" to make "profile_3", and then in the plotting script,

```bash
 plot "profile_2" u 1:REF w linespoint t "FUNCTION 2^x", "profile_3" u 1:REF w lp t "FUNCTION 3·(2^x)"
```

## Parallelization in domains vs states

We can divide up the work among the processors in different ways, by dividing up the points into domains for each processor, or dividing the states into groups for each processor, or a combination of both. Try out different combinations by adding to your input file

{{< code-block >}}
 {{< variable "ParStates" >}} = 2
 {{< variable "ParDomains" >}} = 12
{{< /code-block >}}

and run on 24 cores, with different numbers in the first two fields whose product is the total number of processors (''e.g.'' 6 x 4, 3 x 8, ...).

## PFFT Poisson solver

Another thing you can try is to compare the PFFT (parallel FFT) Poisson solver against the one we were using before (look in the output file to see which one it was). You will need to use this aprun line in your job script instead of the previous one:

```bash
 aprun -n XXX /global/homes/j/joseba/octopus/bin/octopus_mpi &> output_td_XXX
```

in order to use a different octopus compilation that uses that library, and add these lines to your input file:

{{< code-block >}}
 {{< variable "PoissonSolver" >}} = fft
 {{< variable "FFTLibrary" >}} = pfft
{{< /code-block >}}

Compare some runs against one on a similar number of processors that you did previously. How does the time for this solver compare? You can also try ordinary FFTW (not parallel) with

{{< code-block >}}
 {{< variable "FFTLibrary" >}} = fftw
{{< /code-block >}}

## Parallelization of the ground state

We can also try different parameters and algorithms to see their effect on the speed of the ground-state calculation, for 24, 48, or 96 processors. Look each up in the variable reference to see what they mean, and see which of the options you were using in the previous runs.
* parallelization in domains vs states (as above)
* Eigensolver = rmmdiis, plan, cg, cg_new, lobpcg.
* StatesOrthogonalization = cholesky_serial, cholesky_parallel, mgs, qr
* SubspaceDiagonalization = standard, scalapack
* linear-combination of atomic orbitals (LCAO) for initial guess: {{< variable "LCAOAlternative" >}} = yes, no. In this case, add {{< variable "MaximumIter" >}} = 0 to do just LCAO rather than the whole calculations.

---
title: "ARPES"
tags: ["Tutorial", "Advanced", "Time-dependent", "Surface", "User defined potential"]
tutorials: "Periodic Systems"
difficulties: "advanced"
theories: "Independent particles"
system_types: "semi-periodic"
species_types: "User defined potential"
description: "This tutorial aims at introducing the user to the basic concepts needed to calculate ARPES with octopus."
#series: "Tutorial"
Weight: 20
---

This tutorial aims at introducing the user to the basic concepts needed to calculate ARPES with {{< octopus >}}.
We choose as test case a fictitious non-interacting two-dimensional system with semi-periodic boundary conditions.

## Ground state

Before starting with the time propagation we need to obtain the ground state of the system.

To simulate ARPES simulations with the tSURRF implementation in {{<octopus>}}[^footnote-1] one needs to explicitly model the surface of the material. This requirement often leads to heavy simulations. For this reason, in this tutorial, we illustrate the procedure on a model system. Calculating ARPES with TDDFT on an ab-initio model for a surface involves the same steps. 

### Input

The system we choose is a 2D toy model potential simulating an atomic chain. The dimensionality is selected by telling octopus to use on 2 spatial dimensions (xy) 
with {{<code-inline>}}{{<variable "Dimensions">}} = 2{{</code-inline>}} and by imposing periodic boundary conditions along x with 
{{<code-inline>}}{{<variable "PeriodicDimensions">}} = 1{{</code-inline>}}.


{{< code-block >}}
#include_input doc/tutorials/other/arpes/01-gs/inp
{{< /code-block >}}


Octopus can calculate ARPES on a path in reciprocal space. For this reason we have to specify the path with the {{<variable "KPointsPath">}} 
block already at the ground state level. This is needed to generate the KS wave functions that will be evolved in the time propagation. 
It also provides the bands structure that we will use as reference to cross-check the quality of the spectrum.

Note that, since this is a toy model with independent electrons, {{<code-inline>}}{{<variable "TheoryLevel">}} = independent_particles{{</code-inline>}}, 
we do not need to sample the BZ and the grid in reciprocal space is constituted only by the gamma point, 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/01-gs/inp kpoints
{{< /code-block >}}



## Time-dependent run

To calculate the ARPES spectrum we excite the system with a laser field to excite electrons into the vacuum, i.e. on the non-periodic dimension. The photoelectron detection probability constituting ARPES is calculated by analyzing the flux of the ionization current trough a surface positioned at a certain distance from the surface of the system. This current is obtained by propagating in time the KS orbitals under the effect of the laser. 

### Input

This is how the input file should look for the time propagation.

{{% expand "full input file" %}}
{{< code-block >}}
#include_input doc/tutorials/other/arpes/02-td/inp
{{< /code-block >}}
{{% /expand %}}


Where we just changed the {{<variable "CalculationMode">}} to {{<code "td">}} to activate the time propagation. 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp mode
{{< /code-block >}}

We then have to specify the laser field parameters as following.

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp laser_field
{{< /code-block >}}

Since we deal with periodic system we have to specify the light-matter coupling in the so called "velocity gauge" where the field is described by a time-dependent vector potential.
We specify the field as a "vector_potential" polarized along x and with carrier frequency "wpr"

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp vector_potential
{{< /code-block >}}

and a sin^2 envelope function "probe" specified by an analytical expression of time, t.

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp tdfunctions
{{< /code-block >}}

We also specified the total time of the propagation to be exactly synchronized with the pulse switch-off time by setting 
{{<code-inline>}}{{<variable "TDPropagationTime">}} = Tpr{{</code-inline>}}. 
For more details look at the {{<tutorial "basics/time-dependent_propagation" "Time-dependent propagation">}} tutorial.

Photoelectrons ejected from the system by the laser will eventually bounce back from the boundary of the simulation box along the non-periodic dimension (y). We therefore employ absorbing boundary conditions to prevent spurious reflections with the following code block

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp cap
{{< /code-block >}}


Here we employ complex absorbing potential (CAP) boundaries with {{<code-inline>}}{{<variable "AbsorbingBoundaries">}} = cap{{</code-inline>}}, 
specify the CAP parameters {{<variable "ABCapHeight">}} and {{<variable "ABShape">}} to have maximal absorption in the energy region where we expect the photoelectrons.[^footnote-2] 

Finally we specify the parameters for the evaluation of the ARPES spectrum. 

{{< code-block >}}
#include_input_snippet doc/tutorials/other/arpes/02-td/inp arpes
{{< /code-block >}}


Here we tell octopus to use tSURFF with {{<code-inline>}}{{<variable "PhotoElectronSpectrum">}} = pes_flux {{</code-inline >}}, 
specify the position of the analyzing surface at the onset of the absorbing boundaries with 
{{<code-inline>}}{{<variable "PES_Flux_Lsize">}} = Lmin{{</code-inline>}} and define the parameters of the energy grid for the 
final spectrum with the block {{<variable "PES_Flux_EnergyGrid">}}.
The tSURFF method allows to calculate the momentum-resolved photoelectron probability $P({\bf p})$ on an arbitrary grid in momentum space.
The ARPES spectrum is than obtained rewriting $P({\bf p})$ as a function of the total kinetic energy and the momentum parallel to the surface, ${\bf p}\_\parallel$ 
(equivalent to the crystal momentum ${\bf k}$ ) $P({\bf p}\_\parallel={\bf k}, E=({\bf p}\_\parallel + {\bf p}\_\perp)^2/2)$. 
Since a simple cartesian grid in momentum results in a deformed ARPES grid octopus can generate a grid that compensates the deformation such the final ARPES spectrum is in a cartesian grid. 
This grid is generated with {{<code-inline>}}{{<variable "PES_Flux_ARPES_grid">}} = yes{{</code-inline>}}. 
Since this option is true by default for semi-periodic systems we do not have to specify it in the input file.




### Output

If the code runs correctly the standard output should present a Photoelectron section like this:

{{< code-block >}}
*************************** Photoelectron ****************************
Info: Calculating PES using t-surff technique.
Input: [PES_Flux_Shape = pln]
Input: [PES_Flux_Lsize = (8.150,30.00)]
Input: [PES_Flux_Parallelization = pf_none]
Input: [PES_Flux_Momenutum_Grid = cartesian]
Energy grid (Emin, Emax, DE) [H]:  (   1.711,   1.911, 0.10E-01)
Momentum linear grid (Pmin, Pmax, DP) [me*bH(2pi/h)]:  (   1.850,   1.955, 0.14E+00)
Input: [PES_Flux_ARPES_grid = yes]
Number of points with E < p//^2/2 = 0 [of 462]
Input: [PES_Flux_UseSymmetries = no]
Input: [PES_Flux_RuntimeOutput = no]
Info: Total number of surface points = 16
Info: Total number of momentum points = 21
{{< /code-block >}}

## ARPES

The spectral information on the photoelectrons is stored in {{<file "restart/td/pesflux*">}} binary files and can be analyzed in post-processing.

### The {{< command "oct-photoelectron_spectrum" >}} utility

{{<octopus>}} provides an utility called {{<command "oct-photoelectron_spectrum">}} to process the {{<file "pesflux*">}} files and obtain the spectrum.


###  Input  

The utility by default requires no input from the user in calculations designed to obtain an ARPES spectrum. If the td run has been performed with 
{{<code-inline>}}{{<variable "PES_Flux_ARPES_grid">}} = yes{{</code-inline>}} it should have all the information that it needs. 

Just run the utility in the same path where the inp file resides.


### Output

This is what you should get:

{{< code-block >}}
************************** Kpoint selection **************************
Will use a zero-weight path in reciprocal space with the following points
       2    0.000000 |   -0.500000    0.000000 |
       3    0.000000 |   -0.450000    0.000000 |
       4    0.000000 |   -0.400000    0.000000 |
       5    0.000000 |   -0.350000    0.000000 |
       6    0.000000 |   -0.300000    0.000000 |
       7    0.000000 |   -0.250000    0.000000 |
       8    0.000000 |   -0.200000    0.000000 |
       9    0.000000 |   -0.150000    0.000000 |
      10    0.000000 |   -0.100000    0.000000 |
      11    0.000000 |   -0.050000    0.000000 |
      12    0.000000 |    0.000000    0.000000 |
      13    0.000000 |    0.050000    0.000000 |
      14    0.000000 |    0.100000    0.000000 |
      15    0.000000 |    0.150000    0.000000 |
      16    0.000000 |    0.200000    0.000000 |
      17    0.000000 |    0.250000    0.000000 |
      18    0.000000 |    0.300000    0.000000 |
      19    0.000000 |    0.350000    0.000000 |
      20    0.000000 |    0.400000    0.000000 |
      21    0.000000 |    0.450000    0.000000 |
      22    0.000000 |    0.500000    0.000000 |

**********************************************************************

Read PES restart files.
Zenith axis: (      0.00,       0.00,       1.00)
[ 42/ 42] 100%|****************************************|     --:-- ETA

***************** ARPES cut on reciprocal space path *****************
Done
**********************************************************************
{{< /code-block >}}


If all goes well the file is {{<file "PES_ARPES.path">}} will be created containing the ARPES spectrum evaluated on the k-point path.
The spectrum can be visualized with {{<command "gnuplot">}} as a density plot with 

{{< code-block >}}
set pm3d map
sp "PES_ARPES.path" u 1:4:5
{{< /code-block >}}


{{< figure src="images/arpes_2d.png" width="500px" caption="ARPES spectrum for a 2D atomic chain model system." >}}

The ARPES spectrum looks a bit blocky to improve the resolution one has to increase the number of kpoints in "KPointsPath" (requires recomputing the groundstate) and decrease the energy spacing in {{<variable "PES_Flux_EnergyGrid">}}. 
Keep in mind that larger grids implies heavier simulations and longer run times.


Some questions to think about:

* How does the ARPES spectrum compares with the band structure? Plot one on top of the other with {{<command gnuplot>}}.

* How to choose the photoelectron energy range? In the previous example we chose Emin =  wpr - 0.2 and Emax =  wpr in order to see photoelectrons emitted from the valence band. How can we change {{<variable "PES_Flux_EnergyGrid">}} to see the band below? What about conduction band?

* Play around with the laser carrier, "wpr", and pulse envelope, "Tpr". Can you identify the impact of these parameters on the spectrum?



## References

[^footnote-1]: {{< article title="A First-Principles Time-Dependent Density Functional Theory Framework for Spin and Time-Resolved Angular-Resolved Photoelectron Spectroscopy in Periodic Systems" authors="U. De Giovannini, H. Hübener, and A. Rubio" journal="JCTC" volume="13" pages="265" year="2017" doi="10.1021/acs.jctc.6b00897" >}}

[^footnote-2]: {{< article title="Modeling electron dynamics coupled to continuum states in finite volumes with absorbing boundaries" authors="U. De Giovannini, A. H. Larsen, A. Rubio, and A. Rubio" journal="EPJB" volume="88" pages="1" year="2015" doi="10.1140/epjb/e2015-50808-0" >}}


{{< tutorial-footer >}}


---
title: "Wires and slabs"
tags: ["Advanced", "Ground State", "Unoccupied", "Chain", "Slab", "Pseudopotentials", "DFT", "Band Structure"]
tutorials: ["Periodic Systems"]
difficulties: "advanced"
theories: ["DFT"]
calculation_modes: ["Ground State", "Unoccupied"]
system_types: ["Chain", "Slab"]
species_types: "Pseudopotentials"
features: "Band Structure"
description: "Use the flexibility of the real-space grid to treat systems that are periodic in only one or two dimensions."
#series: "Tutorial"
Weight: 15
---


In this tutorial we will explain how to use the flexibility of the real-space grid to treat systems that are periodic in only one or two dimensions. As examples we will use a Na chain and a hexagonal boron nitride (h-BN) monolayer. 

## Introduction

In the [Periodic systems](../Periodic systems) tutorial, we saw that the {{< variable "PeriodicDimensions" >}} input variable controls the number of dimensions to be considered as periodic. In that tutorial we only considered the case <tt>{{< variable "PeriodicDimensions" >}} = 3</tt>. Let's now see in detail the different cases:

* <tt>{{< variable "PeriodicDimensions" >}} = 0</tt> (which is the default) gives a finite system calculation, since Dirichlet zero boundary conditions are used at all the borders of the simulation box;

* <tt>{{< variable "PeriodicDimensions" >}} = 1</tt> means that only the ''x'' axis is periodic, while in all the other directions the system is confined. This value must be used to simulate, for instance, a single infinite wire.

* <tt>{{< variable "PeriodicDimensions" >}} = 2</tt> means that both ''x'' and ''y'' axis are periodic, while zero boundary conditions are imposed at the borders crossed by the ''z'' axis. This value must be used to simulate, for instance, a single infinite slab.

* <tt>{{< variable "PeriodicDimensions" >}} = 3</tt> means that the simulation box is a primitive cell for a fully periodic infinite crystal. Periodic boundary conditions are imposed at all borders.


It is important to understand that performing, for instance, a {{< variable "PeriodicDimensions" >}}<tt> = 1</tt> calculation in {{< octopus >}} is not quite the same as performing a {{< variable "PeriodicDimensions" >}}<tt> = 3</tt> calculation with a large supercell. In the infinite-supercell limit the two approaches reach the same ground state, but this does not hold for the excited states of the system.

Another point worth noting is how the Hartree potential is calculated for periodic systems. In fact the discrete Fourier transform that are used internally in a periodic calculation would always result in a 3D periodic lattice of identical replicas of the simulation box, even if only one or two dimensions are periodic. Fortunately {{< octopus >}} includes a clever system to exactly truncate the long-range part of the Coulomb interaction, in such a way that we can effectively suppress the interactions between replicas of the system along non-periodic axes [^cutoff]
. This is done automatically, since the value of {{< variable "PoissonSolver" >}} in a periodic calculation is chosen according to {{< variable "PeriodicDimensions" >}}. See also the variable documentation for {{< variable "PoissonSolver" >}}.

## Atomic positions

An important point when dealing with semi-periodic systems is that the coordinates of the atoms along the aperiodic directions must be centered around 0, as this is the case for isolated systems.
For slabs, this means that the atoms of the slab must be centered around z=0. If this is not the case, the calculation will lead to wrong results.

## Sodium chain

Let us now calculate some bands for a simple single Na chain (i.e. not a crystal of infinite parallel chains, but just a single infinite chain confined in the other two dimensions).

### Ground-state

First we start we the ground-state calculation using the following input file:

```text
 {{< variable "CalculationMode" >}} = gs
  {{< variable "UnitsOutput" >}} = ev_angstrom
 {{< variable "ExperimentalFeatures" >}} = yes
 
 {{< variable "PeriodicDimensions" >}} = 1
 
 {{< variable "Spacing" >}} = 0.3*angstrom
  
 %{{< variable "Lsize" >}}
  1.99932905*angstrom | 5.29*angstrom | 5.29*angstrom
 %
 
 %{{< variable "Coordinates" >}}
  "Na" | 0.0 | 0.0 | 0.0
 %
 
 %{{< variable "KPointsGrid" >}}
  9 | 1 | 1
 %
 {{< variable "KPointsUseSymmetries" >}} = yes
```

Most of these input variables were already introduced in the [Periodic systems](../Periodic systems) tutorial. Just note that the ''k''-points are all along the first dimension, as that is the only periodic dimension.  

The output should be quite familiar, but the following piece of output confirms that the code is indeed using a cutoff for the calculation of the Hartree potential, as mentioned in the Introduction:
```text

****************************** Hartree *******************************
The chosen Poisson solver is 'fast Fourier transform'
Input: [PoissonFFTKernel = cylindrical]
**********************************************************************

Input: [FFTLibrary = fftw]
Info: FFT grid dimensions       = 13 x 75 x 75
      Total grid size           = 73125 (   0.6 MiB )
      Inefficient FFT grid. A better grid would be: 14 75 75
Info: Poisson Cutoff Radius     =  11.2 A
```

The cutoff used is a cylindrical cutoff and, by comparing the FFT grid dimensions with the size of the simulation box, we see that the Poisson solver is using a supercell doubled in size in the ''y'' and ''z'' directions.

At this point you might want to play around with the number of ''k''-points until you are sure the calculation is converged.

### Band structure

We now modify the input file in the following way:

{{< figure src="/images/Na_chain_bands_1D_3D.png" width="500px" caption="Band structure for a infinite chain of Sodium atoms, calculated for a single chain (purple lines), and a 3D-periodic crystal of chains in a supercell (green lines)." >}}

```text
 {{< variable "CalculationMode" >}} = unocc
 {{< variable "UnitsOutput" >}} = ev_angstrom
 {{< variable "ExperimentalFeatures" >}} = yes
 
 {{< variable "PeriodicDimensions" >}} = 1
 
 {{< variable "Spacing" >}} = 0.3*angstrom
  
 %{{< variable "Lsize" >}}
  1.99932905*angstrom | 5.29*angstrom | 5.29*angstrom
 %
 
 %{{< variable "Coordinates" >}}
  "Na" | 0.0 | 0.0 | 0.0
 %
 
 {{< variable "ExtraStates" >}} = 6
 {{< variable "ExtraStatesToConverge" >}} = 4
 
 %{{< variable "KPointsPath" >}}
  14
  0.0 | 0.0 | 0.0
  0.5 | 0.0 | 0.0
 %
 {{< variable "KPointsUseSymmetries" >}} = no
```

and run the code. You might notice the comments on the the LCAO in the output. What's going on? Why can't a full initialization with LCAO be done?

You can now plot the band structure using the data from the {{< file "static/bandstructure" >}} file, just like in the [Periodic systems](../Periodic systems) tutorial.

In the introduction we mentioned that performing a calculation with <tt>{{< variable "PeriodicDimensions" >}} = 1</tt> is not quite the same as performing a <tt>{{< variable "PeriodicDimensions" >}}= 3</tt> calculation with a large supercell. Lets now check this. Re-run both the ground-state and the unoccupied calculations, but setting <tt>{{< variable "PeriodicDimensions" >}} = 3</tt> in the above input files. Before doing so, make sure you copy the {{< file "static/bandstructure" >}} file to a different place so that it is not overwritten (better yet, run the new calculations in a different folder). You can see the plot of the two band structures on the right. More comments on this in ref. [^cutoff].

## h-BN monolayer

Hexagonal boron nitride (h-BN) is an insulator widely studied which has a similar structure to graphene. Here we will describe how to get the band structure of an h-BN monolayer. 

### Ground-state calculation
We will start by calculating the ground-state using the following input file:

```text
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = ev_angstrom
 {{< variable "ExperimentalFeatures" >}} = yes
 
 {{< variable "PeriodicDimensions" >}} = 2
 
 {{< variable "Spacing" >}} = 0.20*angstrom
 
 BNlength = 1.445*angstrom
 a = sqrt(3)*BNlength
 L = 40
 %{{< variable "LatticeParameters" >}}
  a | a | L
 %
 
 %{{< variable "LatticeVectors" >}}
   1.0 | 0.0       | 0.0
  -1/2 | sqrt(3)/2 | 0.0
   0.0 | 0.0       | 1.0
 %
  
 %{{< variable "ReducedCoordinates" >}}
  'B' | 0.0 | 0.0 | 0.0
  'N' | 1/3 | 2/3 | 0.0
 % 
 
 {{< variable "PseudopotentialSet" >}}=hgh_lda
 
 {{< variable "LCAOStart" >}}=lcao_states 
 
 %{{< variable "KPointsGrid" >}}
   12   | 12   | 1
 %
 {{< variable "KPointsUseSymmetries" >}} = yes
```

Most of the file should be self-explanatory, but here is a more detailed explanation for some of the choices:

* <tt>{{< variable "PeriodicDimensions" >}} = 2</tt>: A layer of h-BN is periodic in the ''x''-''y'' directions, but not in the ''z'' direction, so there are two periodic dimensions.

* <tt>{{< variable "LatticeParameters" >}}</tt>: Here we have set the bond length to 1.445 Å. The box size in the ''z'' direction is ''2 L'' with ''L'' large enough to describe a monolayer in the vacuum. Remember that one should always check the convergence of any quantities of interest with the box length value.

If you now run the code, you will notice that the cutoff used for the calculation of the Hartree potential is different than for the Sodium chain, as is to be expected:
```text

****************************** Hartree *******************************
The chosen Poisson solver is 'fast Fourier transform'
Input: [PoissonFFTKernel = planar]
**********************************************************************

Input: [FFTLibrary = fftw]
Info: FFT grid dimensions       = 13 x 13 x 225
      Total grid size           = 38025 (   0.3 MiB )
      Inefficient FFT grid. A better grid would be: 14 14 225
Info: Poisson Cutoff Radius     =  22.5 A
```

### Band Structure
After the ground-state calculation, we will now calculate the band structure. This is the input for the non-self consistent calculation:

```text
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = ev_angstrom
 {{< variable "ExperimentalFeatures" >}} = yes
 
 {{< variable "PeriodicDimensions" >}} = 2
 
 {{< variable "Spacing" >}} = 0.20*angstrom
 
 BNlength = 1.445*angstrom
 a = sqrt(3)*BNlength
 L = 40
 %{{< variable "LatticeParameters" >}}
  a | a | L
 %
 
 %{{< variable "LatticeVectors" >}}
   1.0 | 0.0       | 0.0
  -1/2 | sqrt(3)/2 | 0.0
   0.0 | 0.0       | 1.0
 %
  
 %{{< variable "ReducedCoordinates" >}}
  'B' | 0.0 | 0.0 | 0.0
  'N' | 1/3 | 2/3 | 0.0
 % 
 
 {{< variable "PseudopotentialSet" >}}=hgh_lda
 
 {{< variable "LCAOStart" >}}=lcao_states 
 
 {{< variable "ExtraStatesToConverge" >}}  = 4
 {{< variable "ExtraStates" >}}  = 8
 
 %{{< variable "KPointsPath" >}}
  12  | 7   | 12   - Number of k point to sample each path
  0.0 | 0.0 | 0.0  - Reduced coordinate of the 'Gamma' k point
  1/3 | 1/3 | 0.0  - Reduced coordinate of the 'K' k point
  1/2 | 0.0 | 0.0  - Reduced coordinate of the 'M' k point
  0.0 | 0.0 | 0.0  - Reduced coordinate of the 'Gamma' k point
 %
 {{< variable "KPointsUseSymmetries" >}} = no
```

{{< figure src="/images/Tutorial_band_structure_HBN.png" width="500px" caption="Band structure of a monolayer h-BN." >}}

In this case, we chose the following path to calculate the band structure: Gamma-K, K-M, M-Gamma, with a sampling of 12-7-12 ''k''-points.

On the right is the resulting plot of the occupied bands and the first four unoccupied bands from the {{< file "static/bandstructure" >}} file.


[^cutoff]: {{< article title="Exact Coulomb cutoff technique for supercell calculations" authors="C. A. Rozzi, D. Varsano, A. Marini, E. K. U. Gross, and A. Rubio" journal="Phys. Rev. B" volume="73" pages="205119" year="2006" doi="10.1103/PhysRevB.73.205119" >}}

{{< tutorial-footer >}}


---
title: "DFT+U"
tags: ["Tutorial", "Advanced", "Ground State", "Bulk", "Pseudopotentials", "DFT+U"]
tutorials: "Periodic Systems"
difficulties: "advanced"
theories: "DFT+U"
system_types: "bulk"
species_types: "Pseudopotentials"
description: "The objective of this tutorial is to give a basic idea of how DFT+U in octopus works"
#series: "Tutorial"
Weight: 20
---


The objective of this tutorial is to give a basic idea of how DFT+U in {{< octopus >}} works.

### Input

As a prototypical example for DFT+U, we will consider bulk NiO in its anti-ferromagnetic configuration. We will neglect the small lattice distortion and consider only its cubic cell.

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=hscv_pbe
  a = 7.8809
  %{{< variable "LatticeParameters" >}}
    a | a | a
  %
  %{{< variable "LatticeVectors" >}}
   0.0 | 1/2 | 1/2
   1/2 | 0.0 | 1/2
   1.0 | 1.0 | 0.0
  %
  
  %{{< variable "Species" >}}
  "Ni" | species_pseudo | hubbard_l | 2 | hubbard_u | 5.0*eV
  %
  
  {{< variable "DFTULevel" >}} = dft_u_empirical
  
  %{{< variable "ReducedCoordinates" >}}
   "Ni" | 0.0 | 0.0 | 0.0
   "Ni" | 0.0 | 0.0 | 0.5
   "O"  | 0.5 | 0.5 | 0.25
   "O"  | 0.5 | 0.5 | 0.75
  %
  {{< variable "Spacing" >}} = 0.5
  %{{< variable "KPointsGrid" >}}
  2 | 2 | 2
  %
  {{< variable "ParDomains" >}} = no
  {{< variable "ParKPoints" >}} = auto
  
  {{< variable "SpinComponents" >}} = polarized
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   8.0
  -8.0
   0.0
   0.0
  %
  
  {{< variable "OutputLDA_U" >}} = occ_matrices
{{< /code-block >}}

As we are interested by the antiferromagnetic order, the primitive cell is doubled along the last lattice vector.
To help the convergence, an initial guess should be added, by adding to the input file the variables {{< variable "GuessMagnetDensity" >}} and {{< variable "AtomsMagnetDirection" >}}.

In order to perform a calculation with a U of 5eV on the 3d orbitals (corresponding to the quantum number l=2),
we define a block Species, where hubbard_l specifies the orbitals (l=0 for s orbitals, l=1 for p orbitals, ...) and hubbard_u is used to set the value of the effective Hubbard U.

In order to activate the DFT+U part, one finally needs to specify the level of DFT+U used . This is done using the variable {{< variable "DFTULevel" >}}. At the moment there is three possible options for this variable, which correspond to no +U correction (dft_u_none), an empirical correction (dft_u_empirical) or the ab initio U correction based on the ACBN0 functional[^footnote-1] (dft_u_acbn0).

Some specific outputs can then be added, such as the density matrix of the selected localized subspaces. 

### Output

## References







---------------------------------------------
[^footnote-1]: {{< article title="Reformulation of $\mathrm{DFT}+U$ as a Pseudohybrid Hubbard Density Functional for Accelerated Materials Discovery" authors="Agapito, Luis A. and Curtarolo, Stefano and Buongiorno Nardelli, Marco" journal="Phys. Rev. X" volume="5" pages="011006" year="2015" doi="10.1103/PhysRevX.5.011006" >}}


{{< tutorial-footer >}}
---
title: "Magnons"
tags: ["Tutorial", "Advanced", "Bulk", "Magnons in real time"]
tutorials: "Periodic Systems"
difficulties: "advanced"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "bulk"
species_types: "Pseudopotentials"
features: "Magnons in real time"
#series: "Tutorial"
description: "This tutorial gives the basic idea of how it is possible to compute magnons and transverse spin susceptibilities from a real-time calculation using octopus."
Weight: 30
---


The objective of this tutorial is to give a basic idea of how it is possible to compute magnons and transverse spin susceptibilities from a real-time calculation using {{< octopus >}}.

{{< notice warning >}}
This tutorial requires to be able to run calculations which have a non-negligible numerical cost and cannot be run locally.
{{< /notice >}}

{{< notice note >}}
This tutorial is still under construction.
{{< /notice >}}


### Magnon kick in supercell

As a first approach, we will investigate how to obtain transverse spin susceptibilities using supercells adapted to host a spin spiral of a specific momentum.

#### Ground-state input file

As a prototypical example for investigating spin-waves, we will consider bulk iron, which is a ferromagnetic material. We start here by computing its ground-state by constructing a supercell corresponding to the cubic cell of iron, doubled along one direction.
For this, we use the following input file:

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | 2*a
    90| 90| 90
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
   "Fe" | 1/2 | 1/2 | 1/4
   "Fe" | 0.0 | 0.0 | 1/2
   "Fe" | 1/2 | 1/2 | 3/4
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 2
  %
  {{< variable "Smearing" >}} = 0.1*eV
  {{< variable "SmearingFunction" >}} = fermi_dirac
  {{< variable "LCAOStart" >}} = lcao_states
  {{< variable "SpinComponents" >}} = spinors
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
   0.0 | 0.0 | 4.0
  %
  {{< variable "EigenSolver" >}} = rmmdiis
  {{< variable "ConvRelDens" >}} = 1e-7
  {{< variable "ExtraStates" >}} = 20
{{< /code-block >}}

Running this input file will give the ground state of bulk iron in a supercell. 
It is important to note here that we are using Pauli spinors to represent the wavefunctions. This is needed because we will later on kick the system with a spiral in spin space, that cannot be represented by collinear spins.

Note that the present input file is not converged in the number of k-points and should not be considered for production runs.
Similarly, the smearing of the occupation is set here to a high temperature of 0.1 eV, which is also not realistic for practical applications. This is done here because the number of k-point is too low in our example to properly sample the energy window close to the Fermi energy. 


#### Time-dependent run

Now that we constructed the supercell, we can investigate a spin spiral. The momemtum q corresponding to this cell has a value of q=(0,0,2 pi/a/2), where a is lattice parameter of iron, which we took as 2.867 angstrom in our example.

Let us now look at how to specify a magnon kick in {{< octopus >}}.
There are three different points that need to be specified:
* The strength of the kick. This determines how strongly we perturb the system.
* The momentum of the kick. This is the momentum of the spin spiral we are imposing to the system's spins.
* The easy axis of the material. This is usually the z direction, but the code allows for defining an arbitrary direction. This direction is used to determine the transverse magnetization from the total magnetization computed in Cartesian coordinates. 

First of all, we need to set that we will use a magnon kick:

{{< code-block >}}
 {{< variable TDDeltaStrengthMode >}} = kick_magnon
{{< /code-block >}}

Then we specify the strength of the perturbation:

{{< code-block >}}
 {{< variable TDDeltaStrength >}} = 0.01
{{< /code-block >}}

The strength of the perturbation should be converged such that the results are guarantied to be valid within linear response. The means that changing the strength of the kick should lead to the same transverse spin susceptibility than the original calculation, as the extracted susceptibility (see next section) does not depend on the kick strength within linear response.
The momentum of spin spiral is set using the block

{{< code-block >}}
  %{{< variable TDMomentumTransfer >}}
   0 | 0 | 2*pi/a/2
  %
{{< /code-block >}}

Finally, the easy axis of the material is defined by

{{< code-block >}}
  %{{< variable TDEasyAxis >}}
    0 | 0 | 1
  %
{{< /code-block >}}

The full input file should read as

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | 2*a
    90| 90| 90
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
   "Fe" | 1/2 | 1/2 | 1/4
   "Fe" | 0.0 | 0.0 | 1/2
   "Fe" | 1/2 | 1/2 | 3/4
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 2
  %
  {{< variable "RestartFixedOccupations" >}} = yes
  {{< variable "SpinComponents" >}} = spinors

  {{< variable "ExtraStates" >}} = 8

  {{< variable "RestartWriteInterval" >}} = 5000
  {{< variable "TDDeltaStrength" >}} = 0.01
  {{< variable "TDDeltaStrengthMode" >}} = kick_magnon
  %{{< variable "TDMomentumTransfer" >}}
   0 | 0 | 2*pi/a/2 
  %
  %{{< variable "TDEasyAxis" >}}
   0 | 0 | 1
  %
  {{< variable "TDTimeStep" >}} = 0.075
  {{< variable "TDPropagator" >}} = aetrs
  {{< variable "TDExponentialMethod" >}} = lanczos
  {{< variable "TDExpOrder" >}} = 16
  {{< variable "TDPropagationTime" >}} = 1800
  {{< variable "TDOutput" >}} = total_magnetization + energy
{{< /code-block >}}

The parameters for time-dependent runs are already explained in other tutorials and are not further detailed here. 
Here the time-propagation is performed only up to 1800 atomic units. This frequency resolution for this time propagation is probably too large for most applications, but we use this value here to maintain the tutorial feasible within a reasonable amount of time.

Looking at the total energy (td.general/energy) of the system, one realizes that it starts to increase, which is not correct. This is due to the poor convergence parameters used for the calculation. In practical calculation, it is strongly advised to check that the energy is conserved throughout the complete simulation.

#### Computing the transverse spin susceptibility

In order to extract the transverse spin susceptibilities from the total magnetization, one needs to use the utility oct-spin_susceptibility.
This utility produces files names td.general/spin_susceptibility_qXXX, where XXX is the index of the q-vector. There is typically only one file produced, except in the multi-q kick mode, see below for some example.

### Magnon kick using the generalized Bloch theorem

If the Hamiltonian does not include any off-diagonal terms in spin space, it is possible to use the so-called generalized Bloch theorem (GBT).
Thanks to the GBT, it is possible to investigate spin spirals using only the primitive cell of the system. This comes at the cost of modifying the periodic boundary conditions by some twisted boundary conditions in which a different phase is applied to each components of the Pauli spinors.

In order to investigate it, we first need to prepare the ground-state of the system in its primitive cell.
For this, we modify the above input file to define only the primitive cell

{{< code-block >}}
  {{< variable "CalculationMode" >}} = gs
  {{< variable "PeriodicDimensions" >}} = 3 
  {{< variable "BoxShape" >}} = parallelepiped
  {{< variable "ExperimentalFeatures" >}} = yes
  {{< variable "PseudopotentialSet" >}}=pseudodojo_lda
  a = 2.867*Angstrom
  %{{< variable "LatticeParameters" >}}
    a | a | a
  %
  %{{< variable "LatticeVectors" >}}
   -0.5 | 0.5 | 0.5
    0.5 |-0.5 | 0.5
    0.5 | 0.5 |-0.5
  %
  %{{< variable "ReducedCoordinates" >}}
   "Fe" | 0.0 | 0.0 | 0.0
  %
  {{< variable "Spacing" >}} = 0.35
  %{{< variable "KPointsGrid" >}}
  4 | 4 | 4
  %
  {{< variable "Smearing" >}} = 0.1*eV
  {{< variable "SmearingFunction" >}} = fermi_dirac
  {{< variable "LCAOStart" >}} = lcao_states
  {{< variable "SpinComponents" >}} = spinors
  {{< variable "GuessMagnetDensity" >}} = user_defined
  %{{< variable "AtomsMagnetDirection" >}}
   0.0 | 0.0 | 4.0
  %
  {{< variable "EigenSolver" >}} = rmmdiis
  {{< variable "ConvRelDens" >}} = 1e-7
  {{< variable "ExtraStates" >}} = 10
{{< /code-block >}}

Once the ground state is converged, we need to run the time-dependent calculation.
This is done using the same variable as used before, only adding the line
{{< code-block >}}
 {{< variable "SpiralBoundaryCondition" >}} = yes 
{{< /code-block >}}

This variable set the use of the GBT and allow to do spin-wave calculations in primitive cell.

{{< tutorial-footer >}}
---
Title: "Periodic Systems"
weight: 40
description: "Tutorial series for systems with periodic boundary conditions"
---


{{< tutorial-series "Periodic Systems" >}}---
title: "BerkeleyGW (2016)"
tags: ["Tutorial", "Bulk", "GW"]
#series: "Tutorial"
hidden: True
---


'''NOTE''': This tutorial page is set up for the [https://www.benasque.org/2016tddft Benasque TDDFT school 2016]. More recent version: [[Tutorial:BerkeleyGW]]

##  Interacting with Cori  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/cori Cori supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainXXX@cori.nersc.gov" >}} in your terminal, substituting the actual name of your training account for {{< code "XXX" >}}.
* Be aware that since this machine is far away, you should not try running X-Windows programs.
* We will work with "debug" jobs which can submitted via the provided scripts in each directory. We are using the SLURM queue manager on Cori.
* To copy files from Cori to your local machine, in a terminal on your local machine, write {{< code "scp trainXXX@cori.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.
* To see if you have jobs running, do {{< code "squeue -u trainXXX" >}}. You should not have more than one in the queue; if you do, cancel them with {{< code "scancel JOBID" >}}, filling in the number for JOBID from the output of {{< code "squeue" >}}.
* Accounts will expire on October 7.
* Cori is actually down for an upgrade since Sept 19, but the files we were using are available if you login to {{< code "edison.nersc.gov" >}} and go to {{< code "$CSCRATCH" >}}.

##  Documentation and resources  

* [https://www.tddft.org/programs/octopus/doc/generated/html/vars.php Octopus variable reference] 
* [https://www.berkeleygw.org/releases/manual_v1.2.0.html BerkeleyGW manual]
* [https://www.benasque.org/2016tddft/talks_contr/154_BerkeleyGW_octopus.pptx.pdf Intro slides from first day]
* [https://www.benasque.org/2016tddft/talks_contr/164_128_Felipe_BSE_Presentation.pdf Intro slides from second day]
* [https://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv
* More extensive [https://sites.google.com/site/berkeleygw2018/about lecture slides] from a longer tutorial devoted solely to BerkeleyGW in January 2018

##  Instructions  

The first time you log in, execute these lines which will help you see color-coding for what is a link, executable, or directory:

```text
 echo 'alias ls="ls --color"' >> ~/.bashrc.ext
 . ~/.bashrc
```

Each time you log in, you should do this:

```text
 - Load modules
 module load berkeleygw/1.2 octopus/6.0
```

```text
 - Go to the scratch directory, where all runs should happen.
 cd $SCRATCH
```

To begin with the examples,

```text
 - List all examples available
 ls /project/projectdirs/mp149/Benasque2016
```

```text
 - Copy 1-silicon example to your directory
 cp -R /project/projectdirs/mp149/Benasque2016/1-boron_nitride .
```

```text
 - Go to your local folder and follow instructions
 cd 1-boron_nitride
 less README
```

##  Schedule  

* Day 1
** 1-boron_nitride, GW
* Day 2
** 2-benzene, GW and [https://en.wikipedia.org/wiki/Hans_Bethe Bethe]-[https://en.wikipedia.org/wiki/Edwin_Ernest_Salpeter Salpeter]
** 3-xct_LiCl, exciton visualization (copy to your local machine according to blackboard or from Cori, unpack the archive, and follow README)
** 4-silicon, Bethe-Salpeter

##  For historical interest  
* [2014 version of this tutorial](../BerkeleyGW (2014) )

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
---
title: "Benasque school"
weight: 100
description: "Basic introduction to run Octopus"
hidden: True
---

---
title: "BerkeleyGW (2014)"
tags: ["Tutorial", "Bulk", "GW"]
#series: "Tutorial"
hidden: True
---


'''NOTE''': This tutorial page was set up for the Benasque TDDFT school 2014. The specific references to the supercomputer used at that time will have to be adapted for others to use this tutorial. More recent version: [[Tutorial:BerkeleyGW]]

##  Interacting with Hopper  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/hopper Hopper supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainX@hopper.nersc.gov" >}} in your terminal, substituting the actual name of your training account.
* Be aware that since this machine is far away, you should not try running X-Windows programs!
* You submit jobs by the {{< code "qsub" >}} command, ''e.g.'' {{< code "qsub job.scr" >}}, which will put them in the queue for execution when there is free space.
* You can see what jobs you currently have in the queue by executing {{< code "qstat -u $USER" >}}, so you can see when your job finishes. A status code will be shown: Q = waiting in the queue, R = running, C = complete.
* You can cancel a job by {{< code "qdel job" >}}, where {{< code "job" >}} is the job id as written by {{< code "qstat" >}}.
* The job script (''e.g.'' {{< code "01-calculate_scf.qsub" >}}) specifies parameters to the PBS/Torque queuing system about how many cores to use, what commands to run, etc.
* To copy files from Hopper to your local machine, in a terminal on your local machine, write {{< code "scp trainX@hopper.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.

##  Getting started in the tutorial  

###  Day 1  
To obtain the files for the boron nitride and benzene examples for the first day of the tutorial:

```text
 cd $SCRATCH
 /project/projectdirs/m1694/BGW-tddft/copy_day_1.sh
```

* In each case, enter your copy of the directory, and look at {{< code "README" >}} and follow instructions given there.
* Start by running {{< code "2-benzene/1-mf/1-scf" >}} and then {{< code "2-benzene/1-mf/2-wfn" >}}. This will take a little while, so while this runs, do the BN example.

###  Day 2  

To obtain the files for the silicon and benzene examples for the second day of the tutorial:

```text
 cd $SCRATCH
 /project/projectdirs/m1694/BGW-tddft/copy_day_2.sh
```

The solution for the benzene example is available at

```text
  /project/projectdirs/m1694/BGW-tddft/2-benzene_run
```

You can copy the necessary files (WFNs, bse*mat, eps*mat, eqp*) from there.

Other instructions:

* We will work on the following directories: 2-benzene and 3-silicon. We will not work on the 1-boron_nitride example!
* Start with the example 3-silicon.
* In each case, enter your copy of the directory, and look at {{< code "README" >}} and follow instructions given there.
* There is additional example for XCrySDen, which is available in the shared folder on imac01 (see instructions on blackboard). If for some reason you are not able to copy it from there, it can also be downloaded [https://civet.berkeley.edu/~jornada/files/xct_LiCl.zip here]. Note: to use XCrySDen on the iMacs, run {{< code "/sw/bin/xcrysden" >}}.

The examples are available for download here: [https://web.mit.edu/~dstrubbe/www/2-benzene.tar.gz 2-benzene.tar.gz], [https://web.mit.edu/~dstrubbe/www/3-silicon.tar.gz 3-silicon.tar.gz].

##  General workflow  

* Finish all basic goals from both the boron nitride and benzene examples before starting any stretch goal.

##  Documentation and resources  

* [https://www.tddft.org/programs/octopus/doc/generated/html/vars.php?page=sections Octopus variable reference] for the pre-release development version. 
* [https://www.berkeleygw.org/releases/manual_v1.0.6.html BerkeleyGW manual].
* [https://benasque.org/2014tddft/talks_contr/115_BerkeleyGW_octopus.pptx.pdf Intro slides from first day]
* [https://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv.
* More extensive [https://www.nersc.gov/users/training/nersc-training-events/berkeleygw2013 lecture slides] from a longer tutorial devoted solely to BerkeleyGW in November 2012.
* The [https://benasque.org/2014tddft/talks_contr/128_Felipe_BSE_Presentation.pdf slides for the second day of tutorial].
Note that we are using the pre-release development version of Octopus in this tutorial, rather than the current release 4.1.2 which lacks full support for BerkeleyGW output. There are some small differences in output from 4.1.2.

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
---
title: "BerkeleyGW"
tags: ["Tutorial", "Bulk", "GW"]
#series: "Tutorial"
hidden: True
---


'''NOTE''': This tutorial page is set up for the [https://www.benasque.org/2018tddft Benasque TDDFT school 2018].

* Your instructors: David Strubbe and Adriel González.

* Test your knowledge with a [https://faculty.ucmerced.edu/dstrubbe/BerkeleyGW_GW_BSE_quiz.pdf quiz]! 

##  Interacting with Cori  

The [[BerkeleyGW]] tutorial is done on the [https://www.nersc.gov/users/computational-systems/cori Cori supercomputer] (at NERSC in California). There are a few key things you need to know about how to interact with the machine:
* To log in, run {{< code "ssh trainXXX@cori.nersc.gov" >}} in your terminal, substituting the actual name of your training account for {{< code "XXX" >}}.
* Be aware that since this machine is far away, running X-Windows programs will be very slow.
* For visualization with XCrySDen or other tools, installing and running on your laptop is recommended.
* We submit jobs using the SLURM queue manager, using a "reservation" for Aug 22 (benasque2018_1), 23 (benasque2018_2), and 24 (benasque2018_3). The reservations are usable from 15:00 to 19:00 each day. If you want to run outside of these times, use the debug queue. The tutorials are designed to be run using an "interactive job", in which you have a node continually available for your use and then executables launched from the command line are run immediately. Start an interactive job like this (do not submit more than one):
```text
 salloc -N 1 -q regular --reservation=benasque2018_2 -t 03:00:00 -C haswell
```
If you are using the debug queue rather than a reservation, use this:
```text
 salloc -N 1 -q debug -t 00:30:00 -C haswell
```
<!-- To use a reservation, you need to add the --reservation=<reservation name> flag into your submission. So for example for the second day you could either put the following line in your batch script:
```text
 -SBATCH --reservation=benasque2018_2
```
or add the flag on the command line when you submit your script:
```text
 sbatch --reservation=benasque2018_2 ./myscript.sl
```
It also works with interactive jobs:
```text
 salloc --reservation=benasque2018_2 <other arguments> -->
```

* Here is an example script using 32 cores. If you are using an interactive job, don't use this.

```text
 -!/bin/bash
 
 -SBATCH -J test_pulpo
 -SBATCH -N 1
 -SBATCH -C haswell
 -SBATCH -p debug
 -SBATCH -t 00:30:00
 -SBATCH --export=ALL
 
 module load octopus/8.2
 srun -n 32 octopus &> output
```

* To copy files from Cori to your local machine, in a terminal on your local machine, write {{< code "scp trainXXX@cori.nersc.gov:FULL_PATH_TO_YOUR_FILE ." >}} (filling in the username and filename) and enter your password when prompted. For very small ASCII files, you may find cut and paste more convenient.
* To see if you have jobs running, do {{< code "squeue -u trainXXX" >}}. You should not have more than one in the queue; if you do, cancel them with {{< code "scancel JOBID" >}}, filling in the number for JOBID from the output of {{< code "squeue" >}}.
* Accounts will expire on September 11. Feel free to copy the files off the machine before that to somewhere else for your future reference.

##  Documentation and resources  

* [https://www.tddft.org/programs/octopus/doc/generated/html/vars.php Octopus variable reference] 
* [https://oldsite.berkeleygw.org/releases/manual_v1.2.0.html BerkeleyGW manual for version 1.2.0]
* [https://www.benasque.org/2018tddft/talks_contr/232_BerkeleyGW_octopus_2018.pdf Lecture: Practical calculations with the GW approximation and Bethe-Salpeter equation in BerkeleyGW]
* [https://faculty.ucmerced.edu/dstrubbe/BerkeleyGW_BSE_2018.pdf Lecture: Practical BSE Calculations with BerkeleyGW + Octopus]
* [https://arxiv.org/abs/1111.4429 BerkeleyGW implementation paper] on arxiv
* More extensive [https://sites.google.com/site/berkeleygw2018/about lecture slides and other examples] from a longer tutorial devoted solely to BerkeleyGW in January 2018
* [https://drive.google.com/file/d/1mLCZrvpZrwsfsnvSBfZs0aKpKFWvSZIU/view?usp=drive_web slides on developing an interface for a new DFT code]

##  Instructions  

The first time you log in, execute these lines which will help you see color-coding for what is a link, executable, or directory:

```text
 echo 'alias ls="ls --color"' >> ~/.bashrc.ext
 . ~/.bashrc
```

Each time you log in, you should do this:

```text
 - Load modules
 module load berkeleygw/1.2
```

```text
 - Go to the scratch directory, where all runs should happen.
 cd $SCRATCH
```

To begin with the examples,

```text
 - List all examples available
 ls /project/projectdirs/mp149/Benasque2018
```

```text
 - Copy 1-boron_nitride example to your directory
 cp -R /project/projectdirs/mp149/Benasque2018/1-boron_nitride .
```

```text
 - Go to your local folder and follow instructions
 cd 1-boron_nitride
 less README
```

##  Schedule  

* Day 1
** 1-boron_nitride, GW
** 2-benzene, GW
* Day 2
** 2-benzene, Bethe-Salpeter
** 3-xct_LiCl, exciton visualization (download [https://faculty.ucmerced.edu/dstrubbe/xctLiCl.zip here] or from Cori, unpack the archive with "unzip", and follow README)
** 4-silicon, Bethe-Salpeter
* Day 3 (optional)
** Stretch goals from 2-benzene
** Any of the plane-wave BerkeleyGW examples, on Cori at {{< code "/project/projectdirs/mp149/BGW-2018/" >}}: silicon with PARATEC; boron nitride sheet, benzene, or sodium with Quantum ESPRESSO.

##  Errata  

* In {{< code "1-boron_nitride/1-mf/2.1-epsilon" >}}, the reference to {{< code "./01-calculate_wfn.qsub" >}} should be {{< code "./02-calculate_wfn.qsub" >}}, and you should run {{< code "01-get_kgrid.sh" >}} first and look at the input and output files.
* In the output of {{< code "kgrid.x" >}}, in the {{< code "KPointsGrid" >}} block, the value for a small q-shift is not quite right, and should be divided by the number of k-points in that direction. The value given in the Octopus input files is the correct one.
* In {{< code "1-boron_nitride/2-bgw/1-epsilon" >}}, the script {{< code "./01-run_epsilon.qsub" >}} should have "32" rather than "48" in the execution line, or you will get some OMP errors.
* In {{< code "1-boron_nitride/2-bgw/2-sigma" >}}, the script {{< code "./01-run_sigma.qsub" >}} should have "32" rather than "48" in the execution line, or you will get some OMP errors.

##  For historical interest  
* [2014 version of this tutorial](../BerkeleyGW (2014) )
* [2016 version of this tutorial](../BerkeleyGW (2016) )

<span class=noprint><hr>
Back to [[Tutorials]]




---------------------------------------------
---
title: "Benasque"
#series: "Tutorial"
hidden: true
---


This is the tutorial of Octopus for the Benasque TDDFT school 2008. 

If you need help you can check the {{< manual "" "manual" >}}. You can also check the documentation for each variable, the easiest way to do it is to run {{< command "oct-help show <variable_name>" >}} from the command line. You can also search among variable names by running {{< command "oct-help search <string>" >}}.

* Day 1: The ground state
** [Hydrogen atom](../Hydrogen atom) - getting started
** [Nitrogen atom](../Nitrogen atom) - basic input variables
** [Methane molecule](../Methane molecule) - converging a ground-state calculation
** [Benzene molecule](../Benzene molecule) - making 3D plots
* Day 2: Excited States
** [Time-dependent run](../Time-dependent run)
** [Optical spectra from time-propagation](../Optical Spectra from TD) - how to obtain the absorption spectrum through the explicit solution of the time-dependent Kohn-Sham equations
** [Optical spectra from Casida's equation](../Optical Spectra from Casida) - how to solve Casida's equation to get an optical spectrum
* Optional
** [Basic QOCT](../Basic QOCT) - getting started with quantum optimal control theory.
---------------------------------------------
---
title: "OpenDX"
tags: ["Obsolete"]
#series: "Tutorial"
---


{{< notice warning >}}
The OpenDX program is no longer supported and the following tutorial should be considered obsolete. We recommend you to use one of the alternatives mentioned in this {{< tutorial "Basics:Visualization" "Visualization" >}}.
{{< /notice >}}
## Input

At this point, the input should be quite familiar to you:

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_Angstrom
 
 {{< variable "Radius" >}} = 5*angstrom
 {{< variable "Spacing" >}} = 0.15*angstrom
 
 {{< variable "Output" >}} = wfs + density + elf + potential
 {{< variable "OutputFormat" >}} = cube + xcrysden + dx + axis_x + plane_z
 
 {{< variable "XYZCoordinates" >}} = "benzene.xyz"
{{< /code-block >}}

Coordinates are in this case given in the file {{< file "benzene.xyz" >}} file. This file should look like:

{{< code-block >}}
 12
    Geometry of benzene (in Angstrom)
 C  0.000  1.396  0.000
 C  1.209  0.698  0.000
 C  1.209 -0.698  0.000
 C  0.000 -1.396  0.000
 C -1.209 -0.698  0.000
 C -1.209  0.698  0.000
 H  0.000  2.479  0.000
 H  2.147  1.240  0.000
 H  2.147 -1.240  0.000
 H  0.000 -2.479  0.000
 H -2.147 -1.240  0.000
 H -2.147  1.240  0.000
{{< /code-block >}}

Note that we asked octopus to output the wavefunctions ({{< code wfs >}}), the density, the electron localization function ({{< code elf >}}) and the Kohn-Sham potential. We ask Octopus to generate this output in several formats, that are requires for the different visualization tools that we mention below. If you want to save some disk space you can just keep the option that corresponds to the program you will use. Take a look at the documentation on variables {{< variable "Output" >}} and {{< variable "OutputFormat" >}} for the full list of possible quantities to output and formats to visualize them in.

### OpenDX

[OpenDX](https://www.opendx.org) is a very powerful, although quite complex, program for the 3D visualization of data. It is completely general-purpose, and uses a visual language to process the data to be plotted. Don't worry, as you will probably not need to learn this language. We have provided a program ({{< file "mf.net" >}} and {{< file "mf.cfg" >}}) that you can get from the {{< file "SHARE/util" >}} directory in the {{< octopus >}} installation (probably it is located in {{< file "/usr/share/octopus/util" >}}). 

Before starting, please make sure you have it installed together with the Chemistry extensions. You may find some instructions [[Releases-OpenDX|here]]. Note that openDX is far from trivial to install and to use. However, the outcomes are really spectacular, so if you are willing to suffer a little at the beginning, we are sure that your efforts will be compensated.

If you run {{< octopus >}} with this input file, afterward you will find in the {{< file "static" >}} directory several files with the extension {{< file ".dx" >}}. These files contain the wave-functions, the density, etc. in a format that can be understood by [https://www.opendx.org/ Open DX] (you will need the chemistry extensions package installed on top of DX, you can find it in our [[releases]] page). 

So, let us start. First copy the files {{< file "mf.net" >}} and {{< file "mf.cfg" >}} to your open directory, and open dx with the command line

```text
  > dx
```

If you followed the instructions given in [[Releases-OpenDX|here]] you should see a window similar to the one shown on the right.
<gallery>
```text
  Image:Tutorial_dx_1.png|Screenshot of the Open DX main menu
  Image:Tutorial_dx_2.png|Main dialog of mf.net
  Image:Tutorial_dx_3.png|Isosurfaces dialog
  Image:Tutorial_dx_4.png|Geometry dialog
```
</gallery>

Now select '''<tt>Run Visual Programs...</tt>''', select the {{< file "mf.net" >}} file, and click on '''<tt>OK</tt>'''. This will open the main dialog of our application. You will no longer need the main DX window, so it is perhaps better to minimize it to get it out of the way. Before continuing, go to the menu '''<tt>Execute > Execute on Change</tt>''', which instructs to open DX to recalculate out plot whenever we make a change in any parameters.

{{< figure src="/images/Tutorial_dx_5.png" width="500px" caption="350px" >}}

Now, let us choose a file to plot. As you can see, the field <tt>File</tt> in the mains dialog is <tt>NULL</tt>. If you click on the three dots next to the field a file selector will open. Just select the file {{< file "static/density-1.dx" >}} and press OK. Now we will instruct open DX to draw a iso-surface. Click on <tt>Isosurfaces</tt> and a new dialog box will open (to close that dialog box, please click again on <tt>Isosurfaces</tt> in the main dialog!) In the end you may have lots of windows lying around, so try to be organized. Then click on <tt>Show Isosurface 1</tt> - an image should appear with the density of benzene. You can zoom, rotate, etc the picture if you follow the menu '''<tt>Options > View control...</tt>''' in the menu of the image panel. You can play also around in the isosurfaces dialog with the <tt>Isosurface value</tt> in order to make the picture more appealing, change the color, the opacity, etc. You can even plot two isosurfaces with different values.

We hope you are not too disappointed with the picture you just plotted. Densities are usually quite boring (even if the Hohenberg-Kohn theorem proves that they know everything)! You can try to plot other quantities like the electronic localization function, the wave-functions, etc, just by choosing the appropriate file in the main dialog. You can try it later: for now we will stick to the "boring density".

Let us now add a structure. Close the isosurfaces dialog '''by clicking again in <tt>Isosurfaces</tt> in the main dialog box''', and open the <tt>Geometry</tt> dialog. In the <tt>File</tt> select the file {{< file "benzene.xyz" >}}, and click on <tt>Show Geometry</tt>. The geometry should appear beneath the isosurface. You can make bonds appear and disappear by using the slide <tt>Bonds length</tt>,  change its colors, etc. in this dialog.

To finish, let us put some "color" in our picture. Close the <tt>Geometry</tt> dialog (you know how), and open again the isosurfaces dialog. Choose <tt>Map isosurface</tt> and in the <tt>File</tt> field choose the file {{< file "static/vh.dx" >}}. There you go, you should have a picture of an isosurface of benzene colored according to the value of the electrostatic potential! Cool, no? If everything went right, you should end with a picture like the one on the right.

As an exercise you can try to plot the wave-functions, etc, make slabs, etc. Play around. As usual some things do not work perfectly, but most of it does!


---------------------------------------------
---
title: "Open boundaries"
tags: ["Obsolete"]
#series: "Tutorial"
---


{{< notice warning >}}
'''WARNING: The open boundaries implementation was removed in version 5.0.0 of {{< octopus >}} and is no longer available.'''
{{< /notice >}}


This tutorial briefly sketches how to perform calculations with open (or transparent) boundaries in {{< octopus >}}.
{{< emph "Important:" >}} Please note that the features described in the following are still subject to changes and further development and certainly not as stable and reliable as required for production runs. Therefore, we very much appreciate feedback and bug reports (or even fixes).

## Geometry 

Open boundary calculations are only possible for the parallelepiped {{< variable "BoxShape" >}}, and only the planes perpendicular to the x-axis are transparent, ''i. e.'' the transport direction is along the ''x''-axis, all other boundaries are rigid walls like for the usual finite-system calculations. As the system is open at the "left" and "right" side we have to specify the character of the outside world here. This is done by giving a unit cell of two semi-infinite leads attached to the central or device region, yielding the following geometry (in 2D):

{{< figure src="/images/Open_boundaries_geometry.png" width="500px" caption="Open boundaries geometry" >}}

The red area is the device region also constituting the simulation box, ''i. e.'' the only part in space we treat explicitly in our simulations, which is attached to semi-infinite periodic leads (shown in green). The blue lines indicate the left and right open boundaries.

## Overview of the method 

Our approach to open-boundary calculations in {{< octopus >}} generally consists of a three-step procedure:

- A periodic ground-state calculation of the lead unit cell provides the eigenstates of the leads.
- We now scatter these "free" eigenstates in the device region to obtain the ground-state wave functions of the system with broken translational symmetry. This is done by solving a Lippmann-Schwinger equation for the central region. Note: this implies symmetric leads.
- We are now in a position to propagate the ground state of step 2 in time. The system may be perturbed by time-dependent electric and magnetic fields in the central region &ndash; just as TDDFT for finite systems &ndash; and, additionally, by different left and right time-dependent but spatially constant potentials in the leads (time-dependent bias). The propagation is performed by a modified Crank-Nicolson propagator that explicitly includes the injection of density from the leads and the hopping in and out of the central region by additional source and memory terms. For details, see 
[^kurth2005transport]


## Obtaining the ground state 

Calculating the ground state for the above geometry requires two runs of {{< octopus >}}:

- A ground-state calculation of an infinite chain of unit cells of the lead.
- A ground-state calculation of the full system, ''i. e.'' with the scattering center in place.

We use the multi-dataset capabilities of {{< octopus >}} to connect the two.

Our first example will calculate the scattering state of an attractive 1D square potential:

{{< figure src="/images/Open_boundaries_square_well.png" width="500px" caption="1D square well potential" >}}

We begin with the calculation of the leads which are, in this case, flat and of zero potential, hence, we know that our result should be plane waves. We write the following input file in which we name our periodic dataset {{< code "lead_" >}} in order to reference it in subsequent calculations:

```text

%CalculationMode
 "lead_"
 gs
%

FromScratch = yes
TheoryLevel = independent_particles

BoxShape = parallelepiped
Dimensions = 1
Spacing = 0.1

%Species
 "flat" | 0 | user_defined | 2.0 | "0"
%


- The lead.
lead_PeriodicDimensions = 1
%lead_Coordinates
 "flat" | 0 | 0 | 0
%

%lead_Lsize
 40*Spacing | 0 | 0
%

%lead_KPoints
 1 | 0.1 | 0 | 0
%

lead_Output = potential
lead_OutputHow = binary
```


We write out the potential in binary format because it is needed in the next run, where we will scatter the outcome of the periodic run at the square well. The following additions to the the input file will do the job:

```text

%CalculationMode
 "lead_" | "well_"
 gs      | gs
%

<...>

V = 0.5
W = 10
%Species
 "flat" | 0 | user_defined | 2.0 | "0"
 "well" | 0 | user_defined | 2.0 | "V*(-step(x+W/2)+step(x-W/2))"
%

<...>

- The extended system with scattering center.
add_ucells = 2
%well_OpenBoundaries
 lead_dataset     | "lead_"
 lead_restart_dir | "lead_restart"
 lead_static_dir  | "lead_static"
 add_unit_cells   | add_ucells
%

%well_Coordinates
 "well" | 0 | 0 | 0
%

%well_Lsize
 20 | 0 | 0
%
```

A few points deserve comments:

* The block {{< variable "OpenBoundaries" >}} enables the entire open-boundaries machinery. In the block, we have to specify which dataset provided the periodic results. Note that only the {{< code "lead_" >}} run is periodic by virtue of <tt>lead_PeriodicDimensions=1</tt>. Furthermore, the algorithm that calculates the ground state of the extended system needs to know where to find the lead potential (line starting with {{< code "lead_static_dir" >}}) and the ground-state wave functions of the lead (line starting with {{< code "lead_restart_dir" >}}).

* The line {{< code "add_unit_cells" >}} is a bit special: it is usually desirable to include the first few unit cells of the lead in the simulation box. To avoid having to put them there by hand {{< octopus >}} attaches as many unit cells as requested in this line to the simulation box of the dataset with open boundaries, {{< code "well_" >}} in this case. So keep in mind that this line may significantly increase your grid sizes.

Processing this input file with {{< octopus >}} gives us the following eigenstate ($|\psi(x)|^2$ is plotted):

{{< figure src="/images/Open_boundaries_square_well_eigenstate.png" width="500px" caption="Extended eigenstate of 1D square well" >}}

Before we turn our attention to a time-propagation, we briefly describe the console output of {{< octopus >}} related to open boundaries:

* In the <tt>Grid</tt> section appear two entries. The first one stating that open boundaries are indeed enabled showing the number of unit cells included in the simulation box. The alert reader will note that it says <tt>3</tt> here but <tt>2</tt> in the input file. The reason is that at least one additional unit cell is required by the propagating algorithm into which we will eventually feed our ground-state wave functions. The second one prints how many grid points are in the interface region which is basically determined by the size of the leads' unit cell.
```text

******************************** Grid ********************************
<...>
Open boundaries in x-direction:
  Additional unit cells left:       3
  Additional unit cells right:      3
  Left lead read from directory:  lead_restart
  Right lead read from directory: lead_restart
<...>
Number of points in left  interface:     80
Number of points in right interface:     80
**********************************************************************
```

* The following section is also related to open boundaries:
```text

************************ Lead Green functions ************************
 st-  Spin  Lead     Energy
   1   --     L     0.005000
   1   --     R     0.005000
**********************************************************************
```

The Green's functions of the leads are required to include the effects of the leads in the Lippmann-Schwinger equation for the simulation region. More information can be found in this [https://www.physik.fu-berlin.de/~lorenzen/physics/transport_nq_yrm08.pdf poster] which was presented at the 5th Nanoquanta Young Researchers Meeting 2008.

## Time propagation 

We give two small examples for time propagations.

{{< notice note >}}
TODO: Add examples with TD potentials in the central region and in the leads. Give an example of a calculation of an I-V-curve (as soon as current calculation is fixed).
{{< /notice >}}

### Propagating an eigenstate 

As a first time-dependent simulation we simply propagate the eigenstate obtained in the previous section in time without any external influences. We expect that $|\psi(x)|^2$ will not change and that $\psi(x)$ picks up a phase which is shown in this movie:

<center>
<flash>file=Open_boundaries_square_well_eigenstate_propagation.swf|width=720|height=504|quality=best</flash>
</center>

In order to obtain the input data for the movie we add a third dataset to our input file starting a <tt>td</tt> calculation:

```text

%CalculationMode
 gs      | gs      | td
 "lead_" | "well_" | "well_"
 1       | 2       | 3
%
```


Note that this dataset shares label with the ground state run. This enables {{< octopus >}} to find the correct initial state and avoids the need to repeat the <tt>OpenBoundaries</tt> block. Apart from the <tt>td</tt> dataset we only have to give the usual TD parameters: timestep and number of steps in the propagation.

```text

TDMaximumIter = 300
TDTimeStep = 0.5
```


For completeness, we mention some of {{< octopus >}}' output related to open boundaries in TD propagations:

```text

************************** Open Boundaries ***************************
Type of memory coefficients:           full
MBytes required for memory term:     59.326
Maximum BiCG iterations:                100
BiCG residual tolerance:            1.0E-12
Included additional terms:              memory source
TD lead potential (L/R):                  0         0
**********************************************************************
```


This block gives some details about the parameters of the propagation algorithm:

* Both source and memory terms will be included in the simulation (this can be changed by the {{< variable "OpenBoundariesAdditionalTerms" >}} variable).
* The memory coefficients are stored as full matrices (by {{< variable "OpenBoundariesMemType" >}}) and require 59 MB of main memory.
* The time-dependent potential of the leads is 0.
* The maximum number of iterations solving the Crank-Nicolson system of linear equations is 100 (by {{< variable "OpenBoundariesBiCGMaxiIter" >}}) and the tolerance for the iterative solver ({{< variable "OpenBoundariesBiCGTol" >}}) is $10^{-12}$.

The calculation of the memory coefficients is accompanied by the oputput

```text

Info: Calculating missing coefficients for memory term of left lead.
[ 52/301]  17%|****                      |     05:06 ETA
```


and takes very long. Indeed, this is one of the bottlenecks in the algorithm.


### Propagation of a Gaussian wavepacket 

Our second example is the propagation of a Gaussian-shaped charge distribution in 2D that is initially completely localized in the central region. This calculation is the "standard" test for all implementations of open boundaries. Therefore, it is worth repeating it here.
As there are no charges outside the central region we have to disable the source term which is repsonsible for injection of density into our simulation region. Furthermore, the initial state of our time propagation is now given as an entry in the {{< variable "UserDefinedStates" >}} block and not as the result of a ground state calculation. Nevertheless, the two ground-state calculations for a completely flat 2D strip have to be performed to initialize certain datastructures.

The part of the input describing the flat 2D strip looks like this:

```text

%CalculationMode
 gs      | gs
 "lead_" | "flat_"
%

FromScratch = yes

TheoryLevel = independent_particles

Dimensions = 2
BoxShape = parallelepiped
Spacing = 0.25

%Species
 "flat" | 0 | user_defined | 2.0 | "0"
%

- The lead.
lead_PeriodicDimensions = 1
%lead_Coordinates
 "flat" | 0 | 0 | 0
%
%lead_Lsize
 4*Spacing | 4 | 0
%
%lead_KPoints
 1 | 0.25 | 0 | 0
%
lead_Output = potential
lead_OutputHow = binary

- The central region.
%flat_OpenBoundaries
 lead_dataset     | "lead_"
 lead_restart_dir | "lead_restart"
 lead_static_dir  | "lead_static"
 add_unit_cells   | 0
%
%flat_Coordinates
 "flat" | 0 | 0 | 0
%
%flat_Lsize
 4 | 4 | 0
%
```


For the propagation we now override the initial state by giving an explicit formula of a 2D Gaussian wavepacket with initial momentum $(1.75, 1)$

```text

i = {0, 1}
kx = 1.75
ky = 1.0
alpha = 0.5
i = {0, 1}
flat_OnlyUserDefinedInitialStates = yes
%flat_UserDefinedStates
 1 | 1 | 1 | formula | "exp(i*(kx*x + ky*y))*exp(-alpha*r*r)"
%
```

and add the TD dataset:

```text

%CalculationMode
 gs      | gs      | td
 "lead_" | "flat_" | "flat_"
%
```


We disable the source term by the line

```text

flat_OpenBoundariesAdditionalTerms = mem_term
```


The default of {{< variable "OpenBoundariesAdditonalTerms" >}} is <tt>mem_term + src_term</tt>, so solely giving <tt>mem_term</tt> disables the source term.

Feeding the input into {{< octopus >}} and generating a movie of $|\psi(x, y, t)|^2$ gives

<center>
<flash>file=Open_boundaries_wavepacket_propagation.swf|width=720|height=504|quality=best</flash>
</center>

Since the wavepacket has an initial momentum in the $x$- as well as in $y$-directions, the open boundary at the right and the rigid wall at the back are easily recognized.

## References 

<references/>


---------------------------------------------
[^kurth2005transport]: {{< article title="Time-dependent quantum transport: A practical scheme using density functional theory" authors="S. Kurth, G. Stefanucci, C.-O. Almbladh, A. Rubio, E. K. U. Gross" journal="Phys. Rev. B" volume="72" pages="35308" year="2005" url="https://journals.aps.org/prb/abstract/10.1103/PhysRevB.72.035308" doi="10.1103/PhysRevB.72.035308" link="APS" >}}

---
title: "Obsolete Material"
weight: 110
description: "Basic introduction to run Octopus"
hidden: True
---

---
title: "1D Helium"
tags: ["Beginner", "Ground State", "Unoccupied", "Model", "User-defined Species", "Independent Particles", "Visualization"]
weight: 3
tutorials: "Model Systems"
difficulties: "beginner"
theories: "Independent particles"
calculation_modes: ["Ground state", "Unoccupied"]
system_types: "Model"
species_types: "User-defined species"
features: "Visualization"
description: "The helium atom in one dimension which also has two electrons."
#series: "Tutorial"
---


The next example will be the helium atom in one dimension which also has two electrons, just as we used for the harmonic oscillator. The main difference is that instead of describing two electrons in one dimension we will describe one electron in two dimensions. The calculation in this case is not a DFT one, but an exact solution of the Schrödinger equation -- not an exact solution of the helium atom, however, since it is a one-dimensional model.

### Equivalence between two 1D electrons and one 2D electron

To show that we can treat two electrons in one dimension as one electron in two dimensions, lets start by calling $x\,$ and $y\,$ the coordinates of the two electrons. The Hamiltonian would be (note that the usual Coulomb interaction between particles is usually substituted, in 1D models, by the *soft Coulomb* potential, 
$u(x)=(1+x^2)^{(-1/2)}\\,$):

$$
  \\hat{H} = -\\frac{1}{2}\\frac{\\partial^2}{\\partial x^2}
            -\\frac{1}{2}\\frac{\\partial^2}{\\partial y^2}
  +\\frac{-2}{\\sqrt{1+x^2}}+\\frac{-2}{\\sqrt{1+y^2}}+\\frac{1}{\\sqrt{1+(x-y)^2}}.
$$

Instead of describing two electrons in one dimension, however, we may very well think of one electron in two dimensions,
subject to a external potential with precisely the shape given by:

$$
  \\frac{-2}{\\sqrt{1+x^2}}+\\frac{-2}{\\sqrt{1+y^2}}+\\frac{1}{\\sqrt{1+(x-y)^2}}
  \\,,
$$

Since the Hamiltonian is identical, we will get the same result. Whether we regard $x\,$ and $y\,$ as the coordinates of two different particles in one dimension or as the coordinates of the same particle along the two axes in two dimensions is entirely up to us. (This idea can actually be generalized to treat two 2D particles via a 4D simulation in {{< octopus >}} too!) Since it is usually easier to treat only one particle, we will solve the one-dimensional helium atom in two dimensions. We will also therefore get a "two-dimensional wave-function". In order to plot this wave-function we specify an output plane instead of an axis.

### Input

With the different potential and one more dimension the new input file looks like the following

{{< code-block >}}
#include_input doc/tutorials/model_systems/1D_Helium/1.gs/inp
{{< /code-block >}}

For more information on how to write a potential formula expression, see {{< manual "Basics:Input file" "Input file" >}}.

We named the species "helium" instead of "He" because "He" is already the name of a pseudopotential for the actual 3D helium atom.

### Running

The calculation should converge within 14 iterations. As usual, the results are summarized in the {{< file "static/info" >}} file, where you can find

{{< code-block >}}
 ...
**************************** Theory Level ****************************
Input: [TheoryLevel = independent_particles]
**********************************************************************

SCF converged in   14 iterations

Some of the states are not fully converged!
Eigenvalues [H]
 #st  Spin   Eigenvalue      Occupation
   1   --    -2.238257       1.000000

Energy [H]:
      Total       =        -2.23825730
      Free        =        -2.23825730
      -----------
      Ion-ion     =         0.00000000
      Eigenvalues =        -2.23825730
      Hartree     =         0.00000000
      Int[n*v_xc] =         0.00000000
      Exchange    =         0.00000000
      Correlation =         0.00000000
      vanderWaals =         0.00000000
      Delta XC    =         0.00000000
      Entropy     =         1.38629436
      -TS         =        -0.00000000
      Kinetic     =         0.28972647
      External    =        -2.52798377
      Non-local   =         0.00000000

Dipole:                 [b]          [Debye]
      <x> =   -2.82217E-07     -7.17323E-07
      <y> =    1.26162E-07      3.20671E-07

Convergence:
      abs_dens =  9.81658651E-06 ( 0.00000000E+00)
      rel_dens =  9.81658651E-06 ( 1.00000000E-05)
      abs_ev =  6.00463679E-10 ( 0.00000000E+00) [H]
{{< /code-block >}}

As we are running with non-interacting electrons, the Hartree, exchange and correlation components of the energy are zero. Also the ion-ion term is zero, as we only have one "ion".

###  Unoccupied States  

Now you can do just the same thing we did for the {{< tutorial "Model:1D_Harmonic_Oscillator" "harmonic oscillator">}} and change the unoccupied calculation mode:

{{< code-block >}}
#include_input doc/tutorials/model_systems/1D_Helium/2.unocc/inp
{{< /code-block >}}


We have added extra states and also restricted the wavefunctions to plot ({{< code-inline >}}{{< variable "OutputWfsNumber" >}} = "1-4,6"{{< /code-inline >}}).

The results of this calculation can be found in the file {{< file "static/eigenvalues" >}}. In this case it looks like

{{< code-block >}}
Some of the states are not fully converged!
Criterion =      0.100000E-05

Eigenvalues [H]
 #st  Spin   Eigenvalue      Occupation     Error
   1   --    -2.238257       1.000000      (2.9E-06)
   2   --    -1.815718       0.000000      (7.9E-07)
   3   --    -1.701549       0.000000      (9.7E-07)
   4   --    -1.629240       0.000000      (9.6E-07)
   5   --    -1.608656       0.000000      (9.3E-07)
   6   --    -1.509599       0.000000      (4.1E-07)
{{< /code-block >}}

Apart from the eigenvalues and occupation numbers we asked {{< octopus >}} to output the wave-functions. To plot them, we will use gnuplot. You start it and type

```bash
#include_file doc/tutorials/model_systems/1D_Helium/2.unocc/plot.gp
```

We plot the ground-state, 1st and 2nd excited-state wave-functions. (If you get this, ignore it: {{< code-inline >}} warning: Cannot contour non grid data. Please use "set dgrid3d".{{< /code-inline >}}) Which correspond to singlet and which to triplet states?

#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0001.eps caption="Ground-state of He in 1D" 
#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0002.eps caption="1st excited-state of He in 1D" 
#include_eps doc/tutorials/model_systems/1D_Helium/2.unocc/Wfs-st0003.eps caption="2nd excited-state of He in 1D" 


### Exercises

* Calculate the helium atom in 1D, assuming that the 2 electrons of helium do not interact (using {{< variable "Dimensions" >}} = 1). Can you justify the differences?

* See how the results change when you change the interaction. Often one models the Coulomb interaction by $1/\sqrt{a^2+r^2}\,$, and fits the parameter $a\,$ to reproduce some required property.

{{< tutorial-footer >}}
---
title: "Model Systems"
weight: 30
description: "Model Systems"
---


{{<octopus>}} has the unusual feature for DFT codes of being able to handle "model systems," i.e. ones that have a user-defined arbitrary potential as opposed to a system of real atoms. This can be useful for simplified test calculations or for computations in the "effective mass" or "envelope function" approximation, e.g. for quantum dots or quantum wells. In this series of linked tutorials, we will show how to run 1D and 2D systems, how to set up a user-defined arbitrary potential, and showcase some more unusual features of {{<octopus>}}. These tutorials assume that you are already familiar with using {{<octopus>}}. If that is not the case, start first with the 
{{< tutorial "basics" "Octopus basics" >}} tutorials.

{{< tutorial-series "Model Systems" >}}
---
title: "Particle in a box"
tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles"]
weight: 2
tutorials: "Model Systems"
difficulties: "beginner"
theories: "Independent particles"
calculation_modes: "Ground state"
system_types: "Model"
species_types: "User-defined species"
description: "The classic quantum problem of a particle in a box."
#series: "Tutorial"
---


In the previous tutorial, we considered applying a user-defined potential. What if we wanted to do the classic quantum problem of a particle in a box, ''i.e.'' an infinite square well?

$$
V(x) = 
\begin{cases}
0, & \frac{L}{2} < x < \frac{L}{2} \cr
\infty, & \text{otherwise}
\end{cases}
$$

There is no meaningful way to set an infinite value of the potential for a numerical calculation. However, we can instead use the boundary conditions to set up this problem. In the locations where the potential is infinite, the wavefunctions must be zero. Therefore, it is equivalent to solve for an electron in the potential above in all space, or to solve for an electron just in the domain $x \in [-\tfrac{L}{2}, \tfrac{L}{2}]$ with zero boundary conditions on the edges. In the following input file, we can accomplish this by setting the "radius" to $\tfrac{L}{2}$, for the default box shape of "sphere" which means a line in 1D.

## Input
As usual, we will need to write an input file describing the system we want to calculate:

{{< code-block >}}
#include_input doc/tutorials/model_systems/particle_in_a_box/inp
{{< /code-block >}}

Run this input file and look at the ground-state energy and the eigenvalue of the single state.

## Exercises

* Calculate unoccupied states and check that they obey the expected relation.
* Vary the box size and check that the energy has the correct dependence.
* Plot the wavefunctions and compare to your expectation.
* Set up a calculation of a ''finite'' square well and compare results to the infinite one as a function of potential step. (Hint: along with the usual arithmetic operations, you may also use logical operators to create a piecewise-defined expression. See {{< manual "Input file" "Input file" >}}).
* Try a 2D or 3D infinite square well.


{{< tutorial-footer >}}
---
title: "Kronig-Penney Model"
#tags: ["Beginner", "Ground State", "Model", "Chain", "Band Structure", "User-defined Species", "Independent Particles"]
section: "/home/lueders/Downloads/Tutorial Kronig-Penney_Model"
tutorials: "Model Systems"
theories: "Independent particles"
calculation_modes: "Ground state"
system_types: "Chain"
species_types: "User-defined species"
features: "Band structure"
difficulties: "beginner"
description: "Calculate the bandstructure for Kronig-Penney Model."
---


The Kronig-Penney model is a 1D system that demonstrates band gaps, which relate to the allowed energies for electrons in a material. In this tutorial we calculate the bandstructure for Kronig-Penney Model. The Kronig-Penney Model has a periodic potential of

$$
V(x) =
\begin{cases}
      V_0 & -b < x < 0 \cr
      0 & 0 < x < a
\end{cases}
$$

Where b is the width of each barrier, and a is the spacing between them. 

## Input
The following input file will be used for the ground state calculation:

{{< code-block >}}
#include_input doc/tutorials/model_systems/kronig_penney_model/1.gs/inp
{{< /code-block >}}

{{< figure src="/images/Kp_wavefunctions.png" width="500px" caption="The first two wavefunctions plotted alongside the potential." >}}


## Bandstructure
To calculate the bandstructure simply change the {{< variable "CalculationMode" >}} to unocc.

{{< code-block >}}
#include_input doc/tutorials/model_systems/kronig_penney_model/2.unocc/inp
{{< /code-block >}}


#include_eps doc/tutorials/model_systems/kronig_penney_model/2.unocc/bandstructure.eps caption="The band structure for Kronig-Penney Model."


To plot the bandstructure, we will use the same command from the {{< tutorial "Periodic_Systems" "Periodic Systems" >}} (assuming you are using gnuplot).

{{< code-block >}}
 plot for [col=5:5+9] 'static/bandstructure' u 1:(column(col)) w l notitle ls 1
{{< /code-block >}}

Reference:
Sidebottom DL. Fundamentals of condensed matter and crystalline physics: an introduction for students of physics and materials science. New York: Cambridge University Press; 2012.



{{< tutorial-footer >}}





---
title: "1D Harmonic Oscillator"
tags:: ["Beginner", "Ground State", "Unoccupied", "Model", "User-defined Species", "Independent Particles"]
weight: 1
tutorials: "Model Systems"
difficulties: "beginner"
theories: "Independent particles"
calculation_modes: ["Ground state", "Unoccupied"]
species_types: "User-defined species"
description: "The standard textbook harmonic oscillator in one dimension with two non-interacting electrons."
#series: "Tutorial"
---


As a first example we use the standard textbook harmonic oscillator in one dimension and fill it with two non-interacting electrons. 

## Input

The first thing to do is to tell {{< octopus >}} what we want it to do. Write the following lines and save the file as {{< file "inp" >}}.

{{< code-block >}}
#include_input doc/tutorials/model_systems/1d_harmonic_oscillator/1.main/inp
{{< /code-block >}}


Most of these input variables should already be familiar. Here is a more detailed explanation for some of the values:
* {{< code-inline >}}{{< variable "Dimensions" >}} = 1{{< /code-inline >}}: This means that the code should run for 1-dimensional systems. Other options are 2, or 3 (the default). You can actually run in 4D too if you have compiled with the configure flag <tt>--max-dim=4</tt>.

* {{< code-inline >}}{{< variable "TheoryLevel" >}} = independent_particles{{< /code-inline >}}: We tell {{< octopus >}} to treat the electrons as non-interacting.

* {{< code-inline >}}{{< variable "Radius" >}} = 10.0{{< /code-inline >}}: The radius of the 1D "sphere," ''i.e.'' a line; therefore domain extends from -10 to +10 bohr.

* {{< code-inline >}}%{{< variable "Species" >}}{{< /code-inline >}}: The species name is "HO", then the potential formula is given, and finally the number of valence electrons. See {{< manual "Input file" "Input file" >}} for a description of what kind of expressions can be given for the potential formula.

* {{< code-inline >}}%{{< variable "Coordinates" >}}{{< /code-inline >}}: add the external potential defined for species "HO" in the position (0,0,0). The coordinates used in the potential formula are relative to this point.

## Output
Now one can execute this file by running {{< octopus >}}. Here are a few things worthy of a closer look in the standard output.

First one finds the listing of the species:

{{< code-block >}}
 ****************************** Species *******************************
 Species "HO" is an user-defined potential.
    Potential = 0.5*x^2
 Number of orbitals:      5
 **********************************************************************
{{< /code-block >}}

The potential is $V(x) = 0.5 x^2$, and 5 Hermite-polynomial orbitals are available for LCAO (the number is based on the valence). The theory level is as we requested:

{{< code-block >}}
 **************************** Theory Level ****************************
 Input: [{{< variable "TheoryLevel" >}} = independent_particles]
 **********************************************************************
{{< /code-block >}}

The electrons are treated as "non-interacting", which means that the Hartree and exchange-correlation terms are not included. This is usually appropriate for model systems, in particular because the standard XC approximations we use for actual electrons are not correct for "effective electrons" with a different mass.

{{< code-block >}}
 Input: [{{< variable "MixField" >}} = none] (what to mix during SCF cycles)
{{< /code-block >}}

Since we are using independent particles (and only one electron) there is no need to use a mixing scheme to accelerate the SCF convergence. 

{{< code-block >}}
 Input: [{{< variable "LCAOStart" >}} = lcao_none]
 Info: Unnormalized total charge =      2.000000
 Info: Renormalized total charge =      2.000000
 Info: Setting up Hamiltonian.
 Orthogonalizing wavefunctions.
{{< /code-block >}}

Starting from scratch means that {{< octopus >}} generates a starting density from the sum of atomic densities. This is then renormalized to integrate to the total number of electrons present in the system. For atomic systems, the default is to find a first estimate for the wave-functions using a {{< manual "Calculations/Ground_State#LCAO" "linear combination of atomic orbitals" >}} (LCAO) technique, using the atomic wavefunctions from the pseudopotentials. However, we do not necessarily have such corresponding wavefunctions for a user-defined potential, so LCAO is turned off by default here. The starting orbitals will then be random but orthogonal.

Now the self-consistent cycle starts and you should see the following output: 

{{< code-block >}}
Info: Starting SCF iteration.
ETA: .......1......2.......3......4......5.......6......7.......8......9......0


*********************** SCF CYCLE ITER #    1 ************************
 etot  =  1.24484521E+00 abs_ev   =  9.49E+01 rel_ev   =  7.62E+01
 ediff =        1.24E+00 abs_dens =  3.27E+00 rel_dens =  1.63E+00
Matrix vector products:     27
Converged eigenvectors:      0

#  State  Eigenvalue [H]  Occupation    Error
      1        0.622423    2.000000   (6.3E+00)

Elapsed time for SCF step     1:          0.00
**********************************************************************

ETA: .......1......2.......3......4......5.......6......7.......8......9......0
{{< /code-block >}}

and after a few iterations it converges:

{{< code-block >}}
*********************** SCF CYCLE ITER #   22 ************************
 etot  =  1.00000001E+00 abs_ev   =  1.95E-08 rel_ev   =  1.95E-08
 ediff =       -1.95E-08 abs_dens =  1.61E-05 rel_dens =  8.07E-06
Matrix vector products:     27
Converged eigenvectors:      0

#  State  Eigenvalue [H]  Occupation    Error
      1        0.500000    2.000000   (1.4E-03)

Elapsed time for SCF step    22:          0.00
**********************************************************************


             Info: Writing states. 2018/08/08 at 14:13:39


        Info: Finished writing states. 2018/08/08 at 14:13:39

Info: SCF converged in   22 iterations
{{< /code-block >}}

Note that, as the electrons are non-interacting, there is actually no self-consistency needed. Nevertheless, {{< octopus >}} takes several SCF iterations to converge. This is because it takes more than one SCF iteration for the eigensolver to converge the wave-functions and eigenvalues. To improve this we can try having a better initial guess for the wave-functions by turning the LCAO on. You can do this by setting <tt>{{< variable "LCAOStart" >}} = lcao_states</tt> -- compare how many iterations and matrix-vector products (in total) are required now. Why? (Hint: what are Hermite polynomials?)

##  Exercises  

The first thing to do is to look at the {{< file "static/info" >}} file. Look at the eigenvalue and total energy. Compare to your expectation from the analytic solution to this problem!

We can also play a little bit with the input file and add some other features. For example, what about the other eigenstates of this system? To obtain them we can calculate some unoccupied states. This can be done by changing the <tt>CalculationMode</tt> to

{{< code-block >}}
 {{< variable "CalculationMode" >}} = unocc
{{< /code-block >}}

in the input file. In this mode {{< octopus >}} will not only give you the occupied states (which contribute to the density) but also the unoccupied ones. Set the number of states with an extra line

{{< code-block >}}
 {{< variable "ExtraStates" >}} = 10
{{< /code-block >}}

that will calculate 10 empty states. A thing to note here is that {{< octopus >}} will need the density for this calculation. (Actually for non-interacting electrons the density is not needed for anything, but since {{< octopus >}} is designed for interacting electrons it will try to read the density anyways.) So if you have already performed a static calculation (<tt>CalculationMode = gs</tt>) it will just use this result. 

Compared to the ground-state calculation we also have to change the convergence criterion. The unoccupied states do not contribute to the density so they might not (and actually will not) converge properly if we use the density for the convergence check. Therefore, {{< octopus >}} now checks whether all calculated states converge separately by just looking at the biggest error in the wave-functions. From the calculation you get another file in the {{< file "static" >}} directory called {{< file "eigenvalues" >}}. The {{< file "info" >}} file will only contain the information about the ground-state; all eigenvalues and occupation numbers will be in the {{< file "eigenvalues" >}} file.

If we also want to plot, say, the wave-function, at the end of the calculation, we have to tell {{< octopus >}} to give us this wave-function and how it should do this. We just include

{{< code-block >}}
 %{{< variable "Output" >}}
    wfs
 %
 {{< variable "OutputFormat" >}} = axis_x
{{< /code-block >}}

The first line tells {{< octopus >}} to give out the wave-functions and the second line says it should do so along the x-axis. We can also select the wave-functions we would like as output, for example the first and second, the fourth and the sixth. (If you don't specify anything {{< octopus >}} will give them all.)

{{< code-block >}}
 {{< variable "OutputWfsNumber" >}} = "1-2,4,6"
{{< /code-block >}}

{{< octopus >}} will store the wave-functions in the same folder {{< file "static" >}} where the {{< file "info" >}} file is, under a meaningful name. They are stored as pairs of the ''x''-coordinate and the value of the wave-function at that position ''x''. One can easily plot them with gnuplot or a similar program.

It is also possible to extract a couple of other things from {{< octopus >}} like the density or the Kohn-Sham potential.

{{< tutorial-footer >}}
---
title: "Particle in an octopus"
tags: ["Beginner", "Ground State", "Model", "User-defined Species", "Independent Particles", "Visualization"]
weight: 4
tutorials: "Model Systems"
difficulties: "beginner"
theories: "Independent particles"
calculation_modes: "Ground state"
species_types: "User-defined species"
features: "Visualization"
description: "2D systems where the shape of the simulation box is defined by an image file."
#series: "Tutorial"
---


{{< octopus >}} can actually run 2D systems where the shape of the simulation box is defined by what is white in an image file. Here is an example of a "particle in an octopus", in which we have a constant potential and an octopus-shaped quantum dot. To run it, you will need to have built the code with the [optional library GDLIB](https://libgd.github.io).

## Input

For this example we will need two files:

#### {{< file "inp" >}}

{{< code-block >}}
#include_input doc/tutorials/model_systems/particle_in_a_octopus/inp 
{{< /code-block >}}


#### {{< file "Gdlib.png" >}} 

Make this file available in the run directory. You can download it by clicking on the image bellow. It is also available in the {{< file "PREFIX/share/octopus" >}} directory from your {{< octopus >}} installation.

{{< figure src="/images/Gdlib.png" >}}

## Plotting

The wavefunction are obtained by using:

{{< code-block >}}
 %{{< variable "Output" >}}
   wfs
 %
 {{< variable "OutputFormat" >}} = plane_z
{{< /code-block >}}

View it in {{< code "gnuplot" >}} with

```bash
 plot 'static/wf-st0001.z=0' u 1:2:3 linetype palette
```

or

```bash
 splot 'static/wf-st0001.z=0' u 1:2:(0):($3*500) with pm3d
```

Where does the wavefunction localize, and why?

## Exercises
* See how the total energy scales with the size of the system (controlled by the {{< code ff >}} parameter in the input file). How does it compare to the formula for a particle in a box?
* Look at the wavefunctions of the unoccupied states.
* Think of a serious application that would use the ability to define the simulation box by an image!

{{< tutorial-footer >}}
---
title: "Use of symmetries in optical spectra from time-propagation"
#tags: ["Advanced", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
tutorials: "Optical Response"
difficulties: "advanced"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
weight: 6
#series: "Tutorial"
description: "Reduce the number of time-propagations required to compute the absorption cross-section."
---


In this tutorial we will see how the spatial symmetries of a molecule can be used to reduce the number of time-propagations required to compute the absorption cross-section.

## Introduction

The dynamic polarizability (related to optical absorption cross-section via $\sigma = \frac{4 \pi \omega}{c} \mathrm{Im} \alpha $) is, in its most general form, a 3x3 tensor. The reason is that we can shine light on the system polarized in any of the three Cartesian axes, and for each of these three cases measure how the dipole of the molecule oscillates along the three Cartesian axes. This usually means that to obtain the full dynamic polarizability of the molecule we usually need to apply 3 different perturbations along $x, y, z\,$, by setting {{< variable "TDPolarizationDirection" >}} to 1, 2, or 3.

However, if the molecule has some symmetries, it is in general possible to reduce the total number of calculations required to obtain the tensor from 3 to 2, or even 1.[^footnote-1]

To use this formalism in {{< octopus >}}, you need to supply some extra information. The most important thing that the code requires is the information about equivalent axes, that is, directions that are related through some symmetry transformation. Using these axes, we construct a reference frame and specify it with the {{< variable "TDPolarization" >}} block. Note that these axes need not be orthogonal, but they must be linearly-independent. The {{< variable "TDPolarizationEquivAxes" >}} tells the code how many equivalent axes there are. Ideally, the reference frame should be chosen to maximize the number of equivalent axes. When using three equivalent axes, an extra input variable, {{< variable "TDPolarizationWprime" >}}, is also required.

Let us give a couple of examples, which should make all these concepts easier to understand.

## Methane

As seen in previous tutorials, the methane molecule has $T_d$ symmetry. This means that it is trivial to find three linearly-independent equivalent axes such that only one time-propagation is needed to obtain the whole tensor. As it happens, we can use the usual $x$, $y$, and $z$ directions, with all of them being equivalent (check the [symmetry operations](https://en.wikipedia.org/wiki/Symmetry_operation) of the $T_d$ [point group](https://en.wikipedia.org/wiki/Molecular_symmetry) if you are not convinced). Therefore we can perform just one propagation with the perturbation along the $x$ direction adding the following to the input file used previously: 

{{< code-block >}}
#include_input_snippet doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/2.td/inp polarization
{{< /code-block >}}

{{% expand "complete input file" %}}
{{< code-block >}}
#include_input doc/tutorials/optical_response/use_of_symmetries_in_optical_spectra_from_time-propagation/2.td/inp
{{< /code-block >}}
{{% /expand %}}


Note that we had omitted the blocks {{< variable "TDPolarization" >}} and {{< variable "TDPolarizationWprime" >}} in the previous tutorials, as these are their default 

Once the time-propagation is finished, you will find, as usual, a {{< file "td.general/multipoles" >}} file. This time the file contains all the necessary information for {{< octopus >}} to compute the full tensor, so running the {{< file "oct-propagation_spectrum" >}} utility will produce two files: {{< file "cross_section_vector.1" >}} and {{< file "cross_section_tensor" >}}. The later should look like this (assuming you used the converged grid parameters from the
[Convergence of the optical spectra](../Convergence of the optical spectra) tutorial):

{{< code-block >}}
# nspin         1
# kick mode    0
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction    1
# Equiv. axes  3
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
#       Energy         (1/3)*Tr[sigma]    Anisotropy[sigma]      sigma(1,1,1)        sigma(1,2,1)        sigma(1,3,1)       ...
#        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]           ...
      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00    ...
      0.10000000E-01     -0.21016843E-08      0.17355148E-11     -0.21016843E-08      0.12207713E-11      0.88665555E-13    ...
      0.20000000E-01     -0.83911822E-08      0.69337823E-11     -0.83911822E-08      0.48772812E-11      0.35411614E-12    ...
...
{{< /code-block >}}

Try comparing the spectrum for each component of the $\sigma$ tensor.

## Linear molecule

Now let us look at a linear molecule. In this case, you might think that we need two calculations to obtain the whole tensor, one for the direction along the axis of the molecule, and another for the axis perpendicular to the molecule. The fact is that we need only one, in a specially chosen direction, so that our field has components both along the axis of the molecule and perpendicular to it. Let us assume that the axis of the molecule is oriented along the $x\,$ axis. Then we can use

{{< code-block >}}
 %{{< variable "TDPolarization" >}}
  1/sqrt(2) | -1/sqrt(2) | 0
  1/sqrt(2) |  1/sqrt(2) | 0
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
 {{< variable "TDPolarizationDirection" >}} = 1
 {{< variable "TDPolarizationEquivAxes" >}} = 3
 %{{< variable "TDPolarizationWprime" >}}
  1/sqrt(2) |  0         | 1/sqrt(2)
 %
{{< /code-block >}}

You should try to convince yourself that the three axes are indeed equivalent and linearly independent. The first and second axes are connected through a simple reflection in the $xz$ plane, transforming the $y$ coordinate from $-1/\sqrt{2}$ into $1/\sqrt{2}$. {{< variable "TDPolarizationWprime" >}} should be set to the result obtained by applying the inverse symmetry operation to the third axis. This actually leaves the third axis unchanged.

## Planar molecule

Finally, let us look at a general planar molecule (in the $xy$ plane). In principle we need only two calculations (that is reduced to one if more symmetries are present like, ''e.g.'', in benzene). In this case we chose one of the polarization axes on the plane and the other two rotated 45 degrees:

{{< code-block >}}
 %{{< variable "TDPolarization" >}}
  1/sqrt(2) | 0 | 1/sqrt(2)
  1/sqrt(2) | 0 |-1/sqrt(2)
  0         | 1 | 0
 %
 
 {{< variable "TDPolarizationEquivAxes" >}} = 2
{{< /code-block >}}

In this case, we need two runs, one for {{< variable "TDPolarizationDirection" >}} equal to 1, and another for it equal to 3. Note that if there are less than 3 equivalent axes, {{< variable "TDPolarizationWprime" >}} is irrelevant.

[^footnote-1]: {{< article title="On the use of Neumann's principle for the calculation of the polarizability tensor of nanostructures" authors="M.J.T. Oliveira, A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Nanoscience and Nanotechnology" volume="8" pages="1-7" year="2008" doi="10.1166/jnn.2008.142" arxiv="0710.2624v1" >}}


{{< tutorial-footer >}}
---
title: "Convergence of the optical spectra"
Weight: 2
#tags: ["Beginner", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
difficulties: "beginner"
description: "Convergence study of the spectrum with respect to the grid parameters."
---


In the previous {{< tutorial "Response/Optical Spectra from time-propagation" "Optical Spectra from TD" >}} we have seen how to obtain the absorption spectrum of a molecule. Here we will perform a convergence study of the spectrum with respect to the grid parameters using again the methane molecule as an example.

## Convergence with the spacing

In the {{< tutorial "Basics/Total energy convergence" "Total energy convergence tutorial" >}} we found out that the total energy of the methane molecule was converged to within 0.1 eV for a spacing of 0.18 Å and a radius of 3.5 Å. We will therefore use these values as a starting point for our convergence study. Like for the total energy, the only way to check if the absorption spectrum is converged with respect to the spacing is by repeating a series of calculations, with identical input files except for the grid spacing. The main difference in this case is that each calculation for a given spacing requires several steps:

- Ground-state calculation,
- Three time-dependent propagations with perturbations along ''x'', ''y'' and ''z'', 
- Calculation of the spectrum using the {{< file "oct-propagation_spectrum" >}} utility.

As we have seen before, the $T_d$ symmetry of methane means that the spectrum is identical for all directions, so in this case we only need to perform the convergence study with a perturbation along the ''x'' direction.

All these calculations can be done by hand, but we would rather use a script for this boring an repetitive task. We will then need the following two files.

#### inp
This file is the same as the one used in the previous {{< tutorial "Response/Optical Spectra from time-propagation" "Optical Spectra from TD" >}}:

{{< code-block >}}
#include_input doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/inp
{{< /code-block >}}


#### spacing.sh

And this is the bash script:

```bash
#include_file doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/spacing.sh
```

This script performs the three steps mentioned above for a given list of spacings. You will notice that we are not including smaller spacings than 0.18 Å. This is because we know from previous experience that the optical spectrum is well converged for spacings larger than the one required to converge the total energy.

Now run the script typing {{< command-line "source spacing.sh" >}}. This will take some time, so you might want to have a break. Once the script is finished running, we need to compare the spectra for the different spacings. Using your favorite plotting software, e.g. {{< command-line "gnuplot" >}}, plot the second column of the {{< file "cross_section_vector" >}} files vs the first column. You can see on the right how this plot should look like. To make the plot easier to read, we have restricted the energy range and omitted some of the spacings. From the plot it is clear that a spacing of 0.24 Å is sufficient to have the position of the peaks converged to within 0.1 eV.

#include_eps doc/tutorials/optical_response/convergence_of_the_optical_spectra/1.spacing/Absorption_spectrum_CH4_spacing.eps caption="Convergence with spacing of methane absorption spectrum."



## Convergence with the radius

We will now study how the spectrum changes with the radius. We will proceed as for the spacing, by performing a series of calculations with identical input files except for the radius. For this we will use the two following files. 

#### inp
This file is the same as the prevous one, with the exception of the spacing and the time-step. As we have just seen, a spacing of 0.24 Å is sufficient, which in turns allows us to use a larger time-step.

{{< code-block >}}
#include_input doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/inp
{{< /code-block >}}


#### radius.sh

And this is the bash script:

```bash
#include_file doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/radius.sh
```

This script works just like the previous one, but changing the radius instead of the spacing. Now run the script ( {{< command "source radius.sh" >}}) and wait until it is finished.

Now plot the spectra from the {{< file "cross_section_vector" >}} files. You should get something similar to the plot on the right. We see that in this case the changes are quite dramatic. To get the first peak close to 10 eV converged within 0.1 eV a radius of 6.5 Å is necessary. The converged transition energy of 9.2 eV agrees quite well with the experimental value of 9.6 eV and with the TDDFT value of 9.25 eV obtained by other authors.[^footnote-1]

#include_eps doc/tutorials/optical_response/convergence_of_the_optical_spectra/2.radius/Absorption_spectrum_CH4_radius.eps caption="Convergence with spacing of methane absorption spectrum."


What about the peaks at higher energies? Since we are running in a box with zero boundary conditions, we will also have the box states in the spectrum. The energy of those states will change with the inverse of the size of the box and the corresponding peaks will keep shifting to lower energies. However, as one is usually only interested in bound-bound transitions in the optical or near-ultraviolet regimes, we will stop our convergence study here. As an exercise you might want to converge the next experimentally observed transition. Just keep in mind that calculations with larger boxes might take quite some time to run.

[^footnote-1]: {{< article title="Time-Dependent Density Functional Theory Calculations of Photoabsorption Spectra in the Vacuum Ultraviolet Region" authors="N. N. Matsuzawa, A. Ishitani, D. A. Dixon, and T. Uda" journal="J. Phys. Chem. A" volume="105" pages="4953–4962" year="2001" doi="10.1021/jp003937v" >}}


{{< tutorial-footer >}}
---
title: "Optical spectra from time-propagation"
Weight: 1
#tags: ["Beginner", "Time-dependent", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-propagation_spectrum"]
difficulties: "beginner"
theories: "DFT"
calculation_modes: "Time-dependent"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical Absorption"
utilities: "oct-propagation_spectrum"
tutorials: "Optical Response"
#series: "Tutorial"
description: "Absorption spectrum of a molecule from the explicit solution of the time-dependent Kohn-Sham equations: methane"
---


In this tutorial we will learn how to obtain the absorption spectrum of a molecule from the explicit solution of the time-dependent Kohn-Sham equations. We choose as a test case methane (CH<sub>4</sub>).

## Ground state

Before starting out time-dependent simulations, we need to obtain the ground state of the system. For this we use basically the same {{< file "inp" >}} file as in the 
{{< tutorial "Basics/Total energy convergence" "total energy convergence tutorial" >}}:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/1.gs/inp
{{< /code-block >}}

After running {{< octopus >}}, we will have the Kohn-Sham wave-functions of the ground-state in the directory {{< file "restart/gs" >}}. As we are going to propagate these wave-functions, they have to be well converged. It is important not only to converge the energy (that is relatively easy to converge), but the density. 
The default {{< code-inline >}}{{< variable "ConvRelDens" >}} = 1e-05{{< /code-inline >}} is usually enough, though.

## Time-dependent run

To calculate absorption, we excite the system with an infinitesimal electric-field pulse, and then propagate the time-dependent Kohn-Sham equations for a certain time ''T''. The spectrum can then be evaluated from the time-dependent dipole moment.

#### Input

This is how the input file should look for the time propagation. It is similar to the one from the {{< tutorial "Basics/Time-dependent propagation" "Time-dependent propagation tutorial" >}}.

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_time-propagation/2.td_x/inp
{{< /code-block >}}

Besides changing the {{< variable "CalculationMode" >}} to {{< code td >}}, we have added {{< code-inline >}}{{< variable "FromScratch" >}} = yes{{< /code-inline >}}. This will be useful if you decide to run the propagation for other polarization directions (see bellow). For the time-evolution we use again the Approximately Enforced Time-Reversal Symmetry (aetrs) propagator. The time-step is chosen such that the propagation remains numerically stable. You should have learned how to set it up in the tutorial {{< tutorial "Baics/Time-dependent propagation" "time-dependent propagation tutorial" >}}. Finally, we set the number of time steps with the variable {{< variable "TDMaxSteps" >}}. To have a maximum propagation time of 10 $\hbar/{\rm eV}$ we will need around 4350 iterations.

We have also introduced two new input variables to define our perturbation:

* {{< variable "TDDeltaStrength" >}}: this is the strength of the perturbation. This number should be small to keep the response linear, but should be sufficiently large to avoid numerical problems.

* {{< variable "TDPolarizationDirection" >}}: this variable sets the polarization of our perturbation to be on the first axis (''x'').

Note that you will be calculating the singlet dipole spectrum. You can also obtain the triplet by using {{< variable "TDDeltaStrengthMode" >}}, and other multipole responses by using {{< variable "TDKickFunction" >}}. For details on triplet calculations see the {{< tutorial "Response/Triplet excitations" "Triplet excitations tutorial" >}}.

You can now start {{< octopus >}} and go for a quick coffee (this should take a few minutes depending on your machine). Propagations are slow, but the good news is that they scale very well with the size of the system. This means that even if methane is very slow, a molecule with 200 atoms can still be calculated without big problems.

#### Output

The output should be very similar to the one from the {{< tutorial "Basics/Time-dependent propagation" "Time-dependent propagation tutorial" >}}. The main difference is the information about the perturbation:

{{< code-block >}}
Info: Applying delta kick: k =    0.005292
Info: kick function: dipole.
Info: Delta kick mode: Density mode
{{< /code-block >}}

You should also get a {{< file "td.general/multipoles" >}} file, which contains the necessary information to calculate the spectrum. The beginning of this file should look like this:

{{< code-block >}}
################################################################################
# HEADER
# nspin         1
# lmax          1
# kick mode     0
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction     1
# Equiv. axes   0
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
# Iter             t          Electronic charge(1)       <x>(1)              <y>(1)              <z>(1)       
#[Iter n.]     [hbar/eV]           Electrons              [A]                 [A]                 [A]         
################################################################################
       0  0.000000000000e+00  8.000000000000e+00 -5.105069841261e-11 -2.458728679159e-11 -6.191777548908e-11
       1  2.300000000000e-03  7.999999999898e+00 -1.441204189916e-03 -2.454483489431e-11 -6.178536414960e-11
       2  4.600000000000e-03  7.999999999794e+00 -2.849596623682e-03 -2.441311503813e-11 -6.139129111463e-11
       3  6.900000000000e-03  7.999999999692e+00 -4.212595751069e-03 -2.420204397245e-11 -6.075589410104e-11
{{< /code-block >}}

Note how the dipole along the ''x'' direction (forth column) changes in response to the perturbation. 

##  Optical spectra  

In order to obtain the spectrum for a general system one would need to perform three time-dependent runs, each for a perturbation along a different Cartesian direction (''x'', ''y'', and ''z''). In practice this is what you would need to do:
- Set the direction of the perturbation along ''x'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 1{{< /code-inline >}}),
- Run the time-propagation,
- Rename the {{< file "td.general/multipoles" >}} file to {{< file "td.general/multipoles.1" >}},
- Repeat the above step for directions ''y'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 2{{< /code-inline >}}) and ''z'' ({{< code-inline >}}{{< variable "TDPolarizationDirection" >}} = 3{{< /code-inline >}}) to obtain files {{< file "td.general/multipoles.2" >}} and {{< file "td.general/multipoles.3" >}}.

Nevertheless, the $T_d$ symmetry of methane means that the response is identical for all directions and the absorption spectrum for ''x''-polarization will in fact be equivalent to the spectrum averaged over the three directions. You can perform the calculations for the ''y'' and ''z'' directions if you need to convince yourself that they are indeed equivalent, but you this will not be necessary to complete this tutorial.

Note that {{< octopus >}} can actually use the knowledge of the symmetries of the system when calculating the spectrum. However, this is fairly complicated to understand for the purposes of this introductory tutorial, so it will be covered in the 
{{< tutorial "Response/Use of symmetries in optical spectra from time-propagation" "Use of symmetries in optical spectra from time-propagation tutorial" >}}.

### The {{< file "oct-propagation_spectrum" >}} utility

{{< octopus >}} provides an utility called {{< manual "external_utilities/oct-propagation_spectrum" "oct-propagation_spectrum" >}} to process the {{< file "multipoles" >}} files and obtain the spectrum. 

#### Input
This utility requires little information from the input file, as most of what it needs is provided in the header of the {{< file "multipoles" >}} files. If you want you can reuse the same input file as for the time-propagation run, but the following input file will also work:

{{< code-block >}}
 {{< variable "UnitsOutput" >}} = eV_angstrom
{{< /code-block >}}

Now run the utility. Don't worry about the warnings generated, we know what we are doing!

#### Output

This is what you should get:

{{< code-block >}}
************************** Spectrum Options **************************
Input: [PropagationSpectrumType = AbsorptionSpectrum]
Input: [SpectrumMethod = fourier]
Input: [PropagationSpectrumDampMode = polynomial]
Input: [PropagationSpectrumTransform = sine]
Input: [PropagationSpectrumStartTime = 0.000 hbar/eV]
Input: [PropagationSpectrumEndTime = -.3675E-01 hbar/eV]
Input: [PropagationSpectrumEnergyStep = 0.1000E-01 eV]
Input: [PropagationSpectrumMaxEnergy = 20.00 eV]
Input: [PropagationSpectrumDampFactor = -27.21 hbar/eV^-1]
Input: [PropagationSpectrumSigmaDiagonalization = no]
**********************************************************************

File "multipoles" found. This will be the only file to be processed.
(If more than one file is to be used, the files should be called
"multipoles.1", "multipoles.2", etc.)


** Warning:
**   The file "multipoles" tells me that the system has no usable symmetry.
**   However, I am only using this file; cannot calculate the full tensor.
**   A file "XXXX_vector" will be generated instead.


Octopus emitted 1 warning.
{{< /code-block >}}

You will notice that the file {{< file "cross_section_vector" >}} is created. If you have all the required information (either a symmetric molecule and one multipole file, or a less symmetric molecule and multiple multipole files), {{< file "oct-propagation_spectrum" >}} will generate the whole polarizability tensor, {{< file "cross_section_tensor" >}}. If not enough multipole files are available, it can only generate the response to the particular polarization you chose in you time-dependent run, {{< file "cross_section_vector" >}}. 

### Cross-section vector

Let us look first at {{< file "cross_section_vector" >}}:

{{< code-block >}}
# nspin         1
# kick mode    0
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction    1
# Equiv. axes  0
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
#%
{{< /code-block >}}

The beginning of the file just repeats some information concerning the run.

{{< code-block >}}
# Number of time steps =     4350
# PropagationSpectrumDampMode   =    2
# PropagationSpectrumDampFactor =   -27.2114
# PropagationSpectrumStartTime  =     0.0000
# PropagationSpectrumEndTime    =    10.0050
# PropagationSpectrumMaxEnergy  =    20.0000
# PropagationSpectrumEnergyStep =     0.0100
#%
{{< /code-block >}}

Now comes a summary of the variables used to calculate the spectra. Note that you can change all these settings just by adding these variables to the input file before running {{< file "oct-propagation_spectrum" >}}.  Of special importance are perhaps {{< variable "PropagationSpectrumMaxEnergy" >}} that sets the maximum energy that will be calculated, and {{< variable "PropagationSpectrumEnergyStep" >}} that determines how many points your spectrum will contain. To have smoother curves, you should reduce this last variable.

{{< code-block >}}
# Electronic sum rule       =         3.682690
# Static polarizability (from sum rule) =         2.064454 A
#%
{{< /code-block >}}

Now comes some information from the sum rules. The first is just the ''f''-sum rule, which should yield the number of active electrons in the calculations. We have 8 valence electrons (4 from carbon and 4 from hydrogen), but the sum rule gives a number that is much smaller than 8! The reason is that we are just summing our spectrum up to 20 eV (see {{< variable "PropagationSpectrumMaxEnergy" >}}), but for the sum rule to be fulfilled, we should go to infinity. Of course infinity is a bit too large, but increasing 20 eV to a somewhat larger number will improve dramatically the ''f''-sum rule. The second number is the static polarizability also calculated from the sum rule. Again, do not forget to converge this number with {{< variable "PropagationSpectrumMaxEnergy" >}}.

{{< code-block >}}
#       Energy        sigma(1, nspin=1)   sigma(2, nspin=1)   sigma(3, nspin=1)   StrengthFunction(1)
#        [eV]               [A^2]               [A^2]               [A^2]               [1/eV]
      0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00
      0.10000000E-01     -0.42061649E-08     -0.71263367E-14      0.21595354E-13     -0.38321132E-08
      0.20000000E-01     -0.16805683E-07     -0.28462613E-13      0.86270383E-13     -0.15311164E-07
...
{{< /code-block >}}

Finally comes the spectrum. The first column is the energy (frequency), the next three columns are a row of the cross-section tensor, and the last one is the strength function for this run.

The dynamic polarizability is related to optical absorption cross-section via $\sigma \left( \omega \right) = \frac{4 \pi \omega}{c} \mathrm{Im}\ \alpha \left( \omega \right) $ in atomic units, or more generally $4 \pi \omega \tilde{\alpha}\ \mathrm{Im}\ \alpha \left( \omega \right) $ (where $\tilde{\alpha}$ is the fine-structure constant) or $\frac{\omega e^2}{\epsilon_0 c} \mathrm{Im}\ \alpha \left( \omega \right) $. The cross-section is related to the strength function by $S \left( \omega \right) = \frac{mc}{2 \pi^2 \hbar^2} \sigma \left( \omega \right)$.

### Cross-section tensor

If all three directions were done, we would have four files: {{< file "cross_section_vector.1" >}}, {{< file "cross_section_vector.2" >}}, {{< file "cross_section_vector.3" >}}, and {{< file "cross_section_tensor" >}}. The latter would be similar to:

{{< code-block >}}
# nspin         1
# kick mode    0
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction    3
# Equiv. axes  0
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
#       Energy         (1/3)*Tr[sigma]    Anisotropy[sigma]      sigma(1,1,1)        sigma(1,2,1)        sigma(1,3,1)        ...
#        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]            ...
      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00      0.00000000E+00     ...
      0.10000000E-01     -0.42062861E-08      0.16626408E-12     -0.42061649E-08      0.18742668E-12      0.18909806E-12     ...
      0.20000000E-01     -0.16806167E-07      0.66452380E-12     -0.16805683E-07      0.74882279E-12      0.75548381E-12     ...
...
{{< /code-block >}}

#include_eps doc/tutorials/optical_response/optical_spectra_from_time-propagation/6.propagation_spectra_tensor/Absorption_spectrum_CH4.eps caption="Absorption spectrum of methane"

The columns are now the energy, the average absorption coefficient (Tr $\sigma/3$), the anisotropy, and 
then all the 9 components of the tensor. The third number is the spin component, just 1 here since it is unpolarized. The anisotropy is defined as

$$
 (\\Delta \\sigma)^2 = \\frac{1}{3} \\{ 3{\\rm Tr}(\\sigma^2) - \[{\\rm Tr}(\\sigma)\]^2 \\}
$$

The reason for this definition is that it is identically equal to:

$$
 (\\Delta \\sigma)^2 = (1/3) \[ (\\sigma\_1-\\sigma\_2)^2 + (\\sigma\_1-\\sigma\_3)^2 + (\\sigma\_2-\\sigma\_3)^2 \]\\,
$$

where $\{\sigma_1, \sigma_2, \sigma_3\}\,$ are the eigenvalues of $\sigma\,$. An "isotropic" tensor is characterized by having three equal eigenvalues, which leads to zero anisotropy. The more different that the eigenvalues are, the larger the anisotropy is.

If you now plot the absorption spectrum (column 5 vs 1 in {{< file "cross_section_vector" >}}, but you can use {{< file "cross_section_tensor" >}} for this exercise in case you did propagate in all three directions), you should obtain the plot shown on the right. Of course, you should now try to converge this spectrum with respect to the calculation parameters. In particular:

* Increasing the total propagation time will reduce the width of the peaks. In fact, the width of the peaks in this methods is absolutely artificial, and is inversely proportional to the total propagation time. Do not forget, however, that the area under the peaks has a physical meaning: it is the oscillator strength of the transition.
* A convergence study with respect to the spacing and box size might be necessary. This is covered in the next tutorial.

Some questions to think about:
* What is the equation for the peak width in terms of the propagation time? How does the observed width compare to your expectation?
* What is the highest energy excitation obtainable with the calculation parameters here? Which is the key one controlling the maximum energy?
* Why does the spectrum go below zero? Is this physical? What calculation parameters might change this?


{{< tutorial-footer >}}
---
title: "Optical spectra from Casida"
#tags: ["Beginner", "Unoccupied", "Casida", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "oct-casida_spectrum"]
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Time-dependent"
difficulties: "beginner"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Optical absorption"
utilities: "oct-casida_spectrum"
Weight: 3
#series: "Tutorial"
description: "Calculate the absorption spectrum of methane using Casida's equations"
---


In this tutorial we will again calculate the absorption spectrum of methane, but this time using Casida's equations.


## Ground-state

Once again our first step will be the calculation of the ground state. We will use the following input file:

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_casida/1.gs/inp
{{< /code-block >}}

Note that we are using the values for the spacing and radius that were found in the {{< tutorial "Response/Convergence of the optical spectra" "Convergence_of_the_optical_spectra tutorial" >}} to converge the absorption spectrum.

## Unoccupied States

The Casida equation is a (pseudo-)eigenvalue equation written in the basis of particle-hole states. This means that we need both the occupied states -- computed in the ground-state calculation -- as well as the unoccupied states, that we will now obtain, via a non-self-consistent calculation using the density computed in {{< code gs >}}. The input file we will use is

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_casida/2.unocc/inp
{{< /code-block >}}

Here we have changed the {{< variable "CalculationMode" >}} to {{< code unocc>}} and added 10 extra states by setting the {{< variable "ExtraStates" >}} input variable.

By running {{< octopus >}}, you will obtain the first 10 unoccupied states (do not forget to run the ground-state calculation first). The solution of the unoccupied states is controlled by the variables in section {{< variable "Eigensolver" >}}. You can take a look at the eigenvalues of the unoccupied states in the file {{< file "static/eigenvalues" >}}:

{{< code-block >}}
All states converged.
Criterion =      0.100000E-05

Eigenvalues [eV]
 #st  Spin   Eigenvalue      Occupation     Error
   1   --   -16.787587       2.000000      (7.6E-07)
   2   --    -9.464233       2.000000      (7.6E-07)
   3   --    -9.464233       2.000000      (7.6E-07)
   4   --    -9.464233       2.000000      (7.6E-07)
   5   --    -0.287342       0.000000      (9.8E-07)
   6   --     0.822599       0.000000      (8.0E-07)
   7   --     0.822599       0.000000      (8.6E-07)
   8   --     0.822599       0.000000      (6.5E-07)
   9   --     1.599977       0.000000      (7.9E-07)
  10   --     1.599977       0.000000      (9.2E-07)
  11   --     1.599977       0.000000      (7.3E-07)
  12   --     1.813297       0.000000      (7.5E-07)
  13   --     1.813297       0.000000      (7.6E-07)
  14   --     2.309432       0.000000      (8.8E-07)
{{< /code-block >}}

##  Casida calculation  

Now modify the {{< variable "CalculationMode" >}} to {{< code casida>}} and rerun {{< octopus >}}. Note that by default {{< octopus >}} will use all occupied and unoccupied states that it has available. 

Sometimes, it is useful not to use all states. For example, if you have a molecule with 200 atoms and 400 occupied states ranging from -50 to -2 eV, and you are interested in looking at excitations in the visible, you can try to use only the states that are within 10 eV from the Fermi energy. You could select the states to use with {{< variable "CasidaKohnShamStates" >}} or {{< variable "CasidaKSEnergyWindow" >}}.

A new directory will appear named {{< file "casida" >}}, where you can find the file {{< file "casida/casida" >}}:

{{< code-block >}}
                E [eV]         <x> [A]         <y> [A]         <z> [A]             <f>
     1  9.27286861E+00  7.88719340E-02  3.02317792E-01 -1.41580187E-01  9.54564687E-02
     2  9.27291282E+00 -3.06114652E-01  5.96207735E-03 -1.57881218E-01  9.62734206E-02
     3  9.27292772E+00 -1.36862383E-01  1.62030534E-01  2.71506435E-01  9.63001401E-02
     4  1.02415035E+01 -3.18483852E-05  1.77939091E-02 -7.67661722E-04  2.84230878E-04
     5  1.02439430E+01 -6.36128746E-05  8.14736582E-03 -5.21336390E-04  5.97390614E-05
     6  1.02564557E+01 -1.79749506E-03 -1.09950087E-01  4.16624712E-03  1.08663407E-02
     7  1.02593552E+01  3.83016766E-03  3.82893874E-03  1.03764102E-01  9.69062218E-03
     8  1.02594070E+01 -1.04409592E-01  2.07398640E-03  3.78813752E-03  9.80169803E-03
...
{{< /code-block >}}

The \<x\>, \<y\> and \<z\> are the transition dipole moments:

$$
  \<x\> = \<\\Phi\_0|x|\\Phi\_I\>
  \\,;\\qquad
  \<y\> = \<\\Phi\_0|y|\\Phi\_I\>
  \\,;\\qquad
  \<z\> = \<\\Phi\_0|z|\\Phi\_I\>
$$

where $\Phi_0$ is the ground state and $\Phi_I$ is the given excitation. The
oscillator strength is given by:

$$
  f\_I = \\frac{2 m\_e}{3 \\hbar^2} \\omega\_I \\sum\_{n\\in x,y,z} |\<\\Phi\_0|n|\\Phi\_I\>|^2\\,
$$

as the average over the three directions. The optical absorption spectrum can be given as the "strength function",
which is

$$
  S(\\omega) = \\sum\_I f\_I \\delta(\\omega-\\omega\_I)\\,
$$

Note that the excitations are degenerate with degeneracy 3. This could already be expected from the $T_d$ symmetry of methane.

''Note that within the degenerate subspaces, there is some arbitrariness (possibly dependent on the details of your compilation and machine) in the linear combinations of transitions. Therefore, you should not be concerned if you do not have the same results for the components of the transition dipole moments (above) and analysis of the excitations (below). However, the energies and resulting spectra should agree.''

Further information concerning the excitations can be found in the directory {{< file "casida/casida_excitations" >}}. For example, the first excitation at 9.27 eV is analyzed in the file {{< file "casida/casida_excitations/00001" >}}:

{{< code-block >}}
# Energy [eV] =    9.27287E+00
# <x> [A] =    7.88719E-02
# <y> [A] =    3.02318E-01
# <z> [A] =   -1.41580E-01
           1           5           1  -1.5104487570024249E-005
           2           5           1  0.46406175367379687     
           3           5           1  0.86441581738495976     
           4           5           1 -0.18515534208137396     
           1           6           1   3.2556357546167634E-003
...
{{< /code-block >}}

These files contain basically the eigenvector of the Casida equation. The first two columns are respectively the index of the occupied and the index of the unoccupied state, the third is the spin index (always 1 when spin-unpolarized), and the fourth is the coefficient of that state in the Casida eigenvector. This eigenvector is normalized to one, so in this case one can say that 74.7% (0.864<sup>2</sup>) of the excitation is from state 3 to state 5 (one of 3 HOMO orbitals->LUMO) with small contribution from some other transitions.

##  Absorption spectrum  

#include_eps doc/tutorials/optical_response/optical_spectra_from_casida/4.spectrum/Absorption_spectrum_CH4_casida.eps caption="Absorption spectrum of methane calculated with time-propagation and with the Casida equation."

To visualize the spectrum, we need to broaden these delta functions with the utility {{< file "oct-casida_spectrum" >}} (run in your working directory, not in {{< file "casida" >}}). It convolves the delta functions with Lorentzian functions. The operation of {{< file "oct-casida_spectrum" >}} is controlled by the variables {{< variable "CasidaSpectrumBroadening" >}} (the width of this Lorentzian), {{< variable "CasidaSpectrumEnergyStep" >}}, {{< variable "CasidaSpectrumMinEnergy" >}}, and {{< variable "CasidaSpectrumMaxEnergy" >}}. If you run  {{< file "oct-casida_spectrum" >}} you obtain the file {{< file "casida/spectrum.casida" >}}. It contains all columns of {{< file "casida" >}} broadened. If you are interested in the total absorption spectrum, then you should plot the first and fifth columns. You should obtain a picture like the one on the right.

Comparing the spectrum obtained with the time-propagation in the {{< tutorial "Response/Convergence of the optical spectra" "Convergence of the optical spectra tutorial" >}} with the one obtained with the Casida approach using the same grid parameters, we can see that

* The peaks of the time-propagation are broader. This can be solved by either increasing the total propagation time, or by increasing the broadening in the Casida approach.
* The first two peaks are nearly the same. Probably also the third is OK, but the low resolution of the time-propagation does not allow to distinguish the two close peaks that compose it.
* For high energies the spectra differ a lot. The reason is that we only used 10 empty states in the Casida approach. In order to describe better this region of the spectrum we would need more. This is why one should always check the convergence of relevant peaks with respect to the number of empty states.

You probably noticed that the Casida calculation took much less time than the time-propagation. This is clearly true for small or medium-sized systems. However, the implementation of Casida in {{< octopus >}} has a much worse scaling with the size of the system than the time-propagation, so for larger systems the situation may be different. Note also that in Casida one needs a fair amount of unoccupied states which are fairly difficult to obtain.

If you are interested, you may also compare the Casida results against the Kohn-Sham eigenvalue differences in {{< file "casida/spectrum.eps_diff" >}} and the Petersilka approximation to Casida in {{< file "casida/spectrum.petersilka" >}}.

{{< tutorial-footer >}}
---
title: "Optical Response"
weight: 20
description: "Calculating the optical response"
---

In this series of linked tutorials we will show how to use TDDFT to calculate the optical absorption of molecules using different methods implemented in Octopus. These tutorials assume that you are already familiar with using Octopus. If that is not the case, start first with the Octopus basics tutorials. 

{{< tutorial-series "Optical Response" >}}

---
title: "Triplet excitations"
#tags: ["Advanced", "Time-dependent", "Casida", "Molecule", "Pseudopotentials", "DFT", "Triplet Excitations", "oct-propagation_spectrum", "oct-casida_spectrum"]
weight: 5
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: ["Time-dependent", "Casida"]
system_types: "Molecule"
species_types: "Pseudopotentials"
features: "Triplet Excitations"
utilities: ["oct-propagation_spectrum", "oct-casida_spectrum"]
difficulties: "advanced"
#series: "Tutorial"
description: "Calculate triplet excitations for methane with time-propagation and Casida methods."
---


In this tutorial, we will calculate triplet excitations for methane with time-propagation and Casida methods.

## Time-propagation

### Ground-state
We begin with a spin-polarized calculation of the ground-state, as before but with {{< code-inline >}}{{< variable "SpinComponents" >}} = spin_polarized{{< /code-inline >}} specified now.

{{< code-block >}}
#include_input doc/tutorials/optical_response/triplet_excitations/1.time_propagation/1.gs_pol/inp
{{< /code-block >}}

You can verify that the results are identical in detail to the non-spin-polarized calculation since this is a non-magnetic system.

### Time-propagation

Next, we perform the time-propagation using the following input file:

{{< code-block >}}
#include_input doc/tutorials/optical_response/triplet_excitations/1.time_propagation/2.td/inp
{{< /code-block >}}

Besides the {{< variable "SpinComponents" >}} variable, the main difference is the type of perturbation that is applied to the system. By setting {{< code-inline >}}
{{< variable "TDDeltaStrengthMode" >}} = kick_spin{{< /code-inline >}}, the kick will have opposite sign for up and down states. Whereas the ordinary kick ({{< code kick_density >}}) yields the response to a homogeneous electric field, ''i.e.'' the electric dipole response, this kick yields the response to a homogeneous magnetic field, ''i.e.'' the magnetic dipole response. Note however that only the spin degree of freedom is coupling to the field; a different calculation would be required to obtain the orbital part of the response. Only singlet excited states contribute to the spectrum with {{< code kick_density >}}, and only triplet excited states contribute with {{< code kick_spin >}}. We will see below how to use symmetry to obtain both at once with {{< code kick_spin_and_density >}}.

### Spectrum

When the propagation completes, run the {{< manual "external_utilities/oct-propagation_spectrum" "oct-propagation_spectrum" >}} utility to obtain the spectrum. This is how the {{< file "cross_section_vector" >}} should look like.

{{< code-block >}}
# nspin         2
# kick mode    1
# kick strength    0.005291772086
# pol(1)           1.000000000000    0.000000000000    0.000000000000
# pol(2)           0.000000000000    1.000000000000    0.000000000000
# pol(3)           0.000000000000    0.000000000000    1.000000000000
# direction    1
# Equiv. axes  0
# wprime           0.000000000000    0.000000000000    1.000000000000
# kick time        0.000000000000
#%
# Number of time steps =     2500
# PropagationSpectrumDampMode   =    2
# PropagationSpectrumDampFactor =     4.0817
# PropagationSpectrumStartTime  =     0.0000
# PropagationSpectrumEndTime    =    10.0000
# PropagationSpectrumMaxEnergy  =    20.0000
# PropagationSpectrumEnergyStep =     0.0100
#%
#       Energy        sigma(1, nspin=1)   sigma(2, nspin=1)   sigma(3, nspin=1)   sigma(1, nspin=2)   sigma(2, nspin=2)   sigma(3, nspin=2)   StrengthFunction(1) StrengthFunction(2)
#        [eV]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]               [A^2]               [1/eV]              [1/eV]       
      0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00     -0.00000000E+00
      0.10000000E-01     -0.24458945E-08      0.16678045E-11      0.19693772E-11      0.53059881E-08      0.20990576E-11      0.66109455E-12     -0.22283826E-08      0.48341299E-08
      0.20000000E-01     -0.97344740E-08      0.66632996E-11      0.78681086E-11      0.21157433E-07      0.83862450E-11      0.26412279E-11     -0.88687932E-08      0.19275916E-07
{{< /code-block >}}

#include_eps doc/tutorials/optical_response/triplet_excitations/1.time_propagation/3.spectrum/Singlet_triplet_spectrum_CH4.eps caption="Comparison of absorption spectrum of methane calculated with time-propagation for singlets and triplets."

You can see that there are now separate columns for cross-section and strength function for each spin. The physically meaningful strength function for the magnetic excitation is given by {{< code "StrengthFunction(1)" >}} - {{< code "StrengthFunction(2)" >}} (since the kick was opposite for the two spins). (If we had obtained {{< file "cross_section_tensor" >}}, then the trace in the second column would be the appropriate cross-section to consider.) We can plot and compare to the singlet results obtained before. You can see how this looks on the right. The first triplet transition is found at 9.05 eV, slightly lower energy than the lowest singlet transition. 

If you are interested, you can also repeat the calculation for {{< code-inline >}}{{< variable "TDDeltaStrengthMode" >}} = kick_density{{< /code-inline >}} (the default) and confirm that the result is the same as for the non-spin-polarized calculation.

### Using symmetries of non-magnetic systems

As said before, methane is a non-magnetic system, that is, the up and down densities are the same and the magnetization density is zero everywhere:

$$
 \\rho^{\\uparrow}(\\mathbf r) = \\rho^{\\downarrow}(\\mathbf r)
$$

Note that it is not enough that the total magnetic moment of the system is zero as the previous condition does not hold for anti-ferromagnetic systems. The symmetry in the spin-densities can actually be exploited in order to obtain both the singlet and triplet spectra with a single calculation. This is done by using a special perturbation that is only applied to the spin up states.[^footnote-1] 
To use this perturbation, we need to set {{< code-inline >}}{{< variable "TDDeltaStrengthMode" >}} = kick_spin_and_density{{</ code-inline >}}. If you repeat the time-propagation with this kick, you should obtain a different {{< file "cross_section_vector" >}} file containing both the singlet and triplet spectra. The singlets are given by {{< code "StrengthFunction(1)" >}} + {{< code "StrengthFunction(2)" >}}, while the triplets are given by {{< code "StrengthFunction(1)" >}} - {{< code "StrengthFunction(2)" >}}.

## Casida equation

The calculation of triplets with the {{< code casida>}} mode for spin-polarized systems is currently not implement in {{< octopus >}}. Nevertheless, just like for the time-propagation, we can take advantage that for non-magnetic systems the two spin are equivalent. In this case it allow us to calculate triplets without the need for a spin-polarized run. The effective kernels in these cases are:

$f_{\rm Hxc}^{\rm singlet} \left[ \rho \right] = f^{\uparrow}_{\rm Hxc} \left[ \rho ^{\uparrow} \right] + f^{\uparrow}_{\rm Hxc} \left[ \rho^{\downarrow} \right] = f_{\rm H} \left[ \rho \right] + f^{\uparrow}_{\rm xc} \left[ \rho ^{\uparrow} \right] + f^{\uparrow}_{\rm xc} \left[ \rho^{\downarrow} \right]$

$f_{\rm Hxc}^{\rm triplet} \left[ \rho \right] = f^{\uparrow}_{\rm Hxc} \left[ \rho ^{\uparrow} \right] - f^{\uparrow}_{\rm Hxc} \left[ \rho^{\downarrow} \right] = f^{\uparrow}_{\rm xc} \left[ \rho ^{\uparrow} \right] - f^{\uparrow}_{\rm xc} \left[ \rho^{\downarrow} \right]$

Therefore, we start by doing a ground-state and unoccupied states runs exactly as was done in the {{< tutorial "Response/Optical spectra from Casida" "Optical spectra from Casida tutorial" >}}. Then, do a Casida run with the following input file:

#include_eps doc/tutorials/optical_response/triplet_excitations/2.casida/4.spectrum/Singlet_triplet_spectrum_Casida_CH4.eps caption="Comparison of absorption spectrum of methane calculated with the Casida equation for singlets and triplets."


{{< code-block >}}
#include_input doc/tutorials/optical_response/triplet_excitations/2.casida/5.extra_states_convergence/inp
{{< /code-block >}}

The only difference with respect to the calculation of the singlets, is the {{< variable "CasidaCalcTriplet" >}} input variable.

Once the Casida calculation is finished, run {{< file "oct-casida_spectrum" >}}, plot the results, and compare to the singlet calculation. On the right you can see how this plot should look like. How do the singlet and triplet energy levels compare? Can you explain a general relation between them? How does the run-time compare between singlet and triplet, and why?

## Comparison

As for the singlet spectrum, we can compare the time-propagation and Casida results. What is the main difference, and what is the reason for it?

#include_eps doc/tutorials/optical_response/triplet_excitations/3.comparison/Triplet_casida_td_CH4.eps caption="Comparison of triplet absorption spectrum of methane calculated with time-propagation and with the Casida equation."


[^footnote-1]: {{< article title="On the use of Neumann's principle for the calculation of the polarizability tensor of nanostructures" authors="M.J.T. Oliveira, A. Castro, M.A.L. Marques, and A. Rubio" journal="J. Nanoscience and Nanotechnology" volume="8" pages="1-7" year="2008" doi="10.1166/jnn.2008.142" arxiv="0710.2624v1" >}}


{{< tutorial-footer >}}
---
title: "Optical spectra from Sternheimer"
#tags: ["Beginner", "Molecule", "Pseudopotentials", "DFT", "Optical Absorption", "Electromagnetic Response", "Sternheimer"]
tutorials: "Optical Response"
theories: "DFT"
calculation_modes: "Electromagnetic response"
difficulties: "beginner"
system_types: "Molecule"
species_types: "Pseudopotentials"
features: ["Optical absorption", "Sternheimer"]
Weight: 4
#series: "Tutorial"
description: "Calculate optical spectra in the frequency domain from linear response."
---


We have just seen how to calculate optical spectra in the time domain with a finite perturbation, and in a frequency-domain, linear-response matrix formulation with the Casida equation. Now we will try a third approach, which is in the frequency domain and linear response but rather than using a pseudo-eigenvalue equation as in Casida, uses a self-consistent linear equation, the Sternheimer equation. This approach is also known as density-functional perturbation theory. It has superior scaling, is more efficient for dense spectra, and is more applicable to nonlinear response. One disadvantage is that one needs to proceed one frequency point at a time, rather than getting the whole spectrum at once. We will find we can obtain equivalent results with this approach for the optical spectra as for time propagation and Casida, by calculating the polarizability and taking the imaginary part.

## Ground state

Before doing linear response, we need to obtain the ground state of the system, for which we can use the same input file as for [Optical spectra from Casida](../Optical spectra from Casida), but we will use a tighter numerical tolerance, which helps the Sternheimer equation to be solved more rapidly. Unlike for Casida, no unoccupied states are required. If they are present, they won't be used anyway.

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_sternheimer/1.gs/inp
{{< /code-block >}}

{{< code-block >}}
 {{< variable "CalculationMode" >}} = gs
 {{< variable "UnitsOutput" >}} = eV_angstrom
 
 {{< variable "Radius" >}} = 6.5*angstrom
 {{< variable "Spacing" >}} = 0.24*angstrom
 
 CH = 1.097*angstrom
 %{{< variable "Coordinates" >}}
   "C" |           0 |          0 |           0
   "H" |  CH/sqrt(3) | CH/sqrt(3) |  CH/sqrt(3)
   "H" | -CH/sqrt(3) |-CH/sqrt(3) |  CH/sqrt(3)
   "H" |  CH/sqrt(3) |-CH/sqrt(3) | -CH/sqrt(3)
   "H" | -CH/sqrt(3) | CH/sqrt(3) | -CH/sqrt(3)
 %
 {{< variable "ConvRelDens" >}} = 1e-7
{{< /code-block >}}

## Linear response

Add the lines below to the input file (replacing the {{< variable "CalculationMode" >}} line).

The frequencies of interest must be specified, and we choose them based on the what we have seen from the Casida spectrum. The block below specifies 5 frequencies spanning the range 0 to 8 eV (below the resonances) and 9 frequencies spanning the range 10 to 12 eV (where there are peaks). If we didn't know where to look, then looking at a coarse frequency grid and then sampling more points in the region that seems to have a peak (including looking for signs of resonances in the real part of the polarizability) would be a reasonable approach. We must add a small imaginary part ($\eta$) to the frequency in order to be able to obtain the imaginary part of the response, and to avoid divergence at resonances. The resonances are broadened into Lorentzians with this width. The larger the $\eta$, the easier the SCF convergence is, but the lower the resolution of the spectrum.

To help in the numerical solution, we turn off the preconditioner (which sometimes causes trouble here), and use a linear solver that is experimental but will give convergence much faster than the default one.

{{< code-block >}}
#include_input doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/inp
{{< /code-block >}}


## Spectrum

This Perl script can extract the needed info out of the files in the {{< file "em_resp" >}} directory:

```text
#include_file doc/tutorials/optical_response/optical_spectra_from_sternheimer/2.em_resp/extract.pl
```

Our result is

```text
 # Energy (eV)    Re alpha         Im alpha    Cross-section (A^2)
    0.000          2.64122000       0.00000000       0.00000000
    2.000          2.69965700       0.00041822       0.00007670
    4.000          2.89886900       0.00101989       0.00037410
    6.000          3.35097500       0.00234919       0.00129254
    8.000          4.70072400       0.01004299       0.00736763
   10.000          3.58189900       0.04386600       0.04022566
   10.250          4.14888500       0.11959579       0.11241259
   10.500          4.86710100       0.05001950       0.04816193
   10.750          6.78277900       0.07330539       0.07226359
   11.000         10.49267200       0.26261161       0.26489990
   11.250          3.33778500       1.06435158       1.09802650
   11.500         -0.47501400       0.20728472       0.21859505
   11.750          4.28837500       0.28012324       0.30182986
   12.000         -1.06249200       0.21176755       0.23303215
```

##  See also  

{{< tutorial "Sternheimer linear response" "Sternheimer linear response" >}}

{{< tutorial-footer >}}
---
Title: "Gaussian-shaped external current density"
tutorials: "Maxwell"
Weight: 14
---

### Gaussian-shaped external current density passed by a Gaussian temporal current pulse

#### No absorbing boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp
{{< /code-block >}}
{{< /expand >}}

Instead of an incoming external plane wave, Octopus can simulate also external current
densities placed inside the simulation box. In this example we place one shape of such
a current density in the simulation box

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp boundaries
{{< /code-block >}}

Since we start with no absorbing boundaries, we reset the box size to 10.0.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp box
{{< /code-block >}}

The external current density is switched on by the corresponding options and
two blocks define its spatial distribution and its temporal behavior. The
spatial distribution of our example external current is a Gaussian distribution
in 3D. The temporal pulse is one Gaussian along the y-axis and one along the
opposite direction but time shifted.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/1.gaussian_current_pulse/inp current
{{< /code-block >}}

{{% expand "gnuplot script" %}}
```
set pm3d
set view map
set palette defined (-0.005 "blue", 0 "white", 0.005"red")
set term png size 1000,500

unset surface
unset key

set output 'plot1.png'

set xlabel 'x-direction'
set ylabel 'y-direction'
set cbrange [-0.005:0.005]

set multiplot

set origin 0.025,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.252788 au)'
sp [-10:10][-10:10][-0.01:0.01] 'Maxwell/output_iter/td.0000120/e_field-z.z=0' u 1:2:3

set origin 0.525,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.379182 au)'
sp [-10:10][-10:10][-0.01:0.01] 'Maxwell/output_iter/td.0000180/e_field-z.z=0' u 1:2:3

unset multiplot
```
{{% /expand %}}

Contour plot of the electric field in z-direction after 120 time steps for
t=0.24 and 180 time steps for t=0.36:
{{< figure src="/images/Maxwell/tutorial_04.1-plot1.png" width="50%" >}}



#### Mask absorbing boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/3.external-current/2.gaussian_current_pulse_with_mask/inp
{{< /code-block >}}
{{< /expand >}}

We can add now mask absorbing boundaries.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/2.gaussian_current_pulse_with_mask/inp boundaries
{{< /code-block >}}

Accordingly to the additional absorbing width, we have to update the simulation
box dimensions.
{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/2.gaussian_current_pulse_with_mask/inp box
{{< /code-block >}}


Contour plot of the electric field in z-direction after 120 time steps for
t=0.24 and 180 time steps for t=0.36:
{{< figure src="/images/Maxwell/tutorial_04.2-plot1.png" width="50%" >}}

#### PML boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/3.external-current/3.gaussian_current_pulse_with_pml/inp
{{< /code-block >}}
{{< /expand >}}


We can repeat the simulation using PML absorbing boundaries.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/3.external-current/3.gaussian_current_pulse_with_pml/inp boundaries
{{< /code-block >}}

Contour plot of the electric field in z-direction after 120 time steps for
t=0.24 and 180 time steps for t=0.36:
{{< figure src="/images/Maxwell/tutorial_04.3-plot1.png" width="50%" >}}


Maxwell fields at the origin and Maxwell energy inside the free Maxwell
propagation region of the simulation box:
{{< figure src="/images/Maxwell/tutorial_04_maxwell_energy_and_fields.png" width="50%" >}}

<!---
Contour plot of the electric field in z-direction after 120 time steps for t=0.24 and 180 time steps for t=0.36:
{{< figure src="/images/Maxwell/tutorial_08_run_electric_field_contour.png" width="720px" >}}
--->

{{< tutorial-footer >}}
---
Title: "Maxwell input file"
tutorials: "Maxwell"
weight: 4
---

### Maxwell Input File

#### Input file variable description

##### Calculation mode and parallelization strategy


At the beginning of the input file, the basic option variable {{< variable
"CalculationMode" >}} selects the run mode of Octopus and has always to be set.
In case of a parallel run, there are some variables to set the proper
parallelization options. For Maxwell propagation, parallelization in domains is
possible, but parallelization in states is not needed, as there are always 3 or
6 states, depending on the Hamiltonian (see below).

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp calc_mode
{{< /code-block >}}


##### Multisystem setup

The Maxwell system is now implemented in Octopus' multisystem framework, which
allows to calculate several different systems, which can interact with each
other.

Currently implemented system types are:

* electronic: An electronic system. (only partly implemented)
* maxwell: A maxwell system.
* classical_particle: A classical particle. Used for testing purposes only.
* charged_particle: A charged classical particle.
* multisystem: A system containing other systems.
* dftbplus: A DFTB+ system
* linear_medium: A linear medium for classical electrodynamics.


##### Definition of the systems

In this tutorial, we will use a pure Maxwell system:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp system
{{< /code-block >}}

Subsequently, variables relating to a named system are prefixed by the system
name. However, this prefix is not mandatory for calculations considering only
one system.

##### Maxwell box variables and parameters

The Maxwell box is only defined as a parallelepiped box.

{{< expand "Example" >}}

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp box
{{< /code-block >}}
{{< /expand >}}


##### Maxwell options

The Maxwell options determine the Maxwell propagation scheme. First, the
Maxwell propagation operator options consist of the type of {{< variable
"MaxwellHamiltonianOperator" >}}, which can be set as a vacuum
Maxwell-Hamiltonian (faraday_ampere) or a Maxwell-Hamiltonian in linear medium
(faraday_ampere_medium).

Currently, the Maxwell system does not define its own propagator, but uses the
system default propagator, defined in {{< variable "TDSystemPropagator" >}}. So
far, only the exponential midpoint is implemented.

The Maxwell Boundary conditions can be chosen for each direction differently
with {{< variable "MaxwellBoundaryConditions" >}}. Possible choices are zero and
absorbing boundary conditions, as well as plane waves and constant field
boundary conditions. Additionally to the boundary condition, the code can
include absorbing boundaries. This means that all outgoing waves are absorbed
while all incoming signals still arise at the boundaries and propagate into
the simulation box. The absorbing boundaries can be achieved by a mask
function or by a perfectly matched layer (PML) calculation with additional
parameters. For more information on the physical meaning of these parameters,
please refer to the [Maxwell-TDDFT paper].

[Maxwell-TDDFT paper]: https://doi.org/10.1080/00018732.2019.1695875

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp maxwell_calc
{{< /code-block >}}

##### Output options

The output option {{< variable "OutputFormat" >}} can be used exactly as in the
case of TDDFT calculations, with the exception that some formats are inconsistent
with a Maxwell calculation (like xyz). It is possible to chose multiple
formats. The Maxwell output options are defined in the {{< variable
"MaxwellOutput" >}} block (written to output_iter every {{< variable
"MaxwellOutputInterval" >}} steps) and through the {{< variable
"MaxwellTDOutput" >}} (written to td.general every step).

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp output2
{{< /code-block >}}


##### Time step variables

The {{< variable "TDTimeStep" >}} option in general defines the time step used
for each system. For the stability of the time propagation, it is recommended for
Maxwell systems to always fulfill the Courant condition. In some cases larger
time steps are possible but this should be checked case by case.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp timestep
{{< /code-block >}}


##### Maxwell field variables

The incident plane waves can be defined at the boundaries using the {{<
variable "MaxwellIncidentWaves" >}} block. This block defines the incoming wave
type and complex electric field amplitude, and a Maxwell envelope function
name, which must be defined in a separate block. Multiple plane waves can be
defined, and their superposition gives the final incident wave result at the
boundary. This plane waves boundaries are evaluated at each timestep, and give
the boundary conditions for the propagation inside the box. To use them {{<
variable "MaxwellBoundaryConditions" >}} must be set to "plane_waves" in the
necessary directions.

The {{< variable "MaxwellFunctions" >}} block reads the additional required
information for the plane wave pulse, for each envelope function name given in
the {{< variable "MaxwellIncidentWaves" >}} block.

Inside the simulation box, the initial electromagnetic field is set up by the
{{< variable "UserDefinedInitialMaxwellStates" >}} block. If the pulse
parameters are such that it is completely outside the box at the initial time,
this variables does not need to be set, as the default initial field inside the
box is zero.

In case of {{< code-inline >}}{{< variable "MaxwellPlaneWavesInBox" >}} = yes
{{< /code-inline >}}, the code uses the initial analytical plane wave values
also inside the simulation box. Hence, there is no Maxwell field propagation
and the values are set analytically.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp field
{{< /code-block >}}

##### External current

An external current can be switched that adds an external current contribution
to the internal current. Such an external current is calculated analytically by
the {{< variable "UserDefinedMaxwellExternalCurrent" >}} block. This block
defines the spatial profile and frequency of each external current
contribution, and it sets an envelope function name, that must be defined in a
separate {{< variable "TDFunctions" >}} block.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/00-maxwell-input-options.inp ext_current
{{< /code-block >}}


##### Linear Medium

Different shapes of linear media can be included in the calculation either
through a simple box, or through a file describing more complex geometries. For
more information check the Linear Medium tutorial.


{{< tutorial-footer >}}

---
Title: "Maxwell overview"
tutorials: "Maxwell"
weight: 1
---


## Tutorial for Maxwell propagation and its coupling with TDDFT in Octopus

#### Authors: Rene Jestädt, Franco Bonafé and Heiko Appel

#### Max Planck Institute for the Structure and Dynamics of Matter - Theory Department


This tutorial is based on the following publication:

[Light-matter interactions within the Ehrenfest–Maxwell–Pauli–Kohn–Sham
framework: fundamentals, implementation, and nano-optical
applications](./Light_matter_interactions_within_the_Ehrenfest_Maxwell_Pauli_Kohn_Sham_framework_fundamentals_implementation_and_nano_optical_applications.pdf)
René Jestädt, Michael Ruggenthaler, Micael J. T. Oliveira, Angel Rubio, and
Heiko Appel https://doi.org/10.1080/00018732.2019.1695875

---

[Slides of Heiko Appel's talk](https://theory.mpsd.mpg.de/talks/heiko.appel/2020-01-21-Uni-Jena)

---

Using the Riemann-Silberstein representation for the Maxwell's equations, the
corresponding time-evolution of Maxwell fields is implemented as quantum
mechanical like propagation in the TDDFT code Octopus.

The program can be run in a free Maxwell propagation mode where the
electromagnetic fields are calculated inside a simulation box. This box can
include various shapes for linear media. An initial electromagnetic field can be
set up in the box as well as an external current density or incoming plane
waves at the box boundaries. The simulation box boundaries can be selected as
zero boundaries, perfectly electric conductor (PEC mirror), perfectly magnetic
conductor (PMC) or absorbing boundaries with different methods.

Combining all features of a free Maxwell propagation with the normal TDDFT
matter propagation in Octopus, the code can solve fully coupled
Maxwell-Kohn-Sham systems by propagating both systems on two separated grids.
In this case, the Kohn-Sham system gives rise to an internal current density that
influences the Maxwell field propagation, and in turn the electromagnetic field
is part of the Kohn-Sham Hamiltonian. Both coupled systems are propagated
self-consistently, but the code can also be used for only forward or only
backward coupling.

In the following we introduce the different input options for the free Maxwell
propagation and the fully coupled Maxwell-Kohn-Sham propagation, and different
examples to show the various features.

<!-- All the input files can be downloaded [here](./maxwelltddft_tutorial_files.tar). -->

---

#### Simulation box and relevant input file variables

* {{< tutorial "maxwell/simulationbox" "Simulation box" >}}
* {{< tutorial "maxwell/maxwellinputfile" "Input file variables" >}}




#### Examples of electromagnetic field propagation without coupling to TDDFT

* {{< tutorial "maxwell/run01" "Cosinoidal plane wave in vacuum" >}}
* {{< tutorial "maxwell/run02" "Interference of two cosinoidal plane waves" >}}
* {{< tutorial "maxwell/run03" "Cosinoidal plane wave hitting a linear medium box with and without absorbing boundaries" >}}
* {{< tutorial "maxwell/run04" "Gaussian-shaped spatial external current distribution density passed by a Gaussian temporal current pulse" >}}


<!--
#### Coupled Maxwell-TDDFT propagation

* {{< tutorial "maxwell/benzene" "Benzene ground state and dynamics coupled to an external EM pulse" >}}


To check what will be the input syntax once the Maxwell propagation branch is merged into the current version of Octopus and how to obtain and compile the code, click on {{< tutorial "multisystem" "Multisystems" >}}

-->
{{< tutorial-footer >}}

---
Title: " Coupled Maxwell-TDDFT propagation"
tutorials: "Maxwell"
Weight: 20
draft: yes
---

### Coupled Maxwell-TDDFT propagation 

{{< notice warning >}}
This needs to be removed as Maxwell-matter coupling is not yet working.{{< /notice >}}

#### Ground state of benzene


To simulate a coupled Maxwell-Kohn-Sham system, we have to obtain first an initial Kohn-Sham system. Here in our case a groundstate of benzene. 

```
# ----- Calculation node and parallelization ----------------------------------------------------

CalculationMode = gs
ParDomains = auto
ParStates = no

# ----- Matter box variables --------------------------------------------------------------------

lsize_ma = 10.0
dx_ma = 0.5

Dimensions = 3
BoxShape = parallelepiped

%Lsize
lsize_ma | lsize_ma | lsize_ma
%

%Spacing
dx_ma | dx_ma | dx_ma
%

# ----- Species variables -----------------------------------------------------------------------

XYZCoordinates = "geometry_z_dir.bohr.xyz"

# ----- Hatter calculation variables

Extrastates = 2
XCFunctional = lda_x + lda_c_gl
EigenSolver = cg
EigenSolverTolerance = 1.00e-12
EigensolverMaxIter = 50
CoanelDens = 1.00e-11
MaximumIter = 1000
ConvForce = 0.0
SmearingFunction = fermi_dirac
Smearing = 0.001
MixingScheme = broyden
Mixing = 0.02
LCAOStart = lcao_states
```
#### Coupled dynamics of benzene and a cosinoidal EM pulse 

** Folder [/02_maxwell_ks_propagation/01_ed_coupling/01_no_backreaction/](./02_maxwell_ks_propagation/01_ed_coupling/01_no_backreaction/) **

A cosinoidal laser pulse hits a molecule in groundstate usind electric dipole approximation without any matter to Maxwell back-reaction. 

```
# ----- Calculation mode and parallelization ----------------------------------------------------

CalculationMode = maxwell_ks
ParDomains = auto
ParStates = no
MaxwellParDomains = auto
MaxwellParStates = no

# ----- Matter box variables --------------------------------------------------------------------

lsize_ma = 10.0
dx_ma = 0.5
Dimensions = 3
BoxShape = parallelepiped

%Lsize
lsize_ma | lsize_ma | lsize_ma
%

%Spacing
dx_ma | dx_ma | dx_ma
%

# ----- Maxwell box variables -------------------------------------------------------------------

# free maxwell box limit of 10.0 plus for the 5.0 absorbing pml boundaries 
# plus 2.0 for the incident wave boundaries with der_order = 4 times dx_mx
lsize_mx = 17.0 
dx_mx = 0.5
MaxwellDimensions = 3
MaxwellBoxShape = parallelepiped

%MaxwellLsize
lsize_mx | lsize_mx | lsize_mx
%

%MaxwellSpacing
dx_mx | dx_mx | dx_mx
%

# ----- Species variables -----------------------------------------------------------------------

XYZCoordinates = "geometry_z_dir.bohr.xyz"

# ----- Matter calculation variables ------------------------------------------------------------

Extrastates = 2
XCFunctional = lda_x + lda_c_gl
EigenSolver = cg
EigenSolverTolerance = 1.00e-12
EigensolverMaxIter = 50
CoanelDens = 1.00e-12
MaximumIter = 1000
DerivativesStencil = stencil_starplus
DerivativesOrder = 4
TDExpOrder = 4
AbsorbingBoundaries = not_absorbing

# ----- Maxwell calculation variables -----------------------------------------------------------

MaxwellHamiltonianOperator = faraday_ampere
MaxwellTDOperatorMethod = maxwell_op_fd

MaxwellTDPropagator = maxwell_etrs
MatterToMaxwellCoupling = no
MaxwellToMatterCoupling = yes

MaxwellCouplingOrder = electric_dipole_coupling
MaxwellTransFieldCalculationMethod = trans_field_poisson
MaxwellPoissonSolver = isf
MaxwellPoissonSolverBoundaries = multipole

MaxwellDerivativesStencil = stencil_starplus
MaxwellDerivativesOrder = 4
MaxwellTDExpOrder = 4

%MaxwellBoundaryConditions
maxwell_plane_waves | maxwell_plane_waves | maxwell_plane_waves
%

%MaxwellAbsorbingBoundaries
cpml | cpml | cpml
%

MaxwellABPMLWidth = 5.0
MaxwellABPMLKappaMax = 1.0
MaxwellABPMLAlphaMax = 1.0 
MaxwellABPMLPower = 2.0
MaxwellABPMLReflectionError = 1.0e-16

# ----- Output variables ------------------------------------------------------------------------

OutputFormat = plane_x + plane_y + plane_z + vtk + xyz + axis_x

# ----- Matter output variables -----------------------------------------------------------------

Output = potential + density + current + geometry + forces + elf

OutputInterval = 50
TDOutput = energy + multipoles + laser + geometry

# ----- Maxwell output variables ----------------------------------------------------------------

MaxwellOutput = maxwell electric field + maxwell magnetic field + maxwell_energy_density + maxwell_trans_§lectric_field #_(has to be written in one line)

MaxwellOutputInterval = 1

MaxwellTDOutput = maxwell_energy + maxwell_fields

%MaxwellFieldsCoordinate
0.00 | 0.00 | 0.00
%

# ----- Time step variables ---------------------------------------------------------------------
TDTimeStep = 0.002
TDMaxSteps = 200

TDEnergyUpdateIter = 1
MaxwellTDIntervalSteps = 1
TDMaxwellTDRelaxationSteps = 0
TDMaxwellKSRelaxationSteps = 0

MaxwellTDETRSApprox = no
CurrentPropagationTest = no

# ----- Maxwell field variables -----------------------------------------------------------------

lambda = 10.0
omega = 2 * pi * c / lambda
kx = omega / c
E2 = 0.05
pw = 10.0
ps = - 5 * 5.0

%UserDefinedMaxwellIncidentWaves
plane_waves_mx_function | 0 | 0 | E2 | "plane_waves_function" | plane_wave
%

%MaxwellFunctions
"plane_waves_function" | mxf_cosinoidal_wave | kx | 0 | 0 | ps | 0 | 0 | pw
%

# cosinoidal pulse
%UserDefinedInitialMaxwellStates
3 | formula | electric_field | “ Ez*cos(kx*(x-ps))*cos(pi/2*(x-ps-2*pw)/pw) * step(pw-abs((kx*x-kx*ps)/kx^2)) “ # (in one line)
2 | formula | magnetic_field | “ -1/c*Ez*cos(kx*(x-ps))*cos(pi/2*(x-ps-2*pw)/pw) * step(pw-abs((kx*x-kx*ps)/kx^2)) " # (in one line)
%
```

To explore the effect of the back reaction of the field (considering the fully coupled Maxwell-matter interaction) and the effect of higher order coupling terms (magnetic dipole and electric quadrupole) run the input files in folders {{< file "02_maxwell_ks_propagation/01_ed_coupling/02_fully_coupled/" >}} , {{< file "2_maxwell_ks_propagation/02_ed_md_eq_coupling/01_no_backreaction/" >}} and {{< file "02_maxwell_ks_propagation/02_ed_md_eq_coupling/01_no_backreaction/" >}} . 

{{< tutorial-footer >}}

---
title: "Maxwell Systems"
weight: 60
description: "Maxwell Systems"
---

Since version 10, {{< octopus >}} can also be used to calculate Maxwell systems, i.e. propagating the electro-magnetic fields.
In future versions, this will be coupled to matter systems, to go beyond the usual forward-only coupling between the fields and the electrons.

{{< tutorial-series "Maxwell" >}}---
Title: "Interference of two cosinoidal plane waves"
tutorials: "Maxwell"
Weight: 12
---

## Interference of two cosinoidal plane waves


Instead of only one plane wave, we simulate two different plane waves with
different wave-vectors entering the simulation box, interfering and leaving the
box again. In addition to the wave from the last tutorial, we add a second
wave with different wave length, and entering the box at an angle, and shifted
by 28 Bohr along the corresponding direction of propagation.

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/1.free-propagation/2.2_pulses_td/inp
{{< /code-block >}}
{{< /expand >}}

Both electric fields are polarized only in z-direction, and the magnetic field
only in y-direction.

We can start from the last input file, and add the second wave, according to
the following excerpt:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/2.2_pulses_td/inp field
{{< /code-block >}}

Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:
{{< figure src="/images/Maxwell/tutorial_02_run_electric_field_contour.png" width="50%" >}}

Maxwell fields at the origin and Maxwell energy inside the free Maxwell
propagation region of the simulation box:
{{< figure src="/images/Maxwell/tutorial_02_run_maxwell_energy_and_fields.png" width="50%" >}}

{{< tutorial-footer >}}
---
Title: "Cosinoidal plane wave in vacuum"
tutorials: "Maxwell"
Weight: 11
---

## Cosinoidal plane wave in vacuum

In this tutorial, we will describe the propagation of a cosinoidal wave in vacuum.

{{% expand "click for full input file" %}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp
{{< /code-block >}}
{{% /expand %}}

{{% notice warning %}}
  Currently, this run only works with 1 or 2 MPI tasks.
{{% /notice %}}

First, define the calculation mode and the parallelization strategy:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp calc_mode
{{< /code-block >}}


The total size of the (physical) simulation box is 10.0 x 10.0 x 10.0 Bohr with
a spacing of 0.5 Bohr in each direction. As discussed in the {{< tutorial
"maxwell/simulationbox" "section on simulation boxes" >}}, also the boundary
points have to be accounted for. Hence the total size of the simulation box is
chosen to be 12.0 x 12.0 x 12.0 Bohr.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp box
{{< /code-block >}}

Since we propagate without a linear medium, we can use the vacuum Maxwell
Hamiltonian operator. Derivatives and exponential expansion order are four.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp calculation
{{< /code-block >}}

The boundary conditions are chosen as plane waves to simulate the incoming wave
without any absorption at the boundaries.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp boundaries
{{< /code-block >}}

The simulation time step is set to 0.002 and the system propagates to a maximum
time step of 200. For each time step, we write an output of the electric and
the magnetic field as well as the energy density, the total Maxwell energy.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp timestep
{{< /code-block >}}

The incident cosinoidal plane wave passes the simulation box, where the
cosinoidal pulse has a width of 10.0 Bohr, a spatial shift of 25.0 Bohr in the
negative x-direction, an amplitude of 0.05 a.u., and a wavelength of 10.0 Bohr.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp field
{{< /code-block >}}


Finally, the output options are set:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/1.free-propagation/1.1_pulse_td/inp output
{{< /code-block >}}

Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:

{{% expand "Example gnuplot script" %}}
```bash
set pm3d
set view map
set palette defined (-0.1 "blue", 0 "white", 0.1 "red")
set term png size 1000,500

unset surface
unset key

set output 'plot.png'

set xlabel 'x-direction'
set ylabel 'y-direction'
set cbrange [-0.1:0.1]

set multiplot

set origin 0.025,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.105328 au)'
sp [-10:10][-10:10][-0.1:0.1] 'Maxwell/output_iter/td.0000050/e_field-z.z=0' u 1:2:3

set origin 0.525,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.210656 au)'
sp [-10:10][-10:10][-0.1:0.1] 'Maxwell/output_iter/td.0000100/e_field-z.z=0' u 1:2:3

unset multiplot
```
If you copy this script into a file called {{< file "plot.gnu" >}}, you can create the plots by:
```bash
gnuplot plot.gnu
```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_01_run_electric_field_contour.png" width="50%" >}}


Maxwell fields at the origin and Maxwell energy inside the free Maxwell
propagation region of the simulation box:

{{% expand "gnuplot script" %}}
```bash
set term png size 1000,400
set output 'plot2.png'

set multiplot
set style data l

set xlabel 'time [a.u.]'

set origin 0.025,0
set size 0.45,0.9
set title 'Electric and magnetic fields'
set ylabel 'Maxwell fields [a.u.]'
p 'Maxwell/td.general/total_e_field_z' u 2:3 t 'E_z', 'Maxwell/td.general/total_b_field_y' u 2:($3*100) t '100 * B_y'

set origin 0.525,0
set size 0.45,0.9
set title 'Maxwell energy'
set ylabel 'Maxwell energy [a.u.]'
p  'Maxwell/td.general/maxwell_energy' u 2:3  t 'Maxwell energy'

unset multiplot

```
{{% /expand %}}

{{< figure src="/images/Maxwell/tutorial_01_run_maxwell_energy_and_fields.png" width="50%" >}}

{{< tutorial-footer >}}
---
Title: "Cosinoidal plane wave hitting a linear medium box"
tutorials: "Maxwell"
Weight: 13
---

## Cosinoidal plane wave hitting a linear medium box

An arbitrary number of linear medium shapes can be placed inside the Maxwell
simulation box.

Linear media are considered a separate system type, for example:

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp systems
{{< /code-block >}}


The box shape can be defined in two ways, which are given by the variable
{{<variable "LinearMediumBoxShape">}}: if set to medium_parallelepiped, the
parallelepiped box will be defined by its center, and size in each dimension
through the {{<variable "LinearMediumBoxSize">}} block; if the box shape is
defined as medium_box_file, the box shape will be read from an external OFF
file, defined through the variable {{<variable "LinearMediumBoxFile">}}. To
produce such files, check the {{<tutorial "Maxwell/details/creating_geometries" 
"Creating geometries" >}} tutorial. In either case, the
electromagnetic properties of each medium must be defined in the {{<variable
"LinearMediumProperties">}} block, which specifies the relative electric
permittivity, relative magnetic permeability, and the electric and magnetic
conductivities (for lossy media). Finally, the spatial profile assumed to
calculate gradients at the box edges is defined through the {{<variable
"LinearMediumEdgeProfile">}} variable, and can be either "edged" or "smooth"
(for linear media shape read from file, it can only be edged).

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp medium_box
{{< /code-block >}}

Note that we also need to change the variables describing the propagation.
{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp timestep
{{< /code-block >}}

{{< notice note >}}
Linear media are static. Therefore only the propagator
{{<code-inline>}}{{<variable "TDSystemPropagator">}} = static
{{</code-inline>}} is permitted.  By moving the variable {{<variable
"TDSystemPropagator">}} into the {{< code "Maxwell" >}} namespace, the default
propagator (which is {{<code "static">}}) is automatically inherited by the
{{<code "Medium">}} system and only the propagator for the {{<code
"Maxwell">}} system is explicitly set to {{<code "exp_mid">}}. As the {{<code
"exp_mid">}} propagator contains two algorithmic steps, it must be clocked
twice as fast as the {{<code "static">}} propagator and hence the time step 
of the medium bust be set to half the Maxwell time step. This is currently a
workaround, and is likely to change in future releases of the code.
{{< /notice >}}

### No absorbing boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp
{{< /code-block >}}
{{< /expand >}}


To illustrate a laser pulse hitting a linear medium box, in addition to
adding the linear_medium system, we need to switch the Maxwell Hamiltonian
operator from the 3x3 vacuum representation to the 6x6 coupled representation
with the six component Riemann-Silberstein vector that includes also the complex
conjugate Riemann-Silberstein vector.


{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp calculation
{{< /code-block >}}

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp box
{{< /code-block >}}

For this run, we use the previous incident plane wave propagating only in the
x-direction but place a medium box inside the Simulation box.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/1.cosinoidal_pulse_td/inp field
{{< /code-block >}}

{{% expand "gnuplot script" %}}
```
set pm3d
set view map
set palette defined (-0.05 "blue", 0 "white", 0.05"red")
set term png size 1000,500

unset surface
unset key

set output 'plot1.png'

set xlabel 'x-direction'
set ylabel 'y-direction'
set cbrange [-0.05:0.05]

set multiplot

set origin 0.025,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.105328 au)'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000050/e_field-z.z=0' u 1:2:3

set origin 0.525,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.210656 au)'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000100/e_field-z.z=0' u 1:2:3

unset multiplot

set output 'plot2.png'

set multiplot

set origin 0.025,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.263321 au)'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000125/e_field-z.z=0' u 1:2:3

set origin 0.525,0
set size 0.45,0.9
set size square
set title 'Electric field E_z (t=0.315985 au)'
sp [-10:10][-10:10] 'Maxwell/output_iter/td.0000150/e_field-z.z=0' u 1:2:3

unset multiplot

```
{{% /expand %}}

Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:
{{< figure src="/images/Maxwell/tutorial_03.1-plot1.png" width="50%" >}}

Contour plot of the electric field in z-direction after 125 time steps for
t=0.26 and 150 time steps for t=0.32:
{{< figure src="/images/Maxwell/tutorial_03.1-plot2.png" width="50%" >}}

In the last panel, it can be seen that there is a significant amount of
scattered waves which, in large parts, are scattered from the box boundaries.
In the following we will use different boundary conditions, in order to reduce
this spurious scattering.


### Mask absorbing boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/2.linear-medium/2.cosinoidal_pulse_td_mask/inp
{{< /code-block >}}
{{< /expand >}}

The previous run was without absorbing boundaries, now we switch on the mask
function that damps the field at the boundary by multiplying a scalar mask
function. The {{<variable "MaxwellAbsorbingBoundaries">}} block is updated to
run with mask absorbing. The {{<variable "MaxwellABMaskWidth">}} is set to
five.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/2.cosinoidal_pulse_td_mask/inp boundaries
{{< /code-block >}}

As a consequence of the additional region for the absorbing boundary condition,
we have to change the box size to still obtain the same free Maxwell
propagation box. Therefore, the {{<code "lsize_mx">}} value is now 17.0, which
is the previously used size that includes the incident wave boundary width plus
the absorbing boundary width. For more details, see the information on 
{{<tutorial "Maxwell/details/simulationbox" "simulation boxes">}}

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/2.cosinoidal_pulse_td_mask/inp box
{{< /code-block >}}

Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:
{{< figure src="/images/Maxwell/tutorial_03.2-plot1.png" width="50%" >}}

Contour plot of the electric field in z-direction after 125 time steps for
t=0.26 and 150 time steps for t=0.32:
{{< figure src="/images/Maxwell/tutorial_03.2-plot2.png" width="50%" >}}

It can be seen that the scattering is slightly reduced, but still noticeable.

### Perfectly matched layer boundaries

{{< expand "click for complete input" >}}
{{< code-block >}}
#include_input doc/tutorials/maxwell/2.linear-medium/3.cosinoidal_pulse_td_pml/inp
{{< /code-block >}}
{{< /expand >}}


The mask absorbing method can be replaced by the more accurate perfectly
matched layer (PML) method. Therefore, the {{<variable
"MaxwellAbsorbingBoundaries">}} by the {{<code "cpml">}} option. The PML
requires some additional parameters to the width for a full definition.

{{< code-block >}}
#include_input_snippet doc/tutorials/maxwell/2.linear-medium/3.cosinoidal_pulse_td_pml/inp boundaries
{{< /code-block >}}


Contour plot of the electric field in z-direction after 50 time steps for
t=0.11 and 100 time steps for t=0.21:
{{< figure src="/images/Maxwell/tutorial_03.3-plot1.png" width="50%" >}}

Contour plot of the electric field in z-direction after 125 time steps for
t=0.26 and 150 time steps for t=0.32:
{{< figure src="/images/Maxwell/tutorial_03.3-plot2.png" width="50%" >}}

The PML boundary conditions in this case are comparable to the mask boundary conditions,
but lead to an increased computational cost.

{{< tutorial-footer >}}
---
Title: "Creating geometries"
tutorials: "Maxwell_details"
Weight: 100
---


### Creating a geometry for the linear medium

Octopus uses the [CGAL](https://www.cgal.org) library, which allows to parse
standard geometry files, as used e.g. for 3D printing. Currently, we support
only the [OFF](https://en.wikipedia.org/wiki/OFF_(file_format)) file format.

There are many software packages around to generate such files, but we
recommend [openSCAD](https://www.openscad.org/index.html) and `ctmconv` which
is part of [OpenCTM](https://openctm.sourceforge.net/).

{{< notice info >}}
    In Debian style systems, these can be installed using
    {{< code-block >}}
    apt-get install openscad
    apt-get install openctm-tools
    {{< /code-block >}}
{{< /notice >}}

OpenSCAD provides its own scripting language for creating geometries, which can
be created in a GUI, but also can be rendered on the command line.  The
documentation of OpenSCAD can be found
[here](https://www.openscad.org/documentation.html).

{{< expand "Example: creating a simple lens with Openscad" >}}
<br>
Here is an example script for a lens, generated from the intersection of two spheres:
{{< code-block >}}
$fs=0.5;
$fa=0.5;
intersection(){
  translate([-22,0,0]) sphere(r=25);
  translate([ 22,0,0]) sphere(r=25);
};
{{< /code-block >}}
For details on the syntax of this, see the OpenSCAD {{< link "https://www.openscad.org/documentation.html" "tutorials" >}}.

This creates the following lens: {{< figure src="/images/Maxwell/Lense.png" width="500px" >}}

{{< /expand >}}

If you copy the above code into a file `lens.scad`, you can generate the required `OFF` file by the command:
```bash
openscad -o lens.off lens.scad
```

{{< notice warning >}}
    The `OFF` files, generated by openSCAD are malformed. If you want to further manipulate them with `ctmconv`, the header needs to be fixed:
    {{< code-block >}}
    OFF 1162 2320 0
    1.81312 4.79828 5.90861
    1.73416 3.80751 6.86969
    1.98506 3.84776 5.90861
    1.72673 7.81373 -0.993913
    1.39732 8.75085 -0.993913
    1.41584 8.75777 0
    ...
    {{< /code-block >}}
    The first line needs to be `OFF` only, and the first line of the data block needs to contain the number of vertices, planes and (optional) connections.
    You can add comments (starting with `#` of empty lines are allowed after the header).
    {{< code-block >}}
    OFF
    # This is a lens (generated with OpenSCAD)

    1162 2320 0
    1.81312 4.79828 5.90861
    1.73416 3.80751 6.86969
    1.98506 3.84776 5.90861
    1.72673 7.81373 -0.993913
    1.39732 8.75085 -0.993913
    1.41584 8.75777 0
    ...
    {{< /code-block >}}
{{< /notice >}}




---
Title: "Simulation box"
tutorials: "Maxwell_details"
Weight: 101
---


### Simulation box

The Maxwell simulation box is divided into mainly two regions, one inner region
for the requested Maxwell propagation and an outer region to simulate the
proper boundary conditions. Later in the input file, the simulation box sizes
refer to the total simulation box which includes the boundary region. The
inner simulation box is defined by the total box size minus the boundary
width. In case of zero boundary condition, there is no boundary region. {{<
figure caption="Boundaries in 2D" src="/images/Maxwell/boundaries_in_2D.png"
width="500px" >}}

The boundary region can be set up by absorbing boundaries, incident plane waves
or a combination of both. In the latter case, the absorbing boundary region is
the inner boundary region, and the plane waves region is the outer boundary
region. Again, the box sizes are determined by the total simulation box size
and the corresponding boundary width. {{< figure caption="Plane Waves and PML"
src="/images/Maxwell/plane_wave_and_pml_boundaries_in_2D.png" width="500px" >}}

{{< notice note >}} For a given size of the propagation region, the boundaries
have to be added to the total box size. The boundary width is given by the
derivative order (default is 4) times the spacing. The width of the absorbing
boundary region is defined by the variable {{< variable "MaxwellABWidth" >}}.
{{< /notice >}}

The matter grid is in general located inside the Maxwell grid. There are
several possible types of grids to describe a coupled Maxwell-matter system.
The following figure illustrates some possible cases overlaying Maxwell and matter
grids. In the Octopus code, currently only the types e), f), and g) are
implemented. So the matter box sizes and Maxwell box sizes can be chosen
independently, whereas the spacings of both grids have to be equal and the grid
points have to lie on the top of each other. The only exception is type g),
where the matter grid is much smaller than the Maxwell grid. In this case, the
matter grid size has to be smaller than the Maxwell grid spacing. {{< figure
caption="Multiscale" src="/images/Maxwell/multiscale_figures.png" width="500px"
>}}

=====
DFTD3
=====

This is a repackaged version of the `DFTD3 program
<http://www.thch.uni-bonn.de/tc/index.php?section=downloads&subsection=getd3>`_
by S. Grimme and his coworkers.

The original program (V3.1 Rev 1) was downloaded at 2016-04-03. It has been
converted to free format and encapsulated into modules. The source has been
split into two parts:

* A library with the core functionality. This can be directly used by third
  party applications wishing to calculate dispersion with the DFT-D3
  approach.
  
* Additional extensions which are necessary for the command line tool DFTD3 and
  the command line tool itself.

* Updated dftd3 code to include refitted/modified zero- and BJ-damped D3
  versions of Sherrill and coworkers (-bjm and -zerom)
  (Functionality corresponds to V3.2 Rev0)

Compilation
===========

Edit the file `make.arch` to reflect your compiler and linker. Then you can
issue one of the following commands:

* ``make lib``: to build the library `libdftd3.a` and the necessary
  module files (`*.mod`) in the directory `lib/`.

* ``make dftd3``: to build the executable `dftd3` in the directory `prg/`.

* ``make testapi``: to build a simple tester for the library (`testapi`) in the
  directory `test/`. The source code of this tester demonstrates how the library
  can be used by third party codes.

If you just issue ``make``, all three targets will be compiled.


Credits
=======

When using the library or the dftd3 tool, please cite:

  S. Grimme, J. Antony, S. Ehrlich and H. Krieg
  J. Chem. Phys, 132 (2010), 154104.
 
If BJ-damping is used 

  S. Grimme, S. Ehrlich and L. Goerigk
  J. Comput. Chem, 32 (2011), 1456-1465.

should be cited as well.


License
=======

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 1, or (at your option) any later version.
==========
Change Log
==========


0.9
===

Added
-----

* Implementation of revised damping parameters as desribed in
  D. G. A. Smith, L. A. Burns, K. Patkowski, and C. D. Sherrill
  J. Phys. Chem. Lett., 2016, 7, pp 2197–2203.
  (Functionality should correspond to V3.2 Rev 0 of original dftd3 code.)

Fixed
-----

* Routine dftd3_pbc_dispersion delivers now correct stress tensor.


0.1
===

Added
-----

* Interface (API) for third party codes.

* Tester for the API.

* Default settings for various compilers.


Changed
-------

* Original source file converted to free format using Metcalf's convert
  tool.

* Code packed into modules.

* Core library separated from dftd3 tools.
spglib
======

|Version Badge| |Downloads Badge| |PyPi downloads| |Build Status|

Python bindings for C library for finding and handling crystal
symmetries

Installation
------------

The package is developed on github. You can get the source for the
released versions from the
`repository <https://github.com/atztogo/spglib/releases>`__.

Using package distribution service
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The easiest way to install python-spglib is to use the pypi package, for
which numpy is required to be installed before the installation. A
command to install spglib is:

::

    pip install spglib

Conda is another choice:

::

    conda install -c atztogo spglib

These packages are made by Paweł T. Jochym.

Building using setup.py
~~~~~~~~~~~~~~~~~~~~~~~

To manually install python-spglib using ``setup.py``, python header
files (python-dev), C-compiler (e.g., gcc, clang), and numpy are
required before the build. The installation steps are shown as follows:

1. Go to the python directory
2. Type the command:

   ::

       python setup.py install --user

.. |Version Badge| image:: https://anaconda.org/atztogo/spglib/badges/version.svg
   :target: https://anaconda.org/atztogo/spglib
.. |Downloads Badge| image:: https://anaconda.org/atztogo/spglib/badges/downloads.svg
   :target: https://anaconda.org/atztogo/spglib
.. |PyPi downloads| image:: https://img.shields.io/pypi/dm/spglib.svg?maxAge=2592000   
   :target: https://pypi.python.org/pypi/spglib
.. |Build Status| image:: https://travis-ci.org/atztogo/spglib.svg?branch=master
   :target: https://travis-ci.org/atztogo/spglib
Fortran and Ruby interfaces
============================

It is expected that Fortran and Ruby interfaces work at a certain
level but not as goog as Python interface.

Ruby interface
---------------

The ruby script ``symPoscar.rb`` in ``ruby`` directory reads a VASP
POSCAR formatted file and finds the space-group type. The ``-l``
option refines a slightly distorted structure within tolerance
spacified with the ``-s`` option. The ``-o`` option gives the symmetry
operations of the input structure.

::

   ruby symPoscar.rb POSCAR -s 0.1 -o -l


Some more options are shown with ``--help`` option::

   ruby symPoscar.rb --help

The way to compile the ruby module is explained in 
`ruby/README <https://github.com/atztogo/spglib/blob/master/ruby/README>`_.

Fortran interface
------------------

Fortran interface `spglib_f08.f90
<https://github.com/atztogo/spglib/blob/master/example/example_f08.f90>`_
is found in `example
<https://github.com/atztogo/spglib/tree/master/example>`_
directory. This fortran interface is the contribution from Dimitar
Pashov. `spglib_f.c
<https://github.com/atztogo/spglib/blob/master/src/spglib_f.c>`_ is
obsolete and is not recommended to use.
.. _python_spglib:

Spglib for Python
==================

.. contents::
   :depth: 2
   :local:

Installation
-------------

The source code is downloaded at
https://github.com/atztogo/spglib/releases . 
But minor updates are not included in this package. If you want the
latest version, you can git-clone the spglib repository::

   % git clone https://github.com/atztogo/spglib.git

It is also possible to install spglib for python via the following
package distribution services or from building using setup.py.

The following pip and conda packages are made and maintained by
`Paweł T. Jochym <https://github.com/jochym>`_, which is of great help
to keeping spglib handy and useful.

Using pip
^^^^^^^^^

Numpy is required before the python-spglib installation. The command to
install spglib is::

   % pip install --user spglib

If you see the error message like below in the installation process::

   _spglib.c:35:20: fatal error: Python.h: No such file or directory

development tools for building python module are additionally
necessary and are installed using OS's package management system,
e.g.,::

   % sudo apt-get install python-dev

If your installation by pip failed, you may need to
upgrade setuptools, e.g., by::

   % pip install --upgrade --user setuptools

Using conda
^^^^^^^^^^^

Conda is another choice::

   % conda install -c conda-forge spglib

Building using setup.py
^^^^^^^^^^^^^^^^^^^^^^^

To manually install python-spglib using ``setup.py``, python header
files (python-dev), C-compiler (e.g., gcc, clang), and numpy are
required before the build. The installation steps are shown as
follows:

1. Go to the :file:`python` directory
2. Type the command::

      % python setup.py install --user

   Document about where spglib is installed is found at the
   links below:
   
   - https://docs.python.org/2/install/#alternate-installation-the-user-scheme
   - https://docs.python.org/3/install/#alternate-installation-the-user-scheme

If your installation by setup.py failed, you may need to upgrade
setuptools, e.g., by::

   % pip install --upgrade --user setuptools


Test
-----

The test script :file:`test_spglib.py` is found in :file:`python/test`
directory. Got to this directory and run this script. It will be like below::

   % python test_spglib.py
   test_find_primitive (__main__.TestSpglib) ... ok
   test_get_symmetry (__main__.TestSpglib) ... ok
   test_get_symmetry_dataset (__main__.TestSpglib) ... ok
   test_refine_cell (__main__.TestSpglib) ... ok
   
   ----------------------------------------------------------------------
   Ran 4 tests in 13.147s
   
   OK

How to import spglib module
---------------------------

**Change in version 1.9.0!**

For versions 1.9.x or later::

   import spglib     

For versions 1.8.x or before::

   from pyspglib import spglib

If the version is not sure::

   try:
       import spglib as spg
   except ImportError:
       from pyspglib import spglib as spg   

Version number
--------------

In version 1.8.3 or later, the version number is obtained by
``spglib.__version__`` or :ref:`method_get_version`.

Example
--------

Examples are found in `examples
<https://github.com/atztogo/spglib/tree/master/python/examples>`_
directory.

.. _py_variables:

Variables
----------

.. _py_variables_crystal_structure:

Crystal structure (``cell``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A crystal structure is given by a **tuple**. This tuple format is
supported at version 1.9.1 or later. Optionally, an **ASE Atoms-like
object** is also supported. An alternative Atoms class (`atoms.py
<https://github.com/atztogo/spglib/blob/master/python/examples/atoms.py>`_)
that contains minimum set of methods is prepared in the `examples
<https://github.com/atztogo/spglib/tree/master/python/examples>`_
directory. When using ASE Atoms-like object, ``get_symmetry`` with
collinear polarizations is not supported.

The tuple format is shown as follows. There are three or four elements
in the tuple: ``cell = (lattice, positions, numbers)`` or ``cell =
(lattice, positions, numbers, magmoms)`` where ``magmoms`` represents
collinear polarizations on atoms and is optional.

Lattice parameters ``lattice`` are given by a 3x3 matrix with floating
point values, where :math:`\mathbf{a}, \mathbf{b}, \mathbf{c}` are
given as rows, which results in the transpose of the definition for
C-API (:ref:`variables_lattice`). Fractional atomic positions
``positions`` are given by a Nx3 matrix with floating point values,
where N is the number of atoms. Numbers to distinguish atomic species
``numbers`` are given by a list of N integers. The collinear polarizations
``magmoms`` only work with ``get_symmetry`` and are given
as a list of N floating point values.

::

   lattice = [[a_x, a_y, a_z],
              [b_x, b_y, b_z],
              [c_x, c_y, c_z]]
   positions = [[a_1, b_1, c_1],
                [a_2, b_2, c_2],
                [a_3, b_3, c_3],
                ...]
   numbers = [n_1, n_2, n_3, ...]
   magmoms = [m_1, m_2, m_3, ...]  # Only works with get_symmetry


**Version 1.9.5 or later**:
When a crystal structure is not properly given, the methods that use
the crsytal strcutre will return ``None``.

Symmetry tolerance (``symprec``)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Distance tolerance in Cartesian coordinates to find crystal symmetry.

Methods
--------

.. _method_get_version:

``get_version``
^^^^^^^^^^^^^^^^

**New in version 1.8.3**

::

    version = get_version()

This returns version number of spglib by tuple with three numbers.

``get_error_message``
^^^^^^^^^^^^^^^^^^^^^^

**New in version 1.9.5**

This method may be used to see why spglib failed though error handling
in spglib is not very sophisticated.

::

   error_message = get_error_message()

``get_spacegroup``
^^^^^^^^^^^^^^^^^^^

::

    spacegroup = get_spacegroup(cell, symprec=1e-5)

International space group short symbol and number are obtained as a
string. With ``symbol_type=1``, Schoenflies symbol is given instead of
international symbol.

.. _py_method_get_symmetry:

``get_symmetry``
^^^^^^^^^^^^^^^^^^

::

    symmetry = get_symmetry(cell, symprec=1e-5)

Symmetry operations are obtained as a dictionary. The key ``rotation``
contains a numpy array of integer, which is "number of symmetry
operations" x "3x3 matrices". The key ``translation`` contains a numpy
array of float, which is "number of symmetry operations" x
"vectors". The orders of the rotation matrices and the translation
vectors correspond with each other, e.g. , the second symmetry
operation is organized by the set of the second rotation matrix and second
translation vector in the respective arrays. Therefore a set of
symmetry operations may obtained by::

   [(r, t) for r, t in zip(dataset['rotations'], dataset['translations'])]

The operations are given with respect to the fractional coordinates
(not for Cartesian coordinates). The rotation matrix and translation
vector are used as follows::

    new_vector[3x1] = rotation[3x3] * vector[3x1] + translation[3x1]

The three values in the vector are given for the a, b, and c axes,
respectively. The key ``equivalent_atoms`` gives a mapping table of
atoms to symmetrically independent atoms. This is used to find
symmetrically equivalent atoms. The numbers contained are the indices
of atoms starting from 0, i.e., the first atom is numbered as 0, and
then 1, 2, 3, ... ``np.unique(equivalent_atoms)`` gives representative
symmetrically independent atoms. A list of atoms that are
symmetrically euivalent to some independent atom (here for example 1
is in ``equivalent_atom``) is found by
``np.where(equivalent_atom=1)[0]``. When the search failed, ``None``
is returned.

If ``cell`` is given as a tuple and collinear polarizations are given
as the fourth element of this tuple, symmetry operations are searched
considering this freedome. In ASE Atoms-class object, this is not supported.

``refine_cell``
^^^^^^^^^^^^^^^^

**Behaviour changed in version 1.8.x**

::

    lattice, scaled_positions, numbers = refine_cell(cell, symprec=1e-5)

Standardized crystal structure is obtained as a tuple of lattice (a 3x3
numpy array), atomic scaled positions (a numpy array of
[number_of_atoms,3]), and atomic numbers (a 1D numpy array) that are
symmetrized following space group type. When the search
failed, ``None`` is returned.

The detailed control of standardization of unit cell is achieved using
``standardize_cell``.

``find_primitive``
^^^^^^^^^^^^^^^^^^^

**Behaviour changed in version 1.8.x**

::

   lattice, scaled_positions, numbers = find_primitive(cell, symprec=1e-5)

When a primitive cell is found, lattice parameters (a 3x3 numpy array),
scaled positions (a numpy array of [number_of_atoms,3]), and atomic
numbers (a 1D numpy array) is returned. When the search failed,
``None`` is returned.

The detailed control of standardization of unit cell can be done using
``standardize_cell``.

``standardize_cell``
^^^^^^^^^^^^^^^^^^^^^

**New in version 1.8.x**

::

   lattice, scaled_positions, numbers = standardize_cell(bulk, to_primitive=False, no_idealize=False, symprec=1e-5)

``to_primitive=True`` is used to create the standardized primitive
cell, and ``no_idealize=True`` disables to idealize lengths and angles
of basis vectors and positions of atoms according to crystal
symmetry. Now ``refine_cell`` and ``find_primitive`` are shorthands of
this method with combinations of these options. When the search
failed, ``None`` is returned.  is returned. More detailed explanation
is shown in the spglib (C-API) document.

.. _py_method_get_symmetry_dataset:

``get_symmetry_dataset``
^^^^^^^^^^^^^^^^^^^^^^^^^^

**At version 1.9.4, the member 'choice' is added.**

::

    dataset = get_symmetry_dataset(cell, symprec=1e-5, angle_tolerance=-1.0, hall_number=0)

The arguments are:

* ``cell`` and ``symprec``: See :ref:`py_variables`.
* ``angle_tolerance``: An experimental argument that controls angle
  tolerance between basis vectors. Normally it is not recommended to use
  this argument. See a bit more detail at
  :ref:`variables_angle_tolerance`.
* ``hall_number`` (see the definition of this number at
  :ref:`api_spg_get_dataset_spacegroup_type`): The argument to
  constrain the space-group-type search only for the Hall symbol
  corresponding to it. The mapping from Hall symbols to a
  space-group-type is the many-to-one mapping. Without specifying this
  option (i.e., in the case of ``hall_number=0``), always the first one
  (the smallest serial number corresponding to the space-group-type in
  `list of space groups (Seto's web site)
  <http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en>`_)
  among possible choices and settings is chosen as default. This
  argument is useful when the other choice (or settting) is
  expected to be hooked. This affects to the obtained values of ``international``,
  ``hall``, ``hall_number``, ``choice``, ``transformation_matrix``,
  ``origin shift``, ``wyckoffs``, ``std_lattice``, ``std_positions``,
  and ``std_types``, but not to ``rotations`` and ``translations`` since
  the later set is defined with respect to the basis vectors of user's
  input (the ``cell`` argument).

``dataset`` is a dictionary. The keys are:

* ``number``: International space group number
* ``international``: International short symbol
* ``hall``: Hall symbol
* ``hall_number``: Hall number. This number is used in
  :ref:`py_method_get_symmetry_from_database` and
  :ref:`py_method_get_spacegroup_type`.
* ``choice``: Centring, origin, basis vector setting
* ``transformation_matrix``: See the detail at :ref:`api_origin_shift_and_transformation`.
* ``origin shift``: See the detail at :ref:`api_origin_shift_and_transformation`.
* ``wyckoffs``: Wyckoff letters
* ``equivalent_atoms``: Mapping table to equivalent atoms
* ``rotations`` and ``translations``: Rotation matrices and
  translation vectors. See :ref:`py_method_get_symmetry` for more
  details.
* ``pointgroup``: Symbol of the crystallographic point group in
  the Hermann–Mauguin notation.
* ``std_lattice``, ``std_positions``, ``std_types``: Standardized
  crystal structure corresponding to a Hall symbol found. These are
  equivalently given in the array formats of ``lattice``,
  ``positions``, and ``numbers`` presented at
  :ref:`py_variables_crystal_structure`, respectively.

..
   * ``pointgrouop_number``: Serial number of the crystallographic point
     group, which refers list of space groups (Seto’s web site)

When the search failed, ``None`` is returned. See more details of the
keys at :ref:`api_struct_spglibdataset`.

.. _py_method_get_symmetry_from_database:

``get_symmetry_from_database``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   symmetry = get_symmetry_from_database(hall_number)

A set of crystallographic symmetry operations corresponding to
``hall_number`` is returned by a dictionary where rotation parts and
translation parts are accessed by the keys ``rotations`` and
``translations``, respectively. The definition of ``hall_number`` is
found at :ref:`api_spg_get_dataset_spacegroup_type`.

When something wrong happened, ``None`` is returned.

.. _py_method_get_spacegroup_type:

``get_spacegroup_type``
^^^^^^^^^^^^^^^^^^^^^^^^

**New at version 1.9.4**

::

   spacegroup_type = get_spacegroup_type(hall_number)

This function allows to directly access to the space-group-type
database in spglib (spg_database.c). A dictionary is returned. To
specify the space group type with a specific choice, ``hall_number``
is used. The definition of ``hall_number`` is found at
:ref:`api_spg_get_dataset_spacegroup_type`. The keys of the returned
dictionary is as follows:

::

   number
   international_short
   international_full
   international
   schoenflies
   hall_symbol
   choice
   pointgroup_schoenflies
   pointgroup_international
   arithmetic_crystal_class_number
   arithmetic_crystal_class_symbol

Here ``spacegroup_type['international_short']`` is equivalent to
``dataset['international']`` of ``get_symmetry_dataset``,
``spacegroup_type['hall_symbol']`` is equivalent to
``dataset['hall']`` of ``get_symmetry_dataset``, and
``spacegroup_type['pointgroup_international']`` is equivalent to
``dataset['pointgroup_symbol']`` of ``get_symmetry_dataset``.

When something wrong happened, ``None`` is returned.

``niggli_reduce``
^^^^^^^^^^^^^^^^^^

**New at version 1.9.4**

::

   niggli_lattice = niggli_reduce(lattice, eps=1e-5)

Niggli reduction is achieved using this method. The algorithm detail
is found at https://atztogo.github.io/niggli/ and the references are
there in. Original basis vectors are stored in ``lattice`` and the
Niggli reduced basis vectors are given in ``niggli_lattice``. The
format of basis vectors are found at
:ref:`py_variables_crystal_structure`. ``esp`` is the tolerance
parameter, but unlike ``symprec`` the unit is not a length. This is
used to check if difference of norms of two basis vectors is close to
zero or not and if two basis vectors are orthogonal by the value of
dot product being close to zero or not.  The detail is shown at
https://atztogo.github.io/niggli/.

When the search failed, ``None`` is returned.

The transformation from original basis vectors :math:`( \mathbf{a}
\; \mathbf{b} \; \mathbf{c} )` to final baiss vectors :math:`(
\mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by :math:`\boldsymbol{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
elements have to be almost integers.

``delaunay_reduce``
^^^^^^^^^^^^^^^^^^^^

**New at version 1.9.4**

::

   delaunay_lattice = delaunay_reduce(lattice, eps=1e-5)

Delaunay reduction is achieved using this method. The algorithm is
found in the international tables for crystallography
volume A. Original basis vectors are stored in ``lattice`` and the
Delaunay reduced basis vectors are given in ``delaunay_lattice``,
where the format of basis vectors are shown in
:ref:`py_variables_crystal_structure`. ``esp`` is the tolerance
parameter, but unlike ``symprec`` the unit is not a length. This is
used as the criterion if volume is close to zero or not and if two
basis vectors are orthogonal by the value of dot product being close
to zero or not.

When the search failed, ``None`` is returned.

The transformation from original basis vectors :math:`( \mathbf{a}
\; \mathbf{b} \; \mathbf{c} )` to final basis vectors :math:`(
\mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by :math:`\boldsymbol{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
elements have to be almost integers.

``get_ir_reciprocal_mesh``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

::

   mapping, grid = get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

Irreducible k-points are obtained from a sampling mesh of k-points.
``mesh`` is given by three integers by array and specifies mesh
numbers along reciprocal primitive axis. ``is_shift`` is given by the
three integers by array. When ``is_shift`` is set for each reciprocal
primitive axis, the mesh is shifted along the axis in half of adjacent
mesh points irrespective of the mesh numbers. When the value is not 0,
``is_shift`` is set.

``mapping`` and ``grid`` are returned. ``grid`` gives the mesh points in
fractional coordinates in reciprocal space. ``mapping`` gives mapping to
the irreducible k-point indices that are obtained by ::

   np.unique(mapping)

Here ``np`` means the numpy module. The grid point is accessed by
``grid[index]``.

When the sesarch failed, ``None`` is returned.

An example is shown below::

   import numpy as np
   import spglib
   
   lattice = np.array([[0.0, 0.5, 0.5],
                       [0.5, 0.0, 0.5],
                       [0.5, 0.5, 0.0]]) * 5.4
   positions = [[0.875, 0.875, 0.875],
                [0.125, 0.125, 0.125]]
   numbers= [1,] * 2
   cell = (lattice, positions, numbers)
   print(spglib.get_spacegroup(cell, symprec=1e-5))
   mesh = [8, 8, 8]
   
   #
   # Gamma centre mesh
   #
   mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
   
   # All k-points and mapping to ir-grid points
   for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
       print("%3d ->%3d %s" % (i, ir_gp_id, gp.astype(float) / mesh))
   
   # Irreducible k-points
   print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
   print(grid[np.unique(mapping)] / np.array(mesh, dtype=float))
   
   #
   # With shift
   #
   mapping, grid = spglib.get_ir_reciprocal_mesh(mesh, cell, is_shift=[1, 1, 1])
   
   # All k-points and mapping to ir-grid points
   for i, (ir_gp_id, gp) in enumerate(zip(mapping, grid)):
       print("%3d ->%3d %s" % (i, ir_gp_id, (gp + [0.5, 0.5, 0.5]) / mesh))
   
   # Irreducible k-points
   print("Number of ir-kpoints: %d" % len(np.unique(mapping)))
   print((grid[np.unique(mapping)] + [0.5, 0.5, 0.5]) / mesh)
Spglib
======

Spglib is a library for finding and handling crystal symmetries
written in C.

Spglib for python is found at :ref:`python_spglib`.

Version history is summarized in `ChangeLog
<https://github.com/atztogo/spglib/blob/master/ChangeLog>`_.

News
-----

* SpglibDataset is slightly modified at version 1.9.4. This may hit
  some users.

* Python module pyspglib is renamed to spglib in version 1.9.0.

  In versions 1.8.x or before, the python wrapper (pyspglib) was
  considered as a special extention of spglib and separately maintained,
  but in version 1.9.x and later, it starts to be a part of spglib as a
  usual extension, and this python wrapper will be maintained to work
  following every update of c-spglib. The package name is changed from
  pyspglib to spglib. Following this change, the way to import python
  spglib module is changed and it is written in :ref:`python_spglib`.

Features
----------

* Find symmetry operations
* Identify space-group type
* Wyckoff position assignment
* Refine crystal structure
* Find a primitive cell
* Search irreducible k-points

Documentation
--------------

.. toctree::
   :maxdepth: 1

   install
   api
   variable
   Python interface <python-spglib>
   interface
   definition

References
----------

General materials

* International tables for crystallography volumes A, B, C
* Online dictionary of crystallography:
  http://reference.iucr.org/dictionary/
* Symmetry relationships between crystal structures by Ulrich Muller,
  Oxford university press (2013)

Specific algorithms

* Space-group type determination: R. W. Grosse-Kunstleve, Acta Cryst.,
  A55, 383-395 (1999)
* Crystal refinement: R. W. Grosse-Kunstleve and P. D. Adams, Acta
  Cryst., A58, 60-65 (2002)

Mailing list
-------------

For questions, bug reports, and comments, please visit following
mailing list:

https://lists.sourceforge.net/lists/listinfo/spglib-users

For more information
---------------------

* Repository: https://github.com/atztogo/spglib 
* License: New BSD after version 1.0.beta-1. The older versions are under GPLv2.
* Contact: atz.togo@gmail.com
* Authour: Atsushi Togo

Acknowledgments
----------------

Spglib project acknowledges Yusuke Seto for the Crystallographic
database, Dimitar Pashov for the fortran interface, Paweł T. Jochym
for Python packages.
How to install spglib C-API
============================

Download
---------

The source code is downloaded at
https://github.com/atztogo/spglib/releases . 
But minor updates are not included in this package. If you want the
latest version, you can git-clone the spglib repository::

   % git clone https://github.com/atztogo/spglib.git

Install
--------

Compiling using cmake
^^^^^^^^^^^^^^^^^^^^^^

1. After expanding source code, go into the source code directory::

     % tar xvfz spglib-1.9.8.tar.gz
     % cd spglib-1.9.8

   The current directory is ``PROJECT_SOURCE_DIR``.

2. Create a directory for the build::

     % mkdir _build; cd _build
     % cmake ..
     % make
     % make install

3. The libraries are installed at ``PROJECT_SOURCE_DIR/lib`` and the
   header file is installed at ``PROJECT_SOURCE_DIR/include``.

Compiling using configure script
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. The configure script is prepared using
   autotools and libtool as follows::

     % aclocal
     % autoheader
     % libtoolize # or glibtoolize with macport etc
     % touch INSTALL NEWS README AUTHORS
     % automake -acf
     % autoconf


2. Run configure script::

     % tar xvfz spglib-1.9.8.tar.gz
     % cd spglib-1.9.8
     % ./configure --prefix=INSTALLATION_LOCATION
     % make
     % make install

3. The libraries are installed at ``INSTALLATION_LOCATION/lib`` or found in
   ``src/.libs`` if you don't run ``make install``.

OpenMP
^^^^^^^

Bottle neck of symmetry operation search may be eased using the OpenMP
threading. In the case of gcc (> 4.2), set the following environment
variables before running configure script::

   % export LIBS='-lgomp'
   % export CFLAGS='-fopenmp'

Overhead of threading is relatively large when the number of atoms is
small. Therefore the threading is activated when the number of atoms
>= 1000.

Usage
------

1. Include ``spglib.h``
2. Link ``libsymspg.a`` or ``libsymspg.so``
3. A compilation example is shown in  `example/README
   <https://github.com/atztogo/spglib/blob/master/example/README>`_.

Example
--------

A few examples are found in `example
<https://github.com/atztogo/spglib/tree/master/example>`_ directory.
C-APIs
======

.. contents:: List of C-APIs
   :depth: 3
   :local:

``spg_get_major_version``, ``spg_get_minor_version``, ``spg_get_micro_version``
--------------------------------------------------------------------------------

**New in version 1.8.3**

Version number of spglib is obtained. These three functions return
integers that correspond to spglib version [major].[minor].[micro].

|

``spg_get_error_code`` and ``spg_get_error_message``
-----------------------------------------------------

**New in version 1.9.5**

These methods may be used to see why spglib failed though error handling
in spglib is not very sophisticated.

::

   SpglibError spg_get_error_code(void);

::

   char * spg_get_error_message(SpglibError spglib_error);

The ``SpglibError`` type is a enum type as shown below.

::

   typedef enum {
     SPGLIB_SUCCESS = 0,
     SPGERR_SPACEGROUP_SEARCH_FAILED,
     SPGERR_CELL_STANDARDIZATION_FAILED,
     SPGERR_SYMMETRY_OPERATION_SEARCH_FAILED,
     SPGERR_ATOMS_TOO_CLOSE,
     SPGERR_POINTGROUP_NOT_FOUND,
     SPGERR_NIGGLI_FAILED,
     SPGERR_DELAUNAY_FAILED,
     SPGERR_ARRAY_SIZE_SHORTAGE,
     SPGERR_NONE,
   } SpglibError;

The usage is as follows::

   SpglibError error;
   error = spg_get_error_code();
   printf("%s\n", spg_get_error_message(error));
   
|

.. _api_spg_get_symmetry:

``spg_get_symmetry``
---------------------

This function finds a set of representative symmetry operations for
primitive cells or its extension with lattice translations for
supercells. 0 is returned if it failed.

::

  int spg_get_symmetry(int rotation[][3][3],
  		       double translation[][3],
  		       const int max_size,
		       const double lattice[3][3],
  		       const double position[][3],
		       const int types[],
  		       const int num_atom,
		       const double symprec);

The operations are stored in ``rotation`` and ``translation``. The
number of operations is return as the return value. Rotations and
translations are given in fractional coordinates, and ``rotation[i]``
and ``translation[i]`` with same index give a symmetry operations,
i.e., these have to be used together.

As an exceptional case, if a supercell has the basis vectors of the
lattice that break crsytallographic point group, the crystallographic
symmetry operations are searched with this broken symmetry, i.e., at
most the crystallographic point group found in this case is the point
group of the lattice. For example, this happens for the :math:`2\times
1\times 1` supercell of a conventional cubic unit cell. This may not
be understandable in crystallographic sense, but is practically useful
treatment for research in computational materials science.

|

``spg_get_international``
--------------------------

Space group type is found and returned in international table symbol
to ``symbol`` and also as a number (return value). 0 is returned if
it failed.

::

  int spg_get_international(char symbol[11],
                            const double lattice[3][3],
                            const double position[][3],
                            const int types[],
			    const int num_atom,
                            const double symprec);


|

``spg_get_schoenflies``
-------------------------

Space group type is found and returned in schoenflies to ``symbol``
and also as a number (return value). 0 is returned if it failed.

::

  int spg_get_schoenflies(char symbol[10],
                          const double lattice[3][3],
                          const double position[][3],
                          const int types[],
                          const int num_atom,
                          const double symprec);



|

``spg_standardize_cell``
-------------------------

The standardized unit cell (see :ref:`def_standardized_unit_cell`) is
generated from an input unit cell structure and its space group type
determined about a symmetry search tolerance. Usually
``to_primitive=0`` and ``no_idealize=0`` are recommended to set and
this setting results in the same behavior as ``spg_refine_cell``. 0 is
returned if it failed.

::

   int spg_standardize_cell(double lattice[3][3],
                            double position[][3],
                            int types[],
                            const int num_atom,
                            const int to_primitive,
                            const int no_idealize,
                            const double symprec);

Number of atoms in the found standardized unit (primitive) cell is
returned.

``to_primitive=1`` is used to create the standardized primitive cell
with the transformation matricies shown at
:ref:`def_standardized_primitive_cell`, otherwise ``to_primitive=0``
must be specified. The found basis vectors and
atomic point coordinates and types are overwritten in ``lattice``,
``position``, and ``types``, respectively. Therefore with
``to_primitive=0``, at a maximum four times larger array size for
``position`` and ``types`` than the those size of the input unit cell
is required to store a standardized unit cell with face centring found
in the case that the input unit cell is a primitive cell.

``no_idealize=1`` disables to idealize lengths and angles of basis
vectors and positions of atoms according to crystal symmetry. The
detail of the idealization (``no_idealize=0``) is written at
:ref:`def_idealize_cell`. ``no_idealize=1`` may be used when we want to
leave basis vectors and atomic positions in Cartesianl coordinates
fixed.

|

``spg_find_primitive``
-----------------------

**Behavior is changed. This function is now a shortcut of**
``spg_standardize_cell`` **with**
``to_primitive=1`` **and** ``no_idealize=0``.

A primitive cell is found from an input unit cell. 0 is returned if it
failed.

::

  int spg_find_primitive(double lattice[3][3],
                         double position[][3],
                         int types[],
			 const int num_atom,
			 const double symprec);

``lattice``, ``position``, and ``types`` are overwritten. Number of
atoms in the found primitive cell is returned.

|

``spg_refine_cell``
--------------------

**This function exists for backward compatibility since it is same as** ``spg_standardize_cell`` **with** ``to_primitive=0`` **and** ``leave_distorted=0``.

The standardized crystal structure is obtained from a non-standard
crystal structure which may be slightly distorted within a symmetry
recognition tolerance, or whose primitive vectors are differently
chosen, etc. 0 is returned if it failed.

::

  int spg_refine_cell(double lattice[3][3],
		      double position[][3],
		      int types[],
		      const int num_atom,
 		      const double symprec);

The calculated standardized lattice and atomic positions overwrites
``lattice``, ``position``, and ``types``. The number of atoms in the
standardized unit cell is returned as the return value. When the input
unit cell is a primitive cell and is the face centring symmetry, the
number of the atoms returned becomes four times large. Since this
function does not have any means of checking the array size (memory
space) of these variables, the array size (memory space) for
``position`` and ``types`` should be prepared **four times more** than
those required for the input unit cell in general.

|

.. _api_spg_get_dataset:

``spg_get_dataset`` and ``spg_get_dataset_with_hall_number``
--------------------------------------------------------------

**Changed in version 1.8.1**

For an input unit cell structure, symmetry operations of the crystal
are searched. Then they are compared with the crsytallographic
database and the space group type is determined. The result is
returned as the ``SpglibDataset`` structure as a dataset. The default
choice of setting of basis vectors in spglib is explained in the
manuscript found at http://arxiv.org/abs/1506.01455.

Usage
^^^^^^

Dataset corresponding to the space group type in the standard setting
is obtained by ``spg_get_dataset``. If this symmetry search fails,
``NULL`` is returned in version 1.8.1 or later (spacegroup_number = 0
is returned in the previous versions). In this function, the other
crystallographic setting is not obtained.

::

   SpglibDataset * spg_get_dataset(const double lattice[3][3],
                                   const double position[][3],
                                   const int types[],
                                   const int num_atom,
                                   const double symprec);

To specify the other crystallographic choice (setting, origin, axis,
or cell choice), ``spg_get_dataset_with_hall_number`` is used.

::

   SpglibDataset * spg_get_dataset_with_hall_number(SPGCONST double lattice[3][3],
						    SPGCONST double position[][3],
						    const int types[],
						    const int num_atom,
						    const int hall_number,
						    const double symprec)

where ``hall_number`` is used to specify the choice. The possible
choices and those serial numbers are found at `list of space groups
(Seto's web site)
<http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en>`_.
The crystal structure has to possess the space-group type of the Hall
symbol. If the symmetry search fails or the specified ``hall_number``
is not in the list of Hall symbols for the space group type of the
crystal structure, ``spacegroup_number`` in the ``SpglibDataset``
structure is set 0.

Finally, its allocated memory space must be freed by calling
``spg_free_dataset``.

.. _api_struct_spglibdataset:

Dataset
^^^^^^^^

**At version 1.9.4, SpglibDataset was modified.** The member name
``setting`` is changed to ``choice`` and ``pointgroup_number`` is removed.

The dataset is accessible through the C-structure given by

::

   typedef struct {
     int spacegroup_number;
     int hall_number;
     char international_symbol[11];
     char hall_symbol[17];
     char choice[6];
     double transformation_matrix[3][3];
     double origin_shift[3];
     int n_operations;
     int (*rotations)[3][3];
     double (*translations)[3];
     int n_atoms;
     int *wyckoffs;
     int *equivalent_atoms;
     int n_std_atoms;             /* n_brv_atoms before version 1.8.1 */
     double std_lattice[3][3];    /* brv_lattice before version 1.8.1 */
     int *std_types;              /* brv_types before version 1.8.1 */
     double (*std_positions)[3];  /* brv_positions before version 1.8.1 */
     char pointgroup_symbol[6];
   } SpglibDataset;

In **versions before 1.8.1**, the member names of ``n_std_atoms``,
``std_lattice``, ``std_types``, and ``std_positions`` were
``n_brv_atoms``, ``brv_lattice``, ``brv_types``, and
``brv_positions``, respectively.

.. _api_spg_get_dataset_spacegroup_type:

Space group type
"""""""""""""""""

``spacegroup_number`` is the space group type number defined in
International Tables for Crystallography (ITA). ``hall_number`` is the
serial number between 1 and 530 which are found at `list of space
groups (Seto's web site)
<http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en>`_.
The (full) Hermann–Mauguin notation of space group type is given by
``international_symbol``. The Hall symbol is stored in
``hall_symbol``. The information on unique axis,
setting or cell choices is found in ``choice``.

Symmetry operations
"""""""""""""""""""""""

The symmetry operations of the input unit cell are stored in
``rotations`` and ``translations``. A crystallographic symmetry
operation :math:`(\boldsymbol{W}, \boldsymbol{w})` is made from a pair
of rotation :math:`\boldsymbol{W}` and translation
:math:`\boldsymbol{w}` parts with the same index. Number of symmetry
operations is given as ``n_operations``. The detailed explanation of
the values is found at :ref:`api_spg_get_symmetry`.

.. _api_spg_get_dataset_site_symmetry:

Site symmetry
""""""""""""""

``n_atoms`` is the number of atoms of the input unit
cell. ``wyckoffs`` gives Wyckoff letters that are assigned to atomic
positions of the input unit cell. The numbers of 0, 1, 2,
:math:`\ldots`, correspond to the a, b, c, :math:`\ldots`,
respectively. Number of elements in ``wyckoffs`` is same as
``n_atoms``. ``equivalent_atoms`` is a list of atomic indices that map
to indices of symmetrically independent atoms, where the list index
corresponds to atomic index of the input crystal structure.

.. _api_origin_shift_and_transformation:

Origin shift and lattice transformation
""""""""""""""""""""""""""""""""""""""""

**Changed in version 1.8.1**

``transformation_matrix`` and ``origin_shift`` are obtained as a
result of space-group-type matching under a set of unique axis,
setting and cell choices. In this matching, basis vectors and atomic
point coordinates have to be standardized to compare with the database
of symmetry operations. The basis vectors are transformed to those of
a standardized unit cell. Atomic point coordinates are shifted so that
symmetry operations have the standard
origin. ``transformation_matrix`` (:math:`\boldsymbol{P}`) is the
matrix to transform the input basis vectors to the standardized basis
vectors, wihch is represented as

.. math::

   ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
   = ( \mathbf{a}_\mathrm{s} \; \mathbf{b}_\mathrm{s} \; \mathbf{c}_\mathrm{s} )  \boldsymbol{P}

where :math:`\mathbf{a}`, :math:`\mathbf{b}`, and :math:`\mathbf{c}`
are the input (original) basis vectors, and
:math:`\mathbf{a}_\mathrm{s}`, :math:`\mathbf{b}_\mathrm{s}`, and
:math:`\mathbf{c}_\mathrm{s}` are the standardized basis vectors. The
``origin_shift`` (:math:`\boldsymbol{p}`) is the vector from the
origin of the standardized coordinate system to the origin of the
input (original) coordinate system measured in the standardized
coordinate system. The atomic point shift is measured from the
standardized unit cell (conventional unit cell) to the original unit
cell measured in the coordinates of the standardized unit cell. An
atomic point in the original unit cell :math:`\boldsymbol{x}` (input
data) is mapped to that in the standardized unit cell
:math:`\boldsymbol{x}_\mathrm{s}` by

.. math::

   \boldsymbol{x}_\mathrm{s} = \boldsymbol{P}\boldsymbol{x} +
   \boldsymbol{p} \;\;(\mathrm{mod}\; \mathbf{1}).

In **versions 1.7.x and 1.8 or before**, ``transformation_matrix`` and
``origin_shift`` are defined as follows:

.. math::

   ( \mathbf{a}_\mathrm{s} \; \mathbf{b}_\mathrm{s} \;
   \mathbf{c}_\mathrm{s} ) = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c}
   ) \boldsymbol{P} \;\; \text{and} \;\; \boldsymbol{x}_\mathrm{s} =
   \boldsymbol{P}^{-1}\boldsymbol{x} - \boldsymbol{p}
   \;\;(\mathrm{mod}\; \mathbf{1}),

respectively.

Standardized crystal structure
"""""""""""""""""""""""""""""""

**Changed in version 1.8.1**

The standardized crystal structure corresponding to a Hall symbol is
stored in ``n_std_atoms``, ``std_lattice``, ``std_types``, and
``std_positions``.

In **versions 1.7.x and 1.8 or before**, the variable names of the
members corresponding to those above are ``n_brv_atoms``,
``brv_lattice``, ``brv_types``, and ``brv_positions``, respectively.

Crystallographic point group
"""""""""""""""""""""""""""""

**New in version 1.8.1**

``pointgroup_number`` is the serial number of the crystallographic
point group, which refers `list of space
groups (Seto's web site)
<http://pmsl.planet.sci.kobe-u.ac.jp/~seto/?page_id=37&lang=en>`_.
``pointgroup_symbol`` is the symbol of the crystallographic point
group in the Hermann–Mauguin notation.

|

``spg_free_dataset``
---------------------

Allocated memoery space of the C-structure of ``SpglibDataset`` is
freed by calling ``spg_free_dataset``.

::

  void spg_free_dataset(SpglibDataset *dataset);

|

``spg_get_spacegroup_type``
-----------------------------

**Changed at version 1.9.4: Some members are added and the member name 'setting' is changed to 'choice'.**

This function allows to directly access to the space-group-type
database in spglib (spg_database.c). To specify the space group type
with a specific choice, ``hall_number`` is used. The definition of
``hall_number`` is found at
:ref:`api_spg_get_dataset_spacegroup_type`.
``number = 0`` is returned when it failed.

::

   SpglibSpacegroupType spg_get_spacegroup_type(const int hall_number)

``SpglibSpacegroupType`` structure is as follows:

::

   typedef struct {
     int number;
     char international_short[11];
     char international_full[20];
     char international[32];
     char schoenflies[7];
     char hall_symbol[17];
     char choice[6];
     char pointgroup_schoenflies[4];
     char pointgroup_international[6];
     int arithmetic_crystal_class_number;
     char arithmetic_crystal_class_symbol[7];
   } SpglibSpacegroupType;

|

``spg_get_symmetry_from_database``
-----------------------------------

This function allows to directly access to the space group operations
in the spglib database (spg_database.c). To specify the space group
type with a specific choice, ``hall_number`` is used. The definition
of ``hall_number`` is found at
:ref:`api_spg_get_dataset_spacegroup_type`. 0 is returned when it
failed.

::

   int spg_get_symmetry_from_database(int rotations[192][3][3],
				      double translations[192][3],
				      const int hall_number);

The returned value is the number of space group operations. The space
group operations are stored in ``rotations`` and ``translations``.

|

``spg_get_multiplicity``
-------------------------

This function returns exact number of symmetry operations. 0 is
returned when it failed.

::

  int spg_get_multiplicity(const double lattice[3][3],
  			   const double position[][3],
  			   const int types[],
			   const int num_atom,
  			   const double symprec);

This function may be used in advance to allocate memoery space for
symmetry operations.

|

``spg_get_symmetry_with_collinear_spin``
-----------------------------------------

This function finds symmetry operations with collinear polarizations
(spins) on atoms. Except for the argument of ``const double spins[]``,
the usage is basically the same as ``spg_get_symmetry``, but as an
output, ``equivalent_atoms`` are obtained. The size of this array is
the same of ``num_atom``. See :ref:`api_spg_get_dataset_site_symmetry`
for the definition ``equivalent_atoms``. 0 is returned when it failed.

::

  int spg_get_symmetry_with_collinear_spin(int rotation[][3][3],
                                           double translation[][3],
					   int equivalent_atoms[],
                                           const int max_size,
                                           SPGCONST double lattice[3][3],
                                           SPGCONST double position[][3],
                                           const int types[],
                                           const double spins[],
                                           const int num_atom,
                                           const double symprec);


|

``spg_niggli_reduce``
----------------------

Niggli reduction is applied to input basis vectors ``lattice`` and the
reduced basis vectors are overwritten to ``lattice``. 0 is returned if
it failed.

::

   int spg_niggli_reduce(double lattice[3][3], const double symprec);

The transformation from original basis vectors :math:`( \mathbf{a}
\; \mathbf{b} \; \mathbf{c} )` to final basis vectors :math:`(
\mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by :math:`\boldsymbol{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
elements have to be almost integers.

|

``spg_delaunay_reduce``
------------------------

Delaunay reduction is applied to input basis vectors ``lattice`` and
the reduced basis vectors are overwritten to ``lattice``. 0 is
returned if it failed.

::

   int spg_delaunay_reduce(double lattice[3][3], const double symprec);

The transformation from original basis vectors :math:`( \mathbf{a}
\; \mathbf{b} \; \mathbf{c} )` to final basis vectors :math:`(
\mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )` is achieved by linear
combination of basis vectors with integer coefficients without
rotating coordinates. Therefore the transformation matrix is obtained
by :math:`\boldsymbol{P} = ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
( \mathbf{a}' \; \mathbf{b}' \; \mathbf{c}' )^{-1}` and the matrix
elements have to be almost integers.


``spg_get_ir_reciprocal_mesh``
-------------------------------

Irreducible reciprocal grid points are searched from uniform mesh grid
points specified by ``mesh`` and ``is_shift``.

::

   int spg_get_ir_reciprocal_mesh(int grid_address[][3],
                                  int map[],
                                  const int mesh[3],
                                  const int is_shift[3],
                                  const int is_time_reversal,
                                  const double lattice[3][3],
                                  const double position[][3],
                                  const int types[],
                                  const int num_atom,
                                  const double symprec)

``mesh`` stores three integers. Reciprocal primitive vectors are
divided by the number stored in ``mesh`` with (0,0,0) point
centering. The center of grid mesh is shifted +1/2 of a grid spacing
along corresponding reciprocal axis by setting 1 to a ``is_shift``
element. No grid mesh shift is made if 0 is set for ``is_shift``.

The reducible uniform grid points are returned in fractional coordinates
as ``grid_address``. A map between reducible and irreducible points are
returned as ``map`` as in the indices of ``grid_address``. The number of
the irreducible k-points are returned as the return value.  The time
reversal symmetry is imposed by setting ``is_time_reversal`` 1.

Grid points are stored in the order that runs left most element
first, e.g. (4x4x4 mesh).::

   [[ 0  0  0]
    [ 1  0  0]
    [ 2  0  0]
    [-1  0  0]
    [ 0  1  0]
    [ 1  1  0]
    [ 2  1  0]
    [-1  1  0]
    ....      ]

where the first index runs first.  k-qpoints are calculated by
``(grid_address + is_shift / 2) / mesh``. A grid point index is
recovered from ``grid_address`` by ``numpy.dot(grid_address % mesh,
[1, mesh[0], mesh[0] * mesh[1]])`` in Python-numpy notation, where
``%`` always returns non-negative integers. The order of
``grid_address`` can be changed so that the last index runs first by
setting the macro ``GRID_ORDER_XYZ`` in ``kpoint.c``. In this case the
grid point index is recovered by ``numpy.dot(grid_address % mesh,
[mesh[2] * mesh[1], mesh[2], 1])``.

|

``spg_get_stabilized_reciprocal_mesh``
---------------------------------------

The irreducible k-points are searched from unique k-point mesh grids
from direct (real space) basis vectors and a set of rotation parts of
symmetry operations in direct space with one or multiple
stabilizers.

::

   int spg_get_stabilized_reciprocal_mesh(int grid_address[][3],
                                          int map[],
                                          const int mesh[3],
                                          const int is_shift[3],
                                          const int is_time_reversal,
                                          const int num_rot,
                                          const int rotations[][3][3],
                                          const int num_q,
                                          const double qpoints[][3])

The stabilizers are written in fractional coordinates. Number of the
stabilizers are given by ``num_q``. Symmetrically equivalent k-points
(stars) in fractional coordinates are stored in ``map`` as indices of
``grid_address``. The number of reduced k-points with the stabilizers
are returned as the return value.

This function can be used to obtain all mesh grid points by setting
``num_rot = 1``, ``rotations = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}``,
``num_q = 1``, and ``qpoints = {0, 0, 0}``.
Variables
==========

Variables frequently used in C-API are explained. In interfaces of
other languages, the definition may be different, e.g., ``lattice``
for Python interface is transposed.

.. contents::
   :depth: 2
   :local:

.. _variables_lattice:

``lattice``
-------------

**Attention**: In Python interface, ``lattice`` is transposed
(:ref:`py_variables_crystal_structure`.)

Basis vectors (in Cartesian)

::

  [ [ a_x, b_x, c_x ],
    [ a_y, b_y, c_y ],
    [ a_z, b_z, c_z ] ]

Cartesian position :math:`\mathbf{x}_\mathrm{Cart}` is obtained by

.. math::

  \mathbf{x}_\mathrm{Cart} = \mathrm{L}\cdot\mathbf{x}_\mathrm{frac}

where :math:`\mathrm{L}` is the basis vectors defined above and is
:math:`\mathrm{L}=(\mathbf{a},\mathbf{b},\mathbf{c})`, and
:math:`\mathbf{x}_\mathrm{frac}` is the atomic position in fractional
coordinates given below.


``position``
--------------

Atomic points (in fractional coordinates with respect to basis vectors)

::

  [ [ x1_a, x1_b, x1_c ], 
    [ x2_a, x2_b, x2_c ], 
    [ x3_a, x3_b, x3_c ], 
    ...                   ]


``types``
----------

Atomic species are differenciated by integers. The nubmers are just
used to distinguishes atoms, and is not releated to atomic numbers.

::

  [ type_1, type_2, type_3, ... ]

The number of elements is same as the number of atoms.

``rotation``
--------------

Rotation matricies (ingeger) of symmetry operations.

::

    [ [ r_aa, r_ab, r_ac ],
      [ r_ba, r_bb, r_bc ],
      [ r_ca, r_cb, r_cc ] ]

``translation``
-----------------

Translation vectors corresponding to symmetry operations in fractional
coordinates.

::

    [ t_a, t_b, t_c ]

``symprec``
------------

Tolerance of distance between atomic positions and between lengths of
lattice vectors to be tolerated in the symmetry finding. The angle
distortion between lattice vectors is converted to a length and
compared with this distance tolerance. If the explicit angle tolerance
is expected, see ``angle_tolerance``.

.. _variables_angle_tolerance:

``angle_tolerance``
--------------------

**Experimental**

Tolerance of angle between lattice vectors in degrees to be tolerated
in the symmetry finding. To use angle tolerance, another set of
functions are prepared as follows

::

   spgat_get_dataset
   spgat_get_symmetry
   spgat_get_symmetry_with_collinear_spin
   spgat_get_multiplicity
   spgat_find_primitive
   spgat_get_international
   spgat_get_schoenflies
   spgat_refine_cell

These functions are called by the same way with an additional argument
of ``const double angle_tolerance`` in degrees. By specifying a negative
value, the behavior becomes the same as usual functions. The default
value of ``angle_tolerance`` is a negative value.
Definitions and conventions
============================

Information in this page is valid for spglib 1.8.1 or later. The
definitions of transformation matrix and origin shift were different
in the previous versions.

.. contents::
   :depth: 2
   :local:


References
-----------

Some references about crystallographic definitions and conventions are
shown below. Though spglib may not follow them fully, it doesn't mean
spglib doesn't respect them, rather it is due to the lack of
understanding by the author of spglib.

* `International Tables for Crystallography <http://it.iucr.org/>`_.
* `Bilbao Crystallographic Server <http://www.cryst.ehu.es/>`_. The
  references of many useful papers are found at
  http://www.cryst.ehu.es/wiki/index.php/Articles.
* Ulrich Müller, "Symmetry Relationships between Crystal Structures"
* E. Parthé, K. Cenzual, and R. E. Gladyshevskii, "Standardization of
  crystal structure data as an aid to the classification of crystal
  structure types", Journal of Alloys and Compounds, **197**, 291-301
  (1993). [`doi2
  <https://dx.doi.org/10.1016/0925-8388(93)90049-S>`_]
* E. Parthé and L. M. Gelato, "The ’best’ unit cell for monoclinic
  structures consistent with b axis unique and cell choice 1
  of international tables for crystallography (1983)", Acta
  Cryst. A **41**, 142-151 (1985) [`doi3
  <https://doi.org/10.1107/S0108767385000289>`_]
* E. Parthé and L. M. Gelato, "The standardization of inorganic
  crystal-structure data", Acta Cryst. A
  **40**, 169-183 (1984) [`doi4
  <https://doi.org/10.1107/S0108767384000416>`_]
* S. Hall, "Space-group notation with an explicit origin", Acta
  Cryst. A **37**, 517-525 (1981) [`doi1
  <https://doi.org/10.1107/S0567739481001228>`_]

Basis vectors :math:`(\mathbf{a}, \mathbf{b}, \mathbf{c})` or :math:`(\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3)`
------------------------------------------------------------------------------------------------------------------

In spglib, basis vectors are represented by three column vectors:

.. math::

   \mathbf{a}= \begin{pmatrix}
   a_x \\
   a_y \\
   a_z \\
   \end{pmatrix},
   \mathbf{b}= \begin{pmatrix}
   b_x \\
   b_y \\
   b_z \\
   \end{pmatrix},
   \mathbf{c}= \begin{pmatrix}
   c_x \\
   c_y \\
   c_z \\
   \end{pmatrix},

in Cartesian coordinates. Depending on the situation,
:math:`(\mathbf{a}_1, \mathbf{a}_2, \mathbf{a}_3)` is used instead of
:math:`(\mathbf{a}, \mathbf{b}, \mathbf{c})`.

Atomic point coordinates :math:`\boldsymbol{x}`
-----------------------------------------------

Coordinates of an atomic point :math:`\boldsymbol{x}` are represented
as three fractional values relative to basis vectors as follows,

.. math::

   \boldsymbol{x}= \begin{pmatrix}
   x_1 \\
   x_2 \\
   x_3 \\
   \end{pmatrix},

where :math:`0 \le x_i < 1`. A position vector :math:`\mathbf{x}` in
Cartesian coordinates is obtained by

.. math::

   \mathbf{x} = (\mathbf{a}, \mathbf{b}, \mathbf{c}) \boldsymbol{x}.

or 

.. math::

   \mathbf{x} = \sum_i x_i \mathbf{a}_i.

Symmetry operation :math:`(\boldsymbol{W}, \boldsymbol{w})`
-----------------------------------------------------------

A symmetry operation consists of a pair of the rotation part
:math:`\boldsymbol{W}` and translation part :math:`\boldsymbol{w}`,
and is represented as :math:`(\boldsymbol{W}, \boldsymbol{w})` in the
spglib document. The symmetry operation transfers :math:`\boldsymbol{x}` to
:math:`\tilde{\boldsymbol{x}}` as follows:

.. math::

  \tilde{\boldsymbol{x}} = \boldsymbol{W}\boldsymbol{x} + \boldsymbol{w}.

Transformation matrix :math:`\boldsymbol{P}` and origin shift :math:`\boldsymbol{p}`
-------------------------------------------------------------------------------------

The transformation matrix :math:`\boldsymbol{P}` changes choice of
basis vectors as follows

.. math::

   ( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )
   = ( \mathbf{a}_\mathrm{s} \; \mathbf{b}_\mathrm{s} \;
   \mathbf{c}_\mathrm{s} )  \boldsymbol{P},

where :math:`( \mathbf{a} \; \mathbf{b} \; \mathbf{c} )` and :math:`(
\mathbf{a}_\mathrm{s} \; \mathbf{b}_\mathrm{s} \;
\mathbf{c}_\mathrm{s} )` are the basis vectors of an arbitrary system
and of a starndardized system, respectively. Transformation matrix
**doesn't** rotate a crystal in Cartesian coordinates, but just
changes the choices of basis vectors.

The origin shift :math:`\boldsymbol{p}` gives the vector from the
origin of the standardized system :math:`\boldsymbol{O}_\mathrm{s}` to
the origin of the arbitrary system :math:`\boldsymbol{O}`,

.. math::

   \boldsymbol{p} = \boldsymbol{O} - \boldsymbol{O}_\mathrm{s}.

Origin shift **doesn't** move a crystal in Cartesian coordinates, but
just changes the origin to measure the coordinates of atomic points.

   
A change of basis is described by the combination of the
transformation matrix and the origin shift denoted by
:math:`(\boldsymbol{P}, \boldsymbol{p})` where first the
transformation matrix is applied and then origin shift. The points in
the standardized system :math:`\boldsymbol{x}_\mathrm{s}` and
arbitrary system :math:`\boldsymbol{x}` are related by

.. math::

  \boldsymbol{x}_\mathrm{s} = \boldsymbol{P}\boldsymbol{x} +
  \boldsymbol{p},

or equivalently,

.. math::

  \boldsymbol{x} = \boldsymbol{P}^{-1}\boldsymbol{x}_\mathrm{s} -
  \boldsymbol{P}^{-1}\boldsymbol{p}.
  

A graphical example is shown below.

.. |cob| image:: change-of-basis.png
         :width: 20%

|cob|

(click the figure to enlarge)

In this example,

.. math::

   \renewcommand*{\arraystretch}{1.4}
   \boldsymbol{P} = \begin{pmatrix}
   \frac{1}{2} & \frac{1}{2} & 0 \\
   \frac{\bar{1}}{2} & \frac{1}{2} & 0 \\
   0 & 0 & 1 
   \end{pmatrix}.


.. _def_standardized_unit_cell:

Conventions of standardized unit cell
--------------------------------------

Choice of basis vectors
^^^^^^^^^^^^^^^^^^^^^^^^

Using the APIs ``spg_get_dataset``,
``spg_get_dataset_with_hall_number``, or ``spg_standardize_cell``, the
starndardized unit cell is obtained. The "starndardized unit cell" in
this document means that the (conventional) unit cell structure is
standardized by the crystal symmetry and lengths of basis vectors.
Crystals are categorized by Hall symbols in 530 different types in
terms of 230 space group types, unique axes, settings, and cell
choices. Moreover in spglib, lengths of basis vectors are used to
choose the order of :math:`(\mathbf{a}, \mathbf{b}, \mathbf{c})` if
the order can not be determined only by the symmetrical conventions.

.. _def_standardized_primitive_cell:

Transformation to the primitive cell
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In the standardized unit cells, there are five different centring
types available, base centrings of A and C, rhombohedral (R), body centred
(I), and face centred (F). The transformation is applied to the
standardized unit cell by

.. math::

   ( \mathbf{a}_\mathrm{p} \; \mathbf{b}_\mathrm{p} \; \mathbf{c}_\mathrm{p} )
   = ( \mathbf{a}_\mathrm{s} \; \mathbf{b}_\mathrm{s} \;
   \mathbf{c}_\mathrm{s} )  \boldsymbol{P}_\mathrm{c},

where :math:`\mathbf{a}_\mathrm{p}`, :math:`\mathbf{b}_\mathrm{p}`,
and :math:`\mathbf{c}_\mathrm{p}` are the basis vectors of the
primitive cell and :math:`\boldsymbol{P}_\mathrm{c}` is the
transformation matrix from the standardized unit cell to the primitive
cell. :math:`\boldsymbol{P}_\mathrm{c}` for centring types are given
as follows:

.. math::

   \renewcommand*{\arraystretch}{1.4}
   \boldsymbol{P}_\mathrm{A} = 
   \begin{pmatrix}
   1 & 0 & 0 \\
   0 & \frac{1}{2} & \frac{\bar{1}}{2} \\
   0 & \frac{1}{2} & \frac{{1}}{2}
   \end{pmatrix},
   \renewcommand*{\arraystretch}{1.4}
   \boldsymbol{P}_\mathrm{C} = 
   \begin{pmatrix}
   \frac{1}{2} & \frac{{1}}{2} & 0 \\
   \frac{\bar{1}}{2} & \frac{1}{2} & 0\\
   0 & 0 & 1
   \end{pmatrix},
   \boldsymbol{P}_\mathrm{R} = 
   \begin{pmatrix}
   \frac{2}{3} & \frac{\bar{1}}{3} & \frac{\bar{1}}{3} \\
   \frac{1}{3} & \frac{{1}}{3} & \frac{\bar{2}}{3} \\
   \frac{1}{3} & \frac{{1}}{3} & \frac{{1}}{3}
   \end{pmatrix},
   \boldsymbol{P}_\mathrm{I} = 
   \begin{pmatrix}
   \frac{\bar{1}}{2} & \frac{{1}}{2} & \frac{{1}}{2} \\
   \frac{{1}}{2} & \frac{\bar{1}}{2} & \frac{{1}}{2} \\
   \frac{{1}}{2} & \frac{{1}}{2} & \frac{\bar{1}}{2}
   \end{pmatrix},
   \boldsymbol{P}_\mathrm{F} = 
   \begin{pmatrix}
   0 & \frac{{1}}{2} & \frac{{1}}{2} \\
   \frac{{1}}{2} & 0 & \frac{{1}}{2} \\
   \frac{{1}}{2} & \frac{{1}}{2} & 0
   \end{pmatrix}.

For rhombohedral lattice systems with the choice of hexagonal axes,
:math:`\boldsymbol{P}_\mathrm{R}` is applied.

.. _def_idealize_cell:

Idealization of unit cell structure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Spglib allows tolerance parameters to match a slightly distorted unit
cell structure to a space group type with some higher symmetry. Using
obtained symmetry operations, the distortion is removed to idealize
the unit cell structure. The coordinates of atomic points are
idealized using respective site-symmetries (Grosse-Kunstleve *et
al*. (2002)). The basis vectors are idealized by forceing them into
respective lattice shapes as follows. In this treatment, except for
triclinic crystals, crystals can be rotated in Cartesian coordinates,
which is the different type of transformation from that of the
change-of-basis transformation explained above.

Triclinic lattice
""""""""""""""""""

- Niggli reduced cell is used for choosing :math:`\mathbf{a}, \mathbf{b}, \mathbf{c}`.
- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set in :math:`x\text{-}y` plane of Cartesian
  coordinates so that :math:`\mathbf{a}\times\mathbf{b}` is along
  :math:`+z` direction of Cartesian coordinates.

Monoclinic lattice
"""""""""""""""""""

- :math:`b` axis is taken as the unique axis.
- :math:`\alpha = 90^\circ` and :math:`\gamma = 90^\circ`
- :math:`90^\circ < \beta < 120^\circ`.

- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set along :math:`+y` direction of Cartesian coordinates.
- :math:`\mathbf{c}` is set in :math:`x\text{-}z` plane of Cartesian coordinates.

Orthorhombic lattice
"""""""""""""""""""""

- :math:`\alpha = \beta = \gamma = 90^\circ`.

- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set along :math:`+y` direction of Cartesian coordinates.
- :math:`\mathbf{c}` is set along :math:`+z` direction of Cartesian coordinates.

Tetragonal lattice
"""""""""""""""""""

- :math:`\alpha = \beta = \gamma = 90^\circ`.
- :math:`a=b`.

- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set along :math:`+y` direction of Cartesian coordinates.
- :math:`\mathbf{c}` is set along :math:`+z` direction of Cartesian coordinates.

Rhombohedral lattice
"""""""""""""""""""""

- :math:`\alpha = \beta = \gamma`.
- :math:`a=b=c`.

- Let :math:`\mathbf{a}`, :math:`\mathbf{b}`, and :math:`\mathbf{c}`
  projected on :math:`x\text{-}y` plane in Cartesian coordinates be
  :math:`\mathbf{a}_{xy}`, :math:`\mathbf{b}_{xy}`, and
  :math:`\mathbf{c}_{xy}`, respectively, and their angles be
  :math:`\alpha_{xy}`, :math:`\beta_{xy}`,
  :math:`\gamma_{xy}`, respectively.
- Let :math:`\mathbf{a}`, :math:`\mathbf{b}`, and :math:`\mathbf{c}`
  projected along :math:`z`-axis in Cartesian coordinates be
  :math:`\mathbf{a}_{z}`, :math:`\mathbf{b}_{z}`, and
  :math:`\mathbf{c}_{z}`, respectively.

- :math:`\mathbf{a}_{xy}` is set along :math:`+x` direction of Cartesian
  coordinates, and :math:`\mathbf{b}_{xy}` and :math:`\mathbf{c}_{xy}`
  are placed by angles :math:`120^\circ` and :math:`240^\circ` from
  :math:`\mathbf{a}_{xy}` counter-clockwise, respectively.
- :math:`\alpha_{xy} = \beta_{xy} = \gamma_{xy} = 120^\circ`.
- :math:`a_{xy} = b_{xy} = c_{xy}`.
- :math:`a_{z} = b_{z} = c_{z}`.


Hexagonal lattice
""""""""""""""""""

- :math:`\alpha = \beta = 90^\circ`.
- :math:`\gamma = 120^\circ`.
- :math:`a=b`.

- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set in :math:`x\text{-}y` plane of Cartesian coordinates.
- :math:`\mathbf{c}` is set along :math:`+z` direction of Cartesian coordinates.

Cubic lattice
""""""""""""""

- :math:`\alpha = \beta = \gamma = 90^\circ`.
- :math:`a=b=c`.

- :math:`\mathbf{a}` is set along :math:`+x` direction of Cartesian coordinates.
- :math:`\mathbf{b}` is set along :math:`+y` direction of Cartesian coordinates.
- :math:`\mathbf{c}` is set along :math:`+z` direction of Cartesian coordinates.
