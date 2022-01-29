<p align="center">
  <img width="256" src="mkdocs/docs/figs/logo_nap_256.png">
</p>

[![DOI](https://zenodo.org/badge/20047602.svg)](https://zenodo.org/badge/latestdoi/20047602)

# News

- **nap** was selected as **2020 Github archive program** and stored on magnetic tape(?) in the [Arctic vault](https://archiveprogram.github.com).

# What's nap
**Nagoya Atomistic-simulation Package (nap)** includes the following programs and utilities:
- parallel molecular dynamics simulation (*pmd*)
- potential parameter fitting (*fitpot* for neural-network potential and *fp.py* for other classical potentials)
- python modules for pre/post-processes (*nappy*)

The program, *pmd*, includes various interatomic potentials for metals and semiconductors,
and uses spatial decomposition technique for the parallelization, and linked-cell method for efficient neighbor search.

# Who made this?
* [Ryo KOBAYASHI](http://ryokbys.web.nitech.ac.jp/index.html)
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology.

# Requirements and dependencies

To compile *pmd* and *fitpot*, the following programs/libraries are required:

- Fortran compiler
- MPI library

The *nappy* and *fp.py* requires the following packages:

- *numpy*
- *scipy*
- *pandas*
- *docopt*
- *ASE*


# Compilation and usage

For the short test, whether or not you can use this program in your environment,

```bash
$ git clone https://github.com/ryokbys/nap.git
$ cd nap/
$ ./configure --prefix=$(pwd)
$ make test
```

If it works, you can use this program in your system.
To install the python package *nappy*,

```shell
$ python setup.py sdist
$ pip install -e .
```

Then you can use the nappy commands, `napsys` and `napopt`, in the terminal and can import `nappy` package in python programs.

For details, please see the [documentation](http://ryokbys.web.nitech.ac.jp/contents/nap_docs) or ask me via e-mail (kobayashi.ryo[at]nitech.ac.jp).

# Acknowledgements
This program was supported in part by ["Materials research by Information Integration" Initiative (MI2I)](http://www.nims.go.jp/MII-I/) project of the Support Program for Starting Up Innovation Hub from Japan Science and Technology Agency (JST).


# LICENSE
This software is released under the MIT License, see the LICENSE.

---
title: 'nap: A molecular dynamics package with parameter-optimization programs for classical and machine-learning potentials'
tags:
  - Fortran
  - Python
  - materials science
  - molecular dynamics
  - interatomic potential
  - neural-network potential
  - meta-heuristics
authors:
  - name: Ryo KOBAYASHI
    orcid: 0000-0001-8244-5844
    affiliation: 1
affiliations:
 - name: Department of Physical Science and Engineering, Nagoya Institute of Technology, Gokiso, Showa, Nagoya 466-8555, Japan
   index: 1
date: 28 August 2020
bibliography: paper.bib
---

# Summary

The **nap** is a package for molecular dynamics (MD) simulation consisting of an MD program (**pmd**) that can perform large-scale simulation using a spatial-decomposition technique and two parameter-optimization programs: one for classical (CL) potentials (**fp.py**) and another for machine-learning (ML) potentials (**fitpot**).
Since the numbers of parameters to be optimized are much different between CL and ML potentials, optimization approaches for them are also different; meta-heuristic global minimum-search algorithms for the CL potentials, in which the numbers of parameters are usually much less than one hundred, and gradient-based methods for the ML potentials.
The parameters of CL potentials can be optimized to any target quantity that can be computed using the potentials since meta-heuristic methods do not require the derivatives of the quantity with respect to parameters. On the other hand, ML-potential parameters can be optimized to only energies, forces on atoms and stress components of reference systems, mainly because gradient-based methods require the derivatives of other quantities with respect to parameters, and the analytical derivatives and the coding of them are usually painful and sometimes impossible.
Potentials can be used in combination with any other potential, such as pair and angular potentials, short-range and long-range potentials, CL and ML potentials.
With using the **nap** package, users can perform MD simulation of solid-state materials with the choice of different levels of complexity (CL or ML) by creating interatomic potentials optimized to quantum-mechanical calculation data even if no potential is available.

# Statement of need

MD simulation is widely used in many research fields such as materials science, chemistry, physics, etc., to study dynamics of atoms or molecules. In order to perform MD simulation of systems including large number of atoms, where quantum-mechanical calculations can not be used due to their computational cost, empirical interatomic potentials between species are required. And the results of MD simulation are strongly dependent on the property or accuracy of the potentials used in the simulation. Hence, there are a lot of CL potential models have been proposed such as Lennard-Jones (LJ) potential for van der Waals interaction, Coulombic potential for ionic interaction, Morse potential for covalent interactions [@Morse1929-ow], angular-dependent models for angles between covalent bonds, bond-order models for more complex systems, etc.
Recently ML potentials have been also actively studied because they are usually more flexible and can reproduce reference data more precisely than CL potentials.
Even though the potential is flexible or suitable to problems considered, the parameters in the potential model still need to be optimized to well reproduce the properties or phenomena that are in focus.

There are already a lot of MD programs such as LAMMPS [@Plimpton1995-az] and IMD [@Stadler1997-wr], and some of them can use both CL and ML potentials.
And also there are some parameter-optimization programs that can produce parameter sets available in the other MD programs such as aenet [@Artrith2016-mu] for ML potentials and potfit [@Brommer2007-kr; @Brommer2015-hw] for CL potentials.
However, there is also a demand of combining ML potentials and simpler CL potentials, e.g., ML potential with Coulomb interactions [@Morawietz2013-qq] and ML potential with core-core repulsions [@Wang2019-py], since creating an ML potential that covers very short-range and/or very long-range interactions is very inefficient. Thus it is beneficial if the programs of parameter-optimization for both CL and ML potentials are in one package and highly connected to one MD program and it will be more efficient than using several different programs to optimize parameters of CL and ML potentials.


# Programs and functionalities

- **pmd**: A Fortran program to perform MD simulation. It can perform large-scale MD simulation efficiently by using linked-cell list [@Allen2017-nw] and spatial decomposition with MPI parallelization. Several widely-used interatomic potentials for solid-state materials are implemented, for example, pair potentials such as LJ, Morse, Coulomb and screened Coulomb; angular-dependent potentials such as Stillinger-Weber [@Stillinger1985-gy] and Tersoff [@Tersoff1988-bo]; multi-body potentials such as Finnis-Sinclair [@Finnis1984-hv] and embeded-atom method [@Daw1984-hj]; machine-learning potentials such as linear-regression [@Seko2014-mh] and neural-network (NN) [@Behler2007-hr].
- **fp.py**: A Python program to optimize the parameters of classical potentials using meta-heuristic algorithms such as cuckoo search [@Yang2009-at]. Meta-heuristic algorithms need to perform MD simulations and to evaluate the quantities to be compared with target values of each individual representing a candidate parameter set. These jobs by individuals are performed as child processes of the program by calling a shell script that describes what to compute using given parameters. This program has a functionality of automatic update of the search range of parameters and it allows to optimize parameters efficiently and less dependent on the initial ranges of parameters [@Kobayashi2020-rz].
- **fitpot**: A Fortran program to optimize the parameters of ML potentials, such as linear regression and NN, using gradient-based methods [@Kobayashi2017-ky]. It can use high-performance computer resources to evaluate energies, forces and stresses of a huge number of sample systems in parallel using the MPI library, which is sometimes crucial since in general ML models require a lot of sampling to avoid over-fitting and make them robust.


# Acknowledgements

This work was supported in part by JSPS KAKENHI Grant Number JP20H05290 (Grant-in-Aid for Scientific Research on Innovative Areas "Interface Ionics") and by the "Materials Research by Information Integration" Initiative (MI2I) project of the "Support Program for Starting Up Innovation Hub" from JST.

# References
# nap/examples/pmd_W

This example shows how to perform an MD simulation of bcc-W system using the Ito potential[1]. The files required to run the MD simulation in this case are:

- `pmdini` -- information about simulation cell and positions of atoms
- `in.pmd` -- information about MD simulation setting

To perform MD,
```bash
$ /path/to/pmd | tee out.pmd
```

Users can check the result by looking at the output `out.pmd` and compare it with `out.pmd.REF` which is the reference output. Users may be compare them by looking at `Potential` and confirm that the results are identical.
```bash
$ grep 'Potential' out.pmd*
out.pmd:   Potential energy=       -463.72014 eV =     -8.587 eV/atom
out.pmd:   Potential energy=       -460.77567 eV =     -8.533 eV/atom
out.pmd.REF:   Potential energy=       -463.72014 eV =     -8.587 eV/atom
out.pmd.REF:   Potential energy=       -460.77567 eV =     -8.533 eV/atom
```


## References

1. A.M. Ito, Y. Yoshimoto, S. Saito, A. Takayama, H. Nakamura, Phys. Scripta T159, 014062 (2014).
# nap/example/fp_LATP

This example shows how to use `fp.py` to optimize the parameters in interatomic potential for Li-Al-Ti-P-O system. The potential forms are screened Coulomb, Morse and angular and the target quantities are RDF, ADF and equilibrium volume, see the Ref. [1] for details. The files required to perform `fp.py` are:

- `in.fitpot` -- information about `fp.py` setting
- `in.vars.fitpot` -- initial values of potential parameters and search ranges of the parameters
- `in.params.Coulomb` -- parameter file for Coulomb potential
- `data.ref.xxx` -- target quantities (RDF, ADF and volume)
- `subjob.sh` -- shell-script that describe what to compute target quantities with the current parameters and
- `in.pmd.NpT`, `pmdini`, and something -- files need to perform MD simulations described in `subjob.sh`

To see command-line options,
```bash
$ python /path/to/fp.py -h
```

To perform `fp.py` using four threads,
```bash
$ /path/to/fp.py --nproc=4 | tee out.fp
```

Since the `fp.py` uses random values to generate trial parameters and the random values would be different every time it is executed, the output obtained will not be identical to that in `out.fp.REF`. But it will be OK if `tail out.fp` looks like the following,
```bash
$ tail out.fp
   iid,losses=        6     1.2219     3.3689     0.0120     4.60281
   iid,losses=        9     0.6474     7.5324     0.1334     8.31316
 step,time,best,vars=      1    173.1    2.2234  1.000  0.748  0.783  1.019  1.071  0.689  0.883  1.629  1.829  1.615  1.951  2.442  1.606  1.739  1.791  3.519
   iid,losses=       11     1.7628     0.2426     2.4705     4.47583
   iid,losses=       12     2.0236     1.5161     0.7651     4.30480
   iid,losses=       10     1.2067     0.6369     0.3003     2.14386
   iid,losses=       13     0.6055     9.3178     0.1672    10.09049
   iid,losses=       14     1.0248     8.7687     0.0168     9.81031
 step,time,best,vars=      2    346.8    2.1439  1.000  0.785  0.785  1.030  1.090  0.687  0.911  1.627  1.825  1.615  1.957  2.442  1.653  1.770  1.790  3.429
elapsed time = 346.771660 sec.
```
And when you use `fp.py`, you had better use more processes than 4 like in this case to efficiently run the program.

## Differences from `examples/fp_LZP`

This example used **whatever mode** for target quantities not like **DF-matching mode** in `examples/fp_LZP`. Differences from `examples/fp_LZP` include:

- The format of `data.ref.xxx` files
- `match` entry in `in.fitpot`
- `subjob.sh` script (`--out4fp` option is specified in this case.)


## References

1. R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, “High-throughput production of force-fields for solid-state electrolyte materials”, APL Materials 8, 081111 (2020). [link](https://aip.scitation.org/doi/10.1063/5.0015373).
# nap/example/fp_LZP

This example shows how to use `fp.py` to optimize the parameters in interatomic potential for Li-Zr-P-O system. The potential forms are screened Coulomb, Morse and angular and the target quantities are RDF, ADF and equilibrium volume, see the Ref. [1] for details. The files required to perform `fp.py` are:

- `in.fitpot` -- information about `fp.py` setting
- `in.vars.fitpot` -- initial values of potential parameters and search ranges of the parameters
- `in.params.Coulomb` -- parameter file for Coulomb potential
- `data.ref.xxx` -- target quantities (RDF, ADF and volume)
- `subjob.sh` -- shell-script that describe what to compute target quantities with the current parameters and
- `in.pmd.NpT`, `pmdini`, and something -- files need to perform MD simulations described in `subjob.sh`

To see command-line options,
```bash
$ python /path/to/fp.py -h
```

To perform `fp.py` using four threads,
```bash
$ /path/to/fp.py --nproc=4 | tee out.fp
```

Since the `fp.py` uses random values to generate trial parameters and the random values would be different every time it is executed, the output obtained will not be identical to that in `out.fp.REF`. But it will be OK if `tail out.fp` looks like the following,
```bash
$ tail out.fp
   iid,Lr,Lth,Lvol,Llat,L=        1    1.3376     0.7049     7.7690     0.0000     9.8116
   iid,Lr,Lth,Lvol,Llat,L=        3    0.8912     0.6027     4.7938     0.0000     6.2877
 step,time,best,vars=      0     40.3    6.2877  1.000  1.142  1.255  1.054  0.868  1.208  2.106  1.778  2.988  2.044  1.921  4.330  2.013  1.536  2.486  1.000
   iid,Lr,Lth,Lvol,Llat,L=        8    3.1944     0.7949    25.6610     0.0000    29.6503
   iid,Lr,Lth,Lvol,Llat,L=        7    1.8878     0.7908    19.5151     0.0000    22.1937
   iid,Lr,Lth,Lvol,Llat,L=        6    1.2444     0.7123     7.4910     0.0000     9.4477
   iid,Lr,Lth,Lvol,Llat,L=        5    0.8285     0.6159     5.3013     0.0000     6.7457
   iid,Lr,Lth,Lvol,Llat,L=        9    2.0158     0.7721    15.0153     0.0000    17.8033
 step,time,best,vars=      1    104.3    6.2877  1.000  1.142  1.255  1.054  0.868  1.208  2.106  1.778  2.988  2.044  1.921  4.330  2.013  1.536  2.486  1.000
elapsed time = 104.290595 sec.
```
And when you use `fp.py`, you had better use more processes than 4 like in this case to efficiently run the program.

## References

1. R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, “High-throughput production of force-fields for solid-state electrolyte materials”, APL Materials 8, 081111 (2020). [link](https://aip.scitation.org/doi/10.1063/5.0015373).
# nap/examples/pmd_DNN_SiO

This example shows how to perform an MD simulation of SiO2 system using a neural-network (NN) potential. The files required to run the MD simulation in this case are:

- `pmdini` -- information about simulation cell and positions of atoms
- `in.pmd` -- information about MD simulation setting
- `in.params.desc`, `in.params.DNN`, `in.params.ZBL` -- parameter files for the NN potential for Si-O.

To perform MD,
```bash
$ /path/to/pmd | tee out.pmd
```

Users can check the result by looking at the output `out.pmd` and compare it with `out.pmd.REF` which is the reference output. Users may be compare them by looking at `Potential` and confirm that the results are identical.
```bash
$ grep 'Potential' out.pmd*
out.pmd:   Potential energy=        -71.13849 eV =     -7.904 eV/atom
out.pmd:   Potential energy=        -70.97027 eV =     -7.886 eV/atom
out.pmd.REF:   Potential energy=        -71.13849 eV =     -7.904 eV/atom
out.pmd.REF:   Potential energy=        -70.97027 eV =     -7.886 eV/atom
```

And there are also some `dump_###` files that contain snapshot configurations of atoms and users can see the movie of the simulation using some visualization software such as [Ovito](https://www.ovito.org/about/).
# nap/example/fitpot_DNN_SiO

This example shows how to use `fitpot` to optimize the parameters in neural-network potential for Si-O system, see the Ref. [1,2] for details of NN potential. The files required to perform `fitpot` are:

- `in.fitpot` -- information about `fitpot` setting.
- `in.params.desc` -- information about descriptors used as inputs for the NN potential.
- `in.params.DNN` -- initial values of potential parameters (NN weights) and search ranges of the parameters
- `in.params.ZBL` -- classical potential added to the NN potential (in this case, ZBL potential)
- `dataset/smpl_XXX/` -- reference dataset

To perform `fitpot` using 2 MPI processes,
```bash
$ mpirun -np 2 /path/to/fitpot | tee out.fitpot
```

If the tail of output shows like the following, at least `fitpot` program finished corectly without errors.
```bash
$ tail out.fitpot
 ENERGY:      100          16.58    0.0000479    0.0000000    0.0000576    0.0000000    1.0000000    0.0000000
 FORCE:       100          16.59    0.0020653    0.0000000    0.0049990    0.0000000    0.9999964    0.0000000
 STRESS:      100          16.59   39.5048633    0.0000000   98.5031577    0.0000000  -99.7924332    0.0000000
 Number of func and grad calls =  203  101
 Memory/proc =         0.420 MB
 Time func =           1.423 sec
 Time grad =          14.367 sec
 Time comm =           0.001 sec
 Time      =          16.588 sec  =   0h00m16s
 Job finished at 10:25:05 on 2020-08-20
```
Please see the nap documentation for more details how to use `fitpot` to obtain the NN potential parameters.

## References

1. Behler, Jörg, and Michele Parrinello. 2007. “Generalized Neural-Network Representation of High-Dimensional Potential-Energy Surfaces.” Physical Review Letters 98 (14): 146401–146401
2. Kobayashi, Ryo, Daniele Giofré, Till Junge, Michele Ceriotti, and William Arthur Curtin. 2017. “Neural Network Potential for Al-Mg-Si Alloys.” Physical Review Materials 1 (5): 53604–11.[link](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.053604)
# Input file: in.pmd

*pmd* starts with reading setting file `in.pmd` and atom-configuration
files `pmdini`. So simulation settings except the atom configuration
must be described in `in.pmd`.

## Example of in.pmd


    #
    #  unit of time  = femto sec
    #  unit of length= Angstrom
    #  unit of mass  = unified atomic mass unit
    #
      io_format         ascii
      print_level       1

      time_interval     2d0
      num_iteration     1000
      min_iteration     10
      num_out_energy    100

      flag_out_pos      1
      num_out_pos       10

      force_type        Morse Coulomb
      cutoff_radius     5.8d0
      cutoff_buffer     0.2d0

      flag_damping      2
      damping_coeff     0.5d0
      converge_eps      1d-4
      converge_num      3

      initial_temperature    -2000d0
      final_temperature      -2000d0
      temperature_control     none
      temperature_target      1  100d0
      temperature_relax_time  1d0

      factor_direction 3 2
        1.000d0  1.000d0  1.000d0
        1.000d0  0.000d0  1.000d0

      stress_control      none
      stress_relax_time   100d0
      stress_target
        0.00d0   0.00d0   0.00d0
        0.00d0   0.00d0   0.00d0
        0.00d0   0.00d0   0.00d0
      pressure_target     1.00

      boundary   ppp

Here, lines begin with `!` or `#` are treated as comment lines.


## Input parameters


### num_nodes_x

- Default: `-1`

Number of division in x, y, or z direction. If one of these is
non-positive (`<=0`), these numbers are automatically estimated from the
system size and the number of MPI processes used. If all of these are
positive, specified values are used. The product of these, $xyz$, should
be the same as the number of divided atom-configuration files and
computer nodes specified when executing *mpirun* or *mpiexec* command.

------------------------------------------------------------------------

### io_format

- Default: `ascii`

You can choose either `ascii` or `binary` format of atom-configuration
files. When you perform large scale simulation, you should choose
`binary` for efficient reading/writing atom-configuration files.

------------------------------------------------------------------------

### print_level

- Default: `1`

How much information is written out during the run.

- `1` -- Normal information for MD simulation run.
- `100` -- Debug info is written out.

------------------------------------------------------------------------

### time_interval

-   Default: 1.0

Time interval in the unit of **femto second**. If negative, it activates
*variable time-step mode* and its absolute value is the maximum time
interval $\Delta t_{\max}$ in the mode. Appropriate range of dt_max
would be 2.0 to 5.0 depending on the minimum mass of ion in the system.

------------------------------------------------------------------------

### vardt_length_scale

-   Default: 0.1

The specific length $L^*$ of the *variable time-step mode* where the
time interval is determined as

$$\begin{equation}
   \Delta t = \min \left( \Delta t_\mathrm{max}, \frac{L^*}{v_\mathrm{max}}\right).
\end{equation}$$

------------------------------------------------------------------------

### num_iteration

- Default: 0
- Alternative: `num_steps`

Number of MD steps. Simulation time equals `time_interval` times
`num_iteration`.

------------------------------------------------------------------------

### min_iteration

- Default: 0
- Alternative: `min_steps`

Minimum number of MD steps. In the case you want the MD simulation at
least *min_iteration*, you can set this parameter.

------------------------------------------------------------------------

### num_out_energy

- Default: 1000

Number of outputs of energies.

------------------------------------------------------------------------

### flag_out_pos

- Default: 1
- Alternative: `flag_out_pmd`

A flag whether or not to write atomic configurations to files at certain
steps.

- `0` -- Not to write.
- `1` -- Write *pmd*-format atomic configurations to files `pmd_####` where `####` indicates sequential number of the files.
- `2` -- Write LAMMPS *dump*-fomrat atomic configurations to files `dump_####`.

------------------------------------------------------------------------

### num_out_pos

- Default: 10
- Alternative: `num_out_pmd`

Number of atom-configuration files to be written.

------------------------------------------------------------------------

### flag_sort

- Default: `1`

A flag whether or not to sort the order of atoms by tag before writing
out atomic configurations to *pmd* files. It might cost some time for
large scale simulation.

- `1` -- Do sorting
- `2` -- Do not sorting

------------------------------------------------------------------------

### force_type

- Default: `None`

Choice of the interatomic potential. Available potentials are listed
below:

- `LJ` : Lennard-Jones potential for Ar system.
- `SW_Si` : Stillinger-Weber potential for Si system.
- `EDIP_Si` : Environment Dependent Interatomic Potential for Si system.
- `Ito3_WHe` : EAM potential for W-He system made by Ito et al. at NIFS.
- `Morse` : Morse potential that requires an input files `in.params.Morse`.
- `Coulomb` : Coulomb potential that requires an input files `in.params.Coulomb`.
- `DNN` : Neural-network potential that requires two input files `in.params.desc` and `in.params.NN2`.

------------------------------------------------------------------------

### flag_damping

- Default: `0`

A flag whether or not damp atom velocities.

- `0` -- No damping.
- `1` -- Simple damped MD using the following **damping_coeff**.
- `2` -- FIRE algorithm. This is usually much faster and stabler than the simple damped MD.

------------------------------------------------------------------------

### damping_coeff

- Default: `0.9`

Damping coefficient.

------------------------------------------------------------------------

### converge_eps

- Default: `1d-4`

Convergence criterion in eV. If it is negative value, not to stop because of convergence.

------------------------------------------------------------------------

### converge_num

- Default: `1`

Convergence is achieved if the convergence criterion is sufficed this times successively.

------------------------------------------------------------------------

### initial_temperature

- Default: `-1.0`

Initial temperature of all atoms.

------------------------------------------------------------------------

### final_temperature

- Default: `-1.0`

Final temperature of all atoms. If it is set, target temperature of all
atoms changes linearly from `initial_temperature` as simulation
proceeds.

------------------------------------------------------------------------

### temperature_control

- Default: `none`

Temperature-control method, `none`, `Berendsen`, and `Langevin` are now
available.

------------------------------------------------------------------------

### temperature_target

- Default: `300.0`

Target temperature (K) of atoms specified by *ifmv*. For example, :

    temperature_target   1  300.0

indicates setting the temperature of *ifmv=1* to 300 K.

------------------------------------------------------------------------

### temperature_relax_time

- Default: `100.0`

Relaxation time of Berendsen thermostat (fs).

------------------------------------------------------------------------

### stress_control

- Default: `none`

Type of barostat. Following methods are available:

- `Berendsen` / `vc-Berendsen`: variable-cell Berendsen method.
- `vv-Berendsen`: variable-volume Berendsen method which keeps the
    cell shape but changes volume by scaling all the cell vectors.

See Berendsen's paper[^Berendsen1984] for the detail.

[^Berendsen1984]: Berendsen, C., Postma, P. M., Van Gunsteren, W. F., Dinola, A., & Haak, R. (1984). Molecular dynamics with coupling to an external bath,
81(8), 3684--3690.

------------------------------------------------------------------------

### pressure_target

Default: `0.0`

Target hydrostatic pressure which only works when *stress_control* is
`vv-Berendsen`.

------------------------------------------------------------------------

### stress_target

Default:

     0.00d0  0.00d0  0.00d0
     0.00d0  0.00d0  0.00d0
     0.00d0  0.00d0  0.00d0

Target stress tensor which only works when *stress_control* is
`Berendsen` or `vc-Berendsen`.

------------------------------------------------------------------------

### stress_relax_time

- Default: `100d0`

Relaxation time (fs) of the *Berendsen* barostat.

------------------------------------------------------------------------

### zload_type

- Default: `no`

How to apply z-direction strain:

- `atoms` -- atoms whose ifmv value are 2 and relative z-position are
    over 0.5 are moved to upward, those of relative z-position under 0.5
    are moved downward.
- `box` -- control z-component of simulation box matrix.
- `no` -- Not to apply z-direction strain.

------------------------------------------------------------------------

### final_strain

- Default: `0.0`

Final strain value (%). Thus strain rate can be given as `final_strain`
/ ( `time_interval` * `num_iteration` ).

------------------------------------------------------------------------

### shear_stress

- Default: `0.0`

Shear stress value applied to the system.

------------------------------------------------------------------------

### cutoff_radius

- Default: `5.0`

Cutoff radius (Angstrom) of the interatomic potential used.

------------------------------------------------------------------------

### flag_temp_dist

- Default: `.false.`

Flag about whether or not writing out temperature distribution to
`out.temp-dist` file.

------------------------------------------------------------------------

### num_temp_dist

- Default: `1`

Number of bins along *x*-direction where the temperature is calculated.
This value must be a multiple of `num_nodes_x`.

------------------------------------------------------------------------

### mass

- Default: masses of given species are set automatically.

If masses of some species need to be set different from those of
elements, masses should be specified as follows. :

    mass   Si  28.0855
    mass   He   4.00

------------------------------------------------------------------------

### boundary

- Default: `ppp`

Boundary conditions for each axis, 123.

- `p`: periodic boundary condition
- `f`: free boundary condition
# pmd usage

## Install

### Download

Download the source code of current version of the whole packange from
the Github site, [https://github.com/ryokbys/nap](https://github.com/ryokbys/nap).


### Compile

You can download `nap-master.zip` file from the site. And you can get
`nap-master` directory when you unzip the zip file. For the ease of
following explanation, change the directory name to `nap` as,

    $ unzip nap-master.zip
    $ mv nap-master nap
    $ cd nap/

Then you can compile the *pmd* program as following,

    $ ./configure --prefix=$(pwd)
    $ cd pmd
    $ make pmd

If you get an error of finding an MPI Fortran compiler when you are
running `configure`, you have to find an MPI Fortran compiler by
yourself or asking the system administrator and do `configure` again by
specifying the compiler position as,

    $ ./configure --prefix=$(pwd) FC=<path/to/fortran-compiler>

The option `--prefx=$(pwd)` is not necessary for *pmd*, but when you use
`fitpot` program, it is required.

!!! Note
    If you get an error related to C preprocessor such as,
    
        configure: error: C preprocessor "/lib/cpp" fails sanity check
    
    you may have to specify true C preprocessor path to `configure` command
    by adding an option like `CPP=/usr/bin/cpp`.


Although the *pmd* command should be available with this compilation,
this *pmd* may not be optimized to the system in which it is compiled.
You might need to add some options relevant for the system in which it
is compiled.

#### gfortran and openmpi

In the case of `gfortran` with `openmpi`,

    $ ./configure --prefix=$(pwd) FCFLAGS="-O2 -g"

The optimization option over `-O3` seems to cause some errors depending
on the version of gfortran, so it is recommended to use `-O2`
optimization option.

When you are debugging, you may had better set some warning and checking
options enabled as follows.

    $ ./configure --prefix=$(pwd) FCFLAGS="-O2 -g -fbounds-check -Wuninitialized -fbacktrace"

!!! Note
    Compilation with LLVM version gcc is not tested. Use Homebrew version of gcc and openmpi.


#### Intel Fortran compiler

If you can use Intel Fortran Compiler`ifort` in your system, the
configure command would be like,

    $ ./configure --prefix=$(pwd) FCFLAGS="-xHOST -O3 -ip -ipo -no-prec-div"

The options `-ip` and `-ipo` have to do with inline expansions and are
relevant to the efficiency of *pmd*.

#### PGI fortran compiler

If the MPI fortran command `mpif90` is linked to PGI fortran compiler,
you can use the compiler by just specifying the compiler path as,

    $ FC=/path/to/mpif90 ./configure --prefix=$(pwd) FCFLAGS='-Minfo -O2 -g'



#### Fujitsu Fortran in FX?

It is easier to compile on the computation node not on the login node.
Since there are some difference about configuring/compiling on those
nodes. To configure and compile the *pmd*, first you need to login to a
computation node by doing `pjsub --interact`.

    $ pjsub --interact
    or
    $ pjsub --interact -L rscgrp=fx-interactive,node=1  <== in case of flow-fx@nagoya-u
    $ ./configure --prefix=$(pwd) FCFLAGS="-O3"
    $ cd pmd
    $ make pmd
    $ exit

!!! Note
    In case that the `configure` returns errors and exit without completing
    the configuration and the error message is related to cross compilation,
    you may need to add an option like `--host=sparc64` to the above command
    line.


#### Fujitsu Fortran in CX400

In the case of Fujitsu Fortran compiler `mpifrt` in CX400,

    $ ./configure --prefix=$(pwd) FCFLAGS="-Kfast,parallel"


#### Helios in Rokkasho-mura

It is Linux OS on Intel CPU, and the compilation seems to be basic one.
But one needs to add specific options as following,

    $ ./configure --prefix=$(pwd) FC=mpiifort FCFLAGS="-xAVX -O3 -ip -ipo -g -CB"

If you don\'t specify the `mpiifor` explicitly, `ifort` is set by
default and the compilation does not work correctly.

------------------------------------------------------------------------

## Run a sample simulation

There are some input files in `example/test-W/` directory. ( `in.pmd`
and `pmdini` ) These input files are for the system of BCC tungsten
crystalline structure including one helium atom.

    $ cd example/test-W
    $ ../../pmd/pmd


When you run the *pmd* command like above, *NVE* -MD simulation of 100
steps is performed. And the total, kinetic, and potential energies are
output in `out.erg` file. So you can look at the evoluation of these
energies using `gnuplot` command as,

    $ gnuplot
    gnuplot> plot 'out.erg' us 1:3 w l, 'out.erg' us 1:4 w l, 'out.erg' us 1:5 w l

In this case, since you are performing *NVE* -MD simulation of bcc-W,
the total energy conserves conpensating the deviations of kinetic and
potential energies.

![image](./figs/graph_energy-steps.png)

!!! Note
    The format ot `out.erg` is a bit changed from that of before 2018-11-01
    versions. The total and potential energies are raw values not being
    subtracted the initial values.

And also configurations of atoms at each 10 steps out of 100 steps are
written in LAMMPS-dump format, *e.g.*, `dump_0`, `dump_10`,\...,
`dump_100`.

------------------------------------------------------------------------

## Input files needed to run pmd

To run *pmd*, the following files are required in the working directory,

- `in.pmd` -- Input file that describes simulation setting.
- `pmdini` -- Cell information and initial positions and velocities of atoms.

And there are some optional files required by the *pmd* if you use
interatomic potentials that require input parameters from files such as
`in.params.xxx`.

![image](./figs/pmd.png)

After running *pmd* , some output files appear in the same directory.

------------------------------------------------------------------------

## Units used in pmd

- Length: Angstrom
- Time: femto second (fs)
- Energy: electron volt (eV)
- Mass: 1/12 of carbon atom

------------------------------------------------------------------------

## Make an initial atom-configuration file

Please [Atom-configuration file](pmd-file.md) for detail.

One has to make an initial atom-configuration file, `pmdini`, to run
*pmd*. There are already some programs that make initial
atom-configuration files of some systems (`mkconf/mkconf_????.F` and/or
`nappy/mkcell/cell_maker.py`). You can make your own initial
atom-configuration file by looking at those program codes.

If there is already a program that makes an atom-configuration file of
your target system, you can make an atom-configuration file as,

    $ cd mkconf
    $ emacs makefile
    (find which mkconf_* will be made)
    $ make mkconf_Si_disl
    $ ./mkconf_Si_disl

or you can use `nappy/mkcell/cell_maker.py` as well,

    $ python /path/to/nappy/mkcell/cell_maker.py -h
    ...
    $ python /path/to/nappy/mkcell/cell_maker.py dia -l 5.427 -s 4,4,4

Then you get an atom-configuration file `pmdini`.

!!! Note
    If you have to make the program which makes an atom-configuration file,
    copy any program like `mkconf_BCC.F` , modify it, add an entry into
    `makefile` , and compile.

------------------------------------------------------------------------

## Make the in.pmd file

Please refer `in-pmd`{.interpreted-text role="ref"} for details of
`in.pmd` file.

For instance, `in.pmd` file for the system of 1000 step MD simulation
using `SW_Si` potential is as follows,

    #
    #  unit of time  = femto sec
    #  unit of length= Angstrom
    #  unit of mass  = unified atomic mass unit
    #

    io_format         ascii
    print_level       1

    time_interval     2d0
    num_iteration     1000
    num_out_energy    100

    flag_out_pmd      1
    num_out_pmd       10

    force_type        SW_Si
    cutoff_radius     3.7712d0
    cutoff_buffer     0.2d0

    flag_damping      2
    damping_coeff     0.5d0
    converge_eps      1d-4
    converge_num      3

    initial_temperature     -2000d0
    final_temperature     -2000d0
    temperature_control     none
    temperature_target      100d0
    temperature_relax_time  1d0

    factor_direction 3 2
      1.000d0  1.000d0  1.000d0
      1.000d0  0.000d0  1.000d0

    stress_control       none
    stress_relax_time   100d0
    stress_target
      0.00d0   0.00d0   0.00d0
      0.00d0   0.00d0   0.00d0
      0.00d0   0.00d0   0.00d0
    pressure_target     1.00

    shear_stress   0.00

Here, the lines begin with `!` or `#` are treated as comment lines and
blanc lines are skipped.

------------------------------------------------------------------------

## Run pmd

### Run pmd on 1-process

It is really easy to run *pmd* on 1-process. On the directory where
`in.pmd` and `pmdini` exist, just execute *pmd* as,

    $ /path/to/pmd/pmd

If you want to perform it background,

    $ /path/to/pmd/pmd > out.pmd 2>&1 &

The following files appear when you perform *pmd*,

- `out.erg` -- Total, kinetic, potential energies, and temperature,
    volume, pressure.
- `dump_##` -- Atom-configurations at a certain MD step is written in
    LAMMPS-dump format by default. `##` means the MD step.

### Run pmd on parallel-nodes

Different from the old version of *pmd* which requires divided atom
configuration files for parallel nodes, in the current version (since
2016-05-05), the parallel simulation can be performed almost the same as
the serial run.

Just you need to describe how many divisions on each direction in
`in.pmd` such as `num_nodes_x`, `num_nodes_y` and `num_nodes_z` ,and run
*pmd* with `mpirun` or `mpiexec` command to run MPI executable.

    $ mpirun -np 8 /path/to/pmd > out.pmd 2>&1 &

Here, *pmd* will be executed on 8-nodes and the standard output is
written into `out.pmd` .

If any job-scheduling system is available on the system you are using,
describe the above command in your job script to be submitted.

------------------------------------------------------------------------

## Notes on performing massively parallel simulation

When you perform parallel simulation with over one million atoms, the
data of atom-configuration files becomes considerably large and
reading/writing data takes long time compared with intrinsic
computation. So *pmd* can read/write binary version of
atom-configuration files that are way more small amount of data. If you
want to read/write binary files, describe following in `in.pmd`,

    io_format   binary

And also you have to write code of writing binary atom-configuration
file in `mkconf_*.F`. In `mkconf_*.F` files, there is a line of
subroutine call `call write_pmd0_ascii` , you have replace it to
`call write_pmd0_bin` and recompile the program.
# Force fields

Force fields (FFs) are specified at `force_type` keyword in `in.pmd`
file. Plural FFs can be specified as a space-separated list,

    force_type    NN Morse Coulomb

Each FF reads one or some parameter files in the working directory,
which are usually named like `in.params.Coulomb` or so, specified by FF.
For example, **DNN** requires two files, `in.params.DNN` (file for
weights in NN) and `in.params.desc` (file for descriptor information).

Available FFs are listed below:

- [force_fields/DNN](force_fields/DNN.md)
- [force_fields/Coulomb](force_fields/Coulomb.md)
- [force_fields/Morse](force_fields/Morse.md)
- [force_fields/SW](force_fields/SW.md)
- Lennard-Jones (Ar)
- Linear regression

# nap documentation

This is a documentation of **Nagoya Atomistic-simulation Package
(nap)**.

-   **nap** is hosted at Github: <https://github.com/ryokbys/nap>
-   **nap** includes the following programs and utilities:
    -   *pmd*: Fortran program of massively parallel molecular dynamics
    -   *fitpot*: Fortran program for fitting the parameters of
        neural-network potential
    -   *fp.py*: Python script for fitting the parameters of classical
        potentials
    -   *nappy*: Python scripts for pre- and post-processing
-   **Nagoya** is the name of a city in Japan, where the project
    started.

------------------------------------------------------------------------

## How to cite nap

- [R. Kobayashi, Y. Miyaji, K. Nakano, M. Nakayama, APL Materials 8, 081111 (2020)](https://aip.scitation.org/doi/10.1063/5.0015373) for *nap* package and/or *fp.py* program.
- [Kobayashi, Ryo, Daniele Giofré, Till Junge, Michele Ceriotti, and William Arthur Curtin. Physical Review Materials 1, 53604--11 (2017)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.053604) for *fitpot* program and/or neural-network (NN) potential.

------------------------------------------------------------------------

## License

The **nap** is distributed under the MIT license, please see the license statement in the GitHub repository.

------------------------------------------------------------------------

## Contact

**Ryo KOBAYASHI**  
kobayashi.ryo [at] nitech.ac.jp  
Department of Physical Science and Engineering, Nagoya Institute of Technology

------------------------------------------------------------------------

## Acknowledgements

This program was supported in part by ["Materials research by Information Integration" Initiative (MI2I)](http://www.nims.go.jp/MII-I/) project of the Support Program for Starting Up Innovation Hub from Japan Science and Technology Agency (JST).


# fp.py -- fit parameters of classical potentials

The python program *fp.py* is another program of fitting potential
parameters. The *fitpot* focuses on the neural-network (NN) potential,
on the other hand, this *fp.py* focuses on classical potentials that
have much less potential parameters compared to NN potential, usually
less than 100. Since the number of parameters to be optimized is small,
*fp.py* employs meta-heuristic or nature-inspired methods, which can
adopt any physical value as a learning target since derivatives of the
target values w.r.t. optimizing parameters are not required.

------------------------------------------------------------------------

## Setup

To use *fp.py*, you should check if *nappy* works correctly, since
*fp.py* is actually a part of *nappy* package (*fp.py* exists at
`nap/nappy/fitpot/fp.py`). See [nappy_setup](./nappy.html#setup) of *nappy* package.

------------------------------------------------------------------------

## What does fp.py do?

In the *fp.py*, the following loss function is minimized by optimizing
potential parameters,

$$\mathcal{L} = \sum_t \mathcal{L}_t$$

where $t$ stands for the type of target and it can be anything if it is
computed from *ab-initio* program and MD program using the potential. If
the target is the volume of the system,

$$\mathcal{L}_\mathrm{V} = \left( \frac{V_\mathrm{FF} -V_\mathrm{ref}}{V_\mathrm{ref}} \right)^2.$$

If the target is the radial distribution function (RDF),

$$\mathcal{L}_\mathrm{R} = \frac{1}{N_p} \sum_p^{N_p} \left( \frac{\sum_i^{N_\mathrm{R}} (g_{p,\mathrm{FF}}(r_i)- g_{p,\mathrm{ref}}(r_i))^2}{\sum_i^{N_\mathrm{R}} g_{p,\mathrm{ref}}^2(r_i)} \right),$$

where $N_p$ is the number of RDF pairs to be considered, $N_\mathrm{R}$
is the number of sampling points of RDF.

To minimize the above loss function, the following metaheuristic methods
are available:

-   Cuckoo search (CS)
-   Differential evolution (DE)

------------------------------------------------------------------------

## Quick trial with examples

There are two examples of *fp.py* in `nap/examples/`,

-   `fp_LZP/` -- fitting of parameters of Morse, Coulomb and SW-like
    angular potentials to RDF, angular distribution function (ADF) and
    equilibrium volume of Li-Zr-P-O system.
-   `fp_LATP/` -- fitting of parameters of Morse, Coulomb and SW-like
    angular potentials to RDF, ADF and equilibrium volume of
    Li-Al-Ti-P-O system.

In either of these two directories, you can try *fp.py* by running the
following command,

    $ python /path/to/nap/nappy/fitpot/fp.py --nproc 4 | tee out.fp

This could take a few minutes using 4 processes. And you can see some
output files written by *fp.py*. See how to discuss [result](#results-and-outputs).

------------------------------------------------------------------------

## Files needed to run fp.py

As you can see in `nap/examples/fp_LZP/` or `nap_examples/fp_LATP/`,
there are several files needed to run the *fp.py* program.

-   `in.fitpot` -- *fp.py* configuration file
-   `in.vars.fitpot` -- optimizing parameter file
-   `in.params.XXX` -- potential parameter files that are not
    actually used during optimization, except that `in.params.Coulomb`
    must exist becuase charge information of each element is read from
    it.
-   `data.ref.XXX` -- Reference data file
-   `subjob.sh` -- Shell-script file to perform sub jobs
-   Files needed to perform sub jobs:
    -   `pmdini` -- atom configuration file (cell info, positions,
        and velocities)
    -   `in.pmd.NpT` -- input file for *pmd*

## Reference data for fp.py

Reference data should be stored in files `data.ref.XXX` where `XXX`
indicates the type of target quantity, such as `rdf`, `adf`, `vol`, etc.

Currently (July 2020), *fp.py* can run in two different modes,

-   **distribution function (DF)-matching mode** -- RDF, ADF,
    equilibrium volume and lattice constants are adopted as targets.
-   **whatever mode** -- any quantity that is computable can be used
    as a target.

In the **DF-matching mode**, the reference data formats of these targets
are different.

On the other hand, in the **whatever mode**, the format of reference
data is fixed.

### Reference data format in DF-matching mode

See the format of each target (RDF, ADF, vol and lat) in
`nap/examples/fp_LZP/`.

### Reference data format in whatever mode

The format of reference data as follows, :

    # comment line begins with `#`
    #
        100    1.0
        0.1234   0.2345  0.3456  0.4567  0.5678  0.6789
        0.7890   0.8901  0.9012  0.0123  0.1234  0.2345
        ...

-   Lines begining with `#` **at the head of the file** are treated as
    comment lines.
-   1st line -- `NDAT` the number of data, `WGT` the weight for this
    target.
-   2nd line and later -- data, no limitation to the number of entries
    in a line, but it is recommended to include 6 entries in a line.

## Input file in.fitpot

Control parameters for *fp.py* are read from `in.fitpot` in the working
directory. There are some diffierences between **DF-matching mode** and
**whatever mode** in `in.fitpot`.

First, in the case of **whatever mode**, the `in.fitpot` in the example
`nap/examples/fp_LATP/` is shown below,

    num_iteration      100
    print_level         1

    fitting_method   cs
    sample_directory "./"
    param_file in.vars.fitpot

    match     rdf adf vol
    potential   BVSx

    cs_num_individuals   20
    cs_fraction          0.25
    update_vrange        10
    fval_upper_limit     100.0

    specorder  Li Al Ti P O

    interactions  7
      Li  O
      Al  O
      Ti  O
      P   O
      Al  O  O
      Ti  O  O
      P   O  O


-   `num_iteration` -- Number of iterations (generations) to be
    computed
-   `print_level` -- Frequency of output \[default: `1`\]
-   `fitting_method` -- Optimization algorithm \[default: `cs`\]
-   `sample_directory` -- Directory where the reference data,
    `data.ref.XXX`, exist.
-   `param_file` -- Parameter file that contains initial values and
    ranges.
-   `match xxx yyy zzz` -- List of quantities used as optimization
    targets
-   `potential` -- Potential type whose parameters to be optimized.
    Currently available potentials are Morse, BVS, and BVSx.
-   `cs_XXXX` -- Parameters related to CS.
    -   `cs_num_individuals` -- Number of individuals (nests) in a
        generation.
    -   `cs_fraction` -- Fraction of abandons in a generation.
-   `update_vrange` --
-   `fval_upper_limit` -- Upper limit of loss function. The loss
    functions above this limit is set to this value.
-   `specorder` -- Order of species used in reference and MD program.
-   `interaction` -- Pairs and triples that are taken into account for
    optimization.

## Parameter file in.vars.fitpot

The parameter file `in.vars.fitpot` contains initial values and ranges
of each parameter to be explored. The file can be specified by
`param_file` in `in.fitpot` file.

    #  hard-limit:   T
    #
      10     6.000   3.000
         1.0000     1.0000     1.0000     1.000    1.000
         0.9858     0.5000     1.5000     0.500    3.000
         0.8000     0.5000     1.5000     0.500    3.000
         0.9160     0.5000     1.5000     0.500    3.000
         1.1822     0.5000     5.0000     0.100   10.000
         2.1302     1.5000     3.0000     0.100   10.000
         1.9400     1.5000     2.5000     0.100   10.000
         4.1963     3.0000     8.0000     0.100   10.000
         2.5823     1.5000     3.0000     0.100   10.000
         1.4407     1.2000     2.0000     0.100   10.000


-   Lines begin with `#` at the head of the file are treated as comment
    lines.
-   `hard-limit:  T` in comment line is a optional setting. The
    `hard-limit` set additional hard limit for parameters for automatic
    update of the search range.
-   1st line -- Number of optimizing parameters `NVAR`, cutoff radius
    for 2-body potential `RCUT2`, and cutoff for 3-body potential
    `RCUT3`, respectively.
-   2nd line and later -- initiall value, soft-limit (lower and upper),
    hard-limit (lower and upper), respectively. If `hard-limit: F`
    (hard-limit is not set), entries for hard-limit are not required in
    a line.

## Subjob script subjob.sh

The `subjob.sh` is used to perform MD runs and extract data for
evaluating the loss function of each nest (individual). :

    #!/bin/bash
    #=======================================================================
    #  Script to be called from fp.py to perfom pmd simulation
    #  and to extract RDF, ADF, and volume data.
    #
    #  Usage:
    #    $ run_pmds.sh
    #=======================================================================

    #...copy filed required for pmd calculation
    cp ../in.pmd* ../pmdini ./

    #...cd to the directory and clean up
    rm -f dump_* out.* data.pmd.*

    #...NpT MD
    cp in.pmd.NpT in.pmd
    pmd 2>&1 > out.pmd.NpT
    head -n166 out.pmd.NpT
    tail -n20 out.pmd.NpT
    echo "NpT-MD done at" `date`
    #...extract rdf, adf, vol and rename files
    python ~/src/nap/nappy/rdf.py -d 0.05 -r 5.0 --gsmear=2 --skip=80 --specorder=La,Li,F --pairs=La-F,Li-F --out4fp -o data.pmd.rdf dump_* 2>&1
    python ~/src/nap/nappy/adf.py --gsmear=2 --triplets=Li-F-F --out4fp --skip=80 -o data.pmd.adf dump_* 2>&1
    python ~/src/nap/nappy/vol_lat.py --out4fp --skip=80 dump_* 2>&1
    echo "post-processing done at" `date`


-   `--pairs` and `--triplets` should be correctly set in `rdf.py` and
    `adf.py` as well as `--specorder` options.
-   `--out4fp` option is required to write **whatever mode** format of
    reference data. On the other hand, in the case of **DF-matching
    mode**, `--out4fp` option should not be used.

## in.pmd file in the subjob


Here is an example of `in.pmd` file used in *subjob* of each individual
(nest), acually named `in.pmd.NpT` in `nap/examples/fp_LATP`.

    max_num_neighbors         200

    time_interval              2.0
    num_iteration            10000
    min_iteration               5
    num_out_energy           1000

    flag_out_pmd                1
    num_out_pmd               100
    flag_sort                   1

    force_type           Morse Coulomb angular
    cutoff_radius                6.0
    cutoff_buffer                0.3

    flag_damping                 0
    damping_coeff                0.99
    converge_eps                 1.0e-05
    converge_num                 3

    initial_temperature        300.0
    temperature_control        Langevin
    temperature_target         1  300.0
    temperature_relax_time     50.0
    remove_translation         1

    factor_direction          3 1
        1.00   1.00   1.00

    stress_control              vc-Berendsen
    pressure_target              0.0
    stress_relax_time           50.0


See [Input file: in.pmd](in-pmd.md) for detailed meaning of the
input file.

In short, this `in.pmd.NpT` is going to perform a MD simulation of
10,000 steps with Morse, Coulomb and angular potentials at 300 K under
NpT condition.

And from the output `pmd_###` files, target quantities are extracted
using some python scripts as described in `subjob.sh`. Those python
scripts create `data.pmd.XXX` files as output and *fp.py* is going to
read those data files to evaluate the loss function of each individual
(nest).

## Run fp.py


    $ python ~/src/nap/fitpot/fp.py --nproc 4 | tee out.fp


-   `--nproc` sets number of processes used for the evaluation of
    individuals.
-   `--subjob-script` option sets which script file is used for to
    perform subjob. [default: `subjob.sh`]
-   `--subdir` option sets the prefix of directories where the subjobs
    are performed. [default: `subdir`]

## Results and outputs


Files and directories created by *fp.py* are,

-   `out.fp` -- Standard output.
-   `out.cs.generations` -- Information of generations.
-   `out.cs.individuals` -- Information of all the individuals.
-   `in.vars.fitpot.####` -- Parameter file that is written whenever
    the best individual is updated.
-   `in.vars.fitpot.best` -- Parameter file of the best individual in
    the run.
-   `subdir_###` -- Directories used for the calculations of
    individuals. You can remove these directories after the run.


### Convert *fp.py* parameter file to *pmd* parameter files

    $ python ~/src/nap/nappy/fitpot/fp2prms.py BVSx in.vars.fitpot.best

This command will create `in.params.Morse`, `in.params.Coulomb` and
`in.params.angular` files (the keyword `BVSx` means that these 3
potentials).

### Visualize the evolution of optimization

One can plot loss function values of all the individuals appeared during
optimization as a function of generation using *gnuplot* as, :

    $ gnuplot
    gnuplot> set ylabel 'Loss function value'
    gnuplot> set xlabel 'Generation'
    gnuplot> p 'out.cs.generations' us 1:3 w p pt 5

-   Check if the loss function converges.
-   Check that the minimum loss function value is sufficiently small
    (below 0.01 per target would be good enough).

### Visualize the distribution of each parameters

You can plot the parameter values of all the individuals using the data
in `out.cs.individuals` as,

    $ gnuplot
    gnuplot> p 'out.cs.individuals' us 7:2 w p t 'D (Li-S)', '' us 8:2 w p t 'alpha (Li-S)', '' us 9:2 w p t 'Rmin (Li-S)'
# Examples

------------------------------------------------------------------------

## Structure relaxation


In general, the atom configuration made systematically using `mkconf_*`
may not be most stable structure and one needs to relax the structure.
To relax the strucure, you can run MD simulation with **FIRE algorithm**[^Bitzek2006]
as,

- without controlling temperature ( `temperature_control` to be
    `none`);
- with applying velocity damping ( `flag_damping` to be `2` to use
    **FIRE**);
- with convergence criterion 1.0e-4 eV ( `converge_eps` to be `1.0d-4`
    ).

Then atoms move towards the direction where the forces on atoms indicate
with reducing their energies, and finally the system becomes the least
energy structure.

After relaxing the structure with some outputs, use the relaxed
structure as the initial structure of the next simulation as,

    $ mv pmdini pmdorig     <== save initial structure to pmdini
    $ mv pmdfin pmdini      <== replace pmdini with pmdfin (relaxed structure)

------------------------------------------------------------------------

## Lattice constant and bulk modulus


In order to calculate lattice constant and bulk modulus, we have to
change lattice size and evaluate potential energy of the system. The
lattice size that minimizes the potential energy is the lattice
constant, and curvature around the lattice constant corresponds to the
bulk modulus.

1.  Prepare a small simulation cell of the system.

2.  Make `pmdini` file.

3.  Set `num_iteration` value in `in.pmd` file to zero. (Because only
    the 1st evaluation of the potential energy and forces are required.)

4.  Perform test run of *pmd* as,

        $ /path/to/pmd/pmd

    And confirm that *pmd* was done correctly.

5.  Run the script as,

        $ /path/to/nap/nappy/energy_vs_size.py 3.1  3.3

    where 3.1 and 3.3 are the min and max value of lattice size. Then
    you get the following output and graph.

        3.1000   804.3570    -455.5304887
        3.1100   812.1662    -457.0570473
        3.1200   820.0259    -458.4111792
        3.1300   827.9360    -459.6000610
        3.1400   835.8969    -460.6305856
        3.1500   843.9086    -461.5093734
        3.1600   851.9714    -462.2420295
        3.1700   860.0854    -462.8289708
        3.1800   868.2507    -463.2688466
        3.1900   876.4675    -463.5607352
        3.2000   884.7360    -463.7041313
        3.2100   893.0563    -463.6989250
        3.2200   901.4287    -463.5453811
        3.2300   909.8532    -463.2441193
        3.2400   918.3300    -462.7960993
        3.2500   926.8594    -462.2027522
        3.2600   935.4414    -461.4659414
        3.2700   944.0761    -460.5877459
        3.2800   952.7639    -459.5704416
        3.2900   961.5048    -458.4164932
        3.3000   970.2990    -457.1285461
        plsq= [   1.8861122     1.58686615  888.41535225 -463.71506645]
        =============================== RESULTS ================================
        Lattice constant =     3.2044 Ang.
        Cohesive energy  =     -8.587 eV
        Bulk modulus     =     302.16 GPa
        ================================ OUTPUT ================================
        * out.Ecoh-vs-size
        * log.Ecoh-vs-size.eps

    Then you can make a graph using some plotting program like `gnuplot`
    with the `out.energy_vs_size` file.

    ![image](./figs/Ecoh-vs-size.png)

------------------------------------------------------------------------

## Elastic constants

By applying cell deformations which correspond to the elastic constants
you want to calculate, you can obtain elastic constants by calculating
potential energy for each deformed structure. (Note that this script
works only for cubic systems.)

And here it is assumed that `lattice-constant`{.interpreted-text
role="ref"} is done, and the lattice constant is alreadly obtained.

1.  Set the lattice constant in `pmdini` file to the value obtained in
    `lattice-constant`{.interpreted-text role="ref"} .
2.  Run the script as follows, then you can get the following outputs
    and graph.

        $ calc-elastic-constants.py
             0.0000    -222.6114952    -222.6114952    -222.6114952
             0.0010    -222.6110355    -222.6111095    -222.6112676
             0.0020    -222.6096187    -222.6099522    -222.6105848
             0.0030    -222.6072423    -222.6080232    -222.6094468
             0.0040    -222.6039040    -222.6053225    -222.6078536
             0.0050    -222.5996012    -222.6018497    -222.6058052
             0.0060    -222.5943316    -222.5976044    -222.6033015
             0.0070    -222.5880928    -222.5925864    -222.6003426
             0.0080    -222.5808824    -222.5867951    -222.5969284
             0.0090    -222.5726978    -222.5802299    -222.5930588
             0.0100    -222.5635369    -222.5728903    -222.5887338
        =============================== RESULTS ================================
         C11     =    244.481 GPa
         C11-C12 =     98.392 GPa
         C12     =    146.089 GPa
         C44     =    116.030 GPa
         Following values maybe only valid for isotropic materials...
         Young's modulus =    215.743 GPa
         Poisson's ratio =      0.203
         shear modulus   =     89.698 GPa
        ================================ OUTPUT ================================
         * out.elastic-constants
         * graph.elastic-constants.eps

    ![image](./figs/graph_elastic-constants.png)

------------------------------------------------------------------------

## Phonon dispersion using phonopy

Phonon dispersion relation can be calculated using
[phonopy](http://phonopy.sourceforge.net) program. So you have to
install *phonopy* before moving forward in this topic.

First, prepare the atom configuration file for pmd `pmdini` which
contains cell structure and atom positions. Generally this should be the
primitive unit cell of the crystall structure you are considering now.

!!! Warning
    Following approach, making **FORCE\_CONSTANTS**, is an old fashion way.
    See [pmd2phonopy](#pmd2phonopy) , you can run only one command to get
    the phonon dispersion relation.


2nd, run the following command to get the **force constansts** via
finite displacement approach as,

    $ python /path/to/nap/nappy/force_constant.py -d 0.0001 -r 3.772 pmdini

You get `POSCAR` and `FORCE_CONSTANTS` to be used by `phonopy` program.
The option with `-d` means the magnitude of displacement in Angstrom,
and `-r` for cutoff radius of the interatomic potential used.
`force_constant.py` will show an output as following,

    displacement =  0.0001  Ang.
    rcut         =  3.772  Ang.
    POSCAR was written.
    vol of unit cell= 40.0456573564
    num of cells in each axis= 5 5 5
    num of atoms in extended system= 250
    sysext.num_atoms()= 250

It means that the `FORCE_CONSTANTS` file contains 5x5x5 cells of the
original primitive unit cell. This values will be passed to phonopy
below.

3rd, prepare a configuration file for `phonopy` (here it is named as
`conf.phonopy`).

    ATOM_NAME = Si
    CELL_FILENAME = POSCAR
    CREATE_DISPLACEMENTS = .FALSE.
    FORCE_CONSTANTS = READ
    DIM = 5 5 5
    BAND = 0 0 0  1/2 1/2 0  1/2 1/2 1/2  1/2 0 0  0 0 0

Here `DIM` should be the same as the values above. And
`FORCE_CONSTANTS = READ` let phonopy read force constants from the file
`FORCE_CONSTANTS`. Running phonopy with specifying this `conf.phonopy`
as input,

    $ phonopy -p conf.phonopy

you can get phonon dispersion graph as following.

![image](./figs/phonon-dispersion-Si.png)

If you specify the phonopy configuration file like,

    ATOM_NAME = Si
    CELL_FILENAME = POSCAR
    CREATE_DISPLACEMENTS = .FALSE.
    FORCE_CONSTANTS = READ
    DIM = 5 5 5
    MP = 5 5 5
    DOS_RANGE = 0 20 0.1
    SIGMA = 0.2

You can get a phonon DOS.

![image](./figs/phonon-dos-Si.png)

### pmd2phonopy.py to get the phonon dispersion directory

Above approach, in which `FORCE_CONSTANTS` are obtained, is old fashion
and lengthy. Now users can get the phonon dispersion relation directory
with one command.

First, prepare `pmdini` and `band.conf` files. `band.conf` file should
be like the following,

    ATOM_NAME = Si
    DIM =  4  4  4
    BAND = 0 0 0  1/2 0 1/2,  1/2 1/2 1  0 0 0  1/2 1/2 1/2

Here `ATOM_NAME` is necessary if you want to get correct frequency
values in phonon dispersion graph. Only phonon frequencies correspond to
the wave numbers given by `BAND` parameters are shown.

Users do not need to consider `DIM`, these values are automatically
determined in the following command.

``` {.bash}
$ python /path/to/pmd2phonopy.py -c 4.0 pmdini
```

Then you can get a band output file `out_band`, so you can see a graph
with `gnuplot` as,

``` {.bash}
$ gnuplot
gnuplot> plot "out_band" us 1:2 w l
```

------------------------------------------------------------------------

## Non-equilibrium molecular dynamics (NEMD) for heat flux

One can perform NEMD heat flux simulation applying different
temperatures at different regions. The figure below shows a setting of
heat-flux NEMD simulation.

![image](./figs/NEMD-setting.png)

Heat flux must be flown along *x*-direction. In this case, at both left
and right side of the system, atoms are fixed during simulation. And
vacuum region is placed in order to avoid interactions between hot and
cool atoms.

The digits at the bottom of above picture are *ifmv* values. The *ifmv*
values of fixed atoms should be 0, hot atoms to be 2, cool atoms 3,
intermediate atoms 1. And the temperature setting in `in.pmd` is like
following,

    initial_temperature     -600d0
    temperature_control     Berendsen
    temperature_target   1  -300d0
    temperature_target   2  350d0
    temperature_target   3  250d0
    temperature_relax_time  100d0
    factor_direction 3 3
      1.000d0  1.000d0  1.000d0
      1.000d0  1.000d0  1.000d0
      1.000d0  1.000d0  1.000d0  

Target temperatures of hot and cool atoms are 350 K and 250 K. Whereas
the target temperature of intermediate atoms are not set and the value
is set negative.

Temperature distribution along *x* is obtained by specifying as follows
in `in.pmd`,

    flag_temp_dist  T
    num_temp_dist  50

which means setting flag `.true.` and the number of bins along *x* is
50. The `num_temp_dist` value must be a multiple of `num_nodes_x`.
Results will be written in `out.temp-dist`. The `out.temp-dist` file
contains temperatures of 50 points of every `num_iteration`/
`num_out_energy` steps.




[^Bitzek2006]: E. Bitzek, P. Koskinen, F. Gähler, M. Moseler, and P. Gumbsch, Phys. Rev. Lett. 97, 170201 (2006).
# Publications

List of publications in which *nap* was used. Let us know if you want to show your publications in which *nap* was used.


-   [Kobayashi, R., Miyaji, Y., Nakano, K. & Nakayama, M.
    High-throughput production of force-fields for solid-state
    electrolyte materials. APL Materials 8,
    081111 (2020)](https://aip.scitation.org/doi/10.1063/5.0015373)
-   [Nakano, K. et al. Exhaustive and informatics-aided search for fast
    Li-ion conductor with NASICON-type structure using material
    simulation and Bayesian optimization Exhaustive and
    informatics-aided search for fast Li-ion conductor with NASICON-type
    structure using material. APL Materials 041112,
    041112 (2020)](https://aip.scitation.org/doi/10.1063/5.0007414)
-   [Kobayashi, R., Giofré, D., Junge, T., Ceriotti, M. & Curtin, W. A.
    Neural network potential for Al-Mg-Si alloys. Physical Review
    Materials 1,
    53604--53611 (2017)](https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.1.053604)
-   [Kobayashi, R., Hattori, T., Tamura, T. & Ogata, S. A molecular
    dynamics study on bubble growth in tungsten under helium
    irradiation. J. Nucl. Mater. 463,
    1071--1074 (2015)](https://www.sciencedirect.com/science/article/abs/pii/S002231151400991X?via%3Dihub)
-   [Kobayashi, R., Ohba, N., Tamura, T. & Ogata, S. A Monte Carlo study
    of host-material deformation effect on Li migration in graphite. J.
    Phys. Soc. Jpn.
    82, (2013)](https://journals.jps.jp/doi/10.7566/JPSJ.82.094603)
-   [Kobayashi, R., Nakamura, T. & Ogata, S. A Coupled Molecular
    Dynamics/Coarse-Grained-Particle Method for Dynamic Simulation of
    Crack Growth at Finite Temperatures. Mater. Trans. 52,
    1603--1610 (2011)](http://dx.doi.org/10.2320/matertrans.M2011116)
# How to analyze simulation results


------------------------------------------------------------------------

## Plot evolution of energies


Energies are basic values in MD. If you perform *NVE* simulation in
which no damping, temperature control, nor external forces, the total
energy has to be conserved whereas kinetic and potential energies
fluctuate. This is often used to check the validity of forces on atoms
when programming new potentials.

When one performs structure relaxation by damping velocities of atoms,
as MD steps increase the kinetic energy has to decrease to zero.

When one performs constant temperture simulation, kinetic energy should
be almost constant during the simulation.

Since total, kinetic, and potential energies are written in `out.erg`
file, users can plot energy evolution using *gnuplot* command as, :

    $ gnuplot
    gnuplot> plot 'out.erg' us 1:3 w l t 'total', 'out.erg' us 1:4 w l t 'kinetic', 'out.erg' us 1:5 w l t 'potential'

Or copy `util/gp.erg` script to the working directory and, :

    gnuplot> load 'gp.erg'

------------------------------------------------------------------------

## Visualize atom configuration

### Convert from pmd-format files

There is a conversion program that changes from pmd format to
visualization software format. When visualizing the atom configuration
by using [Ovito](https://www.ovito.org) , first convert the *pmd*-format
files, `pmd_###` , to *LAMMPS-dump* format files by doing the following,

    $ /path/to/nap/nappy/napsys.py convert --specorder=A,B,C pmd_#### dump_####

where `--specorder=A,B,C` specifies the order of species used in the
pmd-format file. Then you get `dump_####` file where `####` should be a
4-digit number. The [Ovito](https://www.ovito.org) can open the
LAMMPS-dump format file.

If there are sequential `pmd_####` files, one can convert those files by
using bash for-statement as :

    $ for f in pmd_[0-9]*; do /path/to/nap/nappy/napsys.py convert \
        --specorder=A,B,C $f `echo $f | sed 's/pmd/dump/'`; done

### Write LAMMPS-dump file directly

If you set `flag_out_pmd` as `2`, the *pmd* program writes atomic
configurations in *LAMMPS-dump* format with file names `dump_####`. So
you can visualize those files without any conversion using *Ovito* or
some other programs that can visualize *LAMMPS-dump* file.

------------------------------------------------------------------------

## Radial distribution function (RDF)

RDF is commonly used to investigate structure of the system, e.g., solid vs liquid, or crystal vs amorphous.

To get the RDF,

    $ python /path/to/nap/nappy/rdf.py [options] dump_0*

then, you get averaged RDF over atoms in `out.rdf`.

Given atom configuration files, `dump_####`, are read and average over
atoms in those files are taken.

Options are shown below,

    Options:
      -h, --help  Show this help message and exit.
      -d DR       Width of the bin. [default: 0.1]
      -r RMAX     Cutoff radius of radial distribution. [default: 5.0]
      --gsmear=SIGMA
                  Width of Gaussian smearing, zero means no smearing. [default: 0]
      -o OUT      Output file name. [default: out.rdf]
      --specorder=SPECORDER
                  Order of species separated by comma, like, --specorder=W,H. [default: None]
      --skip=NSKIP 
                  Skip first NSKIP steps from the statistics. [default: 0]
      --no-average
                  Not to take average over files.
      --no-normalize
                  Not to normalize by the density.
      --plot      Plot figures. [default: False]

The RDF of each pair of species normalized with the density of the pair.

![image](./figs/graph_rdf.png)

------------------------------------------------------------------------

## Angular distribution function (ADF)

To get ADF, perform `adf.py` something like,

    $ python /path/to/nap/nappy/adf.py --triplets=La-F-F,Ba-F-F dump_0*

The triplets consisting angles must be provided via the option
`--triplets`. Note that the 1st species in the triplet is the central
atom having bonds to the other two atoms, which maybe counter-intuitive.

![image](./figs/graph_adf.png)

Here is some options of `adf.py`,

    Options:
      -h, --help  Show this help message and exit.
      -w DEG      Width of the angular degree. [default: 1.0]
      -r RCUT     Cutoff radius of the bonding pair. [default: 3.0]
      --gsmear=SIGMA
                  Width of Gaussian smearing, zero means no smearing. [default: 0]
      --triplets=TRIPLETS
                  Triplets whose angles are to be computed. Three species should be specified connected by hyphen,
                  and separated by comma, e.g.) P-O-O,Li-O-O. [default: None]
      -o OUT      Output file name [default: out.adf]
      --skip=NSKIP 
                  Skip first NSKIP steps from the statistics. [default: 0]
      --no-average
                  Not to take average over files.
      --plot      Plot figures. [default: False]

------------------------------------------------------------------------

## Voronoi analysis

The other way to analyze the local structure is **Voronoi analysis**
which makes a cell which contains one atom according to the rule; the
minimum cell made of faces that bisect atom bonds.

To perform Voronoi analysis, you need to instal `voro++` first, and then

    $ python /path/to/nap/nappy/voro.py pmd_####

This command will provide an output `pmd_####.voro` and
`pmd_####.voro.vol` which are input and output of `voro++`,
respectively.

------------------------------------------------------------------------

## Velocity autocorrelation and power spectrum

!!! Note
    The file format used in this section, **akr**, is no longer used since
    years before. So the explanation here would not work directly. But, we
    do not remove this section because somebody might be interested in the
    procedure used here...


In order to get power spectrum from the MD simulation result, firstly we
have to think how long the MD simulation has to be run. In case of Si,
its phonon DOS exists up to about 16\~18 THz which is the inverse of
time interval of sampling data. And the frequency resolution is the
inverse of simulation time. So the time interval of sampling data should
be about 20 fs (which corresponds to 25 THz since the half of data will
be omitted because of the symmetry.) And the simulation time should be
10,000 fs which corresponds to the frequency resolution 0.1 THz.
Usually, one has to make about **1,000 akr files** for the power
spectrum calculation.

To get the velocity autocorrelation and power specturm, you can use
`power_spectrum.py` in `nappy` directory.

    $ python /path/to/nap/nappy/power_spectrum.py -t 20.0 --relax 5000.0 akr0???

Here `-t` option specifies the time interval between successive akr
files. `--relax` specifies relaxation time of the decaying factor for
autocorrelation function, if this is omitted no decaying factor is
applied.. The you get `dat.autocorr` and `dat.power` files.
`dat.autocorr` includes velocity autocorrelation functions of *x*, *y*,
*z*, and sum of those. `dat.power` also includes power spectrums of *x*,
*y*, *z*, and sum of those.

If this power spectrum graph seems too spiky, you can smear it by using
`gaussian_smear.py` as,

    $ python /path/to/nap/nappy/gaussian_smear.py -x 1 -y 5 -s 2.0 dat.power

Then you get `dat.power.smeared` file which contains only 2 columns of
blurred data of 1st and 5th columns of `dat.power`.
# nappy -- NAP PYthon utilities

*nappy* is a set of python utilities. It is contained in `nap/nappy/`
directory, but currenly it is not exactly a python module. Users may
have to use those utilities by calling them directly from the shell.

## Setup

*nappy* should work with python-2.7 and python-3 series. But, as the
maintenance of python-2.x series have stopped, it is recommended to use
python-3.

*nappy* requires the following packages, which can be installed using
`pip` command,

- [docopt](http://docopt.org)
- [ASE](https://wiki.fysik.dtu.dk/ase/index.html)

To use *nappy* in python program, it is required to add a path to
`nap/nappy` directory to the environment variable `PYTHONPATH`. In case
of `bash`, you can achieve this by adding the following line to
`~/.bash_profile`,

``` {.bash}
export PYTHONPATH=${PYTHONPATH}:/path/to/nap
```

You can check whether the path to `nappy` is added to `PYTHONPATH` by
the following command, :

    $ python -c 'import nappy; print(nappy.__file__)'


-----

## Quick start


Once *nappy* is installed, do the following code on *ipython* or
*jupyter notebook*,

```python
from nappy.napsys import NAPSystem, analyze
nsys = NAPSystem(fname='/path/to/nap/example/test_W/pmdini')
analyze(nsys)
```

Above code will show the following result. :

    a1 vector = [    10.254,      0.000,      0.000]
    a2 vector = [     0.000,      9.613,      0.000]
    a3 vector = [     0.000,      0.000,      9.613]
    a =     10.254 A
    b =      9.613 A
    c =      9.613 A
    alpha =   90.00 deg.
    beta  =   90.00 deg.
    gamma =   90.00 deg.
    volume=    947.617 A^3
    number of atoms   =  54


---

## Read and write files

To read atomic structures from files of several formats,

```python
nsys = NAPSystem(fname='POSCAR')
nsys = NAPSystem(fname='structure.dat', ffmt='xsf', specorder=['Al','O'])
```

Current available formats are:

- `pmd`: file format specific for `pmd` program
- `POSCAR`: VASP POSCAR file
- `dump`: LAMMPS dump file
- `xsf`
- `akr`: file format for Akira viewer
- other formats readable by [ASE](https://wiki.fysik.dtu.dk/ase/) .

And to write the structure to a file,

```python
nsys.write('POSCAR')
nsys.write('pmdini')
nsys.write_dump()
```

The `write` function will parse file name and choose appropriate file
format from the name.
# Development

This package, **nap** , is an open source project. And you can modify
and use for free of charge and by your own risk.

Usually users have to modify program source codes when they apply the MD
program to specific subjects they are tackling. The source codes in this
package contains a lot of lines, but the algorithms usesd in this
package are mostly basic and users can easily understand most of these
if they are familiar with the molecular dynamics.

*Enjoy coding and your MD simulation ;)*

------------------------------------------------------------------------

## Formulas

For those who want to read and modify program sources, it would be
helpful to show some details and formulas of what the **pmd** is doing
in the program. Basically the **pmd** performs standard MD with the
velocity Verlet algorithm, one should firstly understand the standard MD
formulas using some standard textbooks[^Frenkel].

However, of course, there are a lot code-specific definitions to
implement the standard MD. This section shows some of the code-specific
definitions used in the **pmd**.

------------------------------------------------------------------------

### Simulation cell

The simulation cell consists of one scalar value, $l$, and three
vectors, $\mathbf{a}, \mathbf{b}, \mathbf{c}$, and is shown in the
atom-configuration file as, :

    5.472                     <--- l
    0.500    0.500   0.000    <--- ax, ay, az
    0.500    0.000   0.500    <--- bx, by, bz
    0.000    0.500   0.500    <--- cx, cy, cz

The cell information actually used in the code is the matrix,
$\mathbf{h}$, defined as,

$$\begin{aligned}
\begin{equation}
   \mathbf{h} = l \times (\mathbf{a},\mathbf{b},\mathbf{c})= l \times
   \begin{pmatrix}
   a_x & b_x & c_x \\
   a_y & b_y & c_y \\
   a_z & b_z & c_z
   \end{pmatrix}.
\end{equation}
\end{aligned}$$

------------------------------------------------------------------------

### Positions

Atom positions, which are defined as `ra(1:3,1:natm)` in the code, are
normalized within 0.0 and 1.0 during the MD simulation so that they become
absolute positions after multiplying the cell matrix, $\mathbf{h}$, as,

$$\begin{aligned}
\begin{eqnarray*}
\mathbf{r}' & = & \mathbf{h} \cdot \mathbf{r}, \\
            & = & r_a\mathbf{a} +r_b\mathbf{b} +r_c\mathbf{c},
\end{eqnarray*}
\end{aligned}$$

which is written in the Fortran code as follows,

```fortran
xi(1:3)= h(1:3,1)*ra(1,i) +h(1:3,2)*ra(2,i) +h(1:3,3)*ra(3,i)
```

------------------------------------------------------------------------

### Velocities

Atom velocities, `va(1:3,1:natm)`, are also scaled by the cell matrix.
However, not only that, but also scale by time, which means the
velocities in the code have actually length scale but velocity scale.
Thus, in the velocity Verlet algorithm, velocities are directly added to
the positions without multiplying time interval, $\Delta t$.

------------------------------------------------------------------------

### Accelerations

Atom accelerations, `aa(1:3,1:natm)`, are also scaled by the cell matrix
and are multplied by $\Delta t^2$ so that they become the units of
positions at the end of the force calculation subroutines. Thus, the
accelerations as well can be directly added to the positions without
multiplying $\Delta t^2$ in the velocity Verlet loop.

------------------------------------------------------------------------

### Velocity Verlet {#vv}

Basically MD with the velocity Verlet algorithm is very simple:

1.  Compute initial forces.
2.  Velocity Verlet loop starts.
    1.  Update velocities using the current forces with a half of time
        interval, $\Delta t/2$.
    2.  Update positions using the current velocities with a time
        interval, $\Delta t$.
    3.  Compute forces from the current positions.
    4.  Update velocities using the current forces with a half of time
        interval, $\Delta t/2$.

This is written mathematically as,

$$\begin{aligned}
\begin{eqnarray*}
\mathbf{v}_i^{(n+1)} &=& \mathbf{v}_i^* +\frac{\mathbf{f}_i^{(n)}}{m_i} \frac{\Delta t}{2}, \\
\mathbf{r}_i^{(n+1)} &=& \mathbf{r}_i^{(n)} +\mathbf{v}_i^{(n+1)}\Delta t, \\
\text{Compute}\ &\ & \mathbf{f}_i^{(n+1)}\left(\left\{\mathbf{r}^{(n+1)}\right\}\right), \\
\mathbf{v}_i^* &=& \mathbf{v}_i^{(n+1)} +\frac{\mathbf{f}_i^{(n+1)}}{m_i} \frac{\Delta t}{2},
\end{eqnarray*}
\end{aligned}$$

where subscript *i* is an atomic index and superscript *n* is an
MD-step.

This is implemented in the code as,

```Fortran
call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
&     ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
&     ,mpi_md_world,myid_md,epi,epot0,nismax,acon,avol
&     ,cforce)
...
do istp=1,nstp
   ...
   va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)
   ...
   ra(1:3,1:natm)=ra(1:3,1:natm) +va(1:3,1:natm)
   ...
   call get_force(namax,natm,tag,ra,nnmax,aa,strs,h,hi
&         ,tcom,nb,nbmax,lsb,lsrc,myparity,nn,sv,rc,lspr
&         ,mpi_md_world,myid_md,epi,epot,nismax,acon,avol
&         ,cforce)
   ...
   va(1:3,1:natm)=va(1:3,1:natm) +aa(1:3,1:natm)
   ...
enddo
```

As you can see, the basic construction is very simple. However, the
codes for miscellaneous stuff such as thermostat, isobaric, and
parallelization hidden in the above code-block as `...` are rather
lengthy.

------------------------------------------------------------------------

### Berendsen barostat

Berendsen barostat is similar to Berendsen thermostat, which control the
stress or temperature moderately to the target ones.

------------------------------------------------------------------------

## Force implementation


Implementation of the interatomic force and potential is the core of MD
programming. And if you want to perform simulation that includes a
combination of elements that is not implemented in the **pmd**, you have
to implement the force calculation routine by yourself. Every force
routine is defined in a separated module such as `force_SW_Si.F90`,
which indicates the force routine of Stillinger-Weber type for Si
system, and is called via `get_force` subroutine in `force_common.F`
file. Thus you can write your own force routine by following `get_force`
and `force_SW_Si` subroutines.

The **force-type** used for the simulation is determined by a variable,
`cforce`, which is a character and specified in *in.pmd* file. The
corresponding force routine is called according to the **force-type**
name. Each force subroutine has to have the same arguments and order and
should be provided using `use` in the beginning of the `get_force`
routine.

------------------------------------------------------------------------

## ASE interface

There is a python script that connects pmd to **ASE (atomistic
simulation environment)**. This enable us small calculation of pmd much
easier. The following code shows how to use
`nappy/interface/ase/pmdrun.py` with ase program.

``` {.python}
import os,sys
from ase.io import read
from nappy.interface.ase.pmdrun import PMD

atoms=read('POSCAR',index=0,format='vasp')
os.system('cp /path/to/in.*.NN ./')
calc= PMD(label='pmd',command='/path/to/nap/pmd/pmd > out.pmd',force_type='NN')
atoms.set_calculator(calc)
print atoms.get_potential_energy()
print atoms.get_forces()
```

When `atoms.get_potential_energy()` is called, pmd program is performed
background and ASE gets calculation results from `out.pmd` file.

------------------------------------------------------------------------

## Versioning and tagging

The **nap** uses Git for version controlling. Developers are recommended
to make new branches to modify codes, and merge the change to the master
branch. Tags are named like `rev170605` which means the *revision* of
the date 2017-06-05. If you need to make new tag at the same date, which
would not happen often, you can name the tags like `rev170605_1` or
something like this by adding some suffix after an underscore.

Versioning in **nap** follows (not strictly) the [semantic
versioning](https://semver.org). However, since API is not defined, the
major version will never get one and the version will be like `v0.x,x`.

[^Frenkel]: Dann Frenkel and Berend Smit, *Understanding Molecular Simulation*, Academic Press
# Cluster Manager utility

## What's `clmgr.py`?

This python utility automates the process of submitting jobs on clusters and supercomputers such as

* to make job scripts
* to submit jobs with specifying appropriate resources
* to combine some jobs to submit in one job script to reduce the actual number of submissions


### What is not provided by `clmgr.py`?

* to make input files for calculations
* to send input files from local computer to remote clusters
* to gather results from remote clusters
* to parse and check calculation results

There are some very sophisticated programs that can do these processes as well as the processes that `clmgr.py` provides, such as AiiDA by MARVEL and fireworks by materials project, but to remember how to use these very sophisticated program is also not so easy. Compared to those programs, `clmgr.py` is thought to be much easier to setup and use.


### Why is it included in nappy?

Basically this `clmgr.py` can be provided as a separated module, but it requires some parsers such as VASP parser which is already implemented in nappy and it is just easier for me to make this by utilizing the existing package. And also this will make installation process easier as well.


## Install nappy

```
export PYTHONPATH=/path/to/nap:$PYTHONPATH
```

Then you can test if the installation went well by
```
$ python -c "import nappy; print nappy.__file__"
```
If you get the path to nappy directory, the installation is OK.

## Setup on remote cluster machine

Following files are required to run `/path/to/nappy/clutil/clmgr.py`,
probably on a remote cluster or supercomputer.

* `~/.nappy/clmgr`
* `~/.nappy/clmgr/machine`
* `~/.nappy/vasp.conf`


### `machine` file

The machine file should be written in YAML format and contain information about:

* scheduler-type to be used
* list of queues and corresponding resources
* number of processors/cores in one node
* MPI command to be used

An example of `machine` file for Fujitsu FX100 is as follows,
```
scheduler: fujitsu
queues:
    'fx-debug':  {'num_nodes':  32, 'default_sec':  3600, 'limit_sec':   3600}
    'fx-small':  {'num_nodes':  16, 'default_sec': 86400, 'limit_sec': 604800}
    'fx-middle': {'num_nodes':  96, 'default_sec': 86400, 'limit_sec': 259200}
    'fx-large':  {'num_nodes': 192, 'default_sec': 86400, 'limit_sec': 259200}
    'fx-xlarge': {'num_nodes': 864, 'default_sec': 86400, 'limit_sec':  86400}
nprocs_per_node: 32
mpi_command: "mpiexec --vcoordfile {rankfile} -n {npara} --stdout {out} {exec_path}"
```

#### scheduler

`scheduler` can take values either one of the following,

* `fujitsu`
* `pbs` or `torque`


#### queue

`queues` should be a list of queue names and each queue name should have a dictionary with entries, `num_nodes`, `default_sec`, and `limit_sec`.

#### nprocs_per_node

This value is very machine specific and you must provide it.

#### mpi_command

The `mpi_command` entry specifies how the MPI command is used in the machine.
Currently following parameters can be used in the command,

* `{npara}`: number of parallel processes for the MPI run
* `{out}`: output file path
* `{exec_path}`: executable file path
* `{rankfile}`: path to the rankfile which is used in Fujitsu FX100


### `vasp.conf` file

This file should be written in YAML format and should contains the path to the *VASP* executable file like,

```
exec_path: '~/bin/vasp541fx'
```


## Usage

The usage can be shown by typing `clmgr.py -h`,
```
Usage:
  clmgr.py [options] DIRS [DIRS..]

Options:
  -h, --help  Show this message and exit.
  -c CALCULATOR
              Set calculator type. Available types are: vasp
              [default: vasp]
  -d          Dryrun: run clmgr without actually submitting jobs.
              [default: False]
  -m          Enable multiple jobs in one submission.
              In case that there is a limitation to the number of jobs 
              that can run at the same time in the server, you may 
              had better set this TRUE.
              [default: False]
  -q QUEUE_NAME
              Queue name the jobs are submitted to. [default: default]
```

So, for example, in Linux cluster without any limitation of running jobs per user, one can run `clmgr.py` like,
```
$ /path/to/clmgr.py -q queue_name calc_dir
```

Or in some supercomputer which limits the number of running jobs at a time,
```
$ /path/to/clmgr.py -q queue_name -m calc_dir
```

### log file

`clmgr.py` creates a log file at `~/.nappy/clmgr/` with the name including *PID* of the process that `clmgr.py` was running. In the file, you can find which directories are treated as directories where jobs will be running and to which job id the jobs are assigned.

# fitpot -- fit parameters of neural-network potential

The validity of MD simulations depends strongly on the accuracy of the
interatomic potential used in the simulation. So, when you think of
doing some simulation of specific system, you have to prepare an
interatomic potential that reproduces the phenomena you are thinking of.

Here we indroduce how to make a neural-network potential and fit
potential parameters with the *fitpot* program included in **nap**
package.

!!! Note
    Currently, the *fitpot* program is used only for neural-network
    potential. For other classical potentials, use *fp.py* instead.


## What does fitpot do?

In the *fitpot*, the following loss function is minimized by optimizing
potential parameters $\{ w \}$.

$$\mathcal{L}(\{w\}) = \frac{1}{\eta N_s}\sum_s^{N_s} \left[ \Delta E^2 +\sum_i^{N^s_\mathrm{a}}\left| \Delta \boldsymbol{F}_i\right|^2 +\left| \Delta \sigma \right|^2\right]$$

in the case of fitting energies and forces. Here, $s$ is the sample
number, $N_s$ the number of samples, $N^s_\mathrm{a}$ the number of
atoms in the sample $s$, $\eta$ the parameter corresponding to how many
properties are fitted (in this case $\eta = 3$ because energy, force and
stress are used.)

To minimize the above loss function, the following gradient-based
methods are available in *fitpot*:

-   Steepest descent (SD)
-   Quasi-Newton method (BFGS)
-   Conjugate gradient (CG)
-   Stochastic gradient descent (SGD)

## Compilation

Since some modules in *pmd* program are required for the compilation of
*fitpot*, compile *pmd* before compiling *fitpot*. :

    $ cd /path/to/nap/
    $ ./configure --prefix=$(pwd)
    $ cd pmd
    $ make pmd
    $ cd ../fitpot/
    $ make fitpot

------------------------------------------------------------------------

## Quick trial with an example

There is an example of *fitpot* with minimal dataset to see how it
works. Go to the directory `examples/fitpot_DNN_SiO/`, read `README.md`,
try running *fitpot*, and look at some output files.

------------------------------------------------------------------------

## Fitting procedure

Hereafter, we assume that the reference data are obtained by using an
*ab-initio* calculation program, VASP.

Potential parameters are fitted as the following procedure:

1. [vasp-data](#prepare-reference-data)
2. [prepare-inputs](#prepare-input-files)
3. [exec-fitpot](#run-fitpot-program)

------------------------------------------------------------------------

### Prepare reference data

Assuming that there are some reference data in `dataset/` directory, and
all the data are stored in the directories whose names start with
`smpl_`.

The following files are required in eacy sample directory (`smpl_*`):

-   `pos`
-   `erg.ref`
-   `frc.ref`
-   `strs.ref`

`pos` is a pmd-format atom-configuration file, `erg.erg` contains a
scalar value of energy of the system, `frc.ref` contains the number of
atoms in the system and forces of all atoms shown as, :

    4
    0.1000   0.0000   0.0000
    0.0000   0.1000   0.0000
    0.0000  -0.1000   0.0000
    -0.1000   0.0000   0.0000

In the case of extracting DFT data from *ab-initio* MD runs with VASP,
positions, energy, forces and stress of each MD step can be obtained
from `vasprun.xml` file as follows, :

    $ python path/to/nap/nappy/vasp/vasprun2fp.py /path/to/dir/that/includes/vasprun.xml/

Then you get directories with names like `#####` including `pos`,
`erg.ref`, `frc.ref` and `strs.ref` files in them.

### Prepare input files

Inputs files needed for *fitpot* are the following:

-   `in.fitpot`
-   `in.params.DNN`
-   `in.params.desc`
-   `in.params.Coulomb` in each `smpl_XXX` directory in some special
     cases

You have to specify the `num_samples` in `in.fitpot` file which is the
number of samples in `dataset/` directory. The number of sample
directories can be counted by the following command,

```bash
$ ls /path/to/dataset | grep smpl_ -c
```

### Run fitpot program

In the directory where `dataset/` directory and `in.fitpot` file exist,
you can run the *fitpot* program as, :

    $ ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

Or if you want it to run in parallel mode, :

    $ mpirun -np 10 ~/src/nap/fitpot/fitpot > out.fitpot 2>&1 &

There are some output files:

- `out.erg.trn.fin`, `out.erg.tst.fin` -- These files include reference and *pmd* data of energies. To see whether the fitting went well or not, plot these data by using `gnuplot` as,

        $ gnuplot
        gnuplot> plot 'out.erg.trn.fin' us 1:2 w p t 'training set'
        gnuplot> rep 'out.erg.tst.fin' us 1:2 w p t 'test set'

- `out.frc.trn.fin`, `out.frc.tst.fin` --  These files include reference and *pmd* data of forces.

------------------------------------------------------------------------

Input file for *fitpot*
-----------------------

The following code shows an example of the input file `in.fitpot`.

    num_samples       14
    test_ratio        0.1
    num_iteration     100
    num_iter_eval     1

    fitting_method    bfgs
    sample_directory  "../dataset"
    param_file        in.params.DNN
    normalize_input   none

    energy_match       T
    force_match        T
    stress_match       T
    potential          DNN

    ftol              1.0e-5
    xtol              1.0e-4

    penalty           none
    penalty_weight    1d-3

    # Species order:  1) Al, 2) Mg, 3) Si
    specorder    Al  Mg  Si


------------------------------------------------------------------------

### num_samples

Default: *none*

Number of reference samples to be used for training and test.

------------------------------------------------------------------------

### sample_list

Default: `none`

Path to the file that contains a list of samples to be used for training
and test. The format of the list file should be like, :

    smpl_001
    smpl_002
    smpl_003
    ...

or with specifying which samples are training (`1`) or test (`2`) as, :

    smpl_001  1
    smpl_002  2
    smpl_003  1
    ...

If whether training or test is specified in the list,
[test_ratio](#test_ratio) will be neglected.

------------------------------------------------------------------------

### test_ratio

Default: `0.1`

The ratio of test data set $r$ within whole data set $N$. Thus the
number of test data set is $rN$, and the number of training data set is
$(1-r)N$.

------------------------------------------------------------------------

### num_iteration

Default: `1`

Number of iterations of a minimization method.

------------------------------------------------------------------------

### num_iter_eval

Default: `1`

Test data set will be evaluated every *num_iter_eval* iterations.


------------------------------------------------------------------------

### fitting_method

Default: *test*

The method used to fit parameters to the sample data. Available methods
are the following:

- `sd`/`SD` -- Steepest descent algorithm which requires gradient information.
- `cg`/`CG` -- Conjugate gradient algorithm which requires gradient information.
- `bfgs`/`BFGS` -- Quasi-Newton method with BFGS. This requires gradient information.
- `check_grad` -- Comparison of analytical derivative and numerical derivative. Use this to check the implemented analytical gradient.
- `test`/`TEST` -- Just calculate function L and gradient of L w.r.t. fitting parameters.


------------------------------------------------------------------------

### sample_directory

Default: `dataset`

The directory that includes sample data. We call this `dataset` in the
above instruction.

If you want to use `..` to specify the directory relative to the current
working directory, e.g. `../dataset`, you need to enclose with
double-quotation marks like `"../dataset"`.

------------------------------------------------------------------------

### param_file

Default: *in.params.DNN*

The name of the file that has parameter values in it. This is passed to
`pmd` program.

------------------------------------------------------------------------

### ftol

Default: *1.0e-5*

The tolerance of difference of the loss function value.

------------------------------------------------------------------------

### xtol

Default: *1.0e-4*

The tolerance of the change of variables which are optimized. If either
one of [ftol]{.title-ref} or [xtol]{.title-ref} is achieved, the
optimization stops.

------------------------------------------------------------------------

### energy_match, force_match, stress_match

Default: *True* for energy, *False* for force and stress

Whether or not to match forces. ( *True* or *False* ) It is recommended
to match not only energy but also forces, since forces are important for
molecular dynamics.

------------------------------------------------------------------------

### potential or force_field

Default: `DNN`

The potential whose parameters you are going to fit.
Potentials currently available are:

- `DNN` -- Neural-network potential
- `linreg` -- Linear regression potential

------------------------------------------------------------------------

### random_seed

Default: *12345d0*

Initial random seed for the uniform random numbers used in the *fitpot*.
This is used to change the random choice of training and test sets.

------------------------------------------------------------------------

### regularize

Whether or not regularize bases obtained in *linreg* and *DNN*
potentials. ( *True* or *False* )

Default: *False*

------------------------------------------------------------------------

### penalty

Default: *no*

Type of penalty term, *lasso* which is L1-norm penalty or *ridge* which
is L2-norm penalty, or *no* which means no penalty term.


------------------------------------------------------------------------

### penalty_weight

Default: *1.0*

The weight applied to the penalty term. This value also has to be
determined through cross-validation scoring\...

------------------------------------------------------------------------

### sample_error

Default: *0*

The number of samples whose errors are to be given. These errors appear
at the denominators of energy and force in the evaluation function such
that

$$\left( \frac{E^\mathrm{NN}-E^\mathrm{DFT}}{N_\mathrm{a}\varepsilon_\mathrm{e}}\right)^2 +\sum_i^{N_\mathrm{a}} \sum_\alpha^{xyz} \frac{1}{3N_\mathrm{a}}\left( \frac{F^\mathrm{NN}_{i\alpha} -F^\mathrm{DFT}_{i\alpha}}{\varepsilon_\mathrm{f}}\right)^2$$

If the difference between NN energy and DFT energy/force is lower than
this value, this term becomes less than 1.0, which means the
energy/force of the sample is thought to be converged. The initial
values of the errors are 0.001 (eV/atom) and 0.1 (eV/Ang) for energy and
force, respectively.

There must be the same number of following entry lines as the above
value which determine the errors of energy and force of each sample like
the this, :

    sample_error   2
        Al_fcc    0.001  0.2  1.0
        Al_bcc    0.001  0.2  1.0

The each entry has *entry_name*, *error of energy (eV/atom)*, *error of
forces (eV/Ang)* and *error of stresses (GPa)*. The error values are
applied to all the samples that contain *entry_name* in their directory
names.

------------------------------------------------------------------------

### force_denom_type

`relative` or `absolute`

Default: `relative`

Which type of denominator of force term in the loss function is used. If
`absolute` is specified, the *fitpot* uses an *error of forces*
specified in the [sample_error](#sample_error) for the
denominator of force term. If `relative` is specified, the *fitpot* uses
a magnitude of force on the atom in the denominator of force term.

------------------------------------------------------------------------

### specorder

Default: *none*

The order of species common in fitpot. This must be specified before
`atom_energy` entry and must hold for every samples.

------------------------------------------------------------------------

### init_params

Default: `read`

Whether the paramters to be optimized are read from the file or
initialized.

- `read` -- Read parameters from the file.
- `gaussian` --  Parameters are initialized with Gaussian distribution according to *init_params_sgm* and *init_params_mu*.

------------------------------------------------------------------------

### init_params_sgm

Default: `1d0`

Variance of Gaussian distribution of the initial values for parameters.

------------------------------------------------------------------------

### init_params_mu

Default: `0d0`

Mean value of Gaussian distribution of the initial values for
parameters.

------------------------------------------------------------------------

### init_params_rs

Default: `12345.0`

Random seed for the initialization of parameters. This random seed is
only used for this purpose and does not affect random seed for the
choice of training and test sets, which is affected by
`random_seed`{.interpreted-text role="ref"}.

分子動力学法(MD)の基礎
======================

MDの概要
--------

分子動力学法(molecular dynamics,
MD)は，原子・分子の動きを計算機上で再現することで，

> -   実験では観測できないミクロの世界の現象を観たり，
> -   実験では行うことのできない状況の仮想実験を行ってみたり，
> -   統計データからマクロな物理量を計算したり，

と，様々なことが可能な計算機シミュレーション手法である．

単純に言えば，MDシミュレーションとは，

> 1.  原子間に働く力，つまりは原子に働く加速度を求める．
> 2.  加速度に従い原子を動かす．

の2つを繰り返すシミュレーションであり，
原子間に働く力をいかに正確に求められるかがシミュレーションの精度を決定する．
単純であるが(もしくは単純であるが故に?)，応用範囲は広く，研究レベルでは

> -   工学
> -   化学
> -   生物，創薬
> -   物理

の様々な分野で使われている．
しかし，有限要素法や第一原理計算に比べると，企業レベルではまだまだ知名度は低いと思われる．

分子動力学シミュレーションにおいて重要なことは，

> -   原子間ポテンシャル
> -   モデリング
> -   結果の解釈

である．すでに述べたように，まずは原子間ポテンシャルが正しくないと，原子間に働く力が現実と異なるので，どんなに他を頑張ってもダメなことがある．
しかし，完璧な原子間ポテンシャルなどは存在しないので，その限界や適用範囲などを理解しながらシミュレーションを行い，
結果を解釈できれば，目的は達成できると期待される．
そのため，シミュレーションをするに際し，何を解決するために，どのようなことが分かれば良いのかを明確にしなければならない．
そして初めてどのようなシミュレーションを行えば良いのかが決まる．
この段階が **モデリング** と呼ばれる．

原子間ポテンシャル
------------------

原子に働く力は，原子の周りに存在している電子を介した量子論的な相互作用により初めて求まるものである．
そのため，量子論的な計算を行い，原子間に働く力を求め，それを用いて分子動力学シミュレーションを行う
**第一原理分子動力学法** もある．
ただし，一般的に量子論的な計算は非常に計算量が多いために，多数の原子や，長時間のダイナミクスを行うことは現在でも難しい．
そこで，量子論的な相互作用を何らかの関数系で近似し，原子に働く力を求めてシミュレーションを行うことを(古典)分子動力学法という．
これまでに，様々な関数系が種々の系に対して提案されている．

### Lennard-Joneポテンシャル

希ガス原子の間の相互作用は電子の揺らぎに起因する静電相互作用であり，van
der Waals相互作用と呼ばれる． これは距離 $r$
の-6乗に比例する引力となることが知られている．
この相互作用を表現するのがLennard-Jonesポテンシャルであり，次のような形となっている．

$$\phi (r_{ij}) = \epsilon_{ij} \left[ \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{12} - \left( \frac{\sigma_{ij}}{r_{ij}} \right)^{6} \right]$$

斥力項には物理的な根拠はないが，一般的には-12乗が用いられる．

### Moorseポテンシャル

### Embeded Atom Method (EAM) ポテンシャル

-   EAMポテンシャル
-   FSポテンシャル
-   MEAMポテンシャル
-   GEAMポテンシャル

### Stillinger-Weberポテンシャル

### Tersoffポテンシャル

### 長距離力 (クーロン力)

高速化技法
----------

### Cell list法とBook keeping法

原子間ポテンシャルの計算には近くに存在する原子の情報が必要である．
大抵の原子間ポテンシャルにおいて，隣接原子を検索する段階が計算スピードのボトルネックとなる．
そのため，MDシミュレーションでは，隣接原子の検索をどれだけ速く計算できるかが重要なポイントである．

単純に全ての他の原子との距離を計算すると，N個の原子一つ一つが(N-1)個の他の原子との距離を計算するので，
O(N\^2)の計算量となる． これを減らすための手法が **Cell list法** と
**Book keeping法** である．

**Cell list法**

:   系を小さなセルに分割し，そこに原子を登録していく．
    こうすることで，全ての原子との距離を計算する必要がなくなり，隣のセルに入っている原子とのみ距離を測れば良くなり，
    O(N)の計算量となる．

**Book keeping法**

:   Cell
    list法では，27個のセルの中に入っている原子の数だけ検索すれば良いので，O(N\^2)に比べれば非常に高速だが，
    そのセルの端にある原子はカットオフ半径の外にあるので，毎回検索するとムダが多い．
    そのため，カットオフ半径に入っている原子だけをリストに保存しておいて，ポテンシャルと力の計算に用いることで計算量を減らすことができる．
    その際，カットオフ半径よりも少し大きめの半径を取っておけば，毎ステップこのリストを更新する必要がなくなるため，
    全シミュレーションとしての計算量を減らすことができる．

### 空間分割並列化法

多くの原子間ポテンシャルは短距離相互作用となっているので，遠くまで相互作用は及ばない．
そのことを活かし，大規模な系を空間分割し，各領域を各プロセスとして同時並列的に各領域の計算を行う．
こうすることで，およぼ並列数倍だけ計算が速くなる．

各領域の端にある原子は，隣りの領域（つまりは別プロセス）の原子と相互作用する必要があるため，
毎ステップ原子座標情報を隣りの領域と交換する必要がある．
このデータ通信量が大きいと並列性能がでない．
隣りの領域と情報交換する原子は，領域の端からカットオフ半径程度にある原子であるため，
カットオフ半径が長いポテンシャルは，隣接原子探索にも時間がかかるし，並列性能も悪くなる．

MDで得られる物理量
------------------

### 凝集エネルギー

### 輸送係数

-   グリーン・久保公式
Notes
=====

-   Temperature control with simple velocity scaling is added. So the
    format in input file \'in.pmd\' is changed a little.
-   van der Waals potential for Brenner is not correct, because
    parameters in the vdW potential are fitted with other potential.
-   Check num. of divisions and cutoff length in \'in.pmd\' before the
    simulation run.
-   Use variable array TAG instead of IS, because of the efficient
    parallelization. TAG includes species, index of FMV, total id.

History
-------

-   2014.02.11 Input and output files are seperated to folders of name
    of output number.
-   2014.01.?? Units are changed from atomic unit to eV, Angstrom, and
    mass of 1/12.
-   2012.05.22 Change input format to user readable one, and input file
    name became \'in.pmd\'.
-   2009.05.12 Use TAG instead of IS. Reduce num of MPI message passing
    in BAMOVE and BACOPY.
-   2009.04.28 Add Brenner potential with van der Waals term. Also
    simple velocity scaling is added.
-   2009.03.24 Add smoothing for embedded term, too.
-   2009.03.20 Bugs fixed about copy of IFMV in subroutine BAMOVE.
-   2009.03.20 Add smoothing for 2-body terms in EAM potential.
# Atom-configuration file

Original file format described here is used in *pmd* .

## File format

    1 :!
    2 :!  specorder:  W  H
    3 :!    
    4 :  2.855300000E+000
    5 :  3.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
    6 :  0.00000000000000E+000  3.00000000000000E+000  0.00000000000000E+000
    7 :  0.00000000000000E+000  0.00000000000000E+000  3.00000000000000E+000
    8 :  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
    9 :  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
    10:  0.00000000000000E+000  0.00000000000000E+000  0.00000000000000E+000
    11:        55
    12:  1.10000000000001E+000  1.000E-007  1.000E-007  1.000E-007  3.62E-004  1.60E-004  9.60E-004
    13:  1.10000000000002E+000  1.666E-001  1.666E-001  1.666E-001  3.62E-004  1.60E-004  9.60E-004
    14:  ...
    15:  1.10000000000054E+000  9.00E-001  9.00E-007  9.00E-001  3.62E-004  1.60E-004  9.60E-004
    16:  2.10000000000055E+000  5.33E-001  5.33E-001  5.33E-001  3.62E-004  1.60E-004  9.60E-004

Here, line numbers are shown for the ease of explanation.

- **Line 1-3**: Lines begin with `!` are treated as comment lines. There are some
    keywords that are used to specify some additional data to *pmd* when
    they are at a comment line at the beginning.

    -   `specorder:` specifies the species order used in *pmd*.

    !!! Note
        **specorder** must be specified in the current *pmd* (since
        *rev190515*), as the masses and the interatomic potentials are
        determined using this information.
        
- **Line 4**: Superficial or apparent lattie constant. This value is to be
    multiplied to the cell vectors below to obtain absolute cell
    vectors.
- **Line 5 to 7**: Lattice vectors. The 2nd line is *a1* vector, 3rd line for *a2*, and
    4th line for *a3*. Columns 1, 2, and 3 are *x* , *y* , *z*
    components of each vector.
- **Line 8 to 10**: Velocities of lattice vectors that are used in *NpT* -ensemble
    simulation which involves lattice deformation.
- **Line 11**: Number of atoms in the system or decomposed region.
- **After line 11**: One atom information per one line.
- **1st column after line 11**: Tag of an atom. The digit in the one's place means species of the
    atom. The digit in the tenth's place is *ifmv* value which controls
    the direction of motion of the atom.
- **2-4th column after line 11**: *x* , *y* , and *z* coordinates of the atom normalized by the
    lattice vectors. Thus they should be in `(0:1]`
- **5-7th column after line 11**:  *x* , *y* , and *z* components of the atom velocities that are also normalized.

------------------------------------------------------------------------

## Sample Fortran code

```Fortran
open(ionum,file=cfname,status='replace')
write(ionum,'(es23.14e3)') hunit
write(ionum,'(3es23.14e3)') (((h(ia,ib,l)/hunit,ia=1,3) &
     ,ib=1,3),l=0,1)
write(ionum,'(i10)') natm
do i=1,ntot
  write(ionum,'(7es23.14e3)') tag(i), rtot(1:3,i), vtot(1:3,i)
  close(ionum)
enddo
```

Detail explanations of variables are omitted. Users can write their own
code by following this sample Fortran code.

------------------------------------------------------------------------

## Format conversion


There is a python utility, `napsys.py` that can convert files amoung the
following formats,

- `pmd`: input format for *pmd* program.
- `POSCAR`: input format of *VASP* program.
- `dump`: output format of *LAMMPS*-dump command.

You can use like following, :

    $ python /path/to/nap/nappy/napsys.py convert pmdini POSCAR

Here `pmdini` file will be converted to `POSCAR` file in POSCAR format,
where the file format is determined automatically from file names. Users
can also specify input/output file format by the options `--in-format`
and `--out-format`.

See the help message for more details as, :

    $ python /path/to/nap/nappy/napsys.py -h

------------------------------------------------------------------------

## Make crystalline structures

There is also a python utility, `cell_maker.py`, which makes typical
conventional crystalline structures. You can make a *pmd* format file of
diamond structured cubic system with 8 atoms as,

    $ python /path/to/nap/nappy/mkcell/cell_maker.py dia -l 5.473 -o pmdini

The option `-l` specifies the lattice constant of the lattice. Output
format is automatically detected from the file name. You can also make
*fcc*, *bcc*, *sc (simple cubic)*, and *hcp* structures as well.
# Introduction

*pmd* is an acronym of **parallel molecular dynamics** which means that
molecular dynamics (MD) using spatial decomposition technique on
parallel (ditributed-memory) computers. And *pmd* is a part of the
**Nagoya atomistic-simulation package (nap)**.

The main features of *pmd* are the following:

- several interatomic potentials for solid state systems are available;
- **Deep neural-network (DNN) interatomic potential**;
- **QEq** or **variable charge** Coulombic potential;
- parallel computation using spatial decomposition technique;
- efficient searching of neighbor atoms using linked-list cell
  method;
- structure relaxation using simple velocity damping or **FIRE**
  algorithm;
- thermostats: Berendsen and Langevin;
- barostat: Berendsen;
- **variable-timestep** MD for high-energy ion-bombardment
  simulation;
- **non-equilibrium MD (NEMD)** for heat flux simulation;
- **two-temperature model MD (TTM-MD)** for laser-ablation
  simulation.

Since this program has been developed for the purpose of personal
research tool, there are not so many functionalities. And there are some
(open source) MD programs that can do almost the same thing that *pmd*
can do. But there are some features only *pmd* or the parent packange
**nap** can do. Please feel free to contact me to ask anything if you
want to use *pmd* or **nap** for your specific purpose.

Bug reports and questions about *pmd* and **nap** are welcome, but I am
afraid that I might not be able to respond all the reports or questions.


## Requirements

*pmd* can be executed in Unix/Linux, macOS X, and Windows with using a
Fortran compiler and an MPI library.

For some analysis tools written in Python language, you may also need
*Python 3.* and some python-utilities such as *numpy*, *scipy*,
*pandas*, *docopt* and *ASE*.
# FAQs

------------------------------------------------------------------------

## There are not enough slots available...

If you get the following error message when you run `pmd` with some
number of parallel processes,

    --------------------------------------------------------------------------
    There are not enough slots available in the system to satisfy the 4 slots
    that were requested by the application:
      pmd

    Either request fewer slots for your application, or make more slots available
    for use.
    --------------------------------------------------------------------------

you may have to increase the number of slots available by the
`--oversubscribe` option as, :

    $ mpirun --oversubscribe -np 4 pmd

Reference:

-   <https://stackoverflow.com/questions/35704637/mpirun-not-enough-slots-available>

------------------------------------------------------------------------

## A system call failed during shared memory initialization...

The following error message could appear at the end of `pmd` output,
when the number of MPI processes and the number of automatically
determined spatial decomposition divisions are different.

    --------------------------------------------------------------------------
    A system call failed during shared memory initialization that should
    not have.  It is likely that your MPI job will now either abort or
    experience performance degradation.

      Local host:  mbp-rk-eth.ogt.nitech.ac.jp
      System call: unlink(2) /var/folders/fh/cd90qtsj3d1fr3wmjw954xt40000gn/T//ompi.mbp-rk-eth.501/pid.52045/1/vader_segment.mbp-rk-eth.95f30001.2
      Error:       No such file or directory (errno 2)
    --------------------------------------------------------------------------

But the `pmd` results are correct. You can remove it by setting an
environment variable as, :

    export OMPI_MCA_btl=self,tcp

Reference:

-   <https://github.com/open-mpi/ompi/issues/6518>

-------------
Neural-network (NN) potential {#neural_network}
=============================

NN potential requires following two input files in the working
directory.

-   `in.const.NN`: network structure, types of symmetry functions,
    parameters of the symmetry functions, and interacting pairs.
-   `in.params.NN`: number of weights, cutoff radii for two- and
    three-body interactions, and weight values of the network.

These two files must be consistent such that the number of weights must
correspond to the number of symmetry functions, number of layers, and
number of nodes in each layer. The examples of these files can be found
in `pmd/force_params/NN_??????` directories.

There are some rules in `in.const.NN` as follows:

-   Three-body terms must come after all the two-body terms.

-   Within two-body terms, terms for the specific pair must be written
    in consecutive lines such as :

        1   1   1   10.000     3.0000
        1   1   1   10.000     4.0000
        1   1   2   10.000     3.0000
        1   1   2   10.000     4.0000
        1   1   2   10.000     5.0000
        1   1   3   10.000     3.0000

    Not like :

        1   1   1   10.000     3.0000
        1   1   2   10.000     3.0000
        1   1   3   10.000     3.0000
        1   1   1   10.000     4.0000
        1   1   1   10.000     5.0000
        1   1   2   10.000     4.0000

    where the pairs 1-2 and 1-3 (see 2nd and 3rd columns) appear before
    finishing all the inputs of 1-1 pair.
# Stillinger-Weber (SW) potential


The SW potential requires a parameter file `in.params.SW` of the
following format.

    unit  2.1678  2.0951
    1  1  7.049556277  0.6022245584  4.0  0.0  1.0  1.8
    1  1  1  21.0  1.20

-   Line starting with `unit` defines the units of energy and length.
-   Line with 8 entries is for two-body parameters with the format,

        isp, jsp, Eij, Aij, Bij, pij, qij, cij, rcij

-   Line with 5 entries is for three-body parameters with the format,

        isp, jsp, ksp, sij, tij
Neural-network (NN2) potential {#NN2}
==============================

::: {.note}
::: {.title}
Note
:::

The NN potential was updated to ver 2.0 on 7 June 2018 by the
force-field name [NN2]{.title-ref}. The document here is for the NN2
potential. If one has to use ver 1.0, see
`neural_network`{.interpreted-text role="doc"}.
:::

NN2 potential requires following two input files in the working
directory.

-   `in.params.desc`: types of symmetry functions, parameters of the
    symmetry functions, their cutoff radii, and interacting pairs.
-   `in.params.NN`: NN structures, number of layers, number of nodes for
    each layer, and weight values of the network.

These two files must be consistent such that the number of weights must
correspond to the number of symmetry functions, number of layers, and
number of nodes in each layer. The examples of these files can be found
in `pmd/force_params/NN2_??????` directories.

in.params.desc
--------------

::: {.note}
::: {.title}
Note
:::

The format of `in.params.desc` file has changed from the previous one at
around May 2019, where the species are written directly by their acronym
name such as *Si* not by species-ID (digit 1-9). If you want to use the
previous-format `in.params.desc`, you should modify it by replacing
species-ID with species name.
:::

The format of `in.params.desc` is like the following, :

    2    20
    1  W   W   10.000    2.000
    1  W   W   10.000    3.000
    1  W   W   10.000    4.000
    1  W   H   10.000    2.000
    1  W   H   10.000    3.000
    1  W   H   10.000    4.000
    2  W   W   -0.900
    2  W   W   -0.800
    2  W   W   -0.700
    ...

-   1st line has two entries, *number of speceis* and *number of
    descriptors*.
-   Following lines have each descriptor information, the 1st entry is
    the type of descriptor, 2nd and 3rd are species of interaction pair,
    from the 4th to the end are parameters of the descriptor. The number
    of parameters depend on the type of descriptor.

in.params.NN2
-------------

`in.params.NN2` file should have the following format. :

    1   18   10
    -3.64106023330479E-001 -1.0000E-01  1.0000E-01
    -2.01340565152879E+000 -1.0000E-01  1.0000E-01

where three digits in the 1st line are *number of layers*, *number of
input nodes*, and *number of nodes in the 1st layer*. There should be
190 (= 18\*10 + 10) following lines (in this case), with *NN weight* and
two dummy values.
# Deep neural-network (DNN) potential

!!! Note
    The NN potentials (NN and NN2) were replaced with this `DNN` potential
    in January 2020. It is strongly recommended to use `DNN` potential
    instead of `NN2`, because `DNN` is a super-set of `NN2` and `NN2` is no
    longer maintenanced in the future.


DNN potential requires the following two input files in the working
directory.

-   `in.params.desc`: types of symmetry functions, parameters of the
    symmetry functions, their cutoff radii, and interacting pairs.
-   `in.params.DNN`: NN structures, number of layers, number of nodes
    for each layer, and weight values of the network.

These two files must be consistent such that the number of descriptors
must correspond to the number of inputs in the NN. The examples of these
files can be found in `pmd/force_params/DNN_??????` directories.

## in.params.desc

!!! Note
    The format of `in.params.desc` file has changed from the previous one at
    around May 2019, where the species are written directly by their acronym
    name such as *Si* not by species-ID (digit 1-9). If you want to use the
    previous-format `in.params.desc`, you should modify it by replacing
    species-ID with species name.

The format of `in.params.desc` is like the following,

    2    20
    1    W   W   10.000    2.000
    1    W   W   10.000    3.000
    1    W   W   10.000    4.000
    1    W   H   10.000    2.000
    1    W   H   10.000    3.000
    1    W   H   10.000    4.000
    2    W   W   -0.900
    2    W   W   -0.800
    2    W   W   -0.700
    ...

-   1st line has two entries, *number of speceis* and *number of
    descriptors*.
-   Following lines have each descriptor information, the 1st entry is
    the type of descriptor, 2nd and 3rd are species of interaction pair,
    from the 4th to the end are parameters of the descriptor. The number
    of parameters depend on the type of descriptor.

## in.params.DNN

`in.params.DNN` file should have the following format.

    !  sigtype: 2
    ! 
       3   20   10  10  5
     -3.64106023330479E-001 -1.0000E-01  1.0000E-01
     -2.01340565152879E+000 -1.0000E-01  1.0000E-01

-   Lines starting with `!` are comment lines. If any of special keyword
    is found just after `!`, an option will be passed to the program.
-   The 1st entry of the 1st line following comments is the number of
    hidden layers in the NN.
-   The 2nd entry is the number of nodes in 0-th layer, which is called
    input layer.
-   The following digits are the number of nodes in hidden layers. The
    number of these digits must be the same as the number of hidden
    layers given by the 1st entry.
-   Following lines include weight value, lower and upper bounds, which
    are only used in `fitpot`.
# Morse potential


Potential form of the Morse potential is,

$$\begin{equation*}
E_{ij}(R_{ij}) = D_{ij} \left\{ [\exp (\alpha_{ij} (R_{0,ij} -R_{ij}))-1]^2 -1\right\}
\end{equation*}$$

The potential requires a parameter file `in.params.Morse` that has the
following format.

    #  isp, jsp, D_ij, alpha_ij, R0_ij
    1    1    3.7    2.0    1.68
    1    2    3.2    1.9    2.5
    2    2    2.5    2.1    2.3

In the above case, there are three interactions between 1-1, 1-2 and
2-2. Each interaction has three parameters:

-   *D_ij*: depth of the potential curve (eV)
-   *alpha_ij*: related to the width of the potential well
    (Ang.^{-1}).
-   *R0_ij*: position of the potential minimum (Ang.)
# Coulomb potential (Ewald or short/screend)

Coulomb potential requires a parameter file `in.params.Coulomb` that has
the following format.

    charges  fixed
      Si   1.0
      O   -2.0
    interactions
      Si  Si
      Si  O
      O   O
    terms  short

    sigma  2.5

In this format, black lines are neglected. There are some keywords:

- `charges` : `fixed` or `variable` (or `qeq`) -- Followed by the lines of species charges, e.g., $q_1 = 1.0$ and $q_2 = -2.0$.
- `interactions` *(optional)* -- Followed by the pairs of species. If not specified, all the  interactions are taken into account.
- `terms` : `full`, `long` or `short`/`screened_cut` -- Either full Ewald terms, long-range term only, short-range term only, or short-term with smooth cutoff.
- `sigma` -- Width of the Gaussian charge distribution, which is related to the
    accuracy in case of Ewald method.

## Variable charge or QEq

Coulomb potential can treat **variable charge** or **QEq** by specifying
`variable` or `qeq` to the `charges` keyword as shown below.

    charges  variable
      Si  4.7695    8.7893    0.0   0.0  2.4
      O   7.5405    15.8067   0.0  -1.2  0.0
    interactions
      Si  Si
      Si  O
      O   O
    terms  short
    sigma  2.5
    conv_eps  1.0d-6

Here, `charges variable` requires some following lines that have

    name,  chi,  Jii,  E0,  qlow,  qup

-   `name`: name of the chemical species
-   `chi`: electronegativity of the species (eV)
-   `Jii`: hardness of the species (eV)
-   `E0`: atomic energy (eV)
-   `qlow`: lower limit of the charge of the species
-   `qup`: upper limit of the charge of the species
-   `conv_eps`: Convergence criterion of the charge optimization. Usually it should
    be very small to achieve good energy conservation.
# What's qmcl
qmcl is a program which performs QM/MM simulation using VASP and pmd for QM and CL calculation, respectively.
The package includes various force interfaces such as QM only, CL only, and QM/MM.

# Who made this?
* Ryo KOBAYASHI
* Assistant Professor in the department of mechanical engineering, Nagoya Institute of Technology. (2014-03-25)

# Compilation and usage
See the manual web site below (but Japanese only),
http://locs.bw.nitech.ac.jp/~kobayashi/qmcl_manual

