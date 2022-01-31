# Citing Castro

There are a number of Castro papers that describe parts of the
algorithm.  We ask that you cite all of the appropriate papers
describing the capabilities you used.

## General code use

If you use Castro, we appreciate you citing the most recent code paper:

```
@article{Almgren2020,
  doi = {10.21105/joss.02513},
  url = {https://doi.org/10.21105/joss.02513},
  year = {2020},
  publisher = {The Open Journal},
  volume = {5},
  number = {54},
  pages = {2513},
  author = {Ann Almgren and Maria Barrios Sazo and John Bell and Alice Harpole and Max Katz and Jean Sexton and Donald Willcox and Weiqun Zhang and Michael Zingale},
  title = {CASTRO: A Massively Parallel Compressible Astrophysics Simulation Code},
  journal = {Journal of Open Source Software}
}
```

You are welcome to cite the original code paper as well (which
provides some details on the algorithmic implementations):

```
    @ARTICLE{2010ApJ...715.1221A,
       author = {{Almgren}, A.~S. and {Beckner}, V.~E. and {Bell},
		      J.~B. and {Day}, M.~S. and {Howell}, L.~H. and
		      {Joggerst}, C.~C. and {Lijewski}, M.~J. and
		      {Nonaka}, A. and {Singer}, M. and {Zingale}, M.},
	title = "{CASTRO: A New Compressible Astrophysical
		      Solver. I. Hydrodynamics and Self-gravity}",
      journal = {\apj},
    archivePrefix = "arXiv",
       eprint = {1005.0114},
     primaryClass = "astro-ph.IM",
     keywords = {equation of state, gravitation, hydrodynamics, methods:
		      numerical, nuclear reactions, nucleosynthesis,
		      abundances},
	 year = 2010,
	month = jun,
       volume = 715,
	pages = {1221-1238},
	  doi = {10.1088/0004-637X/715/2/1221},
       adsurl = {http://adsabs.harvard.edu/abs/2010ApJ...715.1221A},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

You should also cite the zenodo DOI for the code release.  A bibtex
entry for the latest release can be found here:

https://doi.org/10.5281/zenodo.2301848

## Radiation hydrodynamics

If you use the radiation hydrodynamics capabilities, please additionally
cite the following:

```
    @ARTICLE{2011ApJS..196...20Z,
       author = {{Zhang}, W. and {Howell}, L. and {Almgren}, A. and
		      {Burrows}, A. and {Bell}, J.},
	title = "{CASTRO: A New Compressible Astrophysical
		      Solver. II. Gray Radiation Hydrodynamics}",
      journal = {\apjs},
    archivePrefix = "arXiv",
       eprint = {1105.2466},
     primaryClass = "astro-ph.IM",
     keywords = {diffusion, hydrodynamics, methods: numerical, radiative
		      transfer},
	 year = 2011,
	month = oct,
       volume = 196,
	  eid = {20},
	pages = {20},
	  doi = {10.1088/0067-0049/196/2/20},
       adsurl = {http://adsabs.harvard.edu/abs/2011ApJS..196...20Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

```
    @ARTICLE{2013ApJS..204....7Z,
       author = {{Zhang}, W. and {Howell}, L. and {Almgren}, A. and
		      {Burrows}, A. and {Dolence}, J. and {Bell}, J.},
	title = "{CASTRO: A New Compressible Astrophysical
		      Solver. III. Multigroup Radiation Hydrodynamics}",
      journal = {\apjs},
    archivePrefix = "arXiv",
       eprint = {1207.3845},
     primaryClass = "astro-ph.IM",
     keywords = {diffusion, hydrodynamics, methods: numerical, radiative
		      transfer },
	 year = 2013,
	month = jan,
       volume = 204,
	  eid = {7},
	pages = {7},
	  doi = {10.1088/0067-0049/204/1/7},
       adsurl = {http://adsabs.harvard.edu/abs/2013ApJS..204....7Z},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## Gravity and rotation

If you use the Poisson gravity solver or rotation, please cite the
following, which included a lot of improvements to the coupling of
hydro and gravity:

```
    @ARTICLE{2016ApJ...819...94K,
       author = {{Katz}, M.~P. and {Zingale}, M. and {Calder}, A.~C. and
		      {Swesty}, F.~D. and {Almgren}, A.~S. and {Zhang},
		      W.},
	title = "{White Dwarf Mergers on Adaptive Meshes. I. Methodology
		      and Code Verification}",
      journal = {\apj},
    archivePrefix = "arXiv",
       eprint = {1512.06099},
     primaryClass = "astro-ph.HE",
     keywords = {hydrodynamics, methods: numerical, supernovae: general,
		      white dwarfs},
	 year = 2016,
	month = mar,
       volume = 819,
	  eid = {94},
	pages = {94},
	  doi = {10.3847/0004-637X/819/2/94},
       adsurl = {http://adsabs.harvard.edu/abs/2016ApJ...819...94K},
      adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## GPUs and scaling

For CPU performance numbers, please cite:

```
    @INPROCEEDINGS{2018JPhCS1031a2024Z,
           author = {{Zingale}, M. and {Almgren}, A.~S. and {Barrios Sazo}, M.~G. and
             {Beckner}, V.~E. and {Bell}, J.~B. and {Friesen}, B. and
             {Jacobs}, A.~M. and {Katz}, M.~P. and {Malone}, C.~M. and
             {Nonaka}, A.~J. and {Willcox}, D.~E. and {Zhang}, W.},
            title = "{Meeting the Challenges of Modeling Astrophysical Thermonuclear Explosions: Castro, Maestro, and the AMReX Astrophysics Suite}",
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics},
        booktitle = {Journal of Physics Conference Series},
             year = 2018,
           series = {Journal of Physics Conference Series},
           volume = {1031},
            month = may,
              eid = {012024},
            pages = {012024},
              doi = {10.1088/1742-6596/1031/1/012024},
    archivePrefix = {arXiv},
           eprint = {1711.06203},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2018JPhCS1031a2024Z},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

For GPU performance, please cite:

```
    @ARTICLE{2020arXiv200705218K,
           author = {{Katz}, Max P. and {Almgren}, Ann and {Barrios Sazo}, Maria and
             {Eiden}, Kiran and {Gott}, Kevin and {Harpole}, Alice and
             {Sexton}, Jean M. and {Willcox}, Don E. and {Zhang}, Weiqun and
             {Zingale}, Michael},
            title = "{Preparing Nuclear Astrophysics for Exascale}",
          journal = {arXiv e-prints},
         keywords = {Astrophysics - Instrumentation and Methods for Astrophysics, Astrophysics - High Energy Astrophysical Phenomena},
             year = 2020,
            month = jul,
              eid = {arXiv:2007.05218},
            pages = {arXiv:2007.05218},
    archivePrefix = {arXiv},
           eprint = {2007.05218},
     primaryClass = {astro-ph.IM},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2020arXiv200705218K},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

## Spectral deferred corrections

For the 2nd and 4th order SDC coupling of hydro and reactions, please cite:

```
    @ARTICLE{2019ApJ...886..105Z,
           author = {{Zingale}, M. and {Katz}, M.~P. and {Bell}, J.~B. and {Minion}, M.~L. and
             {Nonaka}, A.~J. and {Zhang}, W.},
            title = "{Improved Coupling of Hydrodynamics and Nuclear Reactions via Spectral Deferred Corrections}",
          journal = {\apj},
         keywords = {Hydrodynamics, Astrophysical fluid dynamics, Computational methods, Computational astronomy, Astronomy software, Nuclear astrophysics, Nucleosynthesis, Stellar nucleosynthesis, Physics - Computational Physics, Astrophysics - Instrumentation and Methods for Astrophysics},
             year = 2019,
            month = dec,
           volume = {886},
           number = {2},
              eid = {105},
            pages = {105},
              doi = {10.3847/1538-4357/ab4e1d},
    archivePrefix = {arXiv},
           eprint = {1908.03661},
     primaryClass = {physics.comp-ph},
           adsurl = {https://ui.adsabs.harvard.edu/abs/2019ApJ...886..105Z},
          adsnote = {Provided by the SAO/NASA Astrophysics Data System}
    }
```

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2301848.svg)](https://doi.org/10.5281/zenodo.2301848)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02513/status.svg)](https://doi.org/10.21105/joss.02513)
[![AMReX](https://amrex-codes.github.io/badges/powered%20by-AMReX-red.svg)](https://amrex-codes.github.io)
[![yt-project](https://img.shields.io/static/v1?label="works%20with"&message="yt"&color="blueviolet")](https://yt-project.org)
[![github pages](https://github.com/AMReX-Astro/Castro/workflows/github%20pages/badge.svg)](https://github.com/AMReX-Astro/Castro/actions?query=workflow%3A%22github+pages%22)

![Castro](https://github.com/AMReX-Astro/Castro/blob/development/Util/logo/castro_logo_hot_200.png)

*an adaptive mesh, astrophysical radiation hydrodynamics simulation code*

`Castro` is an adaptive-mesh compressible radiation / MHD / hydrodynamics
code for astrophysical flows.  `Castro` supports a general equation of
state, full Poisson gravity, and reactive flows, and is parallelized
with MPI + OpenMP for CPUs and MPI + CUDA for GPUs.

More information on Castro can be found here:

http://amrex-astro.github.io/Castro/


## Getting Started

The "Getting Started" section of the User's Guide walks you
through running your first problem:

https://amrex-astro.github.io/Castro/docs/getting_started.html

This will have you clone Castro and its dependencies (AMReX and
StarKiller Microphysics),

The User's Guide in written in re-structured text using Sphinx, with
the source in `Castro/Docs/`, and is built automatically
from the `development` branch.

## Running at Supercomputer Centers

Documentation for running the AMReX Astrophysics codes at popular
supercomputing centers can be found at:
https://amrex-astro.github.io/workflow/

## Call Graph

A doxygen-generated call graph for `Castro` is available here:

http://bender.astro.sunysb.edu/Castro/staging/Castro/html/


## Development Model:

Development generally follows the following ideas:

  * New features are committed to the `development` branch.

    Nightly regression testing is used to ensure that no answers
    change (or if they do, that the changes were expected).

    If a change is critical, we can cherry-pick the commit from
    `development` to `main`.

  * Contributions are welcomed from anyone in the form of a pull
    request from your fork of Castro, targeting the `development`
    branch. (If you mistakenly target `main`, we can change it
    for you.)

    Please add a line to `CHANGES.md` summarizing your change if it
    is a bug fix or new feature.  Reference the PR or issue as
    appropriate. Additionally, if your change fixes a bug (or if
    you find a bug but do not fix it), and there is no current
    issue describing the bug, please file a separate issue describing
    the bug, regardless of how significant the bug is. If possible,
    in both the `CHANGES.md` file and the issue, please cite the pull
    request numbers or git commit hashes where the problem was
    introduced and fixed, respectively.

    We will squash commits upon merge to have a clean history.
    *Please ensure that your PR title and and the PR summary field are
    descriptive, since these will be used for a squashed commit message.*

  * On the first workday of each month, we perform a merge of
    `development` into `main`, in coordination with `AMReX`,
    `Maestro`, and `Microphysics`.  For this merge to take place, we
    need to be passing the regression tests.

    To accommodate this need, we close the merge window into
    `development` a few days before the merge day.  While the merge
    window is closed, only bug fixes should be pushed into
    `development`.  Once the merge from `development` -> `main` is
    done, the merge window reopens.


## Core Developers

People who make a number of substantive contributions will be named
"core developers" of Castro.  The criteria for becoming a core
developer are flexible, but generally involve one of the following:

  * 10 non-merge commits to `Castro/Source/` or `Castro/Docs/`
    or one of the problems that is not your own science problem *or*

  * addition of a new algorithm / module  *or*

  * substantial input into the code design process or testing

Core developers will be recognized in the following ways:

  * invited to the group's slack team

  * listed in the User's Guide and website as a core developer

  * listed in the author list on the Zenodo DOI for the project
    (as given in the .zenodo.json file)

  * invited to co-author general code papers / proceedings describing
    Castro, its performance, etc.  (Note: science papers will always
    be left to the science leads to determine authorship).

If a core developer is inactive for 3 years, we may reassess their
status as a core developer.



## Getting help

We use Github discussions for asking general questions about the code:

https://github.com/AMReX-Astro/Castro/discussions
# 21.12

   * Tiling was added to main loop in MHD algorithm to enable 
     scaling performance increase when using multiple threads
     in with OpenMP. See issue #2038.

   * `castro.hse_fixed_temp` was added to allow for a fixed temperature
     at an HSE boundary. It can be enabled by setting it to a positive
     value and setting castro.hse_interp_temp=0. (#2042)

# 21.10

   * A new option, `castro.drive_initial_convection` was added that
     uses the temperature interpolated from the initial model instead
     of the value on the grid to call the reactions.  This helps
     prevent a reactive zone from burning in place before a convective
     velocity field is established to carry off the energy.

   * The `burn_weights` are no longer stored by default in the plotfile.
     Instead, they are now enabled by setting
     castro.store_burn_weights=1.  Additionally, they now give a better
     estimate of the cost for the numerical Jacobian (#1946, #1949)

   * `castro.change_max` is now required to be greater than 1.0. To enforce
     a timestep cap but still allow the timestep to decrease, use
     castro.max_dt. (#1976)

   * Gravity was modified to introduce parallel plane gravity with a
     point mass by setting the radius of the star by
     `castro.point_mass_location_offset` and the integer
     `castro.point_mass_offset_is_true` == 1. By default, both
     parameters are 0.0 and 0, respectively.


# 21.09

   * `castro.source_term_predictor` now works for simplified-SDC to
     disable the source predictor for the hydrodynamics states to the
     interface. (#1968)

   * `castro.add_sdc_react_source_to_advection` was added to disable
     react source to advection in simplified-SDC (#1969)

# 21.07

   * The sponge is now applied in a fully implicit manner at the end of
     the CTU advance, rather than using a predictor-corrector approach
     with time centering. This is more consistent with the original form
     of the sponge in Castro. (#1876)

   * Castro can now validate the runtime parameters set in the inputs
     file or on the commandline by setting
     castro.abort_on_invalid_params=1 (#1882)

# 21.06

   * Starting with this release, problem setups written in Fortran are
     no longer supported and will no longer work. Please consult the
     code documentation and example problem setups in Exec/ to understand
     the new problem setup format. If you need help converting a Fortran
     setup to C++, please file an issue. (#1728, #1732)

   * Sponge parameters are now only accepted through the inputs file; the
     &sponge namelist in the probin file is no longer read. (#1731)

   * Ambient parameters are now only accepted through the inputs file; the
     &ambient namelist in the probin file is no longer read. (#1742)

   * The update_sponge_params hook has been removed. (#1716)

   * The Fortran problem-specific source file, ext_src_nd.F90, has been
     removed. Problem-specific sources should be implemented in C++ in
     problem_source.H. (#1856)

   * Support for the legacy tagging scheme based on probin parameters (denerr,
     tempgrad, etc.) has been removed. These can be replaced with equivalent
     tagging criteria constructed in the inputs file; see the docs or examples
     in Exec/ to see how to use `amr.refinement_indicators`. (#1834)

   * The Fortran set_problem_tags hook has been removed. The C++ replacement
     is `problem_tagging()` in `problem_tagging.H`. (#1828)

   * The PrescribedGrav functionality has been removed (not replaced with a C++
     implementation). If you want to obtain the same functionality, you can use
     a problem-defined source term (look for problem_source in the documentation)
     and make the appropriate modification for applying it directly to the state
     (e.g. the momentum source term is rho * g). (#1854)

   * The custom radiation boundary using lo_bcflag and hi_bcflag coupled with
     an implementation of rbndry has been removed. (#1743)

   * We no longer store Reactions_Type in checkpoint files.  This means
     that newer versions of Castro will not restart from old version.
     
# 21.05

   * The parameter use_eos_in_riemann was removed -- we found no
     instances of it being used (#1623)

   * The option castro.apply_sources_consecutively was removed (#1636)

# 21.04

   * For simplified-SDC, we now correctly store only the reactive
     part of the update for rho_enuc, enuc, and rho_omegadot
     plotfile variables (#1602)

# 21.02

   * In axisymmetric geometry, there are additional forces that arise
     due to the changing direction of the unit vectors in the div{rho
     U U} term. The paper by Bernand-Champmartin discusses this. See
     issue #913. This adds those forces.  Note that the Coriolis force
     in 2-d axisymmetry is already consistent with a right-handed
     system despite our internal ordering of the state was r, z,
     theta.  (#923)

   * We can now set any of the Microphysics runtime parameters in the
     inputs file instead of probin.  Each group of parameters has a
     namesapce for the inputs file when set this way
     (e.g. eos.use_coulomb = 1), and the C++ inputs value will take
     precedence over the value set in probin if it is set in both
     places. (#1527)

# 21.01

   * The minimum C++ standard supported by Castro is now C++17. Most modern compilers
     support C++17; the notable exception is RHEL 7 and its derivatives like CentOS 7,
     where the default compiler is gcc 4.8. In that case a newer compiler must be loaded,
     particularly a version of gcc >= 7.0, for example by installing devtoolset-7 or (if
     running on an HPC cluster that provides modules) using a more recent gcc module. (#1506)

   * There can now be multiple _prob_params files throughout the source
     tree.  We read the problem's file last and that takes precedence over
     any other _prob_params files found. (#1500)

   * The timestep limiter dtnuc_T has been removed. dtnuc_e and dtnuc_X
     are still available for controlling the burning timestep. (#1501)

   * A bug was fixed in the 2nd order true SDC (with reactions) that
     was giving the wrong solution and convergence (#1494).  A second
     bug was fixed in defining the weights for the Radau quadrature
     when using true SDC (#1493)

   * Compiling with the PGI compiler is no longer a requirement for the CUDA build of Castro.
     We recommend using COMP=gnu with a version of gcc that is C++17 compliant (gcc >= 7).

# 20.12

   * An issue with incorrect application of HSE boundary conditions on derived quantities
     is now resolved (#1356). Also, at this point the old Fortran implementations hypfill,
     denfill, ext_hypfill, and ext_denfill have been removed; problem-specific boundary
     conditions should be implemented using the new C++ interface in this release from #1289.

   * The minimum supported Hypre version is now 2.19.0. (#1333)

   * We have switched from a Fortran to a C++ implementation of VODE in Microphysics.
     As a result we have also switched the Strang and simplified SDC burners in Castro
     to use this C++ implementation. Most networks used in Castro have already been
     ported to C++. While networks are not required to have a C++ implementation,
     networks implemented only in Fortran  will not be useable on GPUs, and eventually
     we will use C++ only. (#1313)

   * `problem_checkpoint` and `problem_restart` are moved to C++ from Fortran. See
     Exec/science/wdmerger for an example of the new scheme. `Problem.f90` and `Problem_F.H`
     are now deleted from the code; if you were using these to implement problem-specific
     functionality, you can still manually add these files to the `Make.package` for your
     problem setup. (#1311)

   * For setups using Poisson gravity, tagging is now turned off in locations where
     the fine levels would have been adjacent to a physical boundary. (This previously
     led to an abort.) (#1302)

   * An interface for doing problem tagging in C++ has been added. (#1289)

   * Simplified SDC now only supports the C++ integrators (#1294)

   * MHD problems can now do the magnetic field initialization in C++
     (#1298)

# 20.11

   * The minimum C++ standard supported by Castro is now C++14. Most modern compilers
     support C++14; the notable exception is RHEL 7 and its derivatives like CentOS 7,
     where the default compiler is gcc 4.8. In that case a newer compiler must be loaded,
     particularly a version of gcc >= 5.0, for example by installing devtoolset-7 or (if
     running on an HPC cluster that provides modules) using a more recent gcc module. (#1284)

   * A new option, `castro.retry_small_density_cutoff`, has been added. In some
     cases a small or negative density retry may be triggered on an update that
     moves a zone already close to small_dens just below it. This is not uncommon
     for "ambient"/"fluff" material outside a star. Since these zones are not
     dynamically important anyway, triggering a retry is unnecessary (and possibly
     counterproductive, since it may require a very small timestep to avoid). By
     setting this cutoff value appropriately, the retry will be skipped if the
     density of the zone prior to the update was below the cutoff. (#1273)

# 20.10

   * A new refinement scheme using the inputs file rather than the Fortran
     tagging namelist has been added. (#1243, #1246) As an example, consider:

     ```
     amr.refinement_indicators = dens temp

     amr.refine.dens.max_level = 1
     amr.refine.dens.value_greater = 2.0
     amr.refine.dens.field_name = density

     amr.refine.temp.max_level = 2
     amr.refine.temp.value_less = 1.0
     amr.refine.temp.field_name = Temp
     ```

     `amr.refinement_indicators` is a list of user-defined names for refinement
     schemes. For each defined name, amr.refine.<name> accepts predefined fields
     describing when to tag. In the current implementation, these are `max_level`
     (maximum level to refine to), `start_time` (when to start tagging), `end_time`
     (when to stop tagging), `value_greater` (value above which we refine),
     `value_less` (value below which to refine), `gradient` (absolute value of the
     difference between adjacent cells above which we refine), and `field_name`
     (name of the string defining the field in the code). If a refinement indicator
     is added, either `value_greater`, `value_less`, or `gradient` must be provided.

   * Automatic problem parameter configuration is now available to every
     problem by placing a _prob_params file in your problem directory.
     Examples can be found in most of the problems in Castro/Exec, and you
     can look at the "Setting Up Your Own Problem" section of the documentation
     for more information. This functionality is optional, however note that
     a file containing a Fortran module named "probdata_module" is now
     automatically generated, so if you have a legacy probdata.F90 file
     containing a module with that name it should be renamed. (#1210)

   * The variable "center" (in the `problem` namespace) is now part of this
     automatically generated probdata module; at the present time, the only
     valid way to change the problem center to a value other than zero is in
     amrex_probinit(). (#1222)

   * Initialization of these problem parameters is now done automatically for
     you, so a call to probdata_init() is no longer required in amrex_probinit(). (#1226)

   * Problems may now be initialized in C++ instead of Fortran. Instead of implementing
     amrex_probinit() and ca_initdata(), the problem should implement the analogous
     functions initialize_problem() and initialize_problem_state_data(). If you switch to
     the new C++ initialization, be sure to delete your Prob_nd.F90 file. By default both
     implementations do nothing, so you can pick either one but do not pick both. (#1227)

   * The external heat source term routines have been ported to C++
     (#1191).  Any problem using an external heat source should look
     convert the code over to C++.

   * The interpolate_nd.F90 file has been moved to Util/interpolate and
     is only compiled into Castro if you set USE_INTERPOLATE=TRUE

# 20.09

   * Reactions now work with MHD (#1179)

   * MHD now uses the main slope routine (#1058) The order of the
     slope is now controlled by plm_iorder, just as with hydro.  There
     is an additional option, mhd_limit_characteristic, that
     determines if the limiting is done on the primitive or
     characteristic variables (the default).

# 20.08

   * Rotation_Type has been removed from StateData. (#1128)

   * castro.use_post_step_regrid now unconditionally regrids after
     every timestep on every level. (#898)

   * An issue with gravity.max_solve_level resulting in accesses to invalid data
     (#469, #1118) has been resolved. (#1123)

   * If castro.speed_limit is set to a number greater than zero, this
     will now be strictly enforced on the magnitude of the velocity. (#1115)

   * When using AMR and gravity or rotation, the source terms applied after
     a reflux would have been incorrect if the previous timestep had a retry
     (#1020). This has now been fixed. (#1112)

   * We now have the ability to access the problem-specific runtime
     parameters in C++ (#1093)

# 20.07

   * The master branch has been renamed the main branch. If you have an
     existing clone of Castro, then do the following to update for this
     change. First, do `git checkout master` if you're not already on the
     old master branch. Then do `git pull`. This will gather the updates
     to the repo, but will fail with the message `Your configuration specifies
     to merge with the ref 'refs/heads/master' from the remote, but no such ref
     was fetched.` Then you can simply do `git checkout main` and your local
     repo should automatically switch to that branch and track updates from
     the upstream repo on GitHub. If you like, you can then delete the old
     master branch with `git branch -D master`.

   * The CUDA build no longer has a requirement that amr.blocking_factor
     be a multiple of 8. Though this is recommended for performance reasons,
     it was previously required due to correctness reasons because of the
     use of an AMReX Fortran function, amrex_filccn. As noted in #1048, this
     function is no longer required due to recent changes in Castro (problems
     overriding bc_fill_nd.F90 or bc_ext_fill_nd.F90 do not need to provide an
     initial fill of the ghost zone data before implementing their specific
     boundary conditions; this is now done for you). Calling this function
     may now result in race conditions and correctness issues in the CUDA
     build, so it should be removed from any problem setups. (#1049)

   * The functionality that permitted the rotation rate to change as a
     function of time, castro.rotation_include_domegadt and
     castro.rotational_dPdt, has been removed. (#1045)

   * A CUDA illegal memory access error in Poisson gravity and diffusion
     has been fixed (#1039).

   * The parameter castro.track_grid_losses has been removed. (#1035)

   * The parameter castro.print_fortran_warnings, which no longer had any
     effect, has been removed. (#1036)

   * PPM reconstruction has been added to the MHD solver (#1002)

   * The Reactions_Type StateData has been reworked so that its first
     NumSpec components are rho * omegadot rather than omegadot; then,
     the NumAux auxiliary components are stored, if the network has any
     auxiliary variables; then, rho * enuc is stored (enuc itself is
     removed), and finally the burn weights are stored. The checkpoint
     version has been incremented, so this version of the code cannot
     restart from checkpoints generated with earlier versions of the
     code. (#927)

   * A bug where refluxing between AMR levels resulted in incorrect results
     when a retry occurred in the previous timestep has been fixed. (#1018)

# 20.06

   * The parameter castro.density_reset_method has been removed. A density
     reset now unconditionally sets the density to small_dens, the temperature
     to small_temp, and zeros out the velocities. (#989)

   * A constrained-transport corner transport upwind MHD solver has been
     added.  This can be used by compiling with USE_MPI = TRUE.  Presently
     it only works for a single level (no AMR).  (#307)

   * A burning timestep limiter dtnuc_T has been added which restricts the
     burning from updating the temperature by more than the factor
     dtnuc_T * T / dT/dt. (#972)

   * The reaction weights metric implemented in version 20.05 (#863) has been
     added to the simplified SDC reactions driver. (#930)

   * When using the simplified SDC integration scheme, we now save new-time
     Reactions_Type data to plotfiles. (#929)

# 20.05

   * The parameter use_custom_knapsack_weights and its associated
     functionality have been removed. (#877)

   * We've changed how the runtime parameters are stored.  Previously
     they were static members of their respective class, but this
     prevented their use in lambda-capture functions on GPUs.  Now the
     runtime parameters are grouped into namespaces as extern managed
     data. (#873)

   * We currently have a scheme for storing reactions weightings, which
     are a measure of the number of RHS evaluations during the burn and
     therefore a proxy for the difficulty of the burn. These weights were
     added as separate StateData depending on the runtime option
     use_custom_knapsack_weights. Now, instead we place the weights
     directly in the Reactions_Type StateData as a new component.

     The number of ghost zones in Reactions_Type is increased to 4.

     The checkpoint version has now been incremented; this version of the
     code will not be able to restart from a checkpoint generated by earlier
     versions of the code. (#863)

   * The meaning of dt_cutoff has changed: it is now the fraction of the
     current simulation time which dt may be no smaller than, instead of
     being an absolute measure. We now have set a non-zero default
     (1.e-12) as well. (#865)

   * Backwards compatibility in restarting from a checkpoint is no longer
     supported. Checkpoints from older versions of the code (as determined
     by the checkpoint version in the CastroHeader file in the checkpoint
     directory) cannot be restarted from. (#860)

   * Added an option to do CTU reactions in C++.  A compile flag
     USE_CXX_REACTIONS is added which switches to the C++ integrator
     in Microphysics. Since we will be doing a phased implementation
     of the networks in Microphysics, this is opt-in for now.  (#836)

   * More of the core routines have been ported to C++, including the
     hydro and diffusion timestep estimators (#853) and the sponge
     (#857)

   * AMReX provides CpuBndryFuncFab and GpuBndryFuncFab which are very
     similar to what generic_fill and hypfill did. The AMReX
     implementations are now used. We still have a hypfill and denfill
     function, so that existing problems are not broken, but the main
     one in Source/ no longer calls amrex_filcc (it only has the
     ambient code now). The problems that do override bc_fill_nd.F90
     are thus no longer required to call amrex_filcc. (#837)

   * We now always issue a timestep retry if the density after an
     advance is negative (or less than small_dens). The parameter
     castro.retry_neg_dens_factor is removed. The parameter
     castro.retry_tolerance is also removed as it no longer has
     any effect. (#796)

   * The timestep control parameter castro.change_max now also will
     prevent the timestep by shrinking too much in one timestep
     (previously it would only prevent it from growing too much).
     If change_max is violated in a timestep we will do a retry
     to take more graceful steps. (#844)

   * We now check if the problem setup initialized the density or
     temperature to a value near small_dens or small_temp and abort.
     If this happens, the recourse is to adjust small_dens and
     small_temp to a meaningful value for your problem.  (#822)

   * The src_q multifab was removed and instead we convert the
     conserved state sources to primitive state sources FAB by FAB.
     This saves a lot of memory at the expense of an EOS call. (#829)

   * The plm_well_balanced option was removed.  It was essentially the
     same as use_pslope except it was lower order and only worked with
     constant gravity.  use_pslope now works with both CTU+PLM and
     SDC2+PLM.  A new test problem, hse_convergence, was added to look
     at the behavior of the different reconstruction methods with HSE.

# 20.04

   * A potential undefined flux from the HLL solver when using
     hybrid_riemann has been fixed (#823)

   * The parameter castro.allow_small_energy has been removed. The
     code behavior is now similar to what it would have been with
     allow_small_energy == 0 (the internal energy can never be
     smaller than that allowed by small_temp). (#817)

   * The BC interfaces have been merged and converted to a new FAB
     interface as part of the port to C++. (#819)

   * All boundary fill interfaces other than hypfill and denfill have
     been removed. So, we no longer support overriding the boundary
     conditions for data other than State_Type. Radiation still has
     its own set of custom boundary conditions that can be accessed
     through the inputs file, as described in the docs. (#815)

   * The conversion of the CTU hydrodynamics code to C++ continues.
     The Riemann solvers were converted to C++ (#801) and the
     hybrid momentum routines (#805), the PLM reconstruction (#814),
     the conversion of primitive to conserved variables (#804)

   * We've changed how the backup for retries is done.  Presently if
     use_retry is enabled we make a pre-emptive copy of the StateData
     right at the beginning of the timestep.  Now we only backup when
     we detect that a retry is needed (#812)

# 20.03

   * We now depend on the fundamental constants from Microphysics
     instead of keep our own copy in Castro (#787)

   * We removed the ppm_predict_gammae option for the CTU hydro solver.
     This was not used frequently and did not show much difference with
     the default (rho e) reconstruction. (#780)

   * The Microphysics "extern" parameters are now available in C++

   * We've started converting the CTU hydro solver from Fortran to C++
     (#731).  The PPM reconstruction is now done in C++ (#784).

   * The option ppm_temp_fix = 3 was removed.  This used a
     temperature-based eigensystem for characteristic tracing but was
     never used for production science.

   * If a derived variable has multiple components, all components are now
     added to plotfiles. Previously only the first component was used. (#758)

   * We have updated our workflow when it comes to Castro's dependencies.

     Previously Castro shipped with it a minimal set of microphysics that
     allowed basic problem setups like Sedov to compile, and more advanced
     setups (like ones that include nuclear burning) required downloading
     the starkiller-astro Microphysics repository as an additional step.
     Now, that Microphysics repository is a requirement for using Castro.
     If you are a current user of the Microphysics repository and prefer
     the current workflow where you maintain Microphysics as a separate
     installation from Castro, no change in your workflow is necessary:
     if MICROPHYSICS_HOME is set as an environment variable, Castro will
     use the Microphysics installation in that directory. However we have
     also added Microphysics as a git submodule to Castro, which is now
     the required path if you previously were not using the more advanced
     microphysics (but is also a possibility for those previously using a
     standalone Microphysics installation). To obtain this, you can use
     git submodule update --init --recursive from the top-level directory
     of Castro. The developer team ensures that the version of Microphysics
     that you obtain this way is consistent with the current version of Castro.
     Then, you can keep up to date with the code mostly as normal, except now
     using git pull --recurse-submodules instead of git pull.

     Similarly, AMReX is now maintained as a git submodule rather than as an
     external standalone installation. If you use the same git submodule command
     as above, you'll obtain AMReX. As with Microphysics, you may opt to
     rely on your own installation of AMReX by setting the AMREX_HOME
     environment variable. However you are then responsible for keeping it
     in sync with Castro; if you use the submodule, then you'll get the version
     of AMReX that we have tested to ensure compatibility with the current
     version of Castro. (#651, #760, #762, #765)

   * The names of the conserved state variables in C++ (Density, Xmom, etc.)
     have been changed to match the names in Fortran (URHO, UMX, etc.).
     For user code, this will only affect problem-specific setup code
     like Prob.cpp that references specific state variables. For compatibility,
     we have kept a copy of the old names around that redirect to the
     new names, but the old names are now considered deprecated and will
     be removed in a future release. (#757)

# 20.02

   * Fixed a bug in the nuclear burning timestep estimator when on GPUs
     (#745)

   * rewrote the 4th order SDC hydro driver in C++ to allow code reuse
     with other solvers (#742), and simplified the 2nd order SDC code
     to do dimensional sweeps to reduce memory (#749)

   * The option radiation.integrate_planck has been removed; it was only
     used by one test. By default we always do the full integral of the
     Planck function. (#740)

   * Most of the radiation test problems have been moved over to a new
     opacity directory, rad_power_law, and all of the parameters that
     controlled the behavior of the power law opacity have been moved
     to the extern probin module. We now always expect you to pick a
     specific opacity implementation, so the parameter
     radiation.use_opacity_table_module has been removed. The "null"
     opacity implementation has been previously moved, and the code
     will fail to compile if you attempt to use it; you will need to
     update to rad_power_law. (See the documentation for information
     about how to use this new implementation.)

     Additionally, the code for the multigroup solver was effectively
     previously setting the Rosseland opacity, kappa_r, equal to the
     Planck opacity, kappa_p, if the latter was set but the former was
     not. There was similar unintuitive behavior for the behavior of
     the scattering parameter. Now you will get exactly what you ask
     for in the probin file, given the defaults in the _parameters file
     for the rad_power_law opacity. By default the constant coefficients
     for both are negative, which is invalid, so both must be set to a
     non-negative value for the code to work. Problems that were previously
     setting const_kappa_p but not const_kappa_r should set the latter
     equal to the former to maintain the same code behavior. The analogous
     thing should be done for the exponents (kappa_p_exp_m, kappa_p_exp_n,
     and kappa_p_exp_p). (#725)

   * The parameter radiation.do_real_eos = 0 has been removed, and its
     functionality is now enabled with a new equation of state called
     rad_power_law. This new EOS is only compatible with the pure
     radiation-diffusion tests, not with castro.do_hydro = 1. (#722)

   * We now default to use_retry = 1, instructing Castro to retry a
     step with a smaller dt if there is a CFL violation, burning
     failure, or negative timestep.  For the burning failure, we have
     Castro set the Microphysics parameter abort_on_failure to .false.
     at a high priority (so it overrides the Microphysics default).
     We also check to make sure the combination of parameters makes
     sense at runtime. (#724)

   * The parameter castro.hard_cfl_limit has been removed. (#723)

   * Some unnecessary clean_state calls were removed (#721)

   * Support for neutrino radiation diffusion has been removed.

   * A bug was fixed in the hydro CFL timestep estimator for
     simplified-SDC.  The timestep was more restrictive than it needed
     to be. (#727)

   * A bug was fixed in the simplified-SDC nuclear burning timestep
     estimator (#733)

# 20.01

   * A new option castro.limit_fluxes_on_large_vel has been added. It
     is similar to the existing option limit_fluxes_on_small_dens --
     fluxes are limited to prevent the velocity in any zone from
     getting too high. The largest legal speed is set by
     castro.speed_limit. (#712) This is more general than the previous
     solution proposed by castro.riemann_speed_limit, so that
     parameter has been removed. (#714)

   * The AMR parameter amr.compute_new_dt_on_regrid is now on by
     default. This avoids crashes that result from the CFL number
     being too large after regridding, because we update the
     timestep after seeing that larger velocity. You can still opt
     to set this off if you want to in your inputs file. (#720)

   * We have added calls into Hypre that only exist as of version
     2.15.0, so that is the new minimum requirement for Castro
     radiation. Note that Hypre is now hosted on GitHub at
     https://github.com/hypre-space/hypre.

   * A new option castro.limit_fluxes_on_large_vel has been added. It
     is similar to the existing option limit_fluxes_on_small_dens --
     fluxes are limited to prevent the velocity in any zone from
     getting too high. The largest legal speed is set by
     castro.riemann_speed_limit. (#712)

   * A new option castro.apply_sources_consecutively has been
     added. By default we add all source terms together at once. This
     option, if enabled, adds the sources one at a time, so that each
     source sees the effect of the previously added sources. This can
     matter, as an example, for the sponge source term, which may be
     more effective if it is added after source terms such as gravity
     that update the velocity. (#710)

   * A new option castro.ext_src_implicit has been added. The external
     source terms were previously only implemented as an explicit
     predictor-corrector scheme. The new option, if turned on, changes
     the handling of the external source terms to allow an implicit
     solve. This is done by subtracting the full old-time source and
     adding the full new-time source in the corrector, rather than
     -0.5 and +0.5 of each, respectively. It is still up to the
     individual problem to make sure it is consistent with this scheme
     if the option is turned on. (#709)

   * Add option for using monopole BCs in 3D.  By setting
     gravity.max_multipole_order to a negative number, you can use
     monopole gravity to fill the boundary conditions, rather than the
     multiple BCs. This is useful for debugging purposes.  To make the
     behavior consistent, we now use multipole BCs by default in 2D as
     well. (#716)


# 19.12

   * The use_retry mechanism has been enabled for the simplified
     SDC time integration method. (#695)

   * A case where use_retry could result in a very small last
     subcycle has been avoided. (#701)

   * We no longer allocate memory for sources for the species
     in the conserved state unless PRIM_SPECIES_HAVE_SOURCES is set
     (#699)

   * A subroutine eos_on_host has been added to the EOS module.
     This is a wrapper for the EOS that must be used for CUDA
     builds if the EOS is being called in probinit or other
     places that don't run on the GPU. (#693)

   * We now use VODE90 instead of VODE by default. (#677)

   * A new unit test was added, model_burner, which reads in a 1-d
     initial model and calls the reaction network on it.  This can
     be used to test tolerances, etc.

# 19.11

   * The density flux limiter was simplified and fixes a race condition
     (#646)

   * The SDC algorithm can now use Radau quadrature instead of
     Gauss-Lobatto quadrature. (#666)

   * The option castro.ppm_reference_eigenvectors has been removed.  This
     is now used by default with the CTU PPM solver.

# 19.10

   * The SDC algorithm now implements the burning conditionals
     depending on rho and T (react_rho_min, react_rho_max,
     react_T_min, react_T_max) (#598, #654)

   * The SDC/MOL PLM reconstruction now implements reflecting BCs on
     the interface states (#652, #654)

   * A well-balanced scheme has been added to the piecewise linear SDC
     method, enabled with castro.plm_well_balanaced=1.  At the moment
     it only supports constant gravity. (#294, $654))

   * The weighting of the time-node fluxes stored in the flux registers
     for SDC has been fixed (#654, #658)

   * As before, we can choose the reconstruction with PLM using the
     castro.plm_iorder flag: 1 = piecewise constant, 2 = piecewise
     linear slopes.  Now we added a way to specify the limiter used
     with the linear slopes.  castro.plm_limiter = 1 will use the 2nd
     order MC limiter and castro.plm_limiter = 2 will use the default
     4th order MC limiter (previously there was no way to select the
     2nd order limiter). (#654)

   * The Runge-Kutta based method-of-lines integration method has
     been removed in favor of the SDC integration. (#657)

   * A new way of specifying the problem runtime parameters has been
     introduced.  They can now be specified in a plain text file,
     _prob_params, and at compile time, the probdata_module is
     automatically created.  This automates the creation of the
     probdata variables, the namelist for reading them, setting them
     as managed for CUDA, and adds the ability to output the values to
     a file (like job_info).  This feature is opt-in. You need to set
     USE_PROB_PARAMS in your GNUmakefile and then define the problem
     parameters in a file _prob_params in the problem directory.
     (#234, #619, #673)

   * The time to output is now stored in the job_info file (#365)

   * The SDC time advancement method has been documented

   * The job_info file now reports the number of GPUs being used.

# 19.09

   * You can now type ./Castro.gnu.ex  *describe to see the list of
     modules / compilers the code was built with (#660)

   * The reaction quantities are now computed as proper 4th order
     averages for the plotfile, when using sdc_order = 4 (#647)

   * The velerr tagging now takes the abs() of the velocity component
     to ensure we tag large positive and negative velocities.


# 19.08.1

   * Fix CUDA compilation

   * Remove special treatment of 4th order outflow BCs (see #648)

# 19.08

   * We slightly changed how the characteristic tracing is done for
     the CTU PPM hydro solver  * we now use the limit of the parabola
     as the edge state if the wave is not moving toward the interface
     (#632)

   * The CTU PPM solver now uses a lot less memory by computing the
     integrals under the parabolas as needed instead of precomputing
     and storing them (#624)

   * We created a new error wrapper, castro_error(), to replace the
     AMReX amrex_error().  This will allow us to deal with error when
     on the GPU.

   * The new SDC solver has had substantial improvements:

      * Explicit thermal diffusion is now implemented for both 2nd and
        4th order accurate solvers. (#610)

      * There is a new option (castro.sdc_extra) for taking extra SDC
        iterations.

      * The Newton solver for the SDC update can now subcycle.

      * The sdc_solver_relax_factor option was fixed

      * There is now an absolute tolerance on species used in the
        error check for the SDC solve.

      * Some scalings of terms in the Jacobian were fixed.

      * We now use one-sided stencils for the reconstruction at
        physical boundaries (#633)

   * The flame problem setup now can initialize conservatively from a
     pre-existing flame solution.  The analysis routines have seen
     some improvements for working with 4th order accurate
     simulations.

   * the diffusion_test unit test now works for 4th order problems.

# 19.07

   * There is now a separate set of boundary filling routines for
     source terms, source_fill.F90.  Previously this was handled
     by generic_fill.F90 (#627)

   * The tagging routines have been reformulated so that they run on
     the GPU.  Since the tags_and_untags routine in AMReX is not
     GPU-accelerated, we opt to directly fill in the TagBoxArray in
     the tagging routines. We now pass the TagBox to Fortran as an
     int8_t. This means that the interface to problem_tagging_nd.F90
     has been updated to use integer(1).

     Castro_prob_err_list.H and other related files have been deleted
     as they were not actually used anywhere.

     Castro_error.cpp is now removed and there is no further support
     for writing custom tagging routines. The set of variables that we
     check for tagging is hard-coded in Castro and can be expanded as
     needed. Problem-specific tagging should be done through the
     set_problem_tags functionality. (#611)

   * The dimension-specific code for problem initialization, boundary
     conditions, and external heat terms has been removed, as warned
     in the previous release notice.

# 19.06

   * Deprecation notice: as of the 19.06 release, dimension-specific
     problem setups are deprecated. Presently they are opt-in by
     adding DIMENSION_AGNOSTIC = TRUE to your makefile, and using
     a Prob_nd.F90 file instead of a Prob_{1,2,3}d.F90 file. The
     dimension-agnostic Prob_nd.F90 is used to fill initial data
     for all dimensions. There is always a 3D loop over dimensions,
     and in 1D or 2D the unused dimensions have indices (lo, hi) =
     (0, 0) which is valid in Fortran. The current interface is found
     in Source/problems/Prob_nd.F90. Most of the problems have been
     converted to dimension-agnostic form and any remaining ones will
     be done shortly, so you can use e.g. the Sedov or DustCollapse
     problems to model you own problem on. The dimension agnostic
     problem setup also implies dimension agnostic helper routines
     such as in bc_fill_nd.F90  * any user-facing file ending in a
     1d/2d/3d.f90 file extension is deprecated. In the 19.07 release
     support for these will be removed and problem setups will only
     work if they are dimension agnostic. Please file an issue if you
     need assistance converting your problem.

   * Deprecation notice: as of the 19.06 release, problem-specific
     overrides of Castro_error.cpp, and in general custom tagging
     routines (including Castro_prob_err.cpp and associated files),
     are deprecated. The only supported mechanism for problem-specific
     tagging is through the set_problem_tags function in problem_tagging_nd.F90.
     (There are also dimension-specific versions of this file, but these
     are now deprecated as above.) Please file an issue if you need
     assistance converting your tagging setup to use the problem tagging,
     or if you need more data in that interface to be able to implement
     your tagging scheme. Support will be removed in the 19.07 release.

   * Deprecation notice: as of the 19.06 release, the problem_pre_tagging_hook
     and problem_post_tagging_hook are deprecated. These were not actually
     being used in any problem. These will be removed in the 19.07 release.

# 19.05

   * The dimension agnostic version of the external source term in
     ext_src_nd.F90 has been updated to use the ISO C binding interface,
     and two parameters, time and dt, are now passed by value  * see
     Source/sources/ext_src_nd.F90 for the new interface.

   * problem_derive_nd.f90 has been renamed to problem_derive_nd.F90.

   * The velocity calculated for the interface in the Riemann solve for
     the CGF/CG Riemann solvers can no longer exceed the speed of light.
     The parameter castro.riemann_speed_limit can be set to control the
     speed limit applied in the Riemann solver  * this is useful for
     preventing unphysically large velocities from being created at
     shock fronts or near large density gradients.

   * The algorithm for recalculating source terms after an AMR reflux
     did not set some data needed to correctly calculate positions on
     the fine levels, so position-dependent source terms (rotation,
     hybrid momentum) were being calculated incorrectly. This caused
     some strange effects in AMR simulations with rotation, such as
     binary star orbits getting wider with time and drifting relative
     to the system center of mass. This has now been fixed. (#599)

   * density_reset_method == 3 (which, in the event of a density reset,
     reset to the density at the beginning of the timestep) no longer
     exists. (#538)

   * The behavior of use_retry = 1 with retry_neg_dens_factor changed
     slightly since we now choose the retry timestep based on the difference
     between the (incorrect) negative density and the density it was
     reset to, rather than the old density at the beginning of the step.
     It still does a similar type of timestep limiting, but quantitatively
     the timesteps it chooses will be different. (#538)

   * A sign error was fixed in the hybrid_hydro angular momentum
     algorithm This addresses issue #462.  During commmit 0f09693, a
     change in signs was introduced in add_hybrid_momentum_sources,
     which should be analogous to linear_to_hybrid (#594)

# 19.04

   * The runtime parameter castro.fix_mass_flux has been removed: it is not
     clear what the use case is, and it had no test suite coverage. (#572)

   * Fixed a bug introduced in August 2015 that resulted in incorrect
     real bounds being passed to ca_initdata after a restart for problems
     using a grown domain. This would have resulted in incorrect initialization
     for problems using the grown restart capability if their initialization
     depended on the position on the grid. (#566)

   * Using point-mass gravity no longer requires USE_POINTMASS = TRUE
     in your makefile; USE_GRAV = TRUE is sufficient. However, to
     compensate for this, you must now include castro.use_point_mass = 1
     in your inputs file to enable the point mass. This input parameter
     already existed, but was defaulted to 1 since it only mattered
     if the compile flag was enabled. Now the default is 0.

   * Also, a couple bugs in the point-mass gravity have been fixed.
     The algorithm was not correct in 1D and 2D, and this has been
     resolved. And the point mass value was not being preserved across
     restarts, which is an issue if you're using point_mass_fix_solution
     to update the point mass value as mass accretes to the center of
     the domain. This has been fixed as well.

   * fixed a bug in the source term to (rho e) evolution when using
     MOL or the new SDC integration (#543, #545) and also no longer
     recompute the source terms after reflux for these methods (#549)

   * Dimension agnostic problem setups have had the interface to the
     physical boundary conditions changed (hypfill, denfill, etc.).
     If your problem is dimension agnostic, please consult the new
     interfaces in Source/problems/bc_fill_nd.F90 to understand how
     to convert your problem. The changes are that (1) the "ca_" prefixes
     have been removed from the subroutine names, (2) the "time" parameter
     is passed by value, and (3) the (lo, hi) indices that are the target
     region to update the boundaries on are explicitly passed in as
     the first two arguments. (#546)

   * removed the code to extrapolate diffusion terms to ghost cells
     as it is no longer necessary (#532)

   * we remove enthalpy diffusion (there were no known applications of
     this) and species and velocity diffusion (they were 1-d only).
     None of these routines were regularly tested.  (#534)

   * the problem diagnostics in Castro/Diagnostics have been converted to
     C++ to remain compatible with the AMReX build system.

# 19.03

   * Fixed a minor long-standing bug in the simplified SDC implementation
     involving incorrect indexing. This changes results slightly.

   * a number of tests involving reactions have been moved from
     hydro_tests to reacting_tests (#527)

   * The old spectral deferred corrections method has been renamed
     "simplified" SDC. It is accessed with time_integration_method = 3.
     This still requires building with USE_SDC = TRUE, and when building
     this way, the other time integration methods are unavailable.

   * The 4th order hydro was extended to support general equations
     of state (it is still single level only).  Artificial viscosity
     was also added.

   * A framework for a new spectral deferred corrections method was
     added.  This will allow for higher-order time integration. (#310)

   * Fix a bug where the self_heat parameter was not being initialized
     for the burning timestep limiter (#521).

   * By default, we no longer allocation storage for source terms to
     species in the primitive variable state.  This is set via the
     _variables file, parsed by set_variables.py.  To allow for
     species sources, you need to set PRIM_SPECIES_HAVE_SOURCES.  This
     is done currently for SDC. (#519)

   * renamed QVAR to NQSRC to make it clear that this is the number of
     source terms for the primitive variable state.  We also fixed a
     number of places where QVAR was used instead of NQ. (#517)

   * A new runtime parameter, T_guess, was created.  This is used as
     the initial temperature guess when calling the EOS with inputs
     other than rho, T, for the initial Newton iteration. (#509)

   * The CTU hydrodynamics driver has been rewritten in C++
     (#498)

   * The input parameter castro.do_ctu has been renamed
     castro.time_integration_method. The current legal values
     are 0 (CTU, the default) and 1 (MOL).

   * fixed a bug in the ppm_temp_fix = 1 reconstruction  * we were not
     doing the initial reconstruction of temperature

# 19.02

   * The flux limiter used with the options limit_fluxes_on_low_dens
     was not implemented correctly.  This has been fixed.  (PR #493)

   * the CTU hydro solver no longer does any allocation in any of the
     support routines  * it is all done in the top-level
     driver. (#455)

   * The CTU solver now makes explicit the range of cells looped over
     in the transverse routines (#456)

   * The plotfile quantities divu and magvort were fixed in
     axisymmetric coordinates and diff_term was fixed in all
     geometries.  (#446, 448, 449, 472)

   * abar is a new derived variable (#470)

   * the job_info file now stores domain information (#451), and the
     job_info files is also stored in checkpoints now too (#450)

   * we can now refine on enuc  * the nuclear energy generation rate.
     This is controlled by the parameters (in the &tagging probin
     namespace) enucerr, and max_enucerr_lev.  We also moved the dxnuc
     tagging parameters from inputs to probin, where they are now
     named dxnuc_min (formerly castro.dxnuc), dxnuc_max, and
     max_dxnuc_lev.  (#364, #437, #473)

   * The diffusion cutoff now is a linear ramp instead of a
     discontinous cutoff.  For densities less than
     diffuse_cutoff_density_hi, the transport coefficient is scaled
     linearly until the density reaches diffuse_cutoff_density, where
     it is zero.

# 19.01.4

   * fixed the .zenodo.json

# 19.01

   * The User's Guide is now automatically built from the development
     branch using travis.

   * we now store the job_info file in the checkpoints (#450)

   * we now automatically generate Doxygen docs along with the User's
     Guide and have started adding docstrings throughout the
     code. (#461)

   * The MG solver was optimized a bit (#464)

# 18.12

   * fixed a bug in the CUDA version of the MOL integrator  * it was
     setting the interface states incorrectly.

   * removed ppm_type = 2.  This was not used for science simulations
     because it was never shown to be completely stable.  In the near
     future, the full fourth order method will be merged which will be
     a better replacement.

   * a bug was fixed in the 1-d SDC integration  * we were not
     applying the reactive source terms to the hydrodynamics interface
     prediction.

   * we now apply the source terms for PLM before the transverse
     routines.  This is more consistent with the PPM version.

   * the angular momentum in the plotfile is not computed with respect
     to the center array initialized in the probinit

   * fixed a bug in 2-d with rotation  * we were adding the source
     terms to the out-of-plane velocity twice in the prediction of the
     interface states, resulting in an overall first-order
     approximation there.

   * fixed an issue that resulted in the custom knapsack distribution
     map not generating a useful distribution map for the reactions.
     Also, fixed an edge case where this custom distribution map could
     cause a numerical overflow and code crash. (#382)

   * rewrote the reconstruction routines to eliminate temporary arrays
     to make them GPU friendly

   * the documentation was converted to sphinx

   * changed the interface to the conductivity routines so the conductivity
     is not part of eos_t (#431)


# 18.11

   * we've restructured the CTU solver to no longer use a slab
     approach in 3d.  This is in preparation for offloading the solver
     to GPUs.

   * a few minor bugs were fixed in the transverse routines of the CTU
     solver with radiation

   * simplified conductivity interface to match the eos interface by
     moving the conductivity variable in the arguments into the eos
     type.


# 18.10

   * fixed handling of external BCs (#402)

   * unified some CUDA hydro solver code with the method-of-lines code

   * offloaded gravity source terms, rotation, reactions, and sponge
     to GPUs with CUDA

   * merged the different dimensional versions of the CTU consup
     routine into a single routine (#399)

   * removed the unsafe option "allow_negative_energy"


# 18.09

   * we now only trace under sources that are non-zero, to save
     computational expense.  (#381)

   * we now update T for consistency when we reset a small internal
     energy, e (#384)

   * The parameter dual_energy_update_E_from_e has been removed,
     and the default behavior is now that the total energy (E)
     will not be changed when we reset the internal energy (e).
     This will cause changes in simulation output. (#368)

   * The probin parameter eos_input_is_constant is now true by
     default. This means that when calling the EOS in the mode
     eos_input_re, the energy will not be updated after the EOS
     call (e.g. by the Newton iteration scheme in Helmholtz).
     This will cause changes in simulation output. (#368)

   * the problem-specific runtime parameters (probin) are now
     written to the job_info file (#380)

   * we now skip the initial EOS call prior to the burn  * this
     was redundant because we already did a clean_state (#377)

   * we now support recent versions of hypre (#373)

# 18.08

   * the old use_mlmg_solver parameters were removed  * this
     has been the only multigrid solver in Castro for some time,
     so the parameters had no effect.

   * The parameter dual_energy_eta3 was removed. This had been
     introduced mostly for testing purposes and ended up being
     unnecessary. Also, the EOS call at the beginning of the burn was
     removed; this should provide a modest computational gain. Answers
     may change at the level of the EOS tolerance (#377).

   * The functionality that checks for a regrid at the end of
     the timestep will now apply all tagging criteria, not
     just the t_sound / t_enuc criterion.

   * A bug was fixed that occurred when a retry was triggered in
     the middle of a group of subcycles, rather than at the
     beginning of the advance (#358).

   * A bug with the logic of refluxing when using retries was
     fixed (#357).

   * Tagging can now be done on relative gradients, in addition
     to the existing capability for absolute gradients (#354).
     For example, tempgrad_rel is a relative gradient criterion
     that will tag a zone with temperature T if any adjacent zone
     has a temperature that is different by more than tempgrad_rel * T.
     The tagging is enabled up to a given level with the parameter
     max_tempgrad_rel_lev. The corresponding new tagging criteria
     for other fields are named similarly.

   * Retries can now be done in a dynamic fashion (#179). An
     advance is itself a subcycled advance always, and we keep
     subcycling until we are done with the step. By default we
     still do only a single timestep when use_retry is not
     enabled, but this helps in cases where we have to reject
     the timestep for (say) a CFL violation, and then during
     the retry the CFL criterion is violated again. In the past,
     we would simply have to abort the run if this happened. Now
     we can cut the timestep again and keep going. Additionally,
     if you set abort_on_false to F in your probin file's extern
     parameters, then a burn in Microphysics will not cause an
     abort of the run, and Castro now knows how to deal with that
     by doing a retry and taking a shorter timestep (under the
     logic that most burn failures come from taking too large of
     a hydrodynamic timestep for the burner to be able to keep up).

# 18.07

   * A new GPU (CUDA) hydrodynamics solver (based on the
     method-of-lines solver) has been added, based on the work
     initially done in StarLord.  This is a work in progress, and
     requires the "gpu" branch of AMReX.

   * We removed all dependencies on the AMReX F_BoxLib source, in
     preparation for this source being removed in the future.

   * we now set the number of variables at compile time, by parsing
     the _variables file and interpreting which options are set in the
     preprocessor.  Note that a side effect of this change is that the
     number of radiation groups is now set at compile time instead of
     at runtime.  This change is needed for the GPU port.

     To set the number of radiation groups, set NGROUPS=4, e.g. for
     4 groups, in your problem's GNUmakefile.  Similar options exist
     for neutrinos.

     A related change is that it is now possible to set the number
     of advected quantities (that are not species or EOS auxillary
     fields) via NUMADV in your GNUmakefile.

# 18.06

   * The new multilevel multigrid solvers (MLMG) in the AMReX
     framework are now the default for self-gravity and constructing
     the operator for explicit diffusion.

   * A new test problem, the classic double Mach reflection, was
     added.

   * fix an issue in the retry logic that could sometimes result in
     an overflow of the estimated number of subcycle steps.

   * Improved the behavior of retries when they hit CFL violations
     (#334, #335).

# 18.05

   * Gamma_1 is now a derived variable

   * a new diffusion solver is implemented that uses the new muligrid
     framework in AMReX to compute the diffusive operator.  This can be
     enabled with diffusion.use_mlmg_solver = 1.

# 18.04

   * The job_info file now indicates which runtime parameters were
     changed from their default value (#314)

   * Improvements made to the 4th order hydro solver for single-level

# 18.03

   * The option ppm_trace_sources has been removed  * we now
     always trace on source terms with ppm.  Additionally, all
     sources are now traced, not just momentum sources.

   * The method-of-lines integrator has been rewritten.  It now works
     properly with sources.  Note: it is not intended for multilevel
     yet.  (#288, #287, #286, #164, #291, #137)

# 18.02

   * The approximate state Riemann solvers have been split into two
     parts: the computation of the interface state and the evaluation
     of the fluxes from this interface state.  This gives additional
     flexibililty in using these solvers in other methods.

# 18.01

   * The parameter dtnuc_mode has been removed. This was initially used
     for testing various forms of the burning timestep limiter before a
     final form was settled on.

   * Minor inconsistencies in how the external and diffusion source terms
     were constructed when simultaneously using reactions (#268, #269)
     have been fixed (#271).

   * The deprecated parameter castro.use_colglaz is removed. It was
     deprecated in June 2016 because it was obsoleted by the parameter
     castro.riemann_solver, which can be set to 1 to use the Colella
     and Glaz Riemann solver.

   * The state variable indicies in Fortran are now all defined in a
     single file, Source/driver/_variables.  This makes it much
     clearer and consistent and will allow for autodocumentation and
     clean-ups for GPU acceleration in the future.

# 17.12

   * The sponge can now operate based on pressure. The new parameters
     belong in the sponge namelist, and are named
     sponge_lower_pressure and sponge_upper_pressure. It works on the
     same principle as the density sponge.

   * The sponge can now drive the system to a particular velocity (the
     default is still zero velocity). The new parameters belong in the
     sponge namelist in your probin file, and are named
     sponge_target_{x,y,z}_velocity.

   * The SDC_Source_Type StateData was removed, as its purpose is now
     supplanted by the change to always keep the source terms in
     StateData (see below), and it was thus redundant. This does not
     change code output but does mean that old checkpoints generated
     while using SDC are no longer compatible with the current
     code. Normally we strive to maintain forward compatibility of
     checkpoints, but this change was considered justified because SDC
     is still considered an experimental feature under development and
     to our knowledge no production science runs have yet been done
     using SDC.

   * The parameter gravity.max_solve_level was removed. This was added
     to work around convergence issues in the multigrid solve, but
     those convergence issues have since been fixed, so the parameter
     is no longer used.

   * The Source_Type StateData now holds the actual old- and new-time
     sources (previously it held the source term predictor). This
     StateData is used to fill the sources_for_hydro MultiFab which
     provides the source terms used in the hydrodynamic update. Since
     it is StateData, this operation is done with a
     FillPatch. Consequently the sources_for_hydro data has meaningful
     data in both physical domain ghost zones and ghost zones at a
     coarse-fine interface (previously it only had meaningful data on
     fully interior ghost zones). Checkpoints will now have both old
     and new data. This change does result in a difference in the
     simulation output for simulations with source terms, as the
     answer will be different at the boundaries and at coarse-fine
     interfaces. (#116, #253) A related bug when using SDC was fixed
     too. (#56)

   * The parameter castro.keep_sources_until_end has been removed.

   * Source terms (gravity, rotation, etc.) have now been all
     coalesced into a single MultiFab. This reduces the memory
     footprint of the code. This negates the need for the code
     parameters update_state_between_sources and
     coalesce_update_diagnostics, so they have been removed. This will
     cause a change in results: the previous code behavior was to
     update the state with each successive source term as it was
     applied at the new time. Now every source term will be calculated
     using the same new-time state (the one coming out of the hydro)
     and the source terms are all applied to the state in one shot at
     the end. Note that both this and the old approach are
     second-order accurate in time. For problems that combine multiple
     source terms, this will cause changes that are larger than
     roundoff.  In particular, if rotation is used in conjunction with
     other source terms, changes will be relatively large, because the
     Coriolis force depends on the velocity, so the source term is a
     bit different now that it is seeing a different velocity. (#165,
     #249)

   * As of 17.10, there is a new option castro.plot_per_is_exact. If
     this is set to 1, timesteps will be shortened to exactly hit the
     time interval specified by amr.plot_per. An issue with this
     (#242) was fixed (#243) where an incorrect timestep would be
     taken after a restart if the previous step had been shortened.

   * We can now use the new multigrid solver from AMReX (implemented
     in C++) instead of the older Fortran solver (#241).  This is
     enabled with gravity.use_mlmg_solver=1.  Note, this only works in
     3-d currently.  This has several options:

     gravity.mlmg_max_fmg_iter = 0 : This integer parameter determines
     how many FMG cycles will be performed before switching to
     V-cycle.

     gravity.mlmg_agglomeration = 0 : This boolean flag determines if
     AMR level 0 grids will be agglomerated as the grids are coarsen
     in the multi-grid hierarchy.

     gravity.mlmg_consolidation = 0 : This boolean flag determines if
     grids on an AMR fine level that is in a single-level solve or the
     lowest AMR level of a multi-level composite solve will be
     consolidated as the grids are coarsen in the multi-grid
     hierarchy. Here, consolidation means the grids will be moved to
     fewer MPI processes.

     Numerical experiments have show this scales better than the old
     solver.

   * Apply the sources to the state's ghost zones (#255).  This fixes
     #213  * it ensures the advances give a valid update for the ghost
     zones in the State_Type.


# 17.11.1

   * Minor bug fixes from the 17.11 release. There is a corresponding
     17.11.1 release of AMReX.

# 17.11

   * A bug was fixed in non-Cartesian simulations with AMR
     (1D spherical and 2D cylindrical). The bug was introduced
     around version 17.02 and resulted in incorrect synchronization
     of the pressure term in the momentum equations. The effect
     would have manifested as non-conservation of momentum or
     strange effects at coarse-fine interfaces.

   * The sponge is now always time centered. The option to
     do otherwise was introduced in 17.02, and has now been
     removed. Additionally, the form of the energy source
     term has been corrected for the time centered case,
     and brought into line with how we do the energy source
     term for other sources. (Issue #7, Issue #57)

   * Fixed a bug in the fix for #188.

   * Conductivity_dir has been renamed CONDUCTIVITY_DIR to be consistent
     with EOS_DIR and NETWORK_DIR

   * we no longer get the compositional derivatives as part of the EOS
     call.  If you need this functionality, you need to set the
     preprocessor variable (in Fortran), EXTRA_THERMO

   * you can now use a system's BLAS routines, instead of compiling
     the versions from Microphysics by setting USE_SYSTEM_BLAS=TRUE.
     This then looks at BLAS_LIBRARY for the link line.


# 17.10

   * It is sometimes useful to be able to do some sort of initialization
     phase in your simulation, stop with a checkpoint, and then restart
     (possibly with different options) to do the main phase of the run.
     In this case, you may want to reset the simulation time to zero for
     analysis purposes. The new option castro.reset_checkpoint_time allows
     you to do this: by setting it to the time you want, the checkpoint you
     generate will have this new time. Similarly, castro.reset_checkpoint_step
     allows you to reset the timestep number (for example, to 0). Both options
     only work when you're using amr.checkpoint_on_restart=1, which itself
     requires amr.regrid_on_restart=1. This option is only intended to be used
     for the case where you're generating this checkpoint, so you also need
     to temporarily set max_step and stop_time to the target values you're
     resetting them to, to prevent further steps after the restart. After you
     have the new checkpoint, then you can undo those temporary variables
     and continue your run as usual.

   * A minor error in the gravity source terms was fixed (#109).
     This error should not normally have been observable.

   * fixed a bug in the artifical viscosity in 1-d in
     non-Cartesian geometries (issue #175)

   * the README.md now describes the process to become a
     "core developer" of Castro, and what this means.

   * Network_dir has been renamed NETWORK_DIR and EOS_dir has been
     renamed EOS_DIR.  All of the problem GNUmakefiles have been
     updated.  The old names will continue to work in the near future,
     but users are encouraged to change any of their problems to use
     the new all-caps names (PR #184)

   * the density flux limiting functionality now has a small tolerance
     (#185). It has also been documented (#193).

   * the timestep retry is now conservative  * this was accomplished
     by saving the density fluxes to use in the conservative gravity
     update (#178). Also, a bug in the timestep retry for Poisson
     gravity was fixed (#188).

# 17.09

   * the Source/ directory structure has been reorganized,
     putting the source files into directories by physics and
     eliminating the Src_1d, Src_2d, ... subdirectories

   * the Riemann solvers have been merged into a single
     dimensional-agnostic version in Src_nd.  In 2-d there was an
     issue with how the Godunov state for the CG solver was stored on
     interfaces, which would affect the internal enery evolution.

   * the PLM and PPM reconstruction routines were merged into
     a single dimensional-agnostic version in hydro/

   * the characteristic tracing routines were merged into
     dimensional-agnostic versions in hydro/ and radiation/.  This
     change fixed and outstanding issue  * the PLM reconstruction in
     1-d now uses a reference state. (issue #11)

# 17.08

   * the option castro.limit_fluxes_on_small_dens now only limits
     on density as the name suggests. It originally also limited
     fluxes if the internal energy would go negative, but this
     caused problems in runs with MPI, so it was removed. It was
     not strictly needed anyway, as the normal logic for handling
     negative internal energies is reliable.

   * two errors were fixed in the implementation of the triggered
     regrid at the end of a timestep. The method now correctly conserves
     fluid variables at coarse-fine boundaries.

   * the XGRAPH stuff that output xmgr-compatible 1-d ASCII profiles

   * fixed a bug where the gravity runtime parameters were not
     being properly initialized in the Fortran side of the code.

   * the viscosity routine is now separate from conductivity
     in Microphysics/.  Also, Castro can now use the stellar
     conductivity that is part of StarKiller.

   * the StarKiller-astro Microphysics repo now uses a denser table
     for the Helmholtz EOS (thanks to Frank Timmes).  If you are using
     this EOS, the new table will be soft-linked to your build
     directory automatically.  If you have an old copy laying around,
     it might fail to run, with an I/O error.

# 17.07

   * start of some code cleaning for eventual GPU offload support
     merging from StarLord

   * added a framework for method-of-lines integration of the
     hydrodynamics.  This does not do characteristic tracing, but
     instead does reconstruction through multiple stages of an ODE
     integrator.  At the moment, radiation is not supported.

# 17.06

   * we now require the AMReX library instead of the BoxLib library

   * the Microphysics repository that we rely on for the EOS and
     reaction networks is now part of the starkiller-astro github.
     You can change your clone to refer to this via:

     git remote set-url origin ssh://git@github.com/starkiller-astro/Microphysics

   * a new mechanism for using a stability criterion to trigger
     a regrid at the end of a timestep was added (PR #122)

   * some cleaning of the logic for momentum fluxes and
     limit_hydro_fluxes_on_small_dens (issues #130, #131)

# 17.05

   * some protections added in the retry code


# 17.04

   * rewrote the conservative gravity formulation to work off of the
     potential.  This gives the best conservation with AMR.  This is
     equivalent to the description in Appendix B from Katz et
     al. (2016).


# 17.03

   * the new refluxing method introduced in 16.11 has been removed,
     as it was determined to not provide any benefit in accuracy.

   * new derived plot variables are available when using thermal
     diffusion, including the conductivity, diffusion coefficient, and
     the entire diffusion term to the energy equation.  (issue #104)

   * when all derived variables were stored in the plotfile, we were
     storing the mass fractions twice.  E.g. for he4, we were saving
     "he4" and "X(he4)".  Now only the latter is stored.

   * created a post_simulation() function that is called at the end of
     a simulation.  An example is provided by test_diffusion where we
     output the norm of the error against the analytic solution.
     (issue #107, 108)


# 17.02

   * diagnostic information about how source terms update the state
     has been overhauled and made uniform. All source terms, including
     hydro and resets, print their changes to the state in the same
     format. The parameter print_energy_diagnostics has been renamed
     print_update_diagnostics, and a new parameter
     coalesce_update_diagnostics has been added so that you can
     combine all of the old-time and new-time updates into one print.
     (issue #58)

   * to support both single and double precision, all of the floating
     point declarations use the amrex_real type defined in the
     amrex_fort_module  * this is set to single or double precision at
     compile time.  All constants also now use this type.  For
     brevity, we rename it to 'rt' in the use statement.  (issue #34)

   * the sponge is now time-centered by default (issue #7)

   * the ppm_temp_fix stuff has been documented and made
     consistent across dimensions (issue #25)

   * the job info git information is now the output of git describe,
     this gives more information, including the last tag, how far we
     are from the tag, and an abbreviated hash.  It also indicates if
     your working directory is dirty


# 17.01

   * the radiation-specific version of computeTemp has been removed
     and instead everything goes through the main Castro computeTemp.
     This affects, in particular, how we treated small internal
     energies in radiation. (issue #64)

   * the radiation-specific versions of umeth and consup have been
     merged with the pure hydro routines.  This gives round-off level
     differences.  This also now uses all velocity components in the
     kinetic energy correction for radiation. (issues #66, 70)

   * a minor bug was fixed in the 3-d radiation characteristic tracing,
     regarding which gamma_1 (gamc) is used.


# 16.12a

   * fix a restart bug with radiation that was introduced after 16.10
     (this was cherry-picked from development) (issues #76, 78)

# 16.12

   * BoxLib now requires a C++ 11 compiler by default.  As part of
     this transition, PArrays are replaced by C++ Arrays.
     Additionally, changes to the BoxLib build system mean that we
     only need to supple COMP for the compiler.  FCOMP is now
     ignored.

   * The User's Guide has been updated to reflect the current flow of
     the algorithm.

# 16.11

   * we now distinguish between gravity (which can include a
     constant gravitational acceleration) and self-gravity,
     with the GRAVITY and SELF_GRAVITY preprocessor directives

   * some work on the sync between levels was done  * this will
     be described in a forthcoming paper. The main change by default
     is that after a reflux, we recompute the value of the source terms
     on the affected levels so that the new-time source term knows about
     the updated state due to the flux. For gravity, this resembles
     what the original Castro paper described for a sync source, but this
     is now done in a consistent way for all source terms. This should be fairly
     cheap which is why it is enabled by default, but you can disable it
     (see castro.update_sources_after_reflux). An additional optional
     change is a new strategy for refluxing (see castro.reflux_strategy).
     In the existing standard method, we only reflux after all fine timesteps
     over a coarse timestep have completed. In the new method, we do a
     partial reflux at the end of each fine timestep. This means that
     the coarse state used in time interpolation for the fine level
     is slightly more accurate as we go for later fine timesteps. It
     should also be needed for self-consistently conserving energy for gravity.
     At present it is more expensive than the standard method when there
     are gravity sync solves because there are more of them, but the tradeoff
     is that the simulation is more accurate.

   * the order of computing the temperature and reseting internal
     energy was changed in a few spots. This will change results by default.

   * the radiation-specific source was moved into the Radiation/
     subdirectory

# 16.10

   * the parameter first_order_hydro has been moved from the
     radiation namespace to the castro namespace

   * the problem setups have been moved into sub-directory
     categories to make it easier to read (issue #32)

   * the way we we use tolerances in the multigrid solve for Poisson
     gravity has changed. The old behavior is that you would pass in
     gravity.ml_tol as the relative tolerance on each grid level, and
     absolute tolerances would not be used. This suffered from some
     defects, notably that on fine grids you often had to loosen the
     relative tolerance on each higher level to achieve convergence,
     and in certain cases the scheme would fail completely, for
     example if the fine grids were not covering the mass on the
     grid. We now use an input parameter gravity.abs_tol which
     controls the absolute scale of the tolerance. This can either be
     an array of values, one for each level, or a single scalar
     value. If it is the latter, then the absolute tolerance passed
     into the multigrid scheme is the tolerance multiplied by the
     maximum value of the RHS over the entire domain.  On the coarse
     grid, then, the absolute tolerance is 4*pi*G*rho_max*abs_tol, and
     on fine grids this is multiplied by ref_ratio**2. If you do not
     specify gravity.abs_tol, then a reasonable value is selected for
     the coarse level, and the same scheme is used to give it
     reasonable values on the fine levels as well. The parameter
     gravity.ml_tol has been renamed gravity.rel_tol, and has the same
     meaning as before, but it now defaults to zero. gravity.ml_tol is
     now deprecated, and will be removed in a future release. Note
     that the tolerance used in determining convergence is always the
     less restrictive of the relative and absolute tolerance
     requirements.  gravity.delta_tol has been removed. (issue #43)

   * the radiation hydro solver, that used to live in
     CastroRadiation.git has now been completely integrated into the
     main Castro git repo.  The history was preserved in the
     transition It has also been cleaned up a little (issues #24, #31,
     #33, #48)

     The radiation build variable Network_inputs was renamed
     to NETWORK_INPUTS for consistency.

     The EOSes that used to come with CastroRadiation are available
     in Microphysics.git

   * the gravity and diffusion runtime parameters have been moved
     to the master _cpp_parameters file (issue #42)

   * enthalpy, species, and temperature diffusion are now properly
     time-centered (issue #22), and a bug in the hi boundary
     inflow boundary conditions for diffusion was fixed (issue #41)

   * a flux limiter has been added that limits the size of the hydro
     fluxes if they would cause rho or (rho e) to go negative. This
     can be used with castro.limit_hydro_fluxes_on_small_dens = 1.

   * a bug for single-level problems with Poisson gravity has been
     fixed where the multi-grid tolerance was being set to an
     uninitialized value

   * a flaw in the source terms for the primitive variables in the
     hydro update has been fixed, so that source terms like gravity
     should no longer affect the pressure and (rho e) interface states
     (issue #19)

   * the prediction of source terms to the half-time used in the
     hydrodynamics reconstruction was not being done properly.
     This has been fixed (issue #18)

   * the radiation hydro ppm now implements the ppm_predict_gammae
     option

   * we no longer ship VODE or BLAS with Castro  * these are provided
     by the separate Microphysics git repo

   * the documentation of the architecture of Castro has been
     significantly improved (issues #20, #23, #29, #31)

# 16.09:

   * the PPM tracing routine for radiation was synced up with the pure
     hydro version.  In particular, it now supports ppm_trace_sources,
     implements the reference states and fixes an issue with the
     flattening.

   * The 1-d PPM routine was also updated to support tracing,
     predicting gamma_e instead of (rho e), and an inconsistency in
     the flattening was fixed.

   * the parameters ppm_reference and ppm_reference_edge_limit
     have been removed  * there was no reason to use anything other
     than the defaults

   * the parameter ppm_tau_in_tracing has been removed.  The useful
     part of this is preserved in the ppm_predict_gammae = 1
     functionality, which uses a different set of primitive variables
     (tau, u, p, gamma_e) in the prediction of the interface states.

   * the flux correction for axisymmetric and 1-d spherical coords
     has been fixed.  In particular, there is now a separate flux
     register for the pressure term that enters as a gradient in the
     momentum equation.

   * The sign on the gravitational potential has been flipped to be
     consistent with the usual convention in the physics literature,
     i.e. the potential is negative and we solve del**2 phi = 4 * pi *
     G * rho.

   * Castro_advance.cpp has been significantly cleaned up. Each source
     term (gravity, rotation, diffusion, etc.) has a MultiFab
     associated with it through which it affects the state data. This
     has changed results slightly (typically a relative change no
     larger than the 1e-4 level) due to updates in the order of
     operations and consistency in usage on ghost zones.

   * An iterative solver for coupling between reactions and
     hydrodynamics has been introduced, which you can enable with
     USE_SDC = TRUE in the makefile. The number of iterations done for
     each timestep is controlled with castro.sdc_max_iters.

   * We changed the defaults for the gravity and rotation sources.
     now we do grav_source_type and rot_source_type = 4 by default.
     This is a conservative formulation for the energy equation that
     incorporates the source as a flux in the energy equation.  See
     Katz et al. 2016 for details.

     We also do implicit_rotation_update = 1 by default  * this does a
     slightly better coupling of the Coriolis force in the momentum
     equation by doing an implicit velocity update

     We also make ppm_trace_sources = 1 the default  * this does
     parabolic reconstruction of the momentum sources and traces
     under them when constructing the interface states

   * we now set castro.cg_blend = 2 by default.  This has no effect for
     the default CGF Riemann solver, but for the Colella & Glaz solver
     (castro.riemann_solver = 1), this will augment the secant iteration
     for the pstar find with bisection if we fail to converge.  This
     makes the root find for the star state more robust.

   * a new "fake" setup, riemann_test_zone can be use to send a left /
     right hydro state to the CG Riemann solver for testing  * this acts
     as a unit test for that solver.

   * the default for castro.allow_negative_energy is now 0  * this is
     the safer choice.

   * the default for castro.point_mass_fix_solution was changed to 0
      * this is a more expected behavior for new users.


# 16.08

   * A new parameter gravity.max_multipole_moment_level was added.
     This comes into play when using the multipole solver to compute
     the boundary conditions on the domain for isolated mass
     distributions.  The default behavior in Castro when constructing
     boundary conditions for the gravity solve is to do a multipole
     expansion sum over the density on the coarse grid only.  If you
     increase the value of that new parameter from its default of 0 to
     some higher number, it will use the data from those more refined
     levels in constructing the boundary condition values.

   * The file sponge_nd.f90 in Source/Src_nd/ has been renamed to
     sponge_nd.F90, the file extension change indicating that it can
     now be run through the preprocessor. Please update your local
     name for this file if you're overriding it in your problem setup.

   * The sponge update is usually done in an implicit fashion,
     but you can now instead do an explicit update with
     castro.sponge_implicit == 0.

   * the shock variable is now output if we are running with shock
     detection enabled

   * Microphysics/eos is now Microphysics/EOS

   * a number of changes were done in the Microphysics repo  * see
     Microphysics/CHANGES for a log of those


# 16.07

   * For consistency across the BoxLib suite of astro codes, we've
     renamed the main environment variables.  CASTRO_HOME now replaces
     CASTRO_DIR; MICROPHYSICS_HOME now replaces MICROPHYSICS_DIR.

   * The EOS, network, and conductivity routines have been moved to
     sub-directories or Castro/Microphysics/.  This reflects the way
     the layout in the standalone Microphysics repo as well as that in
     Maestro.

   * Some of the routines in Source/Src_nd/ have been renamed from
     .f90 files to .F90 files so that we can use the preprocessor. If
     you were using any of them (Prob_nd.f90, problem_tagging_nd.f90,
     Rotation_frequency.f90, or ext_src_nd.f90) by having local copies
     in your problem directory that overwrote them, please be sure to
     update the file extension so that Castro will recognize them.

   * If you were using allow_negative_energy == 0, the case where
     (rho*E), the total gas energy of the zone, was negative was
     indirectly covered and it would be reset in this case due to the
     way the logic worked for resetting the internal energy and then
     updating the total energy to be consistent with it. However at
     one point we added an option castro.dual_energy_update_E_from_e
     which disabled that second update and also meant that negative
     (rho*E) was again possible. This possibility has now been
     precluded directly, by resetting (rho*E) the same way if we
     detect that it is negative. This should not change results unless
     you were using castro.dual_energy_update_E_from_e = 1. This is
     also a good time to plug the newer option
     castro.allow_small_energy, which if set to 1 will reset when you
     hit a (rho*e) that is less than the smallest possible energy for
     the (rho, small_temp, X) in that zone. Note that it requires an
     extra EOS call.

   * The default interpolation for coarse zones into fine zones is
     piecewise linear.  There is now an option to use piecewise
     constant instead  * set castro.state_interp_order to 0. Note that
     if you use piecewise linear you can set
     castro.lin_limit_state_interp to 1 if you want to preserve linear
     combinations and therefore guarantee that, say, sum(X) = 1.

   * If you set the new option castro.limit_fluxes_on_small_dens = 1,
     the fluxes will be explicitly limited such that a negative
     density is never created.

   * Along similar lines, there are also new options for how to reset
     a negative density if one should arise. Set
     castro.density_reset_method = 2 to use the average of all
     adjacent zones instead of the default, which is the
     characteristics of the adjacent zone with the highest
     density. Set it to 3 if you want to reset it to the original zone
     state before the hydro update.

   * We have fixed an issue where diffusion did not work correctly if
     add_ext_src = 0.  The diffusion source term is now independent of
     whether you have user-defined source terms.

   * ConvertCheckpoint/ now lives under Util/

   * UsersGuide/ is now Docs/  * this is consistent with the other
     BoxLib codes

   * Burning is no longer done in ghost cells for boundaries with
     neighbors on the same level of refinement.  Instead a ghost cell
     fill is done to fill the like-level neighbor cells.  As a
     consequence of this change, if reset_internal_energy() is invoked
     in a cell, to reset the internal energy to E - K, this reset is
     now reflected in the ghost cells (this is a more consistent
     behavior).  Previously, the energy was never reset in the ghost
     cells.
---
title: 'CASTRO: A Massively Parallel Compressible Astrophysics Simulation Code'
tags:
  - C++
  - Fortran90
  - convection
  - hydrodynamics
  - nuclear reactions
  - nucleosynthesis
  - abundances
  - supernovae
authors:
  - name: Ann Almgren
    orcid: 0000-0003-2103-312X
    affiliation: 1
  - name: Maria Barrios Sazo
    orcid: 0000-0002-3185-9809
    affiliation: 2
  - name: John Bell
    orcid: 0000-0002-5749-334X
    affiliation: 1
  - name: Alice Harpole
    orcid: 0000-0002-1530-781X
    affiliation: 2
  - name: Max Katz
    orcid: 0000-0003-0439-4556
    affiliation: 3
  - name: Jean Sexton
    orcid: 0000-0003-2551-1678
    affiliation: 1
  - name: Donald Willcox
    orcid: 0000-0003-2300-5165
    affiliation: 1
  - name: Weiqun Zhang
    orcid: 0000-0001-8092-1974
    affiliation: 1
  - name: Michael Zingale
    orcid: 0000-0001-8401-030X
    affiliation: "2, 4"
affiliations:
  - name: Center for Computational Sciences and Engineering, Lawrence Berkeley National Laboratory
    index: 1
  - name: Department of Physics and Astronomy, Stony Brook University
    index: 2
  - name: NVIDIA Corporation
    index: 3
  - name: Center for Computational Astrophysics, Flatiron Institute
    index: 4
date: 02 July 2020
bibliography: paper.bib
---

# Summary 
Castro is a highly parallel, adaptive mesh, multiphysics
simulation code for compressible astrophysical flows.  It has been
used to simulate different progenitor models of Type Ia supernovae,
X-ray bursts, core-collapse and electron capture supernovae, and
dynamics in exoplanets.  Together, Castro, the low Mach number code
MAESTROeX [@maestroex], and the cosmology code Nyx [@nyx] make up the
AMReX-Astrophysics Suite of open-source, adaptive mesh, performance
portable astrophysical simulation codes.

The core hydrodynamics solver in Castro [@castro] is based on the
directionally unsplit corner transport upwind method of @ctu with
piecewise parabolic reconstruction [@ppm].  Modeling reactive flows in
stellar environments is a core capability of Castro.  Astrophysical
reaction networks are stiff and require implicit integration
techniques for accurate and stable solutions.  In Castro, we have
several modes of coupling the reactions to hydro.  The simplest method
is the traditional operator splitting approach, using Strang splitting
to achieve second-order in time.  However, when the reactions are
energetic this coupling can break down, and we have two different
implementations based on spectral deferred corrections (SDC), a method
that aims to prevent the hydro and reactions from becoming decoupled.  The
simplified SDC method uses the CTU PPM hydro together with an
iterative scheme to fully couple the reactions and hydro, still to
second order [@simple_sdc].  Alternatively, we have implemented a
traditional SDC method that couples hydro and reactions to both second
and fourth-order in space and time [@castro_sdc] (at present, this
method is single-level only).  The Strang splitting and simplified SDC
methods have a retry scheme, where a timestep will be rejected and retried
at a smaller, subcycled timestep if the burning solve fails to meet its
tolerance, negative densities are generated, or we violate one of the
timestepping criteria.

In addition to reactive hydrodynamics, Castro includes full
self-gravity with isolated boundary conditions and rotation, both
implemented in an energy-conserving fashion, explicit thermal
diffusion, and gray [@castroII] and multigroup [@castroIII] flux
limited diffusion radiation hydrodynamics.  A constrained transport
MHD solver based on the CTU algorithm is also available, and can use
the same physics source terms.  Castro can use an arbitrary equation of
state and reaction network, and these microphysics routines are
provided by the StarKiller project [@starkiller].

Castro is built on the AMReX [@AMReX] adaptive mesh refinement (AMR)
library and is largely written in C++ with a few Fortran compute
kernels.  AMR levels are advanced at their own timestep (subcycling)
and jumps by factors of 2 and 4 are supported between levels.  We use
MPI to distribute AMR grids across nodes and use logical tiling with
OpenMP to divide a grid across threads for multi-core CPU machines
(exposing coarse-grained parallelism) or CUDA to spread the work across
GPU threads on GPU-based machines (fine-grained parallelism).  All of
the core physics can run on GPUs and has been shown to scale well to
thousands of GPUs [@castro_2019] and hundreds of thousands of CPU cores
[@castro_2017].  For performance portability, we use the same source code
for both CPUs and GPUs, and implement our parallel loops in an abstraction
layer provided by AMReX. An abstract parallel for loop accepts as arguments
a range of indices and the body of the loop to execute for a given index,
and the AMReX backend dispatches the work appropriately (e.g., one zone per
GPU thread). This strategy is similar to the way the Kokkos [@Kokkos] and
RAJA [@RAJA] abstraction models provide performance portability in C++.

# Statement of Need

While there are a number of astrophysical hydrodynamics simulation codes, Castro
offers a few unique features.  The original motivation for developing
Castro was to build a simulation code based on a modern,
well-supported AMR library (BoxLib which evolved to AMReX), using
unsplit integration techniques and targeting problems in nuclear
astrophysics.  The radiation solver was a key design consideration in
the early development.  The large developer community contributing to AMReX
(representing a large number of application codes across various domains)
results in Castro continually gaining optimizations for new
architectures.  As Castro evolved, we adopted a fully open development
model (as does the Enzo [@enzo] code, for example).  We pride ourselves in
making all of the science problems available in the Castro git repository as
we are developing them, and the infrastructure we use for running our problems
and writing our science papers is publicly available in the AMReX-Astro organization.
Other simulation codes, like Flash [@flash], also work with a general equation of
state and reaction network, but Castro is unique in focusing on
spectral deferred correction techniques for coupling the hydro and
reactions.  Finally, while some astrophysics codes have performance portable forks
(like K-Athena [@kathena], which uses Kokkos), Castro's current design -- which targets both
CPUs and GPUs for all solvers -- achieves performance portability as a core
design principle, avoiding the need for a fractured development model.




# Acknowledgments

The work at Stony Brook was supported by DOE/Office of Nuclear Physics
grant DE-FG02-87ER40317 and NSF award AST-1211563.  MZ acknowledges
support from the Simons Foundation.  This research was supported by
the Exascale Computing Project (17-SC-20-SC), a collaborative effort
of the U.S. Department of Energy Office of Science and the National
Nuclear Security Administration.  The work at LBNL was supported by
U.S. Department of Energy under contract No. DE-AC02-05CH11231.  We
also thank NVIDIA Corporation for the donation of a Titan X Pascal and
Titan V used in this research.  The GPU development of Castro
benefited greatly from numerous GPU hackathons arranged by OLCF.

# References

# Radiation diagnostics

This directory contains diagnostics for Castro's radiation tests.
These are:

- `gaussian_pulse`: Process a 2-d gaussian radiation pulse.
- `lgt_frnt1d`: Process a 1-d (radiative) sedov problem to produce rho, u, and
    p as a function of r, for comparison to the analytic solution.
- `rad_shock`: extract a 1-d slice of the data of a radiative shock test (all
        variables or a single variable) along the specified coordinate direction
        from a plotfile.  The plotfile can be 1-, 2-, or 3-d.
- `rad_source`: Analysis routine for the radiation source test.  Take a list of
    files and print out (rho e) and the total radiation energy density in the
    first zone as a function of time.
- `rad_sphere`: Print out the radiation quantities at a specified distance from
    the origin.  This is written for the 1-d radiating sphere problem.
- `rhd_shocktube`: Analysis routine for RHD_shocktube.

## Building

To build one of the radiation diagnostics, specify the executable name
as the target, e.g.

```
   make DIM=1 rad_sphere.ex
```

Take care that you compile with the correct dimension (i.e. set `DIM`
equal to the same value it had for the code that generated the plotfile
to be investigated).

## Running

Command line arguments must be passed to the executables to provide the plotfile(s)
to be anaylzed and problem parameters:

- `gaussian_pulse`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results, and `--xctr x`, `--yctr y` to provide the coordinates
    of the center of the domain (x,y) (if not provided, both default to 0.0).
- `lgt_frnt1d`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results.
- `rad_shock`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results, and `--idir d` for the direction d along which to take
    the slice (where d is an integer 1-3 corresponding to the x-z directions,
    and defaults to 1).
- `rad_source`: the arguments are assumed to be a list of the plotfiles to be analyzed
- `rad_sphere`:  `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-g groupfile_name` (or `--groupfile`) to provide the name of the group file
    (this is output at runtime by the `RadSphere` problem), and `-r radius` (or
    `--radius`) to provide the radius at which to print out the
    radiation quantities.
- `rhd_shocktube`: `-p plotfile_name` (or `--plotfile`) to provide the plotfile,
    `-g groupfile_name` (or `--groupfile`) to provide the name of the group file
    (this is output at runtime by the `RadSphere` problem), and
    `-s slicefile_name` (or `--slicefile`) to provide the name of the file to
    output the results.
# Sedov diagnostics

Process a Sedov problem to produce rho, u and p as a function of r,
for comparison to the analytic solution.

## Building & running

The n-dimensional diagnostic can be built by executing `make
DIM=n`. This will produce the executable `sedov_nd.exe`. To run, the
executable must be provided with a list of arguments. For all
problems, you must provide the name of the plotfile to be analyzed and
the name of the slicefile where the results shall be output:

```
./sedov_nd.ex -p plotfile_name -s slicefile_name
```

Additional arguments depend on the problem:

- **1d**: no additional arguments are required

- **2d cylindrical**: extra arguents `--xctr x` and `--yctr y` giving
  the coordinates of the domain center (x,y) can be provided (both
  default to 0.0)

- **2d spherical**: the argument `--sphr` *must* be provided to
  indicate a spherical problem, and the `--yctr y` argument can be
  provided to give the y coordinate of the domain center (defaults to
  0.0)

- **3d cylindrical**: extra arguents `--xctr x` and `--yctr y` giving
  the coordinates of the domain center (x,y) can be provided (both
  default to 0.0)

- **3d spherical**: the argument `--sphr` *must* be provided to
  indicate a spherical problem, and extra arguents `--xctr x`, `--yctr
  y` and `--zctr z` giving the coordinates of the domain center
  (x,y,z) can be provided (all default to 0.0)
# Sedov diagnostic routines

The diagnostic routines in `Diagnostics/DustCollapse` are used to construct
the data which appears in Figure 12 in the first CASTRO paper.

The analytic maximum density has been computed from the initial density
and the analytic r(t) using conservation of mass; this is hard-wired
into the `main.cpp` file.

The radius of the star in each plotfile is computed by first computing the
radial average of density, then finding the radius at which
the density equals half of the analytic maximum density.

## Building and running

Typing 'make DIM=n' will build the diagnostic routine for the n-dimensional problem.

To run the 1d problem, run
```
./dustcollapse_1d.exe plotfile1 plotfile2
```
i.e. run the executable followed by a list of the plotfiles to be analyzed.
For the 2d and 3d problems, additional argument(s) can be passed in *before* the
list of plotfiles to give the coordinates of the domain center. If these are not
provided, then they default to 0.0. For the 2d problem, this would be
```
./dustcollapse_2d.exe --xctr x --yctr y plotfile1 plotfile2
```
and similarly for 3d
```
./dustcollapse_3d.exe --xctr x --yctr y --zctr z plotfile1 plotfile2
```
For the 2d and 3d problems, it is also possible to print the profile to file
by providing the argument `--profile`. This will create the profile file
`prof.profile` in the plotfile's directory.

## Analytic solution

To make the analytic profile, type
```
gfortran analytic.f90
```
then
```
a.out > analytic.txt
```
## Plotting

To use gnuplot to make Figure 12, use `DustCollapse.gp` once you have created the files
`analytic.txt`, `results_1d.txt`, `results_2d.txt` and `results_3d.txt`
<!-- Thank you for your PR!  Please provide a descriptive title above
and fill in the following fields as best you can. -->

<!-- Note: your PR should:

    * Target the development branch
    * Follow the style conventions here:
      https://amrex-astro.github.io/Castro/docs/coding_conventions.html  -->

## PR summary

<!-- Please summarize your PR here. We will squash merge this PR, so
     this summary should be suitable as the commit message. If the
     PR addresses any issues, reference them by issue # here as well. -->

## PR motivation

<!-- Please describe here the motivation for the PR, if appropriate.
     This section will not be included in the commit message, and can
     be used to communicate with the developers about why the PR should
     be accepted. This section is optional and can be deleted. -->

## PR checklist

- [ ] test suite needs to be run on this PR
- [ ] this PR will change answers in the test suite to more than roundoff level
- [ ] all newly-added functions have docstrings as per the coding conventions
- [ ] the `CHANGES` file has been updated, if appropriate
- [ ] if appropriate, this change is described in the docs
These are some sample scripts for doing volume rendering with yt.
These were originally written with the WD merger problem in mind.

  -- vol-wd.py :

     this script does basic volume rendering using the perspective
     lens.

  -- vol-wd-spherical.py :

     this script does a spherical projection (all 4 pi steradians)
	 with the camera set at the origin.  This can be used to create
	 360 degree videos for uploading to youtube.

     Note to prepare the video for youtube, follow the instructions
     here:

     https://support.google.com/youtube/answer/6178631?hl=en

     In particular, you need to add metadata to the video, which can be
	 done with the google spatial media metadata injector:

     https://github.com/google/spatial-media/blob/master/spatialmedia/README.md


Construction of the logo:

1. create a bitmap slice of the wdmerger simulation using the script
   logo.py

2. Use gimp to resize the image to just the interactive stars (cutting
   out the axes).  Use "image -> crop to selection".

3. Open inkscape and import the image -- use "file -> import..."

4. Select the image and use "path -> trace bitmap..."

   Set it to "colors" and check the "remove background" box.

5. Select the largest element and resize the canvas using "edit -> resize page to selection"
This problem is the radiation-hydro (Sedov-Taylor like) blast wave from
Castro paper II. The default case is the one from Section 6.9, while the
"kp" variant with a very large Planck mean opacity is the one from
Section 6.10.
# RadFront

This setup is based on the Optically-Thin Streaming of Radiation front
described on the Castro paper II but with some differences:

1. Instead of defining the B.C. by Er, I'm using the corresponding value
   for Fr, assuming that for optically thin Fr=-cEr

2. There is no filled region with radiation

3. I'm using do_hydro=1 instead of having the hyperbolic update off. The medium
   has a very low density and temperature

4. I'm using the same value for Rosseland and Planck coeff.

The 1d version, has a domain size of 100 cm, whereas the 2d has both x and y of ~10^10 cm.
The opacities are changed to have similar optical depths along the domain for both setups.
In both cases it is seen that the radiation front propagates at speed of light. 

# Problem Setup

This is a pure diffusion problem (no hydro).  It uses the explicit
diffusion solver in Castro to diffuse a Gaussian thermal pulse.
Because diffusion in Castro is incorporated in the energy equation, we
are solving:

   d(ρe)/dt = ∇·(k ∇T)

where L is the Laplacian and k is the thermal conductivity.  Here we
assume that k is constant, but Castro does not require that
assumption.  For this problem, we take ρ = constant, and using a
gamma-law EOS, we have e = c_v T, so this can take the form of a
diffusion equation as:

dT/dt = k/(ρ c_v) ∇^2 T = D ∇^2 T

where D is the constant diffusion coefficient.

Because we are doing the diffusion explicitly, there is a constraint
on the timestep,

dt < 0.5 dx^2/D

For constant diffusion coefficient, there is an analytic solution: the
diffusion of a Gaussian remains Gaussian, with the amplitude
descreasing and the width increasing with time.


# Testing

The problem creates a derived variable, `analytic`, which is the
analytic solution at the time of the plotfile output.  This allows you
to compare the current solution to the analytic solution to assess the
error.

It also will use the Castro problem-specific post-simulation hooks (in
`Prob.cpp`) to output the L-inf norm of the error (numerical vs
analyic solution) at the end of the simulation.  This can be used
for convergence testing.


## 1-d spherical with AMR

This uses the 2nd order accurate predictor-corrector formulation of
diffusion that is used with the CTU hydrodynamics solver.  A test of
the diffusion in 1-d spherical coordinates, with 2 levels of
refinement can be run as:

```
./Castro2d.gnu.ex inputs.2d.sph
./Castro2d.gnu.ex inputs.2d.sph amr.n_cell=128 256
./Castro2d.gnu.ex inputs.2d.sph amr.n_cell=256 512
```

At the end, each run will report the norm of the error against the
analytic solution, giving:

```
 base resolution      L-inf error
64                  0.0003707056645
128                 9.414571162e-05
256                 2.437072009e-05
```


## SDC-4 in 1-d

A convergence test of the 4th-order SDC algorithm can be run as:

```
./Castro1d.gnu.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64
./Castro1d.gnu.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=128
./Castro1d.gnu.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256
./Castro1d.gnu.ex inputs.1d castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512
```

Note: this is Cartesian, not spherical (since we don't have spherical
implemented to 4th order).  The norm of the error output at the end as
a function of resolution is:

```
 resolution     L-inf error
 64           8.639542243e-05
128           5.84073812e-06
256           3.725743019e-07
512           2.340631822e-08
```


# Non-constant Conductivity

There is no analytic solution for non-constant conductivity, so we can
only do convergence testing by varying the resolution.  In this
manner, the error reported at the end of a run is meaningless, since
it is comparing against the analytic solution for constant
conductivity.


## Second-order predictor-corrector algorithm

The powerlaw conductivity is simply k = k0 T^ν.  To build the test with this
conductivity we do:

```
 make DIM=1 CONDUCTIVITY_DIR=powerlaw -j 20
```

Tests can be run as:

```
 ./Castro1d.gnu.ex inputs.1d.powerlaw amr.n_cell=64
 mv diffuse_plt00048 diffuse_64
 ./Castro1d.gnu.ex inputs.1d.powerlaw
 mv diffuse_plt00191 diffuse_128
 ./Castro1d.gnu.ex inputs.1d.powerlaw amr.n_cell=256
 mv diffuse_plt00761 diffuse_256
```

Then the error can be measured using the RichardsonConvergenceTest
tool in `amrex/Tools/C_util/Convergence` as:

```
RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_64 mediFile=diffuse_128 fineFile=diffuse_256
```

This gives:

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline 
density&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
xmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
ymom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
zmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
rho_E&    	 3.479414e-04 & 2.012958966 & 8.620750e-05 \\ 
rho_e&    	 3.479414e-04 & 2.012958966 & 8.620750e-05 \\ 
Temp&    	 3.479414e-04 & 2.012958966 & 8.620750e-05 \\ 
rho_X&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
```

(some bits were edited out)

e.g. we see second-order convergence in the temperature


## 4th order SDC

This is built as the previous test:

```
make DIM=1 CONDUCTIVITY_DIR=powerlaw -j 20
```

and then run as:

```
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64
mv diffuse_plt00048 diffuse_64
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4
mv diffuse_plt00190 diffuse_128
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256
mv diffuse_plt00760 diffuse_256
./Castro1d.gnu.ex inputs.1d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512
mv diffuse_plt03040 diffuse_512

RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_64 mediFile=diffuse_128 fineFile=diffuse_256 > convergence_diffusion.1d.lo.sdc4.out
RichardsonConvergenceTest1d.gnu.ex coarFile=diffuse_128 mediFile=diffuse_256 fineFile=diffuse_512 > convergence_diffusion.1d.hi.sdc4.out
```

This gives (for the lower resolution runs):

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline 
density&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
xmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
ymom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
zmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
rho_E&    	 1.111626e-05 & 3.948910124 & 7.198104e-07 \\ 
rho_e&    	 1.111626e-05 & 3.948910124 & 7.198104e-07 \\ 
Temp&    	 1.063477e-05 & 3.952987539 & 6.866892e-07 \\ 
rho_X&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
pressure&    	 7.410837e-06 & 3.948910124 & 4.798736e-07 \\ 
```

e.g. we see fourth-order convergence in the temperature


## 4th order SDC (2-d)

This is built as the previous test:

```
make DIM=2 CONDUCTIVITY_DIR=powerlaw -j 20 USE_MPI=TRUE
```

and then run as:

```
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=64 64
mv diffuse_plt00039 diffuse_2d_64
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4
mv diffuse_plt00157 diffuse_2d_128
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=256 256
mv diffuse_plt00626 diffuse_2d_256
mpiexec -n 16 ./Castro2d.gnu.MPI.ex inputs.2d.powerlaw castro.time_integration_method=2 castro.sdc_order=4 amr.n_cell=512 512
mv diffuse_plt02504 diffuse_2d_512

RichardsonConvergenceTest2d.gnu.ex coarFile=diffuse_2d_64 mediFile=diffuse_2d_128 fineFile=diffuse_2d_256 > convergence_diffusion.2d.lo.sdc4.out
RichardsonConvergenceTest2d.gnu.ex coarFile=diffuse_2d_128 mediFile=diffuse_2d_256 fineFile=diffuse_2d_512 > convergence_diffusion.2d.hi.sdc4.out
```

This gives (for the lower resolution runs):

```
Level  L1 norm of Error in Each Component
-----------------------------------------------
Warning: BoxArray lengths are not the same at level 0
  0    Level  L1 norm of Error in Each Component
-----------------------------------------------
  0    \begin{table}[p]
\begin{center}
\begin{tabular}{|cccc|} \hline
Variable & $e_{4h \rightarrow 2h}$ & Order & $e_{2h \rightarrow h}$\\
\hline 
density&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
xmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
ymom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
zmom&    	 0.000000e+00 & ------------ &0.000000e+00 \\ 
rho_E&    	 1.902161e-06 & 3.957610923 & 1.224299e-07 \\ 
rho_e&    	 1.902161e-06 & 3.957610923 & 1.224299e-07 \\ 
Temp&    	 1.770452e-06 & 3.966033724 & 1.132894e-07 \\ 
```

e.g. we see fourth-order convergence in the temperature




# model_burner

This takes an initial model, reads it via the model_parser, and then simply
calls the burner on all zones in the initial model.
# hse_convergence

This is meant to be a simple 1-d test for assessing the convergence of
hydro + gravity in maintaining HSE.  Convergence can be measure either
with the RichardsonConvergenceTest tool or by looking at the max |U|
in the plotfiles.

To run this problem, use one of the convergence scripts:

  * ``convergence_plm.sh`` :

    this runs CTU + PLM using the default HSE BCs and default
    use_pslope, then with reflect BCs, then without use_pslope, and
    finally runs with reflect instead of HSE BCs.

    These tests show that the best results come from HSE BCs + reflect vel

  * convergence_ppm.sh :

    this runs CTU + PPM in a similar set of configurations as PLM above
    (with one additional one: grav_source_type = 4)

    These tests show that the best results come from HSE BCs + reflect vel

  * convergence_sdc.sh :

    this uses the TRUE_SDC integation, first with SDC-2 + PLM  and reflecting BCs,
    the SDC-2 + PPM and reflecting BCs, then the same but HSE BCs, and finally
    SDC-4 + reflect

    These tests show that the PLM + reflect (which uses the
    well-balanced use_pslope) and the SDC-4 + reflect give the lowest
    errors and expected (or better) convergence:


# Evrard Collapse problem

This is the collapse of an isothermal spherical gas cloud.  This
problem was originally discussed in Evrard 1988.

This implementation of the test comes from section 9.1 of Springel, V,
Monthly Notices of the Royal Astronomical Society, Volume 401, Issue
2, pp. 791-851 (section 9.1)
# hse_convergence_general

This is an HSE convergence test problem based on flame_wave.  It is
meant to be run in 1-d, just to understand how well we preserve HSE
for this problem.



# reacting_convergence

This is a reaction test used for measuring convergence in reacting hydro,
built from the acoustic_pulse_general test.

# convergence testing

* in `Castro/Exec/reacting_tests/reacting_convergence`:
    ```
    make COMP=intel -j 4
    ```

* in output directory:
    ```
    mkdir reacting_convergence
    cd reacting_convergence
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/inputs.* .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/probin .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/helm_table.dat .
    cp ~/Castro/Exec/reacting_tests/reacting_convergence/job_scripts/cori.reacting_convergence.slurm .
    sbatch cori.reacting_convergence.slurm
    ```

* to do the analysis, in the output directory:
    ```
    cp ~/development/Castro//Exec/reacting_tests/reacting_convergence/analysis/check_convergence.sh .
    cp ~/development/Castro//Exec/reacting_tests/reacting_convergence/analysis/create_pretty_tables.py .
    ```

  edit it to use bash and give the proper executable name.  Then:
    ```
    ./check_convergence.sh
    ```
  This outputs 3 sets of convergence data files.  To process these
  into a set of LaTeX table rows::
    ```
    python3 create_pretty_tables.py
    ```
# bubble_convergence

This problem was developed to test the convergence of the 4th order
SDC solver with gravity, reactions, and reflecting boundary
conditions.

This was featured in the Castro SDC paper:
https://ui.adsabs.harvard.edu/abs/2019ApJ...886..105Z


## Measuring convergence

The script `converge_test.sh` will test the convergence for the bubble
rising problem by running tests at 5 different resolutions (the 32**2
resolution was not shown in the paper).

* Build as:
  ```
  make USE_TRUE_SDC=TRUE -j 20
  ```

* Run the tests via: 
  ```
  ./converge_test.sh
   ```

  This relies on having built the RichardsonConvergenceTest too in
  `amrex/Tools/C_util/Convergence/` and putting the executable in your
  path.

  3 output files will be produced: `sdc_converge.lo.out`,
  `sdc_converge.mid.out`, and `sdc_converge.hi.out`, giving the order
  of convergence for the 32-64-128, 64-128-256, and 128-256-512 runs
  respectively.

* Create a summary table using `create_pretty_tables.py`.  E.g., for
  the mid and hi resolution cases, you would do:
  ```
  python create_pretty_tables.py --simple sdc_converge.mid.out sdc_converge.hi.out
  ```

  This gives:
  ```
  $\rho$                           3.590942e+15   3.264       3.739197e+14   3.713       2.852085e+13
  $\rho u$                          1.11992e+24   3.794       8.071755e+22   3.930       5.296302e+21
  $\rho v$                          1.31454e+24   3.544       1.127087e+23   3.839       7.878383e+21
  $\rho E$                         3.701644e+32   2.947       4.801223e+31   3.646       3.834012e+30
  $\rho e$                         3.701199e+32   2.947       4.800721e+31   3.646       3.833886e+30
  $T$                              1.438083e+18   3.508       1.264438e+17   3.829       8.898941e+15
  $\rho X(\isotm{He}{4})$          3.589696e+15   3.266       3.732104e+14   3.711       2.849884e+13
  $\rho X(\isotm{C}{12})$          1.519871e+13   2.544       2.605906e+12   3.797       1.874147e+11
  $\rho X(\isotm{O}{16})$              35896250   3.262            3741870   3.714           285088.2
  $\rho X(\isotm{Fe}{56})$             35908410   3.264            3739049   3.713           285208.5
  ```

  which is essentially the same values shown in table 11 of the Castro SDC paper.

## HSE convergence

This setup can also be used just to test HSE convergence.  

* Build as:
  ```
  make USE_TRUE_SDC=TRUE USE_REACT=FALSE -j 20
  ```
  Then run the tests as:
  ```
  ./converge_test_sdc4_nopert.sh
  ```

  We can get the maximum velocity as:
  ```
  fextrema.gnu.ex bubble_64_plt00667 | grep -i magvel
  fextrema.gnu.ex bubble_128_plt01334 | grep -i magvel
  fextrema.gnu.ex bubble_256_plt02667 | grep -i magvel
  ```

  This gives (for max |U|):
  ```
   64: 0.018842065993
  128: 0.0011912049874
  256: 8.8777896776e-05
  ```
  demonstrating nearly 4th order convergence of the HSE state.


## Original Cori convergence testing

The results in the Castro SDC paper were run on Cori, using
the file in `job_scripts/`.

* in `Castro/Exec/reacting_tests/bubble_convergence`
    ```
    make USE_TRUE_SDC=TRUE COMP=intel -j 4
    ```

* in output directory:
    ```
    mkdir bubble_convergence
    cd bubble_convergence
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/helm_table.dat .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/inputs_2d.* .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/probin.* .
    cp ~/Castro/Exec/reacting_tests/bubble_convergence/job_scripts/cori.bubble_convergence.slurm .
    sbatch cori.bubble_convergence.slurm
    ```

* to measure convergence:

  in the `bubble_convergence` output directory:
    ```
    cp ~/development/Castro//Exec/reacting_tests/bubble_convergence/job_scripts/check_convergence.sh .
    cp ~/development/Castro//Exec/reacting_tests/bubble_convergence/job_scripts/create_pretty_tables.py .
    ```

  edit it to use bash and give the proper executable name.  Then:
    ```
    ./check_convergence.sh
    ```
  This outputs 3 sets of convergence data files.  To process these
  into a set of LaTeX table rows:
    ```
    python3 create_pretty_tables.py
    ```
# nse_test

This is a simple test problem designed to explore how well hydro and
reactions are coupled when a system enters NSE.

This version is based on ``reacting_convergence`` (which is in turn
based on ``acoustic_pulse_general``), but using the ``aprox19``
network with the NSE table enabled.

You can run the Strang convergence test with the script convergence_strang.sh

The script create_pretty_tables.py will take the 2 output files and
make a single LaTeX-formatted table of the results.
# Problem Description

This is a single mode Rayleigh-Taylor instability problem.  The
instability is triggered by perturbing the interface between the dense
and light fluid.

This problem was used in the Castro I paper to compare different PPM
types.


# Particles

This can be run with particles by building as:

```
make USE_PARTICLES=TRUE
```

and running with the `inputs_2d.particles` inputs file.  This will
assign particles according to the ``particle_file`` locations.

These particles are moved with the fluid velocity each time step, and
if they leave the domain then they just "disappear."

This test is primarily designed to test the ability of Castro to
handle particles.
# acoustic_pulse

This is the acoustic pulse problem from McCorquodale & Colella 2011.
We use this to measure convergence.

Note: for DIM > 1, you can still run this in a "1-d" mode by
setting the inputs parameter problem.init_as_1d to 1, 2, or 3, to just
do the pulse in that coordinate direction.


# Convergence testing

Here we detail the procedure used for convergence testing at NERSC and
processing the results for the pretty table included in the Castro SDC
paper.

  in `Castro/Exec/hydro_tests/acoustic_pulse`
    ```
    make COMP=intel -j 4
    ```

  in the output directory:
    ```
    mkdir acoustic_pulse
    cd acoustic_pulse
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/inputs.2d.* .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse/job_scripts/cori.acoustic_pulse_convergence.slurm .
    sbatch cori.acoustic_pulse_convergence.slurm
    ```

  to do the analysis, in the output directory:
    ```
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse/job_scripts/check_convergence.sh .
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse/job_scripts/create_pretty_tables.py .

  you may need to edit `check_convergence.sh` to use bash and have a
  different executable name.  Then
    ```
    ./check_convergence.sh
    ```
  This outputs `convergence.2d.lo.sdc4.out` and
 `convergence.2d.hi.sdc4.out`.

  To process these into a set of LaTeX table rows:
    ```
    python3 ./create_pretty_tables.py
    ```
# Sod

This is a one-dimensional shock tube problem that is used to run the
classic Sod problem.  inputs files are also provided to run
most of the test problems defined in Toro's text.

Plotting against the analytic solution can be done using the scripts
and data in `Verification/`.
The Verification directory contains the analytic solutions for the
Sedov problem, as computed using Frank Timmes' exact Sedov routine,
http://www.cococubed.com/codes/sedov/sedov3.f

To compare Castro output, we need to radially bin the plotfile data to
get density, total velocity, and pressure as a function of radius.
This can be accomplished using fsedov*.f90 routines in
Castro/Diagnostics/Sedov

For the 2-d Cylindrical explosion in Cartesian coordinates:

  * run:
   ```
   ./Castro2d.gnu.MPI.ex inputs.2d.cyl_in_cartcoords
   ```

  * create the angle-averaged profile:

    compile Castro/Diagnostics/Sedov/fsedov2d_cyl_in_cartcoords.f90 as:
    ```
    make programs=fsedov2d_cyl_in_cartcoords
    ```
    and then run it as:
    ```
    ./fsedov2d_cyl_in_cartcoords.Linux.gfortran.exe -p ../plt00144 -s sedov2d.out
    ```

  * plot:

    A Gnuplot script to plot the data overtop the analytic solution is
    provided, simply execute
    ```
    gnuplot sedov2d.gp
    ```
    in the Sedov/Verification/ directory (where the above slices should have been
    output).


When running inputs.2d.sph_in_cylcoords or inputs.3d.sph, use sedov.gp instead
of sedov2d.gp

NOTES:

  * To get good agreement with the analytic solution for the Sedov
    problem, it is necessary to start off with a small timestep.  This
    is accomplished by setting castro.init_shrink = 0.1 in the inputs
    file.

  * Subsampling is critical to ensure that the perturbation is as
    spherical as possible.

  * It is important to use a lot of 'steps' in the sedov3.f routine to
    get a good analytic solution.
# acoustic_pulse

This is an alternate version of the acoustic pulse problem from
McCorquodale & Colella 2011.  In this version, we perturb the pressure
and use a constant entropy via the EOS to find the density.

# Convergence testing

Here we detail the procedure used for convergence testing at NERSC and
processing the results for the pretty table included in the Castro SDC
paper.

  in `Castro/Exec/hydro_tests/acoustic_pulse_general`
    ```
    make COMP=intel -j 4
    ```

  in the output directory:
    ```
    mkdir acoustic_pulse_general
    cd acoustic_pulse_general
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/Castro2d.intel.haswell.MPI.ex .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/inputs.2d.* .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/helm_table.dat .
    cp ~/Castro/Exec/hydro_tests/acoustic_pulse_general/job_scripts/cori.acoustic_pulse_convergence.slurm .
    sbatch cori.acoustic_pulse_convergence.slurm
    ```

  to do the analysis, in the output directory:
    ```
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse_general/job_scripts/check_convergence.sh .
    cp ~/Castro//Exec/hydro_tests/acoustic_pulse_general/job_scripts/create_pretty_tables.py .

  you may need to edit `check_convergence.sh` to use bash and have a
  different executable name.  Then
    ```
    ./check_convergence.sh
    ```
  This outputs `convergence.2d.lo.sdc4.out` and
 `convergence.2d.hi.sdc4.out`.

  To process these into a set of LaTeX table rows:
    ```
    python3 ./create_pretty_tables.py
    ```


# sod_stellar

This problem setup was used in Zingale & Katz (2015) to test the
hydrodynamics with a general equation of state.  These problems are
pure shock tube problems with a general equation of state.  An exact
solution produced for comparison using the exact Riemann solver in
`Castro/Util/exact_riemann/`.

This test checks the initial relaxation phase. We set a relatively large damping
factor and run until the relaxation is turned off (at approximately 12s).
# planet
Model of Hot Jupiter planet to be used in conjunction with Castro(Radiation) code (2-Dimensions)

Actual model parameters defined and calculated in the HotJupiter.cpp file. Outputs to newmodelcpp.hse file for Castro to read in. Planetary model based on paper: 
  Youdin, Mitchell 2010 (The Mechanical Greenhouse: Burial of Heat by Turbulence in Hot Jupiter Atmospheres)

Inputs to Castro are in inputs_2d files

more editable parameters (velocities and vortices) in probin file.

# `nse_detonations`

This directory contains a set of scripts that setup, run, and analyze
a set of detonations, comparing Strang splitting to simplified SDC.

The basic usage is:

#. in a `screen` session of similar do:
   ```
   ./setup_runs.py
   ```
   this will create the run directories, copy the needed files, and
   run the jobs in parallel using the python Pool mechanism.

#. while the jobs are running, you can check the status by doing:
   ```
   ./show_status.py
   ```

#. once the runs are finised, you can make the suite of plots by
   ```
   ./make_plots.py
   ```


# flame_wave

This is the XRB flame setup, which was used in Eiden et al. 2020
(https://ui.adsabs.harvard.edu/abs/2020ApJ...894....6E).  That
original paper used "boosted" flames, but we no longer boost the flame
speed.  You can find the inputs file for boosted flame in Castro
versions 21.01 and earlier.

The current set of inputs files are:

* inputs_He:

  pure He flames.  This includes the inputs files for the follow-on
  paper, Harpole et al. (in prep).

  a 3-d setup is also there

* inputs_H_He:

  mixed H/He flames.  These are intended to be run with the rprox
  network.

These scripts are for Cori.  The script cori.xfer.slurm can be used to
migrate files to HPSS.
The process.xrb script can be run via `screen` on the data transfer
nodes (`dtn`) to archive data to HPSS.# Flame Wave Analysis Scripts

### plot_generator.py

Script for generating plots of a sequence of datasets using yt. The variable to plot, whether to use
logscale, the domain and colorbar bounds, streamline options, etc. are all configurable through
command line arguments. For a usage description and a full list of valid parameters, type
`./plot_generator.py -h`. *TODO*: Add output dpi and figure size settings.

### image_animator.py

Uses ffmpeg directly (through an `os.system` call) to generate an animation from a sequence of
images. Allows the user to specify an output file and offers some support for sorting the images.
Advantages: Fast, low memory footprint, configurable framerate. Restrictions: The width and height
of each image must be divisible by 2.

### mpl_image_animator.py

Uses matplotlib to generate an animation from a sequence of images. Allows the user to specify an
output file and offers some support for sorting the images. Advantages: No height and width
restriction, offers a stack feature for vertically stacking multiple images using separate subplots.
Restrictions: Much slower than the ffmpeg version, and much larger memory footprint. *TODO*: Allow
user to specify framerate.

### front_tracker.py
Script for measuring the location of a flame or shock front over a sequence of snapshots. Allows the
user to specify the metrics (will only track 1 / 1000th the local enuc maximum by default), whether
to use a global or local maximum, domain bounds, and a few other things. For a usage description and
a full list of valid parameters, type `./front_tracker.py -h`. Should work for any dataset, but has
only been tested on flame wave ones. Restrictions: Currently only tracks along one dimension (the
user can tell it how to eliminate the others - either through slicing or averaging), and only tracks
percentages of the maximum of some field. Outputs to space-delimited data file called
front_tracking.dat by default.

### flame_speed.py
Script for reading in the front tracking dataset, plotting it, and fitting a line to some portion of
it. Usage: `./flame_speed.py data_file starting_index`, where starting index is the index of the
first datapoint to consider when fitting the line. The script prints out the slope of the line, the
r-squared value, and the fit error. The plot will appear squashed when it pops up - the window needs
to be enlarged. Uses `scipy` and `pandas`.

### multirays.py
Plot 1D vertical slices of axisymmetric datasets. It generates 3 ortho rays - one at each end of
the domain and one along the center. The variable to plot can be supplied as a command line argument
(e.g. `./multirays.py -v Temp`).

### parallel
Run multiple instances of a particular analysis script in parallel. Can only run one instance of
this at a time on a given directory.

### overview.py
Plots temperature, enuc, and z-velocity in a vertical stack. The fields to plot can be set by
modifying the fields list in script.

### time_series.py
Create a stacked plot of abar at a sequence of time points.

### schlieren.py
Make a Schlieren plot (a plot of ln{(∇^2 ρ) / ρ}) of a dataset.
# Problem Description

This is a simple 1-d flame problem.  It sets up a fuel and ash region
with a smooth transition in temperature and composition.  The density
is computed such that the domain is at constant pressure initially.
Thermal diffusion heats the fuel to the point of ignition and a flame
ignites and propagates to the right.



# Usage



# Publications
These files go together with the boosted flame for the `flame_wave`
problem.  You should compile with the `triple_alpha_plus_cago`
network.
*********
Flowchart
*********

Introduction
============

There are several different time-evolution methods currently
implemented in Castro. As best as possible, they share the same
driver routines and use preprocessor or runtime variables to separate
the different code paths.  These fall into two categories:

.. index:: castro.time_integration_method, USE_SIMPLIFIED_SDC, USE_TRUE_SDC

-  Strang+CTU: the Strang evolution does the burning on the
   state for :math:`\Delta t/2`, then updates the hydrodynamics using the
   burned state, and then does the final :math:`\Delta t/2` burning. No
   explicit coupling of the burning and hydro is done.  This code
   path uses the corner-transport upwind (CTU) method (the unsplit,
   characteristic tracing method of :cite:`colella:1990`).  This is the default method.

   The MHD solver uses this same driver.

-  SDC: a class of iterative methods that couples the advection and reactions
   such that each process explicitly sees the effect of the other.  We have
   two SDC implementations in Castro.

   - The "simplified SDC" method is based on the CTU hydro update.  We
     iterate over the construction of this term, using a lagged
     reaction source as inputs and do the final conservative update by
     integrating the reaction system using an ODE solver with the
     explicit advective source included in a
     piecewise-constant-in-time fastion.

   - The "true SDC" method.  This fully couples the hydro and reactions
     to either 2nd or 4th order.  This approximates the integral in
     time using a simple quadrature rule, and integrates the hydro
     explicitly and reactions implicitly to the next time node.
     Iterations allow each process to see one another and achieve
     high-order in time convergence.  This is described in :cite:`castro-sdc`.


The time-integration method used is controlled by
``castro.time_integration_method``.

  * ``time_integration_method = 0``: this is the original Castro method,
    described in :cite:`castro_I`.  This uses Strang splitting and the CTU
    hydrodynamics scheme.

  * ``time_integration_method = 1``: unused (in Castro 19.08 and
    earlier, this was a method-of-lines integration method with Strang
    splitting for reactions.)

  * ``time_integration_method = 2``: this is a full implementation of
    the spectral deferred corrections formalism, with both 2nd and 4th
    order integration implemented.  At the moment, this does not support
    multilevel domains.  Note: because of differences in the interfaces with the 
    default Strang method, you must compile with ``USE_TRUE_SDC = TRUE`` for this
    method to work (in particular, this defines ``EXTRA_THERMO`` which enables some
    additional EOS derivatives).

  * ``time_integration_method = 3``: this is the simplifed SDC method
    described above that uses the CTU hydro advection and an ODE
    reaction solve.  Note: because this requires a different set of
    state variables, you must compile with ``USE_SIMPLIFIED_SDC = TRUE`` for this
    method to work (in particular, this defines ``PRIM_SPECIES_HAVE_SOURCES``).

.. index:: USE_SIMPLIFIED_SDC, USE_TRUE_SDC

.. note::

   By default, the code is compiled for Strang-split CTU evolution
   (``time_integration_method = 0``).  Because the size of the
   different state arrays differs with the other integration schemes,
   support for them needs to be compiled in, using
   ``USE_SIMPLIFIED_SDC=TRUE`` for the simplified-SDC method
   (``time_integration_method=3``) and ``USE_TRUE_SDC=TRUE`` for the
   true SDC method (``time_integration_method = 2``).

.. note::

   MHD and radiation are currently only supported by the Strang+CTU
   evolution time integration method.

Several helper functions are used throughout:

.. index:: clean_state

-  ``clean_state``:
   There are many ways that the hydrodynamics state may become
   unphysical in the evolution. The ``clean_state()`` routine
   enforces some checks on the state. In particular, it

   #. enforces that the density is above ``castro.small_dens``

   #. enforces that the speeds in the state don't exceed ``castro.speed_limit``

   #. normalizes the species so that the mass fractions sum to 1

   #. syncs up the linear and hybrid momenta (for ``USE_HYBRID_MOMENTUM=TRUE``)

   #. resets the internal energy if necessary (too small or negative)
      and computes the temperature for all zones to be thermodynamically
      consistent with the state.

.. _flow:sec:nosdc:

Main Driver—All Time Integration Methods
========================================

This driver supports the Strang CTU integration.
(``castro.time_integration_method`` = 0)

The main evolution for a single step is contained in
``Castro_advance.cpp``, as ``Castro::advance()``. This does
the following advancement. Note, some parts of this are only done
depending on which preprocessor directives are defined at
compile-time—the relevant directive is noted in the [ ] at the start
of each step.

#. *Initialization* (``initialize_advance()``)

   This sets up the current level for advancement. The following
   actions are performend (note, we omit the actions taken for a retry,
   which we will describe later):

   -  Sync up the level information to the Fortran-side of Castro.

   -  Do any radiation initialization.

   -  Set the maximum density used for Poisson gravity tolerances.

   -  Initialize all of the intermediate storage arrays (like those
      that hold source terms, etc.).

   -  Swap the StateData from the new to old (e.g., ensures that
      the next evolution starts with the result from the previous step).

   -  Call ``clean_state``.

   -  Create the MultiFabs that hold the primitive variable information
      for the hydro solve.

   -  Zero out all of the fluxes.

   -  For true SDC, initialize the data at all time nodes (see :ref:`sec:flow_true_sdc`).

#. *Advancement*

   Call ``do_advance`` to take a single step, incorporating
   hydrodynamics, reactions, and source terms.

   For radiation-hydrodynamics, this step does the
   advective (hyperbolic) portion of the radiation update only.
   Source terms, including gravity, rotation, and diffusion are
   included in this step, and are time-centered to achieve second-order
   accuracy.

   .. index:: retry

   If ``castro.use_retry`` is set, then we subcycle the current
   step if we violated any stability criteria to reach the desired
   :math:`\Delta t`. The idea is the following: if the timestep that you
   took had a timestep that was not sufficient to enforce the stability
   criteria that you would like to achieve, such as the CFL criterion
   for hydrodynamics or the burning stability criterion for reactions,
   you can retry the timestep by setting ``castro.use_retry`` = 1 in
   your inputs file. This will save the current state data at the
   beginning of the level advance, and then if the criteria are not
   satisfied, will reject that advance and start over from the old
   data, with a series of subcycled timesteps that should be small
   enough to satisfy the criteria. Note that this will effectively
   double the memory footprint on each level if you choose to use it.
   See :ref:`ch:retry` for more details on the retry mechanism.

   .. note::

      Only Strang+CTU and simplified-SDC support retries.

#. [AUX_UPDATE] *Auxiliary quantitiy evolution*

   Auxiliary variables in Castro are those that obey a continuity
   equation (with optional sources) that are passed into the EOS, but
   not subjected to the constraint on mass fractions (summing to one).

   The advection and source terms are already dealt with in the
   main hydrodynamics advance (above step). A user-supplied routine
   ca_auxupdate can be provided here to further update these
   quantities.

#. *Radial data and [POINTMASS] point mass*

   If ``castro.spherical_star`` is set, then we average the state data
   over angles here to create a radial profile. This is then used in the
   boundary filling routines to properly set Dirichlet BCs when our domain
   is smaller than the star, so the profile on the boundaries will not
   be uniform.

   If ``castro.point_mass_fix_solution`` is set, then we
   change the mass of the point mass that optionally contributes to the
   gravitational potential by taking mass from the surrounding zones
   (keeping the density in those zones constant).

#. [RADIATION] *Radiation implicit update*

   The ``do_advance()`` routine only handled the hyperbolic
   portion of the radiation update. This step does the implicit solve
   (either gray or multigroup) to advance the radiation energies to the
   new time level. Note that at the moment, this is backward-difference
   implicit (first-order in time) for stability.

   This is handled by ``final_radiation_call()``.

#. [PARTICLES] *Particles*

   If we are including passively-advected particles, they are
   advanced in this step.

#. *Finalize*

   This cleans up at the end of a step:

   -  Update the flux registers to account for mismatches at
      coarse-fine interfaces. This cleans up the memory used during
      the step.

   -  Free any memory allocated for the level advance.


.. _sec:strangctu:

Strang+CTU Evolution
====================

``do_advance_ctu()`` in ``Castro_advance_ctu.cpp`` 

This described the flow using Strang splitting and the CTU
hydrodynamics (or MHD) method, including gravity, rotation, and
diffusion.  This integration is selected via
``castro.time_integration_method = 0``.

The system advancement: reactions, hydrodynamics, diffusion, rotation,
and gravity are all considered here.

Consider our system of equations as:

.. math:: \frac{\partial\Ub}{\partial t} = {\bf A}(\Ub) + \Rb(\Ub) + \Sb,

where :math:`{\bf A}(\Ub) = -\nabla \cdot \Fb(\Ub)`, with :math:`\Fb` the flux vector, :math:`\Rb` are the reaction
source terms, and :math:`\Sb` are the non-reaction source terms, which
includes any user-defined external sources, :math:`\Sb_{\rm ext}`. We use
Strang splitting to discretize the advection-reaction equations. In
summary, for each time step, we update the conservative variables,
:math:`\Ub`, by reacting for half a time step, advecting for a full time
step (ignoring the reaction terms), and reacting for half a time step.
The treatment of source terms complicates this a little. The actual
update, in sequence, looks like:

.. math::
   \begin{aligned}
   \Ub^\star &= \Ub^n + \frac{\dt}{2}\Rb(\Ub^n) \\
   \Ub^{n+1,(a)} &= \Ub^\star + \dt\, \Sb(\Ub^\star) \\
   \Ub^{n+1,(b)} &= \Ub^{n+1,(a)} + \dt\, {\bf A}(\Ub^\star) \\
   \Ub^{n+1,(c)} &= \Ub^{n+1,(b)} + \frac{\dt}{2}\, [\Sb(\Ub^{n+1,(b)}) - \Sb(\Ub^\star)] \\
   \Ub^{n+1}     &= \Ub^{n+1,(c)} + \frac{\dt}{2} \Rb(\Ub^{n+1,(c)})
   \end{aligned}
   :label: eq:source_correct

Note that in the first step, we add a full :math:`\Delta t` of the old-time
source to the state. This prediction ensures consistency when it
comes time to predicting the new-time source at the end of the update.
The construction of the advective terms, :math:`{\bf A(\Ub)}` is purely
explicit, and based on an unsplit second-order Godunov method. We
predict the standard primitive variables, as well as :math:`\rho e`, at
time-centered edges and use an approximate Riemann solver construct
fluxes.

At the beginning of the time step, we assume that :math:`\Ub` and the gravitational potential, :math:`\phi`, are
defined consistently, i.e., :math:`\rho^n` and :math:`\phi^n` satisfy the Poisson equation:

.. math::

   \Delta \phi^n = 4\pi G\rho^n

(see :ref:`ch:gravity` for more details about how the Poisson equation is solved.)  
Note that in
:eq:`eq:source_correct`, we can actually do some
sources implicitly by updating density first, and then momentum,
and then energy. This is done for rotating and gravity, and can
make the update more akin to:

.. math:: \Ub^{n+1,(c)} = \Ub^{n+1,(b)} + \frac{\dt}{2} [\Sb(\Ub^{n+1,(c)}) - \Sb(\Ub^n)]

If we are including radiation, then this part of the update algorithm
only deals with the advective / hyperbolic terms in the radiation update.

Here is the single-level algorithm. The goal here is to update the
``State_Type``  ``StateData`` from the old to new time (see
§ :ref:`soft:sec:statedata`). We will use the following notation
here, consistent with the names used in the code:

-  ``S_old`` is a MultiFab reference to the old-time-level
   ``State_Type`` data.

-  ``Sborder`` is a MultiFab that has ghost cells and is
   initialized from ``S_old``. This is what the hydrodynamic
   reconstruction will work from.

-  ``S_new`` is a MultiFab reference to the new-time-level
   ``State_Type`` data.

- ``old_source`` is a MultiFab reference to the old-time-level ``Source_Type`` data.

- ``new_source`` is a MultiFab reference to the new-time-level ``Source_Type`` data.


Single Step Flowchat
--------------------

In the code, the objective is to evolve the state from the old time,
``S_old``, to the new time, ``S_new``.

#. *Initialize*

   A. In ``initialize_do_advance()``, create ``Sborder``, initialized from ``S_old``

   B. Check for NaNs in the initial state, ``S_old``.


#. *React* :math:`\Delta t/2` [``strang_react_first_half()`` ]

   Update the solution due to the effect of reactions over half a time
   step. The integration method and system of equations used here is
   determined by a host of runtime parameters that are part of the
   Microphysics package. But the basic idea is to evolve the energy
   release from the reactions, the species mass fractions, and
   temperature through :math:`\Delta t/2`.

   Using the notation above, we begin with the time-level :math:`n` state,
   :math:`\Ub^n`, and produce a state that has evolved only due to reactions,
   :math:`\Ub^\star`.

   .. math::

      \begin{aligned}
          (\rho e)^\star &= (\rho e)^n + \frac{\dt}{2} \rho H_\mathrm{nuc} \\
          (\rho E)^\star &= (\rho E)^n + \frac{\dt}{2} \rho H_\mathrm{nuc} \\
          (\rho X_k)^\star &= (\rho X_k)^n + \frac{\dt}{2}(\rho\omegadot_k).
        \end{aligned}

   Here, :math:`H_\mathrm{nuc}` is the energy release (erg/g/s) over the
   burn, and :math:`\omegadot_k` is the creation rate for species :math:`k`.

   After exiting the burner, we call the EOS with :math:`\rho^\star`,
   :math:`e^\star`, and :math:`X_k^\star` to get the new temperature, :math:`T^\star`.

   Note that the density, :math:`\rho`, does not change via reactions in the
   Strang-split formulation.

   The reaction data needs to be valid in the ghost cells, so the reactions
   are applied to the entire patch, including ghost cells.

   After reactions, ``clean_state`` is called.

   At the end of this step, ``Sborder`` sees the effects of the
   reactions.

#. *Construct time-level n sources and apply*
   [``construct_old_gravity()``, ``do_old_sources()`` ]

   The time level :math:`n` sources are computed, and added to the
   StateData ``Source_Type``. 

   The sources that we deal with here are:

   A. sponge : the sponge is a damping term added to
      the momentum equation that is designed to drive the velocities to
      zero over some timescale. Our implementation of the sponge
      follows that of Maestro :cite:`maestro:III`

   B. external sources : users can define problem-specific sources
      in the ``problem_source.H`` file. Sources for the different
      equations in the conservative state vector, :math:`\Ub`, are indexed
      using the integer keys defined in ``state_indices.H``
      (e.g., URHO).

      This is most commonly used for external heat sources (see the
      ``toy_convect`` problem setup) for an example. But most
      problems will not use this.

   C. [``MHD``] thermal source: for the MHD system, we are including
      the "pdV" work for the internal energy equation as a source term
      rather than computing it from the Riemann problem.  This source is
      computed here for the internal energy equation.

   D. [``DIFFUSION``] diffusion : thermal diffusion can be
      added in an explicit formulation. Second-order accuracy is
      achieved by averaging the time-level :math:`n` and :math:`n+1` terms, using
      the same predictor-corrector strategy described here.

      Note: thermal diffusion is distinct from radiation hydrodynamics.

      Also note that incorporating diffusion brings in an additional
      timestep constraint, since the treatment is explicit. See
      Chapter :ref:`ch:diffusion` for more details.

   E. [``HYBRID_MOMENTUM``] angular momentum


   F. [``GRAVITY``] gravity:

      For full Poisson gravity, we solve for for gravity using:

      .. math::

         \gb^n = -\nabla\phi^n, \qquad
               \Delta\phi^n = 4\pi G\rho^n,

      The construction of the form of the gravity source for the
      momentum and energy equation is dependent on the parameter
      ``castro.grav_source_type``. Full details of the gravity
      solver are given in Chapter :ref:`ch:gravity`.


   G. [``ROTATION``] rotation

      We compute the rotational potential (for use in the energy update)
      and the rotational acceleration (for use in the momentum
      equation). This includes the Coriolis and centrifugal terms in a
      constant-angular-velocity co-rotating frame. The form of the
      rotational source that is constructed then depends on the
      parameter ``castro.rot_source_type``. More details are
      given in Chapter :ref:`ch:rotation`.

   The source terms here are evaluated using the post-burn state,
   :math:`\Ub^\star` (``Sborder``), and later corrected by using the
   new state just before the burn, :math:`\Ub^{n+1,(b)}`. This is compatible
   with Strang-splitting, since the hydro and sources takes place
   completely inside of the surrounding burn operations.

   The old-time source terms are stored in ``old_source``.

   The sources are then applied to the state after the burn,
   :math:`\Ub^\star` with a full :math:`\Delta t` weighting (this will
   be corrected later). This produces the intermediate state,
   :math:`\Ub^{n+1,(a)}` (stored in ``S_new``).

#. *Construct the hydro / MHD update* [``construct_ctu_hydro_source()``, ``construct_ctu_mhd_source()``]

   The goal is to advance our system considering only the advective
   terms (which in Cartesian coordinates can be written as the
   divergence of a flux).

   We do the hydro update in two parts—first we construct the
   advective update and store it in the hydro_source
   MultiFab, then we do the conservative update in a separate step. This
   separation allows us to use the advective update separately in more
   complex time-integration schemes.

   In the Strang-split formulation, we start the reconstruction using
   the state after burning, :math:`\Ub^\star` (``Sborder``).  For the
   CTU method, we predict to the half-time (:math:`n+1/2`) to get a
   second-order accurate method. Note: ``Sborder`` does not know of
   any sources except for reactions. 

   The method done here differs depending on whether we are doing hydro or MHD.

   A. hydrodynamics

      The advection step is complicated, and more detail is given in
      Section :ref:`Sec:Advection Step`. Here is the summarized version:

      i. Compute primitive variables.

      ii. Convert the source terms to those acting on primitive variables

      iii. Predict primitive variables to time-centered edges.

      iv. Solve the Riemann problem.

      v. Compute fluxes and advective term.

   B. MHD

      The MHD update is described in :ref:`ch:mhd`.

   To start the hydrodynamics/MHD source construction, we need to know
   the hydrodynamics source terms at time-level :math:`n`, since this
   enters into the prediction to the interface states. This is
   essentially the same vector that was computed in the previous step,
   with a few modifications. The most important is that if we set
   ``castro.source_term_predictor``, then we extrapolate the source
   terms from :math:`n` to :math:`n+1/2`, using the change from the
   previous step.

   Note: we neglect the reaction source terms, since those are already
   accounted for in the state directly, due to the Strang-splitting
   nature of this method.

   The update computed here is then immediately applied to
   ``S_new``.

#. *Clean State* [``clean_state()``]

   This is done on ``S_new``.

   After these checks, we check the state for NaNs.

#. *Update radial data and center of mass for monopole gravity*

   These quantities are computed using ``S_new``.

#. *Correct the source terms with the n+1
   contribution* [``construct_new_gravity()``, ``do_new_sources`` ]

   If we are doing self-gravity, then we first compute the updated gravitational
   potential using the updated density from ``S_new``.

   Now we correct the source terms applied to ``S_new`` so they are time-centered.
   Previously we added :math:`\Delta t\, \Sb(\Ub^\star)` to the state, when
   we really want 
   :math:`(\Delta t/2)[\Sb(\Ub^\star + \Sb(\Ub^{n+1,(b)})]` .

   We start by computing the source term vector :math:`\Sb(\Ub^{n+1,(b)})`
   using the updated state, :math:`\Ub^{n+1,(b)}`. We then compute the
   correction, :math:`(\Delta t/2)[\Sb(\Ub^{n+1,(b)}) - \Sb(\Ub^\star)]` to
   add to :math:`\Ub^{n+1,(b)}` to give us the properly time-centered source,
   and the fully updated state, :math:`\Ub^{n+1,(c)}`. 

   This correction is stored
   in the ``new_sources`` MultiFab [1]_.

   In the process of updating the sources, we update the temperature to
   make it consistent with the new state.

#. *React* :math:`\Delta t/2` [``strang_react_second_half()``]

   We do the final :math:`\dt/2` reacting on the state, begining with :math:`\Ub^{n+1,(c)}` to
   give us the final state on this level, :math:`\Ub^{n+1}`.

   This is largely the same as ``strang_react_first_half()``, but
   it does not currently fill the reactions in the ghost cells.

#. *Finalize* [``finalize_do_advance()``]

   Finalize does the following:

   A. for the momentum sources, we compute :math:`d\Sb/dt`, to use in the
      source term prediction/extrapolation for the hydrodynamic
      interface states during the next step.

   B. If we are doing the hybrid momentum algorithm, then we sync up
      the hybrid and linear momenta

A summary of which state is the input and which is updated for each of
these processes is presented below:

.. table:: update sequence of state arrays for Strang-CTU
   :align: center

   +--------------------+-----------+---------------------+---------------------+
   | *step*             | ``S_old`` | ``Sborder``         | ``S_new``           |
   +====================+===========+=====================+=====================+
   | 1. init            | input     | updated             |                     |
   +--------------------+-----------+---------------------+---------------------+
   | 2. react           |           | input / updated     |                     |
   +--------------------+-----------+---------------------+---------------------+
   | 3. old sources     |           | input               | updated             |
   +--------------------+-----------+---------------------+---------------------+
   | 4. hydro           |           | input               | updated             |
   +--------------------+-----------+---------------------+---------------------+
   | 5. clean           |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+
   | 6. radial / center |           |                     | input               |
   +--------------------+-----------+---------------------+---------------------+
   | 7. correct sources |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+
   | 8. react           |           |                     | input / updated     |
   +--------------------+-----------+---------------------+---------------------+


.. _sec:flow_true_sdc:

SDC Evolution
=============

The SDC evolution is selected by ``castro.time_integration_method = 2``.  It
does away with Strang splitting and instead couples the reactions and hydro
together directly.

.. note::

   At the moment, the SDC solvers do not support multilevel or AMR
   simulation.

.. note::

   The code must be compiled with ``USE_TRUE_SDC = TRUE`` to use this
   evolution type.

The SDC solver follows the algorithm detailed in :cite:`castro-sdc`.
We write our evolution equation as:

.. math::
   \frac{\partial \Ub}{\partial t} = {\bf A}(\Ub) + {\bf R}(\Ub)

where :math:`{\bf A}(\Ub) = -\nabla \cdot {\bf F}(\Ub) + {\bf S}(\Ub)`, with the 
hydrodynamic source terms, :math:`{\bf S}` grouped together with the flux divergence.

The SDC update looks at the solution a several time nodes (the number
depending on the desired temporal order of accuracy), and iteratively
updates the solution from node :math:`m` to :math:`m+1` as:

.. math::
   \begin{align}
   \avg{\Ub}^{m+1,(k+1)} = \avg{\Ub}^{m,(k+1)} &+ \Delta t \left [ \avg{{\bf A}(\Ub)}^{m,(k+1)} - \avg{{\bf A}(\Ub)}^{m,(k)} \right ] \\
                                   &+ \Delta t \left [ \avg{{\bf R}(\Ub)}^{m+1,(k+1)} - \avg{{\bf R}(\Ub)}^{m+1,(k)} \right ] \\
                                   &+ \int_{t^m}^{t^{m+1}} \left [ \avg{{\bf A}(\Ub)}^{(k)} + \avg{{\bf R}(\Ub)}^{(k)} \right ] dt
   \end{align}


.. index:: castro.sdc_order, castro.sdc_quadrature

Where :math:`k` is the iteration index.  In the SDC formalism, each
iteration gains us an order of accuracy in time, up to the order with
which we discretize the integral at the end of the above expression.
We also write the conservative state as :math:`\avg{\Ub}` to remind us
that it is the cell average and not the cell-center.  This distinction
is important when we consider the 4th order method.

In Castro, there are two parameters that together determine the number
and location of the temporal nodes, the accuracy of the integral, and
hence the overall accuracy in time: ``castro.sdc_order`` and
``castro.sdc_quadrature``. 

``castro.sdc_quadrature = 0`` uses
Gauss-Lobatto integration, which includes both the starting and ending
time in the time nodes.  This gives us the trapezoid rule for 2nd
order methods and Simpson's rule for 4th order methods.  Choosing
``castro.sdc_quadrature = 1`` uses Radau IIA integration, which includes
the ending time but not the starting time in the quadrature.


.. table:: SDC quadrature summary
   :align: center

   +--------------+---------------+---------------+-------------------+------------------+
   |``sdc_order`` |``quadrature`` |  # of         |  temporal         |  description     |
   |              |               |  time nodes   |  accuracy         |                  |
   +==============+===============+===============+===================+==================+
   |       2      |         0     |          2    |                2  | trapezoid rule   |
   +--------------+---------------+---------------+-------------------+------------------+
   |       2      |         1     |          3    |                2  | Simpson's rule   |
   +--------------+---------------+---------------+-------------------+------------------+
   |       4      |         0     |          3    |                4  | Radau 2nd order  |
   +--------------+---------------+---------------+-------------------+------------------+
   |       4      |         1     |          4    |                4  | Radau 4th order  |
   +--------------+---------------+---------------+-------------------+------------------+

The overall evolution appears as:

.. index:: k_new, A_old, A_new, R_old

#. *Initialization* (``initialize_advance``)

   Here we create the ``MultiFab`` s that store the needed information
   at the different time nodes.  Each of the quantities below is a
   vector of size ``SDC_NODES``, whose components are the ``MultiFab``
   for that time node:


    * ``k_new`` : the current solution at this time node.

      Note that
      ``k_new[0]`` is aliased to ``S_old``, the solution at the start
      of the step, since this never changes (so long as the 0th time
      node is the start of the timestep).

    * ``A_old`` : the advective term at each time node at the old
      iteration.

    * ``A_new`` : the advective term at each time node at the current
      iteration.
    
    * ``R_old`` : the reactive source term at each time node at the old
      iteration.

#. *Advancement*

   Our iteration loop calls ``do_advance_sdc`` to update the solution through
   all the time nodes for a single iteration.

   The total number of iterations is ``castro.sdc_order`` + ``castro.sdc_extra``.

#. *Finalize*

   This clears the ``MultiFab`` s we allocated.

SDC Single Iteration Flowchart
------------------------------

.. index:: do_advance_sdc

The update through all time nodes for a single iteration is done by
``do_advance_sdc``.  The basic update appears as:

Throughout this driver we use the ``State_Type`` ``StateData`` as
storage for the current node.  In particular, we use the new time slot
in the ``StateData`` (which we refer to as ``S_new``) to allow us to
do ``FillPatch`` operations.

#. *Initialize*

   We allocate ``Sborder``.  Just like with the Strang CTU driver, we
   will use this as input into the hydrodynamics routines.

#. Loop over time nodes

   We'll use ``m`` to denote the current time node and ``sdc_iter`` to
   denote the current (0-based) iteration.  In our loop over time
   nodes, we do the following for each node:

   * Load in the starting data

     * ``S_new`` :math:`\leftarrow` ``k_new[m]``

     * ``clean_state`` on ``S_new``

     * Fill ``Sborder`` using ``S_new``

   * Construct the hydro sources and advective term

     Note: we only do this on the first time node for ``sdc_iter`` = 0, and
     we don't need to do this for the last time node on the last
     iteration.

     * Call ``do_old_sources`` filling the ``Source_Type``
       ``StateData``, ``old_source``.

     * Convert the sources to 4th order averages if needed.

     * Convert the conserved variables to primitive variables

     * Call ``construct_mol_hydro_source`` to get the advective update
       at the current time node, stored in ``A_new[m]``.
 
   * Bootstrap the first iteration.

     For the first iteration, we don't have the old iteration's
     advective and reaction terms needed in the SDC update.  So for
     the first time node (``m = 0``) on the first iteration, we do:

     * ``A_old[n]`` = ``A_old[0]``, where ``n`` loops over all time nodes.

     * Compute the reactive source using the ``m = 0`` node's state and
       store this in ``R_old[0]``.

       Then fill all other time nodes as: ``R_old[n]`` = ``R_old[0]``

    * Do the SDC update from node ``m`` to ``m+1``.

      We call ``do_sdc_update()`` to do the update in time to the next
      node.  This solves the nonlinear system (when we have reactions)
      and stores the solution in ``k_new[m+1]``.

#. Store the advective terms for the next iteration.

   Since we are done with this iteration, we do: ``A_old[n]``
   :math:`\leftarrow` ``A_new[n]``.

   We also store ``R_old`` for the next iteration.  We do this by
   calling the reaction source one last time using the data for each
   time node.

#. Store the new-time solution.

   On the last iteration, we save the solution to the ``State_Type`` ``StateData``:

   ``S_new`` :math:`\leftarrow` ``k_new[SDC_NODES-1]``

#. Call ``finalize_do_advance`` to clean up the memory.
   

Simplified-SDC Evolution
========================

The simplified SDC method uses the CTU advection solver together with
an ODE solution to update the compute advective-reacting system.  This
is selected by ``castro.time_integration_method = 3``.

We use one additional StateData type here, ``Simplified_SDC_React_Type``,
which will hold the reactive source needed by hydrodynamics.

.. note::

   The code must be compiled with ``USE_SIMPLIFIED_SDC = TRUE`` to use this
   evolution type.


We express our system as:

.. math:: \Ub_t = \mathcal{A}(\Ub) + \Rb(\Ub)

here :math:`\mathcal{A}` is the advective source, which includes both the
flux divergence and the hydrodynamic source terms (e.g. gravity):

.. math:: \mathcal{A}(\Ub) = -\nabla \cdot \Fb(\Ub) + \Sb

The simplified-SDC version of the main advance loop looks similar to the Strang CTU
version, but includes an iteration loop over the hydro, gravity, and
reaction update. So the only difference happens in step 2 of the
flowchart outlined in § \ `2 <#flow:sec:nosdc>`__. In particular this
step now proceeds as a loop over ``do_advance_ctu``.  The differences
with the Strang CTU version are highlighted below.


Note that the
radiation implicit update is not done as part of the Simplified-SDC iterations.

Simplified_SDC Hydro Advance
----------------------------

The evolution in ``do_advance`` is substantially different than the
Strang case. In particular, reactions are not evolved. Here we
summarize those differences.

#. *Initialize* [``initialize_do_advance()``]

   This is unchanged from the initialization in the CTU Strang algorithm.

#. *Construct time-level n sources and apply*
   [``construct_old_gravity()``, ``do_old_sources()``]

   Unlike the Strang case, there is no need to extrapolate source
   terms to the half-time for the prediction (the
   ``castro.source_term_predictor`` parameter), since the
   Simplified-SDC provides a natural way to approximate the
   time-centered source—we simply use the iteratively-lagged new-time
   source.  We add the corrector from the previous iteration to the
   source Multifabs before adding the current source.  The corrector
   (stored in ``source_corrector``) has the form:

   .. math::

      \Sb^\mathrm{corr} = \frac{1}{2} \left ( \Sb^{n+1,(k-1)} - S^n \right )

   where :math:`\Sb^n` does not have an iteration subscript, since we always have the
   same old time state.  

   Applying this corrector to the the source at time :math:`n`, will give
   us a source that is time-centered,

   .. math::

      {\bf S}(\Ub)^{n+1/2} = \frac{1}{2} \left ( {\bf S}(\Ub)^n + {\bf S}(\Ub)^{n+1,(k-1)} \right )

   For constructing the time-level :math:`n` source, there are no
   differences compared to the Strang algorithm.

#. *Construct the hydro update* [``construct_hydro_source()``]

   In predicting the interface states, we use an iteratively-lagged
   approximation to the reaction source on the primitive variables,
   :math:`\mathcal{I}_q^{k-1}`.  This addition is done in
   ``construct_ctu_hydro_source()`` after the source terms are
   converted to primitive variables.

   The result of this is an approximation to :math:`- [\nabla \cdot {\bf F}]^{n+1/2}` (not yet the full :math:`\mathcal{A}(\Ub)`)
   stored in ``hydro_sources``.

#. *Clean State* [``clean_state()``]

#. *Update radial data and center of mass for monopole gravity*

#. *Correct the source terms with the n+1 contribution*
   [``construct_new_gravity()``, ``do_new_sources()`` ]

#. *React* :math:`\Delta t` [``react_state()``]

   We first compute :math:`\mathcal{A}(\Ub)` using ``hydro_sources``,
   ``old_source``, and ``new_source`` via the ``sum_of_source()``
   function.  This produces an advective source of the form:
   
   .. math::

      \left [ \mathcal{A}(\Ub) \right ]^{n+1/2} = - [\nabla \cdot {\bf F}]^{n+1/2} + \frac{1}{2} (S^n + S^{n+1})

   We burn for the full :math:`\Delta t` including the advective
   update as a source, integrating

      .. math:: \frac{d\Ub}{dt} = \left [ \mathcal{A}(\Ub) \right ]^{n+1/2} + \Rb(\Ub)

   The result of evolving this equation is stored in ``S_new``.

   Note, if we do not actually burn in a zone (because we don't meet
   the thermodynamic threshold) then this step does nothing, and the
   state updated just via hydrodynamics in ``S_new`` is kept.

#. *Clean state*: This ensures that the thermodynamic state is
   valid and consistent.

#. *Construct reaction source terms*: Construct the change
   in the primitive variables due only to reactions over the
   timestep, :math:`\mathcal{I}_q^{k}`. This will be used in the next
   iteration.

#. *Finalize* [``finalize_do_advance()``]

   This differs from Strang finalization in that we do not construct
   :math:`d\Sb/dt`, but instead store the total hydrodynamical source
   term at the new time. As discussed above, this will be used in the
   next iteration to approximate the time-centered source term.

.. [1]
   The correction for gravity is slightly different since we directly compute the time-centered gravitational source term using the hydrodynamic fluxes.
**************************
Development Best Practices
**************************

Coding Conventions
==================

Castro development should do its best to adhere to the following coding
style guidelines.

General
-------

* Indentation should be spaces, not tabs, with 4 spaces preferred.


C++
---

* Conditional should appear as:

  .. code:: c++

     if (condition)
     {
         ...
     }
     else if (condition)
     {
         ...
     }
     else
     {
         ...
     }


C++ to Fortran
--------------

* All C to Fortran interfaces should use the ISO C binding.  The
  Fortran subroutine should use

  .. code:: fortran

     subroutine subroutine_name(args) bind(C, name="subroutine_name")

* Data passed by reference from C++ should use ``*`` and not ``&``.

* Scalars should be passed by value, using the Fortran ``value`` attribute.

* Fortran routines that are called from C++ should being with ``ca_``.
  This is true even if these routines are also called by other
  Fortran.

* A separate tile ``lo`` and ``hi`` should be passed in to each
  subroutine, with the expectation that the Fortran subroutine will
  act exactly over the domain (lo, hi). If you want to include ghost
  zones in the calculation, include them in your box using
  ``mfi.growntilebox(ng)`` instead of ``mfi.validbox()``.


Fortran
-------

* Put all routines in a module—this ensures that argument lists are
  checked.

* Use the ``only`` clause in module ``use`` statements to explicitly
  make clear what is being accessed.

* In a module, there should be no "top-level" ``use`` statements (with
  the exception of getting access to the ``rt`` type).  Instead each
  function / subroutine in the module  should use what it needs directly.

* New Fortran files should have the .F90 file extension, not the .f90
  file extension, so that they can be preprocessed.


Documentation
-------------

C++ routines should use Doxygen style comments.

* C++ functions should be documented in the header file using the style:

  .. code:: c++

     ///
     /// Description of the function
     ///
     /// @param bar       Brief description of the variable
     ///
     void foo(int bar) { ...

* Member variables can either be documented using the above style of comment block or
  with a brief inline description:

  .. code:: c++

     int var; ///< Brief description after the variable

Fortran functions should use Sphinx style comments

* Fortran functions should be documented by placing a comment block
  immediately after their prototype (i.e. `without` a line in betwen ) using the style:

  .. code:: fortran

     subroutine foo(bar)
       ! Description of the function

       use some_module

       implicit none

       integer, intent(inout) :: bar   ! Brief description of bar
       ...

  Documentation for modules should be similarly formatted, with the comment block again
  coming `immediately` after the module definition.

Castro Releases
===============

This outlines the procedure for doing the monthly Castro release.

Castro uses submodules for dependencies, this means that, at a
minimum, we must update the AMReX and  Microphysics submodules monthly when we
issue new releases. The releases for AMReX and Microphysics must be done
first. Then navigate to each submodule directory, checkout the new
tag, and then from the top-level directory of Castro do a "git add" on
the ``external/`` directory to store the new tags. So, for example, at
the beginning of March 2020 we would first issue the ``20.03`` tag on
Microphysics, and wait for AMReX to release a ``20.03`` tag, then do::

   cd $CASTRO_HOME/external
   cd amrex
   git pull
   git checkout 20.03
   cd ..
   cd Microphysics
   git pull
   git checkout 20.03
   cd ..
   git add -u .
   git commit -m "Update AMReX and Microphysics to release 20.03"

Then we can proceed with issuing our own release.


Each month the ``development`` branch is merged into ``main`` and a
release is tagged.  This needs to be done in coordination with its
submodule dependencies.  The general procedure is:

  * Do 'git pull' in both main and development branches.  (Use `git
    checkout xxx` to change to branch xxx.

  * In main branch, do `git merge development`.  Fix any conflicts
    if there are any.  (There should not be any conflicts unless a
    commit is checked into main directly without going through
    development.)

  * In main branch, commit new release notes (``CHANGES.md``)
    summarizing changes since last major release.

  * Tag the new release: ``git tag -m "Castro YY.MM" YY.MM``

  * ``git push``

  * ``git push --tags``

  * ``git checkout development``

  * ``git merge main``

  * ``git push``


Interim updates
---------------

When breaking changes to Microphysics occur in its development branch
that Castro depends on, we must update the Microphysics submodule on
the Castro development branch in the same way, replacing the git
checkout statement with the latest commit hash on the Microphysics
development branch. (A git submodule always tracks a specific
commit/tag on the target repo -- it is not configured to automatically
track a particular branch.)  Since such breaking changes usually are
accompanied by a Castro change, it is best practice to ensure that
the PRs in both Microphysics and Castro have been approved, then
merge the Microphysics PR, then add the update to the Microphysics
submodule to the Castro PR, then merge. A similar process applies for AMReX.


Continuous Integration
======================

We github actions to run integration tests on the code and to build and deploy the documentation.

Currently, we run the `clang static analyzer <https://clang-analyzer.llvm.org/>`_, which finds potential bugs in the code. It also runs a script to convert any tabs in the code into spaces. Both of these are run on pull requests to the Castro GitHub repo, and are run weekly on the development branch. 

****************
Tracer particles
****************

Tracer particles are to track the Lagrangian evolution of a model
fluid using discrete particles. In hydrodynamical simulations based on
an Eulerian grid (including CASTRO), thermodynamic variables at a
given time are derived by solving the equations of motion of a fluid
between cells. Therefore, in this scheme, the physical quantities that
we can access to are not discretized quantities at any given position,
but rather average values over each cell. However, employing discrete
particles, passively advected with the fluid flow, allows us to obtain
local instantaneous thermodynamic variables, such as the temperature
and the density, at well-defined positions, independent of the spatial
resolution, i.e., the spatial cell size. This means that we can follow
the evolution of the fluid at any given position and time.

CASTRO provides a tracer particle scheme with useful options. In this
scheme, particles are advanced using the midpoint method either with
the cell-centered velocities or the face-centered velocities
(Marker-And-Cell method) [1]_. The number and the initial positions of
particles are flexibly determined according to the purpose of a given
model.

Initializing the Particles
==========================

One must include the tracer particles in the ``GNUmakefile`` by setting::

   USE_PARTICLES = TRUE


And the particles can be initialized via::

   castro.do_tracer_particles = 1

in the ``inputs`` file.

If one wants to investigate the evolution of fluid motions starting from specific positions (or a certain range of area or volume), one should manually specify the positions of particles by providing an input file containing the total number and the initial positions of the particles.
The input file should be in the same directory where your inputs file is located. The name of the input file is determined via::

   particles.particle_init_file = particle_file

Here *particle_file* is the user-specified name of the file. The first
line in this file is assumed to contain the number of particles. Each
line after that contains the positions in a coordinate system adopted
for your model. For 3-D cartesian coordinates, :math:`x ~y ~z` For
example, an input file for a model fluid with 6 particles in 2-D
Cartesian coordinates may look like::

    6
    3.28125e+08 9.9198e+08 
    5.46875e+08 9.9198e+08 
    7.65625e+08 9.9198e+08 
    9.84375e+08 9.9198e+08 
    1.20312e+09 9.9198e+08 
    1.42188e+09 9.9198e+08 

According to this input file, the 6 particles will be positioned at
the same height (same :math:`y` coordinate in the second column),
equally spaced in :math:`x` direction (the first column except for the
particle number on the first line) from :math:`3.28\times10^{8} {\rm
~cm}` to :math:`1.42\times 10^{9} {\rm ~cm}`.

.. _particles:output_file:

Output file
===========

The output files are stored in a directory whose name is determined by
a variable ``particles.timestamp_dir``. For example, if the variable is
set as follows::

  particles.timestamp_dir = particle\_dir

A directory *particle_dir* is automatically made with the directories
for the main CASTRO output file (``pltXXXXX``) once a simulation starts
and the particle output files are stored inside that directory.

The name of the output file consists of ``Timestamp_`` along with a
number at the end. The number increases (typically from 00) as more
processors are involved in following the trajectories of particles. In
parallel computing, a computational domain is divided according to the
number of processors requested. Then each processor only follows the
particles on the domain assigned to that processor and records their
positions and velocities at any given time in a different output
file. Since it is possible for particles to move from one domain to
another during the evolution, its history can be stored in different
files. More output files (with larger numbers at the end of the file
name) can be produced as more processors track the particles.

By default, the output file contains the positions and velocities of
all particles at a given time, meaning [:math:`3+ 2\times`\
dimensionality] columns. For example, for particles in a 3-D domain,
the columns in the output file are,

:math:`{\rm index1}~~{\rm index2}~~x~~ y~~ z~~ t~~ v_{\rm x} ~~v_{\rm y}~~ v_{\rm z}~~ [\rho ~~ T]`

The first two integers correspond to the particle index and the
processor number.  One should use the two numbers in order to identify
a particle and extract its history (i.e., the trajectory in :numref:`fig:particletrajectory`.

.. figure:: fluid_motion.png

   A model atmosphere with the arrows showing the direction of the fluid motion.

.. _fig:particletrajectory:
.. figure:: tracer_trajectory.png 

   The trajectories of 500 particles following the fluid motion on the
   atmosphere. The particles are initially positioned at five
   different heights, :math:`y=13000\mathrm{~km},~11000\mathrm{~km},~
   8000\mathrm{~km},~ 6000\mathrm{~km}, ~38000\mathrm{~km}` (100
   particles at each height).  The solid lines represent the
   trajectories of the particles.

One can also add the last two columns :math:`[\rho ~~ T]`, i.e., the
local density and local temperature of fluid at the position of each
particle by setting the following::

    particles.timestamp_temperature= 1
    particles.timestamp_density = 1

For example, let’s consider 10 particles on a domain. If 4 out 10
particles are initially on a processor and the rest are on another
processor, this means two processors are tracking the particles and
two output files are produced. In the output file written by the
processor with 4 particles, one can find that four lines are stored at
the same time and each line corresponds to each particle info. while
in the other output file for the other 6 particles, 6 lines are stored
at the same time.

If ``particles.write_in_plotfile`` = 1, the particle data are stored
in a binary file along with the main CASTRO output plotfile in
directories ``pltXXXXX/Tracer/``.

Run-time Screen Output
----------------------

The verbosity written to the screen at run-time is turned off by setting::

    particles.v = 0


.. [1]
   One can simplify interpolation with the cell-centered
   velocity. However, this can lead to decoupling of the pressure and
   the velocity components, possibly resulting in instability. This
   can be avoided with the face-centered velocity
**********
Regridding
**********

The details of the regridding strategy are described in
§ :ref:`sec:tagging`; here we cover how the input parameters can
control the gridding.

As described later, the user defines Fortran subroutines which tag
individual cells at a given level if they need refinement. This list
of tagged cells is sent to a grid generation routine, which uses the
Berger-Rigoutsos algorithm :cite:`br-refine` to create rectangular
grids that contain the tagged cells.

The relevant runtime parameters are:

  * ``amr.regrid_file``: name of file from which to read the grids
    (text; default: no file)

    If set to a filename, e.g. ``fixed_girds``, then list of grids at
    each fine level are read in from this file during the gridding
    procedure. These grids must not violate the ``amr.max_grid_size``
    criterion. The rest of the gridding procedure described below will
    not occur if ``amr.regrid_file`` is set.

  * ``amr.n_error_buf``: radius of additional tagging
    around already tagged cells (integer :math:`\geq 0`; default: 1)

  * ``amr.max_grid_size``: maximum size of a grid in any
    direction (integer :math:`> 0`; default: 128 (2-d), 32 (3-d))

    Note: ``amr.max_grid_size`` must be even, and a multiple of
    ``amr.blocking_factor`` at every level.

  * ``amr.blocking_factor``: grid size must be a multiple of this
    (integer :math:`> 0`; default: 2)
    ``amr.blocking_factor`` at every level must be a power of 2
    and the domain size must be a multiple of ``amr.blocking_factor``
    at level 0.

    .. note:: This can be very important for elliptic problems with
       multigrid. A higher blocking factor allows the multigrid
       algorithm to coarsen more at the lowest level, reducing the
       amount of work required by the bottom solver.

  * ``amr.grid_eff``: grid efficiency (Real :math:`>0` and :math:`<1`;
    default: 0.7)

    When creating a refined grid, do we make boxes that only include
    the coarse cells that were explicitly tagged for refinement? or do
    we allow ourselves to encompass nearby, untagged cells in order to
    make larger and more regular boxes? This is the grid efficiency.

    When ``blocking_factor = 1``, *grid efficiency* is exactly the
    fraction of refined cells in the fine ``BoxArray`` which
    correspond to coarse cells which were tagged. For other blocking
    factors, we actually apply ``grid_eff`` at the level which has been
    coarsened by ``blocking_factor``, so it is no longer strictly this
    fraction, but the idea is still the same.

  * ``amr.refine_grid_layout``: refine grids more if # of
    processors :math:`>` # of grids (0 if false, 1 if true; default: 1)

Note also that ``amr.n_error_buf``, ``amr.max_grid_size`` and
``amr.blocking_factor`` can be read in as a single value which is
assigned to every level, or as multiple values, one for each level.

As an example, consider::

    amr.grid_eff = 0.9
    amr.max_grid_size = 64
    amr.blocking_factor} = 32

The grid efficiency, ``amr.grid_eff``, means that during the grid
creation process, at least 90% of the cells in each grid at the level
at which the grid creation occurs must be tagged cells. A higher
grid efficiency means fewer cells at higher levels, but may result
in the production of lots of small grids, which have inefficient cache
and OpenMP performance and higher communication costs.

The ``amr.max_grid_size`` parameter means that the final grids will be
no longer than 64 cells on a side at every level.  Alternately, we
could specify a value for each level of refinement as
``amr.max_grid_size = 64 32 16`` in which case our final grids will be
no longer than 64 cells on a side at level 0, 32 cells on a side at
level 1, and 16 cells on a side at level 2. The
``amr.blocking_factor`` means that all of the final grids will be
multiples of 32 at all levels.  Again, this can be specified on a
level-by-level basis, like ``amr.blocking_factor = 32 16 8``, in which
case the dimensions of all the final grids will be multiples of 32 at
level 0, multiples of 16 at level 1, and multiples of 8 at level 2.

Getting good performance
~~~~~~~~~~~~~~~~~~~~~~~~

These parameters can have a large impact on the performance
of Castro, so taking the time to experiment with is worth the effort.
Having grids that are large enough to coarsen multiple levels in a
V-cycle is essential for good multigrid performance in simulations
that use self-gravity.



How grids are created
~~~~~~~~~~~~~~~~~~~~~

The gridding algorithm proceeds in this order:

#. Grids are created using the Berger-Rigoutsos clustering algorithm
   modified to ensure that all new fine grids are divisible by
   ``amr.blocking_factor``.

#. Next, the grid list is chopped up if any grids are larger than
   ``max_grid_size``.  Note that because ``amr.max_grid_size`` is a
   multiple of ``amr.blocking_factor`` the ``amr.blocking_factor``
   criterion is still satisfied.

#. Next, if ``amr.refine_grid_layout = 1`` and there are more
   processors than grids, and if ``amr.max_grid_size`` / 2 is a
   multiple of ``amr.blocking_factor``, then the grids will be
   redefined, at each level independently, so that the maximum length
   of a grid at level :math:`\ell`, in any dimension, is
   ``amr.max_grid_size`` [:math:`\ell`] / 2.

#. Finally, if ``amr.refine_grid_layout = 1``, and there are still
   more processors than grids, and if ``amr.max_grid_size`` / 4 is a
   multiple of ``amr.blocking_factor``, then the grids will be
   redefined, at each level independently, so that the maximum length
   of a grid at level :math:`\ell`, in any dimension, is
   ``amr.max_grid_size`` [:math:`\ell`] / 4.
**********
References
**********

.. bibliography:: refs.bib
   :style: plain

.. _ch:sdc:

*****************************
Spectral Deferred Corrections
*****************************

The Castro SDC solver couples the hydrodynamics tightly together,
iteratively improving the convergence of the solution.  This is the
basis of the 4th order accurate Castro solver.  The algorithm is described
in :cite:`castro-sdc`.

.. note::

   Here we are referring to the full SDC time integration scheme
   (``castro.time_integration_method = 2``), not the simplified-SDC solver.


The options that describe the quadrature and iterations are:

* ``castro.sdc_order`` : the desired spatial and temporal order.  2 and 4 are supported.

* ``castro.sdc_quadrature`` : the quadrature scheme used for the
  time-integration.  This determines the number and location of the
  temporal nodes.  Supported values are 0 for Gauss-Lobatto and 1 for
  Radau IIA.

* ``castro.sdc_extra`` : the number of extra iterations to take.  By
  default the number of iterations used is equal to the value of
  ``sdc_order``.


The options that affect the nonlinear solve are:

* ``sdc_solver`` : the method we use to do the nonlinear solution of
  the reaction system.  Values are:

  * 1 : pure Newton iteration (we subdivide the time interval if
    needed to get the Newton method to converge).

  * 2 : use VODE to solve the nonlinear system by expressing it as an ODE system.

  * 3 : use VODE for the first iteration and then Newton for the
    subsequent iterations.

* ``sdc_solver_tol_dens`` : the relative error on the density in solving the nonlinear system.

* ``sdc_solver_tol_spec`` : the relative error on the partial densities, :math:`(\rho X_k)`
  in the nonlinear solve.

* ``sdc_solver_tol_ener`` : the relative error of the energy in the nonlinear solve.

* ``sdc_solver_atol`` : the absolute error in the mass fractions during the nonlinear solve.

* ``sdc_solver_relax_factor`` : the factor by which to relax the
  tolerances (i.e. increase them) for earlier iterations.  We reach
  the desired tolerances on the final iteration.

* ``sdc_solve_for_rhoe`` : whether we solve the system in terms of :math:`(\rho e)` or :math:`(\rho E)`.

* ``sdc_newton_use_analytic_jac`` : whether we use the analytic Jacobian when doing Newton iterations for
  the reaction part of the system or compute it numerically.






.. _ch:mpiplusx:

******************************
Running Options: CPUs and GPUs
******************************

Castro uses MPI for coarse parallelization, distributing boxs across
compute nodes.  For fine-grained parallelism, OpenMP is used for
CPU-based computing and CUDA is used for GPUs.

Running on CPUs
===============

The preferred was of running on CPUs is to use MPI+OpenMP, compiling as::

  USE_MPI=TRUE
  USE_OMP=TRUE

Castro uses tiling to divide boxes into smaller tiles and distributes
these tiles to the OpenMP threads.  This is all managed at the MFIter
level -- no OpenMP directives need to be present in the compute
kernels themselves.  See `MFIter with Tiling
<https://amrex-codes.github.io/amrex/docs_html/Basics.html#sec-basics-mfiter-tiling>`_
for more information.

The optimal number of OpenMP threads depends on the computer
architecture, and some experimentation is needed.  Tiling works best
with larger boxes, so increasing ``amr.max_grid_size`` can benefit
performance.


Running on GPUs
===============

Castro's compute kernels can run on GPUs and this is the preferred way
to run on supercomputers with GPUs.  At the moment, offloading is
handled using CUDA and managed memory.  The exact same compute kernels
are used on GPUs as on CPUs.

.. note::

   Almost all of Castro runs on GPUs, with the main exception being
   the true SDC solver (``USE_TRUE_SDC = TRUE``).

To enable GPU computing, compile with::

  USE_MPI = TRUE
  USE_OMP = FALSE
  USE_CUDA = TRUE

When using GPUs, almost all of the computing is done on the GPUs.  In
the MFIter loops over boxes, the loops put a single zone on each GPU
thread, to take advantage of the massive parallelism.  The Microphysics
in StarKiller also takes advantage of GPUs, so entire simulations can
be run on the GPU.

Best performance is obtained with bigger boxes, so setting
``amr.max_grid_size = 128`` and ``amr.blocking_factor = 32`` can give
good performance.


Working at Supercomputing Centers
=================================

Our best practices for running any of the AMReX Astrophysics codes
at different supercomputing centers is produced in our workflow
documentation: https://amrex-astro.github.io/workflow/

.. _ch:buildsystem:

*********************
Build System Overview
*********************


Make Parameters
---------------

These build parameters control the parallelism, included physics,
etc.  In general they can be set to ``TRUE`` or ``FALSE``.  The
Castro-specific ones are interpreted by ``Make.Castro``.

A lot of optional features are enabled at compile time.  This allows
Castro to reduce the memory footprint of the state arrays by not allocating
space for variables that are not used.

General Build Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_ALL_CASTRO, USE_AMR_CORE, USE_HYPRE

These Parameters affect the build (parallelism, performance, etc.)
Most of these are parameters from AMReX.

  * ``USE_ALL_CASTRO``: compile all of the core Castro directories.
    This is the defailt (``TRUE``), and should not be changed for
    general simulations.  The purpose of this flag is for unit tests, which
    do not need all of the Castro directories compiled.  

  * ``USE_AMR_CORE``: compile all of the core AMReX directories, including
    ``Base/``, ``AmrCore/``, ``Amr/``, and ``Boundary/``.  This defaults
    to ``TRUE`` and should be left set for Castro simulations.  The purpose
    of this flag is for unit tests that don't need all of AMReX.

  * ``USE_MLMG``: use the AMReX multi-level multigrid solver for gravity
    and diffusion.  This should always be set to ``TRUE``.

  * ``USE_HYPRE``: compile in the Hypre library.  This will be automatically enabled
    for radiation.  You need to specify the path to the Hypre library via either
    ``HYPRE_DIR`` or ``HYPRE_OMP_DIR``.


Fortran Support
^^^^^^^^^^^^^^^

Many problems can be built without Fortran.  The current exceptions
are MHD, radiation, and anything using a pynucastro-generated network.
These parameters control Fortran support:

  * ``USE_FORT_MICROPHYSICS``: if set to ``TRUE``, then Fortran
    versions of the EOS and burner interface will be compiled.  If you
    are not using a pynucastro network, then you can probably set this
    to ``FALSE``.

  * ``BL_NO_FORT``: if set to ``TRUE``, then no AMReX Fortran source will be built.
    This cannot currently be used for the MHD or radiation solvers.


Parallelization and GPUs
^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_MPI, USE_OMP, USE_CUDA, USE_ACC

The following parameters control how work is divided across nodes, cores, and GPUs.

  * ``USE_CUDA``: compile with GPU support using CUDA. 

  * ``USE_ACC``: compile with OpenACC. Note: this is a work in
    progress and should not be used presently.


  * ``USE_MPI``: compile with the MPI library to allow for distributed parallelism.

  * ``USE_OMP``: compile with OpenMP to allow for shared memory parallelism.



General Physics Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. index:: USE_SIMPLIFIED_SDC, USE_TRUE_SDC

The following parameters control how the coupling between hydro and reactions
is handled.

  * ``USE_SIMPLIFIED_SDC``: use the simplified spectral deferred corrections (SDC)
    solver for coupling hydro and reactions.  At the moment, this
    works with the CTU hydrodynamics solver.  This requires running with
    ``castro.time_integration_method = 3``.

  * ``USE_TRUE_SDC``: use the true SDC method to couple hydro and
    reactions.  This can do 2nd order or 4th order accuracy.  At the
    moment, this works on single level only.  This requires running
    with ``castro.time_integration_method = 2``.



Radiation Parameters
^^^^^^^^^^^^^^^^^^^^

  * ``USE_RAD``: use photon radiation diffusion.  Note, For
    multigroup radiation, you need to set the number of radiation
    groups.  This is controlled by the ``NGROUPS`` parameter.

    .. index:: USE_RAD, NGROUPS



Gravity Parameters
^^^^^^^^^^^^^^^^^^

  * ``USE_GRAV``: use gravity (this could be constant or self-gravity)

    .. index:: USE_GRAV

  * ``USE_GR``: use a post-Newtonian approximation for GR gravity for the monopole
    solver.

    .. index:: USE_GR

  * ``USE_POINTMASS``: include a pointmass source to the gravitational potential.

    .. index:: USE_POINTMASS

Microphysics Parameters
^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_DIFFUSION``: enable thermal diffusion.  The conductivity is
    set via ``CONDUCTIVITY_DIR``, which should be a directory in the
    Microphysics repo.

    .. index:: USE_DIFFUSION, CONDUCTIVITY_DIR

  * ``USE_REACT``: enable reactions.  When reactions are set, we need
    to specify a network and an integrator.  Typically these come from
    the Microphysics repo, but one common exception is the
    ``general_null`` network, which just defines a composition.  The
    parameters that come into play here are:

    * ``NETWORK_DIR``: the network to use.  This is expected to be a subdirectory
      in the Microphysics repo.

    * ``NETWORK_INPUTS``: this is the text file that we read to define the
      composition if we are using the ``general_null`` network (e.g., ``gammalaw.net``).
      The build system will look for this file in the Microphysics repo.

    * ``INTEGRATOR_DIR``: this is the ODE integrator to use to integrate the 
      reaction system.  This is expected to be a subdirectory in the Microphysics
      repo.

    .. index:: USE_REACT, general_null, GENERAL_NET_INPUTS, NETWORK_DIR, INTEGRATOR_DIR

  * ``USE_REACT_SPARSE_JACOBIAN``

  * ``USE_SPARSE_STOP_ON_OOB``

  * ``EOS_DIR``: the equation of state to use.  This will be a subdirectory under the
    Microphysics repo.

    .. index:: EOS_DIR


Hydrodynamics and Source Term Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_ROTATION``: include rotation sources

    .. index:: USE_ROTATION

  * ``USE_HYBRID_MOMENTUM``: have Castro evolve angular momentum in addition to linear
    momentum.

    .. index:: USE_HYBRID_MOMENTUM

  * ``USE_SHOCK_VAR``: include a variable in the State_Type StateData that marks the
    location of a shock.

    .. index:: USE_SHOCK_VAR


Simulation Flow Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_AUX_UPDATE``: some networks define auxillary quantities, which in general
    Castro will advect, but not otherwise change.  If we set ``USE_AUX_UPDATE=TRUE``
    then Castro will call a user-supplied routine ``advance_aux()`` that can
    change the auxillary quantities.

    .. index:: USE_AUX_UPDATE

  * ``USE_POST_SIM``: if this is defined, then Castro will call the user-defined 
    routine ``problem_post_simulation()`` after the full evolution of the problem
    has ended.

    .. index:: USE_POST_SIM

  * ``USE_MAESTRO_INIT``: this enables the code to allow Castro to restart from a 
    Maestro simulation.  This will need to be updated in the future to allow for 
    restarts from MAESTROeX.

    .. index:: USE_MAESTRO_INIT

  * ``USE_HDF5``: compile in support for HDF5.  This is needed for some tables used
    by Microphysics routines.

    .. index:: USE_HDF5

Tracer Particle Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^

  * ``USE_PARTICLES``: compile in support for tracer particles.





Build Process Procedure
-----------------------

.. note::

   At build time, there are a number of source files that are autogenerated based
   on the configuration of the problem.  Most of these files are output into
   ``tmp_build_dir/castro_sources/Nd.COMP.OPTIONS.EXE/``, where ``N`` is the 
   dimensionality, ``COMP`` is the compiler name, and ``OPTIONS`` can be any
   number of options (``MPI``, ``DEBUG``, ...).

This is the current build system process.

* ``set_variables.py`` is called

  .. index:: set_variables.py, _variables, state_indices_nd.F90, state_indices.H

  * This processes the Castro ``_variables`` file and writes
    ``state_indices.H`` (and  ``state_indices_nd.F90`` if Fortran is enabled) into the
    ``tmp_build_dir/castro_sources/`` directory.

    These are used to define the size of the various state arrays and
    the integer keys to index each state variable.

  * The hook for this is in ``Make.auto_source`` in the build rule for ``state_indices_nd.F90``

  * You can test this portion of the build system by doing ``make test_variables``

* (for ``general_null networks``), ``network_properties.H`` (and
  ``actual_network.F90`` if Fortran is enabled) is created

  .. index:: write_network.py

  * This is done by ``write_network.py``

  * The hook for this is in ``$(CASTRO_HOME)/Microphysics/networks/general_null/Make.package``

* Runtime parameter files for the microphysics routines are parsed by ``write_probin.py``

  .. index:: write_probin.py

  * This writes the routines that manage the Microphysics runtime
    parameters: ``extern_parameters.cpp``, ``extern_parameters.H``, and  ``extern.F90``.  This is output in
    ``tmp_build_dir/castro_sources/``.

  * The hook for this is in ``Make.auto_source`` in the rule for ``extern_parameters.H``

* Castro's runtime parameters are parsed by ``parse_castro_params.py``

  .. index:: parse_castro_params.py

  * This writes the C++ header files that manage and read the runtime
    parameters.  These file are output in
    ``tmp_build_dir/castro_sources/``.

  * The hook for this is in ``Make.auto_source`` in the rule for ``castro_params.H``

* Problem-specific runtime parameters are parsed by ``write_probdata.py``

  * If the problem directory defines a ``_prob_params`` then it is parsed
    and used to C++ header and source files ``prob_parameters.H`` and ``prob_parameters.cpp``.
    These handle reading the ``problem.*`` parameters from the inputs file.
    Even without a problem-specific ``_prob_params``, all of the 
    variables in ``Castro/Source/problems/_default_prob_params`` will be included.

  * The script ``Castro/Util/scripts/write_probdata.py`` is used

  * The hook for this is in ``Make.auto_source`` in the ``prob_parameters.H`` rule.

  * These headers are output into ``tmp_build_dir/castro_sources/``.

* The Fortran dependencies file is created

  * This creates the ``f90.depends`` file in the ``tmp_build_dir``

  * The script ``amrex/Tools/F_scripts/dep.py`` is used

  * The hook for this is in ``amrex/Tools/GNUMake/Make.rules`` in the
    ``$(depEXETempDir)/f90.depends`` target

* The C/C++ dependencies file is created

  * This creates the individual ``.d`` files in ``tmp_build_dir``, one for each source file

  * A set of rules in ``Make.rules`` handles this. There is some
    description of what each line does in the comments of the make
    file

* Output to stdout the git version of the sources, via
  ``describe_sources.py``.  This doesn’t affect the build process

For all of this to work, we need the ``tmp_build_dir/s`` directory to
be first in the vpath, so our modified sources are found and used.


*************
Visualization
*************

There are a large number of tools that can be used to read in Castro
or AMReX data and make plots.  These tools all work from Castro
plotfiles.  Here we give an overview of the variables in plotfiles and
controlling their output, as well as some of the tools that can be used
for visualization.





Visualization Tools
===================

amrvis
------

Our favorite visualization tool is amrvis. We heartily encourage you
to build the amrvis2d and amrvis3d executables, and to try using them
to visualize your data. A very useful feature is View/Dataset, which
allows you to actually view the numbers – this can be handy for
debugging. You can modify how many levels of data you want to see,
whether you want to see the grid boxes or not, what palette you use,
etc.

If you like to have amrvis display a certain variable, at a certain
scale, when you first bring up each plotfile (you can always change it
once the amrvis window is open), you can modify the amrvis.defaults
file in your directory to have amrvis default to these settings every
time you run it. The directories CoreCollapse, HSE_test, Sod and
Sedov have amrvis.defaults files in them. If you are working in a new
run directory, simply copy one of these and modify it.

VisIt
-----

VisIt is also a great visualization tool, and it directly handles our
plotfile format (which it calls Boxlib). For more information check
out ``visit.llnl.gov``.

[Useful tip:] To use the Boxlib3D plugin, select it from File
:math:`\rightarrow` Open file :math:`\rightarrow` Open file as type Boxlib, and
then the key is to read the Header file, ``plt00000/Header``, for example,
rather than telling it to read ``plt00000``.

yt
--

yt is a free and open-source software that provides data analysis and
publication-level visualization tools for astrophysical simulation
results such as those Castro produces. As yt is script-based, it’s not
as easy to use as VisIt, and certainly not as easy as amrvis, but the
images can be worth it! Here we do not flesh out yt, but give an
overview intended to get a person started. Full documentation and
explanations from which this section was adapted can be found at
http://yt-project.org/doc/index.html.

Example notebook
^^^^^^^^^^^^^^^^

Using the plotfiles generated in the example in the :doc:`getting_started` section, here we demonstrate how to use ``yt`` to load and visualize data. This section was generated from a Jupyter notebook which can be found in ``Docs/source/yt_example.ipynb`` in the Castro repo. 

.. include:: yt_example.rst
.. _ch:rotation:

********
Rotation
********

Introduction
============

Currently, Castro supports contant, solid-body rotation about a fixed
(in space and time) axis in 2D and 3D by transforming the evolution
equations to the rotating frame of reference.

To include rotation you must set::

    USE_ROTATION = TRUE

in the ``GNUMakefile``. Rotation can then be enabled via::

    castro.do_rotation = 1

in the inputs file. The rotational period must then be set via
``castro.rotational_period``. The rotational period is internally
converted to an angular frequency for use in the source term
equations.

The axis of rotation currently depends on the dimensionality of the
problem and the value of coord_sys; in all cases, however, the
default axis of rotation points from ``problem::center`` in the vertical direction.

.. note:: make sure you have set the ``problem::center()`` variable
   appropriately for you problem.  This can be done by directly
   setting it in the ``problem_initialize()`` function.

The "vertical direction" is defined as follows:

* 2D

  * ``coord_sys = 0``, (x,y): out of the (x,y)-plane along the “z”-axis

  * ``coord_sys = 1``, (r,z): along the z-axis

* 3D

  * ``coord_sys = 0``, (x,y,z): along the z-axis

To change these defaults, modify the omega vector in the
``ca_rotate`` routine found in the ``Rotate_$(DIM)d.f90`` file.

The main parameters that affect rotation are:

-  ``castro.do_rotation`` : include rotation as a forcing
   term (0 or 1; default: 0)

-  ``castro.rotational_period`` : period (s) of rotation
   (default: 0.0)

-  ``castro.rotational_dPdt`` : d(period) / dt for rotation
   (default: 0.0)

-  ``castro.rotation_include_centrifugal`` : whether to
   include the centrifugal forcing (default: 1)

-  ``castro.rotation_include_coriolis`` : whether to
   include the Coriolis forcing (default: 1)

-  ``castro.rotation_include_domegadt`` : whether to
   include the forcing from the time derivative of the rotation
   frequency (default: 1)

-  ``castro.state_in_rotating_frame`` : whether state
   variables are measured in the rotating frame (default: 1)

-  ``castro.rot_source_type`` : method of updating the
   energy during a rotation update (default: 4)

-  ``castro.implicit_rotation_update`` : for the Coriolis
   term, which mixes momenta in the source term, whether we should
   solve for the update implicitly (default: 1)

-  ``castro.rot_axis`` : rotation axis (default: 3
   (Cartesian); 2 (cylindrical))

For completeness, we show below a derivation of the source terms that
appear in the momentum and total energy evolution equations upon
switching to a rotating reference frame.

Coordinate transformation to rotating frame
===========================================

.. figure:: tframes.png
   :alt: inertial vs. corotating frame

   Inertial frame :math:`C` and
   non-inertial frame :math:`\tilde{C}`. We consider a fluid element
   :math:`P`, whose distance in the two frames is related by
   :math:`{\bf r} = \tilde{\bf{r}} + {\bf l}`

Consider an inertial reference frame :math:`C` and a non-inertial
reference frame :math:`\widetilde{C}` whose origins are separated by
the vector :math:`\boldsymbol{l}` (see the figure above). The
non-inertial frame is rotating about the axis :math:`\ob` with a
*constant* angular velocity :math:`\omega`; furthermore, we assume the
*direction* of the rotational axis is fixed. Consider a fluid element
at the point :math:`P` whose location is given by :math:`\rb` in
:math:`C` and by :math:`\rbt` in :math:`\widetilde{C}`:

.. math:: \rb = \rbt + \boldsymbol{l},

or in component notation

.. math::
   r_i\boldsymbol{e_i} = \widetilde{r_i}\widetilde{\boldsymbol{e_i}} + l_i\boldsymbol{e_i},
   :label: eq:r

where :math:`\boldsymbol{e_i}` and :math:`\widetilde{\boldsymbol{e_i}}` are the :math:`i`\ th unit
vectors in the :math:`C` and :math:`\widetilde{C}` coordinate systems,
respectively. The total time rate of change of :eq:`eq:r` is given by

.. math::
   \frac{Dr_i}{Dt}\boldsymbol{e_i} = \frac{D\widetilde{r_i}}{Dt}\widetilde{\boldsymbol{e_i}} + \widetilde{r_i}\frac{D\widetilde{\boldsymbol{e_i}}}{Dt} + \frac{Dl_i}{Dt}\boldsymbol{e_i},
   :label: eq:vcomp


where we have used the fact that the unit vectors of the inertial
frame :math:`C` are not moving (or at least can be considered stationary,
and the change in :math:`\boldsymbol{l}` gives the relative motion of the two
coordinate systems). By definition, a unit vector can not change its
length, and therefore the only change of :math:`\widetilde{\boldsymbol{e_i}}` with
time can come from changing direction. This change is carried out by
a rotation about the :math:`\ob` axis, and the tip of the unit
vector moves circumferentially, that is

.. math::
   \frac{D\widetilde{\boldsymbol{e_i}}}{Dt} = \ob\times\widetilde{\boldsymbol{e_i}}.
   :label: eq:etilde-rot


Plugging :eq:`eq:etilde-rot` into :eq:`eq:vcomp` and switching back to
vector notation, we have

.. math::
   \frac{D\rb}{Dt} = \frac{D\rbt}{Dt} + \ob\times\rbt + \frac{D\boldsymbol{l}}{Dt}.
   :label: eq:r-dot


The left hand side of :eq:`eq:r-dot` is interpretted as the velocity
of the fluid element as seen in the inertial frame; the first term on the
right hand side is the velocity of the fluid element as seen by a
stationary observer in the rotating frame :math:`\widetilde{C}`. The second
and third terms on the right hand side of :eq:`eq:r-dot` describe the
additional velocity due to rotation and translation of the frame
:math:`\widetilde{C}` as seen in :math:`C`. In other words,

.. math::
   \vb = \vbt + \ob\times\rbt + \boldsymbol{v_l},
   :label: eq:v


where we use :math:`\boldsymbol{v_l}` to represent the translational velocity.

Similarly, by taking a second time derivative of :eq:`eq:v` we have

.. math::
   \frac{D\vb}{Dt} = \frac{D\vbt}{Dt} + 2\ob\times\vbt + \ob\times\left[\ob\times\rbt\right] + \frac{D\boldsymbol{v_l}}{Dt}.
   :label: eq:a


Henceforth we will assume the two coordinate systems are not
translating relative to one another, :math:`\boldsymbol{v_l} = 0`. It is
also worth mentioning that derivatives with respect to spatial
coordinates do not involve additional terms due to rotation,
i.e. :math:`\nablab\cdot\vb = \nablab\cdot\vbt`.
Because of this, the continuity equation remains unchanged in the
rotating frame:

.. math::
   \frac{\partial \rho}{\partial t} = -\nablab\cdot\left(\rho\vbt\right),
   :label: eq:cont-rot


or

.. math::
   \frac{D\rho}{Dt} = -\rho\nablab\cdot\vbt.
   :label: eq:cont-rot-total


Momentum equation in rotating frame
===================================

The usual momentum equation applies in an inertial frame:

.. math::
   \frac{D\left(\rho\vb\right)}{Dt} = -\rho\vb\cdot\nablab\vb - \nablab p + \rho\gb.
   :label: eq:mom1


Using the continuity equation, :eq:`eq:cont-rot-total`, and substituting for
the terms in the rotating frame from :eq:`eq:a`, we have from :eq:`eq:mom1`:

.. math::

   \begin{align}
       \rho\left(\frac{D\vbt}{Dt} + 2\ob\times\vbt + \ob\times\left[\ob\times\rbt\right]\right) - \rho\vb\nablab\cdot\vb &= -\rho\vb\cdot\nablab\vb - \nablab p + \rho\gb \nonumber \\
       \rho\left(\frac{\partial\vbt}{\partial t} + \vbt\cdot\nablab\vbt\right) &= -\nablab p + \rho\gb - 2\rho\ob\times\vbt - \rho\ob\times\left[\ob\times\rbt\right] \nonumber \\
     \frac{\partial\left(\rho\vbt\right)}{\partial t} &= -\nablab\cdot\left(\rho\vbt\vbt\right) - \nablab p + \rho\gb - 2\rho\ob\times\vbt \nonumber \\
     &-\ \rho\ob\times\left[\ob\times\rbt\right]\label{eq:mom-rot}
     \end{align}

or

.. math::
   \frac{D\left(\rho\vbt\right)}{Dt} = -\rho\vbt\cdot\nablab\vbt - \nablab p + \rho\gb - 2\rho\ob\times\vbt - \rho\ob\times\left[\ob\times\rbt\right].
   :label: eq:mom-rot-tot


Energy equations in rotating frame
==================================

From :eq:`eq:mom-rot-tot`, we have the velocity evolution equation in
a rotating frame

.. math::
   \frac{D\vbt}{Dt} = -\frac{1}{\rho}\nablab p + \gb - 2\ob\times\vbt - \ob\times\left[\ob\times\rbt\right].
   :label: eq:v-rot


The kinetic energy equation can be obtained from :eq:`eq:v-rot` by
mulitplying by :math:`\rho\vbt`:

.. math::
   \begin{align}
       \rho\vbt\cdot\frac{D\vbt}{Dt} &= -\vbt\cdot\nablab p + \rho\vbt\cdot\gb - 2\rho\vbt\cdot\left[\ob\times\vbt\right] - \rho\vbt\cdot\left\{\ob\times\left[\ob\times\rbt\right]\right\} \nonumber \\
       \frac{1}{2}\frac{D\left(\rho\vbt\cdot\vbt\right)}{Dt} - \frac{1}{2}\vbt\cdot\vbt\frac{D\rho}{Dt} &= -\vbt\cdot\nablab p + \rho\vbt\cdot\gb - \rho\vbt\cdot\left[\left(\ob\cdot\rbt\right)\ob - \rho\omega^2\rbt\right] \nonumber \\
       \frac{1}{2}\frac{D\left(\rho\vbt\cdot\vbt\right)}{Dt} &= -\frac{1}{2}\rho\vbt\cdot\vbt\nablab\cdot\vbt - \vbt\cdot\nablab p + \rho\vbt\cdot\gb - \rho\vbt\cdot\left[\left(\ob\cdot\rbt\right)\ob - \rho\omega^2\rbt\right].
     \end{align}
   :label: eq:ekin-rot-total

The internal energy is simply advected, and, from the first law of
thermodynamics, can change due to :math:`pdV` work:

.. math::
   \frac{D\left(\rho e\right)}{Dt} = -\left(\rho e + p\right)\nablab\cdot\vbt.
   :label: eq:eint-rot-total


Combining :eq:`eq:ekin-rot-total` and :eq:`eq:eint-rot-total` we can
get the evolution of the total specific energy in the rotating frame,
:math:`\rho \widetilde{E} = \rho e + \frac{1}{2}\rho\vbt\cdot\vbt`:

.. math::

   \begin{align}
       \frac{D\left(\rho e\right)}{Dt} + \frac{1}{2}\frac{D\left(\rho\vbt\cdot\vbt\right)}{Dt} &= -\left(\rho e + p + \frac{1}{2}\rho\vbt\cdot\vbt\right)\nablab\cdot\vbt - \vbt\cdot\nablab p \\
                     & + \rho\vbt\cdot\gb -\rho\vbt\cdot\left[\left(\ob\cdot\rbt\right)\ob - \rho\omega^2\rbt\right]\nonumber \\
       \frac{D\left(\rho \widetilde{E}\right)}{Dt} &= -\rho\widetilde{E}\nablab\cdot\vbt - \nablab\cdot\left(p\vbt\right) + \rho\vbt\cdot\gb - \rho\vbt\cdot\left[\left(\ob\cdot\rbt\right)\ob - \rho\omega^2\rbt\right] \label{eq:etot-rot-total}
     \end{align}

or

.. math::

   \label{eq:etot-rot}
       \frac{\partial\left(\rho\widetilde{E}\right)}{\partial t} = -\nablab\cdot\left(\rho\widetilde{E}\vbt + p\vbt\right) + \rho\vbt\cdot\gb - \rho\vbt\cdot\left[\left(\ob\cdot\rbt\right)\ob - \rho\omega^2\rbt\right].

Switching to the rotating reference frame
=========================================

If we choose to be a stationary observer in the rotating reference
frame, we can drop all of the tildes, which indicated terms in the
non-inertial frame :math:`\widetilde{C}`. Doing so, and making sure we
account for the offset, :math:`\boldsymbol{l}`, between the two coordinate systems, we obtain
the following equations for hydrodynamics in a rotating frame of
reference:

.. math::

   \begin{align}
       \frac{\partial\rho}{\partial t} &= -\nablab\cdot\left(\rho\vb\right) \label{eq:cont-rot-switch} \\
       \frac{\partial \left(\rho\vb\right)}{\partial t} &= -\nablab\cdot\left(\rho\vb\vb\right) - \nablab p + \rho\gb - 2\rho\ob\times\vb - \rho\left(\ob\cdot\rb\right)\ob + \rho\omega^2\rb \label{eq:mom-rot-switch} \\
       \frac{\partial\left(\rho E\right)}{\partial t} &= -\nablab\cdot\left(\rho E\vb + p\vb\right) + \rho\vb\cdot\gb - \rho\left(\ob\cdot\rb\right)\left(\ob\cdot\vb\right) + \rho\omega^2\left(\vb\cdot\rb\right). \label{eq:etot-rot-switch}
     \end{align}

Adding the forcing to the hydrodynamics
=======================================

There are several ways to incorporate the effect of the rotation
forcing on the hydrodynamical evolution. We control this through the
use of the runtime parameter castro.rot_source_type. This
is an integer with values currently ranging from 1 through 4, and
these values are all analogous to the way that gravity is used to
update the momentum and energy. For the most part, the differences are
in how the energy update is done:

* ``castro.rot_source_type = 1`` : we use a standard
  predictor-corrector formalism for updating the momentum and
  energy. Specifically, our first update is equal to :math:`\Delta t \times \mathbf{S}^n` ,
  where :math:`\mathbf{S}^n` is the value of
  the source terms at the old-time (which is usually called time-level
  :math:`n`). At the end of the timestep, we do a corrector step where
  we subtract off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on
  :math:`\Delta t / 2 \times \mathbf{S}^{n+1}`, so that at the end of
  the timestep the source term is properly time centered.

* ``castro.rot_source_type = 2`` : we do something very similar
  to 1. The major difference is that when evaluating the energy source
  term at the new time (which is equal to
  :math:`\mathbf{u} \cdot \mathbf{S}^{n+1}_{\rho \mathbf{u}}`, where the latter is the
  momentum source term evaluated at the new time), we first update the
  momentum, rather than using the value of :math:`\mathbf{u}` before
  entering the rotation source terms. This permits a tighter coupling
  between the momentum and energy update and we have seen that it
  usually results in a more accurate evolution.

* ``castro.rot_source_type = 3`` : we do the same momentum update as
  the previous two, but for the energy update, we put all of the work
  into updating the kinetic energy alone. In particular, we explicitly
  ensure that :math:`(rho e)` maintains the same, and update
  :math:`(rho K)` with the work due to rotation, adding the new
  kinetic energy to the old internal energy to determine the final
  total gas energy. The physical motivation is that work should be
  done on the velocity, and should not directly update the temperature
  – only indirectly through things like shocks.

* ``castro.rot_source_type = 4`` : the energy update is done in a
   “conservative” fashion. The previous methods all evaluate the value
   of the source term at the cell center, but this method evaluates
   the change in energy at cell edges, using the hydrodynamical mass
   fluxes, permitting total energy to be conserved (excluding possible
   losses at open domain boundaries). Additionally, the velocity
   update is slightly different—for the corrector step, we note that
   there is an implicit coupling between the velocity components, and
   we directly solve this coupled equation, which results in a
   slightly better coupling and a more accurate evolution.

The other major option is ``castro.implicit_rotation_update``.
This does the update of the Coriolis term in the momentum equation
implicitly (e.g., the velocity in the Coriolis force for the zone
depends on the updated momentum). The energy update is unchanged.

A detailed discussion of these options and some verification
tests is presented in :cite:`katz:2016`.
************************
Timestepping and Retries
************************

Simulation Time
---------------

There are two paramters that can define when a simulation ends:

  * ``max_step``: maximum number of level 0 time steps (integer
    :math:`\geq 0`; default: -1)

  * ``stop_time``: final simulation time (Real :math:`\geq 0`; default:
    -1.0)

To control the number of time steps, you can limit by the maximum
number of level 0 time steps (``max_step``) or by the final
simulation time (``stop_time``), or both. The code will stop at
whichever criterion comes first.

Note that if the code reaches ``stop_time`` then the final time
step will be shortened so as to end exactly at ``stop_time``, not
past it.

As an example::

    max_step  = 1000
    stop_time  = 1.0

will end the calculation when either the simulation time reaches 1.0 or
the number of level 0 steps taken equals 1000, whichever comes first.

Time Step Constraints
---------------------

Hydrodynamics
^^^^^^^^^^^^^

If ``castro.do_hydro = 1``, then typically
the code chooses a time step based on the CFL number:

.. math::
   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\min\left\{\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}},
                                                                  \frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right\}\right]
   :label: eq:cfl

If SDC integration is used instead, then we have

.. math::

   \Delta t = \mathtt{CFL}\, \cdot\, \min_{i,j,k}\left[\left(\frac{\Delta x}{|u|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta y}{|v|_{i,j,k}+c_{i,j,k}}\right)^{-1} +
                                                       \left(\frac{\Delta z}{|w|_{i,j,k}+c_{i,j,k}}\right)^{-1}\right]^{-1}

(If we are simulating in 1D or 2D, the extraneous parts related to :math:`v` and/or :math:`w` are removed.)

Additional Controls
^^^^^^^^^^^^^^^^^^^

The following parameters affect the timestep choice:

  * ``castro.cfl``: CFL number (Real :math:`> 0` and :math:`\leq 1`;
    default: 0.8)

  * ``castro.init_shrink``: factor by which to shrink the initial
    time step (Real :math:`> 0` and :math:`\leq 1`; default: 1.0)

  * ``castro.change_max``: factor by which the time step can
    grow in subsequent steps (Real :math:`\geq 1`; default: 1.1)

  * ``castro.fixed_dt``: level 0 time step regardless of cfl
    or other settings (Real :math:`> 0`; unused if not set)

  * ``castro.initial_dt``: initial level 0 time
    step regardless of other settings (Real :math:`> 0`; unused if not set)

  * ``castro.dt_cutoff``: as a fraction of the current simulation time,
    the time step below which the calculation will abort (Real
    :math:`> 0`; default: 1.e-12); typically not user-defined

As an example, consider::

    castro.cfl = 0.9
    castro.init_shrink = 0.01
    castro.change_max = 1.1

This defines the :math:`\mathtt{cfl}` parameter in :eq:`eq:cfl` to be
0.9, but sets (via ``init_shrink``) the first timestep we take to
be 1% of what it would be otherwise. This allows us to ramp up to
the hydrodynamic timestep at the start of a simulation. The
``change_max`` parameter restricts the timestep from increasing by
more than 10% over a coarse timestep. Note that the time step can
shrink by any factor; this only controls the extent to which it can
grow. The ``dt_cutoff`` parameter will force the code to abort if
the timestep ever drops below :math:`10^{-12}` of the current time. This is a safety
feature—if the code hits such a small value, then something likely
went wrong in the simulation, and by aborting, you won’t burn through
your entire allocation before noticing that there is an issue.

If we know what we are doing, then we can force a particular timestep::

    castro.fixed_dt = 1.e-4

This sets the level 0 time step to be 1.e-4 for the entire simulation,
ignoring the other timestep controls. Note that if
``castro.init_shrink`` :math:`\neq 1` then the first time step will in fact
be ``castro.init_shrink`` :math:`\cdot` ``castro.fixed_dt``.

::

    castro.initial_dt = 1.e-4

sets the *initial* level 0 time step to be :math:`10^{-4}` regardless of
``castro.cfl`` or ``castro.fixed_dt``. The time step can
grow in subsequent steps by a factor of castro.change_max each step.


Diffusion
^^^^^^^^^

If diffusion is enabled, the timestep will also be limited by:

.. math::

   \Delta t = \frac{1}{2}\min_{i,j,k}\left[\min\left\{\frac{\Delta x^2}{D_{i,j,k}},
                                                      \frac{\Delta y^2}{D_{i,j,k}},
                                                      \frac{\Delta z^2}{D_{i,j,k}}\right\}\right]

where :math:`D \equiv k / (\rho c_V)` if we are diffusing temperature,
and :math:`D \equiv k / (\rho c_P)` if we are diffusing enthalpy. No
input parameter is necessary to enable this constraint. See Chapter
:ref:`ch:diffusion` for more details.

Reactions
^^^^^^^^^

If reactions are enabled, the timestep will also
be limited by two constraints:

.. math:: \Delta t = \mathtt{dtnuc\_e}\, \min_{i,j,k} \left\{\frac{e_{i,j,k}}{\dot{e}_{i,j,k}}\right\}

.. math:: \Delta t = \mathtt{dtnuc\_X}\, \min_{i,j,k} \left\{\min_n\frac{X^n_{i,j,k}}{\dot{X}^n_{i,j,k}}\right\}

where :math:`e` is the internal energy, and :math:`X^n` is the mass fraction of
the :math:`n`\ th species. The safety factors correspond to the runtime parameters
``castro.dtnuc_e`` and ``castro.dtnuc_X``. These limiters
say that the timestep must be small enough so that no zone can change
its internal energy by more than the fraction in one
step, and so that no zone can change the abundance of any isotope by
more than the fraction in one step. The time derivatives
:math:`\dot{e}` and :math:`\dot{X}^n` are estimated by calling the right-hand-side
of the nuclear network given the state at the time the timestep limiter
is being calculated. (We use a small number floor to prevent division by zero.)
To prevent the timestep from being dominated by trace species, there is
an additional option ``castro.dtnuc_X_threshold`` which is the
mass fraction threshold below which a species will not be considered in
the timestep constraint. and are set to
a large number by default, effectively disabling them. Typical choices
for these values in the literature are :math:`\sim 0.1`.

Subcycling
----------

Subcycling with AMR means that coarser grids can take a larger timestep
than finer grids.  
Castro supports a number of different modes for subcycling in time,
set via ``amr.subcycling_mode``.

  * ``amr.subcycling_mode`` = ``Auto`` (default): the code will run with
    equal refinement in space and time. In other words, if level
    :math:`n+1` is a factor of 2 refinement above level :math:`n`,
    then :math:`n+1` will take 2 steps of half the duration for every
    level :math:`n` step.

  * If ``amr.subcycling_mode`` = ``None``: the code will not refine in
    time. All levels will advance together with a timestep dictated by
    the level with the strictest :math:`dt`. Note that this is
    identical to the deprecated command ``amr.nosub = 1``.

  * If ``amr.subcycling_mode`` = ``Manual``: the code will subcycle
    according to the values supplied by ``amr.subcycling_iterations``.

In the case of ``amr.subcycling_mode`` = Manual, we subcycle in
manual mode with largest allowable timestep. The number of iterations
at each level is then specified as::

    amr.subcycling_iterations = 1 2 1 2

Here, we take 1 level-0 timestep at a time (required). Take 2 level-1
timesteps for each level-0 step, 1 timestep at level-2 for each
level-1 step, and take 2 timesteps at level-3 for each level-2 step.

Alternately, we could do::

    amr.subcycling_iterations = 2

which will subcycle twice at every level (except level 0).


.. index:: retry

.. _ch:retry:

Retry Mechanism
---------------

Castro's Strang CTU solver has a retry mechanism that can discard a
time step on a level and restart with a smaller timestep, subcycling
within the level to make up the full time step needed for that level.
It is enabled by setting::

   castro.use_retry = 1

.. note::

   The Castro retry mechanism is enabled by default for CTU + Strang
   and Simplified SDC integration.

The number of subcycles to try in the level is controlled via the
``castro.max_subcycles`` parameter.  It is not really suggested to go
beyond ``16``---any more is usually an indication of a bigger problem.

A retry can be triggered by a number of conditions:

  * Exceeding the CFL condition for a level

  * A negative density is encountered

  * Integration failure in the burner

    Note: this requires that the following be set in your ``&extern``
    namelist::

      retry_burn = F
      abort_on_failure = F

    This instructs the integration routine in Microphysics to not
    abort when the integration fails, but instead to tell the calling
    Castro routine that the integration failed so Castro can handle
    the retry itself.

    .. note::

       The combination of ``use_retry = 0`` and ``abort_on_failure = F``
       is unsafe and not supported.

       For true SDC, we disable retry and reset ``abort_on_failure`` to
       always be true, since retry is not supported for that integration.


***************************
Setting Up Your Own Problem
***************************

Castro problems are organized loosely into groups describing their
intent (e.g., science, hydro tests, ...).  These groups are
sub-directories under the ``Castro/Exec/`` directory.  Each problem is
then placed in a sub-directory of the appropriate group (for example,
``Castro/Exec/hydro_tests/Sedov`` holds the Sedov test problem).

To create a new problem, you will create a new directory under one
of the groups and place in it the following files:

  * ``GNUmakefile`` : the makefile for this problem.  This will tell
    Castro what options to use and what network and EOS to build.

  * ``problem_initialize.H`` and
    ``problem_initialize_state_data.H`` : this holds the problem
    initialization routines.  MHD and radiation problems require
    an additional file.

  * ``_prob_params`` (optional) : a list of runtime parameters that
    you problem will read.  These parameters are controlled by the
    set in the inputs file as ``problem.param`` where ``param`` is
    one of the parameters listed in ``_prob_params``.

  * ``Make.package`` : this is a makefile fragment that is included
    during the build process.  It tells the build system about any
    problem-specific files that need to be compiled.

  * ``inputs`` : this is the main inputs file that controls Castro and
    AMReX's behavior.

The best way to get started writing a new problem is to copy an
existing problem setup into a new directory.

.. index:: _prob_params

Runtime Parameters
------------------

The problem-specific runtime parameters are defined in ``_prob_params``.
This has the form::

   name       datatype    default      namelist?     size

Here:

* `name` is the name of the variable to put into ``probdata_module``

* `datatype` is one of ``real``, ``integer``, ``character``, or
  ``logical``.

* `default` is the default value of the runtime parameter.  It may be
  overridden at runtime by reading from the namelist.

* `namelist` indicates if the variable should be able to be set as
  a runtime parameter.  If it is empty or marked with
  ``n``, then the variable will still be creates in the ``problem`` namespace,
  but it will not be able to be set via the commandline or inputs file.
  A common usage of this is to define global variables that might be set
  at problem initialization that are used elsewhere in Castro.

* `size` is for arrays, and gives their size.  It can be any integer
  or variable that is known to the ``probdata_module``.  If you need a
  variable from a different module to specify the size (for example,
  ``nspec`` from ``network``), then you express the size as a tuple:
  ``(nspec, network)`` and the appropriate use statement will be
  added.

The variables will all be initialized for the GPU as well.


Problem Initialization
----------------------

Here we describe the main problem initialization routines. 

.. index:: initialize_problem

* ``initialize_problem()``

  This C++ routine is primarily responsible for doing any one-time
  initialization for the problem (like reading in an
  initial model through the ``model_parser.H`` functionality—see the
  ``toy_convect`` problem setup for an example.

  This routine can also postprocess any of the parameters defined
  in the ``_prob_params`` file, which are defined in ``problem`` namespace.

  .. note:: many problems set the value of the ``problem::center[]`` array
     from the ``problem`` namespace.  This is used to note the
     center of the problem (which does not necessarily need to be
     the center of the domain, e.g., for axisymmetric problems).
     ``center`` is used in source terms (including rotation and
     gravity) and in computing some of the derived variables (like
     angular momentum).

  If you need coordinate information, it can be obtained 
  by constructing a ``Geometry`` object using ``DefaultGeometry()``
  and accessing its ``ProbLo()`` and ``ProbHi()`` methods.


* ``problem_initialize_state_data()``:

  This routine will initialize the state data in a given zone.
  The arguments passed are:

  - ``i``, ``j``, ``k``: the index of the zone to fill the data in

  - ``state``: an array containing the simulation state data

  - ``GeomData``: a ``GeometryData`` object that can be used for obtaining
    ``dx``, ``problo``, ``probhi``, etc.

  Filling data is done by simply writing to ``state(i,j,k,URHO)``, etc.

.. _create:bcs:

Boundary Conditions
-------------------

.. index:: boundary conditions

Standard boundary conditions, including outflow (zero-gradient), periodic,
and symmetry (reflect) are handled by AMReX directly.  Castro has a special
hydrostatic boundary condition that can be used for the lower boundary.  It
is accessed by setting the ``castro.lo_bc`` flag to 1 in the vertical coordinate
direction, e.g., for 2-d as::

   castro.lo_bc       =  0   1

The flag value 1 is traditionally named "inflow" by AMReX, but generally means that
the boundary implementation is left to the user.  To tell Castro to use the
hydrostatic boundary condition here, we set::

   castro.yl_ext_bc_type = 1
   castro.hse_interp_temp = 1
   castro.hse_reflect_vels = 1

The first parameter tells Castro to use the HSE boundary condition for the lower
y direction.
In filling the ghost cells, hydrostatic equilibrum will be integrated
from the last interior zone into the boundary.  We need one more
equation for this integration, so we either interpolate the density or
temperature into the ghost cells, depending on the value of
``castro.hse_interp_temp``.  Finally, ``castro.hse_reflect_vels``
determines how we treat the velocity.  The default is to give is a
zero gradient, but in tests we've found that reflecting the velocity
while integrating the HSE profile can be better.  For modeling a
plane-parallel hydrostatic atmosphere, using the hydrostatic boundary
conditions instead of a simple symmetry boundary is essential when
using the standard CTU PPM solver.

A different special boundary condition, based on outflow, is available at
the upper boundary.  This works together with the ``model_parser``
module to fill the ghost cells at the upper boundary with the initial
model data.  You set this as::

   castro.hi_bc = 2 2

   castro.fill_ambient_bc = 1
   castro.ambient_fill_dir = 1
   castro.ambient_outflow_vel = 1

where ``ambient_fill_dir`` is the 0-based direction to fill using an
ambient state defined by the problem setup.  In this example, we will
override the outflow (2) boundary condition in the y-direction.  That
problem setup needs to fill the ``ambient_state[:]`` array defined in
``ambient.H``.  An example of using this boundary is in the
``flame_wave`` problem.

The implementations of these boundary conditions is found in
``Castro/Source/problems/Castro_bc_fill_nd.cpp``.

Optional Files
--------------

The follow problem-specific files are optional. There are stubs for
each of these in the main source tree.

-  ``problem_checkpoint.H``, ``problem_restart.H`` :

   These provides two routines, respectively ``problem_checkpoint`` and
   ``problem_restart`` that can be used to add information to the
   checkpoint files and read it in upon restart. This is useful for
   some global problem-specific quantities. For instance, the
   ``wdmerger`` problem uses this to store center of mass position and
   velocity information in the checkpoint files that are used for
   runtime diagnostics.

   The name of the checkpoint directory is passed in as an argument.

-  ``problem_tagging.H``

   This implements problem-specific tagging for refinement, through a
   the function ``problem_tagging``. The full hydrodynamic state (State_Type)
   is passed in, and the problem can mark zones for refinement by setting the
   tag variable for a zone to set. An example is provided by the ``toy_convect``
   problem which refines a rectangular region (fuel layer) based on
   a density parameter and the H mass fraction.

-  ``Problem_Derives.H``, ``Problem_Derive.H``, and ``Problem_Derives.cpp``

   Together, these provide a mechanism to create derived quantities
   that can be stored in the plotfile. ``Problem_Derives.H``
   provides the C++ code that defines these new plot variables. It
   does this by adding them to the ``derive_lst``—a list of
   derived variables that Castro knows about. When adding new
   variables, a descriptive name, Fortran routine that does the
   deriving, and component of ``StateData`` are specified.

   The other two files provide the header and implementation of the
   function that computes the derived variable.  A example is provided
   by the ``reacting_bubble`` problem, which derives several new
   quantities (perturbations against a background one-dimensional
   model, in this case).

-  ``Prob.cpp``, ``Problem.H``

   These files provide problem-specific routines for computing global
   diagnostic information through the sum_integrated_quantities
   functionality that is part of the ``Castro`` class.

   An example is provided by ``toy_flame``, where an estimate
   of the flame speed is computed by integrating the mass of fuel on
   the grid.


Model Parser
------------

Many problem setups begin with a 1-d initial model that is mapped onto
the grid.  The ``model_parser.H`` provides the functions that read in
the initial model and map it on the Castro grid.  To enable this, add::

  USE_CXX_MODEL_PARSER = TRUE

to the problem ``GNUmakefile``.  There are 2 other parameters that can
be set in the makefile to control the initial model storage:

  * ``MAX_NPTS_MODEL``: is the maximum number of data points in the
    1-d initial model.  This needs to be known at compile time so we
    can make the data managed for GPUs.

  * ``NUM_MODELS``: this is the number of different initial models we
    want to managed.  Typically we only want 1, but some problems,
    like ``flame_wave`` use 2, applied to different portions of the
    domain.

The general form of the initial model is::

    # npts = 896
    # num of variables = 6
    # density
    # temperature
    # pressure
    # carbon-12
    # oxygen-16
    # magnesium-24
    195312.5000  5437711139.  8805500.952   .4695704813E+28  0.3  0.7  0
    585937.5000  5410152416.  8816689.836  0.4663923963E+28  0.3  0.7  0

The first line gives the number of points in the initial model, the
next gives the number of variables, followed by the variable names
(one per line), and then the data.  The data begins with the
coordinate and then the variables in the model, with one data point
per line.

When the model is read, the variables listed in the file are matched
to the ones that Castro knows about.  If the variable is recognized,
then it is stored in the model data, otherwise, it is ignored.

The data can then be mapped onto the grid using the ``interpolate()``
function, e.g., ::

    Real dens = interpolate(height, model::idens);

This fills ``dens`` with the density at the position ``height``.  In
addition to density, you can specify temperature (``model::itemp``),
pressure (``model::ipres``), species (indexed from ``model::ispec``),
or an auxiliary quantity (indexed from ``model::iaux``).


***********************
Restarting from Maestro
***********************

Overview
========

We can now initialize a Castro simulation using data from a Maestro plotfile. This should not be thought of as a restart mode, but rather
a new simulation with a special initialization. In order to use this
feature, you must make sure the Maestro plotfile has the proper
variables, add some new parameters to your inputs file, and add a few
subroutines to Prob_Xd.f90. You need to build a special executable
with “USE_MAESTRO_INIT=TRUE”, which will add “.MAESTRO” to the
executable string. For multilevel problems, there are a few extra
steps relating to the fact that you have to supply a grids file
consistent with the Maestro grid structure.

MAESTRO Plotfile Requirements
=============================

The Maestro plotfile needs to have the following variables:

-  “x_vel”, “y_vel”, (and “z_vel”, depending on
   dimensionality of the problem)

-  “density” (**castro.MAESTRO_init_type** = 1 and 2 only)

-  Optional species (such as “X(C12)”) - there is an option to
   not read any species from the Maestro plotfile. In this case, you
   must make sure your code manually defines the species cell-by-cell
   in the initial Castro data

-  “tfromp”

-  “pi” (**castro.MAESTRO_init_type** = 2, 3, and 4 only)

-  “entropy” (**castro.MAESTRO_init_type** = 4 only)

Also, model_cc_XXXXX needs to list variables in the following order,
which is the default order found in MAESTRO/Source/base_io.f90: r,
base_r, rho0, p0, gamma1bar, rhoh0, div_coeff, psi, tempbar,
etarho_cc, tempbar_init.

List of Parameters
==================

Here are the additional parameters you must add to your inputs file.

+-----------------------------------+------------------+-----------------+-----------------+
| Parameter                         | Definition       | Type            | Default         |
+===================================+==================+=================+=================+
| **castro.MAESTRO_plotfile**       | name of the      | std::string     | must be set     |
|                                   | Maestro plotfile |                 |                 |
|                                   |                  |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_modelfile**      | name of the      | std::string     | must be set     |
|                                   | Maestro model_cc |                 |                 |
|                                   | file             |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_npts_model**     | number of        | int             | must be set     |
|                                   | points in the    |                 |                 |
|                                   | Maestro model_cc |                 |                 |
|                                   | file             |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_first_species**  | name of the      | std::string     | must be set or  |
|                                   | first species    |                 | else nothing    |
|                                   |                  |                 | will be read in |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_nspec**          | number of        | std::string     | NumSpec in      |
|                                   | species in the   |                 | Castro          |
|                                   | Maestro plotfile |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_cutoff_density** | controls how we  | Real            | must be set     |
|                                   | overwrite data   |                 |                 |
|                                   | at the edge of   |                 |                 |
|                                   | the star         |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_init_type**      | determines how   | int             | must be set     |
|                                   | we initialize    |                 |                 |
|                                   | the              |                 |                 |
|                                   | Castro state     |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+
| **castro.MAESTRO_spherical**      | specifies        | int             | must be set     |
|                                   | planar or        |                 |                 |
|                                   | spherical        |                 |                 |
|                                   | problem          |                 |                 |
+-----------------------------------+------------------+-----------------+-----------------+

Examples of Usage
-----------------

-  **castro.MAESTRO_plotfile** = "wd_384_6.25e8K_norotate_plt120578"

-  **castro.MAESTRO_modelfile** = "./wd_384_6.25e8K_norotate_plt120578/model_cc_120578"

-  | **castro.MAESTRO_npts_model** = 1663
   | This is the number of
     points in **castro.MAESTRO_modelfile**. Note that this is not
     the same thing as “npts_model”, which is the number of points in
     the initial model file used for standard simulations where we do not
     initialize from a Maestro plotfile.

-  **castro.MAESTRO_first_species** = “X(C12)” If you do not
   specify this, no species will be read in. You can always manually
   specify or overwrite the species cell-by-cell later.

-  | **castro.MAESTRO_nspec** = 3
   | If you do not specify this, it
     will default to the number of species in the Castro network,
     “NumSpec”. We have this here because sometimes Maestro and Castro will use different networks with different number of species.

-  | **castro.MAESTRO_cutoff_density** = 1.e6
   | The code will use
     this density to figure out the radial coordinate, r_model_start,
     which is the last radial coordinate before rho0 falls below
     **castro.MAESTRO_cutoff_density**. It is possible to set
     **castro.MAESTRO_cutoff_density** to a tiny value, such that rho0
     never falls below this value, in which case we set r_model_start
     to :math:`\infty`. In INITDATA_MAKEMODEL, we create a new 1D model
     integrating outward starting from r_model_start. Then, in
     INITDATA_OVERWRITE, we overwrite newly initialized Castro data in
     any cell that maps into a radial coordinate greater than
     r_model_start by interpolating from the new 1D model.

-  | **castro.MAESTRO_init_type** = 2
   | Castro will read in data
     from the Maestro plotfile, and then call the EOS to make sure that
     :math:`\rho`, :math:`e`, :math:`T`, and :math:`X_k` are consistent. The inputs to the EOS
     are based on the value of **castro.MAESTRO_init_type**:

   #. :math:`e = e(\rho,T,X_k)`

   #. :math:`e,T = e,T(\rho,p_0+\pi,X_k)`

   #. :math:`\rho,e = \rho,e(p_0+\pi,T,X_k)`

   #. :math:`\rho,T,e = \rho,T,e(p_0+\pi,s,X_k)`

-  | **castro.MAESTRO_spherical** = 1
   | 0 = planar; 1 = spherical.

New Subroutines in Prob_Xd.f90
==============================

There are three routines that need to be added to your local copy of
Prob_Xd.f90. See Castro/Exec/wdconvect/Prob_3d.f90 for
a standard spherical Maestro initialization.

#. | INITDATA_MAESTRO
   | This fills in the Castro state by taking
     the Maestro data, calling the EOS, and making the proper variables
     conserved quantities. Specifically, we need a thermodynamically
     consistent :math:`\rho`, :math:`T`, :math:`e`, and :math:`X_k`, and then algebraically
     compute :math:`\rho{\bf u}`, :math:`\rho e`, :math:`\rho E`, and :math:`\rho X_k`,

#. | INITDATA_MAKEMODEL
   | This creates a user-defined 1D initial model starting from r_model_start.

#. | INITDATA_OVERWRITE
   | This overwrites the initialized Castro data using the new 1D initial model for all cells that map into
     radial coordinates greater than r_model_start.

Additional Notes
================

Note that for both single-level and multilevel Maestro to Castro initialization, the Castro base grid structure does not have to match
the Maestro base grid structure, as long as the problem domain is the
same. For example, if the coarsest level in a Maestro plotfile
contains :math:`64^3` cells divided into 8-\ :math:`32^3` grids, it is ok to use a
Castro base grid structure with 1-\ :math:`64^3` grid, 64-\ :math:`16^3` grids, or
anything else you can imagine - the grids don’t even have to be the
same size. As is normally the case, the Castro base grid structure is
created based on the parameters in the Castro inputs file, such as
**amr.max_grid_size**, **amr.blocking_factor**, etc.

Multilevel Restart
------------------

When initialing from a multilevel Maestro plotfile, there are some
extra steps. First, you need to create a Castro-compatible grids file
from the Maestro plotfile. This can be done with the
BoxLib/Tools/Postprocessing/F_Src/fboxinfo.f90 utility. Compile
and run this using the “``–``\ castro” option, e.g.,
“fboxinfo.Linux.gfortran.exe ``–``\ castro pltxxxxx ``|``
tee gr0.maestro”, to generate the Castro-compatible grids file. Note
that the base grid structure is still controlled by
``amr.max_grid_size``, ``amr.blocking_factor``, etc., since in C++ AMReX, the grids file only indicates the refined grid structure,
whereas in Fortran BoxLib the grids file contains the base grid and
refined grid structures.

Now, when you initialize the Castro simulation, you need to specify
the grid file using **amr.regrid_file = "gr0_3d.128_2levels"**,
for example. You can happily run this now, but note that the
regridding algorithm will never be called (since Castro thinks it’s
started a new simulation from scratch with a grids file, thus
disabling the regridding). If you wish for the grid structure to be
changed, you must do a traditional Castro restart from the
Castro-generated checkpoint file (you can still use the same
“.MAESTRO” executable or an executable built with
USE_MAESTRO_INIT=FALSE), making sure that you **do not** specity
**amr.regrid_file** (or else the grids will stay fixed). You are
free to specify **amr.regrid_on_restart**,
**amr.compute_new_dt_on_regrid**, and
**amr.plotfile_on_restart**.

Sometimes a Maestro plotfile will only have 1 or 2 total levels, but
you ultimately want to run a Castro simulation with many more levels
of refinement. My recommended strategy is the following:

#. Initialize a Castro simulation from the Maestro plotfile
   while preserving the exact same grid structure and run for 10 time
   steps.

#. Do a traditional Castro restart from chk00010, but do not
   increase **amr.max_level**, and run for 10 more time steps. This
   allows a new grid structure with the same effective resolution as
   before settle in using the C BoxLib regridding algorithm.

#. Do a traditional Castro restart from chk00020, but increase
   **amr.max_level** by 1, and run for 10 time steps.

#. Repeat the procedure from the previous step (using the most
   updated checkpoint of course) as many times as desired.
.. _ch:mhd:

***
MHD
***

Introduction
============

Castro implements a constrained transport (CT) corner transport upwind
(CTU) ideal MHD scheme based on the work of Miniati & Martin
:cite:`miniati_martin`.  MHD is enabled by compiling with ``USE_MHD =
TRUE``.  This replaces the pure hydrodynamics solver, but uses the
same driver as the CTU hydrodynamics solver.  This means that all of
the source terms supported by the hydrodynamics solver are also
supported by MHD.

.. note::

   The MHD solver supports 3-d only.

   Currently the MHD solver is single-level only.  AMR support is forthcoming.

Equations and Data Structures
=============================

The ideal MHD equations we solve appear as:

.. math::

   \begin{align}
   \frac{\partial \rho X_k}{\partial t} + \frac{\partial}{\partial x_j} ( \rho U_j X_k) &= \rho \omegadot_k \\
   \frac{\partial \rho U_j}{\partial t} + \frac{\partial}{\partial x_j} (\rho U_i U_j + p \delta_{ij} - B_j B_i) &= S_{\rho U_j} \\
   \frac{\partial \rho E}{\partial t} + \frac{\partial}{\partial x_j} \left [ U_j (\rho E + p) - B_j B_i u_i \right ] &= S_{\rho E} + \rho H_\mathrm{nuc} \\
   \frac{\partial B_i}{\partial t} &= -\frac{\partial}{\partial x_j} (U_j B_i - B_j U_i)
   \end{align}

where

.. math::

   p = p_g + \frac{1}{2} B^2

and the :math:`S` sources represent the hydrodynamical sources and
the remainder are reaction sources.  The MHD solver uses the same
time-advancement driver as the hydrodynamic CTU driver
(:ref:`sec:strangctu`), using Strang splitting for the reactions (by
default).

.. note::

   Note: we are following the convention in the MHD community of setting the permeability to unity.  This is
   why the magnetic pressure has the form :math:`\frac{1}{2} B^2` instead of :math:`\frac{1}{8\pi} B^2`.  In
   effect, we are carrying our magnetic field as :math:`{\bf B}^\prime = {\bf B}/\sqrt{4\pi}`.

   If you wish to express the magnetic field in Gauss, then you will need to multiply the value Castro carries
   by :math:`\sqrt{4\pi}`.


The constrained transport algorithm algebraically ensures that
:math:`\nabla \cdot {\bf B} = 0`.  Throughout the algorithm, the
magnetic fields are face-centered, the electric fields are
edge-centered, and the conserved state is cell-centered (or averages),
unless specified otherwise.  We use the following indexing
conventions:

  * ``U(i,j,k)`` is the cell-centered state, :math:`U_{i,j,k}`

  * ``qleft(i,j,k)`` is the interface value of a state, for example,
    in the x-direction this would be :math:`q_{i-1/2,j,k}`

  * ``Bx(i,j,k)`` is the face-centered x-component of the magnetic field,
    :math:`B_{x,i-1/2,j,k}`

  * ``Ex(i,j,k)`` is the edge-centered x-component of the electric field,
    :math:`E_{x,i,j-1/2,k-1/2}`


Problem Initialization
======================

.. index:: ca_initmag

There is an additional initialization routine for MHD,
``problem_initialize_mhd_data()``,
that is used to initialize the face-centered magnetic field
components.  This is done separately from the main conserved fluid
state.

The conserved fluid state is initialized in ``problem_initialize_state_data()`` just as
with pure hydrodynamics problems. Note that you do not need to include 
the magnetic energy contribution to the total energy density, ``UEDEN``.
After this initialization, the driver handles the addition of the magnetic
contribution.   


Hydrodynamics Update
====================

We use piecewise linear or piecewise parabolic reconstruction with a
characteristic projection and the full 12-Riemann solve corner
transport upwind method.  The HLLD Riemann solver is used.

Within the solver, we use the same indexing into the primitive state
as defined in :ref:`table:primlist`, with the additions of ``QMAGX``,
``QMAGY``, and ``QMAGZ`` for the cell-centered magnetic field
components and ``QPTOT`` for the total pressure (gas + magnetic).

Just like with pure hydrodynamics, the reconstruction type is
controlled by ``castro.ppm_type``, with ``0`` selecting piecewise
linear and ``1`` selecting piecewise parabolic.

For the piecewise linear method, the slope limiting is controlled by
the same ``plm_iorder`` runtime parameter as for hydrodynamics.
Additionally, we have an option ``mhd_limit_characteristic`` that
controls whether you want to do the slope limiting on the
characteristic variables (the default) or the primitive variables.

Electric Update
===============

Coupled to the hydrodynamics update is the electric field update.
Here we update the components of E using the contact upwind scheme
first proposed in :cite:`GS2005`.  The updated electric field then
gives the magnetic field via Faraday's law and the discretization ensures
that :math:`\nabla \cdot {\bf B} = 0`.
.. _ch:gravity:

*******
Gravity
*******


Introduction
============

Gravity Types
--------------------

Castro can incorporate gravity as a constant, monopole approximation,
or a full Poisson solve. To enable gravity in the code, set::

    USE_GRAV = TRUE

in the ``GNUmakefile``, and then turn it on in the inputs file
via ``castro.do_grav`` = 1. If you want to incorporate a point mass
(through ``castro.point_mass``), you must have::

    USE_POINTMASS = TRUE

in the ``GNUmakefile``.

There are currently three options for how gravity is calculated,
controlled by setting ``gravity.gravity_type``. The options are
``ConstantGrav``, ``PoissonGrav``, or ``MonopoleGrav``.
Again, these are only relevant if ``USE_GRAV =
TRUE`` in the ``GNUmakefile`` and ``castro.do_grav`` = 1 in the inputs
file. If both of these are set then the user is required to specify
the gravity type in the inputs file or the program will abort.

.. note:: make sure you have set the ``problem::center[]`` variable
   appropriately for you problem.  This can be done by directly
   setting it in the ``problem_initialize()`` function.


Integration Strategy
--------------------

Castro uses subcycling to integrate levels at different timesteps.
The gravity algorithm needs to respect this to obtain full accuracy.
When self-gravity is computed via a multigrid solve
(``gravity.gravity_type = PoissonGrav``), we correct for this (though
not for other gravity types). At coarse-fine interfaces, the stencil
used in the Laplacian understands the coarse-fine interface and is
different than the stencil used in the interior.

There are two types of solves that we discuss with AMR:

-  *composite solve* : This is a multilevel solve, starting at
   a coarse level (usually level 0) and solving for the potential on
   all levels up to the finest level.

-  *level solve* : This solves for the potential only on
   a particular level. Finer levels are ignored. At coarse-fine
   interfaces, the data from the coarse levels acts as Dirichlet
   boundary conditions for the current-level-solve.

The overall integration strategy is as follows, and is similar to
the discussion in :cite:`castro_I`. Briefly:

-  At the beginning of a simulation, we do a multilevel composite
   solve (if ``gravity.no_composite`` = 0).

   We also do a multilevel composite solve after each regrid.

-  The old-time gravity on the coarse level is defined based on
   this composite solve, but we also do a level solve on the coarse
   level, and use it to determine the difference between the composite
   solve and the level solve, and store that in a MultiFab.

-  After the hydro advance on the coarse level, we do another level
   solve, and use the (level solve - compositive solve) as a lagged
   predictor of how much we need to add back to that level solve to get
   an effective new-time value for phi on the coarse level, and that’s
   what defines the phi used for the new-time gravity

-  Then we do the fine grid timestep(s), each using the same
   strategy

-  At an AMR synchronization step across levels (see Section
   :ref:`sec:amr_synchronization` for a
   description of when these synchronizations occur), if we’re
   choosing to synchronize the gravitational field across levels
   (``gravity.no_sync`` = 0) we then do a solve starting from the coarse
   grid that adjusts for the mismatch between the fine-grid phi and
   the coarse-grid phi, as well as the mismatch between the fine-grid
   density fluxes and the coarse-grid density fluxes, and add the
   resulting sync solve phi to both the coarse and the fine level

   Thus, to within the gravity error tolerance, you get the same final
   result as if you had done a full composite solve at the end of the
   timestep (assuming ``gravity.no_sync`` = 0).

If you do ``gravity.no_composite`` = 1, then you never do a full
multilevel solve, and the gravity on any level is defined only by the
solve on that level. The only time this would be appropriate is if
the fine level(s) cover essentially all of the mass on the grid for
all time.

Controls
--------

-  For the full Poisson solver
   (``gravity.gravity_type`` = ``PoissonGrav``), the behavior
   of the full Poisson solve / multigrid solver is controlled by
   ``gravity.no_sync`` and ``gravity.no_composite``.

-  For isolated boundary conditions, and when
   ``gravity.gravity_type`` = ``PoissonGrav``, the parameters
   ``gravity.max_multipole_order`` and
   ``gravity.direct_sum_bcs`` control the accuracy of
   the Dirichlet boundary conditions. These are described in
   Section `2.3.2 <#sec-poisson-3d-bcs>`__.

-  For ``MonopoleGrav``, in 1D we must have ``coord_sys`` = 2, and in
   2D we must have ``coord_sys`` = 1.

The following parameters apply to gravity
solves:

-  ``gravity.gravity_type`` : how should we calculate gravity?
   Can be ``ConstantGrav``, ``PoissonGrav``, or ``MonopoleGrav``

-  ``gravity.const_grav`` : if ``gravity.gravity_type`` =
   ``ConstantGrav``, set the value of constant gravity (default: 0.0)

-  ``gravity.no_sync`` : ``gravity.gravity_type`` =
   ``PoissonGrav``, do we perform the “sync solve"? (0 or 1; default: 0)

-  ``gravity.no_composite`` : if gravity.gravity_type
   = ``PoissonGrav``, whether to perform a composite solve (0 or 1;
   default: 0)

-  ``gravity.max_solve_level`` : maximum level to solve
   for :math:`\phi` and :math:`\mathbf{g}`; above this level, interpolate from
   below (default: ``MAX_LEV``-1)

-  ``gravity.abs_tol`` : if ``gravity.gravity_type`` = ``PoissonGrav``,
   this is the absolute tolerance for the Poisson solve. You can
   specify a single value for this tolerance (or do nothing, and get a
   reasonable default value), and then the absolute tolerance used by
   the multigrid solve is ``abs_tol`` :math:`\times 4\pi G\,
   \rho_{\text{max}}` where :math:`\rho_{\text{max}}` is the maximum
   value of the density on the domain. On fine levels, this absolute
   tolerance is multiplied by :math:`(\mathtt{ref\_ratio})^2` to account
   for the change in the absolute scale of the Laplacian operator. You
   can also specify an array of values for ``abs_tol``, one for each
   possible level in the simulation, and then the scaling by
   :math:`(\mathtt{ref\_ratio})^2` is not applied.

-  ``gravity.rel_tol`` : if ``gravity.gravity_type`` = ``PoissonGrav``,
   this is the relative tolerance for the Poisson solve. By default it
   is zero. You can specify a single value for this tolerance and it
   will apply on every level, or you can specify an array of values
   for ``rel_tol``, one for each possible level in the
   simulation. This replaces the old parameter ``gravity.ml_tol``.

-  ``gravity.max_multipole_order`` : if ``gravity.gravity_type`` =
   ``PoissonGrav``, this is the max :math:`\ell` value to use for
   multipole BCs (must be :math:`\geq 0`; default: 0)

-  ``gravity.direct_sum_bcs`` : if ``gravity.gravity_type`` =
   ``PoissonGrav``, evaluate BCs using exact sum (0 or 1; default: 0)

-  ``gravity.drdxfac`` : ratio of dr for monopole gravity
   binning to grid resolution

The follow parameters affect the coupling of hydro and gravity:

-  ``castro.do_grav`` : turn on/off gravity

-  ``castro.moving_center`` : do we recompute the center
   used for the multipole gravity solver each step?

-  ``castro.point_mass`` : point mass at the center of the star
   (must be :math:`\geq 0`; default: 0.0)

Note that in the following, ``MAX_LEV`` is a hard-coded parameter
in ``Source/Gravity.cpp`` which is currently set to 15. It
determines how many levels can be tracked by the ``Gravity`` object.

Types of Approximations
=======================

``ConstantGrav``
----------------

Gravity can be defined as constant in direction and magnitude,
defined in the inputs file by::

   gravity.const_grav = -9.8

for example, to set the gravity to have magnitude :math:`9.8` in the
negative :math:`y`-direction if in 2D, negative :math:`z`-direction if in 3-D.
The actual setting is done in Gravity.cpp as::

     grav.setVal(const_grav, AMREX_SPACEDIM-1, 1, ng);

Note that at present we do not fill the gravitational potential
:math:`\phi` in this mode; it will be set to zero.

Note: ``ConstantGrav`` can only be used along a Cartesian direction
(vertical for 2D axisymmetric).

.. _sec-monopole-grav:

``MonopoleGrav``
----------------

``MonopoleGrav`` integrates the mass distribution on the grid
in spherical shells, defining an enclosed mass and uses this
to compute the gravitational potential and acceleration in a
spherically-symmetric fashion.

-  In 1D spherical coordinates we compute

   .. math:: g(r) = -\frac{G M_{\rm enclosed}}{ r^2}

   where :math:`M_{\rm enclosed}` is calculated from the density at
   the time of the call.

   For levels above the coarsest level we define the extent of that
   level’s radial arrays as ranging from the center of the star (:math:`r=0`)
   to the cell at that level farthest away from the origin. If there
   are gaps between fine grids in that range then we interpolate the
   density from a coarser level in order to construct a continuous
   density profile. We note that the location of values in the density
   profile and in the gravitational field exactly match the location of
   data at that level so there is no need to interpolate between points
   when mapping the 1D radial profile of :math:`g` back onto the original
   grid.

-  In 2D or 3D we compute a 1D radial average of density and use
   this to compute gravity as a one-dimensional integral, then
   interpolate the gravity vector back onto the Cartesian grid
   cells. At the coarsest level we define the extent of the 1D arrays
   as ranging from the center of the star to the farthest possible
   point in the grid (plus a few extra cells so that we can fill ghost
   cell values of gravity). At finer levels we first define a single
   box that contains all boxes on that fine level, then we interpolate
   density from coarser levels as needed to fill the value of density
   at every fine cell in that box. The extent of the radial array is
   from the center of the star to the *nearest* cell on one of the
   faces of the single box. This ensures that all cells at that
   maximum radius of the array are contained in this box.

   We then average the density onto a 1D radial array. We note that
   there is a mapping from the Cartesian cells to the radial array and
   back; unlike the 1D case this requires interpolation. We use
   quadratic interpolation with limiting so that the interpolation
   does not create new maxima or minima.

   The default resolution of the radial arrays at a level is the grid
   cell spacing at that level, i.e., :math:`\Delta r = \Delta x`.
   For increased accuracy, one can define ``gravity.drdxfac`` as a number
   greater than :math:`1` (:math:`2` or :math:`4` are recommended) and
   the spacing of the radial array will then satisfy :math:`\Delta x /
   \Delta r =` drdxfac.  Individual Cartesian grid cells are
   subdivided by drdxfac in each coordinate direction for the
   purposing of averaging the density, and the integration that
   creates :math:`g` is done at the finer resolution of the new
   :math:`\Delta r`.

   Note that the center of the star is defined in the subroutine
   ``probinit`` and the radius is computed as the distance from that
   center.

   .. note:: there is an additional correction at the corners in
             ``make_radial_grav`` that accounts for the volume in a shell
             that is not part of the grid.

What about the potential in this case? when does
``make_radial_phi`` come into play?

``PoissonGrav``
---------------

The most general case is a self-induced gravitational field,

.. math:: \mathbf{g}(\mathbf{x},t) = \nabla \phi

where :math:`\phi` is defined by solving

.. math::
   \mathbf{\Delta} \phi = 4 \pi G \rho
   :label: eq:Self Gravity

We only allow ``PoissonGrav`` in 2D or 3D because in 1D, computing
the monopole approximation in spherical coordinates is faster and more
accurate than solving the Poisson equation.

Poisson Boundary Conditions: 2D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In 2D, if boundary conditions are not periodic in both directions, we
use a monopole approximation at the coarsest level. This involves
computing an effective 1D radial density profile (on level = 0 only),
integrating it outwards from the center to get the gravitational
acceleration :math:`\mathbf{g}`, and then integrating :math:`g`
outwards from the center to get :math:`\phi` (using :math:`\phi(0) =
0` as a boundary condition, since no mass is enclosed at :math:`r =
0`). For more details, see Section `2.2 <#sec-monopole-grav>`__.

.. _sec-poisson-3d-bcs:

Poisson Boundary Conditions: 3D
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The following describes methods for doing isolated boundary
conditions. The best reference for Castro’s implementation of this
is :cite:`katz:2016`.

-  **Multipole Expansion**

   In 3D, by default, we use a multipole expansion to estimate the value
   of the boundary conditions. According to, for example, Jackson’s
   *Classical Electrodynamics* (with the corresponding change to
   Poisson’s equation for electric charges and gravitational
   ”charges”), an expansion in spherical harmonics for :math:`\phi` is

   .. math::
      \phi(\mathbf{x}) = -G\sum_{l=0}^{\infty}\sum_{m=-l}^{l} \frac{4\pi}{2l + 1} q_{lm} \frac{Y_{lm}(\theta,\phi)}{r^{l+1}},
      :label: spherical_harmonic_expansion

   The origin of the coordinate system is taken to be the ``center``
   variable, that must be declared and stored in the ``probdata``
   module in your project directory. The validity of the expansion used
   here is based on the assumption that a sphere centered on
   ``center``, of radius approximately equal to the size of half the
   domain, would enclose all of the mass. Furthermore, the lowest order
   terms in the expansion capture further and further departures from
   spherical symmetry. Therefore, it is crucial that ``center`` be
   near the center of mass of the system, for this approach to achieve
   good results.

   The multipole moments :math:`q_{lm}` can be calculated by expanding the
   Green’s function for the Poisson equation as a series of spherical
   harmonics, which yields

   .. math:: q_{lm} = \int Y^*_{lm}(\theta^\prime, \phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3x^\prime. \label{multipole_moments_original}

   Some simplification of :eq:`spherical_harmonic_expansion` can
   be achieved by using the addition theorem for spherical harmonics:

   .. math::

      \begin{aligned}
        &\frac{4\pi}{2l+1} \sum_{m=-l}^{l} Y^*_{lm}(\theta^\prime,\phi^\prime)\, Y_{lm}(\theta, \phi) = P_l(\text{cos}\, \theta) P_l(\text{cos}\, \theta^\prime) \notag \\
        &\ \ + 2 \sum_{m=1}^{l} \frac{(l-m)!}{(l+m)!} P_{l}^{m}(\text{cos}\, \theta)\, P_{l}^{m}(\text{cos}\, \theta^\prime)\, \left[\text{cos}(m\phi)\, \text{cos}(m\phi^\prime) + \text{sin}(m\phi)\, \text{sin}(m\phi^\prime)\right].\end{aligned}

   Here the :math:`P_{l}^{m}` are the associated Legendre polynomials and the
   :math:`P_l` are the Legendre polynomials. After some algebraic
   simplification, the potential outside of the mass distribution can be
   written in the following way:

   .. math:: \phi(\mathbf{x}) \approx -G\sum_{l=0}^{l_{\text{max}}} \left[Q_l^{(0)} \frac{P_l(\text{cos}\, \theta)}{r^{l+1}} + \sum_{m = 1}^{l}\left[ Q_{lm}^{(C)}\, \text{cos}(m\phi) + Q_{lm}^{(S)}\, \text{sin}(m\phi)\right] \frac{P_{l}^{m}(\text{cos}\, \theta)}{r^{l+1}} \right].

   The modified multipole moments are:

   .. math::

      \begin{aligned}
        Q_l^{(0)}   &= \int P_l(\text{cos}\, \theta^\prime)\, {r^{\prime}}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(C)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{cos}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime \\
        Q_{lm}^{(S)} &= 2\frac{(l-m)!}{(l+m)!} \int P_{l}^{m}(\text{cos}\, \theta^\prime)\, \text{sin}(m\phi^\prime)\, {r^\prime}^l \rho(\mathbf{x}^\prime)\, d^3 x^\prime.\end{aligned}

   Our strategy for the multipole boundary conditions, then, is to pick
   some value :math:`l_{\text{max}}` that is of sufficiently high order to
   capture the distribution of mass on the grid, evaluate the discretized
   analog of the modified multipole moments for :math:`0 \leq l \leq
   l_{\text{max}}` and :math:`1 \leq m \leq l`, and then directly compute the
   value of the potential on all of the boundary zones. This is
   ultimately an :math:`\mathcal{O}(N^3)` operation, the same order as the
   monopole approximation, and the wall time required to calculate the
   boundary conditions will depend on the chosen value of
   :math:`l_{\text{max}}`.

   The number of :math:`l` values calculated is controlled by
   ``gravity.max_multipole_order`` in your inputs file. By default, it
   is set to ``0``, which means that a monopole approximation is
   used. There is currently a hard-coded limit of
   :math:`l_{\text{max}} = 50`. This is because the method used to
   generate the Legendre polynomials is not numerically stable for
   arbitrary :math:`l` (because the polynomials get very large, for
   large enough :math:`l`).

-  **Direct Sum**

   Up to truncation error caused by the discretization itself, the
   boundary values for the potential can be computed exactly by a direct
   sum over all cells in the grid. Suppose I consider some ghost cell
   outside of the grid, at location :math:`\mathbf{r}^\prime \equiv (x^\prime,
   y^\prime, z^\prime)`. By the principle of linear superposition as
   applied to the gravitational potential,

   .. math:: \phi(\mathbf{r}^\prime) = \sum_{\text{ijk}} \frac{-G \rho_{\text{ijk}}\, \Delta V_{\text{ijk}}}{\left[(x - x^\prime)^2 + (y - y^\prime)^2 + (z - z^\prime)^2\right]^{1/2}},

   where :math:`x = x(i)`, :math:`y = y(j)` and :math:`z = z(k)` are
   constructed in the usual sense from the zone indices. The sum here
   runs over every cell in the physical domain (that is, the
   calculation is :math:`\mathcal{O}(N^3)` for each boundary
   cell). There are :math:`6N^2` ghost cells needed for the Poisson
   solve (since there are six physical faces of the domain), so the
   total cost of this operation is :math:`\mathcal{O}(N^5)` (this only
   operates on the coarse grid, at present). In practice, we use the
   domain decomposition inherent in the code to implement this solve:
   for the grids living on any MPI task, we create six :math:`N^2`
   arrays representing each of those faces, and then iterate over
   every cell on each of those grids, and compute their respective
   contributions to all of the faces. Then, we do a global reduce to
   add up the contributions from all cells together. Finally, we place
   the boundary condition terms appropriate for each grid onto its
   respective cells.

   This is quite expensive even for reasonable sized domains, so this
   option is recommended only for analysis purposes, to check if the
   other methods are producing accurate results. It can be enabled by
   setting ``gravity.direct_sum_bcs`` = 1 in your inputs file.

Point Mass
----------

Pointmass gravity works with all other forms of gravity, it is not a
separate option. Since the Poisson equation is linear in potential
(and its derivative, the acceleration, is also linear), the point mass
option works by adding the gravitational acceleration of the point
mass onto the acceleration from whatever other gravity type is under
in the simulation.

.. note:: The point mass may have a mass < 0

A useful option is ``point_mass_fix_solution``. If set to 1, then it
takes all zones that are adjacent to the location of the center
variable and keeps their density constant. Any changes in density that
occur after a hydro update in those zones are reset, and the mass
deleted is added to the pointmass. (If there is expansion, and the
density lowers, then the point mass is reduced and the mass is added
back to the grid). This calculation is done in
``pointmass_update()`` in ``Castro_pointmass.cpp``.

GR correction
=============

In the cases of compact objects or very massive stars, the general
relativity (GR) effect starts to play a role [2]_. First, we consider
the hydrostatic equilibrium due to effects of GR then derive
GR-correction term for Newtonian gravity.  The correction term is
applied to the monopole approximation only when ``USE_GR = TRUE`` is
set in the ``GNUmakefile``.

The formulae of GR-correction here are based on
:cite:`grbk1`. For detailed physics, please refer to
:cite:`grbk2`. For describing very strong gravitational
field, we need to use Einstein field equations

.. math::
   R_{ik}-\frac{1}{2}g_{ik}R=\frac{\kappa}{c^{2}}T_{ik} \quad , \quad
   \kappa=\frac{8\pi G}{c^{2}}\quad ,
   :label: field

where :math:`R_{ik}` is the Ricci tensor, :math:`g_{ik}` is the metric
tensor, :math:`R` is the Riemann curvature, :math:`c` is the speed of
light and :math:`G` is gravitational constant. :math:`T_{ik}` is the
energy momentum tensor, which for ideal gas has only the non-vanishing
components :math:`T_{00}` = :math:`\varrho c^2` , :math:`T_{11}` =
:math:`T_{22}` = :math:`T_{33}` = :math:`P` ( contains rest mass and
energy density, :math:`P` is pressure). We are interested in
spherically symmetric mass distribution. Then the line element
:math:`ds` for given spherical coordinate :math:`(r, \vartheta,
\varphi)` has the general form

.. math::

   \label{metric}
     ds^{2} = e^{\nu}c^{2}dt^{2}-e^{\lambda}dr^{2}-r^{2}(d\vartheta^{2}+\sin^{2}
     \vartheta d\varphi) \quad ,

with :math:`\nu = \nu(r)`, :math:`\lambda = \lambda(r)`. Now we can
put the expression of :math:`T_{ik}` and :math:`ds` into :eq:`field`,
then field equations can be reduced to 3 ordinary differential
equations:

.. math::
   \frac{\kappa P}{c^{2}} =
      e^{-\lambda}\left (\frac{\nu^{\prime}}{r}+\frac{1}{r^{2}} \right )-\frac{1}{r^{2}}
      \quad ,
   :label: diff1

.. math::
   \frac{\kappa P}{c^{2}} =
     \frac{1}{2}e^{-\lambda}\left (\nu^{\prime\prime}+\frac{1}{2}{\nu^{\prime}}^{2}+\frac{\nu^
       {\prime}-\lambda^{\prime}}{r}
      -\frac{\nu^{\prime}\lambda^{\prime}}{2} \right ) \quad ,
   :label: diff2

.. math::
   \kappa \varrho =
     e^{-\lambda}\left (\frac{\lambda^{\prime}}{r}-\frac{1}{r^{2}}\right )+\frac{1}{r^{2}} \quad ,
   :label: diff3

where primes means the derivatives with respect to :math:`r`. After
multiplying with :math:`4\pi r^2`, :eq:`diff3` can be
integrated and yields

.. math::

   \label{gmass1}
     \kappa m = 4\pi r (1-e^{-\lambda}) \quad ,

the :math:`m` is called “gravitational mass” inside r defined as

.. math::

   \label{gmass2}
     m = \int_{0}^{r}4\pi r^{2}  \varrho dr\quad .

For the :math:`r = R`, :math:`m` becomes the mass :math:`M` of the
star. :math:`M` contains not only the rest mass but the whole energy
(divided by :math:`c^2`), that includes the internal and gravitational
energy. So the :math:`\varrho = \varrho_0 +U/c^2` contains the whole
energy density :math:`U` and rest-mass density
:math:`\varrho_0`. Differentiation of :eq:`diff1` with respect to
:math:`r` gives :math:`P = P^{\prime}(\lambda,\lambda^{\prime},
\nu,\nu^{\prime},r)`, where
:math:`\lambda,\lambda^{\prime},\nu,\nu^{\prime}` can be eliminated by
:eq:`diff1`, :eq:`diff2`, :eq:`diff3`. Finally we reach
*Tolman-Oppenheinmer-Volkoff (TOV)* equation for hydrostatic
equilibrium in general relativity:

.. math::
   \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \left (1+\frac{P}{\varrho
       c^{2}}\right )\left (1+\frac{4\pi r^3 P}{m c^{2}}\right ) \left (1-\frac{2Gm}{r c^{2}} \right)^{-1} \quad .
   :label: tov

For Newtonian case :math:`c^2 \rightarrow  \infty`, it reverts to usual form

.. math::

   \label{newton}
     \frac{dP}{dr} = -\frac{Gm}{r^{2}}\varrho \quad .

Now we take effective monopole gravity as

.. math::
   \tilde{g} = -\frac{Gm}{r^{2}} (1+\frac{P}{\varrho
     c^{2}})(1+\frac{4\pi r^3 P}{m c^{2}}) (1-\frac{2Gm}{r c^{2}})^{-1}  \quad .
   :label: tov2

For general situations, we neglect the :math:`U/c^2` and potential
energy in m because they are usually much smaller than
:math:`\varrho_0`. Only when :math:`T` reaches :math:`10^{13} K`
(:math:`KT \approx m_{p} c^2`, :math:`m_p` is proton mass) before it
really makes a difference. So :eq:`tov2` can be expressed as

.. math::

   \label{tov3}
     \tilde{g} = -\frac{GM_{\rm enclosed}}{r^{2}} \left (1+\frac{P}{\varrho
       c^{2}} \right )\left (1+\frac{4\pi r^3 P}{M_{\rm enclosed} c^{2}} \right ) \left (1-\frac{2GM_{\rm enclosed}}{r c^{2}} \right )^{-1} \quad ,

where :math:`M_{enclosed}` has the same meaning as with the
``MonopoleGrav`` approximation.

Hydrodynamics Source Terms
==========================

There are several options to incorporate the effects of gravity into
the hydrodynamics system. The main parameter here is
``castro.grav_source_type``.

- ``castro.grav_source_type`` = 1 : we use a standard
  predictor-corrector formalism for updating the momentum and
  energy. Specifically, our first update is equal to :math:`\Delta t
  \times \mathbf{S}^n` , where :math:`\mathbf{S}^n` is the value of
  the source terms at the old-time (which is usually called time-level
  :math:`n`). At the end of the timestep, we do a corrector step where
  we subtract off :math:`\Delta t / 2 \times \mathbf{S}^n` and add on
  :math:`\Delta t / 2 \times \mathbf{S}^{n+1}`, so that at the end of
  the timestep the source term is properly time centered.

- ``castro.grav_source_type`` = 2 : we do something very similar
  to 1. The major difference is that when evaluating the energy source
  term at the new time (which is equal to :math:`\mathbf{u} \cdot
  \mathbf{S}^{n+1}_{\rho \mathbf{u}}`, where the latter is the
  momentum source term evaluated at the new time), we first update the
  momentum, rather than using the value of :math:`\mathbf{u}` before
  entering the gravity source terms. This permits a tighter coupling
  between the momentum and energy update and we have seen that it
  usually results in a more accurate evolution.

- ``castro.grav_source_type`` = 3 : we do the same momentum update as
  the previous two, but for the energy update, we put all of the work
  into updating the kinetic energy alone. In particular, we explicitly
  ensure that :math:`(\rho e)` remains the same, and update
  :math:`(\rho K)` with the work due to gravity, adding the new kinetic
  energy to the old internal energy to determine the final total gas
  energy. The physical motivation is that work should be done on the
  velocity, and should not directly update the temperature—only
  indirectly through things like shocks.

- ``castro.grav_source_type`` = 4 : the energy update is done in a
  “conservative” fashion. The previous methods all evaluate the value
  of the source term at the cell center, but this method evaluates the
  change in energy at cell edges, using the hydrodynamical mass
  fluxes, permitting total energy to be conserved (excluding possible
  losses at open domain boundaries). See
  :cite:`katzthesis` for some more details.

.. [2]
   Note: The GR
   code and text here were contributed by Ken Chen of Univ. of
   Minnesota.
************
Verification
************

Hydrodynamics Test Problems
===========================

Sod’s Problem (and Other Shock Tube Problems)
---------------------------------------------

The Exec/hydro_tests/Sod problem directory sets up a general one-dimensional
shock tube. The left and right primitive-variable states are specified
and the solution evolves until a user-specified end time. For a simple
discontinuity, the exact solution can be found from an exact Riemann
solver. For this problem, the exact solutions were computed with the
exact Riemann solver from Toro :cite:`toro:1997`, Chapter 4.

Sod’s Problem
~~~~~~~~~~~~~

The Sod problem :cite:`sod:1978` is a simple shock tube problem that
exhibits a shock, contact discontinuity, and a rarefaction wave.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = 0 \\
   p_L = 1
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 0.125 \\
   u_R = 0 \\
   p_R = 0.1
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.2` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files setup the
Sod’s problem:


+------------------+-----------------------------------------+
| ``inputs-sod-x`` | Sod’s problem along :math:`x`-direction |
+------------------+-----------------------------------------+
| ``inputs-sod-y`` | Sod’s problem along :math:`y`-direction |
+------------------+-----------------------------------------+
| ``inputs-sod-z`` | Sod’s problem along :math:`z`-direction |
+------------------+-----------------------------------------+

For multi-dimensional runs, the directions transverse to the jump are
kept constant. We use a CFL number of 0.9, an initial timestep shrink
(castro.init_shrink) of 0.1, and the maximum factor by which
the timestep can increase (castro.change_max) of 1.05.

.. _fig:sod:
.. figure:: sod_3d.png
   :alt: sod problem
   :align: center
   :width: 4.75in

   Castro solution for Sod’s problem run in 3-d, with the newest ppm
   limiters, along the :math:`x`, :math:`y`, and :math:`z` axes. A
   coarse grid of 32 zones in the direction of propagation, with 2
   levels of refinement was used. The analytic solution appears as the
   red line.

.. _fig:sod_ppm0:
.. figure:: sod_3d_ppm0.png
   :alt: sod problem
   :align: center
   :width: 4.75in

   Castro solution for Sod’s problem run in 3-d, with the
   piecewise-linear Godunov method with limiters, along the :math:`x`,
   :math:`y`, and :math:`z` axes. A coarse grid of 32 zones in the
   direction of propagation, with 2 levels of refinement was used. The
   analytic solution appears as the red line.

:numref:`fig:sod` shows the Castro solution using the newest PPM limiters
compared to the analytic
solution, showing the density, velocity, pressure, and internal energy.
:numref:`fig:sod_ppm0` is the same as :numref:`fig:sod`,
but with the piecewise-linear Godunov method with limiters,
shown for comparison.

The Verification subdirectory includes the analytic solution for
the Sod problem sod-exact.out, with :math:`\gamma = 1.4`. 1-d slices
can be extracted from the Castro plotfile using the fextract tool
from ``amrex/Tools/Postprocessing/C_Src/``.
The steps to generate this verification plot with Castro are:

#. in ``Exec/hydro_tests/Sod``, build the Castro executable in 3-d

#. run the Sod problem with Castro in the :math:`x`, :math:`y`, and :math:`z` directions::

    ./Castro3d.Linux.Intel.Intel.ex inputs-sod-x
    ./Castro3d.Linux.Intel.Intel.ex inputs-sod-y
    ./Castro3d.Linux.Intel.Intel.ex inputs-sod-z

#. build the fextract tool in ``amrex/Tools/Postprocessing/C_Src/`` .

#. run fextract on the Castro output to generate 1-d slices
   through the output::

    fextract3d.Linux.Intel.exe -d 1 -s sodx.out -p sod_x_plt00034
    fextract3d.Linux.Intel.exe -d 2 -s sody.out -p sod_y_plt00034
    fextract3d.Linux.Intel.exe -d 3 -s sodz.out -p sod_z_plt00034

#. copy the sodx/y/z.out files into the ``Verification/`` directory.

#. in ``Verification/`` run the gnuplot script ``sod_3d.gp`` as::

    gnuplot sod_3d.gp

   This will produce the figure ``sod_3d.eps``.

Double Rarefaction
~~~~~~~~~~~~~~~~~~

The double rarefaction is the “Test 2” problem described by Toro
:cite:`toro:1997`, Chapter 6. In this test, the center of the domain
is evacuated as two rarefaction waves propagate in each direction, outward
from the center. It is difficult to get the internal energy to
behave at the center of the domain because we are creating a vacuum.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = -2 \\
   p_L = 0.4
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 1 \\
   u_R = 2 \\
   p_R = 0.4
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.15` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files setup the
double rarefaction problem:


+-----------------------+-----------------------+
| ``inputs-test2-x``    | Double rarefaction    |
|                       | problem along         |
|                       | :math:`x`-direction   |
+-----------------------+-----------------------+
| ``inputs-test2-y``    | Double rarefaction    |
|                       | problem along         |
|                       | :math:`y`-direction   |
+-----------------------+-----------------------+
| ``inputs-test2-z``    | Double rarefaction    |
|                       | problem along         |
|                       | :math:`z`-direction   |
+-----------------------+-----------------------+


We use a CFL number of 0.8, an initial timestep shrink
(``castro.init_shrink``) of 0.1, and the maximum factor by which the
timestep can increase (``castro.change_max``) of 1.05. The PPM solver with
the new limiters are used.

.. _fig:test2:
.. figure:: test2_3d.png
   :alt: double rarefaction
   :align: center
   :width: 5in

   Castro solution for the double rarefaction problem run in 3-d,
   along the :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid
   of 32 zones in the direction of propagation, with 2 levels of
   refinement was used. The analytic solution appears as the red line.

:numref:`fig:test2` shows the Castro output, run along all 3
coordinate axes in 3-d, compared to the analytic solution.

The comparison to the analytic solution follows the same procedure as
described for the Sod’s problem above. The gnuplot script
``test2_3d.gp`` will generate the figure, from the 1-d slices created by
fextract named test2x.out, test2y.out, and test2z.out.

Strong Shock
~~~~~~~~~~~~

The strong shock test is the “Test 3” problem described by Toro
:cite:`toro:1997`, Chapter 6. In this test, a large pressure jump
at the initial interface creates a very strong rightward moving
shock, followed very closely by a contact discontinuity.
The initial conditions are:

.. math::

   \begin{array}{l}
   \rho_L = 1 \\
   u_L = 0 \\
   p_L = 1000
   \end{array}
   \qquad
   \begin{array}{l}
   \rho_R = 1 \\
   u_R = 0 \\
   p_R = 0.01
   \end{array}

The gamma_law equation of state is used with :math:`\gamma = 1.4`.
The system is evolved until :math:`t = 0.012` s. Setups for 1-, 2-, and 3-d
are provided. The following inputs files and probin files setup the
strong shock problem:


+-----------------------+-----------------------+
| ``inputs-test3-x``    | Strong shock problem  |
|                       | along                 |
|                       | :math:`x`-direction   |
+-----------------------+-----------------------+
| ``inputs-test3-y``    | Strong shock problem  |
|                       | along                 |
|                       | :math:`y`-direction   |
+-----------------------+-----------------------+
| ``inputs-test3-z``    | Strong shock problem  |
|                       | along                 |
|                       | :math:`z`-direction   |
+-----------------------+-----------------------+

We use a CFL number of 0.9, an initial
timestep shrink (``castro.init_shrink``) of 0.1, and the maximum factor by which
the timestep can increase (``castro.change_max``) of 1.05. The PPM
solver with the new limiters are used.

.. _fig:test3:
.. figure:: test3_3d.png
   :alt: strong shock
   :align: center
   :width: 5in

   Castro solution for the strong shock problem run in 3-d, along the
   :math:`x`, :math:`y`, and :math:`z` axes. A coarse grid of 32 zones
   in the direction of propagation, with 2 levels of refinement was
   used. The analytic solution appears as the red line.

:numref:`fig:test3` shows the Castro output, run along all 3
coordinate axes in 3-d, compared to the analytic solution.

The comparison to the analytic solution follows the same procedure as
described for the Sod’s problem above. The gnuplot script
``test3_3d.gp`` will generate the figure, from the 1-d slices created by
fextract named test3x.out, test3y.out, and test3z.out.

Sedov Problem
-------------

The Sedov (or Sedov-Taylor) blast wave is a standard hydrodynamics
test problem. A large amount of energy is placed into a very small
volume, driving a spherical (or cylindrical in 2-d Cartesian
coordinates) blast wave. Analytic solutions were found by Sedov
:cite:`sedov:1959`.

A cylindrical blast wave (e.g. a point explosion in a 2-d plane) can
be modeled in 2-d Cartesian coordinates. A spherical blast wave can
be modeled in 1-d spherical, 2-d axisymmetric (cylindrical :math:`r`-:math:`z`), or 3-d
Cartesian coordinates. This provides a good test on the geometric
factors in the hydrodynamics solver.
We use a publically available code, ``sedov3.f``
:cite:`timmes_sedov_code`, to generate the analytic solutions.

The Castro implementation of the Sedov problem is ``in
Exec/hydro_tests/Sedov/``.  A number of different inputs files
are provided, corresponding to different Sedov/Castro geometries. The
main ones are:


.. _table:sedov_inputs:
.. table:: Sedov problem inputs files

     +---------------------------------+---------------------------------------------+
     | inputs file                     | description                                 |
     +=================================+=============================================+
     | ``inputs.1d.sph``               |  Spherical Sedov explosion modeled          |
     |                                 |  in 1-d spherical coordinates               |
     +---------------------------------+---------------------------------------------+
     | ``inputs.2d.sph_in_cylcoords``  |  Spherical Sedov explosion modeled          |
     |                                 |  in 2-d cylindrical (axisymmetric)          |
     |                                 |  coordinates.                               |
     +---------------------------------+---------------------------------------------+
     | ``inputs.2d.cyl_in_cartcoords`` |  Cylindrical Sedov explosion modeled in     |
     |                                 |  2-d Cartesian coordinates.                 |
     +---------------------------------+---------------------------------------------+
     | ``inputs.3d.sph``               |  Spherical Sedov explosion modeled in       |
     |                                 |  3-d Cartesian coordinates.                 |
     +---------------------------------+---------------------------------------------+

In the Sedov problem, the explosion energy, :math:`E_\mathrm{exp}` (in units
of energy, not energy/mass or energy/volume)
is to be deposited into a single point, in a medium of uniform ambient
density, :math:`\rho_\mathrm{ambient}`, and pressure, :math:`p_\mathrm{ambient}`.
Initializing the problem can be difficult because the small volume is
typically only a cell in extent. This can lead to grid imprinting in
the solution. A standard solution (see for example :cite:`omang:2006`
and the references therein)
is to convert the explosion energy into a pressure contained within a
certain volume, :math:`V_\mathrm{init}`, of radius :math:`r_\mathrm{init}` as

.. math:: p_\mathrm{init} = \frac{(\gamma - 1) E_\mathrm{exp}}{V_\mathrm{init}}

This pressure is then deposited in all of the cells where :math:`r <
r_\mathrm{init}`.

To further minimize any grid effects, we do subsampling
in each zone: each zone is divided it into :math:`N_\mathrm{sub}` subzones in each
coordinate direction, each subzone is initialized independently, and
then the subzones are averaged together (using a volume weighting for
spherical or cylindrical/axisymmetric Castro grids) to determine the
initial state of the full zone.

For these runs, we use :math:`\rho_\mathrm{ambient} = 1`,
:math:`p_\mathrm{ambient} = 10^{-5}`, :math:`E_\mathrm{exp} = 1`, :math:`r_\mathrm{init} = 0.01`,
and :math:`N_\mathrm{sub} = 10`. A base grid with 32 zones in each
coordinate direction plus 3 levels of refinement is used (the finest
mesh would coorespond to 256 zones in a coordinate direction). The
domain runs from 0 to 1 in each coordinate direction.

An analysis routines for the Sedov problem is provided in
``Castro/Diagnostics/Sedov/``.  Typing ``make`` should build it (you
can specify the dimensionality with the ``DIM`` variable in the
build).


A spherical Sedov explosion can be modeled in 1-d spherical, 2-d
cylindrical (axisymmetric), or 3-d Cartesian coordinates, using the
inputs files described in :numref:`table:sedov_inputs`. A 1-d radial
profile can be extracted using the analysis routine. For example, to run and process
the 2-d cylindrical Sedov explosion, one would do:

#. in ``Exec/hydro_tests/Sedov``, build the Castro executable in 2-d

#. run the spherical Sedov problem with Castro in 2-d cylindrical coordinates::

    ./Castro2d.Linux.Intel.Intel.ex inputs.2d.sph_in_cylcoords

#. build the ``sedov_2d_ex`` tool in
   ``Castro/Diagnostics/Sedov``.

#. run the analysis script  on the Castro output to generate 1-d radial
   profiles::

      ./sedov_2d.ex --sphr --yctr 0.5 -s sedov_2d_sph_in_cyl.out \
          -p sedov_2d_sph_in_cyl_plt00246

A similar procedure can be used for the 1-d and 3-d spherical Sedov
explosions (with the output named ``sedov_1d_sph.out`` and
``sedov_3d_sph.out`` respectively). Once this is done, the
``sedov_sph.gp`` gnuplot script can be used to make a plot comparing
the 3 solutions to the analytic solution, ``spherical_sedov.dat``.

:numref:`fig:sedov_sph` shows the comparison of the 3 Castro spherical Sedov explosion simulations to the analytic solution.

.. _fig:sedov_sph:
.. figure:: sedov_sph.png
   :alt: Sedov blast wave
   :align: center
   :width: 5in

   Castro solution for the Sedov blast wave problem run in 1-d
   spherical, 2-d axisymmetric, and 3-d Cartesian coordinates.  Each
   of these geometries produces a spherical Sedov explosion.

Cylindrical Blast Wave
~~~~~~~~~~~~~~~~~~~~~~

.. figure:: sedov_cyl.png
   :alt: Sedov in 2-d
   :align: center
   :width: 5in

   Castro solution for the Sedov blast wave problem run in 2-d
   Cartesian coordinates. This corresponds to a cylindrical Sedov
   explosion.

Rayleigh-Taylor
---------------

2D. Domain size 0.5 by 1.0. 256 by 512 cells, single level
calculation. Periodic in x, solid walls on top and bottom in y.
Gamma law gas with :math:`\gamma=1.4`, no reactions. Zero initial velocity.
Constant :math:`|\gb|=1`. The density profile is essentially :math:`\rho=1` on
bottom, :math:`\rho=2` on top, but with a perturbation. A single-mode
perturbation is constructed as:

.. math:: \tilde y(x) = 0.5 + 0.01 \frac{\cos(4\pi x) + \cos(4\pi(L_x - x))}{2}

We note that the symmetric form of the cosine is done to ensure that
roundoff error does not introduce a left-right asymmetry in the problem.
Without this construction, the R-T instability will lose its symmetry
as it evolves. This then applied to the interface with a tanh profile
to smooth the transition between the high and low density material:

.. math:: \rho(x,y) = 1 + 0.5\left[1+\tanh\left(\frac{y-\tilde y(x)}{0.005}\right)\right]

Hydrostatic pressure with :math:`p=5.0` at bottom of domain, assuming
:math:`\rho=1` on the lower half of the domain, and :math:`\rho=2` on the upper
half and no density perturbation. We run to :math:`t=2.5` with piecewise
linear, old PPM, and new PPM. CFL=0.9. See :numref:`fig:RT`.

.. _fig:RT:
.. figure:: RT_ppm_type.png
   :alt: Rayleigh-Taylor with different PPM types.
   :align: center
   :width: 6.5in

   Rayleigh-Taylor with different PPM types.

Gravity Test Problems
=====================

Radiation Test Problems
=======================

There are two photon radiation solvers in Castro—a gray solver and a
multigroup solver. The gray solver follows the algorithm outlined
in :cite:`howellgreenough:2003`. We use the notation described in that
paper. In particular, the radiation energy equation takes the form
of:

.. math::

   \frac{\partial E_R}{\partial t} =
    \nabla \cdot \left ( \frac{c \lambda(E_R)}{\kappa_R} \nabla E_R \right ) +
    \kappa_P (4 \sigma T^4 - c E_R )

Here, :math:`E_R` is the radiation energy density, :math:`\kappa_R` is the
Roseland-mean opacity, :math:`\kappa_P` is the Planck-mean opaciy, and
:math:`\lambda` is a quantity :math:`\le 1/3` that is subjected to limiting to
keep the radiation field causal. Castro allows for :math:`\kappa_R`
and :math:`\kappa_P` to be set independently as power-laws.

Light Front
-----------

The light front problem tests the ability of the radiation solver to
operate in the free-streaming limit. A radiation front is
estabilished by initializing one end of the computational domain with
a finite radiation field, and zero radiation field everywhere else.
The speed of propagation of the radiation front is keep in check by
the flux-limiters, to prevent it from exceeding :math:`c`.

Diffusion of a Gaussian Pulse
-----------------------------

The diffusion of a Gaussian pulse problem tests the diffusion term in
the radiation energy equation. The radiation energy density is
initialized at time :math:`t = t_0` to a Gaussian distribution:

.. math:: E_R = (E_R)_0 \exp \left \{ - \frac{1}{4 D t_0} |r - r_0|^2 \right \}

As the radiation diffuses, the overall distribution will remain
Gaussian, with the time-dependent solution of:

.. math:: E_R = (E_R)_0 \frac{t_0}{t_0 + t} \exp \left \{ -\frac{1}{4 D (t_0 + t)} |r - r_0|^2 \right \}

Radiation Source Problem
------------------------

The radiation source problem tests the coupling between the radiation
field and the gas energy through the radiation source term. The
problem begins with the radiation field and gas temperature out of
equilibrium. If the gas is too cool, then the radiation field will
heat it. If the gas is too hot, then it will radiate and cool. In
each case, the gas energy and radiation field will evolve until
thermal equilibrium is achieved.

Our implementation of this problem follows that of
:cite:`swestymyra:2009`.

.. figure:: radiating_source.png
   :alt: radiatin source
   :align: center
   :width: 5in

   Castro solution for radiating source test problem. Heating and
   cooling solutions are shown as a function of time, compared to the
   analytic solution. The gray photon solver was used.

Radiating Sphere
----------------

The radiating sphere (RadSphere) is a multigroup radiation
test problem. A hot sphere is centered at the origin in a spherical
geometry. The spectrum from this sphere follows a Planck
distribution. The ambient medium is at a much lower temperature. A
frequency-dependent opacity makes the domain optically thin for high
frequecies and optically thick for low frequency. At long times, the
solution will be a combination of the blackbody radiation from the
ambient medium plus the radiation that propagated from the hot sphere.
An analytic solution exists :cite:`graziani:2008` which gives the
radiation energy as a function of energy group at a specified time and
distance from the radiating sphere.

Our implementation of this problem is in Exec/radiation_tests/RadSphere and
follows that of :cite:`swestymyra:2009`. The routine that computes
the analytic solution is provided as analytic.f90.

.. figure:: radiating_sphere.png
   :alt: radiating sphere
   :width: 5in

   Castro solution for radiating sphere problem, showing the radiation
   energy density as a function of energy group.  This test was run
   with 64 photon energy groups.

Regression Testing
==================

An automated regression test suite for Castro (or any AMReX-based
code) written in Python exists in the AMReX-Codes github organization.

***************
Getting Started
***************

.. note::

   Castro has two source dependencies: `AMReX <https://github.com/AMReX-Codes/amrex>`_, the adaptive mesh
   library, and `StarKiller Microphysics <https://github.com/starkiller-astro/Microphysics>`_, the collection of equations
   of state, reaction networks, and other microphysics.  The
   instructions below describe how to get these dependencies automatically
   with Castro.


The compilation process is managed by AMReX and its build system.  The
general requirements to build Castro are:

 * A C++17 (or later) compiler (e.g. gcc >= 7.0)

 * A Fortran 20xx compiler

 * python (>= 3.6)

 * GNU make (>= 3.82)

GCC is the main compiler suite used by the developers.

For running in parallel, an MPI library is required.  For running on GPUs,
CUDA 11 or later is required.  More information on parallel builds
is given in section :ref:`ch:mpiplusx`.

Downloading the Code
====================

Castro is maintained as a repository on GitHub, and can be obtained
via standard git clone commands. First, make sure that git
is installed on your machine—we recommend version 1.7.x or higher.


#. Clone/fork the Castro repository from the AMReX-Astro GitHub
   organization, using either HTTP access::

       git clone --recursive https://github.com/AMReX-Astro/Castro.git

   or SSH access if you have an SSH key enabled with GitHub::

       git clone --recursive git@github.com:AMReX-Astro/Castro.git

   The ``--recursive`` option to ``git clone`` is used to ensure
   that all of Castro's dependencies are downloaded. Currently this
   requirement is for the AMReX mesh refinement framework, which is
   maintained in the AMReX-Codes organization on GitHub, and the
   Microphysics repository from the starkiller-astro organization.
   AMReX adds the necessary code for the driver code for the simulation,
   while Microphysics adds the equations of state, reaction
   networks, and other microphysics needed to run Castro.

   If you forget to do a recursive clone, you can rectify the
   situation by running the following from the top-level of the Castro
   directory::

       git submodule update --init --recursive

   .. note::

      By default, you will be on the ``main`` branch of the source.
      Development on Castro (and its primary dependencies, AMReX and
      Microphysics) is done in the ``development`` branch, so you
      should work there if you want the latest source::

        git checkout development

      The Castro team runs nightly regression testing on the
      ``development`` branch, so bugs are usually found and fixed
      relatively quickly, but it is generally less stable than staying
      on the ``main`` branch.

#. We recommend setting the ``CASTRO_HOME`` environment
   variable to point to the path name where you have put Castro.
   Add the following to your ``.bashrc``::

       export CASTRO_HOME="/path/to/Castro/"

   (or use the analogous form for a different shell).

#. You can keep the code up to date with::

       git pull --recurse-submodules

   The recommended frequency for doing this is monthly, if you are on the
   stable ``main`` branch of the code; we issue a new release of the code
   at the beginning of each month.

#. *optional, for developers*: If you prefer, you can maintain AMReX and
   Microphysics as standalone repositories rather than as git submodules.
   To do so, you can clone them from GitHub using::

       git clone https://github.com/AMReX-Codes/amrex.git
       git clone https://github.com/starkiller-astro/Microphysics.git

   or via SSH as::

       git clone git@github.com:/AMReX-Codes/amrex.git
       git clone git@github.com:/starkiller-astro/Microphysics.git

   Then, set the ``AMREX_HOME`` environment variable to point to the
   ``amrex/`` directory, and the ``MICROPHYSICS_HOME`` environment
   variable to point to the ``Microphysics/`` directory. Castro will
   look there instead of in its local ``external/`` subdirectory.

Building the Code
=================

In Castro each different problem setup is stored in its own
sub-directory under ``Castro/Exec/``. You build the
Castro executable in the problem sub-directory. Here we’ll
build the Sedov problem:

#. From the directory in which you checked out the Castro git repo,
   type::

       cd Castro/Exec/hydro_tests/Sedov

   This will put you into a directory in which you can run the Sedov
   problem in 1-d, 2-d or 3-d.

#. In ``Sedov/``, edit the ``GNUmakefile``, and set

   * ``DIM = 2``

     This is the dimensionality—here we pick 2-d.

   * ``COMP = gnu``

     This is the set of compilers. GNUu are a good default choice
     (this will use g++ and gfortran). You can also choose ``pgi`` and
     ``intel`` for example.

     If you want to try other compilers than the GNU suite and they
     don’t work, please let us know.

   * ``DEBUG = FALSE``

     This disables debugging checks and results in a more optimized
     executable.

   * ``USE_MPI = FALSE``

     This turns off parallelization via MPI. Set it to ``TRUE`` to build
     with MPI—this requires that you have the MPI library installed on
     your machine. In this case, the build system will need to know
     about your MPI installation. This can be done by editing the
     makefiles in the AMReX tree, but the default fallback is to look
     for the standard MPI wrappers (e.g. ``mpic++`` and ``mpif90``) to do
     the build.

#. Now type ``make``.

   The resulting executable will look something like
   ``Castro2d.gnu.ex``, which means this is a 2-d version
   of the code compiled with ``COMP = gnu``.

More information on the various build options is given in :ref:`ch:buildsystem`.

Running the Code
================

#. Castro takes an input file that overrides the runtime parameter defaults.
   The code is run as::

       ./Castro2d.gnu.ex inputs.2d.cyl_in_cartcoords

   This will run the 2-d cylindrical Sedov problem in Cartesian
   (:math:`x`-:math:`y` coordinates). You can see other possible
   options, which should be clear by the names of the inputs files.

#. You will notice that running the code generates directories that
   look like ``plt00000/``, ``plt00020/``, etc, and ``chk00000/``,
   ``chk00020/``, etc. These are “plotfiles” and “checkpoint”
   files. The plotfiles are used for visualization, the checkpoint
   files are used for restarting the code.

Visualization of the Results
============================

There are several options for visualizing the data. The popular
packages yt and VisIt both support the AMReX file format
natively [1]_. The standard tool used within the AMReX-community is
Amrvis, which we demonstrate here. Amrvis is available on github.


.. _sec:gettingstartedyt:

yt
^^

yt is the primary visualization and analysis tool used by the
developers.  Install yt following their instructions: `Getting yt
<https://yt-project.org/#getyt>`_ .

You should be able to read in your plotfiles using ``yt.load()`` and
do any of the plots described in the `yt Cookbook
<https://yt-project.org/doc/cookbook/index.html>`_ .

Here we do a sample visualization and analysis of the
plotfiles generated.  This section was generated from a
Jupyter notebook which can be found in
``Docs/source/yt_example.ipynb`` in the Castro repo.

.. include:: yt_example.rst


Amrvis
^^^^^^

Amrvis is a tool developed at LBNL to visualize AMReX data.  It
provides a simple GUI that allows you to quickly visualize slices and
the grid structure.

#. Get Amrvis::

       git clone https://github.com/AMReX-Codes/Amrvis

   Then cd into ``Amrvis/``, edit the ``GNUmakefile`` there
   to set ``DIM = 2``, and again set ``COMP`` to compilers that
   you have. Leave ``DEBUG = FALSE``.

   Type ``make`` to build, resulting in an executable that
   looks like ``amrvis2d...ex``.

   If you want to build amrvis with ``DIM = 3``, you must first
   download and build volpack::

       git clone https://ccse.lbl.gov/pub/Downloads/volpack.git

   Then cd into ``volpack/`` and type ``make``.

   Note: Amrvis requires the OSF/Motif libraries and headers. If you
   don’t have these you will need to install the development version
   of motif through your package manager.  On most Linux
   distributions, the motif library is provided by the openmotif
   package, and its header files (like ``Xm.h``) are provided by
   openmotif-devel. If those packages are not installed, then use the
   package management tool to install them, which varies from
   distribution to distribution, but is straightforward.  ``lesstif``
   gives some functionality and will allow you to build the Amrvis
   executable, but Amrvis may not run properly.

   You may then want to create an alias to amrvis2d, for example::

       alias amrvis2d=/tmp/Amrvis/amrvis2d...ex

   where ``/tmp/Amrvis/amrvis2d...ex`` is the full path and name of
   the Amrvis executable.

#. Configure Amrvis:

   Copy the ``amrvis.defaults`` file to your home directory (you can
   rename it to ``.amrvis.defaults`` if you wish). Then edit the
   file, and change the palette line to point to the full
   path/filename of the ``Palette`` file that comes with Amrvis.

#. Visualize:

   Return to the ``Castro/Exec/hydro_tests/Sedov`` directory. You should
   have a number of output files, including some in the form ``pltXXXXX``,
   where XXXXX is a number corresponding to the timestep the file
   was output.

   ``amrvis2d filename`` to see a single plotfile, or ``amrvis2d -a
   plt*``, which will animate the sequence of plotfiles.

   Try playing around with this—you can change which variable you are
   looking at, select a region and click “Dataset” (under View) in
   order to look at the actual numbers, etc. You can also export the
   pictures in several different formats under "File/Export".

   Some users have found that Amrvis does not work properly under X
   with the proprietary Nvidia graphics driver. A fix for this is
   provided in the FAQ (§ :ref:`ch:faq:vis`)—this is due
   to the default behavior of the DAC in mappuing colors.


.. [1]
   Each of these will recognize it as the
   BoxLib format.
Distributed Problem Setups
==========================

There are a number of standard problem setups that come with Castro.
These can be used as a starting point toward writing your own setup.
We organize these into subdirectories by broad type (radiation, hydro,
gravity, etc.): The standard categories and *some* of the included
problems are:

* ``gravity_tests``:

   * ``DustCollapse``: A pressureless cloud collapse that is a
     standard test problem for gravity. An analytic solution that
     describes the radius of the sphere as a function of time is found
     in Colgate and White :cite:`colgwhite`. This problem is also
     found in the FLASH User’s Guide.

   * ``evrard_collapse``: This is the collapse of an isothermal
     spherical gas cloud.  This problem was originally discussed in
     :cite:`Evrard1988`.
     This implementation of the test comes from section 9.1 of
     :cite:`springel:2010`.


   * ``hydrostatic_adjust``: Model a 1-d stellar atmosphere (plane-parallel or
     spherical/self-gravitating) and dump energy in via an analytic
     heat source and watch the atmosphere’s hydrostatic state adjust
     in response. This is the counterpart to the Maestro
     ``test_basestate`` unit test.

   * ``hse_convergence``: This is meant to be a simple 1-d test for assessing the convergence of
     hydro + gravity in maintaining HSE.  Convergence can be measured either
     by looking at the max :math:`|U|` in the plotfiles.

   * ``hydrostatic_adjust``: This is a problem that explores the
     change in a hydrostatic structure due to heating.  This was used
     originally in :cite:`maestro:II`.

   * ``StarGrav``: This problem sets up a single spherical star in
     hydrostatic equilibrium and is used to assess the ability to
     maintain HSE.

   * ``uniform_cube_sphere``: This is used to compute the
     gravitational potential of a perfect cube, for which there is an
     analytic solution.  It tests our isolated boundary conditions.
     This was demonstrated in :cite:`katz:2016`.
 

* ``hydro_tests``:

   * ``acoustic_pulse``: The acoustic pulse problem from
     :cite:`mccorquodalecolella` used to measure convergence of pure
     hydrodynamics problems (as used for Castro in
     :cite:`castro-sdc`).

   * ``acoustic_pulse_general``: a general equation of state version
     of ``acoustic_pulse`` used for measuring convergence in
     :cite:`castro-sdc`.

   * ``double_bubble``: Initialize 1 or 2 bubbles in a stratified
     atmosphere (isothermal or isentropic) and allow for the bubbles
     to have the same or a different :math:`\gamma` from one another /
     the background atmosphere.  This uses the multigamma EOS.
     An analogous problem is implemented in Maestro.

   * ``KH``: A Kelvin-Helmholtz shear instability problem.

   * ``oddeven``: A grid-aligned shock hitting a very small density
     perturbation.  This demonstrates the odd-even decoupling problem
     discussed in :cite:`quirk1997`. This setup serves to test the
     castro.hybrid_riemann option to hydrodynamics.

   * ``RT``: A single-model Rayleigh-Taylor instability problem.

   * ``Sedov``: The standard Sedov-Taylor blast wave problem. This
     setup was used in the first Castro paper :cite:`castro_I`.

   * ``Sod``: A one-dimensional shock tube setup, including the
     classic Sod problem. This setup was used in the original Castro
     paper.

   * ``Sod_stellar``: A version of the Sod shock tube for the general
     stellar equation of state. This setup and the included inputs
     files was used in :cite:`zingalekatz`.

   * ``toy_convect``: A simple nova-like convection problem with an
     external heating source. This problem shows how to use the model
     parser to initialize a 1-d atmosphere on the Castro grid,
     incorporate a custom tagging routine, sponge the fluid above the
     atmosphere, and write a custom diagnostics routine.
     A MAESTROeX version of this problem setup also exists.

* ``mhd_tests``:

   * ``Alfven``: a linearized MHD wave test problem from :cite:`crockett:2005` and :cite:`miniati_martin`.

   * ``BrioWu``: the Brio Wu shock tube problem as described in :cite:`briowu`.  This is a standard
     test problem used in many MHD code papers (e.g. :cite:`athena`).

   * ``DaiWoodward``: a shock tube problem described in :cite:`Dai_1998`

   * ``FastRarefaction``: a shock tube problem dominated by kinetic energy, as described in :cite:`miniati_martin`

   * ``MagnetosonicWaves``: the fast and slow magnetosonic wave problem from :cite:`crockett:2005`

   * ``OrszagTang``: a two-dimensional magnetized vortex problem, following :cite:`athena`

   * ``RT``: a magnetized Rayleigh-Taylor instability problem

   * ``species``: a simple test problem to ensure that species are accurately advected.


* ``radiation_tests``:

   * ``Rad2Tshock``: This sets up a radiating shock that can be
     compared to a semi-analytic solution described in :cite:`lowrieedwards`.
   
   * ``RadFront``: This is the optically-thin streaming of a radiation front problem
     demonstrated originally in Castro in :cite:`CastroII`.

   * ``RadShestakovBolstad``: This is a linear multigroup diffusion test problem first described
     by :cite:`SHESTAKOV2005` and demonstrated in Castro in :cite:`CastroIII`.

   * ``RadSourceTest``: Test the implementation of the source terms in the gray radiation
     solver.  This does the "relaxation to thermal equilibrium" test as
     described in :cite:`swestymyra:2009`  (originally described in :cite:`turnerstone2001`).

   * ``RadSphere``: This is a multigroup radiating sphere test problem with an analytic solution,
     described in :cite:`graziani:2008` and :cite:`swestymyra:2009` and shown in Castro in :cite:`CastroIII`.
 
   * ``RadSuOlson``: This is a non-equlibrium Marshak wave test described in :cite:`suolson:1996` and shown
     in Castro in :cite:`CastroII`.

   * ``RadSuOlsonMG``: This is a multigroup version of ``RadSuOlson`` described in :cite:`suolson:1999`
     and shown in Castro in :cite:`CastroIII`.

   * ``RadThermalWave``: A thermal wave test adapted from :cite:`howellgreenough:2003` and shown in Castro
     in :cite:`CastroII`.
 
* ``reacting_tests``:

   * ``bubble_convergence``: a reacting bubble problem designed for measuring the convergence of
     the reactive hydro algorithms in Castro.  This was used in :cite:`castro-sdc`.

   * ``reacting_bubble``: A reacting bubble in a stratified white
     dwarf atmosphere. This problem was featured in the
     Maestro reaction paper :cite:`maestro:III`.

   * ``reacting_convergence``: a simple reacting hydrodynamics problem for measuring convergence,
     used in :cite:`castro-sdc`.

* ``science``:

  The problems in the science directory are science problems that have
  appeared in papers (or will shortly).  Many of these are being actively used and are shared
  here for reproducibility.

   * ``Detonation``: this sets up a 1-d detonation that propagates through the domain.

   * ``flame``: this sets up a 1-d deflagration that propagates through the domain.  This setup
     was used for the testing in :cite:`eiden:2020`.

   * ``flame_wave``: this is a model of a flame propagating across a neutron star as a model for
     an X-ray burst.  This was presented in :cite:`eiden:2020`.

   * ``nova``: this models convection at the base of an accreted layer
     on a white dwarf as a model of a nova.

   * ``planet``: this is the problem setup from :cite:`ryu:2018` that models shear and turbulence in a
     hot Jupiter atmosphere.

   * ``subchandra``: a model of sub-Chandra Type Ia supernova that initializes a hot spot in a helium
     layer on a low mass carbon-oxygen white dwarf.

   * ``wdmerger``: a problem setup for modeling white dwarf mergers.  This was used in :cite:`katz:2016`.

   * ``xrb_mixed``: a compressible version of the X-ray burst convection problem from :cite:`zingale:2015`.

* ``unit_tests``:

   * ``diffusion_test``: a test of thermal diffusion (without hydro).  This was used to demonstrate convergence
     in both :cite:`castro-sdc` and :cite:`eiden:2020`.

   * ``particles_test``: a test of passive particles.

.. _sponge_section:

******
Sponge
******

Castro uses a sponge source term to damp velocities in regions where
we are not interested in the dynamics, to prevent them from
dominating the timestep constraint.  Often these are buffer regions
between the domain of interest and the boundary conditions.

The sponge parameters are set in the inputs file. The timescale of the
damping is set through ``castro.sponge_timescale``, while factors such as
the radius/density/pressure at which the sponge starts to begin being applied
are described below.

The sponge value, :math:`f_\mathrm{sponge}` is computed as described below
and then the sponge factor, :math:`f` is computed as:

.. math::

   f = - \left [ 1 - \frac{1}{1 + \alpha f_\mathrm{sponge}}\right ]

for an implicit update or

.. math::

   f = -\alpha f_\mathrm{sponge}

.. index:: sponge_implicit, sponge_timescale

for an explicit update.  This choice is controlled by
``sponge_implicit``.  Here, :math:`\alpha` is constructed from the
``sponge_timescale`` parameter, :math:`t_\mathrm{sponge}` as:

.. math::

   \alpha = \frac{\Delta t}{t_\mathrm{sponge}}

The sponge source is then added to the momentum and total energy equations.

The general sponge parameters are:

       ==========================     ========================
         variable                       runtime parameter
       ==========================     ========================
       :math:`f_\mathrm{lower}`       ``sponge_lower_factor``
       :math:`f_\mathrm{upper}`       ``sponge_upper_factor``
       ==========================     ========================

There are three sponges, each controlled by different runtime parameters:

  * **radial sponge** : The radial sponge operates beyond some radius
    (define in terms of the ``center(:)`` position for the problem).
    It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & r < r_\mathrm{lower} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (r - r_\mathrm{lower})}{\Delta r} \right ) \right ]  & r_\mathrm{lower} \le r < r_\mathrm{upper} \\
                     f_\mathrm{upper} & r \ge r_\mathrm{upper} 
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ==========================     ========================
         variable                       runtime parameter
       ==========================     ========================
       :math:`r_\mathrm{lower}`       ``sponge_lower_radius``
       :math:`r_\mathrm{upper}`       ``sponge_upper_radius``
       ==========================     ========================

    and :math:`\Delta r = r_\mathrm{upper} - r_\mathrm{lower}` .


  * **density sponge** : The density sponge turns on based on density
    thresholds.  At high densities, we usually care about the
    dynamics, so we will want the sponge off, but as the density
    lowers, we gradually turn on the sponge.  It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & \rho > \rho_\mathrm{upper} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (\rho - \rho_\mathrm{upper})}{\Delta \rho} \right ) \right ]  & \rho_\mathrm{upper} \ge \rho > \rho_\mathrm{lower} \\
                     f_\mathrm{upper} & \rho < \rho_\mathrm{lower} 
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ============================     ==========================
         variable                          runtime parameter
       ============================     ==========================
       :math:`\rho_\mathrm{lower}`       ``sponge_lower_density``
       :math:`\rho_\mathrm{upper}`       ``sponge_upper_density``
       ============================     ==========================

    and :math:`\Delta \rho = \rho_\mathrm{upper} - \rho_\mathrm{lower}` .


  * **pressure sponge** : The pressure sponge is just like the density sponge,
    except it is keyed on pressure.  It takes the form:

    .. math::

       f_\mathrm{sponge} = \left \{
             \begin{array}{cc}
                     f_\mathrm{lower}   & p > p_\mathrm{upper} \\
                     f_\mathrm{lower} + \frac{f_\mathrm{upper} - f_\mathrm{lower}}{2}
                          \left [ 1 - \cos \left ( \frac{\pi (p - p_\mathrm{upper})}{\Delta p} \right ) \right ]  & p_\mathrm{upper} \ge p \ge p_\mathrm{lower} \\
                     f_\mathrm{upper} & \rho < \rho_\mathrm{lower} 
             \end{array} \right .


    The parameters controlling the various quantities here are:

       ============================     ==========================
         variable                          runtime parameter
       ============================     ==========================
       :math:`p_\mathrm{lower}`         ``sponge_lower_pressure``
       :math:`p_\mathrm{upper}`         ``sponge_upper_pressure``
       ============================     ==========================

    and :math:`\Delta p = p_\mathrm{upper} - p_\mathrm{lower}` .

.. _ch:amr:

************************
Adaptive Mesh Refinement
************************

Our approach to adaptive refinement in Castro uses a nested hierarchy
of logically-rectangular grids with simultaneous refinement of the
grids in both space and time. The integration algorithm on the grid
hierarchy is a recursive procedure in which coarse grids are advanced
in time, fine grids are advanced multiple steps to reach the same time
as the coarse grids and the data at different levels are then
synchronized.

During the regridding step, increasingly finer grids
are recursively embedded in coarse grids until the solution is
sufficiently resolved. An error estimation procedure based on
user-specified criteria (described in Section `1 <#sec:tagging>`__)
evaluates where additional refinement is needed
and grid generation procedures dynamically create or
remove rectangular fine grid patches as resolution requirements change.

A good introduction to the style of AMR used here is in Lecture 1
of the Adaptive Mesh Refinement Short Course at
https://ccse.lbl.gov/people/jbb/shortcourse/lecture1.pdf.

.. _sec:tagging:

Tagging for Refinement
======================

Castro determines what zones should be tagged for refinement at the
next regridding step by using a set of built-in routines that test on
quantities such as the density and pressure and determining whether
the quantities themselves or their gradients pass a user-specified
threshold. This may then be extended if ``amr.n_error_buf`` :math:`> 0`
to a certain number of zones beyond these tagged zones. This section
describes the process by which zones are tagged, and describes how to
add customized tagging criteria.

We also provide a mechanism for defining a limited set of refinement
schemes from the inputs file; for example,

::

   amr.refinement_indicators = dens temp

   amr.refine.dens.max_level = 1
   amr.refine.dens.value_greater = 2.0
   amr.refine.dens.field_name = density

   amr.refine.temp.max_level = 2
   amr.refine.temp.value_less = 1.0
   amr.refine.temp.field_name = Temp

``amr.refinement_indicators`` is a list of user-defined names for refinement
schemes. For each defined name, ``amr.refine.<name>`` accepts predefined fields
describing when to tag. These are:

* ``max_level`` : maximum level to refine to
* ``start_time`` : when to start tagging
* ``end_time`` : when to stop tagging
* ``value_greater`` : value above which we refine
*  ``value_less`` : value below which to refine
* ``gradient`` : absolute value of the difference between adjacent cells above which we refine
* ``relative_gradient`` : relative value of the difference between adjacent cells above which we refine
* ``field_name`` : name of the string defining the field in the code

If a refinement indicator is added, either
``value_greater``, ``value_less``, or ``gradient`` must be provided.

.. note::

   Zones adjacent to a physical boundary cannot be tagged for refinement when
   using the Poisson gravity solver. If your tagging criteria are met in these
   zones, they will be ignored.

.. index:: problem_tagging.H

We provide also the ability for the user to define their own tagging criteria.
This is done through the C++ function ``problem_tagging``
in the file ``problem_tagging.H``. This function is provided the entire
state (including density, temperature, velocity, etc.) and the array
of tagging status for every zone.


.. _sec:amr_synchronization:

Synchronization Algorithm
=========================

Here we present the AMR algorithm for the compressible equations with
self-gravity. The gravity component of the algorithm is closely
related to (but not identical to) that in Miniati and Colella, JCP,
2007. The content here is largely based on the content in the original
Castro paper (:cite:`castro_I`). The most significant difference is the
addition of a different strategy for when to employ the synchronization;
but regardless of whether the original or new strategy is used, the fundamental
synchronization step is identical.

.. _sec:synchronization_methodology:

Synchronization Methodology
---------------------------

Over a coarse grid time step we collect flux register information for
the hyperbolic part of the synchronization:

.. math:: \delta\Fb = -\Delta t_c A^c F^c + \sum \Delta t_f A^f F^f

Analogously, at the end of a coarse grid time step we store the
mismatch in normal gradients of :math:`\phi` at the coarse-fine interface:

.. math::

   \delta F_\phi =  - A^c \frac{\partial \phi^c}{\partial n}
   + \sum A^f \frac{\partial \phi^f}{\partial n}

We want the composite :math:`\phi^{c-f}` to satisfy the multilevel
version of (:eq:`eq:Self Gravity`) at the synchronization time, just
as we want the coarse and fine fluxes at that time to match. So the goal
is to synchronize :math:`\phi` across levels at that time and then zero out
this mismatch register.

At the end of a coarse grid time step we can define
:math:`{\overline{\Ub}}^{c-f}` and :math:`\overline{\phi}^{c-f}` as the composite
of the data from coarse and fine grids as a provisional solution at
time :math:`n+1`. (Assume :math:`\overline{\Ub}` has been averaged down so that
the data on coarse cells underlying fine cells is the average of the
fine cell data above it.)

The synchronization consists of two parts:

-  Step 1: Hyperbolic reflux

   In the hyperbolic reflux step, we update the conserved variables with
   the flux synchronization and adjust the gravitational terms to reflect
   the changes in :math:`\rho` and :math:`\ub`.

   .. math:: {\Ub}^{c, \star} = \overline{\Ub}^{c} + \frac{\delta\Fb}{V},

   where :math:`V` is the volume of the cell and the correction from
   :math:`\delta\Fb` is supported only on coarse cells adjacent to fine grids.

   Note: this can be enabled/disabled via castro.do_reflux. Generally,
   it should be enabled (1).

   Also note that for axisymmetric or 1D spherical coordinates, the
   reflux of the pressure gradient is different, since it cannot be
   expressed as a divergence in those geometries. We use a separate
   flux register in the hydro code to store the pressure term in these
   cases.

-  Step 2: Gravitational synchronization

   In this step we correct for the mismatch in normal derivative in
   :math:`\phi^{c-f}` at the coarse-fine interface, as well as accounting for
   the changes in source terms for :math:`(\rho \ub)` and :math:`(\rho E)` due to the
   change in :math:`\rho.`

   On the coarse grid only, we define

   .. math:: (\delta \rho)^{c} =  \rho^{c, \star} - {\overline{\rho}}^{c}  .

   We then form the composite residual, which is composed of two
   contributions. The first is the degree to which the current :math:`\overline{\phi}^{c-f}` does not satisfy the original equation on a
   composite grid (since we have solved for :math:`\overline{\phi}^{c-f}`
   separately on the coarse and fine levels). The second is the response
   of :math:`\phi` to the change in :math:`\rho.` We define

   .. math::

      R \equiv  4 \pi G \rho^{\star,c-f} - \Delta^{c-f} \; \overline{\phi}^{c-f}
      = - 4 \pi G (\delta \rho)^c - (\nabla \cdot \delta F_\phi ) |_c   .

   Then we solve

   .. math::

      \Delta^{c-f} \; \delta \phi^{c-f} = R
      \label{eq:gravsync}

   as a two level solve at the coarse and fine levels.
   We define the update to gravity,

   .. math:: \delta {\bf g}^{c-f} = \nabla (\delta \phi^{c-f})  .

   Finally, we need to

   -  add :math:`\delta \phi^{c-f}` directly to
      to :math:`\phi^{c}` and :math:`\phi^{f}` and interpolate :math:`\delta \phi^{c-f}` to any finer
      levels and add to the current :math:`\phi` at those levels.

   -  if level :math:`c` is not the coarsest level in the calculation, then we must transmit the
      effect of this change in :math:`\phi` to the coarser levels by updating the flux register between
      level :math:`c` and the next coarser level, :math:`cc.` In particular, we set

      .. math::

         \delta {F_\phi}^{cc-c} = \delta F_\phi^{cc-c}
         + \sum A^c \frac{\partial (\delta \phi)^{c-f}}{\partial n}  .

   The gravity synchronization algorithm can be disabled with
   gravity.no_sync = 1. This should be done with care. Generally,
   it is okay only if he refluxing happens in regions of low density that
   don’t affect the gravity substantially.

.. _sec:synchronization_sources:

Source Terms
------------

After a synchronization has been applied, the state on the coarse grid
has changed, due to the change in fluxes at the coarse-fine boundary as
well as the change in the gravitational field. This poses a problem
regarding the source terms, all of which generally rely either on the
state itself, or on the global variables affected by the synchronization
such as the gravitational field. The new-time sources constructed on the
coarse grid all depended on what the state was after the coarse-grid
hydrodynamic update, but the synchronization and associated flux
correction step retroactively changed that hydrodynamic update. So one
can imagine that in a perfect world, we would have calculated the
hydrodynamic update first, including the coarse-fine mismatch
corrections, and only then computed the source terms at the new time.
Indeed, an algorithm that did not subcycle, but marched every zone along
at the same timestep, could do so – and some codes, like FLASH,
actually do this, where no new-time source terms are computed on any
level until the hydrodynamic update has been fully completed and the
coarse-fine mismatches corrected. But in Castro we cannot do this; in
general we assume the ability to subcycle, so the architecture is set up
to always calculate the new-time source terms on a given level
immediately after the hydrodynamic update on that level. Hence on the
coarse level we calculate the new-time source terms before any fine grid
timesteps occur.

One way to fix this, as suggested by Miniati and Colella for the case
of gravity, is to explicitly compute what the difference in the source
term is as a result of any flux corrections across coarse-fine
boundaries. They work out the form of this update for the case of a
cell-centered gravitational force, which has contributions from both
the density advected across the coarse-fine interfaces
(i.e. :math:`\delta \rho \mathbf{g}`, where :math:`\delta \rho` is the density
change due to the coarse-fine synchronization on the coarse rid), as
well as the global change in the gravitational field due to the
collective mass motion (see Miniati and Colella for the explicit form
of the source term). This has a couple of severe limitations. First,
it means that when the form of the source term is changed, the form of
the corrector term is changed too. For example, it is less easy to
write down the form of this corrector term for the flux-based
gravitational energy source term that is now standard in Castro.
Second, gravity is a relatively easy case due to its linearity in the
density and the gravitational acceleration; other source terms
representing more complicated physics might not have an easily
expressible representation in terms of the reflux contribution. For
example, for a general nuclear reaction network (that does not have an
analytic solution), it is not possible to write down an analytic
expression for the nuclear reactions that occur because of
:math:`\delta \rho`.

Instead we choose a more general approach. On the coarse level, we save
the new-time source terms that were applied until the end of the fine
timesteps. We also save the fine level new-time source terms. Then, when
we do the AMR synchronization after a fine timestep, we first subtract
the previously applied new-time source terms to both the coarse and the
fine level, then do the flux correction and associated gravitational
sync solve, and then re-compute the new-time source terms on both the
coarse and the fine level [1]_. In this way, we get almost
the ideal behavior – if we aren’t subcycling, then we get essentially
the same state at the end of the fine timestep as we would in a code
that explicitly had no subcycling. The cost is re-computing the new-time
source terms that second time on each level. For most common source
terms such as gravity, this is not a serious problem – the cost of
re-computing :math:`\rho \mathbf{g}` (for example, once you already know
:math:`\mathbf{g}`) is negligible compared to the cost of actually computing
:math:`\mathbf{g}` itself (say, for self-gravity). If you believe that the
error in not recomputing the source terms is sufficiently low, or the
computational cost of computing them too high, you can disable this
behavior [2]_ using the
code parameter castro.update_sources_after_reflux.

Note that at present nuclear reactions are not enabled as part of this
scheme, and at present are not automatically updated after an AMR
synchronization. This will be amended in a future release of Castro.

.. _sec:synchronization_timing:

Synchronization Timing
----------------------

The goal of the synchronization step is for the coarse and fine grid to
match at the end of a coarse timesteps, after all subcycled fine grid
timesteps have been completed and the two levels have reached the same
simulation time. If subcycling is disabled, so that the coarse and fine
grid take the same timestep, then this is sufficient. However, in the
general subcycling case, the situation is more complicated. Consider the
discussion about source terms in `2.2 <#sec:synchronization_sources>`__. If
we have a coarse level and one fine level with a refinement ratio of
two, then for normal subcycling the fine grid takes two timesteps for
every one timestep taken by the coarse level. The strategy advocated by
the original Castro paper (and Miniati and Colella) is to only do the
AMR synchronization at the actual synchronization time between coarse
and fine levels, that is, at the end of the second fine timestep.
Consequently, we actually only update the source terms after that second
fine timestep. Thus note that on the fine grid, only the *new-time*
source terms in the *second* fine timestep are updated. But a
moment’s thought should reveal a limitation of this. The first fine grid
timestep was also responsible for modifying the fluxes on the coarse
grid, but the algorithm as presented above didn’t take full account of
this information. So, the gravitational field at the old time in
the second fine timestep is actually missing information that would have
been present if we had updated the coarse grid already. Is there a way
to use this information? For the assumptions we make in Castro, the
answer is actually yes. If we apply the effect of the synchronization
not at the synchronization time but at the end of every fine
timestep, then every fine timestep always has the most up-to-date
information possible about the state of the gravitational field. Now, of
course, in fine timesteps before the last one, we have not actually
reached the synchronization time. But we already know at the end of the
first fine timestep what the synchronization correction will be from
that fine timestep: it will be equal to 1/2 of the coarse contribution
to the flux register and the normal contribution to the flux register
for just that timestep. This is true because in Castro, we assume that
the fluxes provided by the hydrodynamic solver are piecewise-constant
over the timestep, which is all that is needed to be second-order
accurate in time if the fluxes are time centered [3]_. So it is fair to say
that halfway through the coarse timestep, half of the coarse flux has
been advected, and we can mathematically split the flux register into
two contributions that have equal weighting from the coarse flux. (In
general, of course, the coarse flux contribution at each fine timestep
is weighted by :math:`1/R` where :math:`R` is the refinement ratio between the
coarse and fine levels.) So, there is nothing preventing us from
updating the coarse solution at the synchronization time :math:`t^{n+1}_c`
after this first fine timestep; we already know at that point how the
coarse solution will change, so why not use that information? We can
then update the gravitational potential at :math:`t^{n+1/2}_c` that is used to
construct the boundary conditions for the gravitational potential solve
on the fine grid at the beginning of the second fine timestep.

In practice, this just means calling the synchronization routine
described in `2.1 <#sec:synchronization_methodology>`__, with the only
modification being that the flux register contribution from the coarse
grid is appropriately weighted by the fine grid timestep instead of
the coarse grid timestep, and we only include the current fine step:

.. math:: \delta\Fb = -\Delta t_f A^c F^c + \Delta t_f A^f F^f

The form of the :math:`\phi` flux register remains unchanged, because the
intent of the gravity sync solve is to simply instantaneously correct
the mismatch between the fine and coarse grid. The only difference,
then, between the old strategy and this new method is that we call the
synchronization at the end of every fine timestep instead of only the
last subcycled one, and we change the weighting appropriately. This
new method is more expensive as currently implemented because we have
to do :math:`R` gravitational sync solves, refluxes, and source term
recalculations instead of only one. However, it results in maximal
possible accuracy, especially in cases where significant amounts of
material are crossing refinement boundaries. The reflux strategy is
controlled by the parameter castro.reflux_strategy. At present
the old method is still the default.

Note that one does not need to be using self-gravity for this to be
beneficial. Even in pure hydrodynamics this can matter. If a regrid
occurs on the fine level, new zones on the boundaries of the current
fine level are filled by interpolation from the coarse level. In the
old method, that interpolation is not using the most up-to-date data
that accounts for the synchronization.

For multiple levels of refinement, the scheme extends naturally. In
the old method, we always call the synchronization at the
synchronization time between any two levels. So for example with two
jumps in refinement by a factor of two, there is a synchronization at
the end of the first two timesteps on level 2 (between level 1 and
level 2), a synchronization after the next two timesteps on level 2
(again between level 1 and level 2), and then a synchronization
between level 0 and level 1. In the new method, we always call the
synchronization at the end of every timestep *on the finest level
only*, and we simultaneously do the synchronization *on every
level*. The timestep :math:`\Delta t_f` in the flux register is just the
timestep on the finest level. (If this is unclear, give it a sanity
check: when the sum of all flux register totals is added up, the level
0 contribution will have a factor of :math:`\Delta t` equal to the coarse
grid timestep since the sum of the timesteps on the finest level over
the entire advance must equal the level 0 timestep. So, the final
contribution from the flux register is the same as if we had saved up
the flux mismatch until the end of the level 0 timestep.) The
synchronization no longer needs to be called at the end of any coarser
level’s timestep because it will already be up to date as a result of
the synchronizations applied at the end of the fine level timesteps.

.. [1]
   In the absence of a global field like
   the gravitational potential, this would only need to be done on the
   coarse level, as we always assume that the solution on the fine grid is
   correct and average it down to the coarse grid. In Castro we do it by
   default on the fine level too in anticipation of the fact that gravity
   is a common component of many of our production science
   simulations. This could be generalized so that if you aren’t using any
   global force fields, you don’t bother updating the fine level. If this
   is important to the science you want to do, please let the Castro developers know and we can look into it.

.. [2]
   in general it may be desirable for this to be a
   source-term specific setting, so that some source terms that are cheap
   or physically important are re-computed after a synchronization can be
   set to update, while others can be disabled. If this is important for
   your science application, please let the developers know, as this would
   be a straightforward extension of the current architecture.

.. [3]
   If this scheme
   is generalized to higher-order methods, in principle all one would need
   to do is integrate the fluxes until :math:`\Delta t / 2`, which is what we are
   doing here for the constant-in-time flux case.
***********
Input Files
***********

The Castro executable uses an inputs file at runtime to set and
alter the behavior of the algorithm and initial conditions.  

Runtime parameters take the form ``namespace.parameter = value`` ,
where *namespace* identifies the major code component.  Some
parameters take multiple values separated by spaces.

Typically
named ``inputs``, it 

  * Sets the AMReX parameters for gridding, refinement, etc., through the
    ``geometry`` and ``amr`` namespaces.

  * Enables different physics behaviors through the ``castro`` namespace.

  * Sets the problem-specific runtime parameters through the ``problem`` namespace.

  * Sets any Microphysics runtime parameters, through the various namespaces
    defined in Microphysics (like ``eos``, ``integrator``, ``network``, ...).

.. warning:: Because the inputs file is handled by the C++ portion
   of the code, any quantities you specify in scientific notation,
   must take the form ``1.e5`` and not ``1.d5``—the ``d``
   specifier is not recognized.


.. note::

   Additionally, note that in Castro, all quantities are in CGS units.


Common inputs Options
=====================


Problem Geometry
----------------

The ``geometry`` namespace is used by AMReX to define the
computational domain. The main parameters here are:

  * ``geometry.prob_lo``: physical location of low corner of the
    domain (type: ``Real``; must be set)

    Note: a number is needed for each dimension in the problem.

  * ``geometry.prob_hi``: physical location of high corner of the
    domain (type: ``Real``; must be set)

    Note: a number is needed for each dimension in the problem.

  * ``geometry.coord_sys``: coordinate system, 0 = Cartesian,
    1 = :math:`r`-:math:`z` (2-d only), 2 = spherical (1-d only) (must be set)

  * ``geometry.is_periodic``: is the domain periodic in this direction?
    0 if false, 1 if true (default: ``0 0 0``)

    Note: an integer is needed for each dimension in the problem.

  * ``castro.center``: physical location of problem center on the
    domain (type: ``Real``; default: ``0.0 0.0 0.0``). The problem
    center is used for gravity, rotation, and some other quantities.
    This is not necessarily the geometric center of the domain—often
    you should choose it to coincide with the center of mass of your
    system. See § :ref:`soft:prob_params` for more details.

   Note: a number is needed for each dimension in the problem.

As an example, the following::

    geometry.prob_lo = 0 0 0
    geometry.prob_hi = 1.e8 2.e8 2.e8
    geometry.coord_sys = 0
    geometry.is_periodic = 0 1 0
    castro.center = 5.e7 1.e8 1.e8

This defines the domain to run from :math:`(0,0,0)` at the lower left to
:math:`(10^8,\, 2\times 10^8,\, 2\times 10^8)` at the upper right in physical
space, specifies a Cartesian geometry, and makes the domain periodic
in the :math:`y`-direction only. The problem center is set to be halfway in
between the lower left and upper right corners.

Domain Boundary Conditions
--------------------------

Boundary conditions are specified using integer keys that are interpreted
by AMReX. The runtime parameters that we use are:

  * ``castro.lo_bc``: boundary type of each low face (must be set)

  * ``castro.hi_bc``: boundary type of each high face (must be set)

The valid boundary types are:

.. table:: boundary condition types
   :align: center

   +------------------------+-----------------+
   | 0: Interior / Periodic | 3: Symmetry     |
   +------------------------+-----------------+
   | 1: Inflow              | 4: Slip Wall    |
   +------------------------+-----------------+
   | 2: Outflow             | 5: No Slip Wall |
   +------------------------+-----------------+

.. note:: ``castro.lo_bc`` and ``castro.hi_bc`` must be consistent
   with ``geometry.is_periodic``—if the domain is periodic in a
   particular direction then the low and high bc’s must be set to 0
   for that direction.

As an example, the following::

    castro.lo_bc = 1 4 0
    castro.hi_bc = 2 4 0

    geometry.is_periodic = 0 0 1

This defines a problem with inflow (1) in the low-\ :math:`x` direction,
outflow (2) in the high-\ :math:`x` direction, slip wall (4) on
the low and high :math:`y`-faces, and periodic in the :math:`z`-direction.
See § :ref:`soft:phys_bcs`.

Resolution
----------

The grid resolution is specified by defining the resolution at the
coarsest level (level 0) and the number of refinement levels and
factor of refinement between levels. The relevant parameters are:

  * ``amr.n_cell``: number of cells in each direction at the coarsest
    level (integer :math:`> 0`; must be set)

  * ``amr.max_level``: number of levels of refinement above the
    coarsest level (integer :math:`\geq 0`; must be set)

  * ``amr.ref_ratio``: ratio of coarse to fine grid spacing
    between subsequent levels (2 or 4; must be set)

  * ``amr.regrid_int``: how often (in terms of number of steps) to
    regrid (integer; must be set)

  * ``amr.regrid_on_restart``: should we regrid immediately after
    restarting? (0 or 1; default: 0)

.. note:: if ``amr.max_level = 0`` then you do not need to set
   ``amr.ref_ratio`` or ``amr.regrid_int``.

Some examples::

    amr.n_cell = 32 64 64

would define the domain to have 32 cells in the :math:`x`-direction, 64 cells
in the :math:`y`-direction, and 64 cells in the :math:`z`-direction *at the
coarsest level*. (If this line appears in a 2D inputs file then the
final number will be ignored.)

::

    amr.max_level = 2

would allow a maximum of 2 refined levels in addition to the coarse
level. Note that these additional levels will only be created only if
the tagging criteria are such that cells are flagged as needing
refinement. The number of refined levels in a calculation must be
:math:`\leq` ``amr.max_level``, but can change in time and need not
always be equal to ``amr.max_level``.

::

    amr.ref_ratio = 2 4

would set factor of 2 refinement between levels 0 and 1, and factor of 4
refinement between levels 1 and 2. Note that you must have at least
``amr.max_level`` values of ``amr.ref_ratio`` (Additional values
may appear in that line and they will be ignored).

::

    amr.regrid_int = 2 2

tells the code to regrid every 2 steps. Thus in this example, new
level 1 grids will be created every 2 level-0 time steps, and new
level 2 grids will be created every 2 level-1 time steps. If
``amr.regrid_int`` :math:`<` 0 for any level, then regridding starting at that
level will be disabled. If ``amr.regrid_int = -1`` only, then we
never regrid for any level. Note that this is not compatible with
``amr.regrid_on_restart = 1``.


Other parameters
----------------

There are a large number of solver-specific runtime parameters. We describe these
together with the discussion of the physics solvers in later chapters.
******************
Runtime Parameters
******************

Introduction to Runtime Parameters
==================================

Castro runtime parameters are set
in the inputs file and managed by the AMReX ``ParmParse``
class. For Castro-specific parameters, we list the runtime
parameters in a file ``_cpp_parameters`` and generate the
C++ code and headers at compile time.

The behavior of the network, EOS, and other microphysics routines are
controlled by a different set of runtime parameters. These parameters are defined
in plain-text files ``_parameters`` located in the different
directories that hold the microphysics code. At compile time, a
a make function locates all
of the ``_parameters`` files that are needed for the given choice
of network, integrator, and EOS, and creates the ``extern_parameters.H`` and
``extern_parameters.cpp`` files that manage the parameters.  These
are set at runtime via the inputs file.

Castro-specific parameters
--------------------------

The Castro parameters that control the behavior of the code and
physics modules are listed in ``_cpp_parameters`` and take the form of::

    # comment describing the parameter
    name   type   default   need in Fortran?   ifdef

Here,

  * `name` is the name of the parameter that will be looked for
    in the inputs file.

    The name can actually take the form of ``(a, b)``, where ``a`` is
    the name to be used in the inputs file where the parameter is set
    and ``b`` is the name used within the Castro C++ class.  It is not
    recommended to name new parameters with this functionality—this
    was implemented for backwards compatibility.

  * `type` is one of int, Real, or string

  * `default` is the default value of the parameter.

The next columns are optional, but you need to fill in all of the
information up to and including any of the optional columns you need
(e.g., if you are going to provide the fortran name, you also need to
provide "need in Fortran?" and "ifdef".

  * `need in Fortran?` is ``y`` if the runtime parameter should be
    made available in Fortran (through ``meth_params_module``).
    Note: this option is deprecated.

  * `ifdef` provides the name of a preprocessor name that should
    wrap this parameter definition—it will only be compiled in if that
    name is defined to the preprocessor.

Finally, any comment (starting with ``#``) immediately before the
parameter definition will be used to generate the documentation
describing the parameters.

Microphysics/extern parameter format
------------------------------------

The microphysics/extern parameter definitions take the form of::

    # comment describing the parameter
    name              data-type       default-value      priority

Here, the `priority` is simply an integer. When two directories
define the same parameter, but with different defaults, the version of
the parameter with the highest priority takes precedence. This allows
specific implementations to override the general parameter defaults.


Parameters by Namespace
=======================

.. toctree::

   runtime_parameters
**********************
Introduction to Castro
**********************

Castro is a adaptive mesh, radiation/MHD hydrodynamics code that is
designed to model astrophysical reacting flows on massively parallel
computers.

Castro's major capabilities:

  * 1-, 2-, and 3-dimensional unsplit, 2nd-order finite-volume
    hydrodynamics; 4th order hydro for uniform grids.
    (see :ref:`ch:hydro`)

  * 3-dimension constrained transport ideal MHD (single level only currently)
    (see :ref:`ch:mhd`)

  * multigroup flux-limited diffusion radiation hydrodynamics
    (see :ref:`ch:radiation`)

  * generalized retry mechanism for recovering from physical
    violations over a timestep (see :ref:`ch:retry`)

  * adaptive mesh refinement with subcycling; jumps of 2x and 4x
    between levels (see :ref:`ch:amr`)

  * arbitrary equation of state (provided by the companion StarKiller
    Microphysics project)

  * general nuclear reaction networks

  * explicit thermal diffusion (see :ref:`ch:diffusion`)

  * full Poisson gravity (with isolated boundary conditions)
    and a conservative energy formulation (see :ref:`ch:gravity`)

  * rotation (in the co-rotating frame) in 2-d axisymmetric and 3-d
    (see :ref:`ch:rotation`)

  * spectral deferred corrections time integration for coupling hydro
    and reactions (see :ref:`ch:sdc`)

  * parallelization via MPI + OpenMP or MPI + CUDA


Development Model
=================

Castro is developed on github (https://github.com/amrex-astro/Castro
). The ``main`` branch is stable and can be used for day-to-day
science.  New changes are made via pull requests to the
``development`` branch.  This is where the ongoing regression testing
is done (both on CPU and GPU).

At the start of each month, we merge ``development`` → ``main`` and
apply a tag of the form ``YY.MM`` (e.g. ``20.02`` for Feb. 2020).  We
also create a github release and mint a Zenodo DOI using the
information in the ``.zenodo.json`` file at the root level.

Castro "core developers" are those who have made substantial code
contributions (details are in the main ``README.md``).  These
developers are coauthors on the Zenodo DOI and of any papers
describing Castro generally (science papers coauthors are decided by
the science paper lead).


Units and Conventions
=====================

Castro works in CGS units unless otherwise specified.
:numref:`table:units` shows some of the common symbols / names used
throughout the code documentation and papers.

.. _table:units:
  
.. table:: Common quantities and units.

   +-----------------------+-----------------------+-----------------------+
   | name                  | units                 | description           |
   +=======================+=======================+=======================+
   | :math:`t`             | s                     | time                  |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\rho`          | :math:`\gcc`          | mass density          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\ub`           | :math:`\cms`          | velocity vector       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`p`             | :math:`\presunit`     | pressure              |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\gb`           | :math:`\accelunit`    | gravitational         |
   |                       |                       | acceleration          |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\Sb`           | varies                | source term           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`E`             | :math:`\ergg`         | specific total energy |
   +-----------------------+-----------------------+-----------------------+
   | :math:`e`             | :math:`\ergg`         | specific internal     |
   |                       |                       | energy                |
   +-----------------------+-----------------------+-----------------------+
   | :math:`T`             | :math:`K`             | temperature           |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\kth`          | :math:`\mathrm{erg~cm | thermal conductivity  |
   |                       | ^{-1}~s^{-1}~K~{-1}}` |                       |
   +-----------------------+-----------------------+-----------------------+
   | :math:`X_k`           | –                     | mass fraction of      |
   |                       |                       | species :math:`k`     |
   +-----------------------+-----------------------+-----------------------+
   | :math:`\omegadot_k`   | :math:`\mathrm{s^{-1} | species creation rate |
   |                       | }`                    | (from reactions)      |
   +-----------------------+-----------------------+-----------------------+

Physical constants, again using the CGS system are available
in ``Microphysics/constants/``.


************************************
Self-Consistent Field Initialization
************************************

Introduction
============

The Hachisu self-consistent field (SCF) method is a way to generate
equilibrium rotating (or non-rotating) initial configurations. It can
generate single or multiple star systems. The SCF method was originally
developed in the 1960s, but it was a variant proposed by Hachisu in 1986
(:cite:`hachisu:1986a,hachisu:1986b`) that ended up becoming the
most popular technique. The SCF method was originally developed
for rapidly rotating single stars, but was soon extended to apply
to rotating binary systems, and it has been applied to construct
initial conditions by several groups studying binary white dwarf
or neutron star systems using several types of rotation laws
(:cite:`newtohline:1997,swc:2000,motl:2002,dsouza:2006,motl:2007,even:2009,kadam:2018,yoshida:2018`).
It has also been used for constructing toroidal configurations
(:cite:`kim:2016`).
The technique assumes a uniform temperature and composition (more generally,
a barotropic equation of state), self-gravitation represented by the
Poisson equation, and a well-defined rotation law (often rigid-body rotation).
The user is required to specify three quantities: the maximum density of
the equilibrium star(s), and two points on the stellar surface. For a
single star this is usually a point on the equator and a point on a pole.
For a detached binary system this corresponds to the inner and outer points
of the teardrop-shape configuration along the axis joining the binary.

At present we only support generating a single star using the Hachisu
SCF method, but we plan to extend this in a later release to rotating
toroidal configurations and binary star systems.

We note that while the Hachisu SCF method iteratively solves the
integral form of the Euler equations (essentially, it solves the Bernoulli
equation from classical fluid dynamics), other approaches have been used
in the literature. Usually these directly solve the coupled Poisson and
Bernoulli equations rather than relying on indirect iterative coupling
between them. See, for example, :cite:`eriguchi:1985,fujisawa:2015` and
:cite:`clement:1974,aksenov:1994`. These alternative methods are generally
more powerful and promise faster convergence, but are also more difficult
to implement. See also :cite:`jackson:2005` for yet another approach.



Usage and Code Parameters
=========================

To use the SCF initialization technique, you must compile with
``USE_GRAV=TRUE`` and ``USE_ROTATION=TRUE``, and enable the method
with ``castro.do_scf_initial_model`` = ``1``. You are responsible
for initializing the grid with an initial guess at the mass distribution.
This guess does not need to be accurate -- for example, just initializing
with a uniform spherical mass distribution should be sufficient.
(See ``Exec/scf_tests/single_star``.) However, it is important that
the grid is isothermal and has a uniform composition. This temperature
and composition will be retained for the final equilibrium configuration.

Several code parameters are available for controlling problem initialization
with SCF:

- ``castro.scf_maximum_density``: the target maximum density on the domain
- ``castro.scf_equatorial_radius``: the target equatorial radius of the star
- ``castro.scf_polar_radius``: the target polar radius of the star
- ``castro.scf_relax_tol``: tolerance required for SCF convergence

The first three options are required and must be set. One limitation of this
method is that (to our knowledge) there is no known way to specify more natural
parameters such as the total mass of the star.


Single Star Algorithm
=====================

The Bernoulli equation tells us

.. math::
   \Phi + \phi + H = \mathrm{constant}

where :math:`\Phi` is the gravitational potential, :math:`\phi` is the
rotational potential (:math:`\phi = -\omega^2 R^2` for rigid-body rotation,
where :math:`R` is the radius in the equatorial plane), and :math:`H` is the
enthalpy.

We denote the prescribed maximum density as :math:`\rho_0` and the equatorial
and polar radii as :math:`r_A` and :math:`r_B` respectively (the A and B indices
are chosen to be consistent with Hachisu). At these radii the enthalpy should
vanish, so we have:

.. math::
   C = \Phi_A + \phi_A

.. math::
   C = \Phi_B + \phi_B

Note that :math:`C = C_A = C_B` is constant everywhere on the body.

Then, given the value of a Bernoulli constant, we can invert the Bernoulli equation to
obtain the enthalpy:

.. math::
   H = C - \Phi - \phi

Given an enthalpy, we can call the equation of state given the enthalpy (and
composition and temperature) as an input to obtain the density. (Note that
this inversion generally requires a Newton-Raphson iteration in realistic
equations of state.)

.. math::
   \rho = \rho(T, H, X)

As one more housekeeping item, we'll notate the rotational potential as

.. math::
   \phi = \omega^2 \psi

where :math:`\psi = -R^2` for rigid body rotation.

The above ingredients are all that is needed to construct the algorithm for
obtaining an equilibrium rotating single star. Given an initial density distribution,
:math:`\rho^n`, gravitational field, :math:`\Phi^n`, and rotational field,
:math:`\phi^n`, we first calculate an updated guess for the rotation frequency
:math:`\omega`:

.. math::
   \omega^{n+1} = \sqrt{\frac{\Phi_B^n - \Phi_A^n}{\psi_A^n - \psi_B^n}}

which simply involves finding :math:`\Phi` and :math:`\psi` at these vanishing points.

With the updated rotation frequency, we can reconstruct the rotational potential
:math:`\phi`, and then update the enthalpy everywhere on the domain as:

.. math::
   H^{n+1} = C - \Phi - \Phi_R

However, we want to guarantee that the maximum density on the domain is fixed. Given
that this maximum density corresponds to a maximum enthalpy,

.. math::
   H_0 = H(\rho_0, T, X)

we can rescale all of the updated enthalpies such that the maximum is fixed:

.. math::
   H^{n+1} \rightarrow H^{n+1} \left( \frac{H_0}{H^{n+1}_{\mathrm{max}}} \right)

and then invert the EOS to obtain :math:`\rho^{n+1}`. Given the new density
distribution, we can then update the gravitational potential, :math:`\Phi^{n+1}`,
by solving the Poisson equation. This procedure is iterated until no zone
changes its density by more than a factor of ``castro.scf_relax_tol``.
*******
Preface
*******

Welcome to the Castro User’s Guide!

In this User’s Guide we describe how to download and run Castro, a
massively parallel code that solves the multicomponent compressible
hydrodynamic equations for astrophysical flows including self-gravity,
nuclear reactions and radiation. Castro uses an Eulerian grid and
incorporates adaptive mesh refinement (AMR). Our approach to AMR uses
a nested hierarchy of logically-rectangular grids with simultaneous
refinement in both space and time, utilizing the
AMReX library :cite:`amrex-joss`.

The core algorithms in Castro are described in a series of papers:

  * *CASTRO: A New Compressible Astrophysical Solver. I. Hydrodynamics
    and Self-gravity*, A. S. Almgren, V. E. Beckner, J. B. Bell,
    M. S. Day, L. H. Howell, C. C. Joggerst, M. J. Lijewski,
    A. Nonaka, M. Singer, & M. Zingale, 2010, ApJ, 715, 1221
    http://dx.doi.org/10.1088/0004-637X/715/2/1221

  * *CASTRO: A New Compressible Astrophysical Solver. II. Gray
    Radiation
    Hydrodynamics*, W. Zhang, L. Howell, A. Almgren, A. Burrows,
    & J. Bell, 2011, ApJS, 196, 20
    http://dx.doi.org/10.1088/0067-0049/196/2/20

  * *CASTRO: A New Compressible Astrophysical Solver. III. Multigroup
    Radiation
    Hydrodynamics*, W. Zhang, L. Howell, A. Almgren, A. Burrows, J. Dolence,
    & J. Bell, 2013, ApJS, 204, 7
    http://dx.doi.org/10.1088/0067-0049/204/1/7

Improvements to the gravity solver and rotation were described in:

  * *Double White Dwarf Mergers on Adaptive Meshes I. Methodology and
    Code
    Verification*, M. P. Katz, M. Zingale, A. C. Calder, F. D. Swesty,
    A. S. Almgren, & W. Zhang, 2016, ApJ, 819, 94.
    http://dx.doi.org/10.3847/0004-637X/819/2/94

The 4th order reactive hydrodynamics spectral deferred corrections solver
was described in:

  * *Improved Coupling of Hydrodynamics and Nuclear Reactions via Spectral Deferred Corrections*
    M. Zingale, M. P. Katz, J. B. Bell, M. L. Minion, A. J. Nonaka, & W. Zhang,
    2019, ApJ, 886, 105.
    https://ui.adsabs.harvard.edu/link_gateway/2019ApJ...886..105Z/doi:10.3847/1538-4357/ab4e1d

The Castro GPU strategy and performance was described in:

  * *Preparing Nuclear Astrophysics for Exascale*
    M. P. Katz, A. Almgren, M. Barrios Sazo, K. Eiden, K. Gott, A. Harpole, J. M. Sexton, D. E. Willcox, W. Zhang, & M. Zingale
    2020, to appear in Proceedings of SC20
    https://ui.adsabs.harvard.edu/abs/2020arXiv200705218K/abstract


The development of AMReX library is led by the
Center for Computational Sciences and Engineering / Lawrence Berkeley
National Laboratory. Castro development is done collaboratively,
including the CCSE and Stony Brook University.

Castro *core developers* are those who have made substantial
contributions to the code. The process for becoming a core developer
is described in the `README.md <https://github.com/AMReX-Astro/Castro/blob/main/README.md>`_ in the Castro root directory.

Current Castro core developers are listed at https://amrex-astro.github.io/Castro/who.html

All Castro development takes place on the project’s github
page, https://github.com/AMReX-Astro/Castro

External contributions are welcomed. Fork the Castro repo, modify your
local copy, and issue a pull-request to the AMReX-Astro/Castro
project. Further guidelines are given in the `README.md
<https://github.com/AMReX-Astro/Castro/blob/main/README.md>`_ file.

Getting Help
============

We use github discussions to ask questions about the code and get help:

https://github.com/AMReX-Astro/Castro/discussions

You can also post issues on the github page to report bugs.


Acknowledging and Citing Castro
===============================

If you use Castro in your research, we would appreciate it if you
cited the relevant code papers describing its design, features, and
testing. A list of these can be found in the `CITATION
<https://github.com/AMReX-Astro/Castro/blob/main/CITATION>`_ file in
the root ``Castro/`` directory.

The development Castro is supported by the science application
interests of the contributors. There is a lot of effort behind the
scenes: testing, optimization, development of new features, bug
fixing, ..., that is often done under the radar. Nevertheless,
we are happy to volunteer our time to help new users come up to speed
with Castro. When significant new development / debugging for you
application is provided by a member of the Castro development
community, we would appreciate consideration of inviting the
developer(s) for co-authorship on any science paper that results.

**************************
Frequently Asked Questions
**************************

Compiling
=========

#. *Compiling fails giving me a cryptic message about a module not
   being found.*

   This usually indicates that the build system cannot find a source file.
   The source files are specified
   in the various ``Make.package`` files throughout the
   Castro directory hierarchy. make will look through the
   directories in the ``VPATH_LOCATIONS`` to find the files.

   There are 2 things you can do to check what’s happening. First, inspect
   the directories in ``VPATH_LOCATIONS``. This can be done via:

   ::

       make print-VPATH_LOCATIONS

   Next, ask make to tell you where it is finding each of the source
   files. This is done through a script ``find_files_vpath.py``
   that is hooked into Castro’s build system. You can run this as:

   ::

       make file_locations

   At the end of the report, it will list any files it cannot find in
   the vpath. Some of these are to be expected (like ``extern.f90``
   and ``buildInfo.cpp``—these are written at compile-time. But any
   other missing files need to be investigated.

#. *I’m still having trouble compiling. How can I find out what
   all of the make variables are set to?*

   Use:

   ::

       make help

   This will tell you the value of all the compilers and their options.

#. *How can I check to make sure the function signatures defined
   in C are consistent with their implementations in Fortran?*

   Use:

   ::

       make typecheck

   This will compile the code and report on any mismatched function signatures.

.. _debugging_backtrace:

Debugging
=========

#. *Castro crashes with a floating point exception—how can
   I get more information?*

   The best thing to do is to recompile the code with ``TEST=TRUE``
   set in the ``GNUmakefile``. This will have AMReX catch the
   signals raised in both C++ and Fortran functions. Behind the
   scenes, this defines the ``AMREX_TESTING`` preprocessor flag, which
   will initialize memory allocated in fabs or multifabs to
   signaling NaNs (sNaN), and use the ``BLBackTrace::handler()``
   function to handle various signals raised in both C and Fortran
   functions. This is a Linux/UNIX capability. This gives us a chance
   to print out backtrace information. The signals include seg fault,
   floating point exceptions (NaNs, divided by zero and overflow), and
   interruption by the user and system. What signals are handed to
   AMReX are controlled by AMReX(e.g., using interruption by the
   user, this was once used to find an MPI deadlock.) It also includes
   the ``AMREX_ASSERTION`` statements if ``USE_ASSERTION=TRUE`` or
   ``DEBUG=TRUE``.

   The AMReX parameters that affect the behavior are:

   -  ``amrex.fpe_trap_invalid``

   -  ``amrex.fpe_trap_zero``

   -  ``amrex.fpe_trap_overflow``

   For further capabilities, you can get 
   more information than the backtrace of the call stack info by
   instrumenting the code.  Here is an
   example. You know the line ``Real rho = state(cell,0);`` is
   causing a segfault. You could add a print statement before that.
   But it might print out thousands (or even millions) of line before
   it hits the segfault. Instead, you could

   .. code:: c++

             std::ostringstream ss;
             ss << "state.box() = " << state.box() << " cell = " << cell;
             BL_BACKTRACE_PUSH(ss.str()); // PUSH takes std::string

             Real rho = state(cell,0);  // state is a Fab, and cell is an IntVect.

   The "print" prints to a stack of string, not stdout. When it hits
   the segfault, you will only see the last print out in the backtrace
   file (e.g. ``BackTrace.0``).

   You may need to include the header ``AMReX_BLBackTrace.H``.

#. *How can I monitor the state in a zone from the C side
   at various points in the evolution?*

   Given a MultiFab ``mf``, you can dump out the state as:

   ::

           print_state(mf, IntVect(D_DECL(10, 20, 30)));

   Here, the IntVect has the dimension that we were compiled with
   (and this is handled through the preprocessor ``D_DECL``). In
   this case, we are inspecting zone (10, 20, 30), in the global index
   space. Note that since a multifab exists only on a single level, the
   integer indices here refer to the global index space on that level.

#. *What if I want to see all the data in a FArrayBox?*

   You can simply output a FAB to ``std::cout``. Imagine that you
   are in an MFIter loop, with a MultiFab ``mf``:

   ::

           S = FArrayBox& mf[mfi];
           std::cout << S << std::endl;

   This will output the contents on the FAB, one zone per line.

Profiling
=========

#. *How can I get line-by-line profiling information?*

   With the GNU compliers, you can enabling profiling with gprof
   by compiling with

   ::

         USE_GPROF=TRUE

   in your ``GNUmakefile``.

   When you run, a file named ``gmon.out`` will be produced. This can
   be processed with gprof by running:

   ::

         gprof exec-name

   where *exec-name* is the name of the executable. More detailed
   line-by-line information can be obtained by passing the -l
   argument to gprof.

Managing Runs
=============

#. *How can I force the running code to output, even it the plot or
   checkpoint interval parameters don’t require it?*

   Create a file called ``dump_and_continue``, e.g., as:

   ::

       touch dump_and_continue

   This will force the code to output a checkpoint file that can be used
   to restart. Other options are ``plot_and_continue`` to output
   a plotfile, ``dump_and_stop`` to output a checkpoint file
   and halt the code, and ``stop_run`` to simply stop the code.


   .. note::

      The parameter ``amr.message_int`` controls how often the
      existence of these files is checked; by default it is 10, so the
      check will be done at the end of every timestep that is a
      multiple of 10.  Set that to 1 in your inputs file if you’d like
      it to check every timestep.

#. *How can I output plotfiles in single precision?*

   The AMReX runtime parameter:

   ::

       fab.format = NATIVE_32

   controls this (put this in your inputs file). Note: checkpoint files are unaffected
   by this and will always be written out in the native precision (the ‘fab.format‘ parameter
   is overridden in the checkpoint code in AMReX).

#. *How can I check the compilation parameters of a Castro executable?*

   The build information (including git hashes, modules, EoS, network, etc.) can be displayed by running the executable as 

   ::

       ./Castro.exe --describe

.. _ch:faq:vis:

Runtime Errors
==============

#. *When running with retries, Castro requests too many substeps
   and crashes.*

   This can occur due to CFL violations or negative densities.  If
   there are density resets, try running with
   ``castro.limit_fluxes_on_small_dens = 1``.  This will use a flux
   limiter to prevent the density from going negative.

Visualization
=============

#. *When I try to use Amrvis with the Nvidia driver, all I see is
   black—no data. How do I fix this?*

   You need to edit your xorg.conf file (usually found in /etc/X11/
   to enable the Dac8Bit option. The section will look like:

   ::

       Section "Device"
           Identifier     "Device0"
           Driver         "nvidia"
           VendorName     "NVIDIA Corporation"
           Option         "Dac8bit" "True"
       EndSection

   If you don’t already have an ``xorg.conf`` then you can create one
   by running ``nvidia-xconfig`` first.
.. Castro documentation main file, created by
   sphinx-quickstart on Mon Dec 25 18:42:54 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

********************************************************
Castro: an adaptive mesh compressible hydrodynamics code
********************************************************

https://github.com/amrex-astro/Castro

.. toctree::
   :maxdepth: 1
   :caption: Castro basics

   Preface
   Introduction
   getting_started
   inputs
   rp_intro
   mpi_plus_x
   FlowChart
   software
   problem_setups
   timestepping
   creating_a_problem
   io
   regridding
   visualization
   faq

.. toctree::
   :maxdepth: 1
   :caption: Castro reference

   build_system
   debugging
   Hydrodynamics
   mhd
   gravity
   diffusion
   rotation
   sponge
   radiation
   Particles
   EOSNetwork
   sdc
   AMR
   ConvertCheckpoint
   self_consistent_field
   MAESTRO_restart
   Verification
   development

.. toctree::
   :maxdepth: 1
   :caption: API

   filelist
   classlist
   .. namespacelist

.. toctree::
   :caption: References

   zreferences

Indices and tables
==================

* :ref:`genindex`

.. * :ref:`modindex`

* :ref:`search`

.. _ch:io:

**********
Outputting
**********

Restart Capability
------------------

.. index:: amr.check_file, amr.check_int, amr.check_per, amr.restart
.. index:: amr.checkpoint_files_output, amr.check_nfiles, amr.checkpoint_on_restart
.. index:: castro.grown_factor

Castro has a standard sort of checkpointing and restarting capability.
In the inputs file, the following options control the generation of
checkpoint files (which are really directories):

  * ``amr.check_file``: prefix for restart files (text;
    default: chk)

  * ``amr.check_int``: how often (by level 0 time steps) to
    write restart files (integer :math:`> 0`; default: -1)

  * ``amr.check_per``: how often (by simulation time) to
    write restart files (Real :math:`> 0`; default: -1.0)

    Note that ``amr.check_per`` will write a checkpoint at the first
    timestep whose ending time is past an integer multiple of this
    interval.  In particular, the timestep is not modified to match
    this interval, so you won’t get a checkpoint at exactly the time
    you requested.

  * ``amr.restart``: name of the file (directory) from which to
    restart (Text; not used if not set)

  * ``amr.checkpoint_files_output``: should we write
    checkpoint files? (0 or 1; default: 1)

    If you are doing a scaling study then set
    ``amr.checkpoint_files_output`` = 0 so you can test scaling of the
    algorithm without I/O.

  * ``amr.check_nfiles``: how parallel is the writing of
    the checkpoint files? (Integer :math:`\geq 1`; default: 64)

    See the section :ref:`sec:parallel_io` for more details on parallel I/O and the
    ``amr.check_nfiles`` parameter.

  * ``amr.checkpoint_on_restart``: should we write a
    checkpoint immediately after restarting? (0 or 1; default: 0)

  * ``castro.grown_factor``: factor by which domain has been
    grown (Integer :math:`\geq 1`; default: 1)

.. note:: You can specify both ``amr.check_int`` or ``amr.check_per``,
   if you so desire; the code will print a warning in case you did
   this unintentionally. It will work as you would expect – you will
   get checkpoints at integer multiples of ``amr.check_int`` timesteps
   and at integer multiples of ``amr.check_per`` simulation time
   intervals.

   ``amr.plotfile_on_restart`` and ``amr.checkpoint_on_restart``
   require amr.regrid_on_restart to be in effect.

As an example::

    amr.check_file = chk_run
    amr.check_int = 10

means that restart files (really directories) starting with the prefix
“chk_run” will be generated every 10 level-0 time steps. The
directory names will be ``chk_run00000``, ``chk_run00010``,
``chk_run00020``, etc.

If instead you specify::

    amr.check_file = chk_run
    amr.check_per = 0.5

then restart files (really directories) starting with the prefix
“chk_run” will be generated every 0.1 units of
simulation time. The directory names will be ``chk_run00000``,
``chk_run00043``, ``chk_run00061``, etc, where :math:`t = 0.1` after
43 level-0 steps, :math:`t = 0.2` after 61 level-0 steps, etc.

To restart from ``chk_run00061``, for example, then set::

    amr.restart = chk_run00061

.. _sec:PlotFiles:


Plotfile Outputting
-------------------

.. index:: amr.plot_files_output, amr.plotfile_on_restart, amr.write_plotfile_with_checkpoint

Castro has two levels of plotfiles, `regular` plotfiles and `small`
plotfiles.  The idea behind this distinction is that we can output a
small number of variables very frequently in the small plotfiles and
output a large number (or all variables) less frequently.  This helps
keep the data sizes down while allowing for fine-grained temporal
analysis of important quantities.


A few general controls determines whether we want to output plotfiles and when:

  * ``amr.plot_files_output`` : this is set to 1 to output plotfiles

  * ``amr.plotfile_on_restart`` : set this to 1 to dump out a plotfile
    immediately when we restart.

  * ``amr.write_plotfile_with_checkpoint`` : always output a plotfile
    when we dump a checkpoint file.

.. index:: amr.plot_file, amr.plot_per, amr.plot_int

The frequency of outputting and naming of regular plotfiles is
controlled by:

  * ``amr.plot_file`` : this is the base name for the plotfile,
    e.g. ``plt``.

  * ``amr.plot_per`` : this is the amount of simulation time between
    plotfile output

    .. note:: ``amr.plot_per`` will write a plotfile at the first
       timestep whose ending time is past an integer multiple of this
       interval.  In particular, the timestep is not modified to match
       this interval, so you won’t get a checkpoint at exactly the time
       you requested.

  * ``amr.plot_int`` this is the number of timesteps between plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.

.. index:: amr.small_plot_file, amr.small_plot_per, amr.small_plot_int

Similarly, the frequency of outputting and naming of small plotfiles
is controlled by:

  * ``amr.small_plot_file`` : this is the base name for the small plotfile,
    e.g. ``smallplt``.

  * ``amr.small_plot_per`` : this is the amount of simulation time between
    small plotfile output

  * ``amr.small_plot_int`` this is the number of timesteps between small plotfiles.
    Set this to -1 to rely on the simulation-time-based outputting.

Additional output options control how the I/O is done:

  * ``amr.plot_nfiles``: how parallel is the writing of the
    plotfiles? (Integer :math:`\geq 1`; default: 64)

    See the Software Section for more details on parallel I/O and the
    ``amr.plot_nfiles`` parameter.

All the options for ``amr.derive_plot_vars`` are kept in
``derive_lst`` in ``Castro_setup.cpp``. Feel free to look at
it and see what’s there.

.. note:: You can specify both ``amr.plot_int`` or ``amr.plot_per``,
   if you so desire; the code will print a warning in case you did
   this unintentionally. It will work as you would expect – you will
   get plotfiles at integer multiples of amr.plot_int timesteps and at
   integer multiples of amr.plot_per simulation time intervals.

As an example::

    amr.plot_file = plt_run
    amr.plot_int = 10

means that plot files (really directories) starting with the prefix
“plt_run” will be generated every 10 level-0 time steps. The
directory names will be ``plt_run00000``, ``plt_run00010``,
``plt_run00020``, etc.

If instead you specify::

    amr.plot_file = plt_run
    amr.plot_per = 0.5

then restart files (really directories) starting with the prefix
“plt_run” will be generated every 0.1 units of simulation time. The
directory names will be ``plt_run00000``, ``plt_run00043``,
``plt_run00061``, etc, where :math:`t = 0.1` after 43 level-0 steps, :math:`t =
0.2` after 61 level-0 steps, etc.


Controlling What’s in the PlotFile
----------------------------------

.. index:: amr.plot_vars, amr.derive_plot_vars

There are a few options that can be set at runtime to control what
variables appear in the regular plotfile.

  * ``amr.plot_vars``: this controls which of the main
    state variables appear in the plotfile. The default is for all of
    them to be stored. But you can specify a subset by name, e.g.::

        amr.plot_vars = density

    to only store that subset.

  * ``amr.derive_plot_vars``: this controls which of the derived
    variables to be stored in the plotfile. Derived variables are
    created only when the plotfile is being created, using the
    infrastructure provided by AMReX to register variables and the
    associated Fortran routine to do the deriving (``Derive_nd.F90``).

    By default, no derived variables are stored. You can store all
    derived variables that Castro knows about by doing::

       amr.derive_plot_vars = ALL

   or a subset by explicitly listing them, e.g.::

      amr.derive_plot_vars = entropy pressure

   To not output any derived variable,s this is set to ``NONE``.

.. index:: amr.small_plot_vars

For small plotfiles, the controls that lists the variables is:

  * ``amr.small_plot_vars`` : this is a list of which variables
    to include in the small plotfile.

  * ``amr.derive_small_plot_vars`` : this is a list of which derived
    variables to include in the small plotfile.


Plotfile Variables
------------------

Native variables
^^^^^^^^^^^^^^^^

These variables come directly from the ``StateData``, either the
``State_Type`` (for the hydrodynamic variables), ``Reactions_Type``
(for the nuclear energy generation quantities). ``PhiGrav_Type`` and
``Gravity_Type`` (for the gravity quantities), and ``Rad_Type`` (for
radiation quantities).


+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``density``                       | Mass density, :math:`\rho`                        | :math:`\gcc`                         |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``xmom``                          | x-momentum, :math:`(\rho u)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``ymom``                          | y-momentum, :math:`(\rho v)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``zmom``                          | z-momentum, :math:`(\rho w)`                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_E``                         | Total energy density                              | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_e``                         | Internal energy density                           | :math:`{\rm erg~cm^{-3}}`            |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Temp``                          | Temperature                                       | :math:`{\rm K}`                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_X``                         | Mass density of species X                         | :math:`\gcc`                         |
| (where X is any of the species    |                                                   |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``omegadot_X``                    | Creation rate of species X                        | :math:`{\rm s^{-1}}`                 |
| (where X is any of the species    | :math:`\omegadot_k = DX_k/Dt`                     |                                      |
| defined in the network)           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``enuc``                          | Nuclear energy generation rate / gram             | :math:`{\rm erg~g^{-1}~s^{-1}}`      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rho_enuc``                      | Nuclear energy generation rate density            | :math:`{\rm erg~cm^{-3}~s^{-1}}`     |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiGrav``                       | Gravitational potential                           | :math:`{\rm erg~g^{-1}}`             |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``grav_x``, ``grav_y``,           | Gravitational acceleration                        | :math:`{\rm cm~s^{-2}}`              |
| ``grav_z``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rmom``                          | Radial momentum (defined for                      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``lmom``                          | Angular momentum (:math:`\theta`; defined for     | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
|                                   | ``HYBRID_MOMENTUM``)                              |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pmom``                          | z-momentum (defined for ``HYBRID_MOMENTUM``)      | :math:`{\rm g~cm^{-2}~s^{-1}}`       |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Shock``                         | Shock flag (= 1 if a zone has a shock;            | --                                   |
|                                   | defined for ``SHOCK``)                            |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``rad``, ``rad0``, ``rad1``,      | Radiation energy density                          |                                      |
| ...                               | (for multigroup radiation, each group has its     |                                      |
|                                   | own variable)                                     |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



Derived variables
^^^^^^^^^^^^^^^^^

.. index:: castro.domain_is_plane_parallel

+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| variable name                     | description                                       | derive routine              | units                                   |
+===================================+===================================================+=============================+=========================================+
| ``angular_momentum_x``,           | Angular momentum / volume in the x, y, or z dir   | ``derangmomx``,             | :math:`{\rm g~cm^{-1}~s^{-1}}`          |
| ``angular_momentum_y``,           | computed as :math:`[(\rho \ub) \times {\bf r}]_n` | ``derangmomy``,             |                                         |
| ``angular_momentum_z``            | where :math:`{\bf r}` is the distance from        | ``derangmomz``              |                                         |
|                                   | ``center`` and :math:`n` is either x, y, or z     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_coeff``                    | Thermal diffusion coefficient,                    | ``derdiffcoeff``            | :math:`{\rm cm^2~s^{-1}}`               |
|                                   | :math:`\kth/(\rho c_v)`                           |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``diff_term``                     | :math:`\nabla\cdot(\kth\nabla T)`                 | ``derdiffterm``             | :math:`{\rm erg~cm^{-3}~s^{-1}}`        |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``divu``                          | :math:`\nabla \cdot \ub`                          | ``derdivu``                 | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_e``                        | Specific internal energy computed from the        | ``dereint2``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | conserved :math:`(\rho e)` state variable as      |                             |                                         |
|                                   | :math:`e = (\rho e)/\rho`                         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``eint_E``                        | Specific internal energy computed from the        | ``dereint1``                | :math:`{\rm erg~g^{-1}}`                |
|                                   | total energy and momentum conserved state as      |                             |                                         |
|                                   | :math:`e=[(\rho E)-\frac{1}{2}(\rho \ub^2)]/\rho` |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``entropy``                       | Specific entropy, :math:`s`, computed as          | ``derentropy``              | :math:`{\rm erg~g^{-1}~K^{-1}}`         |
|                                   | :math:`s = s(\rho, e, X_k)`, where `e` is         |                             |                                         |
|                                   | computed from :math:`(\rho e)`                    |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Ertot``                         | Total radiation energy density                    | ``derertot``                |                                         |
|                                   | (for multigroup radiation problems)               |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frcomx``, ``Frcomy``,           | Comoving radiation flux                           | ``Radiation.cpp``           |                                         |
| ``Frcomz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Frlabx``, ``Frlaby``,           | Lab-frame radiation flux                          | ``Radiation.cpp``           |                                         |
| ``Frlabz``                        |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``Gamma_1``                       | Adiabatic index,                                  | ``dergamma1``               | --                                      |
|                                   | :math:`d\log p/d\log \rho|_s`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``kineng``                        | Kinetic energy density,                           | ``derkineng``               | :math:`{\rm erg~cm^{-3}}`               |
|                                   | :math:`K = \frac{1}{2} |(\rho \ub)|^2`            |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``lambda``                        | Radiation flux limiter                            |                             | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``logden``                        | :math:`\log_{10} \rho`                            | ``derlogten``               | dimensionless, assuming :math:`\rho`    |
|                                   |                                                   |                             | is in CGS                               |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``MachNumber``                    | Fluid Mach number, :math:`|\ub|/c_s`              | ``dermachnumber``           | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``maggrav``                       | Gravitational acceleration magnitude              | ``dermaggrav``              | :math:`{\rm cm~s^{-2}}`                 |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magmom``                        | Momentum density magnitude,                       | ``dermagmom``               | :math:`{\rm g~cm^{-2}~s^{-1}}`          |
|                                   | :math:`|\rho \ub|`                                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvel``                        | Velocity magnitude, :math:`|\ub|`                 | ``dermagvel``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``magvort``                       | Vorticity magnitude, :math:`|\nabla\times\ub|`    | ``dermagvort``              | :math:`{\rm s^{-1}}`                    |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``pressure``                      | Total pressure, including ions, electrons,        | ``derpres``                 | :math:`{\rm dyn~cm^{-2}}`               |
|                                   | and radiation (for non radhydro problems)         |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``radvel``                        | Radial velocity (measured with respect to         | ``derradialvel``            | :math:`\cms`                            |
|                                   | ``center`` or vertical axis if                    |                             |                                         |
|                                   | ``domain_is_plane_parallel`` is set)              |                             |                                         |
|                                   | :math:`(xu + yv + zw)/r`                          |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``circvel``                       | Circumferential velocity (perpendicular to        | ``derradialvel``            | :math:`\cms`                            |
|                                   | ``radvel``.  If ``domain_is_plane_parallel`` is   |                             |                                         |
|                                   | set, then this is in the x-y plane                |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``soundspeed``                    | Sound speed                                       | ``dersoundspeed``           | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``StateErr``                      |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``thermal_cond``                  | Thermal conductivity, :math:`\kth`                | ``dercond``                 | :math:`{\rm erg~cm^{-1}~s^{-1}~K^{-1}}` |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``t_sound_t_enuc``                |                                                   | ``derenuctimescale``        | --                                      |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uminusc``                       | (only for 1D) x-velocity :math:`-` sound          | ``deruminusc``              | :math:`\cms`                            |
|                                   | speed                                             |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``uplusc``                        | (only for 1D) x-velocity + sound speed            | ``deruplusc``               | :math:`\cms`                            |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``X(q)``                          | Mass fraction of species q                        | ``derspec``                 | --                                      |
|                                   | :math:`X_k = (\rho X_k)/\rho`                     |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+
| ``x_velocity``,                   | Fluid velocity,                                   | ``dervel``                  | :math:`\cms`                            |
| ``y_velocity``,                   | :math:`\ub = (\rho \ub)/\rho`                     |                             |                                         |
| ``z_velocity``                    |                                                   |                             |                                         |
+-----------------------------------+---------------------------------------------------+-----------------------------+-----------------------------------------+


problem-specific plotfile variables
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

+-----------------------------------+---------------------------------------------------+--------------------------------------+
| variable name                     | description                                       | units                                |
+===================================+===================================================+======================================+
| ``analytic``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pi``                            |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``pioverp0``                      |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``primarymask``                   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``secondarymask``                 |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Terror``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``Texact``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_angular_momentum_x``,  |                                                   |                                      |
| ``inertial_angular_momentum_y``,  |                                                   |                                      |
| ``inertial_angular_momentum_z``   |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_momentum_x``,          |                                                   |                                      |
| ``inertial_momentum_y``,          |                                                   |                                      |
| ``inertial_momentum_z``           |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``inertial_radial_momentum_x``,   |                                                   |                                      |
| ``inertial_radial_momentum_y``,   |                                                   |                                      |
| ``inertial_radial_momentum_z``    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEff``                        |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_P``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``phiEffPM_S``                    |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+
| ``tpert``                         |                                                   |                                      |
+-----------------------------------+---------------------------------------------------+--------------------------------------+



Screen Output
-------------

There are several options that set how much output is written to the
screen as Castro runs:

  * ``amr.v``: verbosity of ``Amr.cpp`` (0 or 1; default: 0)

  * ``castro.v``: verbosity of ``Castro.cpp`` (0 or 1; default: 0)

  * ``gravity.v``: verbosity of ``Gravity.cpp`` (0 or 1; default: 0)

  * ``diffusion.v``: verbosity of ``Diffusion.cpp`` (0 or 1;
    default: 0)

  * ``mg.v``: verbosity of multigrid solver (for gravity) (allow
    values: 0, 1, 2, 3, 4; default: 0)

  * ``amr.grid_log``: name of the file to which the grids are
    written (text; not used if not set)

  * ``amr.run_log``: name of the file to which certain output is
    written (text; not used if not set)

  * ``amr.run_log_terse``: name of the file to which certain
    (terser) output is written (text; not used if not set)

  * ``castro.do_special_tagging``: allows the user to set a special
    flag based on user-specified criteria (0 or 1; default: 1)

    ``castro.do_special_tagging`` = 1 can be used, for example, to
    calculate the bounce time in a core collapse simulation; the
    bounce time is defined as the first time at which the maximum
    density in the domain exceeds a user-specified value. This time
    can then be printed into a special file as a useful diagnostic.

As an example::

    amr.grid_log = grdlog
    amr.run_log = runlog

Every time the code regrids it prints a list of grids at all relevant
levels. Here the code will write these grids lists into the file
``grdlog``. Additionally, every time step the code prints certain
statements to the screen (if ``amr.v`` = 1), such as::

    STEP = 1 TIME = 1.91717746 DT = 1.91717746
    PLOTFILE: file = plt00001

The ``run_log`` option will output these statements into
*runlog* as well.

Terser output can be obtained via::

    amr.run_log_terse = runlogterse

This file, ``runlogterse`` differs from ``runlog``, in that it
only contains lines of the form::

    10  0.2  0.005

in which “10” is the number of steps taken, “0.2” is the
simulation time, and “0.005” is the level-0 time step. This file
can be plotted very easily to monitor the time step.



Integral Diagnostics
--------------------

.. index:: castro.sum_interval, integral diagnostics

Castro can calculate integrals of quantities on the grid and other
global quantities and output them to both the screen and to a runtime
file at regular intervals.  By default, this capability is off.  To
enable it, one of the following runtime parameters can be set:

  * ``castro.sum_interval``: if :math:`> 0`, how often (in level-0 time
    steps) to compute and print integral quantities (Integer; default: -1)

    The integral quantities include total mass, momentum and energy in
    the domain every ``castro.sum_interval`` level-0 steps.  The print
    statements have the form::

           TIME= 1.91717746 MASS= 1.792410279e+34

    for example.

  * ``castro.sum_per``: how often in simulation time to output
    integral quantities (this is used as an alternate to
    ``castro.sum_interval``).

By default, 4 output files are created:

  * ``amr_diag.out`` : This includes timestep information, in the
    following columns:

    * timestep
    * time
    * dt
    * finest level
    * coarse timestep walltime

  * ``gravity_diag.out`` : For problems with Poisson gravity, this
    includes the gravitational wave amplitudes

  * ``grid_diag.out`` : This includes integrals of the state data:

    * time
    * mass
    * x-, y-, and z-momentum
    * x-, y-, and z-angular momentum
    * kinetic energy
    * internal energy
    * kinetic + internal energy
    * gravitational potential energy
    * total energy (including gravitational potential energy)

  * ``species_diag.out`` : This contains the mass of each of the nuclear species on the grid.

Some problems have custom versions of the diagnostics with additional information.


.. _sec:parallel_io:

Parallel I/O
------------

Both checkpoint files and plotfiles are really directories containing
subdirectories: one subdirectory for each level of the AMR hierarchy.
The fundamental data structure we read/write to disk is a ``MultiFab``,
which is made up of multiple FAB’s, one FAB per grid. Multiple
``MultiFab`` s may be written to each directory in a checkpoint file.
``MultiFab`` s of course are shared across CPUs; a single ``MultiFab`` may be
shared across thousands of CPUs. Each CPU writes the part of the
``MultiFab`` that it owns to disk, but they don’t each write to their own
distinct file. Instead each MultiFab is written to a runtime
configurable number of files :math:`N` (:math:`N` can be set in the inputs file as the
parameter ``amr.checkpoint_nfiles`` and ``amr.plot_nfiles``; the
default is 64). That is to say, each ``MultiFab`` is written to disk
across at most :math:`N` files, plus a small amount of data that gets written
to a header file describing how the file is laid out in those :math:`N` files.

What happens is :math:`N` CPUs each opens a unique one of the :math:`N` files into
which the ``MultiFab`` is being written, seeks to the end, and writes
their data. The other CPUs are waiting at a barrier for those :math:`N`
writing CPUs to finish. This repeats for another :math:`N` CPUs until all the
data in the ``MultiFab`` is written to disk. All CPUs then pass some data
to CPU 0 which writes a header file describing how the ``MultiFab`` is
laid out on disk.

We also read ``MultiFabs`` from disk in a “chunky” manner, opening only :math:`N`
files for reading at a time. The number :math:`N`, when the ``MultiFab`` s were
written, does not have to match the number :math:`N` when the ``MultiFab`` s are
being read from disk. Nor does the number of CPUs running while
reading in the ``MultiFab`` need to match the number of CPUs running when
the ``MultiFab`` was written to disk.

Think of the number :math:`N` as the number of independent I/O pathways in
your underlying parallel filesystem. Of course a “real” parallel
filesytem should be able to handle any reasonable value of :math:`N`. The
value -1 forces :math:`N` to the number of CPUs on which you’re
running, which means that each CPU writes to a unique file, which can
create a very large number of files, which can lead to inode issues.
*********************
Checkpoint Embiggener
*********************

Within the Castro distribution, there is the capability to “grow” a
checkpoint file so that a calculation can be restarted in a larger
domain covered by grid cells a factor of two or four coarser than the
existing coarsest level. Instructions for how to do so are in the
``Castro/Util/ConvertCheckpoint/README`` file and are included here.
Upon restart the existing data in the checkpoint file will be used to
fill the region of the previous computational domain, and the new
regions will be filled by some other means, typically interpolation
from a 1D model file.

Star in Corner (``star_at_center`` = 0)
=======================================

In this section we consider the case where the star (or feature of interest)
is centered at the lower left corner of the domain, e.g. you are modeling only one
quarter of the star in 2D, or an octant of the star in 3D. Then you only want
to grow the domain in the “high side” directions (e.g., to the upper right).

Converting the Checkpoint File
------------------------------

Let’s say you have a checkpoint file, ``chk00100`` with 5 levels of
refinement and a (real) problem domain size :math:`P` and (integer)
domain size :math:`D` at level 0.  The inputs file that created this
might have contained::

   max_step = 100
   amr.max_level = 5
   amr.n_cell = D D
   geometry.prob_lo = 0 0
   geometry.prob_hi = P P
   amr.ref_ratio = 4 4 4 4

Now let’s suppose that you want to grow the domain by a factor of 8
and cover that new larger domain with a level that is a factor of 2
coarser than the existing level 0 grids.

#. First, set ``DIM =`` in the GNUmakefile, and type ``make`` in the 
   ``Util/ConvertCheckpoint/`` directory.  This will
   make an executable from the ``Embiggen.cpp`` code.

#. Run the embiggening code as follows::

    Embiggen2d.Linux.Intel.Intel.ex checkin=chk00100 checkout=newchk00050 ref_ratio=2 grown_factor=8 star_at_center=0

   (Your executable may have a slightly different name depending on the compilers you
   built it with.)

   This will create a new checkpoint directory, called
   ``newchk00050``, that represents a simulation with *one* additional
   level of refinement *coarser* than the previous level 0 grids by a
   factor of ``ref_ratio`` (in this case, 2).  The new domain will be a
   factor of ``grown_factor`` (in this case, 8) larger than the previous
   domain.

   .. note:: ``ref_ratio`` must be 2 or 4, because those are the only
      acceptable values in Castro.

   ``grown_factor`` can be any reasonable integer; but it’s only been
   tested with 2, 3, 4 and 8. It does not need to be a multiple of 2.

Restarting from a Grown Checkpoint File
---------------------------------------

You should now be able to restart your calculation using ``newchk00050``.
Your inputs file should now contain lines like::

   max_step = 51
   amr.restart = newchk00050
   amr.max_level = 6
   amr.n_cell = 4D 4D
   geometry.prob_lo = 0 0
   geometry.prob_hi = 8P 8P
   castro.grown_factor = 8
   castro.star_at_center = 0
   amr.ref_ratio = 2 4 4 4 4

Important:

 * Unlike earlier, you may now set ``amr.max_level`` to be at most one
   greater than before, but you need not set it that high. For
   example, you could set ``amr.max_level`` the same as before and you
   would lose data at the finest refinement level. You may not set
   ``amr.max_level`` = 0, however, because we have no data at the new
   level 0 until we average down from the new level 1 after the
   restart.

 * You must set ``amr.n_cell`` = (``grown_factor`` / ``ref_ratio``)
   :math:`\times` (the previous value of ``amr.n_cell``). In this case
   ``amr.n_cell`` = (8/2)*D = 4D.

 * You must set ``amr.prob_hi`` to be a factor of ``grown_factor``
   greater than the previous value of ``amr.prob_hi``.

 * You must insert the value of ``ref_ratio`` used in the Embiggen
   call as the first value in the list of ``amr.ref_ratio``, since
   that will now be the refinement ratio between the new level 0 and
   the new level 1.

 * You must set ``castro.grown_factor`` in your inputs file equal to
   the value of ``grown_factor`` you used when you called Embiggen*ex
   so that Castro knows how big the original domain was.

 * Note that if you have run 100 steps at the original level 0, that
   would be equivalent to 50 steps at the new level 0 because you
   coarsened by a factor of 2.  Thus once you re-start from the new
   checkpoint directory, the next step will be 51, not 101. Make sure
   to keep track of your plotfiles accordingly.

 * Don’t forget to adjust ``max_denerr_lev`` and comparable variables
   to control the number of fine levels you now want. If you want to
   have 6 levels of refinement after restart, then make sure
   ``max_denerr_lev``, etc, are set high enough. If you only want to have
   5 levels of refinement (where the new level 5 would now be a factor
   of ``ref_ratio`` coarser than the previous level 5), make sure to
   adjust ``max_denerr_lev`` accordingly as well.

Star at Center of Domain (``star_at_center`` = 1)
=================================================

Now let’s assume that the star (or feature of interest) is centered at
the center of the domain in 2D or 3D Cartesian coordinates. We will
later consider the case of 2D cylindrical (r-z) coordinates in which
the star is centered at the left midpoint.

.. _converting-the-checkpoint-file-1:

Converting the Checkpoint File
------------------------------

Suppose that you want to grow the domain by a factor of 2 and cover
that new larger domain with a level that is a factor of 2 coarser than
the existing level 0 grids.

After you build the Embiggen executable, you type::

  Embiggen2d.Linux.Intel.Intel.ex checkin=chk00100 checkout=newchk00050 ref_ratio=2 grown_factor=2 star_at_center=1

Note that

-  ``ref_ratio`` must still be 2 or 4

-  ``grown_factor`` can only be 2 or 3 in this case.

.. _restarting-from-a-grown-checkpoint-file-1:

Restarting from a Grown Checkpoint File
---------------------------------------

Your inputs file for restarting would now look like::

   max_step = 51
   amr.restart = newchk00050
   amr.max_level = 6
   amr.n_cell = D D
   geometry.prob_lo = -P/2 -P/2
   geometry.prob_hi = 3P/2 3P/2
   castro.grown_factor = 2
   castro.star_at_center = 1
   amr.ref_ratio = 2 4 4 4 4

Cylindrical Coordinates
-----------------------

In the case of 2D cylindrical (r-z) coordinates in which the star is
centered at the left edge but vertical midpoint of the domain, the
embiggening procedure is the same as above (with ``star_at_center`` =
1) but the inputs file for restart is slightly different in that
``geometry.prob_lo`` is modified in the z- but not the r-direction. If
we consider the original inputs file to look like::

   max_step = 100
   amr.max_level = 6
   amr.n_cell = D 2D
   geometry.prob_lo = 0 0
   geometry.prob_hi = P 2P
   amr.ref_ratio = 4 4 4 4

then an inputs file for restart would look like::

   amr.restart = newchk00050
   amr.max_level = 6
   amr.n_cell = D 2D
   geometry.prob_lo = 0 -P
   geometry.prob_hi = 2P 3P
   castro.grown_factor = 2
   castro.star_at_center = 1
   amr.ref_ratio = 2 4 4 4 4


Some results:

.. figure:: corner.png

   Data from checkpoint file before and after the domain has been
   coarsened and grown. This case uses ``star_at_center`` = 0 and
   ``ref_ratio`` = 2. The first grown example has ``grown_factor`` =
   2, the second has ``grown_factor`` = 3. In all figures the level 0
   grids are shown in white, the level 1 grids in red, the level 2
   grids in yellow, and in the grown figures, the level 3 grids are in
   pink.

.. figure:: center.png

   Data from checkpoint file before and after the domain has been
   coarsened and grown. This case uses ``star_at_center`` = 1 and
   ``ref_ratio`` = 2. The first grown example has ``grown_factor`` =
   2, the second has ``grown_factor`` = 3. In all figures the level 0
   grids are shown in white, the level 1 grids in red, the level 2
   grids in yellow, and in the grown figures, the level 3 grids are in
   pink.

************
Microphysics
************

Equation of State
=================

Standard Castro EOSes
---------------------

Castro is written in a modular fashion so that the EOS and network
burning routines can be supplied by the user. However, for the
examples presented later we use several EOS and network routines
that come with the Microphysics distribution.

Castro relies on routines to calculate the equation of state (EOS)
of a fluid, as well as a species network to define the components of
the fluid. The network optionally has the ability to do nuclear burning,
but for this section its main purpose is in defining the species so that
the EOS can calculate fluid properties that depend on composition, such
as electron fraction.

Most of the standard problem setups in Castro (such as the Sedov blast wave)
use the ``gamma_law`` EOS. This represents a gamma law gas, with equation of state:

.. math:: p = (\gamma - 1) \rho e.

The gas is currently assumed to be monatomic and ideal.

Runtime Parameters
------------------

When inverting the EOS (e.g. by using ``eos_input_re``), an initial guess for
the temperature is required. This guess is provided by the runtime parameter
``castro.T_guess``, and should be set to a sensible value for each problem
(it will vary depending on which EOS is used).

EOS Interfaces and Parameters
-----------------------------

.. index:: eos_t

Each EOS should have two main routines by which it interfaces to the
rest of Castro. At the beginning of the simulation, ``eos_init``
will perform any initialization steps and save EOS variables (mainly
``smallt``, the temperature floor, and ``smalld``, the
density floor). Then, whenever you want to call the EOS, use::

 eos (eos_input, eos_state)

The first argument specifies the inputs to the EOS. The options
that are currently available are stored in Microphysics in
``interfaces/eos_type.H``, and are always a combination of two
thermodynamic quantities. For example, ``eos_input_rt`` means
that we call the EOS with :math:`\rho` (density) and :math:`T` (temperature)
and we expect the EOS to return the associated thermodynamic
quantities such as internal energy :math:`e` and entropy :math:`s`.

We note that for real (non-analytic) equations of state
in which :math:`\rho`, :math:`T` and species are the independent variables, such
as the Helmholtz EOS, ``eos_input_rt`` directly calls the EOS
and obtains the other thermodynamic variables. But for other inputs,
e.g. ``eos_input_re``, a Newton-Raphson iteration is performed
to find the density or temperature that corresponds to the given
input.

The eos_state variable is a C struct, ``eos_t``. It stores a complete
set of thermodynamic
variables. When calling the EOS, you should first fill the variables
that are the inputs, for example with

::

      eos_t eos_state;
      ...
      eos_state.rho = state(i,j,k,URHO);
      eos_state.T   = state(i,j,k,UTEMP);
      eos_state.e   = state(i,j,k,UEINT) / state(i,j,k,URHO);
      for (int n = 0; n < NumSpec; ++n) {
          eos_state.xn[n] = state(i,j,k,UFS+n) / state(i,j,k,URHO);
      }
      for (int n = 0; n < NumAux; ++n) {
          eos_state.aux[n] = state(i,j,k,UFX+n) / state(i,j,k,URHO);
      }

Whenever the ``eos_state`` type is initialized, the thermodynamic
state variables are filled with unphysical numbers. If you do not
input the correct arguments to match your input quantities, the EOS
will call an error. This means that it is good practice to fill the
quantities that will be iterated over with an initial guess. Indeed,
this initial guess is typically required for equations of state that
iterate over this variable, as the values they are initialized with
will likely not converge. Usually a prior value of the temperature or
density suffices if it’s available, but if not then use ``T_guess`` or
``small_dens``.

The `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository is the collection of microphysics routines that are compatible with the
AMReX-Astro codes. We refer you to the documentation in that repository for how to set it up
and for information on the equations of state provided. That documentation
also goes into more detail about the details of the EOS code, in case you are interested in
how it works (and in case you want to develop your own EOS).

Required Thermodynamics Quantities
----------------------------------

Three input quantities are required of any EOS:

-  ``eos_input_re``: :math:`\rho`, :math:`e`, and :math:`X_k` are input

-  ``eos_input_rt``: :math:`\rho`, :math:`T`, and :math:`X_k` are input

-  ``eos_input_rp``: :math:`\rho`, :math:`P`, and :math:`X_k` are input

The ``eos_t`` derived type holds a large number of thermodynamics
quantities, but not all of these are needed for basic
Castro operation. The main quantities that any EOS in any mode needs to
supply, if they are not input, are:

-  ``eos_state.T``: the temperature

-  ``eos_state.p``: total pressure

-  ``eos_state.e``: the specific energy

-  ``eos_state.gam1``: the first adiabatic index,
   :math:`\Gamma_1 = d\log P / d\log \rho |_s`

Additionally the ``eos_input_re`` mode also needs to supply:

-  ``eos_state.cs``: the adiabatic sound speed

-  ``eos_state.dpdr_e``: the derivative, :math:`\partial p/\partial \rho |_e`
   — note that the specific internal energy, :math:`e`
   is held constant here.

-  ``eos_state.dpde``: the derivative, :math:`\partial p / \partial e |_\rho`

For radiation hydro, the ``eos_input_rt`` model needs to supply:

-  ``eos_state.cv``: the specific heat capacity.

Other quantities (e.g., entropy) might be needed for the derived
variables that are optional output into the plotfiles.


Composition derivatives
-----------------------

.. index:: eos_xderivs_t

A separate type, ``eos_xderivs_t`` provides access to derivatives with respect to mass fraction.

-  ``eos_xderivs.dhdX[NumSpec]``: the derivative of the
   specific enthalpy with respect to mass fraction at constant
   :math:`T` and :math:`p`:

   .. math:: \xi_k = e_{X_k} + \frac{1}{p_\rho} \left (\frac{p}{\rho^2} - e_\rho \right ) p_{X_k}

-  ``eos_xderivs.dpdX[NumSpec]``: the derivative of the pressure with respect to mass fraction:

   .. math::

      \begin{align}
      p_{X_k} &= \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
                \frac{\partial \bar{A}}{\partial X_k} +
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
                \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &= -\frac{\bar{A}^2}{A_k}
                \left .\frac{\partial p}{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
                \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z} \right )
                \left . \frac{\partial p}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
      \end{align}

-  ``eos_xderivs.dedX[NumSpec]``: the derivative of the specific internal energy with respect to mass fraction:

   .. math::

      \begin{align}
      e_{X_k} &= \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}}
              \frac{\partial \bar{A}}{\partial X_k} +
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
              \frac{\partial \bar{Z}}{\partial X_k} \nonumber \\
              &= -\frac{\bar{A}^2}{A_k}
              \left . \frac{\partial e }{\partial \bar{A}} \right |_{\rho, T, \bar{Z}} +
              \frac{\bar{A}}{A_k} \left (Z_k - \bar{Z}\right )
              \left .\frac{\partial e}{\partial \bar{Z}} \right |_{\rho, T, \bar{A}}
      \end{align}

(see :cite:`maestro:III`, Appendix A).


Nuclear Network
===============

.. index:: burn_t

The nuclear network serves two purposes: it defines the fluid components used
in both the equation of state and the hydrodynamics, and it evolves those
components through a nuclear burning step. Castro comes with a ``general_null``
network (which lives in the ``networks/`` directory). This is a bare interface for a
nuclear reaction network. No reactions are enabled, and no auxiliary variables
are accepted.  It contains several sets of isotopes; for example,
``networks/general_null/triple_alpha_plus_o.net`` would describe the
isotopes needed to represent the triple-\ :math:`\alpha` reaction converting
helium into carbon, as well as oxygen and iron.

The main interface file, ``network.f90``, is a wrapper function. The
actual network details are defined in ``actual_network.f90``, a
file which is automatically generated in your work directory when you compile.
It supplies the number and names of species and auxiliary variables, as
well as other initializing data, such as their mass numbers, proton numbers,
and the binding energies.

The burning front-end interface, ``networks/burner.f90``, accepts a different
derived type called the ``burn_t`` type. Like the ``eos_t``, it has entries
for the basic thermodynamic quantities:

::

      use burn_type_module
      ...
      type (burn_t) :: burn_state
      ...
      burn_state % rho = state(i,j,k,URHO)
      burn_state % T   = state(i,j,k,UTEMP)
      burn_state % e   = state(i,j,k,UEINT) / state(i,j,k,URHO)
      burn_state % xn  = state(i,j,k,UFS:UFS+nspec-1) / state(i,j,k,URHO)

It takes in an input ``burn_t`` and returns an output ``burn_t`` after
the burning has completed. The nuclear energy release can be computed by
taking the difference of ``burn_state_out % e`` and
``burn_state_in % e``. The species change can be computed analogously.
In normal operation in Castro  the integration occurs over a time interval
of :math:`\Delta t/2`, where :math:`\Delta t` is the hydrodynamics timestep.

If you are interested in using actual nuclear burning networks,
you should download the `Microphysics <https://github.com/starkiller-astro/Microphysics>`__
repository. This is a collection of microphysics routines that are compatible with the
AMReX Astro codes. We refer you to the documentation in that repository for how to set it up
and for information on the networks provided. That documentation
also goes into more detail about the details of the network code, in case you are interested in
how it works (and in case you want to develop your own network).


Controlling burning
-------------------

There are a number of reactions-related parameters that can be set at runtime
in the inputs file. Reactions are enabled by setting::

    castro.do_react = 1

(Note: turning reactions off for problems where they're not required can help improve
the efficiency).

It is possible to set the maximum and minimum temperature and density for allowing
reactions to occur in a zone using the parameters ``castro.react_T_min``,
``castro.react_T_max``, ``castro.react_rho_min`` and ``castro.react_rho_max``.

*********
Debugging
*********

There are several methods we typically use to debug issues in Castro.
Some descriptions are given below.

Using the compiler and other tools
==================================

Compiler checks
---------------

Recompile the code in debug-mode by setting::

   DEBUG := TRUE

in the ``GNUmakefile`` (or invoking ``make`` as ``make DEBUG=TRUE``).
This will create an executable with bounds checking and other compiler
options enabled.  Running this sometimes will show issues.


Check for invalid floating point exceptions
-------------------------------------------

Compiling in debug mode also initializes uninitialized variables to
NaN.  For optimized code, the same can be done by setting::

   TEST := TRUE

in the ``GNUmakefile``.  To capture the NaNs, use the runtime parameter::

   amrex.fpe_trap_invalid=1

If the code crashes, it will produce one or more ``Backtrace.*``
files.  Looking over these files should pinpoint where the FPE
occurred allowing you to do further debugging.

You can output more information into the ``Backtrace.*`` files by
pushing it to the backtrace stack as described here:
:ref:`debugging_backtrace`.

Make sure your runtime options are valid
----------------------------------------

Castro can validate the runtime options specified in the inputs file
by running with ``castro.abort_on_invalid_params = 1``.


Valgrind
--------

We frequently run Castro with valgrind to find illegal memory
accesses.  The valgrind documentation can give details on how to use
it.


Clang-tidy
----------

.. index:: clang-tidy

We run `clang-tidy <https://clang.llvm.org/extra/clang-tidy/>`_ on all
pull requests using a `GitHub action
<https://github.com/AMReX-Astro/cpp-linter-action>`_. ``clang-tidy``
analyzes the source code, produces warnings for potential bugs and
offers suggestions for performance improvements.

It can also be run locally. This requires the ``clang-tidy`` and
``bear`` packages, and the python script ``run-clang-tidy.py`` (which
can be downloaded from `here
<https://github.com/AMReX-Astro/cpp-linter-action/blob/main/run-clang-tidy.py>`_). The
analysis is performed by first compiling a problem using the ``bear``
package, then running the python script to analyze the source
files. From within a problem directory, run

.. code-block:: bash

    bear make COMP=llvm -j 20 USE_OMP=FALSE USE_MPI=FALSE DEBUG=TRUE 

    python3 run-clang-tidy.py -header-filter='Castro' -ignore-files='amrex|Microphysics' -j 20 > clang-tidy-report.txt

The compiler flags can be modified to suit the problem to be analyzed, but the ``DEBUG`` flag must be set to ``TRUE``. The ``header-filter`` option for the python script tells the script to only analyze header files containing the given regex pattern, and the ``ignore-files`` flag tells it to ignore any source files containing the given regex pattern. The ``-j`` option tells the script to run a given number of processes in parallel. The output is then redirected to a text file. 

Thread sanitizer
----------------



Instrumenting the Code
======================

Checking for NaNs
-----------------

In the C++ code, you can check whether a FAB contains NaNs using
the ``contains_nan()`` method:

.. code-block:: c++

   for (MFIter mfi(S, True); mfi.isValid(); ++mfi) {

     const Box& bx = mf.tilebox()

     // do some operations on S[mfi]

     if (S[mfi].contains_nan()) {
       amrex::Abort("S has NaNs")
     }
   }

There are other versions of ``contains_nan()`` that can take a Box
to operate over.



Physics issues
==============


***************
Software Design
***************

Code structure
==============

Castro is built upon the AMReX C++ framework. This provides
high-level classes for managing an adaptive mesh refinement
simulation, including the core data structures we will deal with. A
key design pattern in AMReX is that the overall memory management
and parallelization is done in the C++ layer, while the heavy
computational work is done in Fortran kernels. AMReX provides
convenient data structures that allow for this workflow—high level
objects in C++ that communicate with Fortran through pointers to
data regions that appear as multidimensional arrays.

Castro uses a structured-grid approach to hydrodynamics. We work
with square/cubic zones that hold state variables (density, momentum,
etc.) and compute the fluxes of these quantities through the
interfaces of the zones (this is a finite-volume approach).
Parallelization is achieved by domain decomposition. We divide our
domain into many smaller boxes, and distributed these across
processors. When information is needed at the boundaries of the
boxes, messages are exchanged and this data is held in a perimeter of
*ghost cells*. AMReX manages this decompostion and
communication for us. Additionally, AMReX implements adaptive mesh
refinement. In addition to the coarse decomposition of our domain
into zones and boxes, we can refine rectangular regions by adding
finer-gridded boxes on top of the coarser grid. We call the
collection of boxes at the same resolution a *level*.

Castro uses a hybrid MPI + OpenMP approach to parallelism. MPI is
at used to communicate across nodes on a computer and OpenMP is used
within a node, to loop over subregions of a box with different
threads.

The code structure in the Castro/ directory reflects the
division between C++ and Fortran.

-  ``Diagnostics/``: various analysis routines for specific problems

-  ``Docs/``: you’re reading this now!

-  ``Exec/``: various problem implementations, sorted by category:

   -  ``gravity_tests/``: test problems that primarily exercise the gravity solver

   -  ``hydro_tests/``: test problems of the hydrodynamics (with or without reactions)

   -  ``radiation_tests/``: test problems that primarily exercise the radiation hydrodynamics solver

   -  ``reacting_tests/``: test problems that primarily exercise the reactions (and hydro + reaction coupling)

   -  ``scf_tests/``: problem setups that use the self-consistent field initialization

   -  ``science/``: problem setups that were used for scientific investigations

   -  ``unit_tests/``: test problems that exercise primarily a single module

-  ``external/``: if you are using git submodules, the Microphysics and AMReX git
   submodules will be in this directory.

-  ``paper/``: the JOSS paper source

-  ``Source/``: source code. In this main directory is all of the
   code. Sources are mixed C++ and Fortran and are organized by topic
   as:

   -  ``diffusion/`` : thermal diffusion code

   -  ``driver/`` : the main driver, I/O, runtime parameter support

   -  ``gravity/`` : self-gravity code

   -  ``hydro/`` : the compressible hydrodynamics code

   -  ``mhd/`` : the MHD solver code

   -  ``particles/`` : support for particles

   -  ``problems/`` : template code for implementing a problem

   -  ``radiation/`` : the implicit radiation solve code

   -  ``reactions/`` : nuclear reaction code

   -  ``rotation/`` : rotating code

   -  ``scf/`` : the self-consistent field initialization support

   -  ``sdc/``: code specified for the true SDC method

   -  ``sources/`` : hydrodynamics source terms support

-  ``Util/``: a catch-all for additional things you may need

   -  ``ConvertCheckpoint/``: a tool to convert a checkpoint file to
      a larger domain

   -  ...

Major data structures
=====================

The following data structures are the most commonly encountered when
working in the C++ portions of Castro. This are all
AMReX data-structures / classes.

``Amr``
-------

This is the main class that drives the whole simulation. This is
the highest level in Castro.

``AmrLevel`` and Castro classes
-------------------------------

An ``AmrLevel`` is a virtual base class provided by AMReX that
stores all the state data on a single level in the AMR hierarchy and
understands how to advance that data in time.

The most important data managed by the ``AmrLevel`` is an array of
``StateData``, which holds the fluid quantities, etc., in the boxes
that together make up the level.

The ``Castro`` class is derived from the ``AmrLevel``. It provides
the Castro-specific routines to evolve our system of equations. Like
the ``AmrLevel``, there is one ``Castro`` object for each level in the
AMR hierarchry.

A lot of the member data in the ``Castro`` class are static member
variables—this means that they are shared across all instances of
the class. So, in this case, every level will have the same data.
This is done, in particular, for the values of the runtime parameters,
but also for the ``Gravity``, ``Diffusion``, and ``Radiation``
objects. This means that those objects cover all levels and are the
same object in each instantiation of the ``Castro`` class.

Floating point data
-------------------

Floating point data in the C++ AMReX framework is declared as
``Real``. This is typedef to either ``float`` or ``double`` depending
on the make variable ``PRECISION``.

The corresponding type for Fortran is provided by the
``amrex_fort_module`` as ``c_real``. We typically rename
this to rt when using it. An example of a declaration of a
parameter is::

      use amrex_fort_module, only : rt => amrex_real
      real(rt) :: tol = 1.0e-10_rt

The ``amrex_constants_module`` provides common constants that can
be used in the code, like ``ZERO``, ``THIRD``, ``ONE``, etc.

.. note :: single precision support in Castro is not yet complete. In
   particular, a lot of the supporting microphysics has not been updated.

``Box`` and ``FArrayBox``
-------------------------

A ``Box`` is simply a rectangular region in space. It does not hold
data. In AMReX, an AMR level has a global index space, with
:math:`(0,0,0)` being the lower left corner of the domain at that level, and
:math:`(N_x-1, N_y-1, N_z-1)` being the upper right corner of the domain
(for a domain of :math:`N_x \times N_y \times N_z` zones). The location of
any ``Box`` at a level can be uniquely specified with respect to this
global index space by giving the index of its lower-left and
upper-right corners. :numref:`fig:soft:indexspace` shows an
example of three boxes at the same level of refinement.

AMReX provides other data structures that collect Boxes together,
most importantly the ``BoxArray``. We generally do not use these
directly, with the exception of the ``BoxArray`` ``grids``,
which is defined as part of the ``AmrLevel`` class that ``Castro``
inherits. ``grids`` is used when building new ``MultiFabs`` to give
the layout of the boxes at the current level.

.. _fig:soft:indexspace:
.. figure:: index_grid2.png
   :width: 4in

   Three boxes that comprise a single level. At this
   resolution, the domain is 20 :math:`\times` 18 zones. Note that the
   indexing in AMReX starts with :math:`0`.

A ``FArrayBox`` or *FAB*, for *Fortran array box* is a data
structure that contains a ``Box`` locating it in space, as well as a
pointer to a data buffer. The real floating point data are stored as
one-dimensional arrays in ``FArrayBox`` es. The associated ``Box`` can be
used to reshape the 1D array into multi-dimensional arrays to be used
by Fortran subroutines. The key part of the C++ AMReX data
structures is that this data buffer can be sent to Fortran, where it
will appear as a DIM+1 dimensional array (DIM space + 1
component).

.. note:: Castro is complied for a specific dimensionality.

``MultiFab``
------------

At the highest abstraction level, we have the ``MultiFab`` (mulitple
FArrayBoxes). A ``MultiFab`` contains an array of ``Box`` es, including
boxes owned by other processors for the purpose of communication,
an array of MPI ranks specifying which MPI processor owns each ``Box``,
and an array of pointers to ``FArrayBoxes`` owned by this MPI
processor. 

.. note:: a ``MultiFab`` is a collection of the boxes that together
   make up a single level of data in the AMR hierarchy.

A ``MultiFab`` can have multiple components (like density, temperature,
...) as well as a perimeter of ghost cells to exchange data with
neighbors or implement boundary conditions (this is all reflected in
the underlying ``FArrayBox``).

Parallelization in AMReX is done by distributing the FABs across
processors. Each processor knows which FABs are local to it. To loop
over all the boxes local to a processor, an ``MFIter`` is used (more
on this below).

High-level operations exist on ``MultiFab`` s to add, subtract, multiply,
etc., them together or with scalars, so you don’t need to write out
loops over the data directly.

In Castro, ``MultiFab`` s are one of the main data structures you will
interact with in the C++ portions of the code.

.. _soft:sec:statedata:

``StateData``
-------------

``StateData`` is a class that essentially holds a pair of
``MultiFab`` s: one at the old time and one at the new
time. AMReX knows how to interpolate in time between these states to
get data at any intermediate point in time. The main data that we care
about in Castro (the fluid state, gravitational potential, etc.) will
be stored as ``StateData``. Essentially, data is made StateData in
Castro if we need it to be stored in checkpoints / plotfiles, and/or
we want it to be automatically interpolated when we refine.

An ``AmrLevel`` stores an array of ``StateData`` (in a C++ array
called ``state``). We index this array using integer keys (defined
via an enum in ``Castro.H``). The state data is registered
with AMReX in ``Castro_setup.cpp``.

Note that each of the different ``StateData`` carried in the state
array can have different numbers of components, ghost cells, boundary
conditions, etc. This is the main reason we separate all this data
into separate StateData objects collected together in an indexable
array.

The current ``StateData`` names Castro carries are:

-  ``State_Type`` : this is the ``NUM_STATE`` hydrodynamics
   components that make up the conserved hydrodynamics state (usually
   referred to as :math:`\Ub` in these notes. But note that this does
   not include the radiation energy density.

   We access this data using an AMReX ``Array4`` type which is
   of the form ``data(i,j,k,n)``, where ``n`` is the component.
   The integer keys used to index the components are defined
   in ``Source/driver/_variables`` (e.g., ``URHO``, ``UMX``,
   ``UMY``, ...)

   .. note:: regardless of dimensionality, we always carry around all
      three velocity components. The “out-of-plane” components will
      simply be advected, but we will allow rotation (in particular,
      the Coriolis force) to affect them.

   ``State_Type`` ``MultiFab`` s have no ghost cells by default for
   pure hydro and a single ghost cell by default when ``RADIATION``
   is enabled. There is an option to force them to have ghost cells by
   setting the parameter ``castro.state_nghost`` at runtime.

   Note that the prediction of the hydrodynamic state to the interface
   will require 4 ghost cells. This accomodated by creating a separate
   MultiFab, ``Sborder`` that lives at the old-time level and
   has the necessary ghost cells. We will describe this more later.

-  ``Rad_Type`` : this stores the radiation energy density,
   commonly denoted :math:`E_r` in these notes. It has ``nGroups``
   components—the number of energy groups used in the multigroup
   radiation hydrodynamics approximation.

-  ``PhiGrav_Type`` : this is simply the gravitational
   potential, usually denoted :math:`\Phi` in these notes.

-  ``Gravity_Type`` : this is the gravitational
   acceleration. There are always 3 components, regardless of the
   dimensionality (consistent with our choice of always carrying all 3
   velocity components).

-  ``Source_Type`` : this holds the time-rate of change of
   the source terms, :math:`d\Sb/dt`, for each of the ``NUM_STATE``
   ``State_Type`` variables.


   .. note:: we do not make use of the old-time quantity here. In
      fact, we never allocate the ``FArrayBox`` s for the old-time in
      the ``Source_Type`` ``StateData``, so there is not wasted
      memory.

-  ``Reactions_Type`` : this holds the data for the nuclear
   reactions. It has ``NumSpec+2`` components: the species
   creation rates (usually denoted :math:`\omegadot_k` in these notes),
   the specific energy generation rate (:math:`\dot{e}_\mathrm{nuc}`),
   and its density (:math:`\rho \dot{e}_\mathrm{nuc}`).

   These are stored as ``StateData`` so we have access to the reaction terms
   outside of advance, both for diagnostics (like flame speed estimation)
   and for reaction timestep limiting (this in particular needs the
   data stored in checkpoints for continuity of timestepping upon restart).

- ``Mag_Type_x`` : this is defined for MHD and stores the
   face-centered (on x-faces) x-component of the magnetic field.

- ``Mag_Type_y`` : this is defined for MHD and stores the
   face-centered (on y-faces) y-component of the magnetic field.

- ``Mag_Type_z`` : this is defined for MHD and stores the
   face-centered (on z-faces) z-component of the magnetic field.

-  ``Simplified_SDC_React_Type`` : this is used with the SDC
   time-advancement algorithm. This stores the ``NQSRC`` terms
   that describe how the primitive variables change over the timestep
   due only to reactions. These are used when predicting the interface
   states of the primitive variables for the hydrodynamics portion of the
   algorithm.

We access the ``MultiFab`` s that carry the data of interest by interacting
with the ``StateData`` using one of these keys. For instance::

    MultiFab& S_new = get_new_data(State_Type);

gets a pointer to the ``MultiFab`` containing the hydrodynamics state data
at the new time.

``MFIter`` and interacting with Fortran
=======================================

The process of looping over boxes at a given level of refinement and
operating on their data is linked to how Castro achieves
thread-level parallelism. The OpenMP approach in Castro has evolved
considerably since the original paper was written, with the modern
approach, called *tiling*, gearing up to meet the demands of
many-core processors in the next-generation of supercomputers.

Full details of iterating over boxes and calling compute kernels
is given in the AMReX documentation here: https://amrex-codes.github.io/amrex/docs_html/Basics.html#mfiter-and-tiling


Practical Details in Working with Tiling
----------------------------------------

With tiling, the OpenMP is now all in C++, and not in Fortran for all
modules except reactions and ``initdata``.

It is the responsibility of the coder to make sure that the routines
within a tiled region are safe to use with OpenMP. In particular,
note that:

-  tile boxes are non-overlapping

-  the union of tile boxes completely cover the valid region of the
   fab

-  Consider working with a node-centered MultiFab, ``ugdnv``, and
   a cell-centered ``MultiFab`` ``s``:

   -  with ``mfi(s)``, the tiles are based on the cell-centered
      index space. If you have an :math:`8\times 8` box, then and 4 tiles,
      then your tiling boxes will range from :math:`0\rightarrow 3`,
      :math:`4\rightarrow 7`.

   -  with ``mfi(ugdnv)``, the tiles are based on nodal indices,
      so your tiling boxes will range from :math:`0\rightarrow 3`,
      :math:`4\rightarrow 8`.

-  When updating routines to work with tiling, we need to
   understand the distinction between the index-space of the entire box
   (which corresponds to the memory layout) and the index-space of the
   tile.

   -  In the C++ end, we pass (sometimes via the
      ``BL_TO_FORTRAN()`` macro) the ``loVect`` and ``hiVect`` of the
      entire box (including ghost cells). These are then used to
      allocate the array in Fortran as::

            double precision :: a(a_l1:a_h1, a_l2:a_h2, ...)

      When tiling is used, we do not want to loop as do ``a_l1``,
      ``a_h1``, but instead we need to loop over the tiling region. The
      indices of the tiling region need to be passed into the Fortran
      routine separately, and they come from the ``mfi.tilebox()``
      or ``mfi.growntilebox()`` statement.

   -  In Fortran, when initializing an array to 0, do so only
      over the tile region, not for the entire box. For a Fortran array
      a, this means we cannot do::

            a = 0.0
            a(:,:,:,:) = 0.0

      but instead must do::

            a(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3),:) = 0.0

      where ``lo()`` and ``hi()`` are the index-space for the tile box
      returned from ``mfi.tilebox()`` in C++ and passed into the Fortran
      routine.

   -  Look at ``r_old_s`` in ``Exec/gravity_tests/DustCollapse/probdata.f90`` as an
      example of how to declare a ``threadprivate`` variable—this is then used
      in ``sponge_nd.f90``.

Boundaries: ``FillPatch`` and ``FillPatchIterator``
===================================================

AMReX calls the act of filling ghost cells a *fillpatch*
operation. Boundaries between grids are of two types. The first we
call “fine-fine”, which is two grids at the same level. The second
type is "coarse-fine", which needs interpolation from the coarse grid
to fill the fine grid ghost cells. Both of these are part of the
fillpatch operation. Fine-fine fills are just a straight copy from
“valid regions” to ghost cells. Coarse-fine fills are enabled
because the ``StateData`` is not just arrays, they’re “State Data”,
which means that the data knows how to interpolate itself (in an
anthropomorphical sense). The type of interpolation to use is defined
in ``Castro_setup.cpp``—search for
``cell_cons_interp``, for example—that’s “cell conservative
interpolation”, i.e., the data is cell-based (as opposed to
node-based or edge-based) and the interpolation is such that the
average of the fine values created is equal to the coarse value from
which they came. (This wouldn’t be the case with straight linear
interpolation, for example.)

Additionally, since ``StateData`` has an old and new timelevel,
the fill patch operation can interpolate to an intermediate time.

Examples
--------

To illustrate the various ways we fill ghost cells and use the data,
let’s consider the following scenarios:

-  *You have state data that was defined with no ghost cells. You
   want to create a new* ``MultiFab`` *containing a copy of that data with*
   ``NGROW`` *ghost cells.*

   This is the case with ``Sborder`` —the ``MultiFab`` of the
   hydrodynamic state that we use to kick-off the hydrodynamics
   advance.

   ``Sborder`` is declared in ``Castro.H`` simply as:

   .. code:: c++

         Multifab Sborder;

   It is then allocated in ``Castro::initialize_do_advance()``

   .. code:: c++

         Sborder.define(grids, NUM_STATE, NUM_GROW, Fab_allocate);
         const Real prev_time = state[State_Type].prevTime();
         expand_state(Sborder, prev_time, NUM_GROW);

   Note in the call to ``.define()``, we tell AMReX to already
   allocate the data regions for the ``FArrayBox`` s that are part of
   ``Sborder``.

   The actually filling of the ghost cells is done by
   ``Castro::expand_state()``:

   .. code:: c++

         AmrLevel::FillPatch(*this, Sborder, NUM_GROW,
                             prev_time, State_Type, 0, NUM_STATE);

   Here, we are filling the ng ghost cells of ``MultiFab``
   ``Sborder`` at time prev_time. We are using the
   ``StateData`` that is part of the current ``Castro`` object that we
   are part of. Note: ``FillPatch`` takes an object reference as its
   first argument, which is the object that contains the relevant
   ``StateData`` —that is what the this pointer indicates.
   Finally, we are copying the ``State_Type`` data components 0 to
   ``NUM_STATE`` [1]_.

   The result of this operation is that ``Sborder`` will now have
   ``NUM_GROW`` ghost cells consistent with the ``State_Type``
   data at the old time-level.

-  *You have state data that was defined with* ``NGROW`` *ghost
   cells. You want to ensure that the ghost cells are filled
   (including any physical boundaries) with valid data.*

   This is very similar to the procedure shown above. The main
   difference is that for the ``MultiFab`` that will be the target
   of the ghost cell filling, we pass in a reference to the ``StateData`` itself.

   The main thing you need to be careful of here, is that you
   need to ensure that the the time you are at is consistent with
   the ``StateData`` ’s time. Here’s an example from the radiation
   portion of the code ``MGFLDRadSolver.cpp``:

   .. code:: c++

         Real time = castro->get_state_data(Rad_Type).curTime();
         MultiFab& S_new = castro->get_new_data(State_Type);

         AmrLevel::FillPatch(*castro, S_new, ngrow, time, State_Type,
                             0, S_new.nComp(), 0);

   In this example, ``S_new`` is a pointer to the new-time-level
   ``State_Type`` ``MultiFab``. So this operation will use the
   ``State_Type`` data to fill its own ghost cells. we fill the
   ``ngrow`` ghost cells of the new-time-level ``State_Type`` data,
   for all the components.

   Note that in this example, because the ``StateData`` lives in the
   ``castro`` object and we are working from the ``Radiation`` object,
   we need to make reference to the current ``castro`` object
   pointer. If this were all done within the ``castro`` object, then
   the pointer will simply be ``this``, as we saw above.

-  *You have a* ``MultiFab`` *with some derived quantity. You want to
   fill its ghost cells.*

   ``MultiFabs`` have a ``FillBoundary()`` method that will fill all
   the ghost cells between boxes at the same level. It will not fill
   ghost cells at coarse-fine boundaries or at physical boundaries.

-  *You want to loop over the FABs in state data, filling ghost cells
   along the way*

   This is the job of the ``FillPatchIterator``—this iterator is used
   to loop over the grids and fill ghostcells. A key thing to keep in
   mind about the ``FillPatchIterator`` is that you operate on a copy of
   the data—the data is disconnected from the original source. If you
   want to update the data in the source, you need to explicitly copy
   it back. Also note: ``FillPatchIterator`` takes a ``MultiFab``, but this is
   not filled—this is only used to get the grid layout. Finally, the
   way the ``FillPatchIterator`` is implemented is that all the
   communication is done first, and then the iterating over boxes
   commences.

   For example, the loop that calls ``CA_UMDRV`` (all the
   hydrodynamics integration stuff) starts with::

          for (FillPatchIterator fpi(*this, S_new, NUM_GROW,
                                     time, State_Type, strtComp, NUM_STATE);
                fpi.isValid(); ++fpi)
          {
            FArrayBox &state = fpi();
            Box bx(fpi.validbox());

            // work on the state FAB.  The interior (valid) cells will
            // live between bx.loVect() and bx.hiVect()
          }

   Here the ``FillPatchIterator`` is the thing that distributes the
   grids over processors and makes parallel “just work”. This fills the
   single patch ``fpi`` , which has ``NUM_GROW`` ghost cells,
   with data of type ``State_Type`` at time ``time``,
   starting with component strtComp and including a total of
   ``NUM_STATE`` components.

In general, one should never assume that ghostcells are valid, and
instead do a fill patch operation when in doubt. Sometimes we will
use a ``FillPatchIterator`` to fill the ghost cells into a ``MultiFab``
without an explict look. This is done as::

      FillPatchIterator fpi(*this,S_old,1,time,State_Type,0,NUM_STATE);
      MultiFab& state_old = fpi.get_mf();

In this operation, state_old points to the internal
``MultiFab`` in the ``FillPatchIterator``, by getting a reference to it as
``fpi.get_mf()``. This avoids a local copy.

Note that in the examples above, we see that only ``StateData`` can fill
physical boundaries (because these register how to fill the boundaries
when they are defined). There are some advanced operations in
AMReX that can get around this, but we do not use them in Castro.

.. _soft:phys_bcs:

Physical Boundaries
-------------------

.. index:: boundary conditions

Physical boundary conditions are specified by an integer index [2]_ in
the ``inputs`` file, using the ``castro.lo_bc`` and ``castro.hi_bc`` runtime
parameters. The generally supported boundary conditions are, their
corresponding integer key, and the action they take for the normal
velocity, transverse velocity, and generic scalar are shown in 
:numref:`table:castro:bcs`.

The definition of the specific actions are:

-  ``INT_DIR``: data taken from other grids or interpolated

-  ``EXT_DIR``: data specified on EDGE (FACE) of bndry

-  ``HOEXTRAP``: higher order extrapolation to EDGE of bndry

-  ``FOEXTRAP``: first order extrapolation from last cell in interior

-  ``REFLECT_EVEN``: :math:`F(-n) = F(n)` true reflection from interior cells

-  ``REFLECT_ODD``: :math:`F(-n) = -F(n)` true reflection from interior cells

The actual registration of a boundary condition action to a particular
variable is done in ``Castro_setup.cpp``. At the top we define arrays
such as ``scalar_bc``, ``norm_vel_bc``, etc, which say which kind of
bc to use on which kind of physical boundary.  Boundary conditions are
set in functions like ``set_scalar_bc``, which uses the ``scalar_bc``
pre-defined arrays. We also specify the name of the Fortran routine
that is responsible for filling the data there (e.g., ``hypfill``).  These
routines are discussed more below.

If you want to specify a value at a function (like at an inflow
boundary), then you choose an *inflow* boundary at that face of
the domain. You then need to write the implementation code for this.
There is a centralized hydrostatic boundary condition that is implemented
this way—see :ref:`create:bcs`.

.. _table:castro:bcs:
.. table:: Physical boundary conditions supported in Castro.

   +-------------+-------------+-------------+--------------+--------------+
   | **name**    | **integer** | **normal    | **transverse | **scalars**  |
   |             |             | velocity**  | velocity**   |              |
   +=============+=============+=============+==============+==============+
   | interior    | 0           | INT_DIR     | INT_DIR      | INT_DIR      |
   +-------------+-------------+-------------+--------------+--------------+
   | inflow      | 1           | EXT_DIR     | EXT_DIR      | EXT_DIR      |
   +-------------+-------------+-------------+--------------+--------------+
   | outflow     | 2           | FOEXTRAP    | FOEXTRAP     | FOEXTRAP     |
   +-------------+-------------+-------------+--------------+--------------+
   | symmetry    | 3           | REFLECT_ODD | REFLECT_EVEN | REFLECT_EVEN |
   +-------------+-------------+-------------+--------------+--------------+
   | slipwall    | 4           | REFLECT_ODD | REFLECT_EVEN | REFLECT_EVEN |
   +-------------+-------------+-------------+--------------+--------------+
   | noslipwall  | 5           | REFLECT_ODD | REFLECT_EVEN | REFLECT_EVEN |
   +-------------+-------------+-------------+--------------+--------------+

``FluxRegister``
----------------

A ``FluxRegister`` holds face-centered data at the boundaries of a box.
It is composed of a set of ``MultiFab`` s (one for each face, so 6 for
3D). A ``FluxRegister`` stores fluxes at coarse-fine interfaces,
and isused for the flux-correction step.

Other AMReX Concepts
====================

There are a large number of classes that help define the structure of
the grids, metadata associate with the variables, etc. A good way to
get a sense of these is to look at the ``.H`` files in the
``amrex/Src/`` directory.

``Geometry`` class
------------------

There is a ``Geometry`` object, ``geom`` for each level as part of
the ``Castro`` object (this is inhereted through ``AmrLevel``).

``ParmParse`` class
-------------------

Error Estimators
----------------

``Gravity`` class
=================

There is a single ``Gravity`` object, ``gravity``, that is a
static class member of the ``Castro`` object. This means that all
levels refer to the same ``Gravity`` object.

Within the ``Gravity`` object, there are pointers to the ``Amr``
object (as ``parent``), and all of the ``AmrLevels`` (as a ``PArray``,
``LevelData``). The ``gravity`` object gets the geometry
information at each level through the parent ``Amr`` class.

The main job of the ``gravity`` object is to provide the potential
and gravitation acceleration for use in the hydrodynamic sources.
Depending on the approximation used for gravity, this could mean
calling the AMReX multigrid solvers to solve the Poisson equation.

Fortran Helper Modules
======================

There are a number of modules that make data available to the Fortran
side of Castroor perform other useful tasks.

``amrex_constants_module``
--------------------------

This provides double precision constants as Fortran parameters, like
``ZERO``, ``HALF``, and ``ONE``.


``fundamental_constants_module``
--------------------------------

This provides the CGS values of many physical constants.


``meth_params_module``
----------------------

This module provides the integer keys used to access the state arrays
for both the conserved variables (``URHO``, ``UMX``, :math:`\ldots`)
and primitive variables (``QRHO``, ``QU``, :math:`\ldots`), as well as
the number of scalar variables.

It also provides the values of most of the ``castro.*xxxx*``
runtime parameters.


.. [1]
   for clarity and continuity in this
   documentation, some of the variable names have been changed
   compared to the actual code

.. [2]
   the integer values are defined in ``BC_TYPES.H``

.. _ch:diffusion:

*********
Diffusion
*********


Thermal Diffusion
=================

Castro incorporates explicit thermal diffusion into the energy equation.
In terms of the specific internal energy, :math:`e`, this appears as:

.. math:: \rho \frac{De}{Dt} + p \nabla \cdot \ub = \nabla \cdot \kth \nabla T

where :math:`\kth` is the thermal conductivity, with units
:math:`\mathrm{erg~cm^{-1}~s^{-1}~K^{-1}}`.

To see the similarity to the thermal diffusion equation, consider the
special case of constant conductivity, :math:`\kth`, and density, and
assume an ideal gas, so :math:`e = c_v T`, where :math:`c_v` is the
specific heat at constant volume.  Finally, ignore hydrodynamics, so
:math:`\ub = 0`. This gives:

.. math:: \frac{\partial T}{\partial t} = D \nabla^2 T

where :math:`D \equiv \kth/(\rho c_v)`. Solving this equation
explicitly requires a timestep limiter of

.. math:: \Delta t_\mathrm{diff} \le \frac{1}{2} \frac{\Delta x^2}{D}

(this is implemented in ``ca_estdt_temp_diffusion`` in
``Castro/Source/driver/timestep.F90``).

Support for diffusion must be compiled into the code by setting
``USE_DIFFUSION = TRUE`` in your ``GNUmakefile``. It is treated
explicitly, by constructing the contribution to the evolution as a
source term. This is time-centered to achieve second-order accuracy
in time.

The following parameter affects diffusion:

-  ``castro.diffuse_temp``: enable thermal diffusion (0 or 1; default 0)

A pure diffusion problem (with no hydrodynamics) can be run by setting::

    castro.diffuse_temp = 1
    castro.do_hydro = 0

To complete the setup, a thermal conductivity must be specified. The
interface for the conductivity is::

      subroutine actual_conductivity(eos_state)

        type (eos_t), intent(inout) :: eos_state

The density, temperature, and mass fractions come in through the
``eos_state`` type. An EOS call is done in Castro just before the call to
``thermal_conductivity``, so you can assume that the entire state is
consistent.  The conductivity is filled in ``eos_state % conductivity``.

.. index:: CONDUCTIVITY_DIR

There are two conductivity routines provided with Castro by default:

-  ``constant`` : A simple constant thermal conductivity. This can be
   selected by setting::

       CONDUCTIVITY_DIR := constant

   in your ``GNUmakefile``. To set the value of the conductivity (e.g., to
   :math:`100`), you add to your input file::

       conductivity.const_conductivity = 100.0

-  ``constant_opacity`` : A simple constant opacity. This is
   converted to an opacity as:

   .. math:: \kth = \frac{16 \sigma_B T^3}{3 \kappa_\mathrm{const} \rho}

   where :math:`\kappa_\mathrm{const}` is the opacity, with units :math:`\mathrm{cm^2~g^{-1}}`.
   This is selected by setting::

       CONDUCTIVITY_DIR := constant_opacity

   in your ``GNUmakefile``. To set the value of the opacity, e.g., to
   0.2 (for electron scattering), set::

       conductivity.const_opacity = 0.2

   in the inputs file.

.. index:: castro.diffusion_cutoff_density, castro.diffusion_cutoff_density_hi

The diffusion approximation breaks down at the surface of stars,
where the density rapidly drops and the mean free path becomes
large. In those instances, you should use the flux limited diffusion
module in Castro to evolve a radiation field. However, if your
interest is only on the diffusion in the interior, you can use
the parameters:

 * ``castro.diffuse_cutoff_density``

 * ``castro.diffuse_cutoff_density_hi``

to specify a density,
below which, diffusion is not modeled. This is implemented in the
code by linearly scaling the conductivity to zero between these limits, e.g.,

.. math::

   \kth = \kth \cdot \frac{\rho - \mathtt{castro.diffuse\_cutoff\_density}}{\mathtt{castro.diffuse\_cutoff\_density\_hi} - \mathtt{castro.diffuse\_cutoff\_density}}


A simple test problem that sets up a Gaussian temperature profile
and does pure diffusion is provided as ``diffusion_test``.
.. _ch:radiation:

*********
Radiation
*********

Introduction
============

Castro has three radiation solvers:

-  SingleGroupSolver: this solver does not have radiation
   pressre. It is pure hydro plus radiation diffusion. This is only
   applicable when the medium is optically thick and the pressure is small.

-  SGFLDSolver: this is the gray flux-limited diffusion
   radiation hydrodymamics solver. Here the radiation pressure is
   separate from the gas pressure, and both pressures participate in
   the Riemann solver.

-  MGFLDSolver: this is the multigroup flux-limited diffusion
   radiation hydrodynamics solver. As with the gray solver, radiation
   pressure contributes to the pressure in the Riemann solver. Here a
   number of energy groups are used to represent the radiation field,
   and the opacity can be frequency-dependent.

The gray solver has a comoving frame mode and a mixed frame mode,
whereas the MG solver uses the comoving frame approach. More details
about the formulation and algorithm can be found in the series of
Castro papers.

Getting Started
===============

Getting the Code
----------------

The Castro radiation solver is part of the main Castro git repo,
so you already have all the Castro code and problem setups
to exercise radiation. The only other requirement is a copy
of the Hypre library. Hypre provides the algebraic multigrid
solvers used by the implicit radiation update. You can get
a copy at https://github.com/hypre-space/hypre (the minimum
supported release version is 2.23.0). Their install
instructions describe what to do; we recommend using the autotools
and GNU Make build. On HPC clusters, you typically want to build
with the same compiler you're using to build Castro, and you also
want to make sure that the options you're using for Castro are
compatible with the Hypre options, in particular when it comes to
``USE_MPI``, ``USE_OMP``, and ``USE_CUDA``.

As an example, to build Hypre on Summit with MPI and CUDA, you
should load the ``gcc/7.4.0`` and ``spectrum-mpi`` modules and
then do the following from the Hypre ``src/`` directory,
replacing ``/path/to/Hypre/install`` with the target location
where you want the Hypre files to be installed.
::

   HYPRE_CUDA_SM=70 CXX=mpicxx CC=mpicc FC=mpifort ./configure --prefix=/path/to/Hypre/install --with-MPI --with-cuda --enable-unified-memory
   make install

Then, when you are building Castro, you would build with
``USE_MPI=TRUE`` and ``USE_CUDA=TRUE``.

Castro looks for Hypre in the environment variable ``HYPRE_DIR``,
which you should point to the install directory you chose above.
Other than that, the only difference for builds with radiation
is that you must set ``USE_RAD=TRUE``.

Microphysics: EOS, Network, and Opacity
=======================================

EOS
---

Castro provides several types of equation of state (EOS), including
gamma-law and Helmholtz. To use the gamma-law EOS, set

::

    EOS_DIR := gamma_law

in the GNUmakefile.

The original Helmholtz EOS for stellar interiors includes a radiation
contribution. However, for radiation hydrodynamics calculations, the
radition contribution should be taken out of EOS because radiation has
been treated in other places. To use Helmholtz EOS, we will use the
version in Microphysics, as with the pure hydrodynamics code, but
this will interpret the RADIATION preprocessor variable and
disable the radiation portion of the EOS. If you have your own EOS, you
can put it in Microphysics.

There is also an artificial EOS that is used for several test cases called

::

   EOS_DIR := rad_power_law

This EOS should only be used for pure radiation-diffusion tests (i.e.
``castro.do_hydro = 0``). It defines the specific heat as a power law,

   .. math:: c_v = \mathrm{const}\ \rho^m T^{-n}

Network
-------

The radiation solver uses the same networks as we saw for pure hydro,
so nothing needs to change here. Again, if you are not modeling
reactions, then the general_null network can be used to define
the appropriate composition for your problem.

Opacity
-------

The most commonly used opacity setup is

.. math::
   \kappa_\nu = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p} ,
   :label: eq:kappa

where :math:`\kappa` is either the Planck or Rosseland absorption
coefficients, :math:`\rho` is density, :math:`T` is temperature, :math:`\nu` is
frequency, and :math:`m`, :math:`n` and :math:`p` are constants. For the gray solver,
:math:`p = 0`. If :eq:`eq:kappa` is sufficient, set

::

    Opacity_dir := rad_power_law

in your GNUmakefile. See § \ `3.3.1 <#sec:opacpars>`__ for instructions on how
to configure the parameters used for this opacity setup. If you would prefer a different
opacity mechanism, you will need to create your own opacity module by creating a new
directory in the Microphysics/opacity directory, creating the same set of subroutines
that the others have.

Some notes:

-  Here, :math:`\kappa` has units of :math:`\mathrm{cm}^{-1}`. Some papers or
   texts may instead have an implicit density factor in :math:`\kappa`,
   yielding units :math:`\mathrm{cm}^2~\mathrm{g}^{-1}`.

-  Castro allows for two temperatures (different radiation and gas
   temperature, so :math:`E_\mathrm{r} \ne a T_\mathrm{gas}^4`).
   Correspondingly, Castro cares about both the Planck mean,
   :math:`\kappa_P`, and Rosseland mean, :math:`\kappa_R`, opacities—these have
   different weightings.

   If we set :math:`\kappa_P \Delta x \gg 1` (:math:`\kappa_P` is really large),
   then the two temperatures become the same.

   If we set :math:`\kappa_P = \kappa_R`, then we can see how different the
   two temperatures are.

   In an optically thick medium, we would not expect the two temperatures
   to be very different.

.. _sec:opacpars:

Opacity Parameters
~~~~~~~~~~~~~~~~~~

The parameters describing the opacity include:

-  For the Planck opacity of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.
   ``opacity.const_kappa_p`` must be set positive to be used.

   -  ``opacity.const_kappa_p = -1.0``

   -  ``opacity.kappa_p_exp_m = 0.0``

   -  ``opacity.kappa_p_exp_n = 0.0``

   -  ``opacity.kappa_p_exp_p = 0.0``

-  For the Rosseland opacity of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.
   ``opacity.const_kappa_r`` must be set positive to be used.

   -  ``opacity.const_kappa_r = -1.0``

   -  ``opacity.kappa_r_exp_m = 0.0``

   -  ``opacity.kappa_r_exp_n = 0.0``

   -  ``opacity.kappa_r_exp_p = 0.0``

-  For the scattering coefficient of the form in :eq:`eq:kappa`,
   the following parameters set the coefficient and exponents.
   These are set via the ``opacity`` namespace in the inputs file.

   -  ``opacity.const_scatter = 0.0``

   -  ``opacity.scatter_exp_m = 0.0``

   -  ``opacity.scatter_exp_n = 0.0``

   -  ``opacity.scatter_exp_p = 0.0``

-  Since the formula above, :eq:`eq:kappa`, is non-physical and
   singular, we must set some floors in practice to prevent
   numerical issues. We have one floor for the opacity, which is
   applied to both the Planck and Rosseland opacities, and we
   also have a temperature floor. 

   -  ``opacity.kappa_floor = 1.d-50``

   -  ``opacity.rad_temp_floor = 0.0``

-  ``radiation.do_kappa_stm_emission = 0``

   If it is 1, correction for stimulated emission is applied to Planck mean as
   follows

   .. math::

      \kappa = \mathrm{const}\ \rho^{m} T^{-n} \nu^{p}
          \left [1-\exp{\left (-\frac{h\nu}{k T} \right )} \right ].

Note that the unit for opacities is :math:`\mathrm{cm}^{-1}`. For
the gray solver, the total opacity in the diffusion coefficient is the sum
of kappa_r and scattering, whereas for the MG solver,
there are two possibilities. If const_kappa_r is greater than
0, then the total opacity is set by kappa_r alone, otherwise
the total opacity is the sum of kappa_p and scattering.

Radiation Solver Physics
========================

In this section, we list some radiation related parameters that you
can set in an inputs file. Here are some important parameters:

-  radiation.SolverType:

   Set it to 5 for the gray solver, and 6 for the MG solver.

-  castro.do_hydro

   Usually you want to set it to 1. If it is set to 0,
   hydro will be turned off, and the calculation will only solve
   radiation diffusion equation.

-  castro.do_radiation

   If it is 0, the calculation will be pure hydro.

Below are more parameters. For each parameter, the default value is
on the right-hand side of the equal sign.

.. _sec:bothpar:

Verbosity and I/O
-----------------

-  radiation.v = 0

   Verbosity

-  radiation.verbose = 0

   Verbosity

-  radiation.plot_lambda = 0

   If 1, save flux limiter in plotfiles.

-  radiation.plot_kappa_p = 0

   If 1, save Planck mean opacity in plotfiles.

-  radiation.plot_kappa_r = 0

   If 1, save Rosseland mean opacity in plotfiles.

-  radiation.plot_lab_Er = 0

   If 1, save lab frame radiation energy density in plotfiles.
   This flag is ignored when the mixed-frame gray solver is used.

-  radiation.plot_com_flux = 0

   If 1, save comoving frame radiation flux in plotfiles.

-  radiation.plot_lab_flux = 0

   If 1, save lab frame radiation flux in plotfiles.

.. _sec:fluxlimiter:

Flux Limiter and Closure
------------------------

-  radiation.limiter = 2

   Possible values are:

   -   0: No flux limiter

   -   2: Approximate limiter of Levermore & Pomraning

   -  12: Bruenn’s limiter

   -  22: Larsen’s square root limiter

   -  32: Minerbo’s limiter

-  radiation.closure = 3

   Possible values are:

   -  0: :math:`f = \lambda`, where :math:`f` is the scalar Eddington factor
      and :math:`\lambda` is the flux limiter.

   -  1: :math:`f = \frac{1}{3}`

   -  2: :math:`f = 1 - 2 \lambda`

   -  3: :math:`f = \lambda + (\lambda R)^2`, where :math:`R` is the radiation
      Knudsen number.

   -  4: :math:`f = \frac{1}{3} + \frac{2}{3} (\frac{F}{cE})^2`, where
      :math:`F` is the radiation flux, :math:`E` is the radiation energy density,
      and :math:`c` is the speed of light.

Note the behavior of the radiative flux in the optically thin and
optically thick limits. The flux limiter, :math:`\lambda = \lambda(R)`,
where

.. math:: R = \frac{|\nabla E_r^{(0)}|}{\chi_R E_r^{(0)}}

Regardless of the limiter chosen, when we are optically thick,
:math:`\chi_R \rightarrow \infty`, :math:`R \rightarrow 0`, and :math:`\lambda \rightarrow 1/3`.
The radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     \frac{1}{3} \frac{c}{\chi_R} \nabla E_r^{(0)}

And when we are optically thin, :math:`\chi_R \rightarrow 0`, :math:`R \rightarrow \infty`,
and :math:`\lambda \rightarrow 1/R = \chi_R E_r^{(0)}/{|\nabla E_r^{0}|}`, and
the radiative flux then becomes

.. math::

   F_r^{(0)} = -\frac{c\lambda}{\chi_R} \nabla E_r^{(0)} \rightarrow
     -\frac{c}{\chi_R}\frac{\chi_R E_r^{(0)}}{|\nabla E_r^{0}|}
       \nabla E_r^{(0)} = -c E_r^{0}

See Krumholz et al. 2007 for some discussion on this.

Boundary Conditions
-------------------

The following parameters are for the radiation boundary in the diffusion
equation. They do not affect hydrodynamic boundaries.

-  radiation.lo_bc

   This sets the action to take at the lower edge of the domain in
   each coordinate direction. Possible values are:

   -  101 *Dirichlet*:

      Specify the radiation energy density on the boundary.
      For gray radiation, this could be :math:`E_r = a T^4`.

      For multigroup radiation, Castro stores the energy density as
      :math:`\mathrm{erg}~\mathrm{cm}^{-3}`, so the total radiation energy
      can be found by simply summing over the groups. So if you want
      to set the radiation BCs using the Planck function, you simply
      multiply by the group width—see Exec/radiation_tests/RadSphere/Tools/radbc.f90
      for an example.

   -  102 *Neumann*:

      Here, you specify the radiation flux on the boundary. For gray
      radiation, this is the expression given in the gray Castro paper
      (Eq. 7, 8),

      .. math:: F_r = - \frac{c\lambda}{\kappa_R} \nabla E_r

      where :math:`\lambda` is the flux limiter.

      Note that if your boundary represents an incoming flux through
      a vacuum (like stellar irradiation), then :math:`\kappa \rightarrow 0`, leaving

      .. math:: F_r = -c E_r

      (see § \ `4.2 <#sec:fluxlimiter>`__) in that case.

   -  104 *Marshak* (vacuum):

      Here, you specify the incident flux and the outside is a vacuum.
      This differs from the Neumann condition because there is also a
      flux coming from inside, for the net flux across the boundary is
      different than the incident flux.

   -  105 *Sanchez-Pomraning*:

      This is a modified form of the Marshak boundary condition that works with FLD.
      This is like the Marshak condition, but :math:`\lambda = 1/3` is not assumed inside
      the boundary (optical thickness).

-  radiation.hi_bc

   See radiation.lo_bc.

-  radiation.lo_bcflag = 0 0 0

   If it is 0, bcval is used for that dimension, otherwise
   subroutine rbndry in RadBndry_1d.f90 is called to set
   boundary conditions.

-  radiation.hi_bcflag = 0 0 0

   See radiation.lo_bcflag

-  radiation.lo_bcval = 0.0 0.0 0.0

   The actual value to impose for the boundary condition type set by
   radiation.lo_bc. This parameter is interpreted differently
   depending on the boundary condition:

   -  Dirchlet: Dirichlet value of rad energy density

   -  Neumann: inward flux of rad energy

   -  Marshak: incident flux

   -  Sanchez-Pomraning: incident flux

-  radiation.hi_bcval = 0.0 0.0 0.0

   See radiation.lo_bcval

Convergence
-----------

For the gray solver, there is only one iteration in the scheme,
whereas for the MG solver, there are two iterations with an inner
iteration embedded inside an outer iteration. In the following, the
iteration in the gray solver will also be referred as the outer
iteration for convenience. The parameters for the inner iteration are
irrelevant to the gray solver.

radiation.maxiter = 50
    |
    | Maximal number of outer iteration steps.

radiation.miniter = 1
    |
    | Minimal number of outer iteration steps.

radiation.reltol = 1.e-6
    |
    | Relative tolerance for the outer iteration.

radiation.abstol = 0.0
    |
    | Absolute tolerance for the outer iteration.

radiation.maxInIter = 30
    |
    | Maximal number of inner iteration steps.

radiation.minInIter = 1
    |
    | Minimal number of inner iteration steps.

radiation.relInTol = 1.e-4
    |
    | Relative tolerance for the inner iteration.

radiation.absInTol = 0.0
    |
    | Absolute tolerance for the inner iteration.

radiation.convergence_check_type = 0
    |
    | For the MG solver only. This specifiy the way of checking the
      convergence of an outer iteration. Possible values are

    -  0: Check :math:`T`, :math:`Y_e`, and the residues of the equations for
       :math:`\rho e` and :math:`\rho Y_e`

    -  1: Check :math:`\rho e`

    -  2: Check the residues of the equations for :math:`\rho e` and :math:`\rho Y_e`

    -  3: Check :math:`T` and :math:`Y_e`

.. _sec:graypar:

Parameters for Gray Solver
--------------------------

radiation.comoving = 1
    |
    | Do we use the comoving frame approach?

radiation.Er_Lorentz_term = 1
    |
    | If the mixed-frame approach is taken, this parameter decides whether
      Lorentz transformation terms are retained.

radiation.delta_temp = 1.0
    |
    | This is used in computing numerical derivativas with respect to :math:`T`.
      So it should be a small number compared with :math:`T`, but not too small.

radiation.update_limiter = 1000
    |
    | Stop updating flux limiter after update_limiter iteration steps.

radiation.update_planck = 1000
    |
    | Stop updating Planck mean opacity after update_planck iteration steps.

radiation.update_rosseland = 1000
    |
    | Stop updating Rosseland mean opacity after update_rosseland iteration steps.

Grouping in the MG Solver
-------------------------

We provide two methods of setting up groups based upon logarithmic
spacing. In both methods, you must provide:

radiation.nGroups
    |
    | Number of groups.

radiation.lowestGroupHz
    |
    | Frequency of the lower bound for the first group.

In addition, if the parameter groupGrowFactor is provided, then
the first method will be used, otherwise the second method will be
used. In the first way, you must also provide firstGroupWidthHz
(the width of the first group). The width of other groups is set to
be groupGrowFactor times the width of its immediately preceding
group. In the second way, you must provide highestGroupHz as
the upper bound of the last group. It should be noted that
lowestGroupHz can be 0 in the first method, but not the second
method. However, when we compute the group-integrated Planck
function, the lower bound for the first group and the upper bound for
the last group are assumed to be 0 and :math:`\infty`, respectively.

.. _sec:mgpar:

Parameters for MG Solver
------------------------

radiation.delta_e_rat_dt_tol = 100.0
    |
    | Maximally allowed relative change in :math:`e` during one time step.

radiation.delta_T_rat_dt_tol = 100.0
    |
    | Maximally allowed relative change in :math:`T` during one time step.

radiation.delta_Ye_dt_tol = 100.0
    |
    | Maximally allowed absolute change in :math:`Y_e` during one tim estep.

radiation.fspace_advection_type = 2
    |
    | Possible value is 1 or 2. The latter is better.

radiation.integrate_Planck = 1
    |
    | If 1, integrate Planck function for each group. For the first
      group, the lower bound in the integration is assumed to be 0 no
      matter what the grouping is. For the last group, the upper bound in
      the integration is assumed to be :math:`\infty`.

radiation.matter_update_type = 0
    |
    | How to update matter. 0 is proabaly the best.

radiation.accelerate = 2
    |
    | The inner iteration of the MG solver usually requires an
      acceleration scheme. Choices are

    -  0: No acceleration

    -  1: Local acceleration

    -  2: Gray acceleration

radiation.skipAccelAllowed = 0
    |
    | If it is set to 1, skip acceleration if it does not help.

radiation.n_bisect = 1000
    |
    | Do bisection for the outer iteration after n_bisec iteration steps.

radiation.use_dkdT = 1
    |
    | If it is 1, :math:`\frac{\partial \kappa}{\partial T}` is retained in the
      Jacobi matrix for the outer (Newton) iteration.

radiation.update_opacity = 1000
    |
    | Stop updating opacities after update_opacity outer iteration steps.

radiation.inner_update_limiter = 0
    |
    | Stop updating flux limiter after inner_update_limiter inner
      iteration steps. If it is 0, the limiter is lagged by one outer
      iteration. If it is -1, the limiter is lagged by one time step. If
      the inner iteration has difficulty in converging, setting this
      parameter it to -1 can help. Since the flux limiter is only a
      kludge, it is justified to lag it.

.. _sec:hypre:

Linear System Solver
--------------------

There are a number of choices for the linear system solver. The
performance of the solvers usually depends on problems and the
computer. So it is worth trying a few solvers to find out which one
is best for your problem and computer.

radsolve.level_solver_flag: the linear solver
in Hypre to use. The available choices are:

-  0: SMG

-  1: PFMG (:math:`\ge` 2-d only)

-  100: AMG using ParCSR ObjectType

-  102: GMRES using ParCSR ObjectType

-  103: GMRES using SStruct ObjectType

-  104: GMRES using AMG as preconditioner

-  109: GMRES using Struct SMG/PFMG as preconditioner

-  150: AMG using ParCSR ObjectType

-  1002: PCG using ParCSR ObjectType

-  1003: PCG using SStruct ObjectType

As a general rule, the SMG is the most stable solver, but is usually
the slowest. The asymmetry in the linear system comes from the
adaptive mesh, so the PFMG should be your first choice. Note: in
you cannot use PFMG.

Setting this to 109 (GMRES using Struct SMG/PFMG as preconditioner)
should work reasonably well for most problems.

radsolve.maxiter (default: 40):
Maximal number of iteration in Hypre.

radsolve.reltol (default: 1.e-10):
Relative tolerance in Hypre

radsolve.abstol (default: 0):
Absolute tolerance in Hypre

radsolve.v (default: 0):
Verbosity

radsolve.verbos (default: 0):
Verbosity

habec.verbose (default: 0):
Verbosity for level_solver_flag :math:`<` 100

hmabec.verbose (default: 0):
Verbosity for level_solver_flag :math:`>=` 100

Output
======

Gray Solver
-----------

For the gray radiation solver, the radiation energy density is stored in plotfiles
as rad. Note that this quantity has units of :math:`\mathrm{erg~cm^{-3}}`, which
is different that the specify internal energy of the gas :math:`\mathrm{erg~g^{-1}}`.
.. _ch:hydro:

*************
Hydrodynamics
*************

Introduction
============

The hydrodynamics scheme in Castro implements an unsplit
second-order Godunov method. Characteristic tracing is used to
time-center the input states to the Riemann solver. The same
hydrodynamics routines are used for pure hydro and radiation
hydrodynamics.

Some general notes:

-  Regardless of the dimensionality, we always carry around all 3
   components of velocity/momentum—this allows for rotation sources easily.

-  When radiation is enabled (via RADIATION), we discuss
   the gas and radiation quantities separately. This generally applies
   to the temperature, pressure, internal energy, various adiabatic
   indices, and sound speed. When we refer to the “total” value of
   one of these, it means that both gas and radiation contributions
   are included. When we refer to the “gas” quantity, this is what
   the equation of state would return.

   For continuity, we continue to use the “gas” qualifier even if we
   are not solving the radiation hydrodynamics equations. In this
   case, it still means that it comes through the equation of state,
   but note some of our equations of state (like the helmeos) include a
   radiation pressure contribution when we are running without
   radiation hydrodynamics enabled. In this case, we still refer to
   this as the “gas”.

Hydrodynamics Data Structures
=============================

Within the Fortran routines that implement the hydrodynamics, there are
several main data structures that hold the state.

-  conserved state: these arrays generally begin with ``u``,
   e.g., ``uin``, ``uout``. The ``NVAR``
   components for the state data in the array are accessed using
   integer keys defined in :numref:`table:consints`.

   .. _table:consints:
   .. table:: The integer variables to index the conservative state array

      +-----------------------+-----------------------+-------------------------+
      | **variable**          | **quantity**          | **note**                |
      +=======================+=======================+=========================+
      | ``URHO``              | :math:`\rho`          |                         |
      +-----------------------+-----------------------+-------------------------+
      | ``UMX``               | :math:`\rho u`        |                         |
      +-----------------------+-----------------------+-------------------------+
      | ``UMY``               | :math:`\rho v`        |                         |
      +-----------------------+-----------------------+-------------------------+
      | ``UMZ``               | :math:`\rho w`        |                         |
      +-----------------------+-----------------------+-------------------------+
      | ``UEDEN``             | :math:`\rho E`        |                         |
      +-----------------------+-----------------------+-------------------------+
      | ``UEINT``             | :math:`\rho e`        | this is computed from   |
      |                       |                       | the other quantities    |
      |                       |                       | using                   |
      |                       |                       | :math:`\rho e = \rho    |
      |                       |                       | E - \rho |\ub|^2        |
      |                       |                       | / 2`                    |
      +-----------------------+-----------------------+-------------------------+
      | ``UTEMP``             | :math:`T`             | this is computed from   |
      |                       |                       | the other quantities    |
      |                       |                       | using the EOS           |
      +-----------------------+-----------------------+-------------------------+
      | ``UFA``               | :math:`\rho A_1`      | the first advected      |
      |                       |                       | quantity                |
      +-----------------------+-----------------------+-------------------------+
      | ``UFS``               | :math:`\rho X_1`      | the first species       |
      +-----------------------+-----------------------+-------------------------+
      | ``UFX``               | :math:`\rho Y_1`      | the first auxiliary     |
      |                       |                       | variable                |
      +-----------------------+-----------------------+-------------------------+
      | ``USHK``              | a shock flag          | (used for shock         |
      |                       |                       | detection)              |
      +-----------------------+-----------------------+-------------------------+
      | ``UMR``               | radial momentum       | (if ``HYBRID_MOMENTUM`` |
      |                       |                       | is defined)             |
      +-----------------------+-----------------------+-------------------------+
      | ``UML``               | angular momentum      | (if ``HYBRID_MOMENTUM`` |
      |                       |                       | is defined)             |
      +-----------------------+-----------------------+-------------------------+
      | ``UMP``               | vertical momentum     | (if ``HYBRID_MOMENTUM`` |
      |                       |                       | is defined)             |
      +-----------------------+-----------------------+-------------------------+

-  primitive variable state: these arrays generally simply called
   ``q``, and has ``NQ`` components.

   A related quantity is ``NQSRC`` which is the number of primitive variable
   source terms.  ``NQSRC`` ≤ ``NQ``.

   .. note:: if ``RADIATION`` is defined, then only the gas/hydro terms are
      present in ``NQSRC``.  

   :numref:`table:primlist` gives the names of the primitive variable integer
   keys for accessing these arrays. Note, unless otherwise specified the quantities without a subscript
   are “gas” only and those with the “tot” subscript are “gas + radiation”.

   .. _table:primlist:
   .. table:: Primitive State Vector Integer Keys

      +-----------------------+------------------------+-----------------------+
      | **variable**          | **quantity**           | **note**              |
      +=======================+========================+=======================+
      | ``QRHO``              | :math:`\rho`           |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QU``                | :math:`u`              |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QV``                | :math:`v`              |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QW``                | :math:`w`              |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QPRES``             | :math:`p`              |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QREINT``            | :math:`\rho e`         |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QTEMP``             | :math:`T`              |                       |
      +-----------------------+------------------------+-----------------------+
      | ``QFA``               | :math:`A_1`            | the first advected    |
      |                       |                        | quantity              |
      +-----------------------+------------------------+-----------------------+
      | ``QFS``               | :math:`X_1`            | the first species     |
      +-----------------------+------------------------+-----------------------+
      | ``QFX``               | :math:`Y_1`            | the first auxiliary   |
      |                       |                        | variable              |
      +-----------------------+------------------------+-----------------------+
      | ``QPTOT``             | :math:`p_\mathrm{tot}` | the total pressure,   |
      |                       |                        | gas + radiation       |
      +-----------------------+------------------------+-----------------------+
      | ``QREITOT``           | :math:`e_\mathrm{tot}` | the total specific    |
      |                       |                        | internal energy, gas  |
      |                       |                        | + radiation           |
      +-----------------------+------------------------+-----------------------+
      | ``QRAD``              | :math:`E_r`            | the radiation energy  |
      |                       |                        | (there are ngroups of |
      |                       |                        | these)                |
      +-----------------------+------------------------+-----------------------+

-  auxiliary primitive variables: these arrays are generally called
   qaux. The main difference between these and the regular
   primitive variables is that we do not attempt to do any
   reconstruction on their profiles. There are ``NQAUX`` quantities, indexed
   by the integer keys listed in :numref:`table:qauxlist`.
   Note, unless otherwise specified the quantities without a subscript are “gas”
   only and those with the “tot” subscript are “gas + radiation”.

   .. _table:qauxlist:
   .. table:: The integer variable keys for accessing the auxiliary primitive state vector, quax.

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | ``QGAMC``             | :math:`\gamma_1`      | the first adiabatic   |
      |                       |                       | exponent, as returned |
      |                       |                       | from the EOS          |
      +-----------------------+-----------------------+-----------------------+
      | ``QC``                | :math:`c_s`           | the sound speed, as   |
      |                       |                       | returned from the EOS |
      +-----------------------+-----------------------+-----------------------+
      | ``QGAMCG``            | :math:`{\Gamma_1      | includes radiation    |
      |                       | }_\mathrm{tot}`       | components (defined   |
      |                       |                       | only if ``RADIATION`` |
      |                       |                       | is defined)           |
      +-----------------------+-----------------------+-----------------------+
      | ``QCG``               | :math:`{c_s           | total sound speed     |
      |                       | }_\mathrm{tot}`       | including radiation   |
      |                       |                       | (defined only if      |
      |                       |                       | ``RADIATION`` is      |
      |                       |                       | defined)              |
      +-----------------------+-----------------------+-----------------------+
      | ``QLAMS``             | :math:`\lambda_f`     | the ``ngroups`` flux  |
      |                       |                       | limiters (defined     |
      |                       |                       | only if ``RADIATION`` |
      |                       |                       | is defined)           |
      +-----------------------+-----------------------+-----------------------+

-  interface variables: these are the time-centered interface states
   returned by the Riemann solver. They are used to discretize some
   non-conservative terms in the equations. These arrays are generally
   called ``q1``, ``q2``, and ``q3`` for the x, y, and z
   interfaces respectively. There are ``NGDNV`` components accessed with
   the integer keys defined in :numref:`table:gdlist`
   Note, unless otherwise specified the quantities without a subscript are
   “gas” only and those with the “tot” subscript are “gas + radiation”.

   .. _table:gdlist:
   .. table:: The integer variable keys for accessing the Godunov interface state vectors.

      +-----------------------+-----------------------+-----------------------+
      | **variable**          | **quantity**          | **note**              |
      +=======================+=======================+=======================+
      | ``QGDRHO``            | :math:`\rho`          |                       |
      +-----------------------+-----------------------+-----------------------+
      | ``QDU``               | :math:`u`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | ``QDV``               | :math:`v`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | ``QDW``               | :math:`w`             |                       |
      +-----------------------+-----------------------+-----------------------+
      | ``QDPRES``            | :math:`p`             | regardless of whether |
      |                       |                       | ``RADIATION`` is      |
      |                       |                       | defined,              |
      |                       |                       | this is always just   |
      |                       |                       | the gas pressure      |
      +-----------------------+-----------------------+-----------------------+
      | ``QDLAMS``            | :math:`{\lambda_f}`   | the starting index    |
      |                       |                       | for the flux          |
      |                       |                       | limiter—there are     |
      |                       |                       | ngroups components    |
      |                       |                       | (defined only if      |
      |                       |                       | ``RADIATION`` is      |
      |                       |                       | defined)              |
      +-----------------------+-----------------------+-----------------------+
      | ``QDERADS``           | :math:`E_r`           | the starting index    |
      |                       |                       | for the radiation     |
      |                       |                       | energy—there are      |
      |                       |                       | ngroups components    |
      |                       |                       | (defined only if      |
      |                       |                       | ``RADIATION`` is      |
      |                       |                       | defined)              |
      +-----------------------+-----------------------+-----------------------+

Conservation Forms
==================

We begin with the fully compressible equations for the conserved state vector,
:math:`\Ub = (\rho, \rho \ub, \rho E, \rho A_k, \rho X_k, \rho Y_k):`

.. math::

   \begin{align}
   \frac{\partial \rho}{\partial t} &= - \nabla \cdot (\rho \ub) + S_{{\rm ext},\rho}, \\
   \frac{\partial (\rho \ub)}{\partial t} &= - \nabla \cdot (\rho \ub \ub) - \nabla p +\rho \gb + \Sb_{{\rm ext},\rho\ub}, \\
   \frac{\partial (\rho E)}{\partial t} &= - \nabla \cdot (\rho \ub E + p \ub) + \rho \ub \cdot \gb - \sum_k {\rho q_k \dot\omega_k} + \nabla\cdot\kth\nabla T + S_{{\rm ext},\rho E}, \\
   \frac{\partial (\rho A_k)}{\partial t} &= - \nabla \cdot (\rho \ub A_k) + S_{{\rm ext},\rho A_k}, \\
   \frac{\partial (\rho X_k)}{\partial t} &= - \nabla \cdot (\rho \ub X_k) + \rho \dot\omega_k + S_{{\rm ext},\rho X_k}, \\
   \frac{\partial (\rho Y_k)}{\partial t} &= - \nabla \cdot (\rho \ub Y_k) + S_{{\rm ext},\rho Y_k}.\label{eq:compressible-equations}
   \end{align}

Here :math:`\rho, \ub, T, p`, and :math:`\kth` are the density,
velocity, temperature, pressure, and thermal conductivity,
respectively, and :math:`E = e + \ub \cdot \ub / 2` is the total
energy with :math:`e` representing the internal energy. In addition,
:math:`X_k` is the abundance of the :math:`k^{\rm th}` isotope, with
associated production rate, :math:`\dot\omega_k`, and energy release,
:math:`q_k`. Here :math:`\gb` is the gravitational vector, and
:math:`S_{{\rm ext},\rho}, \Sb_{{\rm ext}\rho\ub}`, etc., are
user-specified source terms. :math:`A_k` is an advected quantity,
i.e., a tracer. We also carry around auxiliary variables, :math:`Y_k`,
which have a user-defined evolution equation, but by default are
treated as advected quantities.  These are meant to be defined in the network.

In the code we also carry around :math:`T` and :math:`\rho e` in the conservative
state vector even though they are derived from the other conserved
quantities. The ordering of the elements within :math:`\Ub` is defined
by integer variables into the array—see
:numref:`table:consints`.

Some notes:

-  Regardless of the dimensionality of the problem, we always carry
   all 3 components of the velocity. This allows for, e.g., 2.5-d
   rotation (advecting the component of velocity out of the plane in
   axisymmetric coordinates).

   You should always initialize all velocity components to zero, and
   always construct the kinetic energy with all three velocity components.

-  There are ``NADV`` advected quantities, which range from
   ``UFA: UFA+nadv-1``. The advected quantities have no effect at all on
   the rest of the solution but can be useful as tracer quantities.

-  There are ``NSPEC`` species (defined in the network
   directory), which range from ``UFS: UFS+nspec-1``.

-  There are ``NAUX`` auxiliary variables, from ``UFX:UFX+naux-1``. The
   auxiliary variables are passed into the equation of state routines
   along with the species. An example of an auxiliary variable is the
   electron fraction, :math:`Y_e`, in core collapse simulations.  The
   number and names of the auxiliary variables are defined in the
   network.


Source Terms
============

We now compute explicit source terms for each variable in :math:`\Qb` and
:math:`\Ub`. The primitive variable source terms will be used to construct
time-centered fluxes. The conserved variable source will be used to
advance the solution. We neglect reaction source terms since they are
accounted for in **Steps 1** and **6**. The source terms are:

.. math::

   \Sb_{\Qb}^n =
   \left(\begin{array}{c}
   S_\rho \\
   \Sb_{\ub} \\
   S_p \\
   S_{\rho e} \\
   S_{A_k} \\
   S_{X_k} \\
   S_{Y_k}
   \end{array}\right)^n
   =
   \left(\begin{array}{c}
   S_{{\rm ext},\rho} \\
   \gb + \frac{1}{\rho}\Sb_{{\rm ext},\rho\ub} \\
   \frac{1}{\rho}\frac{\partial p}{\partial e}S_{{\rm ext},\rho E} + \frac{\partial p}{\partial\rho}S_{{\rm ext}\rho} \\
   \nabla\cdot\kth\nabla T + S_{{\rm ext},\rho E} \\
   \frac{1}{\rho}S_{{\rm ext},\rho A_k} \\
   \frac{1}{\rho}S_{{\rm ext},\rho X_k} \\
   \frac{1}{\rho}S_{{\rm ext},\rho Y_k}
   \end{array}\right)^n,

.. math::

   \Sb_{\Ub}^n =
   \left(\begin{array}{c}
   \Sb_{\rho\ub} \\
   S_{\rho E} \\
   S_{\rho A_k} \\
   S_{\rho X_k} \\
   S_{\rho Y_k}
   \end{array}\right)^n
   =
   \left(\begin{array}{c}
   \rho \gb + \Sb_{{\rm ext},\rho\ub} \\
   \rho \ub \cdot \gb + \nabla\cdot\kth\nabla T + S_{{\rm ext},\rho E} \\
   S_{{\rm ext},\rho A_k} \\
   S_{{\rm ext},\rho X_k} \\
   S_{{\rm ext},\rho Y_k}
   \end{array}\right)^n.

Primitive Forms
===============

Castro uses the primitive form of the fluid equations, defined in terms of
the state :math:`\Qb = (\rho, \ub, p, \rho e, A_k, X_k, Y_k)`, to construct the
interface states that are input to the Riemann problem.

The primitive variable equations for density, velocity, and pressure are:

.. math::

   \begin{align}
     \frac{\partial\rho}{\partial t} &= -\ub\cdot\nabla\rho - \rho\nabla\cdot\ub + S_{{\rm ext},\rho} \\
   %
     \frac{\partial\ub}{\partial t} &= -\ub\cdot\nabla\ub - \frac{1}{\rho}\nabla p + \gb + 
   \frac{1}{\rho} (\Sb_{{\rm ext},\rho\ub} - \ub \; S_{{\rm ext},\rho}) \\
   \frac{\partial p}{\partial t} &= -\ub\cdot\nabla p - \rho c^2\nabla\cdot\ub +
   \left(\frac{\partial p}{\partial \rho}\right)_{e,X}S_{{\rm ext},\rho}\nonumber\\
   &+\  \frac{1}{\rho}\sum_k\left(\frac{\partial p}{\partial X_k}\right)_{\rho,e,X_j,j\neq k}\left(\rho\dot\omega_k + S_{{\rm ext},\rho X_k} - X_kS_{{\rm ext},\rho}\right)\nonumber\\
   & +\  \frac{1}{\rho}\left(\frac{\partial p}{\partial e}\right)_{\rho,X}\left[-eS_{{\rm ext},\rho} - \sum_k\rho q_k\dot\omega_k + \nabla\cdot\kth\nabla T \right.\nonumber\\
   & \quad\qquad\qquad\qquad+\ S_{{\rm ext},\rho E} - \ub\cdot\left(\Sb_{{\rm ext},\rho\ub} - \frac{\ub}{2}S_{{\rm ext},\rho}\right)\Biggr] 
   \end{align}

The advected quantities appear as:

.. math::

   \begin{align}
   \frac{\partial A_k}{\partial t} &= -\ub\cdot\nabla A_k + \frac{1}{\rho}
                                        ( S_{{\rm ext},\rho A_k} - A_k S_{{\rm ext},\rho} ), \\
   \frac{\partial X_k}{\partial t} &= -\ub\cdot\nabla X_k + \dot\omega_k + \frac{1}{\rho}
                                        ( S_{{\rm ext},\rho X_k}  - X_k S_{{\rm ext},\rho} ), \\
   \frac{\partial Y_k}{\partial t} &= -\ub\cdot\nabla Y_k + \frac{1}{\rho} 
                                        ( S_{{\rm ext},\rho Y_k}  - Y_k S_{{\rm ext},\rho} ).
   \end{align}

All of the primitive variables are derived from the conservative state
vector, as described in Section `6.1 <#Sec:Compute Primitive Variables>`__.
When accessing the primitive variable state vector, the integer variable
keys for the different quantities are listed in :numref:`table:primlist`.

Internal energy and temperature
-------------------------------

We augment the above system with an internal energy equation:

.. math::

   \begin{align}
   \frac{\partial(\rho e)}{\partial t} &= - \ub\cdot\nabla(\rho e) - (\rho e+p)\nabla\cdot\ub - \sum_k \rho q_k\dot\omega_k 
                                           + \nabla\cdot\kth\nabla T + S_{{\rm ext},\rho E} \nonumber\\
   & -\  \ub\cdot\left(\Sb_{{\rm ext},\rho\ub}-\frac{1}{2}S_{{\rm ext},\rho}\ub\right),
   \end{align}

This has two benefits. First, for a general equation of state,
carrying around an additional thermodynamic quantity allows us to
avoid equation of state calls (in particular, in the Riemann solver,
see e.g. :cite:`colglaz`). Second, it is sometimes the case that the
internal energy calculated as

.. math:: e_T \equiv E - \frac{1}{2} \mathbf{v}^2

is
unreliable. This has two usual causes: one, for high Mach number
flows, the kinetic energy can dominate the total gas energy, making
the subtraction numerically unreliable; two, if you use gravity or
other source terms, these can indirectly alter the value of the
internal energy if obtained from the total energy.

To provide a more reasonable internal energy for defining the
thermodynamic state, we have implemented the dual energy formalism
from ENZO :cite:`bryan:1995`, :cite:`bryan:2014`, where we switch
between :math:`(\rho e)` and :math:`(\rho e_T)` depending on the local
state of the fluid. To do so, we define parameters :math:`\eta_1`,
:math:`\eta_2`, and :math:`\eta_3`, corresponding to the code
parameters castro.dual_energy_eta1, castro.dual_energy_eta2, and
castro.dual_energy_eta3. We then consider the ratio :math:`e_T / E`,
the ratio of the internal energy (derived from the total energy) to
the total energy. These parameters are used as follows:

-  :math:`\eta_1`: If :math:`e_T > \eta_1 E`, then we use :math:`e_T` for the purpose
   of calculating the pressure in the hydrodynamics update. Otherwise,
   we use the :math:`e` from the internal energy equation in our EOS call to
   get the pressure.

-  :math:`\eta_2`: At the end of each hydro advance, we examine whether
   :math:`e_T > \eta_2 E`. If so, we reset :math:`e` to be equal to :math:`e_T`,
   discarding the results of the internal energy equation. Otherwise,
   we keep :math:`e` as it is.

-  :math:`\eta_3`: Similar to :math:`\eta_1`, if :math:`e_T > \eta_3 E`, we use
   :math:`e_T` for the purposes of our nuclear reactions, otherwise, we use
   :math:`e`.

Note that our version of the internal energy equation does not require
an artificial viscosity, as used in some other hydrodynamics
codes. The update for :math:`(\rho e)` uses information from the Riemann
solve to calculate the fluxes, which contains the information
intrinsic to the shock-capturing part of the scheme.

In the code we also carry around :math:`T` in the primitive state vector.

Primitive Variable System
-------------------------

The full primitive variable form (without the advected or auxiliary
quantities) is

.. math:: \frac{\partial\Qb}{\partial t} + \sum_d \Ab_d\frac{\partial\Qb}{\partial x_d} = \Sb_{\Qb}.

For example, in 2D:

.. math::

   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_t
   +
   \left(\begin{array}{cccccc}
   u & \rho & 0 & 0 & 0 & 0 \\
   0 & u & 0 & \frac{1}{\rho} & 0 & 0 \\
   0 & 0 & u & 0 & 0 & 0 \\
   0 & \rho c^2 & 0 & u & 0 & 0 \\
   0 & \rho e + p & 0 & 0 & u & 0 \\
   0 & 0 & 0 & 0 & 0 & u
   \end{array}\right)
   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_x
   +
   \left(\begin{array}{cccccc}
   v & 0 & \rho & 0 & 0 & 0 \\
   0 & v & 0 & 0 & 0 & 0 \\
   0 & 0 & v & \frac{1}{\rho} & 0 & 0 \\
   0 & 0 & \rho c^2 & v & 0 & 0 \\
   0 & 0 & \rho e + p & 0 & v & 0 \\
   0 & 0 & 0 & 0 & 0 & v
   \end{array}\right)
   \left(\begin{array}{c}
   \rho \\
   u \\
   v \\
   p \\
   \rho e \\
   X_k
   \end{array}\right)_y
   =
   \Sb_\Qb

The eigenvalues are:

.. math:: {\bf \Lambda}(\Ab_x) = \{u-c,u,u,u,u,u+c\}, \qquad {\bf \Lambda}(\Ab_y) = \{v-c,v,v,v,v,v+c\} .

The right column eigenvectors are:

.. math::

   \Rb(\Ab_x) =
   \left(\begin{array}{cccccc}
   1 & 1 & 0 & 0 & 0 & 1 \\
   -\frac{c}{\rho} & 0 & 0 & 0 & 0 & \frac{c}{\rho} \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   c^2 & 0 & 0 & 0 & 0 & c^2 \\
   h & 0 & 0 & 1 & 0 & h \\
   0 & 0 & 0 & 0 & 1 & 0 \\
   \end{array}\right),
   \qquad
   \Rb(\Ab_y) =
   \left(\begin{array}{cccccc}
   1 & 1 & 0 & 0 & 0 & 1 \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   -\frac{c}{\rho} & 0 & 0 & 0 & 0 & \frac{c}{\rho} \\
   c^2 & 0 & 0 & 0 & 0 & c^2 \\
   h & 0 & 0 & 1 & 0 & h \\
   0 & 0 & 0 & 0 & 1 & 0 \\
   \end{array}\right).

The left row eigenvectors, normalized so that :math:`\Rb_d\cdot\Lb_d = \Ib` are:

.. math::

   \Lb_x =
   \left(\begin{array}{cccccc}
   0 & -\frac{\rho}{2c} & 0 & \frac{1}{2c^2} & 0 & 0 \\
   1 & 0 & 0 & -\frac{1}{c^2} & 0 & 0 \\
   0 & 0 & 1 & 0 & 0 & 0 \\
   0 & 0 & 0 & -\frac{h}{c^2} & 1 & 0 \\
   0 & 0 & 0 & 0 & 0 & 1 \\
   0 & \frac{\rho}{2c} & 0 & \frac{1}{2c^2} & 0 & 0
   \end{array}\right),
   \qquad
   \Lb_y =
   \left(\begin{array}{cccccc}
   0 & 0 & -\frac{\rho}{2c} & \frac{1}{2c^2} & 0 & 0 \\
   1 & 0 & 0 & -\frac{1}{c^2} & 0 & 0 \\
   0 & 1 & 0 & 0 & 0 & 0 \\
   0 & 0 & 0 & -\frac{h}{c^2} & 1 & 0 \\
   0 & 0 & 0 & 0 & 0 & 1 \\
   0 & 0 & \frac{\rho}{2c} & \frac{1}{2c^2} & 0 & 0
   \end{array}\right).

.. _Sec:Advection Step:

Hydrodynamics Update
====================

There are four major steps in the hydrodynamics update:

#. Converting to primitive variables

#. Construction the edge states

#. Solving the Riemann problem

#. Doing the conservative update

.. index:: castro.do_hydro, castro.add_ext_src, castro.do_sponge, castro.normalize_species, castro.spherical_star, castro.show_center_of_mass

Each of these steps has a variety of runtime parameters that
affect their behavior. Additionally, there are some general
runtime parameters for hydrodynamics:

-  ``castro.do_hydro``: time-advance the fluid dynamical
   equations (0 or 1; must be set)

-  ``castro.add_ext_src``: include additional user-specified
   source term (0 or 1; default 0)

-  ``castro.do_sponge``: call the sponge routine
   after the solution update (0 or 1; default: 0)

   See :ref:`sponge_section` for more details on the sponge.

-  ``castro.normalize_species``: enforce that :math:`\sum_i X_i = 1`
   (0 or 1; default: 0)

-  ``castro.spherical_star``: this is used to set the boundary
   conditions by assuming the star is spherically symmetric in
   the outer regions (0 or 1; default: 0)

   When used, Castro averages the values at a given radius over the
   cells that are inside the domain to define a radial function. This
   function is then used to set the values outside the domain in
   implementing the boundary conditions.

-  ``castro.show_center_of_mass``: (0 or 1; default: 0)

.. index:: castro.small_dens, castro.small_temp, castro.small_pres

Several floors are imposed on the thermodynamic quantities to prevet unphysical
behavior:

-  ``castro.small_dens``: (Real; default: -1.e20)

-  ``castro.small_temp``: (Real; default: -1.e20)

-  ``castro.small_pres``: (Real; default: -1.e20)

.. _Sec:Compute Primitive Variables:

Compute Primitive Variables
---------------------------

We compute the primtive variables from the conserved variables.

-  :math:`\rho, \rho e`: directly copy these from the conserved state
   vector

-  :math:`\ub, A_k, X_k, Y_k`: copy these from the conserved state
   vector, dividing by :math:`\rho`

-  :math:`p,T`: use the EOS.

   First, we use the EOS to ensure :math:`e` is no smaller than :math:`e(\rho,T_{\rm small},X_k)`.
   Then we use the EOS to compute :math:`p,T = p,T(\rho,e,X_k)`.

We also compute the flattening coefficient, :math:`\chi\in[0,1]`, used in
the edge state prediction to further limit slopes near strong shocks.
We use the same flattening procedure described in the the the original
PPM paper :cite:`ppm` and the Flash paper :cite:`flash`.
A flattening coefficient of 1 indicates that no additional limiting
takes place; a flattening coefficient of 0 means we effectively drop
order to a first-order Godunov scheme (this convention is opposite of
that used in the Flash paper). For each cell, we compute the
flattening coefficient for each spatial direction, and choose the
minimum value over all directions. As an example, to compute the
flattening for the x-direction, here are the steps:

#. Define :math:`\zeta`

   .. math:: \zeta_i = \frac{p_{i+1}-p_{i-1}}{\max\left(p_{\rm small},|p_{i+2}-p_{i-2}|\right)}.

#. Define :math:`\tilde\chi`

   .. math:: \tilde\chi_i = \min\left\{1,\max[0,a(\zeta_i - b)]\right\},

   where :math:`a=10` and :math:`b=0.75` are tunable parameters. We are essentially
   setting :math:`\tilde\chi_i=a(\zeta_i-b)`, and then constraining
   :math:`\tilde\chi_i` to lie in the range :math:`[0,1]`. Then, if either
   :math:`u_{i+1}-u_{i-1}<0` or

   .. math:: \frac{p_{i+1}-p_{i-1}}{\min(p_{i+1},p_{i-1})} \le c,

   where :math:`c=1/3` is a tunable parameter, then set :math:`\tilde\chi_i=0`.

#. Define :math:`\chi`

   .. math::

      \chi_i =
      \begin{cases}
      1 - \max(\tilde\chi_i,\tilde\chi_{i-1}) & p_{i+1}-p_{i-1} > 0 \\
      1 - \max(\tilde\chi_i,\tilde\chi_{i+1}) & \text{otherwise}
      \end{cases}.

The following runtime parameters affect the behavior here:

-  castro.use_flattening turns on/off the flattening of parabola
   near shocks (0 or 1; default 1)

Edge State Prediction
---------------------

We wish to compute a left and right state of primitive variables at
each edge to be used as inputs to the Riemann problem. There
are several reconstruction techniques, a piecewise
linear method that follows the description in :cite:`colella:1990`,
the classic PPM limiters :cite:`ppm`, and the new PPM limiters introduced
in :cite:`colellasekora`. The choice of
limiters is determined by castro.ppm_type.

For the new PPM limiters, we have further modified the method
of :cite:`colellasekora` to eliminate sensitivity due to roundoff error
(modifications via personal communication with Colella).

We also use characteristic tracing with corner coupling in 3D, as
described in Miller & Colella (2002) :cite:`millercolella:2002`. We
give full details of the new PPM algorithm, as it has not appeared before
in the literature, and summarize the developments from Miller &
Colella.

The PPM algorithm is used to compute time-centered edge states by
extrapolating the base-time data in space and time. The edge states
are dual-valued, i.e., at each face, there is a left state and a right
state estimate. The spatial extrapolation is one-dimensional, i.e.,
transverse derivatives are ignored. We also use a flattening
procedure to further limit the edge state values. The Miller &
Colella algorithm, which we describe later, incorporates the
transverse terms, and also describes the modifications required for
equations with additional characteristics besides the fluid velocity.
There are four steps to compute these dual-valued edge states (here,
we use :math:`s` to denote an arbitrary scalar from :math:`\Qb`, and we write the
equations in 1D, for simplicity):

-  **Step 1**: Compute :math:`s_{i,+}` and :math:`s_{i,-}`, which are spatial
   interpolations of :math:`s` to the hi and lo side of the face with special
   limiters, respectively. Begin by interpolating :math:`s` to edges using a
   4th-order interpolation in space:

   .. math:: s_{i+\myhalf} = \frac{7}{12}\left(s_{i+1}+s_i\right) - \frac{1}{12}\left(s_{i+2}+s_{i-1}\right).

   Then, if :math:`(s_{i+\myhalf}-s_i)(s_{i+1}-s_{i+\myhalf}) < 0`, we limit
   :math:`s_{i+\myhalf}` a nonlinear combination of approximations to the
   second derivative. The steps are as follows:

   #. Define:

      .. math::

         \begin{align}
         (D^2s)_{i+\myhalf} &= 3\left(s_{i}-2s_{i+\myhalf}+s_{i+1}\right) \\
         (D^2s)_{i+\myhalf,L} &= s_{i-1}-2s_{i}+s_{i+1} \\
         (D^2s)_{i+\myhalf,R} &= s_{i}-2s_{i+1}+s_{i+2}
         \end{align}

   #. Define

      .. math:: s = \text{sign}\left[(D^2s)_{i+\myhalf}\right],

      .. math:: (D^2s)_{i+\myhalf,\text{lim}} = s\max\left\{\min\left[Cs\left|(D^2s)_{i+\myhalf,L}\right|,Cs\left|(D^2s)_{i+\myhalf,R}\right|,s\left|(D^2s)_{i+\myhalf}\right|\right],0\right\},

      where :math:`C=1.25` as used in Colella and Sekora 2009. The limited value
      of :math:`s_{i+\myhalf}` is

      .. math:: s_{i+\myhalf} = \frac{1}{2}\left(s_{i}+s_{i+1}\right) - \frac{1}{6}(D^2s)_{i+\myhalf,\text{lim}}.

   Now we implement an updated implementation of the Colella & Sekora
   algorithm which eliminates sensitivity to roundoff. First we
   need to detect whether a particular cell corresponds to an
   “extremum”. There are two tests.

   -  For the first test, define

      .. math:: \alpha_{i,\pm} = s_{i\pm\myhalf} - s_i.

      If :math:`\alpha_{i,+}\alpha_{i,-} \ge 0`, then we are at an extremum.

   -  We only apply the second test if either
      :math:`|\alpha_{i,\pm}| > 2|\alpha_{i,\mp}|`. If so, we define:

      .. math::

         \begin{align}
         (Ds)_{i,{\rm face},-} &= s_{i-1/2} - s_{i-3/2} \\
         (Ds)_{i,{\rm face},+} &= s_{i+3/2} - s_{i-1/2}
         \end{align}

      .. math:: (Ds)_{i,{\rm face,min}} = \min\left[\left|(Ds)_{i,{\rm face},-}\right|,\left|(Ds)_{i,{\rm face},+}\right|\right].

      .. math::

         \begin{align}
         (Ds)_{i,{\rm cc},-} &= s_{i} - s_{i-1} \\
         (Ds)_{i,{\rm cc},+} &= s_{i+1} - s_{i}
         \end{align}

      .. math:: (Ds)_{i,{\rm cc,min}} = \min\left[\left|(Ds)_{i,{\rm cc},-}\right|,\left|(Ds)_{i,{\rm cc},+}\right|\right].

      If :math:`(Ds)_{i,{\rm face,min}} \ge (Ds)_{i,{\rm cc,min}}`, set
      :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm face},\pm}`. Otherwise, set
      :math:`(Ds)_{i,\pm} = (Ds)_{i,{\rm cc},\pm}`. Finally, we are at an extreumum if
      :math:`(Ds)_{i,+}(Ds)_{i,-} \le 0`.

   Thus concludes the extremum tests. The remaining limiters depend on
   whether we are at an extremum.

   -  If we are at an extremum, we modify :math:`\alpha_{i,\pm}`. First, we
      define

      .. math::

         \begin{align}
         (D^2s)_{i} &= 6(\alpha_{i,+}+\alpha_{i,-}) \\
         (D^2s)_{i,L} &= s_{i-2}-2s_{i-1}+s_{i} \\
         (D^2s)_{i,R} &= s_{i}-2s_{i+1}+s_{i+2} \\
         (D^2s)_{i,C} &= s_{i-1}-2s_{i}+s_{i+1}
         \end{align}

      Then, define

      .. math:: s = \text{sign}\left[(D^2s)_{i}\right],

      .. math:: (D^2s)_{i,\text{lim}} = \max\left\{\min\left[s(D^2s)_{i},Cs\left|(D^2s)_{i,L}\right|,Cs\left|(D^2s)_{i,R}\right|,Cs\left|(D^2s)_{i,C}\right|\right],0\right\}.

      Then,

      .. math:: \alpha_{i,\pm} = \frac{\alpha_{i,\pm}(D^2s)_{i,\text{lim}}}{\max\left[(D^2s)_{i},1\times 10^{-10}\right]}

   -  If we are not at an extremum and 
      :math:`|\alpha_{i,\pm}| > 2|\alpha_{i,\mp}|`, then define

      .. math:: s = \text{sign}(\alpha_{i,\mp})

      .. math:: \delta\mathcal{I}_{\text{ext}} = \frac{-\alpha_{i,\pm}^2}{4\left(\alpha_{j,+}+\alpha_{j,-}\right)},

      .. math:: \delta s = s_{i\mp 1} - s_i,

      If :math:`s\delta\mathcal{I}_{\text{ext}} \ge s\delta s`, then we perform
      the following test. If :math:`s\delta s - \alpha_{i,\mp} \ge 1\times
      10^{-10}`, then

      .. math:: \alpha_{i,\pm} =  -2\delta s - 2s\left[(\delta s)^2 - \delta s \alpha_{i,\mp}\right]^{\myhalf}

      otherwise,

      .. math:: \alpha_{i,\pm} =  -2\alpha_{i,\mp}

   Finally, :math:`s_{i,\pm} = s_i + \alpha_{i,\pm}`.

-  **Step 2**: Construct a quadratic profile using :math:`s_{i,-},s_i`,
   and :math:`s_{i,+}`.

   .. math::
      s_i^I(x) = s_{i,-} + \xi\left[s_{i,+} - s_{i,-} + s_{6,i}(1-\xi)\right],
      :label: Quadratic Interp

   .. math:: s_6 = 6s_{i} - 3\left(s_{i,-}+s_{i,+}\right),

   .. math:: \xi = \frac{x - ih}{h}, ~ 0 \le \xi \le 1.

-  | **Step 3:** Integrate quadratic profiles. We are essentially
     computing the average value swept out by the quadratic profile
     across the face assuming the profile is moving at a speed
     :math:`\lambda_k`.
   | Define the following integrals, where :math:`\sigma_k =
       |\lambda_k|\Delta t/h`:

     .. math::

        \begin{align}
        \mathcal{I}^{(k)}_{+}(s_i) &= \frac{1}{\sigma_k h}\int_{(i+\myhalf)h-\sigma_k h}^{(i+\myhalf)h}s_i^I(x)dx \\
        \mathcal{I}^{(k)}_{-}(s_i) &= \frac{1}{\sigma_k h}\int_{(i-\myhalf)h}^{(i-\myhalf)h+\sigma_k h}s_i^I(x)dx
        \end{align}

     Plugging in :eq:`Quadratic Interp` gives:

     .. math::

        \begin{align}
        \mathcal{I}^{(k)}_{+}(s_i) &= s_{i,+} - \frac{\sigma_k}{2}\left[s_{i,+}-s_{i,-}-\left(1-\frac{2}{3}\sigma_k\right)s_{6,i}\right], \\
        \mathcal{I}^{(k)}_{-}(s_i) &= s_{i,-} + \frac{\sigma_k}{2}\left[s_{i,+}-s_{i,-}+\left(1-\frac{2}{3}\sigma_k\right)s_{6,i}\right].
        \end{align}

-  **Step 4:** Obtain 1D edge states by performing a 1D
   extrapolation to get left and right edge states. Note that we
   include an explicit source term contribution.

   .. math::

      \begin{align}
      s_{L,i+\myhalf} &= s_i - \chi_i\sum_{k:\lambda_k \ge 0}\lb_k\cdot\left[s_i-\mathcal{I}^{(k)}_{+}(s_i)\right]\rb_k + \frac{\dt}{2}S_i^n, \\
      s_{R,i-\myhalf} &= s_i - \chi_i\sum_{k:\lambda_k < 0}\lb_k\cdot\left[s_i-\mathcal{I}^{(k)}_{-}(s_i)\right]\rb_k + \frac{\dt}{2}S_i^n.
      \end{align}

   Here, :math:`\rb_k` is the :math:`k^{\rm th}` right column eigenvector of
   :math:`\Rb(\Ab_d)` and :math:`\lb_k` is the :math:`k^{\rm th}` left row eigenvector lf
   :math:`\Lb(\Ab_d)`. The flattening coefficient is :math:`\chi_i`.

In order to add the transverse terms in an spatial operator unsplit
framework, the details follow exactly as given in Section 4.2.1 in
Miller & Colella, except for the details of the Riemann solver,
which are given below.

.. index:: castro.ppm_type

For the reconstruction of the interface states, the following apply:

-  ``castro.ppm_type`` : use piecewise linear vs PPM algorithm (0 or 1;
   default: 1).  A value of 1 is the standard piecewise parabolic
   reconstruction.

-  ``castro.ppm_temp_fix`` does various attempts to use the
   temperature in the reconstruction of the interface states.
   See :ref:`sec-ppm_temp_fix` for an explanation of the allowed options.

The interface states are corrected with information from the
transverse directions to make this a second-order update. These
transverse directions involve separate Riemann solves. Sometimes, the
update to the interface state from the transverse directions can make
the state ill-posed. There are several parameters that help fix this:

-  ``castro.transverse_use_eos`` : If this is 1, then we call
   the equation of state on the interface, using :math:`\rho`, :math:`e`, and
   :math:`X_k`, to get the interface pressure. This should result in a
   thermodynamically consistent interface state.

-  ``castro.transverse_reset_density`` : If the transverse
   corrections result in a negative density on the interface, then we
   reset all of the interface states to their values before the
   transverse corrections.

-  ``castro.transverse_reset_rhoe`` : The transverse updates operate
   on the conserved state. Usually, we construct the interface
   :math:`(\rho e)` in the transverse update from total energy and the
   kinetic energy, however, if the interface :math:`(rho e)` is negative,
   and ``transverse_reset_rhoe`` = 1, then we explicitly
   discretize an equation for the evolution of :math:`(\rho e)`, including
   its transverse update.

Riemann Problem
---------------

Castro has three main options for the Riemann solver—the
Colella & Glaz solver :cite:`colglaz` (the same solver used
by Flash), a simpler solver described in an unpublished
manuscript by Colella, Glaz, & Ferguson, and an HLLC
solver. The first two are both
two-shock approximate solvers, but differ in how they approximate
the thermodynamics in the “star” region.

.. index:: castro.riemann_speed_limit

.. note::

   These Riemann solvers are for Newtonian hydrodynamics, however, we enforce
   that the interface velocity cannot exceed the speed of light in both the
   Colella & Glaz and Colella, Glaz, & Ferguson solvers.  This excessive speed
   usually is a sign of low density regions and density resets or the flux limiter
   kicking in.  This behavior can be changed with the ``castro.riemann_speed_limit``
   parameter.

Inputs from the edge state prediction are :math:`\rho_{L/R}, u_{L/R},
v_{L/R}, p_{L/R}`, and :math:`(\rho e)_{L/R}` (:math:`v` represents all of the
transverse velocity components). We also compute :math:`\Gamma \equiv d\log
p / d\log \rho |_s` at cell centers and copy these to edges directly
to get the left and right states, :math:`\Gamma_{L/R}`. We also define
:math:`c_{\rm avg}` as a face-centered value that is the average of the
neighboring cell-centered values of :math:`c`. We have also computed
:math:`\rho_{\rm small}, p_{\rm small}`, and :math:`c_{\rm small}` using
cell-centered data.

Here are the steps. First, define 
:math:`(\rho c)_{\rm small} = \rho_{\rm small}c_{\rm small}`. Then, define:

.. math:: (\rho c)_{L/R} = \max\left[(\rho c)_{\rm small},\left|\Gamma_{L/R},p_{L/R},\rho_{L/R}\right|\right].

Define star states:

.. math:: p^* = \max\left[p_{\rm small},\frac{\left[(\rho c)_L p_R + (\rho c)_R p_L\right] + (\rho c)_L(\rho c)_R(u_L-u_R)}{(\rho c)_L + (\rho c)_R}\right],

.. math:: u^* = \frac{\left[(\rho c)_L u_L + (\rho c)_R u_R\right]+ (p_L - p_R)}{(\rho c)_L + (\rho c)_R}.

If :math:`u^* \ge 0` then define :math:`\rho_0, u_0, p_0, (\rho e)_0` and :math:`\Gamma_0` to be the left state. Otherwise, define them to be the right state. Then, set

.. math:: \rho_0 = \max(\rho_{\rm small},\rho_0),

and define

.. math:: c_0 = \max\left(c_{\rm small},\sqrt{\frac{\Gamma_0 p_0}{\rho_0}}\right),

.. math:: \rho^* = \rho_0 + \frac{p^* - p_0}{c_0^2},

.. math:: (\rho e)^* = (\rho e)_0 + (p^* - p_0)\frac{(\rho e)_0 + p_0}{\rho_0 c_0^2},

.. math:: c^* = \max\left(c_{\rm small},\sqrt{\left|\frac{\Gamma_0 p^*}{\rho^*}\right|}\right)

Then,

.. math::

   \begin{align}
   c_{\rm out} &= c_0 - {\rm sign}(u^*)u_0, \\
   c_{\rm in} &= c^* - {\rm sign}(u^*)u^*, \\
   c_{\rm shock} &= \frac{c_{\rm in} + c_{\rm out}}{2}.
   \end{align}

If :math:`p^* - p_0 \ge 0`, then :math:`c_{\rm in} = c_{\rm out} = c_{\rm shock}`.
Then, if :math:`c_{\rm out} = c_{\rm in}`, we define :math:`c_{\rm temp} =
\epsilon c_{\rm avg}`. Otherwise, :math:`c_{\rm temp} = c_{\rm out} -
c_{\rm in}`. We define the fraction

.. math:: f = \half\left[1 + \frac{c_{\rm out} + c_{\rm in}}{c_{\rm temp}}\right],

and constrain :math:`f` to lie in the range :math:`f\in[0,1]`.

To get the final “Godunov” state, for the transverse velocity, we
upwind based on :math:`u^*`.

.. math::

   v_{\rm gdnv} =
   \begin{cases}
   v_L, & u^* \ge 0 \\
   v_R, & {\rm otherwise}
   \end{cases}.

Then, define

.. math::

   \begin{align}
   \rho_{\rm gdnv} &= f\rho^* + (1-f)\rho_0, \\
   u_{\rm gdnv} &= f u^* + (1-f)u_0, \\
   p_{\rm gdnv} &= f p^* + (1-f)p_0, \\
   (\rho e)_{\rm gdnv} &=& f(\rho e)^* + (1-f)(\rho e)_0.
   \end{align}

Finally, if :math:`c_{\rm out} < 0`, set 
:math:`\rho_{\rm gdnv}=\rho_0, u_{\rm gdnv}=u_0, p_{\rm gdnv}=p_0`, and 
:math:`(\rho e)_{\rm gdnv}=(\rho e)_0`.
If :math:`c_{\rm in}\ge 0`, set :math:`\rho_{\rm gdnv}=\rho^*, u_{\rm gdnv}=u^*,
p_{\rm gdnv}=p^*`, and :math:`(\rho e)_{\rm gdnv}=(\rho e)^*`.

If instead the Colella & Glaz solver is used, then we define

.. math:: \gamma \equiv \frac{p}{\rho e} + 1

on each side of the interface and follow the rest of the algorithm as
described in the original paper.

For the construction of the fluxes in the Riemann solver, the following
parameters apply:

-  ``castro.riemann_solver``: this can be one of the following values:

   -  0: the Colella, Glaz, & Ferguson solver.

   -  1: the Colella & Glaz solver

   -  2: the HLLC solver. Note: this should only be used with Cartesian
      geometries because it relies on the pressure term being part of the flux
      in the momentum equation.

   The default is to use the solver based on an unpublished Colella,
   Glaz, & Ferguson manuscript (it also appears in :cite:`pember:1996`),
   as described in the original Castro paper :cite:`castro_I`.

   The Colella & Glaz solver is iterative, and two runtime parameters are used
   to control its behavior:

   -  ``castro.cg_maxiter`` : number of iterations for CG algorithm
      (Integer; default: 12)

   -  ``castro.cg_tol`` : tolerance for CG solver when solving
      for the “star” state (Real; default: 1.0e-5)

   -  ``castro.cg_blend`` : this controls what happens if the root
      finding in the CG solver fails. There is a nonlinear equation to find
      the pressure in the *star* region from the jump conditions for a
      shock (this is the two-shock approximation—the left and right states
      are linked to the star region each by a shock). The default root
      finding algorithm is a secant method, but this can sometimes fail.

      The options here are:

      -  0 : do nothing. The pressure from each iteration is
         printed and the code aborts with a failure

      -  1 : revert to the original guess for p-star and carry
         through on the remainder of the Riemann solve. This is almost like
         dropping down to the CGF solver. The p-star used is very approximate.

      -  2 : switch to bisection and do an additional cg_maxiter
         iterations to find the root. Sometimes this can work where the
         secant method fails.

-  ``castro.hybrid_riemann`` : switch to an HLL Riemann solver when we are
   in a zone with a shock (0 or 1; default 0)

   This eliminates an odd-even decoupling issue (see the oddeven
   problem). Note, this cannot be used with the HLLC solver.

Compute Fluxes and Update
-------------------------

Compute the fluxes as a function of the primitive variables, and then
advance the solution:

.. math:: \Ub^{n+1} = \Ub^n - \dt\nabla\cdot\Fb^\nph + \dt\Sb^n.

Again, note that since the source term is not time centered, this is
not a second-order method. After the advective update, we correct the
solution, effectively time-centering the source term.

.. _sec-ppm_temp_fix:

Temperature Fixes
=================

.. index:: castro.ppm_temp_fix

There are a number of experimental options for improving the behavior
of the temperature in the reconstruction and interface state
prediction. The options are controlled by ``castro.ppm_temp_fix``,
which takes values:

  * 0: the default method—temperature is not considered, and we do
    reconstruction and characteristic tracing on :math:`\rho, u, p,
    (\rho e)`.

  * 1: do parabolic reconstruction on :math:`T`, giving
    :math:`\mathcal{I}_{+}^{(k)}(T_i)`. We then derive the pressure and
    internal energy (gas portion) via the equation of state as:

    .. math::

      \begin{align}
            \mathcal{I}_{+}^{(k)}(p_i) &= p(\mathcal{I}_{+}^{(k)}(\rho_i), \mathcal{I}_{+}^{(k)}(T_i)) \\
            \mathcal{I}_{+}^{(k)}((\rho e)_i) &= (\rho e)(\mathcal{I}_{+}^{(k)}(\rho_i), \mathcal{I}_{+}^{(k)}(T_i))
          \end{align}

    The remainder of the hydrodynamics algorithm then proceeds unchanged.

  * 2: on entering the Riemann solver, we recompute the thermodynamics
    on the interfaces to ensure that they are all consistent. This is
    done by taking the interface values of :math:`\rho`, :math:`e`,
    :math:`X_k`, and computing the corresponding pressure, :math:`p`
    from this.


Resets
======

Density Resets
--------------

Need to document density_reset_method

.. _app:hydro:flux_limiting:

Flux Limiting
-------------

Multi-dimensional hydrodynamic simulations often have numerical
artifacts that result from the sharp density gradients. A somewhat
common issue, especially at low resolution, is negative densities that
occur as a result of a hydro update. Castro contains a prescription
for dealing with negative densities, that resets the negative density
to be similar to nearby zones. Various choices exist for how to do
this, such as resetting it to the original zone density before the
update or resetting it to some linear combination of the density of
nearby zones. The reset is problematic because the strategy is not
unique and no choice is clearly better than the rest in all
cases. Additionally, it is not specified at all how to reset momenta
in such a case. Consequently, we desired to improve the situation by
limiting fluxes such that negative densities could not occur, so that
such a reset would in practice always be avoided. Our solution
implements the positivity-preserving method of :cite:`hu:2013`. This
behavior is controlled by
castro.limit_fluxes_on_small_dens.

A hydrodynamical update to a zone can be broken down into an update
over every face of the zone where a flux crosses the face over the
timestep. The central insight of the positivity-preserving method is
that if the update over every face is positivity-preserving, then the
total update must be positivity-preserving as well. To guarantee
positivity preservation at the zone edge :math:`{\rm i}+1/2`, the flux
:math:`\mathbf{F}^{n+1/2}_{{\rm i}+1/2}` at that face is modified to become:

.. math:: \mathbf{F}^{n+1/2}_{{\rm i}+1/2} \rightarrow \theta_{{\rm i}+1/2} \mathbf{F}^{n+1/2}_{{\rm i}+1/2} + (1 - \theta_{{\rm i}+1/2}) \mathbf{F}^{LF}_{{\rm i}+1/2}, \label{eq:limited_flux}

where :math:`0 \leq \theta_{{\rm i}+1/2} \leq 1` is a scalar, and :math:`\mathbf{F}^{LF}_{{\rm i}+1/2}` is the Lax-Friedrichs flux,

.. math:: \mathbf{F}^{LF}_{{\rm i}+1/2} = \frac{1}{2}\left[\mathbf{F}^{n}_{{\rm i}} + \mathbf{F}^{n}_{{\rm i}+1} + \text{CFL}\frac{\Delta x}{\Delta t} \frac{1}{\alpha}\left(\mathbf{U}^{n}_{{\rm i}} - \mathbf{U}^{n}_{{\rm i}+1}\right)\right],

where :math:`0 < \text{CFL} < 1` is the CFL safety factor (the method is
guaranteed to preserve positivity as long as :math:`\text{CFL} < 1/2`), and
:math:`\alpha` is a scalar that ensures multi-dimensional correctness
(:math:`\alpha = 1` in 1D, :math:`1/2` in 2D, :math:`1/3` in 3D). 
:math:`\mathbf{F}_{{\rm i}}` is the flux of material evaluated at the zone center 
:math:`{\rm i}` using the cell-centered quantities :math:`\mathbf{U}`. The scalar
:math:`\theta_{{\rm i}+1/2}` is chosen at every interface by calculating the
update that would be obtained from , setting
the density component equal to a value just larger than the density floor,
castro.small_dens, and solving
for the value of :math:`\theta` at the interface that makes the equality
hold. In regions where the density is not at risk of going negative,
:math:`\theta \approx 1` and the original hydrodynamic update is recovered.
Further discussion, including a proof of the method, a description of
multi-dimensional effects, and test verification problems, can be
found in :cite:`hu:2013`.


Hybrid Momentum
===============

Castro implements the hybrid momentum scheme of :cite:`byerly:2014`.
In particular, this switches from using the Cartesian momenta,
:math:`(\rho u)`, :math:`(\rho v)`, and :math:`(\rho w)`, to a
cylindrical momentum set, :math:`(\rho v_R)`, :math:`(\rho R v_\phi)`,
and :math:`(\rho v_z)`.  This latter component is identical to the
Cartesian value.  We translate between these sets of momentum throughout the code,
ultimately doing the conservative update in terms of the cylindrical momentum.  Additional
source terms appear in this formulation, which are written out in :cite:`byerly:2014`.

The ``rotating_torus`` problem gives a good test for this.  This problem
originated with :cite:`papaloizoupringle`.  The
problem is initialized as a torus with constant specific angular
momentum, as shown below:

.. figure:: rotating_torus_00000_density.png
   :alt: rotating torus initial density

   Initial density (log scale) for the ``rotating_torus`` problem with
   :math:`64^3` zones.

For the standard hydrodynamics algorithm, the torus gets disrupted and
spreads out into a disk:

.. figure:: rotating_torus_00200_density.png
   :alt: rotating torus normal hydro

   Density (log scale) for the ``rotating_torus`` problem after 200
   timesteps, using :math:`64^3` zones.  Notice that the initial torus
   has become disrupted into a disk.

The hybrid momentum algorithm is enabled by setting::

   USE_HYBRID_MOMENTUM = TRUE

in your ``GNUmakefile``.  With this enabled, we see that the torus remains intact:

.. figure:: rotating_torus_hybrid_00200_density.png
   :alt: rotating torus with hybrid momentum

   Density (log scale) for the ``rotating_torus`` problem after 200
   timesteps with the hybrid momentum algorithm, using :math:`64^3`
   zones.  With this angular-momentum preserving scheme we see that
   the initial torus is largely intact.

