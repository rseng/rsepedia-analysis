Dymos:  Open Source Optimization of Dynamic Multidisciplinary Systems
=====================================================================

[![Dymos Tests](https://github.com/OpenMDAO/dymos/actions/workflows/dymos_tests_workflow.yml/badge.svg)](https://github.com/OpenMDAO/dymos/actions/workflows/dymos_tests_workflow.yml) [![Coverage Status](https://coveralls.io/repos/github/OpenMDAO/dymos/badge.svg?branch=master&t=dJxu2Q)](https://coveralls.io/github/OpenMDAO/dymos?branch=master)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.02809/status.svg)](https://doi.org/10.21105/joss.02809)



Dymos is a framework for the simulation and optimization of dynamical systems within the OpenMDAO Multidisciplinary Analysis and Optimization environment.
Dymos leverages implicit and explicit simulation techniques to simulate generic dynamic systems of arbitary complexity.

The software has two primary objectives:
-   Provide a generic ODE integration interface that allows for the analysis of dynamical systems.
-   Allow the user to solve optimal control problems involving dynamical multidisciplinary systems.

Installation
------------

The default installation of the developmental version of Dymos will install the minimum number of prerequisites:

```
python -m pip install dymos
```

More advanced installation instructions are available [here](https://openmdao.github.io/dymos/installation.html).

Citation
--------

See our [overview paper](https://joss.theoj.org/papers/10.21105/joss.02809) in the Journal of Open Source Software

If you use Dymos in your work, please cite: 
```
@article{Falck2021,
  doi = {10.21105/joss.02809},
  url = {https://doi.org/10.21105/joss.02809},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {59},
  pages = {2809},
  author = {Robert Falck and Justin S. Gray and Kaushik Ponnapalli and Ted Wright},
  title = {dymos: A Python package for optimal control of multidisciplinary systems},
  journal = {Journal of Open Source Software}
}
```



Documentation
-------------

Online documentation is available at [https://openmdao.github.io/dymos/](https://openmdao.github.io/dymos/)

Defining Ordinary Differential Equations
----------------------------------------

The first step in simulating or optimizing a dynamical system is to define the ordinary
differential equations to be integrated.  The user first builds an OpenMDAO model which has outputs
that provide the rates of the state variables.  This model can be an OpenMDAO model of arbitrary
complexity, including nested groups and components, layers of nonlinear solvers, etc.

Dymos solutions are constructed of one or more _Phases_.
When setting up a phase, we add state variables, dynamic controls, and parameters,
tell Dymos how the value of each should be connected to the ODE system, and tell Dymos
the variable paths in the system that contain the rates of our state variables that are to be
integrated.

Integrating Ordinary Differential Equations
-------------------------------------------

Dymos's solver-based pseudspectral transcriptions
provide the ability to numerically integrate the ODE system it is given.
Used in an optimal control context, these provide a shooting method in
which each iteration provides a physically viable trajectory.

Pseudospectral Methods
----------------------

Dymos currently supports the Radau Pseudospectral Method and high-order
Gauss-Lobatto transcriptions.  These implicit techniques rely on the
optimizer to impose "defect" constraints which enforce the physical
accuracy of the resulting trajectories.  To verify the physical
accuracy of the solutions, Dymos can explicitly integrate them using
variable-step methods.

Solving Optimal Control Problems
--------------------------------

Dymos uses the concept of _Phases_ to support optimal control of dynamical systems.
Users connect one or more Phases to construct trajectories.
Each Phase can have its own:

-   Optimal Control Transcription (Gauss-Lobatto or Radau Pseudospectral)
-   Equations of motion
-   Boundary and path constraints

Dymos Phases and Trajectories are ultimately just OpenMDAO Groups that can exist in
a problem along with numerous other models, allowing for the simultaneous
optimization of systems and dynamics.

```python
import numpy as np
import openmdao.api as om
import dymos as dm
import matplotlib.pyplot as plt


class BrachistochroneEOM(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int)

    def setup(self):
        nn = self.options['num_nodes']

        # Inputs
        self.add_input('v',
                       shape=(nn,),
                       desc='velocity',
                       units='m/s')

        self.add_input('theta',
                       shape=(nn,),
                       desc='angle of wire',
                       units='rad')

        self.add_output('xdot',
                        shape=(nn,),
                        desc='velocity component in x',
                        units='m/s')

        self.add_output('ydot',
                        shape=(nn,),
                        desc='velocity component in y',
                        units='m/s')

        self.add_output('vdot',
                        val=np.zeros(nn),
                        desc='acceleration magnitude',
                        units='m/s**2')

        self.add_output('check',
                        val=np.zeros(nn),
                        desc='A check on the solution: v/sin(theta) = constant',
                        units='m/s')

        # Use OpenMDAO's ability to automatically determine a sparse "coloring" of the jacobian
        # for this ODE component.
        self.declare_coloring(wrt='*', method='cs')

    def compute(self, inputs, outputs):
        theta = inputs['theta']
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        v = inputs['v']

        outputs['vdot'] = 9.80665 * cos_theta
        outputs['xdot'] = v * sin_theta
        outputs['ydot'] = -v * cos_theta
        outputs['check'] = v / sin_theta


# Define the OpenMDAO problem
p = om.Problem(model=om.Group())

# Define a Trajectory object
traj = dm.Trajectory()
p.model.add_subsystem('traj', subsys=traj)

# Define a Dymos Phase object with GaussLobatto Transcription
phase = dm.Phase(ode_class=BrachistochroneEOM,
                 transcription=dm.Radau(num_segments=10, order=3))
traj.add_phase(name='phase0', phase=phase)

# Set the time options
phase.set_time_options(fix_initial=True,
                       duration_bounds=(0.5, 10.0))

# Set the state options
phase.add_state('x', rate_source='xdot',
                fix_initial=True, fix_final=True)
phase.add_state('y', rate_source='ydot',
                fix_initial=True, fix_final=True)
phase.add_state('v', rate_source='vdot',
                fix_initial=True, fix_final=False)

# Define theta as a control.
phase.add_control(name='theta', units='rad',
                  lower=0, upper=np.pi)

# Minimize final time.
phase.add_objective('time', loc='final')

# Set the driver.
p.driver = om.ScipyOptimizeDriver()

# Allow OpenMDAO to automatically determine our sparsity pattern
# for TOTAL derivatives. Doing so can significantly speed up the
# execution of dymos.
p.driver.declare_coloring()

# Setup the problem
p.setup()

# Now that the OpenMDAO problem is setup, we can guess the
# values of time, states, and controls.
p.set_val('traj.phase0.t_duration', 2.0)

# States and controls here use a linearly interpolated
# initial guess along the trajectory.
p.set_val('traj.phase0.states:x',
          phase.interpolate(ys=[0, 10], nodes='state_input'),
          units='m')

p.set_val('traj.phase0.states:y',
          phase.interpolate(ys=[10, 5], nodes='state_input'),
          units='m')

p.set_val('traj.phase0.states:v',
          phase.interpolate(ys=[0, 5], nodes='state_input'),
          units='m/s')

p.set_val('traj.phase0.controls:theta',
          phase.interpolate(ys=[5, 45], nodes='control_input'),
          units='deg')

# Use Dymos' run_problem method to run the driver, simulate the results,
# and record the results to 'dymos_solution.db' and 'dymos_simulation.db'.
dm.run_problem(p, simulate=True)

# Load the solution and simulation files.
sol_case = om.CaseReader('dymos_solution.db').get_case('final')
sim_case = om.CaseReader('dymos_simulation.db').get_case('final')

# Plot the results
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12, 4.5))

axes[0].plot(sol_case.get_val('traj.phase0.timeseries.states:x'),
             sol_case.get_val('traj.phase0.timeseries.states:y'),
             'ro', label='solution')

axes[0].plot(sim_case.get_val('traj.phase0.timeseries.states:x'),
             sim_case.get_val('traj.phase0.timeseries.states:y'),
             'b-', label='simulation')

axes[0].set_xlabel('x (m)')
axes[0].set_ylabel('y (m/s)')
axes[0].legend()
axes[0].grid()

axes[1].plot(sol_case.get_val('traj.phase0.timeseries.time'),
             sol_case.get_val('traj.phase0.timeseries.controls:theta',
                              units='deg'),
             'ro', label='solution')

axes[1].plot(sim_case.get_val('traj.phase0.timeseries.time'),
             sim_case.get_val('traj.phase0.timeseries.controls:theta',
                              units='deg'),
             'b-', label='simulation')

axes[1].set_xlabel('time (s)')
axes[1].set_ylabel(r'$\theta$ (deg)')
axes[1].legend()
axes[1].grid()

plt.show()
```

![Brachistochrone Solution](brachistochroneSolution.png "Brachistochrone Solution")
*******************************
# Release Notes for Dymos 1.4.0

January 05, 2022

This is version 1.4.0 of Dymos.
It includes a fix for grid refinement and simulation with parameters, some minor performance improvements, and various documentation updates.

## Enhancements

* Disabled check_partials of Dymos core components by default. [#686](https://github.com/OpenMDAO/dymos/pull/686)
* Added performance improvements for PseudospectralTimeseriesOutputComp. [#688](https://github.com/OpenMDAO/dymos/pull/688)
* Added a warning if bounds are applied to a state during solve_segments. [#691](https://github.com/OpenMDAO/dymos/pull/691)
* Expanded the included 1976 standard atmosphere model to cover -15000 ft to 250000 ft and doc cleanup. [#699](https://github.com/OpenMDAO/dymos/pull/699)
* Added the Bryson-Denham example problem and other various documentation improvements. [#700](https://github.com/OpenMDAO/dymos/pull/700) [#702](https://github.com/OpenMDAO/dymos/pull/702) [#704](https://github.com/OpenMDAO/dymos/pull/704)
* Changed deprecated 'value' metadata usages to 'val'. [#706](https://github.com/OpenMDAO/dymos/pull/706)

## Bug Fixes

* Fixed a bug that could result in incorrect parameter values in simulation. [#684](https://github.com/OpenMDAO/dymos/pull/684)
* Fixed example that doesn't work in google collab. [#690](https://github.com/OpenMDAO/dymos/pull/690)
* Fixed a bug in grid refinement error estimation for state rates not from ODE. [#703](https://github.com/OpenMDAO/dymos/pull/703)

## Miscellaneous

* Added [notebooks] spec when installing with pip. [#695](https://github.com/OpenMDAO/dymos/pull/695)
* Added CI matrix entry for testing without pyoptsparse. [#697](https://github.com/OpenMDAO/dymos/pull/697)

*******************************
# Release Notes for Dymos 1.3.0

November 19, 2021

This is version 1.3.0 of Dymos.

This release of Dymos introduces an ExplicitShooting transcription that provides an explicit Runge-Kutta integration of the ODE across a phase.
This transcription is currently limited to fixed-step RK methods (RK4 being the default).
Timeseries outputs are provided at the start/end of each segment in the phase.
This is similar to the solve-segments capability in the collocation transcriptions, but fixed-step will provide _an_ answer (albeit inaccurate) across the integration rather than failing to converge if the dynamics become highly nonlinear.

## Enhancements

* Added path constraints to the explicit shooting transcription. [#659](https://github.com/OpenMDAO/dymos/pull/659)
* Added control continuity enforcement to ExplicitShooting transcription, and refactored continuity components in general. [#660](https://github.com/OpenMDAO/dymos/pull/660)
* Added indication of fixed variables to linkage report. [#662](https://github.com/OpenMDAO/dymos/pull/662)
* Replaced the tensordot in the compute method of timeseries_output_comp with a regular dot product to remove a performance bottleneck.  [#665](https://github.com/OpenMDAO/dymos/pull/665)
* Added constraint report to summarize boundary and path constraints for each phase of a trajectory. [#666](https://github.com/OpenMDAO/dymos/pull/666)
* Added ExplicitShooting to transcriptions [#669](https://github.com/OpenMDAO/dymos/pull/669)
* Significantly improved speed of ExplicitShooting [#670](https://github.com/OpenMDAO/dymos/pull/670)
* Added Radau, BDF and LSODA as options for scipy's integration method when using simulate [#675](https://github.com/OpenMDAO/dymos/pull/675)

## Bug Fixes

* Removed solver for connected linkages. Its only needed for solve_segments. [#668](https://github.com/OpenMDAO/dymos/pull/668)
* Changed default value of units in Trajectory.add_parameter to _unspecified. [#673](https://github.com/OpenMDAO/dymos/pull/673)
* Added fix to allow parameters with static_targets=True to work with ExplicitShooting [#679](https://github.com/OpenMDAO/dymos/pull/679)
* Fixed formatting in the constraint report [#680](https://github.com/OpenMDAO/dymos/pull/680)

## Miscellaneous

* None

*******************************
# Release Notes for Dymos 1.2.0

October 12, 2021

This is version 1.2.0 of Dymos.

The release provides compatibility with OpenMDAO >=3.13.0 and adds
some performance improvements.

While we are beginning to bring a true explicit shooting capability
to Dymos, those capabilities are not fully filled out as of this release.

## Backwards Incompatible API Changes & Deprecations

* Dymos 1.2.0 requires OpenMDAO >= 3.13.0, due to changes in the way indices are specified in OpenMDAO.

## Enhancements

* Update run_problem.py to return success state [#634](https://github.com/OpenMDAO/dymos/pull/634)
* Added an experimental explicit shooting transcription to dymos [#637](https://github.com/OpenMDAO/dymos/pull/637)
* Added control rates and their derivatives when using ExplicitShooting. [#645](https://github.com/OpenMDAO/dymos/pull/645)
* Rewrite of the USatm1976Comp to use pre-computed akima coefficients for interpolation. It is now complex-safe and considerably faster. [#652](https://github.com/OpenMDAO/dymos/pull/652)
* Allow addition of ODE outputs to ExplicitShooting timeseries [#654](https://github.com/OpenMDAO/dymos/pull/654)

## Bug Fixes

* Fixed an issue where simulation was not working when running under MPI in run_problem [#628](https://github.com/OpenMDAO/dymos/pull/628)
* Added a better error message when simulate fails due to the inability to find a good step size. [#630](https://github.com/OpenMDAO/dymos/pull/630)
* Fixes a bug where t_initial_targets and t_duration_targets would not work if input_initial or input_duration were True, respectively. [#656](https://github.com/OpenMDAO/dymos/pull/656)
* Fix to eliminate warning messages related to the recent indexing update to OpenMDAO. [Requires OpenMDAO >= 3.12.0] [#636](https://github.com/OpenMDAO/dymos/pull/636)
* Removed exceptions introduced in OpenMDAO PR [#2279](https://github.com/OpenMDAO/OpenMDAO/pull/2279). [#653](https://github.com/OpenMDAO/dymos/pull/653)

## Miscellaneous

* Fixed issue in executable notebooks. [#631](https://github.com/OpenMDAO/dymos/pull/631)
* Updated CI matrix to test against latest release and development versions of OpenMDAO. [#638](https://github.com/OpenMDAO/dymos/pull/638)

*******************************
# Release Notes for Dymos 1.1.0

July 22, 2021

This is version 1.1.0 of Dymos.
The release provides compatibility with OpenMDAO >=3.10.0, updates the
documentation to JupyterBook, and adds a few new features.

## Backwards Incompatible API Changes & Deprecations

* Removed vectorize_derivs option from phase objectives due to OpenMDAO update. [#605](https://github.com/OpenMDAO/dymos/pull/605)
* The dynamic argument on add_parameter has been removed. A new argument static_target has been added which has opposite meaning of dynamic. [#591](https://github.com/OpenMDAO/dymos/pull/591)
* Updated phase.interpolate to automatically detect variable type, renamed to phase.interp.  Old version is deprecated. [#592](https://github.com/OpenMDAO/dymos/pull/592)

## Enhancements

* Documentation updated to JupyterBook. [#611](https://github.com/OpenMDAO/dymos/pull/611) [#613](https://github.com/OpenMDAO/dymos/pull/613) [#614](https://github.com/OpenMDAO/dymos/pull/614) [#615](https://github.com/OpenMDAO/dymos/pull/615) [#616](https://github.com/OpenMDAO/dymos/pull/616) [#618](https://github.com/OpenMDAO/dymos/pull/618)
* simulate_options are now stored within each Phase. [#610](https://github.com/OpenMDAO/dymos/pull/610)
* Removed vectorize_derivs option from phase objectives due to OpenMDAO update. [#605](https://github.com/OpenMDAO/dymos/pull/605)
* Updated dymos to handle the new OpenMDAO distributed I/O approach [#597](https://github.com/OpenMDAO/dymos/pull/597)
* The dynamic argument on add_parameter has been removed. A new argument static_target has been added which has opposite meaning of dynamic. [#591](https://github.com/OpenMDAO/dymos/pull/591)
* Updated phase.interpolate to automatically detect variable type, renamed to phase.interp.  Old version is deprecated. [#592](https://github.com/OpenMDAO/dymos/pull/592)

## Bug Fixes

* Fixed an issue with units in linkage constraints. [#620](https://github.com/OpenMDAO/dymos/pull/620)
* Fix for key error when performing order reduction under hp adaptive refinement. [#590](https://github.com/OpenMDAO/dymos/pull/590)
* Parameters in the ODE system now respect both dynamic=True and False. [#581](https://github.com/OpenMDAO/dymos/pull/581)

## Miscellaneous

* Added missing example Brachistochrone with upstream initial and duration states. [#623](https://github.com/OpenMDAO/dymos/pull/623)
* Added require_pyoptsparse to all tests that use pyOptSparseDriver [#624](https://github.com/OpenMDAO/dymos/pull/624)
* Added a couple of fixes to examples docs [#622](https://github.com/OpenMDAO/dymos/pull/622)
* Fixed some typos in 'Getting Started' section of the docs [#621](https://github.com/OpenMDAO/dymos/pull/621)
* Install coveralls from pypi in github workflow [#601](https://github.com/OpenMDAO/dymos/pull/601)
* Change base_dir arg to coveralls [#600](https://github.com/OpenMDAO/dymos/pull/600)
* Fixed minor typos in docs. [#599](https://github.com/OpenMDAO/dymos/pull/599)
* Removed require_pyoptsparse and moved it to OpenMDAO. [#595](https://github.com/OpenMDAO/dymos/pull/595)
* Added publishing mkdocs to gh-pages. [#594](https://github.com/OpenMDAO/dymos/pull/594)
* Added github actions workflow for CI. [#589](https://github.com/OpenMDAO/dymos/pull/589)
* Readme updated to point to JOSS paper. [#578](https://github.com/OpenMDAO/dymos/pull/578)
* Updated JOSS bibliography. [#573](https://github.com/OpenMDAO/dymos/pull/573) [#574](https://github.com/OpenMDAO/dymos/pull/574) [#576](https://github.com/OpenMDAO/dymos/pull/576) [#577](https://github.com/OpenMDAO/dymos/pull/577)
* Fixed some references in JOSS paper. [#572](https://github.com/OpenMDAO/dymos/pull/572)
* Minor grammar and consistency edits for JOSS paper. [#571](https://github.com/OpenMDAO/dymos/pull/571)

********************************
# Release Notes for Dymos 1.0.0

March 25, 2021

This is version 1.0.0 of Dymos.
This release primarily removes some deprecated experimental features, along with implementing a few bug fixes.

## Backwards Incompatible API Changes & Deprecations

* The RungeKutta Transcription is removed. [#550](https://github.com/OpenMDAO/dymos/pull/550)
* Two-character location specifiers ('++', '--', '+-', '-+') are removed in favor of 'initial' and 'final' [#556](https://github.com/OpenMDAO/dymos/pull/556)
* User must now specify `solve_segments='forward'` or `solve_segments='backward'` when using solve_segments capability (True is no longer valid). [#557](https://github.com/OpenMDAO/dymos/pull/557)
* `add_input_parameter` and `add_design_parameter` dropped in favor of `add_parameter`. [#558](https://github.com/OpenMDAO/dymos/pull/558) [#561](https://github.com/OpenMDAO/dymos/pull/561)
* Removed the deprecated command line interface due to it being a somewhat hackish abuse of OpenMDAO's hooks. [#563](https://github.com/OpenMDAO/dymos/pull/563)
* Removed the deprecated 'disc' node subset in favor of the more correct 'state_disc' [#565](https://github.com/OpenMDAO/dymos/pull/565)
* Removed the deprecated 'custom_targets' option for parameters.' [#565](https://github.com/OpenMDAO/dymos/pull/565)

## Enhancements

* User will now be warned when multiple timeseries outputs were attempted with the same name. [#567](https://github.com/OpenMDAO/dymos/pull/567)

## Bug Fixes

* The `ode_class` option of phase is now marked non-recordable.  [#555](https://github.com/OpenMDAO/dymos/pull/555)
* Fixed doc building due to a change in OpenMDAO. [#564](https://github.com/OpenMDAO/dymos/pull/564)

## Miscellaneous

* Modified examples to use the preferred `fix_initial=True` for time and states instead of "pinched bounds", e.g. `initial_bounds=(0, 0)`. [#560](https://github.com/OpenMDAO/dymos/pull/560)


********************************
# Release Notes for Dymos 0.18.1

February 18, 2021

This release of Dymos adds several examples demonstrating various capabilities of the code.
Per user request, the ODE of a system can now be provided via a callable function that
takes `num_nodes` and any other potential initialization keywod arguments.
This allows OpenMDAO's ExecComp to be used as an ODE if wrapped in a lambda, for instance.

Some bugs were fixed as part of introducing these examples.  For instance, default values for
states and times were being ignored - this is now fixed.  In addition, Trajectory parameters
are now saved to the simulation database file.

This is expected to be the final release of Dymos before v1.0.0, when several existing
deprecated features will be removed from the code.

## Backwards Incompatible API Changes & Deprecations

None

## Enhancements

* Added the racecar example from @pwmdebuck, demonstrating cycling constraints and integration about arclength instead of time. [#535](https://github.com/OpenMDAO/dymos/pull/535)
* Simplified the SSTO example to use a single ODE component with complex-step. [#534](https://github.com/OpenMDAO/dymos/pull/534)
* Added a balanced field length example. [#533](https://github.com/OpenMDAO/dymos/pull/533)
* Added the ability to use a callable that returns a System as an ODE. [#528](https://github.com/OpenMDAO/dymos/pull/528)

## Bug Fixes

* Removed duplication of inputs to timeseries when multiple outputs may use the same source data. [#543](https://github.com/OpenMDAO/dymos/pull/543)
* Fixed a bug where timeseries outputs of non-dynamic ODE outputs would cause an exception. [#521](https://github.com/OpenMDAO/dymos/pull/521)

## Miscellaneous

* Reworked the cannonball example to make it more simple [#545](https://github.com/OpenMDAO/dymos/pull/545)
* Placed some more tests under the use_tempdirs decorator to clean up output. [#530](https://github.com/OpenMDAO/dymos/pull/530)
* Added a test that checks all docstrings vs the NumpyDoc standard. [#526](https://github.com/OpenMDAO/dymos/pull/526)
* Clean up implementation of wildcard units when adding multiple timeseries at once. [#523](https://github.com/OpenMDAO/dymos/pull/523)

********************************
# Release Notes for Dymos 0.18.0

January 21, 2021

This release of Dymos brings a few improvements and bug fixes.

We've implemented introspection for boundary and path constraints such that units and shapes of the constrained quantities can be determined during the setup process and no longer are required to be specified.
The use of solve_segments to have a solver converge the dynamics and thus provide a shooting method at the optimizer has been improved.  Option solve_segments can now take on values 'forward' or 'backward' to provide forward or backward shooting.
Linkage constraints will now result in an error if the quantity on each side of the linkage is governed by the optimizer (is a design variable) but is fixed in value (not allowed to be varied).

## Backwards Incompatible API Changes & Deprecations

* The command line interface is deprecated due to inconsistent behavior. [#505](https://github.com/OpenMDAO/dymos/pull/505)

## Enhancements

* Changed CoerceDesvar to handle vector values for ref/ref0/scaler/adder. Carried over changes to for defect ref/defect scaler. [#464](https://github.com/OpenMDAO/dymos/pull/464)
* Changed solve_segments to allow it to be used when neither fix_initial nor fix_final are True. [#490](https://github.com/OpenMDAO/dymos/pull/490)
* Added detection of unsatisfiable linkage constraints (where quantities on both sides cannot be changed by the optimizer). [#502](https://github.com/OpenMDAO/dymos/pull/502)
* Added introspection for boundary and path constraints. [#506](https://github.com/OpenMDAO/dymos/pull/506)

## Bug Fixes

* Fixed a bug with using grid refinement if one or more of the states are vector-valued. [#492](https://github.com/OpenMDAO/dymos/pull/492)
* Fixed issue where show_plots argument to run_problem is permanently changing matplotlib backend. [#497](https://github.com/OpenMDAO/dymos/pull/497)
* Fixed a bug where trajectory.simulate() would crash if trajectory were not named the common 'traj' name. [#500](https://github.com/OpenMDAO/dymos/pull/500)
* Fixed a bug that prevented restart from working with multiphase trajectories, and now include static parameters in the solution file. [#510](https://github.com/OpenMDAO/dymos/pull/510)
* Fixed a bug associated with polynomial control rate sources, and added coverage of more cases. [#513](https://github.com/OpenMDAO/dymos/pull/513)
* Fixed exception when simulating a trajectory with a parameter for which None was declared as the target for one of its phases. [#515](https://github.com/OpenMDAO/dymos/pull/515)

## Miscellaneous

None

****************************************************************
# Release Notes for Dymos 0.17.0

December 14, 2020

This release of Dymos adds a few important features:
- Dymos can now automatically generate plots of all timeseries outputs
- State rates may now be tagged in the ODE

There are also several bug fixes.

## Backwards Incompatible API Changes & Deprecations

* None

## Enhancements

* Dymos can now automatically generate plots of all timeseries outputs. [#469](https://github.com/OpenMDAO/dymos/pull/469)
* States can now be discovered automatically by tagging their rate variables in the ODE. [#477](https://github.com/OpenMDAO/dymos/pull/477) [#481](https://github.com/OpenMDAO/dymos/pull/481) [#484](https://github.com/OpenMDAO/dymos/pull/484)

## Bug Fixes

* Fixed issue where polynomial controls were not able to be linked across phases. [#462](https://github.com/OpenMDAO/dymos/pull/462)
* Fixed a bug that was preventing parameters from being state rate sources. [#466](https://github.com/OpenMDAO/dymos/pull/466)
* Users may now specify parameter shapes as integers or iterables other than tuples. [#467](https://github.com/OpenMDAO/dymos/pull/467)

## Miscellaneous

* Fix for apache license string in setup.py classifiers since it was not recognized by PyPI. [#458](https://github.com/OpenMDAO/dymos/pull/458)
* Added a long description for PyPI. [#460](https://github.com/OpenMDAO/dymos/pull/460)
* Some documentation cleanup for the JOSS review  [#474](https://github.com/OpenMDAO/dymos/pull/474) [#488](https://github.com/OpenMDAO/dymos/pull/488)
* Dropped the dependency on the parameterized package with the intent to utilize subTest in the future.  [#479](https://github.com/OpenMDAO/dymos/pull/479)


********************************
# Release Notes for Dymos 0.16.1

November 16, 2020

This release of Dymos fixes an issue that now allows portions of an array
output to be connected to a parameter.

This version works with OpenMDAO 3.3.0 but version 3.4.1 offers some
improved handling of parameters.

## Backwards Incompatible API Changes & Deprecations

* Parameter shapes in 0.16.1 were stored as (1,) + the shape of the parameter.  Now they are shaped as expected. [#444](https://github.com/OpenMDAO/dymos/pull/444)

## Enhancements

* State shapes and units, if not explicitly given, are now pulled from targets (if present and uniquely defined), or from the rate source variable. [#449](https://github.com/OpenMDAO/dymos/pull/449)

## Bug Fixes

* Fixes a bug where user-defined shapes of states were colliding with those found during introspection, and other state introspection updates. [#449](https://github.com/OpenMDAO/dymos/pull/449)

* Fixes a bug that prevented the use of numpy arrays as a boundary constraints [#450](https://github.com/OpenMDAO/dymos/pull/450)

## Miscellaneous

* Added test to verify functionality of OpenMDAO 3.4.1 that allows the final value of a control (or a partial portion of any output) to be connected to a parameter. [#445](https://github.com/OpenMDAO/dymos/pull/445)

********************************
# Release Notes for Dymos 0.16.0

October 23, 2020

This is a long-overdue release of Dymos that features several useful new features.
By taking advantage of recent changes to the OpenMDAO setup process, we can use introspection to determine a lot of things automatically now.
For instance, the user will no longer be required to specify the shape or default units of timeseries outputs or constraints (although they can override units if they choose).

Also on the topic of timeseries, dymos now supports glob patterns when adding ODE outputs to the timeseries.
Adding all 100 outputs of a large ODE to the timeseries output can now be accomplished in a single line:

```
phase.add_timeseries_output('*')
```

Introspection is also used for states, time, and controls.
States will get their default units and shape from the associated state rate variable now.
Controls, Polynomial Controls, and Parameters will now get their default units and shapes from their targets.

With the advent of automatic IndepVarComps in OpenMDAO, there no longer needs to be a distinction between InputParameters and DesignParameters.
From here on out, users just add `parameters` to a Trajectory or Phase.
A user can choose to make them design variables (in which case they perform as DesignParameters), or not.
Connecting an external variable to the parameter makes it function as a InputParameter.

From here on, we plan on making releases of Dymos roughly once per month.

## Backwards Incompatible API Changes & Deprecations:

- `add_input_parameter` and `add_design_parameter` are **deprecated** and replaced with `add_parameter` [#365](https://github.com/OpenMDAO/dymos/pull/365)
- The use of two-character location strings ('--' and '++') are replaced by 'initial' and 'final' when linking phases. [#427](https://github.com/OpenMDAO/dymos/pull/427)
- The endpoint conditions component, used to pull out initial or final values of a variable in a phase, have been removed.  Values should instead by pulled from the timeseries. [#427](https://github.com/OpenMDAO/dymos/pull/427)
- The results of simulate are now automatically recorded to a database recorder 'dymos_simulation.db' and stored in a case named 'final'.  Use `CaseReader('dymos_simulation.db').get_case('final')` to access the case. instead of the previous behavior `CaseReader('dymos_simulation.db').get_case(-1)` [#437](https://github.com/OpenMDAO/dymos/pull/437)
- The RungeKutta transcription is **deprecated**. It is functionally equivalent to the GaussLobatto transcription with an order of 3, but since it's signficantly different under the hood it was being a drag on development.

One note on another backwards incompatibility.
Since we're using some newer OpenMDAO capabilities for the way in which parameters are handled, it is currently not possible to connect an upstream output to a parameter in dymos _when src_indices are used_.
We plan on fixing this issue with an upcoming OpenMDAO release.
If you rely on this capability, we recommend waiting for that update before using this version of Dymos, or using an intermediate "pass-through" component to pull the correct index from the original output and then passing _that_ as the value of a dymos parameter.

## Enhancements:

* Dymos components are no longer listed in check_partials output by default.  Use `dymos.options['include_check_partials'] = True` to override this for debugging.  [#438](https://github.com/OpenMDAO/dymos/pull/438)
* Any outputs can now be linked across phases.  Added a more general and powerful `add_linkage_constraint` method, but the simpler `link_phases` method remains in place for simple continuity. [#422](https://github.com/OpenMDAO/dymos/pull/422)
* Time, controls, and polynomial control units are now determined from targets. [#412](https://github.com/OpenMDAO/dymos/pull/412)
* States now pull their default units and shapes from the rate source. [#398](https://github.com/OpenMDAO/dymos/pull/398)
* Phase method `add_timeseries_output` now supports wildcards to allow multiple timeseries outputs to be added at once. [#387](https://github.com/OpenMDAO/dymos/pull/387), #387
* Timeseries outputs automatically detect shape and default units. [#380](https://github.com/OpenMDAO/dymos/pull/380)
* Added ph-adaptive refinement method that is capable of shrinking the grid. [#379](https://github.com/OpenMDAO/dymos/pull/379)
* Control targets will automatically be set to a top-level input of the ODE, if present. [#373](https://github.com/OpenMDAO/dymos/pull/373)
* State targets will automatically be set to a top-level input of the ODE, if present. [#356](https://github.com/OpenMDAO/dymos/pull/356)
* Water-powered rocket MDO example added. [#343](https://github.com/OpenMDAO/dymos/pull/343)
* State rates are included in timeseries outputs by default. [#329](https://github.com/OpenMDAO/dymos/pull/329)
* By default, run_problem will use a problem recorder to record only the 'final' solution after an optimization. [#328](https://github.com/OpenMDAO/dymos/pull/328)
* Move most setup functionality to configure to allow more introspection changes. [#327](https://github.com/OpenMDAO/dymos/pull/327)
* Include examples using IPOPT via pyoptsparse. [#311](https://github.com/OpenMDAO/dymos/pull/311)
* Parameters may now be selectively omitted from the timeseries outputs, but are included by default. [#305](https://github.com/OpenMDAO/dymos/pull/305)
* Added a test of the distributed ODE capability. [#304](https://github.com/OpenMDAO/dymos/pull/304)
* Replace `load_case` with `reinterpolate_solution`. [#296](https://github.com/OpenMDAO/dymos/pull/296)

## Bug Fixes:

* Fixed a bug that was causing matrix-shaped states to have connection errors. [#441](https://github.com/OpenMDAO/dymos/pull/441)
* Fixed a bug where simulate failed when both standard and polynomial controls were present. [#434](https://github.com/OpenMDAO/dymos/pull/434)
* Fixed an issue where parameters are promoted twice in some cases. This will cause errors in an upcoming versions of OpenMDAO. [#429](https://github.com/OpenMDAO/dymos/pull/429)
* Fixed a bug that was causing matrix-shaped states to have connection errors. [#441](https://github.com/OpenMDAO/dymos/pull/441)
* AnalysisError is now raised if scipy.integrate.solve_ivp fails during simulation. [#400](https://github.com/OpenMDAO/dymos/pull/400)
* Fixed a bug involving GaussLobatto transcriptions `get_rate_source` method. [#334](https://github.com/OpenMDAO/dymos/pull/334)
* Fix for simulate encountering errors when time options `input_initial` or `input_duration` were True. [#317](https://github.com/OpenMDAO/dymos/pull/317)
* Fixed a bug in which shaped input parameters were breaking simulate. [#301](https://github.com/OpenMDAO/dymos/pull/301)
* Sort linkage orders to make convergence more repeatable. [#298](https://github.com/OpenMDAO/dymos/pull/298)

## Miscellaneous:

* Switch to use preferred `assert_near_equal` method instead of `assert_rel_error` from OpenMDAO. [#320](https://github.com/OpenMDAO/dymos/pull/320)
* Documentation is now handled via mkdocs instead of sphinx. [#337](https://github.com/OpenMDAO/dymos/pull/337), [#332](https://github.com/OpenMDAO/dymos/pull/332)

## 0.15.0
2020-02-12

* [__Bug Fix__] Phase Linkage units now checks all optimal control variables from the first phase in the linkage for units. Previously units for things that were neither states, time, nor controls were not being assigned.
* [__Enhancement__] Removed Python 2 support.
* [__Enhancement__] Removed deprecated ODE Decorators.
* [__Enhancement__] Added ability to subclass Phase with an assigned ODE and default state variable options.
* [__Docs__] Added docs on _subclassing phases_.
* [__Enhancement__] Automated grid refinement is now available via the `dymos.run_problem` function.
* [__Enhancement__] Grid data is now available upon instantiation of the Transcription instead of being deferred to setup time.
* [__Bug Fix__] The user now gets a meaningful message if Phase.interpolate is called too soon.
* [__Bug Fix__] State rates are now correctly passed through the interleave component that provides timeseries outputs for Gauss-Lobatto transcription.
* [__Enhancement__] Added hypersensitive example problem.
* [__Docs__] Documentation added for grid refinement.
* [__Enhancement__] Deprecated the use of ODE decorators
* [__Enhancement__] Added shuttle reentry example problem

## 0.13.0
2019-07-18

* [__Enhancement__] Phase methods like `set_state_options` and `add_control` no longer use **kwargs in order to make them more IDE friendly.
* [__Enhancement__] Additional timeseries can for outputs can be added to a phase, with interpolation onto a new set of grid points. This enables the concept of tandem phases, where two different ODE's operating over the same time interval can be integrated on different grids. Fast state variables can be integrated on a dense grid while slower state variables are integrated on a more sparse grid, for performance. For an example see the tandem phase documentation in the feature docs.
* [__Enhancement__] ODE options can be specified at the phase level rather than in the ODE. This feature is experimental, but it allows one way of programmatically defining an ODE.
* [__Enhancement__] Changed the use of OpenMDAO to do `import openmdao.api as om` for consistency with the OpenMDAO documentation.
* [__Enhancement__] Dymos is now imported in the examples as import `dymos as dm`
### Issue Type

- [x] Bug
- [ ] Enhancement
- [ ] Docs
- [ ] Miscellaneous

### Description

For bugs, explain what is happening and how it differs from the expected behavior. For enhancements, explain the desired functionality.

### Example

For bugs, if at all possible, provide a short snipped which reproduces the issue,
or link to a file in another repository where the issue is demonstrated.

For enhancements this is not as critical, but if possible provide an example
of the proposed usage of the enhancement.  If a reference implementation has
been implemented and you desire to have it merged into Dymos, site that reference
implementation here before issuing a pull request. The pull request should then  
cite this issue number in the "Related issues" section.

### Environment

Operating System: <i.e. OS X 10.14.6, Windows 10, Ubuntu 16.04>
Python environment: <i.e. Anaconda Python 3.7.1>
Packages: <versions of Dymos, OpenMDAO, and other pertinent Python packages. You can just paste the output of `pip freeze` here.>
### Summary

Summary of PR.

### Related Issues

- Resolves #

### Backwards incompatibilities

None

### New Dependencies

None
---
title: 'dymos: A Python package for optimal control of multidisciplinary systems'
tags:
  - Python
  - OpenMDAO
  - optimal control
  - trajectory optimization
  - multidisciplinary optimization
  - NASA
authors:
  - name: Robert Falck
    orcid: 0000-0001-9864-4928
    affiliation: 1
  - name: Justin S. Gray
    orcid: 0000-0002-7506-7360
    affiliation: 1
  - name: Kaushik Ponnapalli
    affiliation: 2
  - name: Ted Wright
    affiliation: 1
affiliations:
  - name: NASA Glenn Research Center
    index: 1
  - name: HX5 LLC
    index: 2
date: 23 March 2021
bibliography: paper.bib
---

# Summary

Dymos is a library for optimizing control schedules for dynamic systems --- sometimes referred to as  optimal control or trajectory optimization.
There are a number of other optimal control libraries that tackle similar kinds of problems, such as OTIS4 [@Paris2006], GPOPS-II [@Patterson2014GPOPSII],and CASADI [@Andersson2018].
These tools all rely on gradient-based optimization to solve optimal control problems, though their methods of computing the gradients vary. 
Dymos is built on top of the OpenMDAO framework [@Gray2019a] and supports its modular derivative system which allows users to mix-and-match from finite-differencing, complex-step, hand-differentiated, and algorithmic differentiation.
This flexibility allows Dymos to efficiently solve optimal control problems constructed with both ordinary differential equations (ODE) and differential-algebraic equations (DAE). 

Dymos can also help solve more general optimization problems where dynamics are only one part in a larger system-level model with additional --- potentially computationally expensive --- calculations that come before and after the dynamic calculations.
These broader problems are commonly referred to as co-design, controls-co-design, and multidisciplinary design optimization.
Dymos provides specific APIs and features that make it possible to integrate traditional optimal-control models into a co-design context, while still supporting analytic derivatives that are necessary for computational efficiency in these complex use cases.
An example of a co-design problem that was solved with Dymos is the coupled trajectory-thermal design of an electric vertical takeoff and landing aircraft where the thermal management and propulsion systems were designed simultaneously with the flight trajectories to ensure no components overheated [@Hariton2020a].

# Difference between optimal-control and co-design

Optimal-control and co-design problems deal with dynamic systems.
The evolution of the states over time is governed by an ordinary differential equation (ODE) or differential-algebraic equation (DAE):
\begin{align*}
  \dot{\bar{x}} = f_{ode}(\bar{x},t,\bar{u},\bar{d})
\end{align*}
Here, $\bar{x}$ is a vector of time-varying state variables whose behavior is affected by time ($t$), a vector of dynamic controls ($\bar{u}$), and a vector of static design parameters ($\bar{d}$).


To optimize a dynamic system we also need to account for the objective function ($J$):
\begin{align*}
  \mathrm{J} = f_{obj}(\bar{x},t,\bar{u},\bar{d})
\end{align*}
In addition, there are constraints that typically need to be enforced: 
\begin{align*}
  \mathrm{Time:}& \qquad {t}_{lb} \leq t \leq {t}_{ub} \\
  \mathrm{State \, Variables:}& \qquad \bar{x}_{lb} \leq \bar{x} \leq \bar{x}_{ub} \\
  \mathrm{Dynamic \, Controls:}& \qquad \bar{u}_{lb} \leq \bar{u} \leq \bar{u}_{ub} \\
  \mathrm{Design \, Parameters:}& \qquad \bar{d}_{lb} \leq \bar{d} \leq \bar{d}_{ub} \\
  \mathrm{Initial \, Boundary \, Constraints:}& \qquad \bar{g}_{0,lb} \leq g_{0}(\bar{x}_0,t_0,\bar{u}_0, \bar{d}) \leq \bar{g}_{0,ub} \\
  \mathrm{Final \, Boundary \, Constraints:}& \qquad \bar{g}_{f,lb} \leq g_{f}(\bar{x}_f,t_f,\bar{u}_f, \bar{d}) \leq \bar{g}_{f,ub} \\
  \mathrm{Path \, Constraints:}& \qquad \bar{p}_{lb} \leq p(\bar{x},t,\bar{u},\bar{d}) \leq \bar{p}_{ub} \\
\end{align*}


In the mathematical sense what distinguishes optimal control from co-design is the particulars of which design variables and constraints are actually considered by the optimization.
Pure optimal-control problems deal with a system of fixed design and seek to maximize performance by adjusting dynamic quantities ($t, \bar{x}, \bar{u}$) such as position, speed, fuel-burned, and battery state-of-charge.
Co-design problems simultaneously vary the static design parameters of a system ($\bar{d}$) and its dynamic behavior ($t, \bar{x}, \bar{u}$) to reach maximum performance. 

In practice, the mathematical distinction is too rigid and a more practical distinction is made based on where the static and dynamic calculations are implemented and how complex each of them is. 
For very simple physical design parameters (e.g. the radius of a cannon ball, spring constants, linkage lengths, etc) it is common to integrate the design calculations directly into the ODE.
Even though the calculations are static in nature, they can easily be coded as part of the ODE and still fit well into the optimal-control paradigm.
The optimization structure thus looks like this: 

![Model structure for a traditional optimal control problem](flow_charts/opt_control.png){width=45%}

However, not all problems can be handled with such a compact implementation.
For example if the physical design problem included shaping of an airfoil using expensive numerical solutions of partial differential equations (PDE) to predict drag, then one would not want to embed that PDE solver into the dynamic model.
Instead the user could set up a coupled model with the PDE solver going first, and passing a table of data to be interpolated to the dynamic model.
This effectively splits calculations up into static and dynamic components.
This implementation structure is called co-design. 

Traditionally, this co-design implementation would be done via sequential optimization with a manual outer design iteration between the static and dynamic models, potentially with different teams of people working on each one.
One team would come up with a physical design using their own internal optimization setup. 
A second team takes the design and generates optimal-control profiles for it. 
Of course, the iterations do not need to be manual.
It is also possible to set up an iterative loop around static and dynamic models to converge the problem numerically. 
A sequential co-design implementation looks like this:

![Model structure for a sequential co-design problem](flow_charts/sequential_co_design.png){width=100%}

Dymos can support sequential co-design, but its unique value is that it also enables a more tightly-coupled
co-design process with a single top level optimizer handling both parts of the problem simultaneously. 

![Model structure for a coupled co-design problem](flow_charts/coupled_co_design.png){width=75%}

Compared to sequential co-design, coupled co-design offers the potential to find better designs with much lower computational cost. 
However, it is also more challenging to implement because the top-level optimizer requires derivatives to be propagated between the static and dynamic parts of the model. 
Dymos overcomes this difficulty by providing APIs to exploit OpenMDAO's analytic derivative functionality at the model level. 
Data can be passed from the static model to the dynamic model and vice versa, allowing the construction of the coupled model for optimization. 

# ODE versus DAE

Optimal-control software typically requires that the dynamics of the system be defined as a set of ordinary differential equations (ODE) that use explicit functions to compute the rates of the state variables to be time-integrated.
Sometimes the dynamics are instead posed as a set of differential-algebraic equations (DAE), where some residual equations need to be satisfied implicitly in order to solve for the state rates. 
From the perspective of an optimal-control or co-design problem both ODE and DAE formulations provide state rates that need to be integrated over time. 
The difference is that ODEs are explicit functions which are relatively easy to differentiate, but DAEs are implicit functions which are much more difficult to differentiate. 
Since the derivatives are needed to perform optimization, DAEs are more challenging to optimize. 

One relatively common use case for DAEs is differential inclusions, in which the state trajectory is posed as a dynamic control and the traditional control schedule needed to achieve that trajectory is found using a nonlinear solver within the dynamic model [@Seywald1994].
For some problems this method provides a more natural and numerically-beneficial design space for the optimizer to traverse,
but the nonlinear solver poses numerical challenges for computing derivatives for the optimizer.
A simple approach to this is to just use finite-differences across the nonlinear solver, but this has been shown to be expensive and numerically unstable [@gray2014derivatives].
Another option, taken by some optimal control libraries, is to apply monolithic algorithmic differentiation [@griewank2003mathematical] across the nonlinear solver.
While it does provide accurate derivatives, the monolithic approach is expensive and uses a lot of memory [@mader2008adjoint; @kenway2019effective].
The most efficient approach is to use a pair of analytic derivative approaches called the direct and adjoint methods, which were generalized in a single unified derivative equation (UDE) by Hwang and Martins [@hwang2018b].

Dymos adopts the UDE approach, which uses a linear solver to compute total derivatives needed by the optimizer using only partial derivatives of the residual equations in the DAE.
This approach offers two key advantages. 
First, partial derivatives of the DAE residual equations are much less computationally challenging to compute. 
Second, by using the OpenMDAO underpinnings of Dymos, users can construct their DAE in a modular fashion and combine various methods of computing the partial derivatives via finite-difference, complex-step [@Martins2003CS], algorithmic differentiation, or hand differentiation as needed. 

## The Dymos perspective on optimal control

Dymos breaks the trajectory into portions of time called _phases_.
Breaking the trajectory into phases provides several capabilities.
Intermediate constraints along a trajectory can be enforced by applying a boundary constraint to a phase that begins or ends at the time of interest.
For instance, the optimal trajectory of a launch vehicle may be required to ascend vertically to clear a launch tower before it pitches over on its way to orbit.
Path constraints can be applied within each phase to bound some performance parameter within that phase.
For example, reentry vehicles may need to adjust their trajectory to limit aerodynamic heating.

Each phase in a trajectory can use its own separate ODE.
For instance, an aircraft with vertical takeoff and landing capability may use different ODEs for vertical flight and horizontal flight.
ODEs are implemented as standard OpenMDAO models which are passed to phases at instantiation time with some additional annotations to identify the states, state-rates, and control inputs.

Every phase uses its own specific time discretization tailored to the dynamics in that portion of the trajectory.
If one part of a trajectory has fast dynamics and another has slow dynamics, the time evolution can be broken into two phases with separate time discretizations.

In the optimal-control community there are a number of different techniques for discretizing the continuous optimal control problem into a form that can be solved by a nonlinear optimization algorithm; each one is called a transcription.
Dymos supports two different collocation transcriptions: high-order Gauss-Lobatto [@Herman1996] and Radau [@Garg2009].
Both of these represent state and control trajectories as piece-wise polynomials of at least 3rd order and are formulated in a way that makes it possible to efficiently compute the needed quantities to perform integration in a numerically rigorous fashion.

Dymos also allows the user to choose whether the optimization problem is solved using an explicit or implicit approach.
Some caution with terminology must be taken here because the term "implicit" is often used to describe time integration schemes (e.g. backwards Euler),
but that is not what is meant in an optimal-control context.
Here, explicit propagation is one where the full state trajectory is computed starting from the initial value and propagating forward or from the final value and propagating backward.
From the optimizer's perspective it will set values for the initial or final state ($\bar{x}$), the design parameters ($\bar{d}$), and the controls ($\bar{u}$) and can expect to be given a physically valid time evolution of the states as the output.
Wrapping an optimizer around an explicit propagation gives what is traditionally called a "shooting method" in the 
optimal-control world.
In contrast, implicit propagation used within an optimization does not provide valid trajectories on its own.
Instead, implicit methods add a discretized time-evolution of the state vector ($\bar{x}$) as an additional design variable to the optimizer and add an associated set of defect constraints that must be driven to zero to enforce physics at some set of discrete points in time where the ODE is evaluated.
The net effect is that the full state trajectory is only known once the optimization is fully converged.
In the context of the multidisciplinary design optimization field, explicit phases are similar to the multidisciplinary design feasible (MDF) optimization architecture and implicit phases are similar to the simultaneous analysis and design (SAND) optimization architecture [@Martins2013].

Both implicit and explicit phases are useful in different circumstances. 
Explicit propagation can seem to many like a more natural way to formulate the problem because it matches the way one would use time integration without optimization.
However, when used with optimization explicit propagation is more computationally expensive,
sensitive to singularities in the ODE, 
and potentially unable to converge to a valid solution. 
Implicit propagation tends to be less intuitive computationally, since it does not provide valid state histories without a converged optimization.
The advantages of implicit propagation are that it tends to be faster, more numerically stable, and more scalable --- though also highly sensitive to initial conditions and optimization scaling.

Dymos supports both explicit and implicit propagation for both its transcriptions,
and even allows mixtures of implicitly and explicitly propagated states within a phase.
This flexibility is valuable because it allows users to tailor their optimization to suit their needs. 
Switching transcriptions and changing from implicit to explicit requires very minor code changes --- typically a single line in the run-script.
Examples of how to swap between them are given in the code sample below. 

# Choice of optimization algorithm 

Dymos is not distributed with an optimizer, but relies on the optimizers that are available in the OpenMDAO installation.
OpenMDAO ships with an interface to the optimizers in SciPy [@2020SciPy-NMeth], 
and an additional wrapper for the pyoptsparse library [@Wu_pyoptsparse_2020] which has more powerful optimizer options such as SNOPT [@GilMS05] and IPOPT [@wachter2006].
OpenMDAO also allows users to integrate their own optimizer of choice, which Dymos can then seamlessly use with without any additional modifications.
For simple problems, Scipy's SLSQP optimizer generally works fine.
On more challenging optimal-control problems higher-quality optimizers are important for getting good performance.

Though one could technically choose any optimization algorithm, Dymos is designed to work primarily with gradient-based algorithms.
In general, optimal-control and co-design problems will have both a very large number of design variables and a very large number of constraints. 
Both of these issues make gradient-based methods the strongly-preferred choice. 
Gradient-free methods could potentially be used in certain narrow circumstances with problems built using purely explicit phases and set up intentionally to have a small set of design variables and constraints. 


## Statement of Need

When dealing with the design of complex systems that include transient behavior, co-design becomes critical [@garciasans2019].
Broadly there are two approaches: sequential co-design or coupled co-design [@Fathy2001;  @Peters2009].
The best choice depends on the degree of interaction, or coupling, between various sub-systems. 
If the coupling is strong a coupled co-design approach is necessary to achieve the best performance.

Though there are a number of effective optimal-control libraries, they tend to assume that they are on top of the modeling stack.
They frame every optimization problem as if it were a pure optimal-control problem, and hence are best suited to be used in a sequential co-design style. 
This poses large challenges when expanding to tightly-coupled problems, where the interactions between the static and dynamic systems are very strong. 

Dymos provides a set of unique capabilities that make coupled co-design possible via efficient gradient-based optimization methods.
It provides differentiated time-integration schemes that can generate transient models from user provided ODEs, 
along with APIs that enable users to couple these transient models with other models to form the co-design system while carrying the differentiation through that coupling.
It also supports efficient differentiation of DAEs that include implicit relationships, which allows for a much broader set of possible ways to pose transient models. 
These two features combined make Dymos capable of handling coupled co-design problems in a manner that is more efficient than a pure optimal-control approach. 


## Selected applications of Dymos

Dymos has been used to demonstrate the coupling of flight dynamics and subsystem thermal constraints in electrical aircraft applications [@Falck2017a; @Hariton2020a].
NASA's X-57 "Maxwell" is using Dymos for mission planning to maximize data collection while abiding the limits of battery storage capacity and subsystem temperatures [@Schnulo2018a; @Schnulo2019a].
Other authors have used Dymos to perform studies of aircraft acoustics [@Ingraham2020a] and the design of supersonic aircraft with thermal fuel management systems [@Jasa2018a].

## Optimal-control example: Brachistochrone

As a simple use case of Dymos, consider the classic brachistochrone optimal-control problem.
There is a bead sliding along a frictionless wire strung between two points of different heights, 
and we seek the shape of the wire such that the bead travels from start to finish in the shortest time. 
The time-varying control is the angle of the wire at each point in time and there are no design parameters, 
which makes this a pure optimal-control problem. 

\small
```python
import numpy as np
import openmdao.api as om
import dymos as dm
import matplotlib.pyplot as plt

# First define a system which computes the equations of motion
class BrachistochroneEOM(om.ExplicitComponent):
    def initialize(self):
        self.options.declare('num_nodes', types=int)

    def setup(self):
        nn = self.options['num_nodes']

        # Inputs
        self.add_input('v', val=np.zeros(nn), units='m/s', desc='velocity')
        self.add_input('theta', val=np.zeros(nn), units='rad',
                       desc='angle of wire')
        self.add_output('xdot', val=np.zeros(nn), units='m/s',
                        desc='x rate of change')
        self.add_output('ydot', val=np.zeros(nn), units='m/s',
                        desc='y rate of change')
        self.add_output('vdot', val=np.zeros(nn), units='m/s**2',
                        desc='v rate of change')

        # Ask OpenMDAO to compute the partial derivatives using complex-step
        # with a partial coloring algorithm for improved performance
        self.declare_partials(of='*', wrt='*', method='cs')
        self.declare_coloring(wrt='*', method='cs', show_summary=True)

    def compute(self, inputs, outputs):
        v, theta = inputs.values()
        outputs['vdot'] = 9.80665 * np.cos(theta)
        outputs['xdot'] = v * np.sin(theta)
        outputs['ydot'] = -v * np.cos(theta)

p = om.Problem()

# Define a Trajectory object
traj = p.model.add_subsystem('traj', dm.Trajectory())

# Define a Dymos Phase object with GaussLobatto Transcription
tx = dm.GaussLobatto(num_segments=10, order=3)
phase = dm.Phase(ode_class=BrachistochroneEOM, transcription=tx)

traj.add_phase(name='phase0', phase=phase)

# Set the time options
phase.set_time_options(fix_initial=True,
                       duration_bounds=(0.5, 10.0))
# Set the state options
phase.add_state('x', rate_source='xdot',
                fix_initial=True, fix_final=True)
phase.add_state('y', rate_source='ydot',
                fix_initial=True, fix_final=True)
phase.add_state('v', rate_source='vdot',
                fix_initial=True, fix_final=False)
# Define theta as a control.
phase.add_control(name='theta', units='rad',
                  lower=0, upper=np.pi)
# Minimize final time.
phase.add_objective('time', loc='final')

# Set the driver.
p.driver = om.ScipyOptimizeDriver()

# Allow OpenMDAO to automatically determine total
# derivative sparsity pattern.
# This works in conjunction with partial derivative
# coloring to give a large speedup
p.driver.declare_coloring()

# Setup the problem
p.setup()

# Now that the OpenMDAO problem is setup, we can guess the
# values of time, states, and controls.
p.set_val('traj.phase0.t_duration', 2.0)

# States and controls here use a linearly interpolated
# initial guess along the trajectory.
p.set_val('traj.phase0.states:x',
          phase.interpolate(ys=[0, 10], nodes='state_input'),
          units='m')
p.set_val('traj.phase0.states:y',
          phase.interpolate(ys=[10, 5], nodes='state_input'),
          units='m')
p.set_val('traj.phase0.states:v',
          phase.interpolate(ys=[0, 5], nodes='state_input'),
          units='m/s')
# constant initial guess for control
p.set_val('traj.phase0.controls:theta', 90, units='deg')

# Run the driver to solve the problem and generate default plots of
# state and control values vs time
dm.run_problem(p, make_plots=True, simulate=True)

# Additional custom plot of y vs x to show the actual wire shape
fig, ax = plt.subplots(figsize=(6.4, 3.2))
x = p.get_val('traj.phase0.timeseries.states:x', units='m')
y = p.get_val('traj.phase0.timeseries.states:y', units='m')
ax.plot(x,y, marker='o')
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')
fig.savefig('brachistochone_yx.png', bbox_inches='tight')
```
\normalsize

The built-in plotting utility in Dymos will plot all relevant quantities vs time:

![Brachistochrone Solution: y state trajectory](brachistochrone_states_y.png){width=50%}
![Brachistochrone Solution: x state trajectory](brachistochrone_states_x.png){width=50%}


The more traditional way to view the brachistochrone solution is to view the actual shape of the wire (i.e. y vs x):

![Brachistochrone Solution: y as a function of x](brachistochone_yx.png){width=50%}


## Coupled co-design example: Designing a cannonball

This co-design example seeks to find the best size cannonball to maximize range, considering aerodynamic drag subject to a limit on initial kinetic energy. 
Given the same kinetic energy, a lighter ball will go faster, and hence farther, if aerodynamic drag is ignored. 
Heavier cannonballs will have more inertia to counteract drag.
There is a balance between these two effects, which the optimizer seeks to find.

Here the static calculations are to find the mass and frontal area of the cannonball, given its radius. 
Then the ODE takes the mass and area as inputs and via Dymos can compute the total range. 
For demonstration purposes the trajectory is broken up into an ascent and descent phase, with the break being set up exactly at the apogee of the flight path. 

\small
```python
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import openmdao.api as om

import dymos as dm
from dymos.models.atmosphere.atmos_1976 import _USatm1976Data as USatm1976Data

# CREATE an atmosphere interpolant
english_to_metric_rho = om.unit_conversion('slug/ft**3', 'kg/m**3')[0]
english_to_metric_alt = om.unit_conversion('ft', 'm')[0]
rho_interp = interp1d(np.array(USatm1976Data.alt*english_to_metric_alt, dtype=complex), 
                      np.array(USatm1976Data.rho*english_to_metric_rho, dtype=complex), kind='linear')


class CannonballSize(om.ExplicitComponent):
    """
    Static calculations performed before the dynamic model
    """

    def setup(self):
        self.add_input(name='radius', val=1.0, 
                       desc='cannonball radius', units='m')
        self.add_input(name='density', val=7870., 
                       desc='cannonball density', units='kg/m**3')

        self.add_output(name='mass', shape=(1,), 
                       desc='cannonball mass', units='kg')
        self.add_output(name='area', shape=(1,), 
                       desc='aerodynamic reference area', units='m**2')

        self.declare_partials(of='*', wrt='*', method='cs')

    def compute(self, inputs, outputs):
        radius = inputs['radius']
        rho = inputs['density']

        outputs['mass'] = (4/3.) * rho * np.pi * radius ** 3
        outputs['area'] = np.pi * radius ** 2


class CannonballODE(om.ExplicitComponent): 
    """
    Cannonball ODE assuming flat earth and accounting for air resistance
    """

    def initialize(self): 
        self.options.declare('num_nodes', types=int)

    def setup(self): 
        nn = self.options['num_nodes']

        # static parameters
        self.add_input('mass', units='kg')
        self.add_input('area', units='m**2')

        # time varying inputs 
        self.add_input('alt', units='m', shape=nn)
        self.add_input('v', units='m/s', shape=nn)
        self.add_input('gam', units='rad', shape=nn)

        # state rates
        self.add_output('v_dot', shape=nn, units='m/s**2')
        self.add_output('gam_dot', shape=nn, units='rad/s')
        self.add_output('h_dot', shape=nn, units='m/s')
        self.add_output('r_dot', shape=nn, units='m/s')
        self.add_output('ke', shape=nn, units='J')

        # Ask OpenMDAO to compute the partial derivatives using complex-step 
        # with a partial coloring algorithm for improved performance
        self.declare_partials('*', '*', method='cs')
        self.declare_coloring(wrt='*', method='cs')

    def compute(self, inputs, outputs): 

        gam = inputs['gam']
        v = inputs['v']
        alt = inputs['alt']
        m = inputs['mass']
        S = inputs['area']

        CD = 0.5 # good assumption for a sphere
        GRAVITY = 9.80665 # m/s**2

        # handle complex-step gracefully from the interpolant
        if np.iscomplexobj(alt): 
            rho = rho_interp(inputs['alt'])
        else: 
            rho = rho_interp(inputs['alt']).real

        q = 0.5*rho*inputs['v']**2
        qS = q * S
        D = qS * CD
        cgam = np.cos(gam)
        sgam = np.sin(gam)
        outputs['v_dot'] = - D/m-GRAVITY*sgam
        outputs['gam_dot'] = -(GRAVITY/v)*cgam
        outputs['h_dot'] = v*sgam
        outputs['r_dot'] = v*cgam
        outputs['ke'] = 0.5*m*v**2

if __name__ == "__main__": 

    p = om.Problem()

    ###################################
    # Co-design part of the model, 
    # static analysis outside of Dymos
    ###################################
    static_calcs = p.model.add_subsystem('static_calcs', CannonballSize())
    static_calcs.add_design_var('radius', lower=0.01, upper=0.10, 
                                ref0=0.01, ref=0.10)

    p.model.connect('static_calcs.mass', 'traj.parameters:mass')
    p.model.connect('static_calcs.area', 'traj.parameters:area')

    traj = p.model.add_subsystem('traj', dm.Trajectory())
    # Declare parameters that will be constant across 
    # the two phases of the trajectory, so we can connect to it only once
    traj.add_parameter('mass', units='kg', val=0.01, dynamic=False)
    traj.add_parameter('area', units='m**2', dynamic=False)

    tx = dm.Radau(num_segments=5, order=3, compressed=True)
    ascent = dm.Phase(transcription=tx, ode_class=CannonballODE)
    traj.add_phase('ascent', ascent)

    ###################################
    # first phase: ascent
    ###################################
    # All initial states except flight path angle are fixed
    ascent.add_state('r', units='m', rate_source='r_dot', 
                     fix_initial=True, fix_final=False)
    ascent.add_state('h', units='m', rate_source='h_dot', 
                     fix_initial=True, fix_final=False)
    ascent.add_state('v', units='m/s', rate_source='v_dot', 
                     fix_initial=False, fix_final=False)
    # Final flight path angle is fixed (
    #     we will set it to zero so that the phase ends at apogee)
    ascent.add_state('gam', units='rad', rate_source='gam_dot', 
                     fix_initial=False, fix_final=True)    
    ascent.set_time_options(fix_initial=True, duration_bounds=(1, 100), 
                            duration_ref=100, units='s')

    ascent.add_parameter('mass', units='kg', val=0.01, dynamic=False)
    ascent.add_parameter('area', units='m**2', dynamic=False)

    # Limit the initial muzzle energy to create a well posed problem 
    # with respect to cannonball size and initial velocity
    ascent.add_boundary_constraint('ke', loc='initial', units='J',
                                   upper=400000, lower=0, ref=100000)

    ###################################
    # second phase: descent
    ###################################
    tx = dm.GaussLobatto(num_segments=5, order=3, compressed=True)
    descent = dm.Phase(transcription=tx, ode_class=CannonballODE)
    traj.add_phase('descent', descent )

    # All initial states and time are free so their 
    #    values can be linked to the final ascent values
    # Final altitude is fixed to 0 to ensure final impact on the ground
    descent.add_state('r', units='m', rate_source='r_dot', 
                      fix_initial=False, fix_final=False)
    descent.add_state('h', units='m', rate_source='h_dot', 
                      fix_initial=False, fix_final=True)
    descent.add_state('gam', units='rad', rate_source='gam_dot', 
                      fix_initial=False, fix_final=False)
    descent.add_state('v', units='m/s', rate_source='v_dot',
                      fix_initial=False, fix_final=False)
    descent.set_time_options(initial_bounds=(.5, 100), duration_bounds=(.5, 100),
                             duration_ref=100, units='s')
    
    descent.add_parameter('mass', units='kg', val=0.01, dynamic=False)
    descent.add_parameter('area', units='m**2', dynamic=False)

    # Link Phases (link time and all state variables)
    traj.link_phases(phases=['ascent', 'descent'], vars=['*'])

    # maximize range
    descent.add_objective('r', loc='final', ref=-1.0)

    p.driver = om.pyOptSparseDriver()
    p.driver.options['optimizer'] = 'SLSQP'
    p.driver.declare_coloring()

    # Finish Problem Setup
    p.setup()

    # Set Initial guesses for static dvs and ascent
    p.set_val('static_calcs.radius', 0.05, units='m')
    p.set_val('traj.ascent.t_duration', 10.0)

    p.set_val('traj.ascent.states:r', 
              ascent.interpolate(ys=[0, 100], nodes='state_input'))
    p.set_val('traj.ascent.states:h', 
              ascent.interpolate(ys=[0, 100], nodes='state_input'))
    p.set_val('traj.ascent.states:v', 
              ascent.interpolate(ys=[200, 150], nodes='state_input'))
    p.set_val('traj.ascent.states:gam', 
              ascent.interpolate(ys=[25, 0], nodes='state_input'), units='deg')

    # more intitial guesses for descent
    p.set_val('traj.descent.t_initial', 10.0)
    p.set_val('traj.descent.t_duration', 10.0)

    p.set_val('traj.descent.states:r', 
               descent.interpolate(ys=[100, 200], nodes='state_input'))
    p.set_val('traj.descent.states:h', 
              descent.interpolate(ys=[100, 0], nodes='state_input'))
    p.set_val('traj.descent.states:v', 
              descent.interpolate(ys=[150, 200], nodes='state_input'))
    p.set_val('traj.descent.states:gam', 
              descent.interpolate(ys=[0, -45], nodes='state_input'), units='deg')

    dm.run_problem(p, simulate=True, make_plots=True)

    fig, ax = plt.subplots()
    x0 = p.get_val('traj.ascent.timeseries.states:r', units='m')
    y0 = p.get_val('traj.ascent.timeseries.states:h', units='m')
    x1 = p.get_val('traj.descent.timeseries.states:r', units='m')
    y1 = p.get_val('traj.descent.timeseries.states:h', units='m')
    tab20 = plt.cm.get_cmap('tab20').colors
    ax.plot(x0,y0, marker='o', label='ascent', color=tab20[0])
    ax.plot(x1,y1, marker='o', label='descent', color=tab20[1])
    ax.legend(loc='best')
    ax.set_xlabel('range (m)')
    ax.set_ylabel('height (m)')
    fig.savefig('cannonball_hr.png', bbox_inches='tight')
```
\normalsize

The built-in plotting in Dymos will give time evolutions of all the time varying quantities.
For example, these are the trajectories for the range and height:

![Cannonball Solution: height vs time](cannonball_states_h.png){width=50%}
![Cannonball Solution: range vs time](cannonball_states_r.png){width=50%}

A more natural way to view the solution is to consider height vs range: 

![Cannonball Solution: height vs time](cannonball_hr.png){width=50%}

The parabolic trajectory is slightly skewed due to the effect of air resistance slowing down the cannonball so it is moving slower during the descent than the ascent. 


# Acknowledgements

Dymos was developed with funding from NASA's Transformational Tools and Technologies ($T^3$) Project.

# References

# Multidisciplinary Optimal Control Library

The goal of Dymos is to enable the optimization of subsystem designs which are tightly connected with each other as well as the operation of the overall system.

Dymos is a library for the optimal control of dynamic multidisciplinary systems.
While it can optimize typical optimal control problems, its key feature is the ability to optimize _systems_ in which a trajectory is just one part of the overall optimization.
Other optimization software frequently relies on the parameterization of the hardware models to, for instance, approximate the mass of an engine as a function of its thrust level.
Instead, Dymos allows you to impose higher-fidelity design considerations on these subsystems - from simple parameterized models to high-fidelity CFD models, and apply the resulting subsystem designs to the trajectory profile.

To do this, Dymos relies on ability of OpenMDAO to compute accurate derivatives very efficiently.
This capability enables users of Dymos to embed iterative procedures within their system dynamics.
While normally this would significantly impair performance, Dymos can optimize such systems with minimal performance degradation, freeing the user from reformulating their design specifically for the purposes of the optimal control implementation.

## Key Features

-   Employ Gauss-Lobatto collocation {cite}`herman1996direct` or the Radau Pseudospectral method {cite}`garg2011direct` to find optimal control for a dynamic system.
-   Find the optimal design of a system that can satisfy a variety of different trajectories.
-   Embed nonlinear solvers within the system dynamics.
-   Transform typical state variables into control variables (differential inclusion {cite}`Seywald1994`).
-   Use nonlinear solvers to satisfy the collocation constraints.
-   Single and multiple shooting within the same interface.
-   Leverage multiprocessing capabilities to improve performance

## Why Dymos?

There is no shortage of optimal control software based on pseudospectral approaches.
There are a number of other optimal control libraries that tackle similar kinds of problems, such as OTIS4 {cite}`paris_riehl_sjauw_2006`, GPOPS-II {cite}`patterson2014gpops`,and CASADI {cite}`Andersson2018`.

Given the amount of software existing in this space, why did we develop Dymos?

We wanted to use optimal control in the context of multidisciplinary optimization.
With existing state-of-the-art tools, this would mean wrapping one of these existing tools and passing between different disciplines.
But for performance reasons, we also wanted to be able to pass derivatives between our different disciplines, including trajectory design.
This approach drives us towards a "monolithic" optimization problem, where a single objective is optimized that encompasses information from all of the subsystems.
The "collaborative" optimization approach that optimizes each subsystem with a minimal amount of information passing our subsystem boundaries is generally slower and fails to find solutions as good as the monolithic approach.
The first objective, therefore, was to develop optimal control software that has the capability to provide accurate, analytic derivatives.

Many state-of-the-art optimal control software packages rely on finite-differencing to estimate derivatives needed by the optimizer to solve the problem, although that is beginning to change.
This inherently couples the accuracy of the derivatives to the scaling of the problem.
We'd like to use analytic derivative calculations to better decouple this interaction.
Even those software packages which employ analytic derivatives generally use a forward-differentiation approach.
The work of Hwang and Martins {cite}`hwang2018b` demonstrated how to develop a framework for accurate derivative calculation, including analytic derivatives and both complex-step and finite-difference approximations as fallbacks.
Their approach gives us a few key capabilities:

- Adjoint differentiation which can be more efficient for pure shooting-methods in optimal control, where the number of constraints/objectives is far fewer than the number of design variables.
- The ability to compute derivatives across complex iterative systems _without the need to reconverge the system_.
- The ability to provide a general, bidirectional derivative coloring system {cite}`gray2019coloring` which can minimize the computational effort required to compute the sensitivies of the outputs with respect to the inputs.

In addition to making optimal control more performant for use in multidisciplinary optimization, we were keen to study what sort of work these capabilities could enable.
Embedding iterative systems within the optimization, in particular, is generally avoided for performance reasons.
But with the state-of-the-art differentiation approach of OpenMDAO, built upon the work of Martins and Hwang, we can embed complex implicit systems with minimal impact on performance.
This enables more efficient optimization via differential inclusion {cite}`Seywald1994`, and allows us to employ shooting methods within the pseudospectral framework.

Some developers involved in Dymos are involved with NASA's legacy optimal control software, OTIS.
The general approach used by Dymos is similar to that of OTIS (trajectories divided into time portions called Phases, dynamic controls and static parameters, and both bound constraints as well as nonlinear boundary constraints and path constraints are all notions carried over from OTIS.
OTIS was pioneering software and offers some great capabilities, but it also lacks a lot of desirable modern features that have been developed by the community since its inception over thirty years ago.
Dymos features a more modular way for users to define their dynamics, additional pseudospectral methods, and better differentiation approaches for more reliable convergence.

## Citation

See our [overview paper](https://joss.theoj.org/papers/10.21105/joss.02809) in the Journal of Open Source Software

If you use Dymos in your work, please cite:
```
@article{Falck2021,
  doi = {10.21105/joss.02809},
  url = {https://doi.org/10.21105/joss.02809},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {59},
  pages = {2809},
  author = {Robert Falck and Justin S. Gray and Kaushik Ponnapalli and Ted Wright},
  title = {dymos: A Python package for optimal control of multidisciplinary systems},
  journal = {Journal of Open Source Software}
}
```

## References

```{bibliography}
:filter: docname in docnames
```
# Installing Dymos

## Basic installation

Dymos can be installed from pip using

~~~bash
python -m pip install dymos
~~~

The cutting-edge development version of Dymos may be installed with

~~~bash
python -m pip install git+https://github.com/OpenMDAO/dymos.git
~~~

## Installation (developer mode)

If you intend to modify Dymos, build the documentation, or run the tests locally, you should add `[all]` to the PIP installation to include some documentation and testing-specific installation dependencies.

After cloning out Dymos to a local directory, installing using the `-e` flag installs Dymos in "developer mode."
This removes the need to reinstall Dymos after changes are made.

~~~bash
git clone https://github.com/OpenMDAO/dymos.git ./dymos.git
python -m pip install -e dymos.git[all]
~~~

## Installing Advanced Dependencies

Many of the simple dependencies for Dymos are installed via pip, however a few more dependencies are required in some more advanced use-cases.

### pyoptsparse (recommended)

Dymos will work "out-of-the-box" using OpenMDAO's `ScipyOptimizeDriver` for many of the more simple problems.
For more complex problems, the implicit optimization techniques used by Dymos can easily generate problems with hundreds or even thousands of design variables and constraints.
When the nonlinear optimization problems generated by Dymos become that complicated, the free optimizers available via `ScipyOptimizeDriver` can struggle to converge.

For this reason, the authors tend to use the open-source [IPOPT](https://coin-or.github.io/Ipopt/) optimizer or the proprietary [SNOPT](https://ccom.ucsd.edu/~optimizers/solvers/snopt/) optimizer.
Both of these optimizers are able to capitalize on the _sparse_ nature of the optimal control problems generated by Dymos, significantly reducing the time and memory required to solve problems.

For uniform access to a variety of optimizers, including `SLSQP`, `IPOPT`, and `SNOPT`, Dymos uses OpenMDAO's `pyOptSparseDriver`, which interfaces to these external optimizers via [pyoptsparse](https://github.com/mdolab/pyoptsparse), produced by MDOLab at the University of Michigan.

Users on OS X and Linux systems can use our script, available at [(https://github.com/OpenMDAO/build_pyoptsparse)](https://github.com/OpenMDAO/build_pyoptsparse) to build and install pyoptsparse with support for IPOPT.
The script will also handle support for SNOPT if the user has access to the SNOPT source code.

Installing these dependencies on Windows-based systems requires [Intel Fortran](https://software.intel.com/content/www/us/en/develop/tools/oneapi/components/fortran-compiler.html).
Free Fortran compilers on Windows, at the time of this writing, are not compatible with the Microsoft ABI used by Python on Windows.

### mpi4py and petsc4py

Together, `mpi4py` and `petsc4py` enable parallel processing in Dymos via MPI.
The easiest way to install `mpi4py` is to use the [anaconda](https://www.anaconda.com/products/individual) python environment.

```
    conda create -n py38 python=3.8;
    conda activate py38;
    conda config --add channels conda-forge;
    conda install mpi4py
    conda install petsc4py
```

MPI capability is generally only necessary in some very niche use-cases.
# Phases

- [Phases](phases.ipynb)
- [Segments](segments.ipynb)
- [Variables](variables.ipynb)
- [Constraints](constraints.ipynb)
- [Objective](objective.ipynb)
- [Time Series](timeseries.ipynb)
# Dymos by Example

The goal of these examples is to walk users through the process of formulating an optimal control problem and solving it using Dymos.

In working through these examples, we'll try to emphasize the following process:

1.  Formulate the optimal control problem
2.  Identify state and control variables, and the ordinary differential equations (ODE) which govern the dynamics.
3.  Build the ODE as an OpenMDAO system.
4.  Test the evaluation of the ODE.
5.  Define the partial derivatives for the ODE.
6.  Test the partial derivatives of the ODE against finite-difference or (preferably) complex-step approximations.


## Prerequisites

These examples assume that the user has a working knowledge of the following:

-   Python
-   The numpy package for numerical computing with Python
-   The matplotlib package for plotting results

Some examples may require additional optional packages that are noted in those notebooks.
# Dymos Reference API

- [The Dymos run_problem function](run_problem)
- [The Phase API](phase_api)
- [The Transcriptions API](transcriptions_api)
- [The Trajectory API](trajectory_api)
# Optimal Control Transcriptions

Solving optimal control problems using classical techniques is typically a matter of finding a continuous control function that moves a system to the desired state while extremizing some objective (minimizing time, fuel, cost, etc.).

To solve these problems directly using nonlinear optimizers, we have to find a way to discretize this continuous problem into a form that can be solved by a nonlinear optimizer. This process is called _transcription_.

Dymos supports two forms of direct transcription: [collocation](getting_started:collocation:what_is_collocation) and explicit shooting.

Collocation-based optimal control methods, also known as pseudospectral methods, are implicit techniques.
The entire state and control history of the system is proposed at some discrete set of points in time.
We can then calculate how accurately the proposed state and control history obeyed the ODE governing the system dynamics - this is a residual we call the _defects_.
By iteratively varying the states, controls, and elapsed time of the trajectory we can reduce the defects to zero, meaning that the proposed state history is a solution to the ODE of the system, with the given control history and time of propagation.

In explicit shooting, an explicit numerical integration technique is used to propagate the initial state over the given time duration subject to the given controls. The process of integrating the trajectory itself eliminates the defects.

## Differences between collocation and explicit shooting

### Design Variables

For implicit collocation techniques, the design variables include:
- the states at some set of points along the trajectory (the initial and final ones may not be included if they're fixed)
- the controls at discrete points in time
- the elapsed time of the trajectory.

For explicit shooting techniques, the design variables include:
- the initial state (assuming it's not fixed, and assuming a single integration interval (_single shooting_))
- the initial state in each integration interval (for _multiple shooting_)
- the control values at some set of discrete points in time
- the elapsed time of the trajectory.

The design variable vector can be considerably larger for the implicit collocation techniques, depending on the size of the state vector and the number of discrete points in time.

### Constraints

For implicit collocation techniques, the constraints include:
- the defect constraints for each state at various discrete points in time throughout the trajectory
- continuity constraints that ensure the state values remain continuous (for `compressed=False`)
- continuity constraints on the control values (for `compressed=False`) and (optionally) their rates
- any other path or boundary constraints imposed

For explicit shooting techniques, the constraints include:
- continuity constraints that ensure the state values remain continuous (assuming multiple integration intervals (_multiple shooting_))
- continuity constraints on the control values (for `compressed=False`) and (optionally) their rates
- any other path or boundary constraints imposed

The constraint vector can be considerably larger for implicit collocation techniques.
Again this depends on the size of the state vector and the number of points into which the trajectory has been discretized.

We should also note that for the implicit collocation techniques, the state values at the start and the end of the trajectory are design variables.
This means that these values can be trivially _fixed_ (removed from the design variable vector) or bound using simple bounds on the corresponding design variable.

**The only way to impose a constraint on the path or the final value of a state in a shooting method is with a nonlinear path or boundary constraint.**

### Performance

Due to the fact that the design variable and constraint vectors are considerably larger for implicit collocation problems, it might be logical to conclude they they're slower than explicit shooting for solving a corresponding problem, but this is not the case.

In a single shooting method, the value of some state later in time is (at least potentially) a function of everything that happened before it.
The corresponding _jacobian matrix_ of derivatives for these state values along a trajectory is lower triangular.

For the implicit collocation techniques, the state in one integration _segment_ is not a function of anything that happens outside of that integration segment.
The jacobian matrix is much more sparse.
This increased sparsity helps in two ways.

First, we can determine the derivatives much more efficiently.
Consider a finite differencing technique.
If we know that a state value in the first integration segment of a trajectory only impacts that segment, and there are 100 segments, then we can simultaneously perturb a state value in each of the 100 segments and compare the resulting defect vector to the nominal defect value to compute the sensitivity via finite difference.
This is a massive performance gain that's not possible with a single shooting technique, though multiple shooting can help this.

Second, some optimizers can capitalize on the sparsity of jacobian matrices to provide significantly improved performance vs those which operate only on dense matrices.

Because of these factors, **collocation techniques can be orders of magnitude faster than their explicit shooting counterparts**.

### Robustness

While faster, implicit techniques impose far more defect constraints on the problem.
In some cases, scaling the optimization problem can be a challenge and convergence is poor.
For instance, if the thrust of a rocket engine is taken from experimental data, it can be extremely noisy.
Without first smoothing the data, an implicit simulation of the resulting trajectory may be difficult.
In this case, explicit intgration methods can power through and provide _an_ answer, even if it is subject to some amount of error.

**Explicit integration, at least using a fixed-step form, doesn't fail to provide an answer - but it's important that the user verify its accuracy.**

On the other hand, explicit integration is subservient to the ODE, and sometimes this feature is problematic.
When following the prescribed control profile can take the system into a region where the dynamics become singular, the resulting errors will cause the optimization to fail.
The [minimum time-to-climb](examples:minimum_time_climb) problem is a classic example of this.
The equations of motion used in this problem are singular in vertical flight - there is division by zero in the derivatives of the equations of motion in this case.
It is relatively easy to prescribe a profile of the angle-of-attack history (the control) that sends the aircraft into vertical flight when the integration is required to follow it.

Conversely, the collocation techniques decouple the proposed control history and the trajectory.
Rather than being governed by the control, the flight path angle at various times throughout the trajectory is itself a design variable, and can be bound to values such that avoid the singularities.
The optimizer will then work to make the control history compatible with the corresponding state trajectory, but the two are only compatible when the defect constraints are satisfied.
# Optimal Control


Optimal control implies the optimization of a dynamical system.  Typically this takes the form
of a trajectory in which the *states* of the system evolve with time.  The evolution of the states
$\left(\bar{x}\right)$ are typically governed by an ordinary differential equation (ODE) or
a differential algebraic equation (DAE).  In Dymos we characterize all dynamics as an (ODE),
although solving systems of DAEs is also possible.

\begin{align}
  \dot{\bar{x}} = f_{ode}(\bar{x},t,\bar{u},\bar{d})
\end{align}

# Controls and Parameters

In order to impact the behavior of such systems we need *controls*.
Controls may be allowed to vary with time, such as the angle-of-attack of an aircraft during flight.
We refer to these as _dynamic_ controls $\left(\bar{u}\right)$.
Other controllable parameters might be fixed with time, such as the wingspan of an aircraft.
We refer to these as _parameters_ $\left(\bar{d}\right)$, although in the literature they may also be referred to as static controls.
The endpoints in time, state values, control values, and design parameter values define the independent variables for our optimization problem.
In Dymos, we discretize these variables in time to convert a continuous-time optimal control problem into a nonlinear programming (NLP) problem.

# Constraints

When optimizing problems there will always be constraints on the system.
Some of these constraints can be characterized as bounds on our independent variables.
For instance, fixing the initial conditions of a trajectory can be accomplished by bounding the initial states and time to the desired value.

\begin{align}
  \mathrm{Time:}& \qquad {t}_{lb} \leq t \leq {t}_{ub} \\
  \mathrm{State \, Variables:}& \qquad \bar{x}_{lb} \leq \bar{x} \leq \bar{x}_{ub} \\
  \mathrm{Dynamic \, Controls:}& \qquad \bar{u}_{lb} \leq \bar{u} \leq \bar{u}_{ub} \\
  \mathrm{Design \, Parameters:}& \qquad \bar{d}_{lb} \leq \bar{d} \leq \bar{d}_{ub} \\
\end{align}

Other times we may want to constrain an output of our system at a specific point along the trajectory (boundary constraints) or along the entire trajectory (path constraints).
These are so-called nonlinear constraints because they *may* be nonlinear functions of our independent variables.

\begin{align}
\mathrm{Initial \, Boundary \, Constraints:}& \qquad \bar{g}_{0,lb} \leq g_{0}(\bar{x}_0,t_0,\bar{u}_0, \bar{d}) \leq \bar{g}_{0,ub} \\
\mathrm{Final \, Boundary \, Constraints:}& \qquad \bar{g}_{f,lb} \leq g_{f}(\bar{x}_f,t_f,\bar{u}_f, \bar{d}) \leq \bar{g}_{f,ub} \\
\mathrm{Path \, Constraints:}& \qquad \bar{p}_{f,lb} \leq p_{f}(\bar{x},t,\bar{u},\bar{d}) \leq \bar{p}_{f,ub} \\
\end{align}

In practice, bound constraints typically govern the limits which the optimizer must observe when providing values of our independent variables.  
Nonlinear constraints, on the other hand, are constraints applied to the outputs of a system.  
To put another way, bound constraints are always observed, while nonlinear constraints are not necessarily observed until the optimizer has finished solving the optimization problem.

# The Objective

Dymos may be used to both simulate and optimize dynamical systems.
The phase construct is generally used in optimization contexts.
Within each phase, the user can set the objective:

\begin{align*}
  \mathrm{J} = f_{obj}(\bar{x},t,\bar{u},\bar{d})
\end{align*}

As with constraints, the objective may be any output within the Phase.  Phases can also be
incorporated into larger models wherein the objective is defined in some subsystem outside of the
phase.  In this case, the standard OpenMDAO method `add_objective` can be used.

# The Overall Optimization Problem

The optimization problem as defined by Dymos can thus be stated as:

\begin{align*}
    \mathrm{Minimize}& \qquad \mathrm{J} = f_{obj}(\bar{x},t,\bar{u},\bar{d}) \\
    \mathrm{Subject \, to:}& \\
    \mathrm{Dynamic \, Constraints:}& \qquad \dot{\bar{x}} = f_{ode}(\bar{x},t,\bar{u},\bar{d}) \\
    \mathrm{Time:}& \qquad {t}_{lb} \leq t \leq {t}_{ub} \\
    \mathrm{State \, Variables:}& \qquad \bar{x}_{lb} \leq \bar{x} \leq \bar{x}_{ub} \\
    \mathrm{Dynamic \, Controls:}& \qquad \bar{u}_{lb} \leq \bar{u} \leq \bar{u}_{ub} \\
    \mathrm{Design \, Parameters:}& \qquad \bar{d}_{lb} \leq \bar{d} \leq \bar{d}_{ub} \\
    \mathrm{Initial \, Boundary \, Constraints:}& \qquad \bar{g}_{0,lb} \leq g_{0}(\bar{x}_0,t_0,\bar{u}_0, \bar{d}) \leq \bar{g}_{0,ub} \\
    \mathrm{Final \, Boundary \, Constraints:}& \qquad \bar{g}_{f,lb} \leq g_{f}(\bar{x}_f,t_f,\bar{u}_f, \bar{d}) \leq \bar{g}_{f,ub} \\
    \mathrm{Path \, Constraints:}& \qquad \bar{p}_{f,lb} \leq p_{f}(\bar{x},t,\bar{u},\bar{d}) \leq \bar{p}_{f,ub}
\end{align*}
# Frequently Asked Questions

- [How do I add an ODE output to the timeseries outputs?](add_ode_output_to_timeseries.md)
- [How do I connect a scalar input to the ODE?](connect_scalar_parameters_to_ode.md)
- [How do connect the outputs of an upstream analysis as inputs to Dymos?](upstream_analysis.md)
- [How do connect the outputs of Dymos to a downstream analysis?](downstream_analysis.md)
- [How do I run two phases parallel-in-time?](tandem_phases.md)
- [How can I more efficiently use finite-differenced components in the ODE?](use_partial_coloring.md)
- [How can I debug models when things go wrong?](debugging.md)