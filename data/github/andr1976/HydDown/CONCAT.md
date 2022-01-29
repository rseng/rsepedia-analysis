---
title: HydDown
subtitle: User guide and technical reference
author: Anders Andreasen
titlepage: true
toc-own-page: true
book: true
reference-section-title: References
bibliography: references.bib
listings: True
---

# Introduction
HydDown is an open source python3 tool for calculation of hydrogen (or other pure gas phase species) vessel/container depressurization and filling.
The HydDown logo shown in [@Fig:logo] visualizes the key parameters and transport phenomena during gas vessel filling or discharging.
The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P).
This is caused by change in fluid inventory (density) due to flow of gas either in or out of the vessel.
Further, heat is transferred from or to the surroundings via convective heat transfer on the in- and outside of the vessel with heat being conducted through the vessel wall.

![HydDown logo](docs/img/Sketch.png){#fig:logo}

Running the code is as simple as: 

    python main.py input.yml

where `main.py` is the main script and `input.yml` is the input file in Yaml syntax. 

## Background
HydDown started as a small spare-time project for calculation of vessel filling and depressurization behaviour.
This was mainly to demonstrate that although perceived as a very tedious and difficult task to write your own code for such an apparently complex problem, actually only a fairly limited amount of code is necessary if you have a good thermodynamic backend.
The code has evolved to a point where additional features will increase the complexity to the point where it is no longer a simple tool. 

A few choices has been made to keep things simple:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Only pure substances are considered (limited multi-component capabilities are included)
- Gas phase only
- No temperature stratification in the gas phase
- No temperature gradient through vessel wall 
- Heat transfer is modelled as constant or simplified using empirical
  correlations

These choices make the problem much simpler to solve:
First of all, the the pure substance Helmholtz energy based equation of state (HEOS) in `coolprop` offers a lot of convenience in terms of the property pairs/state variables that can be set independently.
Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required.
That being said, the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.

## Getting the software
The source code can be obtained either from GitHub (via `git` or via the latest tar-ball release) or via **pip** . No packaged releases have currently been planned for **conda**.  

The main branch is located here:

[`https://github.com/andr1976/HydDown`](https://github.com/andr1976/HydDown)

Clone the repo by: 

    git clone https://github.com/andr1976/HydDown.git

Running from source via `git` the dependencies must be installed manually from the repo root dir:

    pip install -r requirements.txt

Installation of latest release via **pip** also installe dependecies automatically:

    pip install hyddown

In case `pip` links to a v2.7 of python you will get an error. If so try the following:

    python3 -m pip install hyddown

where python3 is the symlink or full path to the python3 executable installed on your system. 

## Requirements

- [Python](http://www.python.org) (3.8 - at least python3, not 3.9, see below)
- [Numpy](https://numpy.org/)
- [matplotlib](https://matplotlib.org/)
- [Coolprop (6.4.1)](http://www.coolprop.org/)
- [cerberus](https://docs.python-cerberus.org/en/latest/)
- [PyYaml](https://pypi.org/project/PyYAML/)
- [pandas](https://pandas.pydata.org/)
- [Scipy](https://www.scipy.org/)
- [tqdm](https://tqdm.github.io/)

It has been reported that the CoolProp `pip` package does not run out of the box with python 3.9. Hence, python 3.8 is recommended. 

The script is running on Windows 10 x64, with stock python installation from python.org and packages installed using pip.
It should also run on Linux (it does on an Ubuntu image on GitHub) or in any conda environment as well, but I haven't checked.

## Testing 
Although testing is mainly intended for automated testing (CI) during development using github actions, testing of the installed package can be done for source install by:

    python -m pytest 

run from the root folder. In writing 33 test should pass.

For the package installed with pip navigate to the install directory `../site-packages/hyddown` and run:

    python -m pytest 

## Units of measure
The SI unit system is adapted for this project.
The following common units are used in the present project and this also applies to the units used in the input files:

Property | Unit | Comment 
----    | ----  | ----
Temperature | K | $^\circ$ C is used in plots
Pressure | Pa    | bar is used in plots
Mass | kg |
Volume | m$^3$ |
Time | s | 
Energy | J |
Duty/power | W 
Length | m
Area | m$^2$
Heat flux | W/m$^2$
Heat transfer coefficient | W/(m$^2$ K)
Density | kg/m$^3$
Heat capacity | J/(kg K)

: Unit system {#tbl:units}

As will be noted when presenting the equations implemented in the code, some of the equations utilise different units than the ones listed in [@tbl:units].
However, it is important to note that unit conversions are built in to the methods implemented, so the user shall not worry about unit conversion.  

## Credit 
In the making of this document a great deal of material has been sourced (and modified) from a good colleague's M.Sc. thesis [@iskov], from co-published papers [@Bjerre2017][@safety4010011] and from on-line material published under permissive licenses (with proper citation).
Further, the making of this project would not have possible without the awesome [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999].
I am also thankful for enlightening discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark)  and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurization, nozzle flow and heat transfer considerations.

The present document is typeset using Markdown + [pandoc](https://pandoc.org/) with the [Eisvogel](https://github.com/Wandmalfarbe/pandoc-latex-template) template. 

## License

MIT License

Copyright (c) 2021 Anders Andreasen

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

# Usage 
## Basic usage
From the source repository, running the code is as simple as: 

    python src/main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax. 

The Yaml input file is edited to reflect the system of interest. 

Further, a usable copy of a main script is installed in the python installation Scripts/ folder:

    /python3x/Scripts/hyddown_main.py

or from GitHub/in the source folder:

    HydDown/Scripts/hyddown_main.py

Various example files are included for inspiration under site-packages:

    /site-packages/hyddown/examples/

or in the git/source directory tree

    HydDown/src/hyddown/examples/

## Demos
A few demonstrations of the codes capability are available as [`streamlit`](https://streamlit.io/) apps:

- Vessel gas pressurisation/depressurization for various gases at [`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py)
- Response to external fire of vessel equipped with a PSV using the Stefan-Boltzmann fire equation at[`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_sbapp.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_sbapp.py)
- Vessel depressurization including slow opening of valve at [`https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_cvapp.py`](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_cvapp.py)

The various streamlit apps are also included in the scripts folder for running the scripts on a local machine. 

## Calculation methods 
The following methods are implemented:

- Isothermal: constant temperature of the fluid during depressurization (for a very slow process with a large heat reservoir)
- Isenthalpic: constant enthalpy of the fluid, no heat transfer with surroundings, no work performed by the expanding fluid for depressurization, no work added from inlet stream
- Isentropic: constant entropy of the fluid, no heat transfer with surroundings, PV work performed by the expanding fluid
- Isenergetic: constant internal energy of the fluid
- Energy balance: this is the most general case and is based on the first law of thermodynamics applied to a flow process.

For `isothermal`/`isenthalpic`/`isentropic`/`isenergetic` calculations the minimal input required are:

- Initial conditions (pressure, temperature)
- vessel dimensions (ID, length)
- valve parameters (Cd, diameter, backpressure)
- Calculation setup (time step, end time)
- Type of gas

If heat transfer is to be considered the calculation type `energybalance` is required. A few options are possible:

- Fixed U: U-value, i.e. overall heat transfer coefficient, is required and ambient temperature required.
Based on the ambient (external) temperature, the fluid temperature, and the U-value, the heat flux Q is calculated for each time step. 
In this calculation method the temperature of the vessel wall is not calculated:
 
$$ Q = UA(T_{ambient} - T_{fluid})  $$

- Fixed Q: The external heat flux, Q, applied to the fluid is specified and constant.
The ambient temperature is not required.
Further, the vessel wall temperature is redundant and not calculated.
- Specified h: The external heat transfer coefficient must be specified and either the internal heat transfer coefficient is provided or calculated from the assumption of natural convection from a vertical cylinder at high Gr number.
Ambient temperature is required. Using this method the wall temperature is calculated from an energy balance over the vessel wall taking in and out flux to/from the external ambient plenum as well as heat flux to/from the fluid inventory into account.
- Fire: The Stefan-Boltzmann equation is applied for estimating the external heat duty.
The fire heat flux depends on the vessel wall temperature, and the wall temperature is continuously updated as in the method with specified h.

More elaborate description of the required input for the different calculation types are provided in [@Sec:input].

## Script 
HydDown comes with a script which can be used as a command-line tool to start calculations.
If an input filename (path) is given as the first argument, this input file will be used.
If no arguments are passed, the script will look for an input file with the name `input.yml`.
The content of the main script is shown below:


~~~ {.Python}
import yaml
import sys
from hyddown import HydDown

if __name__ == "__main__":
    if len(sys.argv) > 1:
        input_filename = sys.argv[1]
    else:
        input_filename = "input.yml"

    with open(input_filename) as infile:
        input = yaml.load(infile, Loader=yaml.FullLoader)


    hdown=HydDown(input)
    hdown.run()
    hdown.verbose=1
    hdown.plot()
~~~

## Module import 
To use HydDown simple import the main calculation class `HydDown`.

~~~ {.Python}
from hyddown import HydDown
~~~

## Input file examples
When using HydDown a dictionary holding all relevant input in order for HydDown to do vessel calculations shall be provided when the class is initialized. One way is to read an input file. For HydDown a Yaml format is chosen, but if JSON is a preference this should also work with the Yaml parser being substituted with a JSON parser.

An example of a minimal file for an isentropic vessel depressurization (no heat transfer) is shown below

~~~ {.Yaml}
vessel:
  length: 1.524
  diameter: 0.273
initial:
  temperature: 388.0
  pressure: 15000000.
  fluid: "N2" 
calculation:
  type: "isentropic"
  time_step: 0.05
  end_time: 100.
valve:
  flow: "discharge"
  type: "orifice"
  diameter: 0.00635
  discharge_coef: 0.8
  back_pressure: 101300.
~~~

A more elaborate example which includes heat transfer and with validation data (some data points dropped for simplicity) included: 

~~~ {.Yaml}
vessel:
  length: 1.524
  diameter: 0.273
  thickness: 0.025
  heat_capacity: 500
  density: 7800. 
  orientation: "vertical"
initial:
  temperature: 288.0
  pressure: 15000000.
  fluid: "N2" 
calculation:
  type: "energybalance"
  time_step: 0.05
  end_time: 100.
valve:
  flow: "discharge"
  type: "orifice"
  diameter: 0.00635
  discharge_coef: 0.8
  back_pressure: 101300.
heat_transfer:
  type: "specified_h"
  temp_ambient: 288.
  h_outer: 5 
  h_inner: 'calc' 
validation:
  temperature:
    gas_high:
      time: [0.050285, ... , 99.994]
      temp: [288.93, ... ,241.29] 
    gas_low:
      time: [0.32393, ..., 100.11]
      temp: [288.67, ... ,215.28]
    wall_low:
      time: [0.32276, ... , 100.08]
      temp: [288.93, ... ,281.72]
    wall_high:
      time: [0.049115, ... ,100.06]
      temp: [289.18, ... ,286.09]
  pressure:
    time: [0.28869, ... , 98.367]
    pres: [150.02, ... ,1.7204]
~~~

## Input fields and hierarchy {#sec:input}
In the following the full hierarchy of input for the different calculation types is summarised.

At the top level the following fields are accepted, with the last being optional and the second last dependant on calculation type:

~~~ {.Yaml}
initial: mandatory
vessel: mandatory
calculation: mandatory
valve: mandatory
heat_transfer: depends on calculation type
validation: optional
~~~

### Calculation

The subfields under `calculation`, with value formats and options are:

~~~ {.Yaml}
calculation:
  type: "isothermal", "isentropic", "isenthalpic", "constantU", "energybalance"
  time_step: number
  end_time: number
~~~

The simulation end time is specified as well as the fixed time step used in the integration of the differential equations to be solved.
The four main calculation types are shown as well.

### Vessel

~~~ {.Yaml}
vessel:
  length: number, mandatory
  diameter: number, mandatory
  thickness: number, required when heat transfer is calculated
  heat_capacity: number, required when heat transfer is calculated
  density: number, required when heat transfer is calculated 
  orientation: string, required when heat transfer is calculated
~~~

### Initial
~~~ {.Yaml}
initial:
  temperature: number, mandatory
  pressure: number, mandatory
  fluid: string, mandatory , e.g. "N2", "H2" 
~~~

### Valve
The `valve` field determines the mode of the process is it for depressurization/discharge from the vessel using the value `discharge` or if mass flow is entering the vessel, filling it up/increasing pressure with the value `filling`.

Different types of mass flow devices can be specified:

- Restriction `orifice`
- Relief valve/pressure safety valve (only for `discharge` not `filling`)
- Control valve (`controlvalve`)
- A specified mass flow (`mdot`)

For the physical devices a `back_pressure` is required for the flow calculations.
The value of the `back_pressure` **is also used to specify the reservoir pressure when the vessel is filled**.
See also [@Sec:flow] for details about the calculation of flow rates. 

~~~ {.Yaml}
valve:
  flow: string, mandatory "discharge" or "filling"
  type: string, mandatory "orifice", "controlvalve", "psv", "mdot"
  back_pressure: number, required for type "orifice", "controlvalve" and "psv"
  diameter: number, required for "orifice" and "psv"
  discharge_coef: number, required for "orifice" and "psv"
  Cv: number, required for "control_valve"
  characteristic:  string, optional for "control_valve"
  time_constant: number, optional for "control_valve"
~~~

For the calculations using a control valve as mass transfer device a linear rate for the actuator can be specified using the `time_constant` field. An associated valve characteristic is associated in order to translate the linear actuator rate to an opening Cv (of max. Cv). Three caracteristics can be specified:

- Linear
- Equal percentage (exponential)
- Quick opening/fast (square root)

### Heat transfer
For more information about the actual estimation of heat transfer see also [@Sec:heat]. For the case of `fire` heat input predetermined parameters for the Stefan-Boltzmann fire equation are used to calculate different background heat loads cf. [@Tbl:fire]

| Source          |  Fire type      |       Heat load ($kW/m^2$) |
|-----------------|-----------------|----------------------------|
| API521          | Pool            |                   60       |
| API521          | Jet             |                  100       |
| Scandpower      | Pool            |                  100       |
| Scandpower      | Jet             |                  100       |

: Fire heat loads {#tbl:fire}

~~~ {.Yaml}
heat_transfer:
  type: string, mandatory, "specified_h", "specified_Q", "specified_U", "s-b"
  temp_ambient: number, required for type "specified_h", "specified_U"
  h_outer: number, required for type "specified_h"
  h_inner: number or 'calc', required for type "specified_h"
  U_fix: number, required for type "specified_U"
  Q_fix: number, required for type "specified_Q"
  fire: string, required for type "s-b", 'api_pool','api_jet','scandpower_pool','scandpower_jet'
  D_thoat: number, required for flow type "filling", set to vessel ID as a starting point
~~~

### Validation
In order to plot measured data against simulated data the field `validation` is included.

The following arrays (one or more) are supported in the sub-field `temperature`:

- `gas_high`: highest measured values of the bulk gas 
- `gas_low`: lowest measured values of the bulk gas 
- `gas_mean`: average measured values of the bulk gas
- `wall_mean`: average measured values of the vessel (inner) wall
- `wall_high`  i.e. highest measured values of the vessel (inner) wall
- `wall_low`  i.e. lowest measured values of the vessel (inner) wall

For each of the above fields arrays for `time`and `temp`shall be supplied with matching length.

A field for measured vessel pressure is also possible, where `time`and `pres` shall be identical length arrays. See also example below.

~~~ {.Yaml}
validation:
  temperature:
    gas_high:
      time: [0.050285, ... , 99.994]
      temp: [288.93, ... ,241.29] 
    gas_low:
      time: [0.32393, ..., 100.11]
      temp: [288.67, ... ,215.28]
    wall_low:
      time: [0.32276, ... , 100.08]
      temp: [288.93, ... ,281.72]
    wall_high:
      time: [0.049115, ... ,100.06]
      temp: [289.18, ... ,286.09]
  pressure:
    time: [0.28869, ... , 98.367]
    pres: [150.02, ... ,1.7204]
~~~


# Theory
In this chapter the basic theory and governing equations for the model implementation in HydDown is presented.
The following main topics are covered: 

- thermodynamics
- mass transfer
- heat transfer

## Thermodynamics

### Equation of state
The equation of state used by HydDown is the Helmholtz energy formulation as implemented in [CoolProp](http://www.coolprop.org/) [@doi:10.1021/ie4033999].
Most of the text in the present section has been sourced from the CoolProp documentation to be as accurate and true to the source as possible.
The Helmholtz energy formulation is a convenient construction of the equation of state because all the thermodynamic properties of interest can be obtained directly from partial derivatives of the Helmholtz energy.

It should be noted that the EOS are typically valid over the entire range of the fluid, from subcooled liquid to superheated vapor, to supercritical fluid.
In general, the EOS are based on non-dimensional terms $\delta$ and $\tau$, where these terms are defined by

$$ \delta=\rho/\rho_c $$
$$ \tau=T_c/T $$
    
where $\rho_c$ and $T_c$ are the critical density of the fluid if it is a pure fluid.
For pseudo-pure mixtures, the critical point is typically not used as the reducing state point, and often the maximum condensing temperature on the saturation curve is used instead.

The non-dimensional Helmholtz energy of the fluid is given by

$$ \alpha=\alpha^0+\alpha^r $$
    
where $\alpha^0$ is the ideal-gas contribution to the Helmholtz energy, and $\alpha^r$ is the residual Helmholtz energy contribution which accounts for non-ideal behavior.
For a given set of $\delta$ and $\tau$, each of the terms $\alpha^0$ and $\alpha^r$ are known.
The exact form of the Helmholtz energy terms is fluid dependent, but a relatively simple example is that of Nitrogen, which has the ideal-gas Helmholtz energy of

$$  \alpha^0=\ln\delta+a_1\ln\tau+a_2+a_3\tau+a_4\tau^{-1}+a_5\tau^{-2}+a_6\tau^{-3}+a_7\ln[1-\exp(-a_8\tau)] $$
    
and the non-dimensional residual Helmholtz energy of

$$    \alpha^r=\sum_{k=1}^{6}{N_k\delta^{i_k}\tau^{j_k}}+\sum_{k=7}^{32}{N_k\delta^{i_k}\tau^{j_k}\exp(-\delta^{l_k})}+\sum_{k=33}^{36}{N_k\delta^{i_k}\tau^{j_k}\exp(-\phi_k(\delta-1)^2-\beta_k(\tau-\gamma_k)^2)} $$
    
and all the terms other than $\delta$ and $\tau$ are fluid-dependent correlation parameters.

The other thermodynamic parameters can then be obtained through analytic derivatives of the Helmholtz energy terms.
For instance, the pressure is given by

$$    p=\rho RT\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau} \right] $$
    
and the specific internal energy by

$$    \frac{u}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] $$

and the specific enthalpy by


$$    \frac{h}{RT}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right] +\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+1 $$

which can also be written as

$$    \frac{h}{RT}=\frac{u}{RT}+\frac{p}{\rho RT} $$
    
The specific entropy is given by

$$    \frac{s}{R}=\tau \left[\left(\frac{\partial \alpha^0}{\partial \tau}\right)_{\delta}+ \left(\frac{\partial \alpha^r}{\partial \tau}\right)_{\delta} \right]-\alpha^0-\alpha^r $$
    
and the specific heats at constant volume and constant pressure respectively are given by

$$    \frac{c_v}{R}=-\tau^2 \left[\left(\frac{\partial^2 \alpha^0}{\partial \tau^2}\right)_{\delta}+ \left(\frac{\partial^2 \alpha^r}{\partial \tau^2}\right)_{\delta} \right] $$
    
$$ \frac{c_p}{R}=\frac{c_v}{R}+\dfrac{\left[1+\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}-\delta\tau\left(\frac{\partial^2 \alpha^r}{\partial \delta\partial\tau}\right)\right]^2}{\left[1+2\delta\left(\frac{\partial \alpha^r}{\partial \delta}\right)_{\tau}+\delta^2\left(\frac{\partial^2 \alpha^r}{\partial \delta^2}\right)_{\tau}\right]} $$
    
The EOS is set up with temperature and density as the two independent properties, but often other inputs are known, most often temperature and pressure because they can be directly measured.
As a result, if the density is desired for a known temperature and pressure, it can be obtained iteratively.

### First law for flow process {#sec:firstlaw}
The control volume sketched in [@Fig:firstlaw], separated from the surrounding by a control surface, is used as a basis for the analysis of an open thermodynamic system with flowing streams (fs) in and out, according to [@sva]

A general mass balance or continuity equation can be written: 

![Control volume with one entrance and one exit. The image has been sourced from [@firstlaw].](docs/img/First_law_open_system.png){#fig:firstlaw}

$$ \frac{m_{cv}}{dt} + \Delta \left( \dot{m} \right)_{fs}= 0 $$ {#eq:continuity}

The first term is the accumulation term i.e. the rate of change of the mass inside the control volume, $m_cv$, and the $\Delta$ in the second term represents the difference between the outflow and the inflow

$$ \Delta \left( \dot{m} \right) _{fs} = \dot{m}_2 - \dot{m}_1 $$

An energy balance for the control volume, with the first law of thermodynamics applied, needs to account for all the energy modes with can cross the control surface. Each stream has a total energy 

$$ U + \frac{1}{2}u^2 + zg $$

where the first term is the specific internal energy, the second term is the kinetic energy and the last term is the potential energy.
The rate at which each stream transports energy in or out of the control volume is given by 

$$ \dot{m} (U + \frac{1}{2}u^2 + zg) $$

and in total 

$$  \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs}$$

Furthermore, work (not to be confused with shaft work) is also associated with each stream in order to move the stream into or out from the control volume (one can think of a hypothetical piston pushing the fluid at constant pressure), and the work is equal to $PV$ on the basis of the specific fluid volume.
The work rate for each stream is 

$$ \dot{m}(PV) $$

and in total 

$$ \Delta\left[ \dot{m}(PV) \right]_{fs} $$

Further, heat may be transferred to (or from) the control volume at a rate $\dot{Q}$ and shaft work may be applied, $\dot{W}_{shaft}$.
Combining all this with the accumulation term given by the change in total internal energy the following general energy balance can be written: 

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} (U + \frac{1}{2}u^2 + zg) \right]_{fs} + \Delta \left[ \dot{m}(PV) \right]_{fs} = \dot{Q} +\dot{W}_{shaft}   $$

Applying the relation $H = U + PV$, setting $\dot{W}_{shaft} = 0$ since no shaft work is applied to or extracted from the vessel, and assuming that kinetic and potential energy changes can be omitted the energy balance simplifies to:

$$ \frac{d(mU)_{cv}}{dt} + \Delta \left[ \dot{m} H \right]_{fs} = \dot{Q}  $$

The equation can be further simplified if only a single port acting as either inlet or outlet is present:

$$ \frac{d(mU)_{cv}}{dt} + \dot{m} H  = \dot{Q}  $$ {#eq:energybalance}

where the sign of $\dot{m}$ determines if the control volume is either emptied or filled.
The continuity equation [@Eq:continuity] and the energy balance [@Eq:energybalance] combined with the equation of state are the key equations that shall be solved/integrated in order to calculate the change in temperature and pressure as a function of time.

## Flow devices {#sec:flow}
### Restriction Orifice
When a fluid flows through a restriction or opening such as an orifice, the velocity will be affected by conditions upstream and downstream.
If the upstream pressure is high enough, relative to the downstream pressure, the velocity will reach the speed of sound (Ma = 1) and the flow rate obtained will be the critical flow rate.
This condition is referred to as choked flow.
The maximum downstream pressure for the flow to still be sonic (Ma = 1), is when $P_d = P_c$.
The ratio of the critical and upstream pressure is defined by equation [@Eq:P_critical].

$$ \frac{P_{c}}{P_u}=\left (\frac{2}{k+1} \right)^\frac{k}{k-1} $$ {#eq:P_critical}

- $P_c$ is the critical pressure. [kPa]
- $P_d$ is the downstream pressure. [kPa]
- $P_u$ is the upstream pressure. [kPa]
- k is the isentropic coefficient, approximated by the ideal gas heat capacity ratio $C_p/C_v$. [-]

In order to calculate the mass flow rate through an orifice equations are used based on literature from the Committee for the Prevention of Disasters [@yellowbook]. 

To account for the difference in choked and non-choked flow a set limit pressure is introduced as in equation [@Eq:plimit].
If the downstream pressure, $P_{down}$, is below the pressure limit, $ P_{limit}$, then the flow is choked, and the pressure used, $P_{used}$, in equation [@Eq:massfloworifice] should be the pressure limit, $P_{limit}$.
Otherwise if the downstream pressure, $P_{down}$, is greater than or equal to the pressure limit, $P_{limit}$, the flow is no longer choked and the pressure used should be the downstream pressure, $P_{down}$ [@yellowbook].

$$ P_{limit}=P_{up} \cdot \left ( \frac{2}{k+1} \right ) ^{\frac{k}{k-1}} $$ {#eq:plimit}

$$ \dot{m}_{flow}= C_d  \cdot A \cdot\sqrt{\left ( \frac{2 k}{k-1}\right )  \cdot  P_{up} \cdot \rho \cdot \left (  \frac{P_{used}}{P_{up}} \right )^{\frac{2}{k}} \left (1-\left ( \frac{P_{used}}{P_{up}} \right ) ^{\frac{k-1}{k}} \right )} $$ {#eq:massfloworifice}

- $\rho$ is the density of the gas upstream. $[kg/m^3]$
- $P_{limit}$ is the pressure limit of the upstream absolute pressure. $[bara]$
- $P_{up}$ is the absolute pressure upstream of the orifice. $[bara]$
- $k$ is the ratio of the heat capacities at constant pressure, $C_p$, and at constant volume, $C_v$.
- $P_{down}$ is the absolute pressure downstream of the orifice. $[bara]$
- $P_{used}$ is the pressure used in the mass flow equation based on choked or non-choked conditions. $[bara]$
- $\dot{m}_{flow}$ is the mass flow through the orifice. $[kg/s]$
- $C_d$ is the discharge coefficient of the orifice opening. $[-]$
- $A$ is the cross sectional area of the orifice. $[m^2]$

### Pressure safety valve / Relief valve
A PSV / relief valve is a mechanical device actuated by the static pressure in the vessel and a conventional PSV is often used for gas/vapor systems.
A conventional PSV is a spring-loaded device which will activate at a predetermined opening pressure, and relieve the vessel pressure until a given reseat pressure has been reached.
Both the opening pressure and the reset pressure is above the vessel operating pressure and the PSV will remain closed until the pressure inside the vessel increases to the opening pressure.
The operation of a conventional spring-loaded PSV is based on a force balance.
A conventional PSV can be seen in [@Fig:psv].
A spring exerts a force on a disc blocking the inlet of the PSV.
When the pressure inside the vessels reaches the opening pressure, the force exerted on the disc by the gas will be larger than the force exerted by the spring and the PSV will open and allow the gas to flow out of the vessel.
The flow of gas out of the vessels will lower the pressure and thereby also the force exerted on the disc.
When the pressure in the vessels is reduced to the reset pressure, the PSV will close and the disc will again hinder the gas flow.

![Conventional/pop action PSV adapted from [@iskov] and [@API520]](docs/img/PSV.pdf){#fig:psv}

The relief valve model implemented in HydDown is the API 520 equations [@API520] for gas relief for both sonic/critical as well as subcritical flow. No corrections factors are implemented in HydDown.

For sonic flow (critical flow), as indicated in equation, the  mass flow through the PSV can be determined by equation [@eq:Stationary_sizing_sonic].

$$ W = \frac{A {C \cdot K_d \cdot  K_b \cdot  K_c \cdot  P_1}}{\sqrt{\frac{T\cdot Z}{M}}} $$ {#eq:Stationary_sizing_sonic}

- A is the effective discharge area. [mm$^2$]
- W is the mass flow through the device. [kg/h]
- C is a coefficient, as a function of k, as defined in equation [@eq:Stationary_sizing_C_function].
- $K_d$, $K_b$, and $K_c$ are correction factors.
- $P_1$ is the allowable upstream absolute pressure. [kPa]
- T is the temperature of the inlet gas at relieving conditions. [K]
- M is the molecular mass of the gas at relieving conditions. [kg/kmol]
- Z is the compressibility factor for the gas.

$K_d$ is the effective coefficient of discharge, with a typical value of 0.975, for an installed PSV.
$K_b$ is a back-pressure correction factor between 0 and 1, assumed to be 1.
$K_c$ is a correction factor used when a rupture disk is installed upstream, otherwise it is 1.
In the present implementation a value of 1 is assumed. 

$$ C=0.03948 \sqrt{k \left (\frac{2}{k+1}\right)^{\left (\frac{k+1}{k-1}\right)}} $$ {#eq:Stationary_sizing_C_function}

For subsonic flow (subcritical flow), the effective discharge area of the PSV is determined by equation [@eq:Stationary_sizing_subsonic].

$$ W = \frac{A \cdot {F_2 \cdot K_d  \cdot  K_c }}{17.9 \sqrt{\frac{T\cdot Z}{M\cdot P_1\cdot (P_1-P_2)}}} $$ {#eq:Stationary_sizing_subsonic}

F$_2$ is the coefficient of subcritical flow which can be determined from [@eq:Stationary_sizing_F2].

$$ F_2=\sqrt{\left ( \frac{k}{k-1} \right ) r^{\left (\frac{2}{k} \right )} \left (  \frac{1-r^{\left ( \frac{k-1}{k} \right )}}{1-r}\right )} $$ {#eq:Stationary_sizing_F2}

where $r$ is the ratio of backpressure to upstream relieving pressure, $P_2 / P_1$.

When modelling a pop action PSV/relief valve under dynamic conditions, the valve will go from closed to fully open in a short period of time when the set pressure, $P_{set}$, is reached.
The pop action is illustrated in [@Fig:psv_hyst] which shows the opening and closing hysteresis of the PSV as a function of pressure.
In order to close the pressure shall be reduced below the reseat pressure. 

![Relief valve hysteresis adapted from [@iskov]](docs/img/hysteresis.pdf){#fig:psv_hyst}

When specifying PSVs it common to use standard API sizes as shown in [@tbl:psv_sizes]

Size    | Area [in$^2$]     |   Area [m$^2$] 
--------|-------------------|----------------------------
D       |   0.110           |   7.09676 $\cdot$ 10$^{-5}$
E       |   0.196           |   1.26451 $\cdot$ 10$^{-4}$
F       |   0.307           |   1.98064 $\cdot$ 10$^{-4}$
G       |   0.503           |   3.24515 $\cdot$ 10$^{-4}$
H       |   0.785           |   5.06450 $\cdot$ 10$^{-4}$
J       |   1.287           |   8.30320 $\cdot$ 10$^{-4}$
K       |   1.838           |   1.18580 $\cdot$ 10$^{-3}$
L       |   2.853           |   1.84064 $\cdot$ 10$^{-3}$
M       |   3.600           |   2.32257 $\cdot$ 10$^{-3}$
N       |   4.340           |   2.79999 $\cdot$ 10$^{-3}$
P       |   6.380           |   4.11612 $\cdot$ 10$^{-3}$
Q       |   11.050          |   7.12901 $\cdot$ 10$^{-3}$
R       |   16.000          |   1.03225 $\cdot$ 10$^{-2}$
T       |   26.000          |   1.67741 $\cdot$ 10$^{-2}$

: Standard PSV orifice sizes according to API {#tbl:psv_sizes}

### Control Valve
For calculating the mass flow through a control valve, the ANSI/ISA [@borden][@ISA] methodology also described in IEC 60534 [@IEC60534] is applied. 

The flow model for a compressible fluid in the turbulent regime is: 

$$ W = C N_6 F_P Y \sqrt{x_{sizing} p_1 \rho_1} $$ 

or equivalently:

$$ W = C N_8 F_P Y \sqrt{\frac{x_{sizing}M}{T_1 Z_1}} $$ 

- C is the flow coefficient ($C_v$ or $K_v$)
- $N_8$ is a unit specific constant, 94.8 for $C_v$ and bar as pressure unit
- $F_P$ is a piping geometry factor [-]
- Y is the expansion factor [-]
- $x_{sizing}$ is is the pressure drop used for sizing[-]
- $p_1$ is the upstream pressure [bar]
- $\rho_1$ is the upstream density [kg/m$^3$]
- M is the molecular weight [kg/kmol]
- $T_1$ is the upstream temperature [K]
- $Z_1$ is the upstream compressibility [-]

In HydDown, the piping geometry factor is not yet implemented and assumed to be 1.
The pressure drop ratio $x_{sizing}$ used for sizing is determined as the lesser of the actual pressure drop ratio, $x$, and the choked pressure drop ratio $x_{choked}$.
The actual pressure drop ratio is given by:

$$ x \frac{\Delta p}{p_1}$$

The pressure drop ratio at which flow no longer increases with increased value in pressure drop ratio is the choked pressure drop ratio, given by the following equation:

$$ x_{choked} = F_\gamma x_{TP} $$

The factor $x_T$ is based on air near atmospheric pressure as the flowing fluid with a specific heat ratio of 1.40.
If the specific heat ratio for the flowing fluid is not 1.40, the factor $F_\gamma$ is used to adjust $x_T$.
Use the following equation to calculate the specific heat ratio factor:

$$ F_\gamma = \frac{\gamma}{1.4} $$

where $\gamma$ is the ideal gas $C_p/C_v$.
It should be noted that the above equation has been derived from perfect gas behaviour and extension of an orifice model with $\gamma$ in the range of 1.08 to 1.65.
If used outside the assumptions flow calculations may become inaccurate. 

The expansion factor Y accounts for the change in density as the fluid passes from the valve inlet to the vena contracta.
It also accounts for the change in the vena contracta area as the pressure differential is varied.

$$ Y = 1 -  \frac{x_{sizing}}{3x_{choked}}$$


## Heat transfer {#sec:heat}

### Natural convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to depressurization can be well approximated by that of natural convection as found from measured Nusselt numbers being well correlated with Rayleigh number, with no apparent improvement in model performance by inclusion of the Reynolds number in the model [@woodfield].

To determine the heat transfer for the gas-wall interface, the following is applied cf. equation [@Eq:newton]:

$$  \frac{dQ}{dt} = \dot{Q} = h A ( T_{s} - T_{gas} ) $$  {#eq:newton}

- $d Q$ is the change in thermal energy due to convective heat transfer. [J]
- $d t$ is the change in time during the heat transfer. [s]
- $h$ is the convective heat transfer. [W/m$^2$ $\cdot$ K ]
- $A$ is the area normal to the direction of the heat transfer. [m$^2$]
- $T_{s}$ is the inner surface temperature of the geometry. [K]
- $T_{gas}$ is the temperature of the bulk gas inside the vessel. [K]

The convective heat transfer will need to be estimated for the the gas-wall interface, by the use of empirical relations for the Nusselt number.
The Nusselt number describes the ratio of convective heat transfer to conductive heat transfer, normal to a surface area, as given in equation [@eq:Nu].

$$ Nu=\frac{hL}{k} $$ {#eq:Nu}

- $Nu$ is the Nusselt number. [-]
- $h$ is the convective heat transfer. [W/m$^2$$\cdot$K]
- $L$ is a characteristic length of the geometry. [m]
- $k$ is the thermal conductivity of the gas. [W/m$\cdot$K]

The characteristic length $L$ used is the height of the gas volume i.e. the length of a vertical vessel or the diameter of a horizontal vessel.

The empirical correlations used to calculate the Nusselt number of the gas-wall interface is a function of the Rayleigh number, which can be defined by the Grashof number and Prandtl number, as in equation [@Eq:rayleigh_gas]:

$$ Ra=Gr \cdot Pr $$ {#eq:rayleigh_gas}

- $Ra$ is the Rayleigh number. [-]
- $Gr$ is the Grashof  number. [-]
- $Pr$ is the Prandtl number. [-]

The Grashof number is a dimensionless number which approximates the ratio of the buoyancy forces to viscous forces, as given in equation [@Eq:grashof_gas}:

$$ Gr=\frac{\beta g\rho^2 L^3 \Delta T }{\mu^2} $$ {#eq:grashof_gas}

The Prandtl number is a dimensionless number defined as the ratio of the momentum diffusivity to thermal diffusivity, as given in equation [@Eq:prandtl_gas]:

$$ Pr=\frac{c_p \mu}{k} $$ {#eq:prandtl_gas}

- $\beta$ is the coefficient of volume expansion. [1/K]
- $g$ is the standard acceleration of gravity. [m/s$^2$]
- $\rho$ is the gas density. [kg/m$^3$]
- $L$ is the characteristic length. [m]
- $\Delta T$ is the temperature difference of the surface and gas. [K] 
- $\mu$ is the dynamic viscosity. [kg/m$\cdot$s]
- $c_p$ is the heat capacity of gas. [J/kg$\cdot$K]
- $k$ is the thermal conductivity of gas. [J/m$\cdot$K]

It is important to note that the properties in the above equations shall be evaluated at the fluid film temperature which can be approximated by the average of the the fluid bulk temperature and the vessel wall temperature [@geankoplis].

The mechanism of heat transfer on the outside of the vessel at ambient conditions is similar to the above.
In HydDown the external heat transfer is not modelled currently.
A heat transfer coefficient shall be provided.

### Mixed convection
Experiments have indicated that the internal heat transfer mechanism for a vessel subject to filling can be well approximated by that of combined natural convection and forced convection as found from measured Nusselt numbers being well correlated with Rayleigh and Reynolds number [@woodfield].

For mixed convection the effective Nusselt number, $Nu$, can be approximated by 

$$ Nu = (Nu_{forced}^n + Nu_{natural}^n)^{\frac{1}{n}} $$

During charging with different gases (H$_2$, N$_2$ and Argon), Woodfield *et al.* demonstrated that in order to provide a good fit to the experimentally determined Nusselt number a correlation based on both Reynolds and Rayleigh number was necessary.
The following formula was found to provide a good fit (n=1):

$$ Nu = Nu_{forced} +  Nu_{natural} = 0.56Re_d^{0.67} + 0.104Ra_H^{0.352} $$

### Conduction
For accurate prediction of the outer and especially the inner wall temperature for correct estimation of internal convective heat transfer and the average material temperature, the general equation of 1-D unsteady heat transfer shall be solved:

$$ \frac{\delta T}{\delta t} = \frac{k}{C_p} \frac{\delta^2 T}{\delta x^2} $$

- T is temperature
- x is the spatial (1-D) coordinate
- k is the thermal conductivity 
- $C_p$ is the heat capacity 

Here it is written in Cartesian coordinates, but for most applications to pressure equipment, cylindrical coordinates are applicable, at least for the shell.
To be solved, the initial values and boundary values must be specified.
In its present state, HydDown does not include the unsteady heat transfer model, i.e. the assumption is that the temperature from outer to inner surface is uniform and equal to the average temperature.
This is obviously a crude approximation, but might be justified depending in the Biot number:

$$ Bi = \frac{hL}{k} $$

The Biot number is a simple measure of the ratio of the heat transfer resistances at the surface of a body to the inside of a body.
The ratio gives an indication to which extent the temperature will vary in space (gradient) when the body is subject to a displacement in temperature at the surface boundary layer.
Striednig *et al.* [@STRIEDNIG] concluded that for a type I (steel) cylinder the Biot number was approx. 0.03 and hence the error in assuming a uniform temperature in the vessel wall was low. 

With a typical thermal conductivity of 45 $W/m K$ for steel and a heat transfer coefficient up to 600 $W/m^2 K$ [@woodfield] the Biot number for a vessel with a wall thickness of 2 cm is 0.27.
This is significantly higher than that approximated by [@STRIEDNIG].
However, the Biot number is lower than 1, and the assumption of a uniform temperature is reasonable.
However, for increased wall thickness, and/or for different materials with lower thermal conductivity, the error may grow to an unacceptable level. 

### Fire heat loads
The heat transfer from the flame to the shell is modelled using the recommended approach from Scandpower [@scandpower] and API [@API521].
The heat transfer from the flame to the vessel shell is divided into radiation, convection, and reradiation as seen in equation [@Eq:flame]:

$$ q_f=\underbrace{{\alpha}_s \cdot {\varepsilon}_f \cdot \sigma \cdot T_f^4}_\text{Radiation}+\underbrace{h_f \cdot (T_f-T_s(t))}_\text{Convection}-\underbrace{{\varepsilon}_s \cdot \sigma \cdot T_s(t)^4 }_\text{Reradiation} $$ {#eq:flame}

- $q_f$ is the flame heat flux. [W/m$^2$]
- ${\alpha}_s$ is the vessel surface absorptivity. [-]
- ${\varepsilon}_f$ is the flame emissivity. [-]
- $\sigma$ is the Stefan-Boltzmann constant, $\sigma$ = $5.67 \cdot 10 ^ {-8}$  [W/m$^2 \cdot$ K$^4$]
- $T_f$ is the flame temperature. [K]
- $h_f$ is the convection heat transfer coefficient between the flame and the surface. [W/m$^2 \cdot$ K]
- $T_s(t)$ is the time dependent surface temperature. [K]
- ${\varepsilon}_s$ is the surface emissivity. [-]

This model assumes that the pressure vessel is fully engulfed by the flame.
This means that the view factor for the radiation is unity and is therefore not taken into consideration.
The convective heat transfer coefficients for a jet fire and a pool fire, and recommended values for the emissivity and absorptivity, are given by Scandpower as [@scandpower]:

- $h_{jet~fire}$ = 100 [W/m$^2 \cdot$K]
- $h_{pool~fire}$ = 30 [W/m$^2 \cdot$K]
- ${\alpha}_s$ = 0.85
- ${\varepsilon}_s$ = 0.85
- ${\varepsilon}_f$ = 1.0 (optical thick flames, thickness > 1 m)

The flame temperature is found by solving equation [@Eq:flame2] for the incident heat flux in relation to the ambient conditions.
The flame temperature is kept constant throughout the simulation:

$$ q_{total}=\sigma \cdot T_f^4 + h_f \cdot (T_f-T_{amb})$$ {#eq:flame2}

- $q_{total}$ is the incident flame heat flux as given in table [@Tbl:heatfluxes1]. [W/m$^2$]
- $T_{amb}$ is the ambient temperature $\approx$ 293 K (20$^\circ$ C)

The heat flux used to calculate the flame temperature is given in table [@tbl:heatfluxes1].

|                        | Small jet fire  [kW/m$^2$]  |  Large jet fire  [kW/m$^2$] |  Pool fire  [kW/m$^2$]
| ----                   |  ----           |  ----           | ----
| Peak heat load         |  250            |  350            | 150         
| Background heat load   |   0             |   100           |  100         

: Incident heat fluxes for various fire scenarios given by Scandpower [@scandpower] {#tbl:heatfluxes1}


## Model implementation
A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation as described in [@Sec:flow]:

$$ m_{cv}(i+1) =  m_{cv}(i) + \dot{m}(i) \Delta t  $$ {#eq:euler_mass}

$$ \dot{m}(i) = f(P,T,) $$

For each step, the mass relief/ left in the vessel is known.
Since the volume is fixed the mass density is directly given. 

For the calculation methods (isentropic, isenthalpic, isenergetic, etc.), CoolProp allows specifying density and either H, S, or U directly - this is very handy and normally only TP, PH, or TS property pairs are implemented and you would need to code a second loop to make it into am UV, VH, or SV calculation.

In the following, the high-level steps in the solution procedure are outlined for the different calculation types.
To illustrate calls to the equation of state in CoolProp, the notation $EOS(D,T)$ corresponds to calculation of the state at known density and temperature, for example.

### Isothermal process
For an isothermal process, the solution procedure for each calculation step is the following with the mass in the control volume being calculated from [@Eq:euler_mass]:

$$ T(i+1) = T(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ P(i+1) = EOS(D(i+1),T(i+1)) $$

### Isenthalpic process
For an isenthalpic process: 

$$ H(i+1) = H(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ P(i+1) = EOS(D(i+1),H(i+1)) $$
$$ T(i+1) = EOS(D(i+1),H(i+1)) $$

### Isentropic process
For an isentropic process:

$$ S(i+1) = S(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V}$$
$$ P(i+1) = EOS(D(i+1),S(i+1)) $$
$$ T(i+1) = EOS(D(i+1),S(i+1)) $$

### Isenergetic process
For an isenergetic process: 

$$ U(i+1) = U(i) $$
$$ D(i+1) = \frac{m_{cv}(i+1)}{V}$$
$$ P(i+1) = EOS(D(i+1),U(i+1)) $$
$$ T(i+1) = EOS(D(i+1),U(i+1)) $$

### Energy balance
The general first law applied to a flow process as outlined in [@Sec:firstlaw] subject to an explicit Euler scheme is: 

$$ D(i+1) = \frac{m_{cv}(i+1)}{V} $$
$$ U_{cv}(i+1) = \frac{m_{cv}(i)U_{cv}(i) - \left( \dot{m}(i) H (i) +  \dot{Q}(i) \right) \Delta t}{m_{cv}(i+1)}  $$ {#eq:firstlaw_euler}

The above assumes that mass flow is positive when leaving the control volume and heat rate is positive when transferred to the control volume.
$H(i)$ is the specific enthalpy of the fluid in the control volume for a discharging process and it is equal to the energy of the entering stream for a filling process.
The heat rate is calculated as outlined in [#Sec:heat]. 

For the vessel wall a simple heat balance is also made:

$$ T_{wall}(i+1) = T_{wall}(i)  + \frac{\dot{Q}_{outer} - \dot{Q}_{inner} } {c_p m_{vessel}} \Delta t $$

where $\dot{Q}_{outer}$ is the convective heat transfer to or from the ambient surroundings from the outer surface of the vessel, with positive values indicating that heat is transferred from the surroundings to the vessel wall.
This is either a fixed heat transfer coefficient (for now) with a specified ambient temperature or a calculated fire heat load.
$\dot{Q}_{inner}$ is the internal heat transfer, either a fixed number or calculated as natural convection (for discharge) or mixed natural and forced convection (filling).

If a fixed $\dot{Q}$ is provided or a fixed overall heat transfer coefficient is provided the vessel wall temperature heat balance is not solved. 

The remaining steps are update of temperature and pressure:

$$ P(i+1) = EOS(D(i+1),U(i+1)) $$
$$ T(i+1) = EOS(D(i+1),U(i+1)) $$
 
## Multicomponent mixtures
Although not an initial requirement, the code can handle multi-component gas mixtures.
However, no validation has been performed.
Furthermore, when calculations are done for multicomponent mixtures, the code runs significantly slower.
Please note, that in case that liquid condensate is formed during discharge calculations and even if the calculations does not stop, the results cannot be considered reliable.
This is because component balances are not made and it is always assumed that the discharge composition is the same as the global composition inside the vessel.
When liquid condensate is formed the composition of the vapour phase will differ from the global composition. 

There are a few examples of multicomponent mixtures included with HydDown.
In order to specify multicomponent mixtures the below example can be used as guidance: 

    "Methane[9.01290129e-01]&Ethane[6.35063506e-02]&N2[7.80078008e-03]&CO2[2.34023402e-02]&Propane[3.50035004e-03]&Butane[5.00050005e-04]"

For component names please refer to the [Coolprop manual](http://www.coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids).

![Gas discharge of multicomponent mixture](docs/img/ng_all.png){#fig:ng_all}

An example calculation for the above mixture is shown in [@Fig:ng_all]. 

![Multicomponent mixture phase envelope and pressure and temperature trajectory of the vessel inventory during discharge](docs/img/ng_pe.png){#fig:ng_pe}

The pressure and temperature trajectory is visualised along with the fluid mixture phase envelope in [@Fig:ng_pe].
As seen from the figure, this case is borderline and the pressure/temperature trajectory just coincides with the dew line on the phase envelope.
This plot is included in the example main script that comes with HydDown and serves as an important quality control. 

# Validation
The code is provided as-is.
However, comparisons have been made to a few experiments from the literature.

The following gases and modes are considered:

- High pressure hydrogen discharge
- High pressure hydrogen filling
- High pressure nitrogen discharge

## Hydrogen discharge
Calculations with HydDown for high pressure hydrogen vessel discharge have been compared to three experiments from [@byrnes].
Byrnes *et al.* [@byrnes] have reported a number of rapid depressurization experiments using both nitrogen and hydrogen, with blowdown times being varied between 14 seconds to 33 minutes.
The vessel used was a type "K" gas cylinder with a volume of 0.0499 m$^3$, an ID of 8.56 inches and shell length of 55 inches.
The vessel is placed in a vertical position.
A control valve is mounted on the neck of the cylinder with a fixed opening percentage during each experiment, set to provide the required depressurization rate.
It is mentioned that the gas cylinder is placed inside a weather proof cabinet with perlite insulation in order to emulate adiabatic conditions.
However the simulations run here are with a finite heat transfer coefficient on the outside of the cylinder.
In their paper, detailed pressure and temperature traces are provided for three experimental runs: 7, 8, and 9, with depressurization from 2000 psia in 30 seconds, 480 seconds and 14 seconds, respectively.
These experiments are simulated in [@Fig:byrnes_run7], [@Fig:byrnes_run8] and [@Fig:byrnes_run9]. 
  
![Calculations of hydrogen discharge experiment run 7 from [@byrnes]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](docs/img/Byrnes_run7.png){#fig:byrnes_run7}

The figure layout is general for all the validation studies conducted and will be described in brief for eased readability.
Each experiment figure is divided in 4 subplots: the upper left shows the simulated bulk gas temperature and vessel wall (average) temperature with measured values included with stipulated lines (if available), the upper right shows the simulated and measured pressure (if available, shown with points), the lower left shows the specific enthalpy, specific internal energy, and the specific entropy of the bulk gas content as a function of time, the lower right shows the simulated mass rate leaving or entering the vessel.

As seen from [@Fig:byrnes_run7], [@Fig:byrnes_run8], and [@Fig:byrnes_run9], the calculations generally compare well with the experimental results.
It appears as though the faster depressurization times are predicted the best, as seen from run 7 and 9.
The much slower run 8 does not have as good a fit for the bulk gas temperature.
However, on an absolute scale the results are still within 5$^\circ$C of the experimentally determined gas temperature.
The trend in *goodness of fit* could indicate that the heat transfer from the outside of the vessel to the gas lacks some details, which become detectable at the very slow depressurization.
For the fast depressurization the heat transfer from the outside is of less importance and the inside heat transfer coefficient may be predicted better.

![Simulation of run 8 from [@byrnes].](docs/img/Byrnes_run8.png){#fig:byrnes_run8}

![Simulation of run 9 from [@byrnes].](docs/img/Byrnes_run9.png){#fig:byrnes_run9}

## Hydrogen filling 
In order to compare HydDown to experimental values of hydrogen filling operation the experiments reported by Striednig *et al.* [@STRIEDNIG] are used.
Experimental results were obtained using a type I tank (steel) with a volume of 0.0235 m$^3$.
The filling of hydrogen of gas into the vessel was carried out with an upstream reservoir kept at 350 bar and the flow was controlled with an electronically controlled dispenser.
Vessel details are provided below: 

| Parameter       | Value       |
|-----------------|-------------|
| Mass            | 69..5 kg    |
| Length          | 0.61 m      |
| OD              | 0.280 m     |
| Wall thickness  | 0.0129 m    |
| Density         | 7740 kg/m$^3$ |
| Heat capacity   | 470 J/(kg K) |

: Key hydrogen vessel data [@STRIEDNIG] {#tbl:h2filling}

[@STRIEDNIG] reports many different experiments.
In this comparison we will use experiments applying different pressurisation rates / durations using identical initial conditions.
The experiments are reported in Fig. 6 and Fig. 7 in [@STRIEDNIG].

In its current stage HydDown does not directly support specifying a fixed pressurisation rate expressed as a constant rate of change in pressure / pressure ramp rate (PRR).
In order to emulate the experiments, an orifice model is used with an artificially high reservoir pressure in order to ensure that the flow is critical during the entire pressurisation.
This gives a constant mass rate and translates to a fairly constant pressure ramp rate, although with real gas effects being clearly observable.

The results of the simulations with HydDown and the experiments from [@STRIEDNIG] are shown in [@Fig:striednig5], [@Fig:striednig10], [@Fig:striednig30].
As seen from the figures the experimental gas temperature appears to be simulated very accurately.
The simulation with the fastest pressurisation under-predicts the experimental temperature slightly.
The main features of the pressurisation temperature profile is a rapid increase in temperature, at the beginning of the pressurisation, where the effects of filling dominate.
At approximately 20-30 s, the temperature increase levels off to a plateau, where the length of the plateau depends on the pressurisation rate.
The slower the depressurization, the longer the plateau.
During this temperature plateau, the heat transfer from the gas to the vessel and surroundings balance the temperature increase that would occur if performed at adiabatic (isolated) conditions.
When the gas flow has ceased upon reaching the final pressure, the temperature starts to drop due to cooling of the gas by the colder vessel wall driven by convective heat transfer.
The temperature stabilises at a temperature which is higher than the initial ambient temperature due to the heating of the steel vessel.
The cooling of the vessel by ambient air is slower and would require a longer run time to be clearly visible. 

![Simulation of $H_2$ pressurization using 5 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_5MPa_min.png){#fig:striednig5}

![Simulation of $H_2$ pressurization using 10 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_10MPa_min.png){#fig:striednig10}

![Simulation of $H_2$ pressurization using 30 MPa/min [@STRIEDNIG].](docs/img/Striednig_fillingH2_30MPa_min.png){#fig:striednig30}

## Nitrogen discharge
Calculations with HydDown are compared to experiment I1 from ref. [@Haque1992b].
The experiment is a blowdown of a vertically oriented cylindrical vessel with flat ends.
The vessel length is 1.524 m, the inside diameter is 0.273 m, and the wall thickness is 25 mm.
The vessel is filled with N$_2$ at 150 bar and  15$^\circ$C.
Ambient temperature is 15$^\circ$C.
The blowdown orifice diameter is 6.35 mm.
The results are shown in [@Fig:N2val].
The discharge coefficient of the orifice has been set to 0.8 in order to match the vessel pressure profile.
The back pressure is set to atmospheric conditions.

![Calculations of nitrogen discharge emulating experiment I1 from [@Haque1992b]. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right).](docs/img/N2_filling.png){#fig:N2val}

As seen from [@Fig:N2val], the calculations compare well with the experimental results.
The calculated temperature of the bulk vapor is within the experimental range of measured temperature at all times during the simulation.
It is also noted that the minimum temperature is reached at approximately the same time as in the experiments.
The calculated vessel inner wall temperature does not decline as rapidly as the experiments, but from around a calculation time of 60 s, the temperature is within the experimentally observed inner wall temperature.
The main reason for the inability to match the vessel wall temperature is that the model ignores the temperature gradient from the outer to the inner wall surface and uses an average material temperature.
Especially at the beginning of the discharge it is considered likely that a significant temperature gradient will exist. 
[![DOI](https://zenodo.org/badge/353152239.svg)](https://zenodo.org/badge/latestdoi/353152239) ![license](https://img.shields.io/github/license/andr1976/HydDown) ![buil](https://github.com/andr1976/HydDown/actions/workflows/python-app.yml/badge.svg) [![codecov](https://codecov.io/gh/andr1976/HydDown/branch/main/graph/badge.svg)](https://codecov.io/gh/andr1976/HydDown) [![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py)
[![CodeQL](https://github.com/andr1976/HydDown/actions/workflows/codeql-analysis.yml/badge.svg)](https://github.com/andr1976/HydDown/actions/workflows/codeql-analysis.yml)
[![status](https://joss.theoj.org/papers/0eed2a25a99589ed8dcdc785c890fb25/status.svg)](https://joss.theoj.org/papers/0eed2a25a99589ed8dcdc785c890fb25)
 
# HydDown
Hydrogen (or other pure gas phase species) depressurization calculations 

This code is published under an MIT license.

Install as simple as:

    pip install hyddown
    
In the case of an error in installing relation to python2, for example:
```
DEPRECATION: Python 2.7 reached the end of its life on January 1st, 2020. Please upgrade your Python as Python 2.7 is no longer maintained. pip 21.0 will drop support for Python 2.7 in January 2021. More details about Python 2 support in pip can be found at https://pip.pypa.io/en/latest/development/release-process/#python-2-support pip 21.0 will remove support for this functionality.
Defaulting to user installation because normal site-packages is not writeable
ERROR: Could not find a version that satisfies the requirement hyddown (from versions: none)
ERROR: No matching distribution found for hyddown
```
please instead install with

    python3 -m pip install hyddown

or try
    
    pip3 install hyddown


Run the code as simple as: 

    python main.py input.yml

where main.py is the main script and input.yml is the input file in Yaml syntax. 

Consult the [manual](https://github.com/andr1976/HydDown/raw/main/docs/MANUAL.pdf) for a more rigorous explanation of the software, the implemented methods, and its usage. Further, the manual also contains a few validation studies. 

## IMPORTANT NOTICE:
This package, or CoolProp to be more precise, runs on python 3.8 and apparently the CoolProp `pip` package has issues with python 3.9. Maybe compiling from source works, but this has not been tested. Hence, it is highly recommended that python 3.8 is used. 

## Demonstration 
The easiest way to explore the capability of HydDown is the [streamlit app](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py). This version allows calculation of:

- filling of vessel with gas (pressurisation)
- discharge of gas (depressurisation)
- various gases (H2, N2, CH4, He, Air)
- variable size pressure cylinder/vessel
- heat transfer between gas and vessel wall can be turned on/off

## Background
This is a small spare time project for calculation of vessel filling/depressurisation behaviour. This is mainly to demonstrate, that although perceived as a very tedious/difficult/complex problem to solve, actually a fairly limited amount of code is necessary if you have a good thermodynamic backend. 

A few choices is made to keep things simple to begin with:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Mainly pure substances are considered (mixtures can be handled - but calculations can be slow)
- Gas phase only 
- No temperture stratification in the gas phase
- No temperture gradient through vessel wall

The code is as simple as possible. The above choices makes the problem a lot more simple to solve, First of all the pure substance Helmholtz energy based equation of state (HEOS) in coolprop offers a lot of convenience in terms of the property pairs/state variables that can be set independently. Using only a single gas phase species also means that component balances is redundant and 2 or 3-phase flash calculations are not required. That being said the principle used for a single component is more or less the same, even for multicomponent mixtures with potentially more than one phase.

## Description
The following methods are implemented:

- Isothermal i.e. constant temperature of the fluid during depressurisation (for a very slow process with a large heat reservoir)
- Isenthalpic/Adiabatic (no heat transfer with surroundings, no work performed by the expanding fluid)
- Isentropic (no heat transfer with surroundings, PV work performed by the expanding fluid)
- Constant internal energy
- Energy balance. This is the most general case and includes the ability to transfer heat with surroundings

Various mass flow equations are enabled: 

- Orifice 
- Control valve 
- Relief valve (discharge only)
- Constant mass flow

A simple (naive) explicit Euler scheme is implemented to integrate the mass balance over time, with the mass rate being calculated from an orifice/valve equation. For each step, the mass relief/ left in the vessel is known. Since the volume is fixed the mass density is directly given. For the simple methods (isentropic,isenthalpic,isenergetic etc), Coolprop allows specifying density and either H,S or U directly - this is very handy and normally only TP, PH, TS property pairs are implemented, and you would need to code a second loop to make it into am UV, VH or SV calculation. Coolprop is very convenient for this, however for a cubic EOS and for multicomponent Helmholtz energy EOS coolprop only supports a subset of state variables to be specified directly (T,P,quality). For this reason single component HEOS is the main target of this small project.  In case the "Energy balance" method is applied, the heat added from convection and work is accounted for. 

## Basic usage
The Yaml input file is edited to reflect the system of interest. For isothermal/isenthalpic/isentropic/isenergetic calculations the minimal input required are:

- Initial conditions (pressure, temperature)
- vessel dimensions (ID/length)
- valve parameters (Cd, diameter, backpressure)
- Calculation setup (time step, end time)
- Type of gas

If heat transfer is to be considered the calculation type "energybalance" is required. A few options are possible:

- Fixed U (U-value required, and ambient temperature)
- Fixed Q (Q to be applied to the fluid is requried)
- Specified h, the external heat transfer coefficient is provided and either the internal is provided or calculated from assumption of natural convection from a vertical cylinder at high Gr number. Ambient temperature is required.
- Detailed 
- Fire with heat load calculated from the Stefan-Boltzmann equation
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity
and orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
  and learning from the experience
* Focusing on what is best not just for us as individuals, but for the
  overall community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
  advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
  address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards of
acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies when
an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail address,
posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at
.
All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series
of actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or
permanent ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior,  harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within
the community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0, available at
https://www.contributor-covenant.org/version/2/0/code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at
https://www.contributor-covenant.org/translations.
# Contributing to HydDowm
HydDown loves your input! We want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Support/questions
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

When contributing to this repository, please first discuss the change you wish to make via filling an issue, email, or any other method with the owners of this repository before making a change. 

Please note we have a [code of conduct](https://github.com/andr1976/HydDown/blob/4665472fc1a1baf60dcdc8485ed2a049ff7273fe/CODE_OF_CONDUCT.md), please follow it in all your interactions with the project.

## Support
For questions, assistance using the code or special considerations please file an [issues](https://github.com/andr1976/HydDown/issues) and add the label 'question' or 'help wanted'

## We Develop with Github
HydDown uses github to host code, to track issues and feature requests, as well as accept pull requests.

Pull requests are the best way to propose changes to the codebase. We actively welcome your pull requests:

1. Fork the repository from 'main' to your own Github account
2. Clone the project to your machine
3. Create a branch locally with a succinct but descriptive name
3. Commit changes to the branch
4. Push changes to your fork
5. If you've added code that should be tested, add tests.
6. If you've changed APIs added functionality, update the documentation.
8. Ensure the test suite passes (handled by github action).
9. Make sure your code lints (handled by github action).
10. Open a PR in our repository

## Any contributions you make will be under the MIT Software License
In short, when you submit code changes, your submissions are understood to be under the same [MIT License](http://choosealicense.com/licenses/mit/) that covers the project. 

## Report bugs using Github's [issues](https://github.com/andr1976/HydDown/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/andr1976/HydDown/issues/new/choose); it's that easy!
Add the label 'bug' to the filed issue - which may be removed by the maintainer if this is actually not a bug, but e.g. enhancement/feature request. 

## Write bug reports with detail, background, and sample code
**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can. 
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

## Use a Consistent Coding Style

## License
By contributing, you agree that your contributions will be licensed under its MIT License.

## References
This document was adapted from the open-source contribution guidelines for [Facebook's Draft](https://github.com/facebook/draft-js/blob/a9316a723f9e918afde44dea68b5f9f39b7d9b00/CONTRIBUTING.md) as well as from [Auth0](https://github.com/auth0/open-source-template/blob/master/GENERAL-CONTRIBUTING.md) and [briandk](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62) 
---
title: 'HydDown: A Python package for calculation of hydrogen (or other gas) pressure vessel filling and discharge'
tags:
  - Python
  - Gas storage
  - Depressurisation
  - Blow-down
  - Pressure cylinder
  - Energy storage
  - Hydrogen
authors:
  - name: Anders Andreasen
    orcid: 0000-0003-0475-323X
    affiliation: 1
affiliations:
 - name: Ramboll Energy, Field Development, Bavnehøjvej 5, DK-6700 Esbjerg, Denmark
   index: 1
date: 04 August 2021
bibliography: paper.bib

---

# Summary
HydDown [@anders_andreasen_2021_5154096] is a Python package for calculation of pressure vessel behaviour during filling (pressurisation) or discharge (depressurisation/blow-down). More specifically, the software allows calculation of vessel pressure, fluid inventory temperature as well as vessel wall temperature as a function of time during either filling or discharge operations. The applications are manifold and some examples are: 

* Rapid filling of vehicle high-pressure hydrogen cylinders
* Rapid emergency discharge from high-pressure gas cylinders
* Dynamic pressure safety valve behaviour for gas-filled pressure vessels subject to fire heat load

A typical system modelled is shown in \autoref{fig:sketch} and it visualizes the key parameters and transport phenomena during gas vessel filling or discharging. The thermodynamic state inside the vessel changes over time as seen from immediately observable variables temperature (T) and pressure (P). This is caused by a change in fluid inventory (density) due to flow of gas either in or out of the vessel. Furthermore, heat is transfered from or to the surroundings via convective heat transfer on the in- and outside of the vessel - with heat being conducted through the vessel wall. The easiest way to explore the basic capabilities fo the code is via the [HydDown Streamlit app](https://share.streamlit.io/andr1976/hyddown/main/scripts/streamlit_app.py).  

![Gas-filled pressure vessel subject to gas discharge and heat transfer between vessel and gas inventory. \label{fig:sketch}](../docs/img/Sketch.png){ width=70% }

In its essence, the code solves the mass and energy balances with gas thermodynamics calculated using the [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999] including both real gas equation of state and gas transport properties. The energy balance is the first law of thermodynamics for an open system exchanging both heat and mass with the surroundings [@sva]. Heat transfer between gas inventory and vessel wall is accounted for using either natural convection or mixed forced convection/natural convection [@woodfield;@geankoplis]. The mass balance is closed using an applicable flow equation e.g. orifice  [@yellowbook], relief valve [@API520],  control valve [@borden;@ISA;@IEC60534], a fixed mass rate or a predefined mass rate in or out of the pressure vessel. 
The code also allows an external heat load to be applied using the Stefan-Boltzmann equation for fire heat load mimicking background heat load from both pool and jet fire [@scandpower;@API521].

A few choices have been made to keep things simple:

- [Coolprop](http://www.coolprop.org/) is used as thermodynamic backend
- Gas phase only
- No temperature stratification in vessel inventory
- No temperature gradient through vessel wall (applicable for high heat conductivity / thin-walled vessels)

Still, the code allows a variety of mass flow devices and heat transfer both with ambient air and also considering fire heat loads, with the vessel fluid inventory being rigorously described by a Helmholtz energy formulation of a real gas equation of state. Typical calculation output is shown in \autoref{fig:N2discharge} and \autoref{fig:H2filling} with experimental data included for comparison. 

![Calculations of nitrogen discharge emulating experiment I1 from @Haque1992b. The figure shows calculated gas and wall temperature (full lines) compared to experiments (upper left), calculated and experimental pressure (upper right), specific thermodynamic state variables (lower left), and the calculated vent rate (lower right). \label{fig:N2discharge}](../docs/img/N2_filling.png)

![Simulation of hydrogen cylinder pressurisation using a pressurisation rate of 10 MPa/min. Comparison between calculated (full line) and measured gas temperature from @STRIEDNIG (stipulated line) is shown in the upper left graph. \label{fig:H2filling}](../docs/img/Striednig_fillingH2_10MPa_min.png)

# Statement of need
With an increasing demand of clean(er) energy and associated storage such as e.g. on-board hydrogen storage for hydrogen powered vehicles, compressed air energy storage (CAES), compressed biogas, compressed natural gas (CNG) etc., the need for tools, which can simulate the filling and emptying of the pressure containers holding the fluid is indeed present both for operational reasons, but even more importantly for safety reasons. 

Apparently, no free (as in speach, as in beer) tool exists today which accomplishes the same tasks as HydDown. It has been one of the main goals in producing this software to provide a free tool, which can calculate the pressure vessel response during filling and discharge of gas. There are, however, many tools available which can do a subset of HydDown capabilities, the same as HydDown, or even more. Most of these tools are commercial/proprietary (and closed source) and come with a significant license cost. A single tool has been identified, which is free (as in beer), but still closed source [@h2fills]. Below is a partial list of tools, which have comparable capabilities as HydDown - for review of more codes, please refer to @SHAFIQ2020104.  

| Software                      | Cost                  |  Filling/discharge    | Fire heat laod    |
|-------------------------------|-----------------------|-----------------------|-------------------|
| [Aspen HYSYS](https://www.aspentech.com/en/products/engineering/aspen-hysys)                  |  Fee           |  Both                 | Yes               |
| [Unisim Design](https://www.honeywellprocess.com/en-US/online_campaigns/connected_plant/Pages/process-simulation.html)                 |  Fee           |  Both                 | Yes               |
| VBsim [@DALESSANDRO2015719]   |  N/A             |  Discharge            | No                |
| BLOWSIM [@blowsim]            |  N/A             |  Discharge            | Yes               |
| H2FillS [@h2fills]            |  Free            |  Filling              | No                |

Compared to other software codes, HydDown has a number of advantages: It is open source, computations are fast - especially for single component gas, the code is simple to install, the code is easy to run using a simple input file format, both filling and discharge can be handled and fire heat load can be handled for calculations of emergency response. In this way, the code has advantages over existing codes on one or more parameters. However, the code is still limited to gas phase behaviour and heat transfer through the vessel wall is not rigorously modelled. Hence composite cylinders with different layers of material may give inaccurate results, especially for low heat conductivity materials. Some of the above codes can handle this.   

It is the intent that the present software can be of use for a number of tasks including (but not limited to):

* Design aid for experimental facilities for e.g. mass flow device sizing e.g. orifice, control valve etc. 
* Benchmarking of other software codes. 
* Safety device sizing and review (blow-down valve, pressure safety valve), for instance how fast shall blow-down be? How large should an installed pressure safety valve be?
* Thermal response of vessel during filling/discharge. Is the design (lower/upper) temperatures exceeded? Should the rate be lower? Is any pre-cooling/heating required? 

# Acknowledgements
The making of this project would not have been possible without the great [CoolProp](http://www.coolprop.org/) library [@doi:10.1021/ie4033999]. The author is also thankful for enlightning discussions with colleague Jacob Gram Iskov Eriksen (Ramboll Energy, Denmark) and former Ramboll Energy colleague Carsten Stegelmann (ORS Consulting) in relation to vessel depressurisation, nozzle flow and heat transfer considerations. This work has also benefitted significantly from the understanding of thermo-mechanical pressure vessel behaviour obtained in previous works [@Bjerre2017;@safety4010011;@iskov]. Furthermore, this project relies on high quality open source Python packages: NumPy [@Walt:2011:NAS:1957373.1957466;@harris2020array], SciPy [@virtanen2020scipy], matplotlib [@Hunter:2007], pandas [@mckinney-proc-scipy-2010], [PyYaml](https://pyyaml.org/wiki/PyYAMLDocumentation), [cerberus](https://docs.python-cerberus.org/en/stable/), [Streamlit](https://streamlit.io/), tqdm [@tqdmref] and pytest [@pytest]. Sussanne Tolstrup, language secretary at Ramboll Energy, Field Development, is acknowledged for proofreading this manuscript. 

# References
