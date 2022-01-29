
o80 (pronounced 'oh eighty') is a tool for synchronizing processes while organizing exchange of information between them.
The information exchanged are commands for computing (robotic) desired states and observations, where state and observation are user class (and o80 classes are templated over them).
o80 is in c++, with python wrappers.

Here is a basic usage of o80, using python:

```python

# starts a process running robot control at
# 1000Hz
frequency = 1000
start_standalone("my_robot",frequency)

# creating a frontend for communicating with the robot
frontend = FrontEnd("my_robot")

# creating a command requesting the desired state of
# the 0th actuator to reach the value of 100 in 3 seconds
frontend.add_command(0,State(100),Duration.seconds(3),QUEUE)
# same for the 1th actuator
frontend.add_command(1,State(100),Duration.seconds(3),QUEUE)

# send the commands. The desired states of 0th and 1st dof
# of the robot interpolates from their initial values to 100.
# Once this is achieved, an observation is returned
observation = frontend.pulse_and_wait()

# reading the state of the robot 
state0 = observation.get_observed_states().get(0).get()
state1 = observation.get_observed_states().get(1).get()

# stopping the robot
stop_standalone("my_robot",frequency)

```

A useful complement to this webpage is the documentation of the [o80_example](http://people.tuebingen.mpg.de/mpi-is-software/o80/docs/o80_example/index.html) package, which provides a concrete implementation of an o80 project.

The github source code of o80 is [here](https://github.com/intelligent-soft-robots/o80).


[![DOI](https://joss.theoj.org/papers/10.21105/joss.02752/status.svg)](https://doi.org/10.21105/joss.02752)
[![build](https://raw.githubusercontent.com/MPI-IS-BambooAgent/sw_badges/master/badges/plans/o80/build.svg?sanitize=true)](url)
[![documentation](https://raw.githubusercontent.com/MPI-IS-BambooAgent/sw_badges/master/badges/plans/o80/doc.svg?sanitize=true)](url)
[![unit_tests](https://raw.githubusercontent.com/MPI-IS-BambooAgent/sw_badges/master/badges/plans/o80/unit%20tests.svg?sanitize=true)](url)


# o80
Synchronization of c++ processes via a custom python API.

See : the [documentation](http://people.tuebingen.mpg.de/mpi-is-software/o80/docs/o80/index.html), an [example project](https://github.com/intelligent-soft-robots/o80_example) and the [installation instructions](http://people.tuebingen.mpg.de/mpi-is-software/o80/docs/o80/doc/02.installation.html)

See also: 
- [o80 for the frank panda robot](https://github.com/Data-Science-in-Mechanical-Engineering/franka_o80)
- [o80 for the Pneumatical Artificial Muscle Robot](https://intelligent-soft-robots.github.io/pam_documentation/index.html)

# Citing o80

Berenz et al., (2021). The o80 C++ templated toolbox: Designing customized Python APIs for synchronizing realtime processes. Journal of Open Source Software, 6(66), 2752, https://doi.org/10.21105/joss.02752

# Author
Vincent Berenz, Max Planck Institute for Intelligent Systems, Empirical Inference Department

# Unit-tests
The unit tests for o80 are with the [o80_example](https://github.com/intelligent-soft-robots/o80_example) package
---
title: 'The o80 C++ templated toolbox: Designing customized Python APIs for synchronizing realtime processes'
tags:
  - Python
  - C++
  - processes
  - robotics
  - shared memory
authors:
  - name: Vincent Berenz^[corresponding author]
    affiliation: 1 
  - name: Maximilien Naveau
    affiliation: 1
  - name: Felix Widmaier
    affiliation: 1
  - name: Manuel Wüthrich
    affiliation: 1
  - name: Jean-Claude Passy
    affiliation: 1
  - name: Simon Guist
    affiliation: 1
  - name: Dieter Büchler
    affiliation: 1
affiliations:
 - name: Max Planck Institute for Intelligent Systems, Tübingen, Germany
   index: 1
date: 05 July 2020
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
#aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
#aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Statement of need

o80 (pronounced "oh-eighty") is software for synchronizing and organizing message exchange between (realtime) processes via simple customized Python APIs. Its target domain is robotics and machine learning. Our motivation for developing o80 is to ease the setup of robotics experiments (i.e., integration of various hardware and software) by machine learning scientists. Such setup typically requires time and technical effort, especially when realtime processes are involved. Ideally, scientists should have access to a simple Python API that hides the lower level communication details and simply allows the sending of actions and receiving of observations. o80 is a tool box for creating such API.

o80 is, in some aspects, similar to ROS's actionlib [@actionlib], as both allow a process to create and monitor execution of commands executed by another process, with Python and C++ interoperability. A particular feature of o80 is its support for queues of command and the ability to request the server to automatically perform linear interpolation between them. o80 also introduces new synchronization methods (see "bursting mode" in the next section). Unlike actionlib, o80 does not support network communication as it expects the processes it orchestrates to run on the same computer.

# Overview

For implementing synchronization, o80 organizes two types of processes:

- A server process encapsulates an instance of the o80 back-end and an instance of a driver. At each iteration, the back-end instance computes for each actuator the desired state to be transmitted to the hardware via the driver. The back-end also reads sensory information from the driver. Typically, the server process is programmed in C++ in consideration of realtime needs.
- A client process encapsulates instances of the o80 front-end, which provides: 1) an interface to send commands to the back-end (e.g., requests to compute desired states), 2) methods for querying sensory information, and 3) methods for synchronizing the client with the server process.

In the background, back-end and front-end(s) communicate by exchanging serialized object instances via a interprocess shared memory. Serialization is based on the cereal library [@cereal], and the shared memory is based on the Boost interprocess library [@boost].

o80 is templated over actuator states and driver, and may therefore support a wide range of systems.

o80's core API supports:

- methods for specifying via commands either full desired state trajectories or partial trajectories relying on interpolation.
- interpolation methods based on specified duration, speed, or number of server iterations.
- methods for either queuing or interrupting trajectories.
- frontend methods for setting commands that are either blocking or non blocking.
- frontend methods for retrieving sensory data that request the latest available information, sensor information corresponding to a specific past server iteration, or the full history of sensory information.
- client processes and the server processes that run asynchronously or synchronously.

Synchronization methods may be based on a fixed desired frequency set for the server, or may be set up by the client process ("bursting mode").

The o80 library may be considered complex, as it is versatile and can be used in the context of a wide range of control strategies. Yet, the objective of o80 is to provide robot users with a simple Python API. For this purpose, o80 provides convenience functions for generating application tailored Python bindings. The expected usage is that an expert developer uses o80 to design the Python API that will hide, from end users, the complexities related to interprocess communication and synchronization. Scientists may then use this simplified API to design new robotic experiments. In this sense, o80 aims to be a toolbox for creating customized Python APIs. Generation of Python bindings via the o80 API is based on pybind11 [@pybind11].

# Modular Implementation

o80 is based on open source packages that are maintained by the Max Planck Institute for Intelligent Systems, and that may be reused in other contexts. These packages are available on the GitHub organizations "intelligent-soft-robots", "machines-in-motion", "mpi-is", and "open dynamic robot initiative". Examples of such packages are:

- synchronizer: a library for synchronizing processes
- shared memory: a wrapper over the Boost interprocess library that makes exchange of serialized data over an interprocess shared memory trivial
- time series: a templated circular buffer with time stamps, supporting multiprocess access and synchronization

The complete list, the sources, the binaries, as well as the documentation of theses packages can be found online [@corerobotics].


# Examples of usage

## Integration with SL

An instance of the o80 backend has been integrated into the SL realtime library [@sl] used for the control of the Apollo manipulator [@apollo]. This allows scientists to program robot behavior using a simple Python interface running at a low non-realtime frequency that synchronizes with the realtime higher frequency C++ control loop of the robot. This interface was used, for example, for the experiments described in [@icsds].

## HYSR training

o80 has been used in the context of reinforcement learning applied to real robotic hardware. [@pam] describes a setup in which a robot arm driven by pneumatic artificial muscles autonomously learns to play table tennis using an hybrid sim and real training approach (HYSR), i.e, performing real hardware motions to interact with simulated balls. o80 was used in this context to:

- provide a Python API that has been integrated into a Gym environment to provide an interface to reinforcement algorithms [@gym]
- setting up the synchronization between the real robot control and the MuJoCo simulator [@mujoco] used for HYSR
- setting up asynchronous processes for live plotting of the robot state

The code and documentation of this project are available as open source online [@pamsource].

## Franka Emica Panda Robot System

o80 drivers for the Franka Emica Panda Robot System are also available [@o80_franka].

# References



# Bursting mode

## Overview

It is possible to start the backend in bursting mode.

```python
import o80_robot

segment_id = "id1"
frequency = -1 # in Hz
bursting_mode = True # !
driver_arg = [1,1]

o80_robot.start_standalone(segment_id,
                           frequency,
                           bursting_mode,
			   *driver_args)
```

When *not* using the bursting mode, the start_standalone function spaws a process that will iterate (i.e interact with the driver) at the specified frequency.

When using the bursting mode, the specified frequency is ignored. Instead, the process will iterate only when a frontend requests it to do so.

```python
frontend = o80_robot.FrontEnd(segment_id)
frontend.burst(10)
```

The code above requests the standalone process to perform exactly 10 iterations, and then wait. The burst function returns after the 10 iterations have been perfomed. The 10 iterations are performed as fast as the driver software allows.

Similarly to the functionos latest, pulse and pulse_and_wait; the burst method returns an observation.

```python
current_iteration = frontend.latest().get_iteration()
observation = frontend.burst(10)
# new_iteration is current_iteration+10
new_iteration = observation.get_iteration()
```

The burst function also send the current list of commands.

```python
target_value = 1000
frontend.add_command(0,o80_robot.State(target_value),
                     o80.Iteration(2000,True,True),
                     o80.Mode.OVERWRITE)
observation = frontend.burst(2000)
# value will be 1000
value = observation.get_desired_states().get(0).get()
```

## Usage

The typical usage of the bursting mode is interaction with a simulator. Via the frontend, it is possible to create commands (or queues of commands), and then having them executed as fast the simulator allows.

Another usage is synchronization between process:

```python
obs1 = frontend1.burst(3)
frontend2.add_command(0,
                      obs1.get_observed_value().get(0),
		      o80.Mode.OVERWRITE)
frontend2.burst(1)
```


# Support

## mail

For support and question, please contact [Vincent Berenz](https://is.mpg.de/person/vberenz).

## issues

If you have ideas for improvement, or would like to report a bug, you may create issues on [github](https://github.com/intelligent-soft-robots/o80).# Developers

This guide is here for those who would like to upgrade/extend o80, or one of its dependencies.

For updating the code, the colcon workspace installation is required, as described [here](file:///home/vberenz/Workspaces/pam_colcon/workspace/install/o80/share/o80/docs/sphinx/html/doc/02.installation.html#via-colcon-workspace).

The workflow is based on these tools: github, treep, ament, colcon, pybind11 and gtest. This page will give a brief usage guide for them, before describing the guidelines to follow.


## Tools

### Github

The source repository is hosted in the "Intelligent Soft Robots" [domain](https://github.com/intelligent-soft-robots) of github.
As treep (see below) uses ssh, you need to register your ssh key to github (see [here](https://docs.github.com/en/github/authenticating-to-github/adding-a-new-ssh-key-to-your-github-account))

### Treep

To get o80 compiled and running, several git repositories must be cloned.
As it would be annoying to clone them one by one, we use the project manager [treep](https://pypi.org/project/treep/), which allows to clone several repositories at once.

#### installation

```bash
pip install treep
```

#### configuration repository

treep needs a configuration folder. We use [treep_isr](https://github.com/intelligent-soft-robots/treep_isr.git):

```bash
mkdir Software # or another folder name
cd Software
git clone https://github.com/intelligent-soft-robots/treep_isr.git
```

In the folder treep_isr, you will find the file *repositories.yaml* which provides repositories names, url and default branch (see also configuration.yaml). Once treep_isr has been cloned, you may see for example the list of repositories:

```bash
treep # treep without argument prints some help
treep --repos
```

Projects (i.e. repositories grouped together) are described in *projects.yaml*.
To see the list of projects:

```bash
treep --projects
```

To see information about a specific projects:

```bash
treep --project O80
```

#### cloning

You can use treep to clone a project (i.e. cloning all the repositories of the project):

```bash
treep --clone O80
```

#### status

The status argument is very useful to see the (git) status of each cloned repos (i.e. if they are behind origin, beyond origin, have modified files, etc)

```bash
treep --status
```

#### pull 

You can use treep to update all the repositories

```bash
treep --pull
```


#### calling from anywhere

Note that you do not need to be in the Software folder to call treep, you may call treep for any of its subfolder, e.g.

```bash
cd Software/workspace/src/pam_interface
treep --status # still provides the status of all repositories.
```

### Ament packages

(Almost) all repositories contain the code for an ament package. An ament package is a "extended" CMake package, i.e. ament provides some cmake functionalities meant to make developers life easier. The documentation can be found [here](https://docs.ros.org/en/foxy/Guides/Ament-CMake-Documentation.html).

An ament package contains typically:

- C++ code, in the include/package_name and the source folders

- python wrappers over this C++ code in the srcpy folder

- native python code in the python folder

- unit tests in the tests folder

- configuration files

- a demos folder, containing source code providing example of usage of the code provided by the package

- a package.xml file listing dependencies

- a CMakeLists.txt file containing the CMake commands required for compilation.

The [context](https://github.com/intelligent-soft-robots/context) and [pam_interface](https://github.com/intelligent-soft-robots/pam_interface) repositories provide examples of a packages containing most of the above.

### Colcon

#### Overview

After cloning a workspace using treep, the workspace/src folder contains a collection of ament packages.
Colcon is the tool used for compiling all these packages:

```bash
cd Software/workspace
colcon build # will attempt to compile the whole workspace
```

If the compilation succeed, the workspace folder will contain the "install" folder which contains the binaries and python packages that result from the compilation.

To be able to "use" the code (i.e. adding executables to path, adding the python package to the python path, adding the libraries to the ld path ...), one must "source" the file Software/workspace/install/setup.bash:

```bash
cd Software/workspace/install
source setup.bash
```

This operation must be done in all new terminal. It is therefore recommanded to add to the ~/.bashrc file (which is sourced everytime a new terminal is opened):

```bash
echo "- sourcing workspace"
source /path/to/Software/workspace/install/setup.bash
```

#### Compiling only one package

To compile only one specific package:

```bash
cd Software/workspace
colcon list # list all packages in the workspace
colcon build --packages-up-to packagename
```

This will compile only the selected package and its required dependencies.


### Pybind11

Most of the code is in C++, but made available as python package.

We use [pybind11](https://pybind11.readthedocs.io/en/master/?badge=master) for creating wrappers over c++ code.

By convention, the c++ binding code should be in the srcpy folder of the package, and in the file wrappers.cpp.

As you can see in this [example](https://github.com/intelligent-soft-robots/context/blob/master/srcpy/wrappers.cpp), creating binders for c++ classes can be trivial.

For example, this wrapper code:

```cpp
   pybind11::class_<Ball>(m, "Ball")
        .def(pybind11::init<int>())
        .def("update", &Ball::update)
        .def("get", &Ball::get);
```

binds this c++ class [ball](https://github.com/intelligent-soft-robots/context/blob/master/include/context/ball.hpp).

The following calls has to be added in the CMakeLists.txt file:

```cmake
pybind11_add_module(context_py srcpy/wrappers.cpp)
target_link_libraries(context_py PRIVATE context ${catkin_LIBRARIES})
set_target_properties(context_py PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY ${CATKIN_DEVEL_PREFIX}/${CATKIN_GLOBAL_PYTHON_DESTINATION}
  OUTPUT_NAME context)
install(TARGETS context_py DESTINATION ${CATKIN_GLOBAL_PYTHON_DESTINATION})
```

(replace "context" by the name of the catkin package).

Note that in the [wrappers.cpp](https://github.com/intelligent-soft-robots/context/blob/master/srcpy/wrappers.cpp#L12
) file, the package name must also be used as first argument:

```cmake
PYBIND11_MODULE(context, m)
```

The above results in the creation of the python package "context", which can be used after compilation and sourcing of the workspace, simply in a terminal:

```bash
python
```

then:

```python
import context
```

### Unit-tests

We use google tests for unit tests. See the [tests folder](https://github.com/intelligent-soft-robots/context/tree/master/tests) of the context package for a simple example to follow.

To run all the unit tests included in the package, run:

```bash
cd Software/workspace
colcon test --event-handlers console_direct+
```

To run the unit test of a selected package:

```bash
colcon test --packages-select package_name --event-handlers console_direct+	
```

## Repositories

Apart of the o80 ament package, the treep project O80 also clone its dependencies, which are listed here:

- real_time_tools: API for spawning thread (realtime threads if compiled on a rt-preempt computer)

- shared_memory: wrappers over the [boost interprocess library](https://www.boost.org/doc/libs/1_63_0/doc/html/interprocess/sharedmemorybetweenprocesses.html), aiming a making easier the usage of it.

- signal_handler: detection of ctrl+c in programs

- time_series: interprocess implementation of a circular buffer

- synchronizer: synchronization between processes made easier

- serialization_utils: convenience wrappers over the [cereal serialization library](https://uscilab.github.io/cereal/).

- o80_example: o80's canonical example

Documentation for these packages may be found [here](http://people.tuebingen.mpg.de/mpi-is-software/corerobotics/). Note that all these packages are installed along o80 as explained [here](http://people.tuebingen.mpg.de/mpi-is-software/o80/docs/o80/index.html).

##  Pull requests on Github

We accept pull requests on github.

- Before starting to work on a repository, create a branch: your_name/what_you_will_work_on.

- Please follow these [guidelines](https://machines-in-motion.github.io/code_documentation/ci_example_cpp/docs/doxygen/html/coding_guidelines_1.html)

- Format your c++ code by running the executable *mpi_cpp_format* (installed along o80)

- Format your python code using [black](https://github.com/psf/black)

- Update the unit-tests. Note that the o80 unit tests are in the [o80_example](https://github.com/intelligent-soft-robots/o80_example) package.

In your pull-request:

- describe what functionality you added or updated, and why

- describe how you know the code you wrote works (e.g. running unit tests and/or demos)

- confirm you followed the guidelines and used the formatting tools






# Integrating your hardware

 [o80 example](https://github.com/intelligent-soft-robots/o80_example) provides an example of (pseudo toy) hardware integration. Here are the main steps to follow.

## Step1: Robot driver

To integrate your hardware, you need to code its corresponding Driver, which should be sublcass of [o80::Driver](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/driver.hpp).

A driver simply provides functions to communicate with the hardware, i.e. to send "action" (arbitrary user templated class IN)  and read sensors (arbitrary user templated class OUT).

Here is a toy examples of a driver:

- [driver.hpp](https://github.com/intelligent-soft-robots/o80_example/blob/master/include/o80_example/driver.hpp)
- [driver.cpp](https://github.com/intelligent-soft-robots/o80_example/blob/master/src/driver.cpp)


## Step2: o80 Standalone

An o80 [Standalone](https://github.com/intelligent-soft-robots/o80/blob/vberenz/doc/include/o80/driver.hpp) is an object that will wrap an instance of Driver so to make it compatible with the o80 API. 

A Standalone is templated over IN and OUT (the templates of the Driver). It is also templated over State and ExtendedState. Explanations:

### State

While IN is the input to a robot driver, State represents the desired state of *one* actuator of the robot.  For example, the State can be the position of a joint. 

State is a user developed class which should inheritate [o80::State](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/state.hpp).

Standalone is templated over State because o80 allows to send desired states command to each actuator independantly. The Standalone will receive commands related to desired States, and will convert them to instance of IN before forwarding them to the driver. The user implemented Standalone class will need to implement a [convert method](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/standalone.hpp#L126) for converting instances of State into instances of IN.

A State class does not only declare the desired state of a joint, but also [methods](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/state.hpp#L43) defining how desired states interpolates.

 This will allow for example o80 user API to support duration command method, which will request the desired state of an actuator to reach, starting from its current desired state, a target desired state of a specified duration (this API is presented in a later section, but [here](https://github.com/intelligent-soft-robots/o80_example/blob/master/demos/duration_commands.py#L30) a preview). 

The o80::State base class implements simple [linear interpolation methods](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/state.hxx) for basic types (int, double, float). If the user State encapsulate more complex data, and/or if the interpolation should not be linear, the "intermediate_state" methods should be overridden. 

If an actuator state consists of a boolean, [o80::BoolState](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/bool_state.hpp)  can be used.

If the actuator state is empty (i.e. the hardware is a sensor that takes no input), [o80::VoidState](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/void_state.hpp) can be used.

For example, the toy Joint class is a valid State class:

- [joint.hpp](https://github.com/intelligent-soft-robots/o80_example/blob/master/include/o80_example/joint.hpp)
- [joint.cpp](https://github.com/intelligent-soft-robots/o80_example/blob/master/src/joint.cpp)

### ExtendedState

The template OUT of a Driver is the arbitrary output of a robot.

A o80 [Observation](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/observation.hpp) output of an o80 Standalone is more structured, and has three parts:

- An instance of [States](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/states.hpp), which encapsulates the current desired state of each actuator
- An instance of [States](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/states.hpp), which encapsulates the current observed state of each actuator
- the extended state (an arbitrary user class)

During runtime, the o80 Standalone will call the [get](https://github.com/intelligent-soft-robots/o80/blob/vberenz/doc/include/o80/driver.hpp#L16) method of the driver and retrieve an instance of OUT. It will then need to "convert" this OUT instance into an instance of [o80::States](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/states.hpp). This instance of o80:States will encapsulate the current observed state of each actuator of the robot. The Standalone user class should implement this convert function.

It may be that the instance of OUT returned by the get function of the driver encapsulate data other than the observed state of the robot actuator. This information may be encapsulated in o80 Observation via the extended state. During runtime, the standalone instance will call its ["enrich_extended_state"](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/state.hpp) function to convert the instances of OUT into instances of ExtendedState. As OUT, EXTENDED_STATE is an arbitrary user class that templates Standalone. Note: OUT and EXTENDED_STATE may be the same class, in which case the "enrich_extended_state" function becomes trivial. If there are no sensory information beyond the actuators state, a [void extended state](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/void_extended_state.hpp) can be used.

### Summary

To implement an o80 Standalone, one must implement:

- a Driver class inherating [o80::Driver](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/driver.hpp). This class imlements the communication with the hardware.
- An OUT class and a out class that templates this driver class. Instances of these class are used by the driver to communicate with the hardware.
- A State class inherating [o80::State](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/state.hpp). An instance of State represente the state of an actuator.
- An ExtendedClass, which is an arbitrary class which will encapsulate robot's information other than states. Often, ExtendedClass can also be the OUT class.
- A Standalone, inherating [o80::Standalone](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/standalone.hpp) which is templated over the Driver class, the State class and the ExtendedClass.

## Step3 : library

The documentation above request to create 5 classes : In, Out, State, Driver and Standalone. Once these class declared, they shall be compiled in a new library that links with robot_interfaces and o80, as shown [here](https://github.com/intelligent-soft-robots/o80_example/blob/master/CMakeLists.txt#L30).

## Step4 : python bindings

[pybind11](https://pybind11.readthedocs.io/en/stable/) can be used to create python bindings over the library above.

o80 provides [helper functions](https://github.com/intelligent-soft-robots/o80/blob/master/include/o80/pybind11_helper.hpp) which will generate the bindings automatically. See the [example](https://github.com/intelligent-soft-robots/o80_example/blob/master/srcpy/wrappers.cpp).   
The macro for the creation of the pybind11 module need to be added to the [CMakeLists.txt](https://github.com/intelligent-soft-robots/o80_example/blob/master/CMakeLists.txt#L42).

The [o80::create_python_bindings](https://github.com/intelligent-soft-robots/o80/blob/e050f1ae16b47c4000f85fb237ded7835b3b3daa/include/o80/pybind11_helper.hpp#L73) function is templated over the Driver and the Standalone class, but other arguments (... DriverArgs). These correspond to the signature of the Driver's constructor. For example, this driver's constructor takes two double as arguments. Thus, the corresponding [create_python_bindings function](https://github.com/intelligent-soft-robots/o80_example/blob/c83398d19cd8bb69930fabe42b7cabf1a1d8fa38/srcpy/wrappers.cpp#L6) is templated over two doubles.

Note: important detail. The name of the pybind11 library compiled and binded in the CMakeLists.txt has to match the name of the library passed in the PYBIND11_MODULE macro (as in [here](https://github.com/intelligent-soft-robots/o80_example/blob/c83398d19cd8bb69930fabe42b7cabf1a1d8fa38/CMakeLists.txt#L45) and [here](https://github.com/intelligent-soft-robots/o80_example/blob/c83398d19cd8bb69930fabe42b7cabf1a1d8fa38/srcpy/wrappers.cpp#L4).

## Step5: compilation

For the moment, o80 and robot_interfaces supports only catkin and catkin tools, and the installation instruction should have resulting in the update of your catkin workspace (or the creation of a new catkin workspace), and therefore "catkin_make" or "catkin build" can be used for compilation.

Note that only python3 is supported for the bindings, and therefore the path to the correct python executable must be provided for compilation, for example:

```bash
catkin_make -DPYTHON_EXECUTABLE=/usr/bin/python3 install
```

# Installation

o80 has been tested exclusively on ubuntu 20.04 and ubuntu 18.04.

## From binaries

o80's binaries are available only for ubuntu18.04/python3.6 and ubuntu20.04/python3.8.

To install, first select a version by visiting : [http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/](http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/) or [http://people.tuebingen.mpg.de/mpi-is-software/o80/older/](http://people.tuebingen.mpg.de/mpi-is-software/o80/older/). Then, for example:

```bash
wget http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/o80_ubuntu20.04_py3.8_1.0.tar.gz
tar -zxvf ./o80_ubuntu20.04_py3.8_1.0.tar.gz
sudo ./apt-dependencies
# for global installation
#sudo ./pip3-dependencies
# for user installation or currently activated virtual environment installation
./pip3-dependencies
# configuration. By default, this will
# install libraries and includes to /usr/local
# and python packages to the active python (including currently
# activated virtual environment, if any). To overwrite the default,
# for example: 
# ./configure --prefix=/usr/ --pythondir=/usr/local/lib/python3.8/dist-packages 
./configure
sudo make install
sudo ldconfig
```

## From source

### From tar ball

Select a version by visiting : [http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/](http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/) or [http://people.tuebingen.mpg.de/mpi-is-software/o80/older/](http://people.tuebingen.mpg.de/mpi-is-software/o80/older/). Then, for example:

```bash
wget http://people.tuebingen.mpg.de/mpi-is-software/o80/latest/o80_source.tar.gz
tar -zxvf ./o80_source.tar.gz
sudo ./apt-dependencies
# for global installation
#sudo ./pip3-dependencies
# for user installation or currently activated virtual environment installation
./pip3-dependencies
# configuration. By default, this will
# install libraries and includes to /usr/local
# and python packages to the active python (including currently
# activated virtual environment, if any). To overwrite the default,
# for example: 
# ./configure --prefix=/usr/ --pythondir=/usr/local/lib/python3.8/dist-packages 
./configure
make
sudo make install
sudo ldconfig
```

### Via colcon workspace

[Colcon](https://colcon.readthedocs.io/en/released/) is the built system of [ROS2](https://docs.ros.org/en/foxy/index.html).
The instructions below will result in the setup of a colcon workspace. Possibly, if you would like to use o80 in a ROS2 project, you may copy the cloned packages to an existing workspace.

#### Adding your ssh key to github

See: [github documentation](https://help.github.com/en/github/authenticating-to-github/connecting-to-github-with-ssh).

All the following instructions assume your ssh key has been activated, i.e.:

```bash
# replace id_rsa by your key file name
ssh-add ~/.ssh/id_rsa
```

#### Installing the dependencies:

```bash
apt install -y python3-pip cmake libcereal-dev \
               libboost-all-dev libgtest-dev \
	       libeigen3-dev libedit-dev \
	       libncurses5-dev freeglut3-dev \
	       libxmu-dev doxygen
```

```bash
pip3 install colcon-common-extensions treep \
             fyplot pyside2 empy \
	     catkin-pkg sphinx breathe
```

#### Cloning the repositories

Creating a folder and cloning the treep configuration:

```bash
mkdir Software # you may use another folder name
cd Software
git clone git@github.com:intelligent-soft-robots/treep_isr.git
```

Cloning all the required repositories:

```bash
treep --clone O80
```

#### Compilation

```bash
cd /path/to/Software
cd workspace
colcon build
```

This will result in a "install" folder containing the compiled binaries

### Activating the workspace

In each new terminal, the workspace needs to be sourced:

```bash
source /path/to/Software/install/setup.bash
```

Possibly, you may want to add the line above to the ~/.bashrc file (so that each new terminal source the workspace automatically).


## Checking installation

In a python3 terminal:

```python
import o80
```

## Running the demos

The o80_example package provides usage examples of o80 running on a dummy robot, see [here](http://people.tuebingen.mpg.de/mpi-is-software/o80/docs/o80_example/doc/02_demos.html).


# Troubleshooting

## Processes hanging

### clearing the shared memory

A recurring issue when using o80 is that processes may hang if trying to connect to a deprecated shared memory. A deprecated shared memory occurs when a process is not cleanly exited (i.e. the destructor methods are not properly called).

This can be solved by called the o80::clear_shared_memory method using the suitable segment_id as argument.

o80::clear_shared_memory has not down side, and it is even good manner to call this function before a new standalone is instantiated.

```python

segment_id = "id0"
o80.clear_shared_memory(segment_id)
o80_robot.start_standalone(segment_id, ...)

```

If the clear_shared_memory does not work, you may try to delete the content of the folder /dev/shm. Note that this may also delete shared memory files used by process not related to o80.

### instantiating in the right order

A FrontEnd targeting a segment_id must always be started *after* a corresponding standalone or BackEnd has been instantiated. If the standalone or BackEnd is destroyed before the FrontEnd, the instance of FrontEnd will hang.




# Overview

o80 (pronounced 'oh eighty') is a tool for synchronizing processes while organizing exchange of information between them.
The information exchanged are commands for computing (robotic) desired states and observations, where state and observation are (almost) arbitrary user declared and defined classes (o80 classes are templated over them).

Here is a basic usage of o80, using python:

```python

# starts a process running robot control at
# 1000Hz
frequency = 1000
start_standalone("my_robot",frequency)

# creating a frontend for communicating with the robot
frontend = FrontEnd("my_robot")

# creating a command requesting the desired state of
# the 0th actuator to reach the value of 100 in 3 seconds
frontend.add_command(0,State(100),Duration.seconds(3),QUEUE)
# same for the 1th actuator
frontend.add_command(1,State(100),Duration.seconds(3),QUEUE)

# send the commands. The desired states of 0th and 1st dof
# of the robot interpolates from their initial values to 100.
# Once this is achieved, an observation is returned
observation = frontend.pulse_and_wait()

# reading the state of the robot 
state0 = observation.get_observed_states().get(0).get()
state1 = observation.get_observed_states().get(1).get()

# stopping the robot
stop_standalone("my_robot",frequency)

```

In the example above, the standalone process spawns a realtime process corresponding
to c++ (realtime) code. The frontend communicates under the hood with this process via
a realtime shared memory.

The frontend provides the user with an "add_command" which is flexibe (i.e. has many convenient overloads). For Example:

```python

# the desired state will reach the value of 100 at 3 units per second
frontend.add_command(0,State(100),Speed.per_second(3),QUEUE)

# after this, the desired state will reach the value of 150 at the 5000th 
# robot's control iteration
frontend.add_command(0,State(150),Iteration(5000),QUEUE)

# the desired states will increase by 5 for each of the following iteration
frontend.add_command(0,State(155),QUEUE)
frontend.add_command(0,State(160),QUEUE)
frontend.add_command(0,State(170),QUEUE)

# sending the commands to the robot and returning immediately
frontend.pulse()

# letting the robot run for 1 second, then cancelling all running
# commands and requesting the desired state to decrease to 0 in
# 5 seconds
time.sleep(1.0)
frontend.add_command(0,State(0),Duration.seconds(5),OVERWRITE)
frontend.pulse()

```

o80 classes are templated over State. In the above, State is a python wrapper over an arbitrary
c++ "State" class. For example, it can be for an actuator of your robot the desired position and velocity.

o80 can also be used to synchronize processes. Here an example that has a simulated robot mirroring a
a real robot:

```python

robot = Frontend("real_robot")
simulation = Frontend("simulation")

# sending a trajectory to one of the
# robot actuator
starting_iteration = robot.read().get_iteration()
end_iteration = starting_iteration+len(values)
for value in values:
    robot.add_command(0,State(value),QUEUE)
robot.pulse()

# having the simulated robot mirroring the real robot
while True:
    observation = robot.wait_for_next()
    observed_states = observation.get_observed_states()
    for actuator in range(nb_actuators):
        state = observed_states.get(actuator)
        simulation.add_command(actuator,state,OVERWRITE)
    simulation.pulse()
    if observation.get_iteration()>=end_iteration:
        break

```

This is the gist of o80. For a fully commented concrete example / tutorial, you may check: [demos](https://github.com/intelligent-soft-robots/o80_example/tree/master/demos).
# Usage

Once the catkin workspace sourced, o80 may be used via its python API. Concrete examples of usage of this API is provided in the [o80 example](https://github.com/intelligent-soft-robots/o80_example) package.

The code snippets below assumes the python binded library is called "o80_robot", and that the State class encapsulates only one integer. 

## Starting a standalone

```python
import o80_robot

segment_id = "id1"
frequency = 1000 # in Hz
bursting_mode = False
driver_arg = [1,1]

o80_robot.start_standalone(segment_id,
                           frequency,
                           bursting_mode,
			   *driver_args)
```

The above starts an o80 Standalone that will run at 1000Hz, i.e. it will run 1000 iterations per second, and at each iteration will call the set() and get() method of the encapsulated driver.
The argument is an arbitrary id (segment_id), the frequency, the bursting mode (later explained, set to False in doubt) and the arguments required for instantiating the driver.
The standalone is spawned in a separated (c++) thread.

## Starting a frontend

```python
import o80
import o80_robot

segment_id = "id1"

frontend = o80_robot.FrontEnd(segment_id)
```

This starts a frontend, i.e. the instance of an object able to communicate with the standalone (i.e. send commands, read observation). The segment_id passed should be the same as the one used when starting the standalone.

**The frontend and the standalone can be instantiated in two different scripts**. Just: the standalone must be started first. The frontend and the standalone will communicate under the hood via a realtime shared memory. 

## Buffering commands

```python
target_value = 1000
# actuator 0
frontend.add_command(0,o80_robot.State(target_value),o80.Mode.OVERWRITE)
# actuator 1
frontend.add_command(1,o80_robot.State(target_value),o80.Mode.OVERWRITE)
```

The above creates a command requesting the desired state of actuator 0 to take the value 1000 as soon as possible (o80.Mode.OVERWRITE explained later). 
It also creates a command requesting the same desired state value of actuator 1.
These commands are not applied yet, they are just buffered in the frontend.

## Sending commands

```python
frontend.pulse()
```

The buffer of command is purged, and all the commands are sent to the Standalone which execute them. Consequently, the desired state of actuators 0 and 1 get to the value 1000.

## Time, Speed and Iteration commands

The commands passed above requested the desired states of the actuators to change immediately upon execution of the command.
One may request instead for the desired states to interpolate from their current value to their target value over several iterations.
Several overload of the add_command function may be used:

### Duration commands

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
```

The command above requests the desired state value to interpolate from its current value to the target value over two seconds.

### Speed commands

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),o80.Speed.per_millisecond(1),o80.Mode.OVERWRITE)
```

The command above requests the desired state value to interpolate from its current value to the target value at the speed of 1 unit per milliseconds.

### Iteration commands

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Iteration(5000),o80.Mode.OVERWRITE)
```

The command above request to interpolate from its current value so that reaching the target value at the 5000th iteration of the standalone. 
Note that it does **not** mean reaching the target value over 5000 iterations. When the standalone is started, it starts iterating at its specified frequency, and keeps count of its iteration number (starting at 0). If the current iteration number is 10 when the command is started, the desired state will iterate over 4990 iteration. If the iteration number is 4950, it will interpolate over 50 iterations.

To specify the number of iteration to interpolate over:

```python
actuator = 0
target_value = 1000
relative = True
reset = True
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Iteration(5000,relative,reset),o80.Mode.OVERWRITE)
```

The above request to reach to target value at the 5000th *relative* iteration number, i.e. the iteration number relative to the last command for which "reset" was True was started.
In this case, this command reset the iteration count, and then considers iteration number relative to this resetted number. It will thus request to interpolate over 5000 iterations.

## Interrupting command

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
frontend.pulse()
time.sleep(1)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
frontend.pulse()
```

The above sends a duration command to the standalone, which starts to execute it. 
But after 1 seconds, the frontend sends another command (note: "pulse" returns immediately, not upon finalized execution of the command). 

Because this second command is sent using the mode "o80.Mode.OVERWRITE", it requests the cancelation of the current command. The standalone interrupt its running command and starts execution of the second command.

This has a different behavior:

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
frontend.pulse()
time.sleep(1)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.QUEUE)
frontend.pulse()
```

The 'QUEUE' mode is used, thus the standalone finishes the first command, then starts execution of the second command.

Note that this is equivalent to:

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),
                     o80.Duration.milliseconds(2000),o80.Mode.QUEUE)
frontend.pulse()
```

## Blocking pulse

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.QUEUE)
frontend.pulse()
```

In the above, the pulse method returns immediately, which allows for the frontend to buffer new commands while other commands are executed.

Alternatively, to block the script until these two commands are completed:

```python
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.QUEUE)
frontend.pulse_and_wait() # !
```

Instead of returning once the commands have been finalized, one may return at a specific standalone iteration number

```python
frontend.pulse(o80.Iteration(5000)) 
```

This function will block until the 5000th iteration of the standalone has been reached.

## Reading observations

### latest observation

The pulse and pulse_and_wait method do not only purge the command buffer, but they also return an observation.

```python
observation = frontend.pulse()
```

Alternatively, the latest method also returns an observation (but does not purge the command buffer)

```python
observation = frontend.latest()
```

Observation encapsulate these data:

```python
observed_states = observation.get_observed_states()
desired_states = observation.get_observed_state()
iteration_number = observation.iteration()
frequency = observation.get_frequency()
extended = observation.get_extended_state()
```

observed_states and desired_states are instances of States. To get the State of the 0th actuator:

```python
state = observed_states.get(0)
value = state.get()
```

*iteration number* is a integer corresponding to the standalone iteration number at the time this observation was created.

*frequency* is the frequency of the standalone as observed at the corresponding iteration. 

*extended* is an instance of the "extended class" templating the Standalone. 


### next observation

```python
frontend.reset_next_index()
observation1 = frontend.wait_for_next()
iteration1 = observation1.get_iteration()
observation2 = frontend.wait_for_next()
# iteration2 will be iteration1+1
iteration2 = observation2.get_iteration()
```

The method wait_for_next returns the "next" observation since the previous call to wait_for_next, or wait for such observation to be generated (by the standalone process).

### observation history

```python
history = frontend.get_latest_observations(1000)
```

History is a list of instances of Observation of size 1000, corresponding to the last 1000 observations writen by the standalone (one per iteration).

Alternatively:
 
```python
starting_iteration = frontend.latest().get_iteration()
actuator = 0
target_value = 1000
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.OVERWRITE)
target_value = 500
frontend.add_command(actuator,o80_robot.State(target_value),o80.Duration.milliseconds(2000),o80.Mode.QUEUE)
frontend.pulse_and_wait() 
history = frontend.get_observation_since(starting_iteration)
```


## Using several frontends

Several frontends may connect simultaneously to the same segment_id. This is useful, for example, to create a logging process. For example, in one python executable one may send commands, while in another independant execuble, one may log information related to observations. For example:

```python
while True:
    observation = frontend.wait_for_next()
    print(observation.display())
    time.sleep(0.01)
```

*Important* : while several instances of FrontEnd may run concurrently, only one of them should by used to send commands. Sending commands via several FrontEnds may have unexpected effects. 

## Putting things together

Using the API described above, it is possible for example:

- to generate in python a full trajectory whic specifies a desired state for each Standalone iteration ([code](https://github.com/intelligent-soft-robots/o80_example/blob/master/demos/full_trajectory.py#L20))

- to have a robot replaying another robot actions ([code](https://github.com/intelligent-soft-robots/o80_example/blob/master/demos/delay.py))

- to have a robot mirroring another robot ([code](https://github.com/intelligent-soft-robots/o80_example/blob/master/demos/mirroring.py))
# Integration in user software

So far, we documented the usage of o80 via a standalone, i.e. a python method which spaws a process running the (c++) interaction with the driver.

But one may desire to use o80's frontend API while using their own software architecture.

It is possible to do so by encapsulating an instance of *o80::BackEnd* to the architecture. This could be done either in python and in c++.

## in python

### backend

Similarly than for the FrontEnd and the start_standalone method, the python package wraps the BackEnd class.

```python
import o80_robot
backend = o80_robot.BackEnd(segment_id)
```

A backend instance provides a single method, pulse, which takes the current observed states of the joints as arguments, and returns the current desired states to apply to the robot (as computed based on the queues of commands sent by an eventual frontend).

A call to BackEnd::pulse performs:

- read the queue of commands as generated by an eventual frontend.
- compute for each joint its desired states, according to the current queue of commands.
- creates an instance of o80::Observation, and writes it to the shared memory (i.e. makes it available to FrontEnd::pulse, FrontEnd::latest, etc)
- returns an instance of o80::States, which encapsulate for each joint its desired state.

An instance of BackEnd may be used in the user's code:

```python

import o80_robot

#
# user initialization code
#
# ...
# 
backend = o80_robot.BackEnd(segment_id)

while running:

     # user code compute and extract from sensors
     # the current joint states

     # calling the pulse method of the backend
     # (current_states: list of size number of joints)
     desired_states = backend.pulse(current_states)

     # user code for applying the desired state
     # to the robot
     
     time.sleep(0.001)

```

The pulse function has several overloads. The most complete is:

```python


time_stamp = o80::time_now()
extended_state = MyUserClass()
iteration_update = False
current_iteration = 1000

frontend.pulse(time_stamp,
               current_states,
	       extended_state,
	       iteration_update,
	       1000)

```

- time_stamp will be the time stamp encapsulated in the Observation the BackEnd will write in the shared memory.
- extended_state is an instance of ExtendedState
- iteration_update to False means the BackEnd should not use its internal counter to set the iteration number of the Observation, but rather use the current_iteration argument (iteration_update to False and current_iteration to -1 to use the internal counter)


Other overloads are:

```python

# will use current time as time_stamp, default constructor of ExtendedState, internal iteration counter, and default constructors for desired states.
frontend.pulse()

# same as above, but will use the passed argument for extended state
frontend.pulse(extended_state)

# same as above, but will use the passed argument for desired states
frontend.pulse(desired_states)

```

### bursting mode

To use the bursting mode in a user software, you may update the control loop:

```python

import o80_robot
import o80

#
# user initialization code
#
# ...
# 
backend = o80_robot.BackEnd(segment_id)
burster = o80.Burster(segment_id)

while running:

     # user code compute and extract from sensors
     # the current joint states

     # calling the pulse method of the backend
     # (current_states: list of size number of joints)
     desired_states = backend.pulse(current_states)

     # user code for applying the desired state
     # to the robot

     # ! commented, frequency will be managed by
     # the burster instead
     # time.sleep(0.001)

     burster.pulse()

```

In the above, when its method pulse is called, the loop will hold until a frontend calls its burst function. 

## in c++

An equivalent API is provided in c++.

